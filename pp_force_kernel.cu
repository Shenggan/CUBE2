#include <cuda_runtime.h>
#include <math.h>
#include <math_constants.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <thrust/detail/caching_allocator.h>
#include <thrust/device_ptr.h>
#include <thrust/functional.h>
#include <thrust/gather.h>
#include <thrust/iterator/reverse_iterator.h>
#include <thrust/replace.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <time.h>

#include <vector>

// -------------------------
// Global configuration
// -------------------------
const int ng = 512;
const int nnt = 4;
const int nns = 4;
const int ratio_cs = 4;
const int nc = ng / ratio_cs;
const int np_nc = ratio_cs;
const int nt = nc / nnt;
const int ntt = nt / nns;
const int ngb = 16;
const int ncb = ngb / ratio_cs;
const int nte = nt + 2 * ncb;
const int ns3 = (nnt * nns) * (nnt * nns) * (nnt * nns);

const bool body_centered_cubic = false;
const int image_buffer = 2;
const int tile_buffer = 3;

const int64_t np_image = (nc * np_nc) * (nc * np_nc) * (nc * np_nc) * (1 + body_centered_cubic);
const int64_t np_image_max = np_image * ((int64_t)nte * 1.0 / (int64_t)nt) *
                             ((int64_t)nte * 1.0 / (int64_t)nt) *
                             ((int64_t)nte * 1.0 / (int64_t)nt) * image_buffer;
const int64_t np_tile_max = np_image / (nnt * nnt * nnt) *
                            ((nte * 1.0 / nt) * (nte * 1.0 / nt) * (nte * 1.0 / nt)) * tile_buffer;
const size_t max_vp = (size_t)np_image_max * 3;
const size_t VP_CHUNK_ELEMS = 3 * 80 * 1024 * 1024;

const int64_t ishift = -((int64_t)1 << 15);  // 2^15
const double rshift = 0.5 - ishift;
const double x_resolution = 1.0 / ((int64_t)1 << 16);

const double app = 0.06;

const int ncore = 32;

const int kCellBuildThreads = 256;
const int kUpdateThreads = 256;
const int kPPThreads = 128;
const int kD2HBufferCount = 2;

#define CUDA_CHECK(call)                                                                        \
    do {                                                                                        \
        cudaError_t err = call;                                                                 \
        if (err != cudaSuccess) {                                                               \
            fprintf(stderr, "CUDA error in %s (%s:%d): %s\n", __FUNCTION__, __FILE__, __LINE__, \
                    cudaGetErrorString(err));                                                   \
            exit(EXIT_FAILURE);                                                                 \
        }                                                                                       \
    } while (0)

static inline double elapsed_sec(const struct timespec& t1, const struct timespec& t2) {
    return (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec) * 1.0e-9;
}

__device__ __forceinline__ float F_ra_dev(float r, float apm) {
    float ep = 2.0f * r / apm;
    float f_ra = 0.0f;

    if (apm == 0.0f || ep > 2.0f) {
        f_ra = 1.0f / (r * r);
    } else if (ep > 1.0f) {
        float ep2 = ep * ep;
        float ep3 = ep2 * ep;
        float ep4 = ep3 * ep;
        float ep5 = ep4 * ep;
        float ep6 = ep5 * ep;
        f_ra =
            (1.0f / 35.0f / (apm * apm)) * (12.0f / ep2 - 224.0f + 896.0f * ep - 840.0f * ep2 +
                                            224.0f * ep3 + 70.0f * ep4 - 48.0f * ep5 + 7.0f * ep6);
    } else if (ep > 0.0f) {
        float ep3 = ep * ep * ep;
        float ep4 = ep3 * ep;
        float ep5 = ep4 * ep;
        float ep6 = ep5 * ep;
        f_ra = (1.0f / 35.0f / (apm * apm)) *
               (224.0f * ep - 224.0f * ep3 + 70.0f * ep4 + 48.0f * ep5 - 21.0f * ep6);
    } else {
        f_ra = 0.0f;
    }
    return f_ra;
}

__global__ void build_cell_start_end_kernel(const uint32_t* __restrict__ keys_sorted, int N,
                                            int* __restrict__ cell_start,
                                            int* __restrict__ cell_end) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    uint32_t k = keys_sorted[i];

    if (i == 0 || keys_sorted[i - 1] != k) {
        cell_start[k] = i;
    }
    if (i == N - 1 || keys_sorted[i + 1] != k) {
        cell_end[k] = i + 1;
    }
}

__device__ __forceinline__ void atomicMaxFloatNonNeg(float* addr, float val) {
    // Only valid when val >= 0 and *addr >= 0.
    int* addr_i = reinterpret_cast<int*>(addr);
    int old = *addr_i;
    int v = __float_as_int(val);

    // If current value is already >= val, return directly.
    while (old < v) {
        int assumed = old;
        old = atomicCAS(addr_i, assumed, v);
        if (old == assumed) break;
    }
}

__global__ void pp_force_kernel(const float* __restrict__ xf0, const float* __restrict__ xf1,
                                const float* __restrict__ xf2, int npgrid, int rcp, float pp_range,
                                float app_local, const int* __restrict__ cell_start,
                                const int* __restrict__ cell_end, float* __restrict__ af0,
                                float* __restrict__ af1, float* __restrict__ af2) {
    int cell = blockIdx.x;  // 0..npgrid^3-1
    int cells2 = npgrid * npgrid;
    int k = cell / cells2;
    int j = (cell / npgrid) % npgrid;
    int i = cell % npgrid;

    int stride = (npgrid + 2 * rcp);
    int idx = (k + rcp) * stride * stride + (j + rcp) * stride + (i + rcp);

    int o0 = cell_start[idx];
    int o1 = cell_end[idx];

    float pp2 = pp_range * pp_range;

    for (int outer = o0 + threadIdx.x; outer < o1; outer += blockDim.x) {
        float xo = xf0[outer], yo = xf1[outer], zo = xf2[outer];
        float ax = 0.f, ay = 0.f, az = 0.f;

        for (int kk = -1; kk <= 1; kk++) {
            // #pragma unroll
            for (int jj = -1; jj <= 1; jj++) {
                // #pragma unroll
                for (int ii = -1; ii <= 1; ii++) {
                    int nidx =
                        (k + kk + rcp) * stride * stride + (j + jj + rcp) * stride + (i + ii + rcp);

                    int in0 = cell_start[nidx];
                    int in1 = cell_end[nidx];

                    for (int inner = in0; inner < in1; inner++) {
                        float dx = xf0[inner] - xo;
                        float dy = xf1[inner] - yo;
                        float dz = xf2[inner] - zo;

                        float r2 = dx * dx + dy * dy + dz * dz;
                        if (r2 > 0.0f && r2 < pp2) {
                            float invr = rsqrtf(r2);
                            float r = r2 * invr;
                            float fpp = F_ra_dev(r, app_local) - F_ra_dev(r, pp_range);
                            ax += fpp * dx * invr;
                            ay += fpp * dy * invr;
                            az += fpp * dz * invr;
                        }
                    }
                }
            }
        }

        af0[outer] = ax;
        af1[outer] = ay;
        af2[outer] = az;
    }
}

__global__ void update_v_kernel(const float* __restrict__ af0, const float* __restrict__ af1,
                                const float* __restrict__ af2, const int64_t* __restrict__ ip_local,
                                const uint32_t* __restrict__ values_sorted,
                                float* __restrict__ vp_flat, int N, float mass2, float a_mid,
                                float dt, float* __restrict__ f2_global) {
    extern __shared__ float s[];
    float* sf2 = s;

    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + tid;

    float f2 = 0.f;

    if (idx < N) {
        float ax = af0[idx] * mass2;
        float ay = af1[idx] * mass2;
        float az = af2[idx] * mass2;

        f2 = ax * ax + ay * ay + az * az;

        uint32_t ori = values_sorted[idx];
        int64_t ip = ip_local[ori];

        const float coef = a_mid * dt / (6.0f * CUDART_PI_F);

        float vx = vp_flat[ip * 3 + 0] + ax * coef;
        float vy = vp_flat[ip * 3 + 1] + ay * coef;
        float vz = vp_flat[ip * 3 + 2] + az * coef;

        vp_flat[ip * 3 + 0] = vx;
        vp_flat[ip * 3 + 1] = vy;
        vp_flat[ip * 3 + 2] = vz;
    }

    sf2[tid] = f2;
    __syncthreads();

    for (int offset = blockDim.x / 2; offset > 0; offset >>= 1) {
        if (tid < offset) {
            sf2[tid] = fmaxf(sf2[tid], sf2[tid + offset]);
        }
        __syncthreads();
    }

    // atomically accumulate each block's maxima into 4 global scalars.
    if (tid == 0) {
        atomicMaxFloatNonNeg(f2_global, sf2[0]);
    }
}

struct PPForceWorkspace {
    uint32_t* d_keys;
    uint32_t* d_values;
    float* d_xf0;
    float* d_xf1;
    float* d_xf2;
    float* d_xf0_sorted;
    float* d_xf1_sorted;
    float* d_xf2_sorted;
    float* d_af0;
    float* d_af1;
    float* d_af2;
    int64_t* d_ip_local;
    int* d_cell_start;
    int* d_cell_end;
    float* d_f2_out;
    float* d_vp;

    float* xf0_h;
    float* xf1_h;
    float* xf2_h;
    int64_t* ip_local_h;
    uint32_t* keys_h;
    uint32_t* values_h;

    float* f2max_teams;
    float* vmax_teams;
    float* vp_d2h_buffer[2];

    int max_rcp;
    int n_cap;
    int grid_num_cap;
    int upd_blocks_cap;

    cudaStream_t h2d_stream;
    cudaStream_t d2h_stream;
    cudaStream_t compute_stream;

    cudaEvent_t event_h2d_done;
    cudaEvent_t event_reduce_done;
};

static void init_pp_workspace(PPForceWorkspace* ws, int max_rcp) {
    memset(ws, 0, sizeof(*ws));

    ws->max_rcp = max_rcp;
    ws->n_cap = (int)np_tile_max;
    ws->upd_blocks_cap = (ws->n_cap + kUpdateThreads - 1) / kUpdateThreads;

    int stride = nt * max_rcp + 2 * max_rcp;
    ws->grid_num_cap = stride * stride * stride;

    CUDA_CHECK(cudaMallocHost(&ws->xf0_h, ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMallocHost(&ws->xf1_h, ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMallocHost(&ws->xf2_h, ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMallocHost(&ws->ip_local_h, ws->n_cap * sizeof(int64_t)));
    CUDA_CHECK(cudaMallocHost(&ws->keys_h, ws->n_cap * sizeof(uint32_t)));
    CUDA_CHECK(cudaMallocHost(&ws->values_h, ws->n_cap * sizeof(uint32_t)));

    CUDA_CHECK(cudaMallocHost(&ws->f2max_teams, (size_t)nnt * nnt * nnt * sizeof(float)));
    CUDA_CHECK(cudaMallocHost(&ws->vmax_teams, (size_t)nnt * nnt * nnt * 3 * sizeof(float)));
    CUDA_CHECK(cudaMallocHost(&ws->vp_d2h_buffer[0], VP_CHUNK_ELEMS * sizeof(float)));
    CUDA_CHECK(cudaMallocHost(&ws->vp_d2h_buffer[1], VP_CHUNK_ELEMS * sizeof(float)));

    CUDA_CHECK(cudaMalloc(&ws->d_keys, (size_t)ws->n_cap * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&ws->d_values, (size_t)ws->n_cap * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&ws->d_xf0, (size_t)ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_xf1, (size_t)ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_xf2, (size_t)ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_xf0_sorted, (size_t)ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_xf1_sorted, (size_t)ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_xf2_sorted, (size_t)ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_af0, (size_t)ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_af1, (size_t)ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_af2, (size_t)ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_ip_local, (size_t)ws->n_cap * sizeof(int64_t)));
    CUDA_CHECK(cudaMalloc(&ws->d_cell_start, (size_t)ws->grid_num_cap * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&ws->d_cell_end, (size_t)ws->grid_num_cap * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&ws->d_f2_out, sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_vp, max_vp * sizeof(float)));

    CUDA_CHECK(cudaStreamCreate(&ws->h2d_stream));
    CUDA_CHECK(cudaStreamCreate(&ws->d2h_stream));
    CUDA_CHECK(cudaStreamCreate(&ws->compute_stream));
    CUDA_CHECK(cudaEventCreate(&ws->event_h2d_done));
    CUDA_CHECK(cudaEventCreate(&ws->event_reduce_done));

    // Pre-record events to avoid first-time recording overhead during the main loop.
    cudaEventRecord(ws->event_reduce_done, ws->compute_stream);

    // warmup thrust sort
    auto exec = thrust::cuda::par(thrust::detail::single_device_tls_caching_allocator())
                    .on(ws->compute_stream);
    thrust::device_ptr<uint32_t> kptr(ws->d_keys);
    thrust::device_ptr<uint32_t> vptr(ws->d_values);
    thrust::sort_by_key(exec, kptr, kptr + ws->n_cap, vptr);
}

static void destroy_pp_workspace(PPForceWorkspace* ws) {
    CUDA_CHECK(cudaEventDestroy(ws->event_h2d_done));
    CUDA_CHECK(cudaEventDestroy(ws->event_reduce_done));
    CUDA_CHECK(cudaStreamDestroy(ws->h2d_stream));
    CUDA_CHECK(cudaStreamDestroy(ws->d2h_stream));
    CUDA_CHECK(cudaStreamDestroy(ws->compute_stream));

    CUDA_CHECK(cudaFreeHost(ws->xf0_h));
    CUDA_CHECK(cudaFreeHost(ws->xf1_h));
    CUDA_CHECK(cudaFreeHost(ws->xf2_h));
    CUDA_CHECK(cudaFreeHost(ws->ip_local_h));
    CUDA_CHECK(cudaFreeHost(ws->keys_h));
    CUDA_CHECK(cudaFreeHost(ws->values_h));

    CUDA_CHECK(cudaFreeHost(ws->f2max_teams));
    CUDA_CHECK(cudaFreeHost(ws->vmax_teams));
    CUDA_CHECK(cudaFreeHost(ws->vp_d2h_buffer[0]));
    CUDA_CHECK(cudaFreeHost(ws->vp_d2h_buffer[1]));

    CUDA_CHECK(cudaFree(ws->d_keys));
    CUDA_CHECK(cudaFree(ws->d_values));
    CUDA_CHECK(cudaFree(ws->d_xf0));
    CUDA_CHECK(cudaFree(ws->d_xf1));
    CUDA_CHECK(cudaFree(ws->d_xf2));
    CUDA_CHECK(cudaFree(ws->d_xf0_sorted));
    CUDA_CHECK(cudaFree(ws->d_xf1_sorted));
    CUDA_CHECK(cudaFree(ws->d_xf2_sorted));
    CUDA_CHECK(cudaFree(ws->d_af0));
    CUDA_CHECK(cudaFree(ws->d_af1));
    CUDA_CHECK(cudaFree(ws->d_af2));
    CUDA_CHECK(cudaFree(ws->d_ip_local));
    CUDA_CHECK(cudaFree(ws->d_cell_start));
    CUDA_CHECK(cudaFree(ws->d_cell_end));
    CUDA_CHECK(cudaFree(ws->d_f2_out));
    CUDA_CHECK(cudaFree(ws->d_vp));
}

static PPForceWorkspace g_pp_ws;
static int g_pp_ws_inited = 0;

static void pp_ws_atexit_cleanup() {
    if (g_pp_ws_inited) {
        destroy_pp_workspace(&g_pp_ws);
        g_pp_ws_inited = 0;
        // Optional: zero-initialize the struct for safety.
        memset(&g_pp_ws, 0, sizeof(g_pp_ws));
    }
}

static inline PPForceWorkspace* get_pp_workspace_static(int max_rcp) {
    if (!g_pp_ws_inited) {
        init_pp_workspace(&g_pp_ws, max_rcp);
        g_pp_ws_inited = 1;
        atexit(pp_ws_atexit_cleanup);
    } else {
        if (g_pp_ws.max_rcp != max_rcp) {
            fprintf(stderr, "[PP] max_rcp changed: old=%d new=%d (static ws, no resize)\n",
                    g_pp_ws.max_rcp, max_rcp);
            exit(EXIT_FAILURE);
        }
    }
    return &g_pp_ws;
}

static inline void accumulate_vp_from_device(PPForceWorkspace* ws, float* vp_host_out, size_t vp_elems_used, float vmax[3]) {

    vmax[0] = 0.0f;
    vmax[1] = 0.0f;
    vmax[2] = 0.0f;

    const size_t vp_elems = vp_elems_used;
    const size_t chunk_elems = VP_CHUNK_ELEMS;
    const size_t num_chunks = (vp_elems + chunk_elems - 1) / chunk_elems;
    if (num_chunks == 0) return;

    cudaStream_t copy_streams[kD2HBufferCount] = {ws->d2h_stream, ws->h2d_stream};
    size_t offset_prev = 0;
    size_t elems_prev = 0;

    // Prime pipeline: launch chunk 0 copy.
    const size_t elems0 = (vp_elems < chunk_elems) ? vp_elems : chunk_elems;
    CUDA_CHECK(cudaMemcpyAsync(ws->vp_d2h_buffer[0], ws->d_vp, elems0 * sizeof(float),
                               cudaMemcpyDeviceToHost, copy_streams[0]));
    offset_prev = 0;
    elems_prev = elems0;

    auto add_vp_chunk = [&vmax](float* __restrict__ dst, const float* __restrict__ src, size_t count) {
        const int n = (int)count;
        float _vmax_0 = 0.0f;
        float _vmax_1 = 0.0f;
        float _vmax_2 = 0.0f;
#pragma omp parallel for simd num_threads(ncore) schedule(static) reduction(max : _vmax_0, _vmax_1, _vmax_2)
        for (int i = 0; i < n; i++) {
            dst[i] += src[i];
            if (src[i] != 0.0f) {
                float v = fabsf(dst[i]);
                if (i % 3 == 0) {
                    _vmax_0 = fmaxf(_vmax_0, v);
                } else if (i % 3 == 1) {
                    _vmax_1 = fmaxf(_vmax_1, v);
                } else {
                    _vmax_2 = fmaxf(_vmax_2, v);
                }
            }
        }
        vmax[0] = fmaxf(vmax[0], _vmax_0);
        vmax[1] = fmaxf(vmax[1], _vmax_1);
        vmax[2] = fmaxf(vmax[2], _vmax_2);
    };

    for (size_t c = 1; c < num_chunks; c++) {
        const int cur = (int)(c & 1);
        const int prev = cur ^ 1;

        const size_t offset_cur = c * chunk_elems;
        const size_t remain = vp_elems - offset_cur;
        const size_t elems_cur = (remain < chunk_elems) ? remain : chunk_elems;

        // Launch current copy first, so GPU copy can overlap with host accumulation of previous
        // chunk.
        CUDA_CHECK(cudaMemcpyAsync(ws->vp_d2h_buffer[cur], ws->d_vp + offset_cur,
                                   elems_cur * sizeof(float), cudaMemcpyDeviceToHost,
                                   copy_streams[cur]));

        CUDA_CHECK(cudaStreamSynchronize(copy_streams[prev]));
        add_vp_chunk(vp_host_out + offset_prev, ws->vp_d2h_buffer[prev], elems_prev);

        offset_prev = offset_cur;
        elems_prev = elems_cur;
    }

    CUDA_CHECK(cudaStreamSynchronize(copy_streams[(num_chunks - 1) & 1]));
    add_vp_chunk(vp_host_out + offset_prev, ws->vp_d2h_buffer[(num_chunks - 1) & 1], elems_prev);
}

extern "C" void c_pp_force_kernel_tile(int itile, int iapm, int tile1[3], int nc1[3], int nc2[3],
                                       int utile_shift[3], int ratio_sf[], float apm3[], int* _rhoc,
                                       int64_t* _idx_b_r, int16_t* _xp, float* d_vp_flat,
                                       float mass_p_cdm, float a_mid, float dt,
                                       PPForceWorkspace* ws, float* f2max_t_,
                                       double* ptotal, double* ttotal) {
    struct timespec t_total0, t_total1;
    clock_gettime(CLOCK_MONOTONIC, &t_total0);

    float pp_range = apm3[iapm - 1];
    int rcp = ratio_sf[iapm - 1];
    int npgrid = nt * rcp;

    int (*rhoc)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb][nt + 2 * ncb] =
        (int (*)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb][nt + 2 * ncb]) _rhoc;
    int64_t (*idx_b_r)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb] =
        (int64_t (*)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb]) _idx_b_r;

    int16_t (*xp)[3] = (int16_t (*)[3])_xp;

    int nc1_offset[3];
    int nc2_offset[3];
    for (int d = 0; d < 3; d++) {
        nc1_offset[d] = nc1[d] + ncb - 1;
        nc2_offset[d] = nc2[d] + ncb - 1;
    }

    int tile1_offset[3];
    for (int d = 0; d < 3; d++) {
        tile1_offset[d] = tile1[d] - 1;
    }

    int nptile = 0;
    for (int k = nc1_offset[2] - 1; k <= nc2_offset[2] + 1; k++) {
        for (int j = nc1_offset[1] - 1; j <= nc2_offset[1] + 1; j++) {
            for (int i = nc1_offset[0] - 1; i <= nc2_offset[0] + 1; i++) {
                nptile += rhoc[tile1_offset[2]][tile1_offset[1]][tile1_offset[0]][k][j][i];
            }
        }
    }

    int nptile_real = 0;
    for (int k = nc1_offset[2]; k <= nc2_offset[2]; k++) {
        for (int j = nc1_offset[1]; j <= nc2_offset[1]; j++) {
            for (int i = nc1_offset[0]; i <= nc2_offset[0]; i++) {
                nptile_real += rhoc[tile1_offset[2]][tile1_offset[1]][tile1_offset[0]][k][j][i];
            }
        }
    }

    int nk = (nc2_offset[2] + 1) - (nc1_offset[2] - 1) + 1;
    int nj = (nc2_offset[1] + 1) - (nc1_offset[1] - 1) + 1;
    int ni = (nc2_offset[0] + 1) - (nc1_offset[0] - 1) + 1;
    int ncell = nk * nj * ni;

    auto cell_id = [&](int k, int j, int i) {
        int kk = k - (nc1_offset[2] - 1);
        int jj = j - (nc1_offset[1] - 1);
        int ii = i - (nc1_offset[0] - 1);
        return (kk * nj + jj) * ni + ii;
    };

    std::vector<int> count(ncell);
    std::vector<int64_t> offset(ncell);

#pragma omp parallel for collapse(3) num_threads(ncore)
    for (int k = nc1_offset[2] - 1; k <= nc2_offset[2] + 1; k++) {
        for (int j = nc1_offset[1] - 1; j <= nc2_offset[1] + 1; j++) {
            for (int i = nc1_offset[0] - 1; i <= nc2_offset[0] + 1; i++) {
                int np = rhoc[tile1_offset[2]][tile1_offset[1]][tile1_offset[0]][k][j][i];
                count[cell_id(k, j, i)] = np;
            }
        }
    }

    offset[0] = 0;
    for (int c = 1; c < ncell; c++) offset[c] = offset[c - 1] + count[c - 1];
    int64_t total = offset[ncell - 1] + count[ncell - 1];

#pragma omp parallel for collapse(3) num_threads(ncore)
    for (int k = nc1_offset[2] - 1; k <= nc2_offset[2] + 1; k++) {
        for (int j = nc1_offset[1] - 1; j <= nc2_offset[1] + 1; j++) {
            for (int i = nc1_offset[0] - 1; i <= nc2_offset[0] + 1; i++) {
                int cid = cell_id(k, j, i);
                int np = count[cid];
                int64_t base = offset[cid];

                int64_t rhoc_i_sum = 0;
                for (int index = i; index < nt + 2 * ncb; index++) {
                    rhoc_i_sum +=
                        rhoc[tile1_offset[2]][tile1_offset[1]][tile1_offset[0]][k][j][index];
                }
                int64_t nzero =
                    idx_b_r[tile1_offset[2]][tile1_offset[1]][tile1_offset[0]][k][j] - rhoc_i_sum;

                for (int l = 0; l < np; l++) {
                    int64_t ip1 = base + l;
                    int64_t ip = nzero + l;

                    ws->ip_local_h[ip1] = ip;

                    float x0 = (float)((i - ncb) +
                                       ((int16_t)(xp[ip][0] + ishift) + rshift) * x_resolution);
                    float x1 = (float)((j - ncb) +
                                       ((int16_t)(xp[ip][1] + ishift) + rshift) * x_resolution);
                    float x2 = (float)((k - ncb) +
                                       ((int16_t)(xp[ip][2] + ishift) + rshift) * x_resolution);

                    ws->xf0_h[ip1] = ratio_cs * x0;
                    ws->xf1_h[ip1] = ratio_cs * x1;
                    ws->xf2_h[ip1] = ratio_cs * x2;

                    int idx_0 = (int)floorf((float)rcp * (x0 - (utile_shift[0] * ntt))) + rcp;
                    int idx_1 = (int)floorf((float)rcp * (x1 - (utile_shift[1] * ntt))) + rcp;
                    int idx_2 = (int)floorf((float)rcp * (x2 - (utile_shift[2] * ntt))) + rcp;

                    uint32_t idx_ = (uint32_t)idx_2 * (uint32_t)(npgrid + 2 * rcp) *
                                        (uint32_t)(npgrid + 2 * rcp) +
                                    (uint32_t)idx_1 * (uint32_t)(npgrid + 2 * rcp) +
                                    (uint32_t)idx_0;

                    ws->keys_h[ip1] = idx_;
                    ws->values_h[ip1] = (uint32_t)ip1;
                }
            }
        }
    }

    int N = (int)total;
    if (N > ws->n_cap) {
        fprintf(stderr, "N=%d exceeds workspace n_cap=%d (np_tile_max)\n", N, ws->n_cap);
        exit(EXIT_FAILURE);
    }

    CUDA_CHECK(cudaMemcpyAsync(ws->d_keys, ws->keys_h, (size_t)N * sizeof(uint32_t),
                               cudaMemcpyHostToDevice, ws->h2d_stream));
    CUDA_CHECK(cudaMemcpyAsync(ws->d_xf0, ws->xf0_h, (size_t)N * sizeof(float),
                               cudaMemcpyHostToDevice, ws->h2d_stream));
    CUDA_CHECK(cudaMemcpyAsync(ws->d_xf1, ws->xf1_h, (size_t)N * sizeof(float),
                               cudaMemcpyHostToDevice, ws->h2d_stream));
    CUDA_CHECK(cudaMemcpyAsync(ws->d_xf2, ws->xf2_h, (size_t)N * sizeof(float),
                               cudaMemcpyHostToDevice, ws->h2d_stream));

    cudaStreamWaitEvent(ws->h2d_stream, ws->event_reduce_done, 0);

    CUDA_CHECK(cudaMemcpyAsync(ws->d_values, ws->values_h, (size_t)N * sizeof(uint32_t),
                               cudaMemcpyHostToDevice, ws->h2d_stream));
    CUDA_CHECK(cudaMemcpyAsync(ws->d_ip_local, ws->ip_local_h, (size_t)N * sizeof(int64_t),
                               cudaMemcpyHostToDevice, ws->h2d_stream));

    CUDA_CHECK(cudaMemsetAsync(ws->d_af0, 0, (size_t)N * sizeof(float), ws->h2d_stream));
    CUDA_CHECK(cudaMemsetAsync(ws->d_af1, 0, (size_t)N * sizeof(float), ws->h2d_stream));
    CUDA_CHECK(cudaMemsetAsync(ws->d_af2, 0, (size_t)N * sizeof(float), ws->h2d_stream));

    cudaEventRecord(ws->event_h2d_done, ws->h2d_stream);
    cudaStreamWaitEvent(ws->compute_stream, ws->event_h2d_done, 0);

    auto exec = thrust::cuda::par(thrust::detail::single_device_tls_caching_allocator())
                    .on(ws->compute_stream);

    // -------------------------
    // thrust sort
    // -------------------------
    {
        thrust::device_ptr<uint32_t> kptr(ws->d_keys);
        thrust::device_ptr<uint32_t> vptr(ws->d_values);
        thrust::sort_by_key(exec, kptr, kptr + N, vptr);
        thrust::device_ptr<uint32_t> map(ws->d_values);
        thrust::device_ptr<float> in0(ws->d_xf0), in1(ws->d_xf1), in2(ws->d_xf2);
        thrust::device_ptr<float> out0(ws->d_xf0_sorted), out1(ws->d_xf1_sorted),
            out2(ws->d_xf2_sorted);
        thrust::gather(exec, map, map + N, in0, out0);
        thrust::gather(exec, map, map + N, in1, out1);
        thrust::gather(exec, map, map + N, in2, out2);
    }

    // -------------------------
    // build cell ranges
    // -------------------------
    {
        int grid_num = (npgrid + 2 * rcp) * (npgrid + 2 * rcp) * (npgrid + 2 * rcp);
        if (grid_num > ws->grid_num_cap) {
            fprintf(stderr, "grid_num=%d exceeds workspace grid_num_cap=%d\n", grid_num,
                    ws->grid_num_cap);
            exit(EXIT_FAILURE);
        }
        CUDA_CHECK(cudaMemsetAsync(ws->d_cell_start, 0xFF, (size_t)grid_num * sizeof(int),
                                ws->compute_stream));
        CUDA_CHECK(
            cudaMemsetAsync(ws->d_cell_end, 0xFF, (size_t)grid_num * sizeof(int), ws->compute_stream));
        const int blocks = (N + kCellBuildThreads - 1) / kCellBuildThreads;
        build_cell_start_end_kernel<<<blocks, kCellBuildThreads, 0, ws->compute_stream>>>(
            ws->d_keys, N, ws->d_cell_start, ws->d_cell_end);
        CUDA_CHECK(cudaGetLastError());

        thrust::device_ptr<int> start_ptr(ws->d_cell_start);
        thrust::device_ptr<int> end_ptr(ws->d_cell_end);

        thrust::replace(exec, start_ptr, start_ptr + grid_num, -1, N);
        thrust::replace(exec, end_ptr, end_ptr + grid_num, -1, N);

        auto rstart_begin = thrust::make_reverse_iterator(start_ptr + grid_num);
        auto rstart_end = thrust::make_reverse_iterator(start_ptr);
        auto rend_begin = thrust::make_reverse_iterator(end_ptr + grid_num);
        auto rend_end = thrust::make_reverse_iterator(end_ptr);

        thrust::inclusive_scan(exec, rstart_begin, rstart_end, rstart_begin,
                               thrust::minimum<int>());
        thrust::inclusive_scan(exec, rend_begin, rend_end, rend_begin, thrust::minimum<int>());
    }

    // -------------------------
    // pp force kernel
    // -------------------------
    int real_cells = npgrid * npgrid * npgrid;
    {
        pp_force_kernel<<<real_cells, kPPThreads, 0, ws->compute_stream>>>(
            ws->d_xf0_sorted, ws->d_xf1_sorted, ws->d_xf2_sorted, npgrid, rcp, pp_range, (float)app,
            ws->d_cell_start, ws->d_cell_end, ws->d_af0, ws->d_af1, ws->d_af2);
        CUDA_CHECK(cudaGetLastError());
    }

    // -------------------------
    // update vp + reduce max
    // -------------------------
    {
        int upd_blocks = (N + kUpdateThreads - 1) / kUpdateThreads;
        if (upd_blocks > ws->upd_blocks_cap) {
            fprintf(stderr, "upd_blocks=%d exceeds workspace upd_blocks_cap=%d\n", upd_blocks,
                    ws->upd_blocks_cap);
            exit(EXIT_FAILURE);
        }
        CUDA_CHECK(cudaMemsetAsync(ws->d_f2_out, 0, sizeof(float), ws->compute_stream));
        size_t shmem = (size_t)kUpdateThreads * sizeof(float);
        update_v_kernel<<<upd_blocks, kUpdateThreads, shmem, ws->compute_stream>>>(
            ws->d_af0, ws->d_af1, ws->d_af2, ws->d_ip_local, ws->d_values, d_vp_flat, N,
            mass_p_cdm * mass_p_cdm, a_mid, dt, ws->d_f2_out);
        CUDA_CHECK(cudaGetLastError());
        cudaEventRecord(ws->event_reduce_done, ws->compute_stream);
    }

    cudaStreamWaitEvent(ws->d2h_stream, ws->event_reduce_done, 0);

    CUDA_CHECK(cudaMemcpyAsync(f2max_t_, ws->d_f2_out, sizeof(float), cudaMemcpyDeviceToHost,
                               ws->d2h_stream));

    ptotal[itile] = (double)nptile_real;
    clock_gettime(CLOCK_MONOTONIC, &t_total1);
    ttotal[itile] = elapsed_sec(t_total0, t_total1);
}

extern "C" void c_pp_force_kernel(int isort[ns3], int ires[ns3], int ixyz3[ns3][6], float apm3[7],
                                  int ratio_sf[7], int* _rhoc, int64_t* _idx_b_r, int16_t* _xp,
                                  float* _vp, float mass_p_cdm, float a_mid, float dt, float* f2max,
                                  float vmax[3], double* ptotal, double* ttotal) {
    int max_rcp = ratio_sf[0];
    for (int i = 1; i < 7; i++) {
        if (ratio_sf[i] > max_rcp) max_rcp = ratio_sf[i];
    }
    PPForceWorkspace* ws = get_pp_workspace_static(max_rcp);

    CUDA_CHECK(cudaMemset(ws->d_vp, 0, max_vp * sizeof(float)));

    for (int it = 0; it < nnt * nnt * nnt; it++) {
        int itile = it;
        int iapm = ires[itile];

        int tile1[3] = {it / (nnt * nnt) + 1, (it / nnt) % nnt + 1, it % nnt + 1};
        int utile_shift[3] = {0, 0, 0};
        int nc1[3] = {1, 1, 1};
        int nc2[3] = {nt, nt, nt};

        c_pp_force_kernel_tile(itile, iapm, tile1, nc1, nc2, utile_shift, ratio_sf, apm3, _rhoc,
                               _idx_b_r, _xp, ws->d_vp, mass_p_cdm, a_mid, dt, ws,
                               &(ws->f2max_teams[it]), ptotal, ttotal);
    }

    float _f2max = 0.0f;

    CUDA_CHECK(cudaDeviceSynchronize());

    for (int it = 0; it < nnt * nnt * nnt; it++) {
        _f2max = fmaxf(_f2max, ws->f2max_teams[it]);
    }

    *f2max = _f2max;

    int64_t (*idx_b_r)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb] =
        (int64_t (*)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb]) _idx_b_r;
    const int j_last = nt + 2 * ncb - 1;
    const int k_last = nt + 2 * ncb - 1;
    int64_t n_used = idx_b_r[nnt - 1][nnt - 1][nnt - 1][k_last][j_last]; // particle count
    size_t vp_elems_used = (size_t)n_used * 3;
    if (vp_elems_used > max_vp) vp_elems_used = max_vp;

    accumulate_vp_from_device(ws, _vp, vp_elems_used, vmax);
}
