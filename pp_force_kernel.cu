#include <cuda_runtime.h>
#include "cuda_parameters_generated.h"
#include <limits.h>
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

#include <algorithm>
#include <vector>

// -------------------------
// Global configuration
// -------------------------
const int ng = CUBE_NG;
const int nnt = CUBE_NNT;
const int nns = CUBE_NNS;
const int ratio_cs = CUBE_RATIO_CS;
const int nc = ng / ratio_cs;
const int np_nc = CUBE_NP_NC;
const int nt = nc / nnt;
const int ntt = nt / nns;
const int ngb = CUBE_NGB;
const int ncb = ngb / ratio_cs;
const int nte = nt + 2 * ncb;
const int ns3 = (nnt * nns) * (nnt * nns) * (nnt * nns);

const bool body_centered_cubic = CUBE_BODY_CENTERED_CUBIC;
const double image_buffer = CUBE_IMAGE_BUFFER;
const double tile_buffer = CUBE_TILE_BUFFER;

const int64_t np_image =
    (int64_t)(nc * np_nc) * (nc * np_nc) * (nc * np_nc) *
    (1 + body_centered_cubic);
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

const double app = CUBE_APP;

const int ncore = CUBE_NCORE;

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
                                const float* __restrict__ xf2, const float* __restrict__ pp_range,
                                const uint8_t* __restrict__ active, const uint32_t* __restrict__ keys,
                                int N, int npgrid, int rcp,
                                float app_local, const int* __restrict__ cell_start,
                                const int* __restrict__ cell_end, float* __restrict__ af0,
                                float* __restrict__ af1, float* __restrict__ af2) {
    int stride = (npgrid + 2 * rcp);
    int outer = blockIdx.x * blockDim.x + threadIdx.x;
    if (outer < N) {
        if (!active[outer]) return;
        const int key = (int)keys[outer];
        const int k = key / (stride * stride) - rcp;
        const int j = (key / stride) % stride - rcp;
        const int i = key % stride - rcp;
        float xo = xf0[outer], yo = xf1[outer], zo = xf2[outer];
        float ax = 0.f, ay = 0.f, az = 0.f;
        float outer_range = pp_range[outer];
        float pp2 = outer_range * outer_range;
        const int outer_cell_radius = (int)ceilf(outer_range * (float)rcp / (float)ratio_cs);

        for (int kk = -outer_cell_radius; kk <= outer_cell_radius; kk++) {
            // #pragma unroll
            for (int jj = -outer_cell_radius; jj <= outer_cell_radius; jj++) {
                // #pragma unroll
                for (int ii = -outer_cell_radius; ii <= outer_cell_radius; ii++) {
                    if (k + kk < -rcp || k + kk >= npgrid + rcp || j + jj < -rcp ||
                        j + jj >= npgrid + rcp || i + ii < -rcp || i + ii >= npgrid + rcp)
                        continue;
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
                            float fpp = F_ra_dev(r, app_local) - F_ra_dev(r, outer_range);
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
                                const float* __restrict__ af2, const uint32_t* __restrict__ ip_local,
                                const uint32_t* __restrict__ values_sorted,
                                float* __restrict__ vp_flat, int N, float mass2, float a_mid,
                                float dt, const uint8_t* __restrict__ active,
                                float* __restrict__ f2_global) {
    extern __shared__ float s[];
    float* sf2 = s;

    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + tid;

    float f2 = 0.f;

    if (idx < N && active[idx]) {
        float ax = af0[idx] * mass2;
        float ay = af1[idx] * mass2;
        float az = af2[idx] * mass2;

        f2 = ax * ax + ay * ay + az * az;

        uint32_t ori = values_sorted[idx];
        uint32_t ip = ip_local[ori];

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
    float* d_pp_range;
    float* d_pp_range_sorted;
    uint32_t* d_ip_local;
    uint8_t* d_active;
    uint8_t* d_active_sorted;
    int* d_cell_start;
    int* d_cell_end;
    float* d_f2_out;
    float* d_vp;

    float* xf0_h;
    float* xf1_h;
    float* xf2_h;
    uint32_t* ip_local_h;
    uint32_t* keys_h;
    uint32_t* values_h;
    float* pp_range_h;
    uint8_t* active_h;

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
    CUDA_CHECK(cudaMallocHost(&ws->ip_local_h, ws->n_cap * sizeof(uint32_t)));
    CUDA_CHECK(cudaMallocHost(&ws->keys_h, ws->n_cap * sizeof(uint32_t)));
    CUDA_CHECK(cudaMallocHost(&ws->values_h, ws->n_cap * sizeof(uint32_t)));
    CUDA_CHECK(cudaMallocHost(&ws->pp_range_h, ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMallocHost(&ws->active_h, ws->n_cap * sizeof(uint8_t)));

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
    CUDA_CHECK(cudaMalloc(&ws->d_pp_range, (size_t)ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_pp_range_sorted, (size_t)ws->n_cap * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_ip_local, (size_t)ws->n_cap * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&ws->d_active, (size_t)ws->n_cap * sizeof(uint8_t)));
    CUDA_CHECK(cudaMalloc(&ws->d_active_sorted, (size_t)ws->n_cap * sizeof(uint8_t)));
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
    CUDA_CHECK(cudaFreeHost(ws->pp_range_h));
    CUDA_CHECK(cudaFreeHost(ws->active_h));

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
    CUDA_CHECK(cudaFree(ws->d_pp_range));
    CUDA_CHECK(cudaFree(ws->d_pp_range_sorted));
    CUDA_CHECK(cudaFree(ws->d_ip_local));
    CUDA_CHECK(cudaFree(ws->d_active));
    CUDA_CHECK(cudaFree(ws->d_active_sorted));
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

static inline void accumulate_vp_from_device(PPForceWorkspace* ws, float* vp_host_out,
                                             const uint8_t* active_particles, size_t vp_elems_used,
                                             float vmax[3]) {

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

    auto add_vp_chunk = [&vmax, active_particles, vp_host_out](float* __restrict__ dst,
                                                                 const float* __restrict__ src,
                                                                 size_t count) {
        const int n = (int)count;
        float _vmax_0 = 0.0f;
        float _vmax_1 = 0.0f;
        float _vmax_2 = 0.0f;
#pragma omp parallel for simd num_threads(ncore) schedule(static) reduction(max : _vmax_0, _vmax_1, _vmax_2)
        for (int i = 0; i < n; i++) {
            dst[i] += src[i];
            const size_t global_component = (size_t)(dst - vp_host_out) + (size_t)i;
            if (active_particles[global_component / 3]) {
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

extern "C" void c_pp_force_kernel_tile(const int ires[ns3], const int subtile_owner[ns3],
                                       int tile1[3], int nc1[3], int nc2[3],
                                       int utile_shift[3], int ratio_sf[], float apm3[], int* _rhoc,
                                       int64_t* _idx_b_r, int16_t* _xp, float* d_vp_flat,
                                       float mass_p_cdm, float a_mid, float dt,
                                       PPForceWorkspace* ws, uint8_t* active_particles,
                                       size_t active_particle_count, float* f2max_t_,
                                       double* ptotal, double* ttotal) {
    // The pinned staging buffers are shared by all tiles.  Do not overwrite
    // them while the previous tile's asynchronous H2D copies may still read
    // from them.
    CUDA_CHECK(cudaStreamSynchronize(ws->h2d_stream));

    struct timespec t_total0, t_total1;
    clock_gettime(CLOCK_MONOTONIC, &t_total0);

    // A tile contains nns^3 PP subtiles.  Use one fine-grid launch for the tile;
    // each physical particle retains the PP range selected by its own ires entry.
    int rcp = ratio_sf[1];
    for (int sz = 0; sz < nns; ++sz) {
        for (int sy = 0; sy < nns; ++sy) {
            for (int sx = 0; sx < nns; ++sx) {
                const int tx = tile1[0] - 1;
                const int ty = tile1[1] - 1;
                const int tz = tile1[2] - 1;
                const int subtile = (((tz * nnt + ty) * nnt + tx) * nns + sz) * nns * nns +
                                    sy * nns + sx;
                const int iapm = ires[subtile];
                if (iapm < 2 || iapm > 6) {
                    fprintf(stderr, "invalid PP ires=%d for subtile %d\n", iapm, subtile);
                    exit(EXIT_FAILURE);
                }
                rcp = std::max(rcp, ratio_sf[iapm - 1]);
            }
        }
    }
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
                const int sx = (i - nc1_offset[0]) / ntt;
                const int sy = (j - nc1_offset[1]) / ntt;
                const int sz = (k - nc1_offset[2]) / ntt;
                const int subtile = (((tile1_offset[2] * nnt + tile1_offset[1]) * nnt +
                                      tile1_offset[0]) * nns + sz) * nns * nns + sy * nns + sx;
                const int np = rhoc[tile1_offset[2]][tile1_offset[1]][tile1_offset[0]][k][j][i];
                nptile_real += np;
                ptotal[subtile_owner[subtile]] += (double)np;
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
    std::vector<int> offset(ncell);

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
    int total = offset[ncell - 1] + count[ncell - 1];

#pragma omp parallel for collapse(3) num_threads(ncore)
    for (int k = nc1_offset[2] - 1; k <= nc2_offset[2] + 1; k++) {
        for (int j = nc1_offset[1] - 1; j <= nc2_offset[1] + 1; j++) {
            for (int i = nc1_offset[0] - 1; i <= nc2_offset[0] + 1; i++) {
                int cid = cell_id(k, j, i);
                int np = count[cid];
                int base = offset[cid];

                int rhoc_i_sum = 0;
                for (int index = i; index < nt + 2 * ncb; index++) {
                    rhoc_i_sum +=
                        rhoc[tile1_offset[2]][tile1_offset[1]][tile1_offset[0]][k][j][index];
                }
                const int64_t nzero64 =
                    idx_b_r[tile1_offset[2]][tile1_offset[1]][tile1_offset[0]][k][j] - rhoc_i_sum;
                if (nzero64 < 0 || nzero64 > INT_MAX) {
                    fprintf(stderr, "PP particle base index %lld exceeds int32 range\n",
                            (long long)nzero64);
                    exit(EXIT_FAILURE);
                }
                const int nzero = (int)nzero64;

                for (int l = 0; l < np; l++) {
                    const int ip1 = base + l;
                    const int ip = nzero + l;

                    ws->ip_local_h[ip1] = (uint32_t)ip;

                    float x0 = (float)((i - ncb) +
                                       ((int16_t)(xp[ip][0] + ishift) + rshift) * x_resolution);
                    float x1 = (float)((j - ncb) +
                                       ((int16_t)(xp[ip][1] + ishift) + rshift) * x_resolution);
                    float x2 = (float)((k - ncb) +
                                       ((int16_t)(xp[ip][2] + ishift) + rshift) * x_resolution);

                    ws->xf0_h[ip1] = ratio_cs * x0;
                    ws->xf1_h[ip1] = ratio_cs * x1;
                    ws->xf2_h[ip1] = ratio_cs * x2;

                    const bool is_real = i >= nc1_offset[0] && i <= nc2_offset[0] &&
                                         j >= nc1_offset[1] && j <= nc2_offset[1] &&
                                         k >= nc1_offset[2] && k <= nc2_offset[2];
                    ws->active_h[ip1] = is_real;
                    if (is_real) {
                        if ((size_t)ip >= active_particle_count) {
                            fprintf(stderr, "PP particle index %lld exceeds active-particle capacity %zu\n",
                                    (long long)ip, active_particle_count);
                            exit(EXIT_FAILURE);
                        }
                        active_particles[ip] = 1;
                        const int sx = (i - nc1_offset[0]) / ntt;
                        const int sy = (j - nc1_offset[1]) / ntt;
                        const int sz = (k - nc1_offset[2]) / ntt;
                        const int tx = tile1[0] - 1;
                        const int ty = tile1[1] - 1;
                        const int tz = tile1[2] - 1;
                        const int subtile = (((tz * nnt + ty) * nnt + tx) * nns + sz) * nns * nns +
                                            sy * nns + sx;
                        const int iapm = ires[subtile];
                        ws->pp_range_h[ip1] = apm3[iapm - 1];
                    } else {
                        // Halo particles are interaction sources only.
                        ws->pp_range_h[ip1] = 0.0f;
                    }

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

    int N = total;
    if (N > ws->n_cap) {
        fprintf(stderr, "N=%d exceeds workspace n_cap=%d (np_tile_max)\n", N, ws->n_cap);
        exit(EXIT_FAILURE);
    }

    // Reuse of the single workspace is safe only after the preceding tile's
    // reduction (and therefore all reads of these buffers) has completed.
    cudaStreamWaitEvent(ws->h2d_stream, ws->event_reduce_done, 0);

    CUDA_CHECK(cudaMemcpyAsync(ws->d_keys, ws->keys_h, (size_t)N * sizeof(uint32_t),
                               cudaMemcpyHostToDevice, ws->h2d_stream));
    CUDA_CHECK(cudaMemcpyAsync(ws->d_xf0, ws->xf0_h, (size_t)N * sizeof(float),
                               cudaMemcpyHostToDevice, ws->h2d_stream));
    CUDA_CHECK(cudaMemcpyAsync(ws->d_xf1, ws->xf1_h, (size_t)N * sizeof(float),
                               cudaMemcpyHostToDevice, ws->h2d_stream));
    CUDA_CHECK(cudaMemcpyAsync(ws->d_xf2, ws->xf2_h, (size_t)N * sizeof(float),
                               cudaMemcpyHostToDevice, ws->h2d_stream));
    CUDA_CHECK(cudaMemcpyAsync(ws->d_pp_range, ws->pp_range_h, (size_t)N * sizeof(float),
                               cudaMemcpyHostToDevice, ws->h2d_stream));
    CUDA_CHECK(cudaMemcpyAsync(ws->d_active, ws->active_h, (size_t)N * sizeof(uint8_t),
                               cudaMemcpyHostToDevice, ws->h2d_stream));

    CUDA_CHECK(cudaMemcpyAsync(ws->d_values, ws->values_h, (size_t)N * sizeof(uint32_t),
                               cudaMemcpyHostToDevice, ws->h2d_stream));
    CUDA_CHECK(cudaMemcpyAsync(ws->d_ip_local, ws->ip_local_h, (size_t)N * sizeof(uint32_t),
                               cudaMemcpyHostToDevice, ws->h2d_stream));

    // CUDA_CHECK(cudaMemsetAsync(ws->d_af0, 0, (size_t)N * sizeof(float), ws->h2d_stream));
    // CUDA_CHECK(cudaMemsetAsync(ws->d_af1, 0, (size_t)N * sizeof(float), ws->h2d_stream));
    // CUDA_CHECK(cudaMemsetAsync(ws->d_af2, 0, (size_t)N * sizeof(float), ws->h2d_stream));

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
        thrust::device_ptr<float> in_range(ws->d_pp_range), out_range(ws->d_pp_range_sorted);
        thrust::device_ptr<uint8_t> in_active(ws->d_active), out_active(ws->d_active_sorted);
        thrust::gather(exec, map, map + N, in0, out0);
        thrust::gather(exec, map, map + N, in1, out1);
        thrust::gather(exec, map, map + N, in2, out2);
        thrust::gather(exec, map, map + N, in_range, out_range);
        thrust::gather(exec, map, map + N, in_active, out_active);
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
        // `0` is the insertion point before the first populated cell.  A
        // forward maximum scan then propagates each populated cell's end to
        // the empty cells immediately following it.
        thrust::replace(exec, end_ptr, end_ptr + grid_num, -1, 0);

        auto rstart_begin = thrust::make_reverse_iterator(start_ptr + grid_num);
        auto rstart_end = thrust::make_reverse_iterator(start_ptr);
        thrust::inclusive_scan(exec, rstart_begin, rstart_end, rstart_begin,
                               thrust::minimum<int>());
        // For an empty cell, end must be the previous populated cell's end
        // (and therefore equal its lower bound), not the next cell's end.
        thrust::inclusive_scan(exec, end_ptr, end_ptr + grid_num, end_ptr,
                               thrust::maximum<int>());
    }

    // -------------------------
    // pp force kernel
    // -------------------------
    int pp_blocks = (N + kPPThreads - 1) / kPPThreads;
    {
        // Each particle uses the cell radius implied by its own ires.  This
        // preserves PP_Force's per-subtile cutoff without wasting work on the
        // largest-radius stencil for high-resolution particles.
        pp_force_kernel<<<pp_blocks, kPPThreads, 0, ws->compute_stream>>>(
            ws->d_xf0_sorted, ws->d_xf1_sorted, ws->d_xf2_sorted, ws->d_pp_range_sorted,
            ws->d_active_sorted, ws->d_keys, N, npgrid, rcp, (float)app, ws->d_cell_start,
            ws->d_cell_end, ws->d_af0, ws->d_af1, ws->d_af2);
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
            mass_p_cdm * mass_p_cdm, a_mid, dt, ws->d_active_sorted, ws->d_f2_out);
        CUDA_CHECK(cudaGetLastError());
        cudaEventRecord(ws->event_reduce_done, ws->compute_stream);
    }

    cudaStreamWaitEvent(ws->d2h_stream, ws->event_reduce_done, 0);

    CUDA_CHECK(cudaMemcpyAsync(f2max_t_, ws->d_f2_out, sizeof(float), cudaMemcpyDeviceToHost,
                               ws->d2h_stream));

    clock_gettime(CLOCK_MONOTONIC, &t_total1);
    // Timing is aggregated per parent tile; particle counts above retain the
    // original PP-subtile indexing used by kick.f90.
    const int tile_first_subtile = (((tile1_offset[2] * nnt + tile1_offset[1]) * nnt +
                                     tile1_offset[0]) * nns) * nns * nns;
    ttotal[subtile_owner[tile_first_subtile]] = elapsed_sec(t_total0, t_total1);
}

extern "C" void c_pp_force_kernel(int isort[ns3], int ires[ns3], int ixyz3[ns3][6], float apm3[7],
                                  int ratio_sf[7], int* _rhoc, int64_t* _idx_b_r, int16_t* _xp,
                                  float* _vp, float mass_p_cdm, float a_mid, float dt, float* f2max,
                                  float vmax[3], double* ptotal, double* ttotal) {
    if (ng > 512 || np_image_max > INT_MAX || np_tile_max > INT_MAX) {
        fprintf(stderr,
                "PP CUDA kernel requires ng <= 512 and int32-safe particle capacities "
                "(ng=%d, image_max=%lld, tile_max=%lld)\n",
                ng, (long long)np_image_max, (long long)np_tile_max);
        exit(EXIT_FAILURE);
    }

    // ires follows the Fortran PP-subtile ordering.  Convert it through ixyz3
    // rather than assuming the array happens to be laid out in that order.
    int ires_by_location[ns3];
    int subtile_owner[ns3];
    // Keep the static workspace sized for every PP level that kick.f90 can
    // select (ires=2..6); the actual tile grid still uses its local maximum.
    int max_rcp = ratio_sf[5];
    for (int it = 0; it < ns3; ++it) {
        const int sx = ixyz3[it][0] - 1, sy = ixyz3[it][1] - 1, sz = ixyz3[it][2] - 1;
        const int tx = ixyz3[it][3] - 1, ty = ixyz3[it][4] - 1, tz = ixyz3[it][5] - 1;
        if (sx < 0 || sx >= nns || sy < 0 || sy >= nns || sz < 0 || sz >= nns || tx < 0 ||
            tx >= nnt || ty < 0 || ty >= nnt || tz < 0 || tz >= nnt || ires[it] < 2 ||
            ires[it] > 6) {
            fprintf(stderr, "invalid PP subtile metadata at %d (ires=%d)\n", it, ires[it]);
            exit(EXIT_FAILURE);
        }
        const int location = (((tz * nnt + ty) * nnt + tx) * nns + sz) * nns * nns + sy * nns + sx;
        ires_by_location[location] = ires[it];
        subtile_owner[location] = it;
    }
    int64_t (*idx_b_r)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb] =
        (int64_t (*)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb])_idx_b_r;
    const int j_last = nt + 2 * ncb - 1;
    const int k_last = nt + 2 * ncb - 1;
    const int64_t n_used64 = idx_b_r[nnt - 1][nnt - 1][nnt - 1][k_last][j_last];
    if (n_used64 < 0 || n_used64 > INT_MAX || (size_t)n_used64 * 3 > max_vp) {
        fprintf(stderr, "invalid PP particle count %lld (max=%zu)\n", (long long)n_used64,
                max_vp / 3);
        exit(EXIT_FAILURE);
    }
    const int n_used = (int)n_used64;
    std::vector<uint8_t> active_particles((size_t)n_used, 0);

    PPForceWorkspace* ws = get_pp_workspace_static(max_rcp);

    CUDA_CHECK(cudaMemset(ws->d_vp, 0, (size_t)n_used * 3 * sizeof(float)));
    for (int it = 0; it < ns3; ++it) {
        ptotal[it] = 0.0;
        ttotal[it] = 0.0;
    }

    // ires is indexed by PP subtile (ns3), not by the nnt^3 parent tiles.
    // Build a parent tile at a time so the CUDA work remains tile-granular.
    for (int it = 0; it < nnt * nnt * nnt; it++) {
        int tile1[3] = {it % nnt + 1, (it / nnt) % nnt + 1, it / (nnt * nnt) + 1};
        int utile_shift[3] = {0, 0, 0};
        int nc1[3] = {1, 1, 1};
        int nc2[3] = {nt, nt, nt};

        c_pp_force_kernel_tile(ires_by_location, subtile_owner, tile1, nc1, nc2, utile_shift, ratio_sf, apm3, _rhoc,
                               _idx_b_r, _xp, ws->d_vp, mass_p_cdm, a_mid, dt, ws,
                               active_particles.data(), active_particles.size(),
                               &(ws->f2max_teams[it]), ptotal, ttotal);
    }

    float _f2max = 0.0f;

    CUDA_CHECK(cudaDeviceSynchronize());

    for (int it = 0; it < nnt * nnt * nnt; it++) {
        _f2max = fmaxf(_f2max, ws->f2max_teams[it]);
    }

    *f2max = _f2max;

    size_t vp_elems_used = (size_t)n_used * 3;

    accumulate_vp_from_device(ws, _vp, active_particles.data(), vp_elems_used, vmax);

    (void)isort;  // PP updates are disjoint by owner subtile; ordering only affected CPU scheduling.
}
