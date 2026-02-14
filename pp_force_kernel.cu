#include <cuda_runtime.h>
#include <math.h>
#include <math_constants.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/gather.h>
#include <thrust/iterator/reverse_iterator.h>
#include <thrust/replace.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <time.h>

#include <vector>

const int ng = 512;
const int nnt = 4;
const int nns = 4;
const int ratio_cs = 4;
const int nc = ng / ratio_cs;
const int nt = nc / nnt;
const int ntt = nt / nns;
const int ngb = 16;
const int ncb = ngb / ratio_cs;
const int nte = nt + 2 * ncb;
const int ns3 = (nnt * nns) * (nnt * nns) * (nnt * nns);

const int ncore = 32;
const int nteam = 4;
const int nnest = 8;

const bool body_centered_cubic = false;
const int image_buffer = 2;
const int tile_buffer = 3;
const int np_nc = ratio_cs;
const int64_t np_image = (nc * np_nc) * (nc * np_nc) * (nc * np_nc) * (1 + body_centered_cubic);
const int64_t np_image_max = np_image * ((int64_t)nte * 1.0 / (int64_t)nt) *
                             ((int64_t)nte * 1.0 / (int64_t)nt) *
                             ((int64_t)nte * 1.0 / (int64_t)nt) * image_buffer;
//
const int64_t np_tile_max = np_image / (nnt * nnt * nnt) *
                            ((nte * 1.0 / nt) * (nte * 1.0 / nt) * (nte * 1.0 / nt)) * tile_buffer;

const int64_t ishift = -((int64_t)1 << 15);  // 2^15
const double rshift = 0.5 - ishift;
const double x_resolution = 1.0 / ((int64_t)1 << 16);

const double app = 0.06;

#define CUDA_CHECK(call)                                                                        \
    do {                                                                                        \
        cudaError_t err = call;                                                                 \
        if (err != cudaSuccess) {                                                               \
            fprintf(stderr, "CUDA error in %s (%s:%d): %s\n", __FUNCTION__, __FILE__, __LINE__, \
                    cudaGetErrorString(err));                                                   \
            exit(EXIT_FAILURE);                                                                 \
        }                                                                                       \
    } while (0)

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

#pragma unroll
        for (int kk = -1; kk <= 1; kk++) {
#pragma unroll
            for (int jj = -1; jj <= 1; jj++) {
#pragma unroll
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
                            float r = sqrtf(r2);
                            float fpp = F_ra_dev(r, app_local) - F_ra_dev(r, pp_range);
                            float invr = 1.0f / r;
                            ax += fpp * dx * invr;
                            ay += fpp * dy * invr;
                            az += fpp * dz * invr;
                        }
                    }
                }
            }
        }

        af0[outer] += ax;
        af1[outer] += ay;
        af2[outer] += az;
    }
}

__global__ void update_v_kernel(const float* __restrict__ af0, const float* __restrict__ af1,
                                const float* __restrict__ af2, const int64_t* __restrict__ ip_local,
                                const uint32_t* __restrict__ values_sorted,
                                float* __restrict__ vp_flat, int N, float mass2, float a_mid,
                                float dt, float* __restrict__ f2_block,
                                float* __restrict__ vmax_block0, float* __restrict__ vmax_block1,
                                float* __restrict__ vmax_block2) {
    extern __shared__ float s[];
    float* sf2 = s;
    float* sv0 = s + blockDim.x;
    float* sv1 = s + 2 * blockDim.x;
    float* sv2 = s + 3 * blockDim.x;

    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + tid;

    float f2 = 0.f, mv0 = 0.f, mv1 = 0.f, mv2 = 0.f;

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

        mv0 = fabsf(vx);
        mv1 = fabsf(vy);
        mv2 = fabsf(vz);
    }

    sf2[tid] = f2;
    sv0[tid] = mv0;
    sv1[tid] = mv1;
    sv2[tid] = mv2;
    __syncthreads();

    for (int offset = blockDim.x / 2; offset > 0; offset >>= 1) {
        if (tid < offset) {
            sf2[tid] = fmaxf(sf2[tid], sf2[tid + offset]);
            sv0[tid] = fmaxf(sv0[tid], sv0[tid + offset]);
            sv1[tid] = fmaxf(sv1[tid], sv1[tid + offset]);
            sv2[tid] = fmaxf(sv2[tid], sv2[tid + offset]);
        }
        __syncthreads();
    }

    if (tid == 0) {
        f2_block[blockIdx.x] = sf2[0];
        vmax_block0[blockIdx.x] = sv0[0];
        vmax_block1[blockIdx.x] = sv1[0];
        vmax_block2[blockIdx.x] = sv2[0];
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
    float* d_f2_block;
    float* d_v0_block;
    float* d_v1_block;
    float* d_v2_block;
    int n_cap;
    int grid_num_cap;
    int upd_blocks_cap;
};

static inline double elapsed_sec(const struct timespec& t1, const struct timespec& t2) {
    return (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec) * 1.0e-9;
}

static int pp_timing_level() {
    static int level = -1;
    if (level >= 0) return level;
    const char* env = getenv("CUBE_PP_TIMING");
    level = env ? atoi(env) : 0;
    return level;
}

static void init_pp_workspace(PPForceWorkspace* ws, int max_rcp) {
    memset(ws, 0, sizeof(*ws));

    const int upd_threads = 256;
    ws->n_cap = (int)np_tile_max;
    ws->upd_blocks_cap = (ws->n_cap + upd_threads - 1) / upd_threads;

    int stride = ntt * max_rcp + 2 * max_rcp;
    ws->grid_num_cap = stride * stride * stride;

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
    CUDA_CHECK(cudaMalloc(&ws->d_f2_block, (size_t)ws->upd_blocks_cap * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_v0_block, (size_t)ws->upd_blocks_cap * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_v1_block, (size_t)ws->upd_blocks_cap * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&ws->d_v2_block, (size_t)ws->upd_blocks_cap * sizeof(float)));
}

static void destroy_pp_workspace(PPForceWorkspace* ws) {
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
    CUDA_CHECK(cudaFree(ws->d_f2_block));
    CUDA_CHECK(cudaFree(ws->d_v0_block));
    CUDA_CHECK(cudaFree(ws->d_v1_block));
    CUDA_CHECK(cudaFree(ws->d_v2_block));
}

extern "C" void c_pp_force_kernel_tile(int itile, int iapm, int tile1[3], int nc1[3], int nc2[3],
                                       int utile_shift[3], int ratio_sf[], float apm3[], int* _rhoc,
                                       int64_t* _idx_b_r, int16_t* _xp, float* d_vp_flat,
                                       float mass_p_cdm, float a_mid, float dt,
                                       PPForceWorkspace* ws, float* f2max_t_, float vmax_t_[3],
                                       double* ptotal, double* ttotal) {
    struct timespec t_total0, t_total1;
    struct timespec t_reorder0, t_reorder1, t_h2d0, t_h2d1, t_sort0, t_sort1;
    struct timespec t_cell0, t_cell1, t_pp0, t_pp1, t_upd0, t_upd1, t_d2h0, t_d2h1;
    double dt_reorder = 0.0, dt_h2d = 0.0, dt_sort = 0.0;
    double dt_cell = 0.0, dt_pp = 0.0, dt_upd = 0.0, dt_d2h = 0.0;
    const int timing = pp_timing_level();
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

    float* xf0_h = (float*)malloc((size_t)nptile * sizeof(float));
    float* xf1_h = (float*)malloc((size_t)nptile * sizeof(float));
    float* xf2_h = (float*)malloc((size_t)nptile * sizeof(float));
    int64_t* ip_local_h = (int64_t*)malloc((size_t)nptile * sizeof(int64_t));
    uint32_t* keys_h = (uint32_t*)malloc((size_t)nptile * sizeof(uint32_t));
    uint32_t* values_h = (uint32_t*)malloc((size_t)nptile * sizeof(uint32_t));

    int64_t ip1 = 0;
    float xvec[3];

    if (timing >= 2) clock_gettime(CLOCK_MONOTONIC, &t_reorder0);

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

                // 你原来的 nzero 计算（每个cell独立，线程安全）
                int64_t rhoc_i_sum = 0;
                for (int index = i; index < nt + 2 * ncb; index++) {
                    rhoc_i_sum +=
                        rhoc[tile1_offset[2]][tile1_offset[1]][tile1_offset[0]][k][j][index];
                }
                int64_t nzero =
                    idx_b_r[tile1_offset[2]][tile1_offset[1]][tile1_offset[0]][k][j] - rhoc_i_sum;

                for (int l = 0; l < np; l++) {
                    int64_t ip1 = base + l;  // 代替原来的全局 ip1++
                    int64_t ip = nzero + l;  // 读xp用的索引

                    ip_local_h[ip1] = ip;

                    float x0 = (float)((i - ncb) +
                                       ((int16_t)(xp[ip][0] + ishift) + rshift) * x_resolution);
                    float x1 = (float)((j - ncb) +
                                       ((int16_t)(xp[ip][1] + ishift) + rshift) * x_resolution);
                    float x2 = (float)((k - ncb) +
                                       ((int16_t)(xp[ip][2] + ishift) + rshift) * x_resolution);

                    xf0_h[ip1] = ratio_cs * x0;
                    xf1_h[ip1] = ratio_cs * x1;
                    xf2_h[ip1] = ratio_cs * x2;

                    int idx_0 = (int)floorf((float)rcp * (x0 - (utile_shift[0] * ntt))) + rcp;
                    int idx_1 = (int)floorf((float)rcp * (x1 - (utile_shift[1] * ntt))) + rcp;
                    int idx_2 = (int)floorf((float)rcp * (x2 - (utile_shift[2] * ntt))) + rcp;

                    uint32_t idx_ = (uint32_t)idx_2 * (uint32_t)(npgrid + 2 * rcp) *
                                        (uint32_t)(npgrid + 2 * rcp) +
                                    (uint32_t)idx_1 * (uint32_t)(npgrid + 2 * rcp) +
                                    (uint32_t)idx_0;

                    keys_h[ip1] = idx_;
                    values_h[ip1] = (uint32_t)ip1;
                }
            }
        }
    }

    int N = (int)total;
    if (timing >= 2) {
        clock_gettime(CLOCK_MONOTONIC, &t_reorder1);
        dt_reorder = elapsed_sec(t_reorder0, t_reorder1);
        clock_gettime(CLOCK_MONOTONIC, &t_h2d0);
    }
    if (N > ws->n_cap) {
        fprintf(stderr, "N=%d exceeds workspace n_cap=%d (np_tile_max)\n", N, ws->n_cap);
        exit(EXIT_FAILURE);
    }

    CUDA_CHECK(
        cudaMemcpy(ws->d_keys, keys_h, (size_t)N * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(
        cudaMemcpy(ws->d_values, values_h, (size_t)N * sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ws->d_xf0, xf0_h, (size_t)N * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ws->d_xf1, xf1_h, (size_t)N * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ws->d_xf2, xf2_h, (size_t)N * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ws->d_ip_local, ip_local_h, (size_t)N * sizeof(int64_t),
                          cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMemset(ws->d_af0, 0, (size_t)N * sizeof(float)));
    CUDA_CHECK(cudaMemset(ws->d_af1, 0, (size_t)N * sizeof(float)));
    CUDA_CHECK(cudaMemset(ws->d_af2, 0, (size_t)N * sizeof(float)));
    if (timing >= 2) {
        CUDA_CHECK(cudaDeviceSynchronize());
        clock_gettime(CLOCK_MONOTONIC, &t_h2d1);
        dt_h2d = elapsed_sec(t_h2d0, t_h2d1);
        clock_gettime(CLOCK_MONOTONIC, &t_sort0);
    }

    // -------------------------
    // thrust sort_by_key(keys, values)
    // -------------------------
    {
        thrust::device_ptr<uint32_t> kptr(ws->d_keys);
        thrust::device_ptr<uint32_t> vptr(ws->d_values);
        thrust::sort_by_key(kptr, kptr + N, vptr);
    }

    // -------------------------
    // gather xf -> xf_sorted
    // -------------------------
    {
        thrust::device_ptr<uint32_t> map(ws->d_values);
        thrust::device_ptr<float> in0(ws->d_xf0), in1(ws->d_xf1), in2(ws->d_xf2);
        thrust::device_ptr<float> out0(ws->d_xf0_sorted), out1(ws->d_xf1_sorted),
            out2(ws->d_xf2_sorted);
        thrust::gather(map, map + N, in0, out0);
        thrust::gather(map, map + N, in1, out1);
        thrust::gather(map, map + N, in2, out2);
    }
    if (timing >= 2) {
        CUDA_CHECK(cudaDeviceSynchronize());
        clock_gettime(CLOCK_MONOTONIC, &t_sort1);
        dt_sort = elapsed_sec(t_sort0, t_sort1);
        clock_gettime(CLOCK_MONOTONIC, &t_cell0);
    }

    // -------------------------
    // build cell ranges
    // -------------------------
    int grid_num = (npgrid + 2 * rcp) * (npgrid + 2 * rcp) * (npgrid + 2 * rcp);
    if (grid_num > ws->grid_num_cap) {
        fprintf(stderr, "grid_num=%d exceeds workspace grid_num_cap=%d\n", grid_num,
                ws->grid_num_cap);
        exit(EXIT_FAILURE);
    }
    CUDA_CHECK(cudaMemset(ws->d_cell_start, 0xFF, (size_t)grid_num * sizeof(int)));  // -1
    CUDA_CHECK(cudaMemset(ws->d_cell_end, 0xFF, (size_t)grid_num * sizeof(int)));    // -1

    {
        int threads = 256;
        int blocks = (N + threads - 1) / threads;
        build_cell_start_end_kernel<<<blocks, threads>>>(ws->d_keys, N, ws->d_cell_start,
                                                         ws->d_cell_end);
        CUDA_CHECK(cudaGetLastError());

        thrust::device_ptr<int> start_ptr(ws->d_cell_start);
        thrust::device_ptr<int> end_ptr(ws->d_cell_end);

        thrust::replace(start_ptr, start_ptr + grid_num, -1, N);
        thrust::replace(end_ptr, end_ptr + grid_num, -1, N);

        auto rstart_begin = thrust::make_reverse_iterator(start_ptr + grid_num);
        auto rstart_end = thrust::make_reverse_iterator(start_ptr);
        auto rend_begin = thrust::make_reverse_iterator(end_ptr + grid_num);
        auto rend_end = thrust::make_reverse_iterator(end_ptr);

        thrust::inclusive_scan(rstart_begin, rstart_end, rstart_begin, thrust::minimum<int>());
        thrust::inclusive_scan(rend_begin, rend_end, rend_begin, thrust::minimum<int>());
    }
    if (timing >= 2) {
        CUDA_CHECK(cudaDeviceSynchronize());
        clock_gettime(CLOCK_MONOTONIC, &t_cell1);
        dt_cell = elapsed_sec(t_cell0, t_cell1);
        clock_gettime(CLOCK_MONOTONIC, &t_pp0);
    }

    // -------------------------
    // pp force kernel
    // -------------------------
    int real_cells = npgrid * npgrid * npgrid;
    {
        int threads = 128;  // 可调
        pp_force_kernel<<<real_cells, threads>>>(
            ws->d_xf0_sorted, ws->d_xf1_sorted, ws->d_xf2_sorted, npgrid, rcp, pp_range, (float)app,
            ws->d_cell_start, ws->d_cell_end, ws->d_af0, ws->d_af1, ws->d_af2);
        CUDA_CHECK(cudaGetLastError());
    }
    if (timing >= 2) {
        CUDA_CHECK(cudaDeviceSynchronize());
        clock_gettime(CLOCK_MONOTONIC, &t_pp1);
        dt_pp = elapsed_sec(t_pp0, t_pp1);
        clock_gettime(CLOCK_MONOTONIC, &t_upd0);
    }

    // -------------------------
    // update vp + reduce block max
    // -------------------------
    int upd_threads = 256;
    int upd_blocks = (N + upd_threads - 1) / upd_threads;
    if (upd_blocks > ws->upd_blocks_cap) {
        fprintf(stderr, "upd_blocks=%d exceeds workspace upd_blocks_cap=%d\n", upd_blocks,
                ws->upd_blocks_cap);
        exit(EXIT_FAILURE);
    }

    {
        size_t shmem = (size_t)upd_threads * 4 * sizeof(float);
        update_v_kernel<<<upd_blocks, upd_threads, shmem>>>(
            ws->d_af0, ws->d_af1, ws->d_af2, ws->d_ip_local,
            ws->d_values,  // values 已经是 sorted 后的（排序把 values 跟着动了）
            d_vp_flat, N, mass_p_cdm * mass_p_cdm, a_mid, dt, ws->d_f2_block, ws->d_v0_block,
            ws->d_v1_block, ws->d_v2_block);
        CUDA_CHECK(cudaGetLastError());
    }
    if (timing >= 2) {
        CUDA_CHECK(cudaDeviceSynchronize());
        clock_gettime(CLOCK_MONOTONIC, &t_upd1);
        dt_upd = elapsed_sec(t_upd0, t_upd1);
        clock_gettime(CLOCK_MONOTONIC, &t_d2h0);
    }

    CUDA_CHECK(cudaDeviceSynchronize());

    // CPU reduce f2/vmax
    float* f2_h = (float*)malloc((size_t)upd_blocks * sizeof(float));
    float* v0_h = (float*)malloc((size_t)upd_blocks * sizeof(float));
    float* v1_h = (float*)malloc((size_t)upd_blocks * sizeof(float));
    float* v2_h = (float*)malloc((size_t)upd_blocks * sizeof(float));

    CUDA_CHECK(cudaMemcpy(f2_h, ws->d_f2_block, (size_t)upd_blocks * sizeof(float),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(v0_h, ws->d_v0_block, (size_t)upd_blocks * sizeof(float),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(v1_h, ws->d_v1_block, (size_t)upd_blocks * sizeof(float),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(v2_h, ws->d_v2_block, (size_t)upd_blocks * sizeof(float),
                          cudaMemcpyDeviceToHost));
    if (timing >= 2) {
        clock_gettime(CLOCK_MONOTONIC, &t_d2h1);
        dt_d2h = elapsed_sec(t_d2h0, t_d2h1);
    }

    float f2max = 0.f, vmax0 = 0.f, vmax1 = 0.f, vmax2 = 0.f;
    for (int b = 0; b < upd_blocks; b++) {
        f2max = fmaxf(f2max, f2_h[b]);
        vmax0 = fmaxf(vmax0, v0_h[b]);
        vmax1 = fmaxf(vmax1, v1_h[b]);
        vmax2 = fmaxf(vmax2, v2_h[b]);
    }

    *f2max_t_ = f2max;
    vmax_t_[0] = vmax0;
    vmax_t_[1] = vmax1;
    vmax_t_[2] = vmax2;

    // -------------------------
    // free
    // -------------------------
    free(xf0_h);
    free(xf1_h);
    free(xf2_h);
    free(ip_local_h);
    free(keys_h);
    free(values_h);
    free(f2_h);
    free(v0_h);
    free(v1_h);
    free(v2_h);

    ptotal[itile] = (double)nptile_real;

    clock_gettime(CLOCK_MONOTONIC, &t_total1);
    ttotal[itile] = elapsed_sec(t_total0, t_total1);
    if (timing >= 2) {
        fprintf(
            stderr,
            "[PP-TILE] itile=%d iapm=%d N=%d nptile_real=%d total=%.6fs reorder=%.6fs h2d=%.6fs "
            "sort=%.6fs cell=%.6fs pp=%.6fs update=%.6fs d2h=%.6fs\n",
            itile, iapm, N, nptile_real, ttotal[itile], dt_reorder, dt_h2d, dt_sort, dt_cell, dt_pp,
            dt_upd, dt_d2h);
    }
}

extern "C" void c_pp_force_kernel(int isort[ns3], int ires[ns3], int ixyz3[ns3][6], float apm3[7],
                                  int ratio_sf[7], int* _rhoc, int64_t* _idx_b_r, int16_t* _xp,
                                  float* _vp, float mass_p_cdm, float a_mid, float dt, float* f2max,
                                  float vmax[3], double* ptotal, double* ttotal) {
    struct timespec t_total0, t_total1, t_h2d0, t_h2d1, t_tile0, t_tile1, t_d2h0, t_d2h1;
    double dt_h2d = 0.0, dt_tiles = 0.0, dt_d2h = 0.0;
    const int timing = pp_timing_level();
    if (timing >= 1) clock_gettime(CLOCK_MONOTONIC, &t_total0);
    float _f2max = 0.0f;
    float _vmax0 = 0.0f, _vmax1 = 0.0f, _vmax2 = 0.0f;

    // move vp to device memory once
    float* d_vp = nullptr;
    size_t N_vp = (size_t)np_image_max * 3 * sizeof(float);
    CUDA_CHECK(cudaMalloc(&d_vp, N_vp));
    if (timing >= 1) clock_gettime(CLOCK_MONOTONIC, &t_h2d0);
    CUDA_CHECK(cudaMemcpy(d_vp, _vp, N_vp, cudaMemcpyHostToDevice));
    if (timing >= 1) {
        CUDA_CHECK(cudaDeviceSynchronize());
        clock_gettime(CLOCK_MONOTONIC, &t_h2d1);
        dt_h2d = elapsed_sec(t_h2d0, t_h2d1);
        clock_gettime(CLOCK_MONOTONIC, &t_tile0);
    }

    int max_rcp = ratio_sf[0];
    for (int i = 1; i < 7; i++) {
        if (ratio_sf[i] > max_rcp) max_rcp = ratio_sf[i];
    }
    PPForceWorkspace ws;
    init_pp_workspace(&ws, max_rcp);

    for (int it = 0; it < nnt * nnt * nnt; it++) {
        int itile = it;
        int iapm = ires[itile];

        int tile1[3] = {it / (nnt * nnt) + 1, (it / nnt) % nnt + 1, it % nnt + 1};
        int utile_shift[3] = {0, 0, 0};
        int nc1[3] = {1, 1, 1};
        int nc2[3] = {nt, nt, nt};

        float f2max_team = 0.f;
        float vmax_team[3] = {0.f, 0.f, 0.f};

        c_pp_force_kernel_tile(itile, iapm, tile1, nc1, nc2, utile_shift, ratio_sf, apm3, _rhoc,
                               _idx_b_r, _xp, d_vp, mass_p_cdm, a_mid, dt, &ws, &f2max_team,
                               vmax_team, ptotal, ttotal);

        _f2max = fmaxf(_f2max, f2max_team);
        _vmax0 = fmaxf(_vmax0, vmax_team[0]);
        _vmax1 = fmaxf(_vmax1, vmax_team[1]);
        _vmax2 = fmaxf(_vmax2, vmax_team[2]);
    }
    if (timing >= 1) {
        clock_gettime(CLOCK_MONOTONIC, &t_tile1);
        dt_tiles = elapsed_sec(t_tile0, t_tile1);
        clock_gettime(CLOCK_MONOTONIC, &t_d2h0);
    }

    *f2max = _f2max;
    vmax[0] = _vmax0;
    vmax[1] = _vmax1;
    vmax[2] = _vmax2;

    CUDA_CHECK(cudaMemcpy(_vp, d_vp, N_vp, cudaMemcpyDeviceToHost));
    if (timing >= 1) {
        CUDA_CHECK(cudaDeviceSynchronize());
        clock_gettime(CLOCK_MONOTONIC, &t_d2h1);
        dt_d2h = elapsed_sec(t_d2h0, t_d2h1);
        clock_gettime(CLOCK_MONOTONIC, &t_total1);
        fprintf(stderr,
                "[PP] total=%.6fs vp_h2d=%.6fs tiles=%.6fs vp_d2h=%.6fs f2max=%.6e vmax=(%.6e, "
                "%.6e, %.6e)\n",
                elapsed_sec(t_total0, t_total1), dt_h2d, dt_tiles, dt_d2h, (double)(*f2max),
                (double)vmax[0], (double)vmax[1], (double)vmax[2]);
    }
    destroy_pp_workspace(&ws);
    CUDA_CHECK(cudaFree(d_vp));
}
