#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int ng = 128;
const int nnt = 4;
const int nns = 4;
const int ratio_cs = 4;
const int nc = ng / ratio_cs;
const int nt = nc / nnt;
const int ntt = nt / nns;
const int ngb = 16;
const int ncb = ngb / ratio_cs;
const int ns3 = (nnt * nns) * (nnt * nns) * (nnt * nns);

const int ncore = 32;
const int nteam = 8;
const int nnest = 4;

const double pi = 3.14159265358979323846;

const int64_t ishift = -((int64_t)1 << 15);  // 2^15
const double rshift = 0.5 - ishift;
const double x_resolution = 1.0 / ((int64_t)1 << 16);

const double app = 0.06;

float F_ra(float r, float apm) {
    float ep;
    float f_ra = 0.0f;

    ep = 2.0f * r / apm;

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

void radix_sort_u32_pair(uint32_t* keys, uint32_t* vals, int N) {
    uint32_t* keys_tmp = (uint32_t*)malloc(N * sizeof(uint32_t));
    uint32_t* vals_tmp = (uint32_t*)malloc(N * sizeof(uint32_t));

    const int BITS = 8;
    const int BUCKETS = 1 << BITS;  // 256
    const int PASSES = 32 / BITS;   // 4

    uint32_t count[BUCKETS];

    for (int pass = 0; pass < PASSES; pass++) {
        int shift = pass * BITS;

        // ----------------
        // 1) histogram
        // ----------------
        memset(count, 0, sizeof(count));

        for (int i = 0; i < N; i++) {
            uint32_t bucket = (keys[i] >> shift) & (BUCKETS - 1);
            count[bucket]++;
        }

        // ----------------
        // 2) prefix sum
        // ----------------
        uint32_t sum = 0;
        for (int b = 0; b < BUCKETS; b++) {
            uint32_t c = count[b];
            count[b] = sum;
            sum += c;
        }

        // ----------------
        // 3) scatter (stable)
        // ----------------
        for (int i = 0; i < N; i++) {
            uint32_t k = keys[i];
            uint32_t v = vals[i];

            uint32_t bucket = (k >> shift) & (BUCKETS - 1);
            uint32_t pos = count[bucket]++;

            keys_tmp[pos] = k;
            vals_tmp[pos] = v;
        }

        // ----------------
        // 4) swap buffers
        // ----------------
        memcpy(keys, keys_tmp, N * sizeof(uint32_t));
        memcpy(vals, vals_tmp, N * sizeof(uint32_t));
    }

    free(keys_tmp);
    free(vals_tmp);
}

void c_pp_force_kernel_subtile(int itile, int iapm, int tile1[3], int nc1[3], int nc2[3],
                               int utile_shift[3], int ratio_sf[], float apm3[],
                               int* _rhoc,         // [nnt, nnt, nnt, nt+2*ncb, nt+2*ncb, nt+2*ncb]
                               int64_t* _idx_b_r,  // [nnt, nnt, nnt, nt+2*ncb, nt+2*ncb]
                               int16_t* _xp,       // [nplocal, 3]
                               float* _vp,         // [nplocal, 3]
                               float mass_p_cdm, float a_mid, float dt, float* f2max_t_,
                               float vmax_t_[3], double* ptotal, double* ttotal) {
    struct timespec t1, t2;
    clock_gettime(CLOCK_MONOTONIC, &t1);

    float pp_range = apm3[iapm - 1];
    int rcp = ratio_sf[iapm - 1];
    int npgrid = ntt * rcp;

    int (*rhoc)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb][nt + 2 * ncb] =
        (int (*)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb][nt + 2 * ncb]) _rhoc;
    int64_t (*idx_b_r)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb] =
        (int64_t (*)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb]) _idx_b_r;

    int16_t (*xp)[3] = (int16_t (*)[3])_xp;
    float (*vp)[3] = (float (*)[3])_vp;

    int nc1_offset[3];  // nc1 + ncb - 1
    int nc2_offset[3];  // nc2 + ncb - 1
    for (int d = 0; d < 3; d++) {
        nc1_offset[d] = nc1[d] + ncb - 1;
        nc2_offset[d] = nc2[d] + ncb - 1;
    }

    int tile1_offset[3];  // tile1 - 1
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

    float* _xf = (float*)malloc(nptile * 3 * sizeof(float));
    float* _vf = (float*)malloc(nptile * 3 * sizeof(float));
    float (*xf)[3] = (float (*)[3])_xf;
    float (*vf)[3] = (float (*)[3])_vf;

    float* _af = (float*)malloc(3 * nptile * sizeof(float));
    int64_t* ip_local = (int64_t*)malloc(nptile * sizeof(int64_t));

    memset(_af, 0, 3 * nptile * sizeof(float));
    float (*af)[nptile] = (float (*)[nptile])_af;

    uint32_t* keys = (uint32_t*)malloc(nptile * sizeof(uint32_t));
    uint32_t* values = (uint32_t*)malloc(nptile * sizeof(uint32_t));

    int64_t ip1 = 0;
    float xvec[3];

    for (int k = nc1_offset[2] - 1; k <= nc2_offset[2] + 1; k++) {
        for (int j = nc1_offset[1] - 1; j <= nc2_offset[1] + 1; j++) {
            for (int i = nc1_offset[0] - 1; i <= nc2_offset[0] + 1; i++) {
                int np = rhoc[tile1_offset[2]][tile1_offset[1]][tile1_offset[0]][k][j][i];
                int64_t rhoc_i_sum = 0;
                for (int index = i; index < nt + 2 * ncb; index++) {
                    rhoc_i_sum +=
                        rhoc[tile1_offset[2]][tile1_offset[1]][tile1_offset[0]][k][j][index];
                }
                int64_t nzero =
                    idx_b_r[tile1_offset[2]][tile1_offset[1]][tile1_offset[0]][k][j] - rhoc_i_sum;
                for (int l = 0; l < np; l++) {
                    int64_t ip = nzero + l;
                    ip_local[ip1] = ip;
                    xvec[0] = (float)((i - ncb) +
                                      ((int16_t)(xp[ip][0] + ishift) + rshift) * x_resolution);
                    xvec[1] = (float)((j - ncb) +
                                      ((int16_t)(xp[ip][1] + ishift) + rshift) * x_resolution);
                    xvec[2] = (float)((k - ncb) +
                                      ((int16_t)(xp[ip][2] + ishift) + rshift) * x_resolution);
                    xf[ip1][0] = ratio_cs * xvec[0];
                    xf[ip1][1] = ratio_cs * xvec[1];
                    xf[ip1][2] = ratio_cs * xvec[2];

                    vf[ip1][0] = vp[ip][0];
                    vf[ip1][1] = vp[ip][1];
                    vf[ip1][2] = vp[ip][2];
                    int idx_0 = floor(rcp * (xvec[0] - (utile_shift[0] * ntt))) + rcp;
                    int idx_1 = floor(rcp * (xvec[1] - (utile_shift[1] * ntt))) + rcp;
                    int idx_2 = floor(rcp * (xvec[2] - (utile_shift[2] * ntt))) + rcp;

                    uint32_t idx_ = (idx_2) * (npgrid + 2 * rcp) * (npgrid + 2 * rcp) +
                                    (idx_1) * (npgrid + 2 * rcp) + (idx_0);

                    keys[ip1] = idx_;
                    values[ip1] = ip1;

                    ip1 += 1;
                }
            }
        }
    }

    int N = (int)ip1;

    radix_sort_u32_pair(keys, values, N);

    float* _xf_sorted = (float*)malloc(3 * N * sizeof(float));
    float* _vf_sorted = (float*)malloc(3 * N * sizeof(float));

    float (*xf_sorted)[N] = (float (*)[N])_xf_sorted;
    float (*vf_sorted)[N] = (float (*)[N])_vf_sorted;

    // gather xf vf
    for (int i = 0; i < N; i++) {
        int ori_idx = values[i];
        xf_sorted[0][i] = xf[ori_idx][0];
        xf_sorted[1][i] = xf[ori_idx][1];
        xf_sorted[2][i] = xf[ori_idx][2];

        vf_sorted[0][i] = vf[ori_idx][0];
        vf_sorted[1][i] = vf[ori_idx][1];
        vf_sorted[2][i] = vf[ori_idx][2];
    }

    // build cell start index
    int grid_num = (npgrid + 2 * rcp) * (npgrid + 2 * rcp) * (npgrid + 2 * rcp);
    int* cell_start = (int*)malloc((grid_num + 1) * sizeof(int));
    memset(cell_start, -1, (grid_num + 1) * sizeof(int));
    for (int i = 0; i < N; i++) {
        uint32_t idx_ = keys[i];
        if (cell_start[idx_] == -1) {
            cell_start[idx_] = i;
        }
    }
    cell_start[grid_num] = N;

    // Fix no particle cell start index to the next cell with particles
    for (int i = grid_num - 1; i >= 0; i--) {
        if (cell_start[i] == -1) {
            cell_start[i] = cell_start[i + 1];
        }
    }

// compute pp force
#pragma omp parallel for default(shared) num_threads(nnest)
    for (int kji = 0; kji < npgrid * npgrid * npgrid; kji++) {
        int k = kji / (npgrid * npgrid);
        int j = (kji / npgrid) % npgrid;
        int i = kji % npgrid;

        int idx_ = (k + rcp) * (npgrid + 2 * rcp) * (npgrid + 2 * rcp) +
                   (j + rcp) * (npgrid + 2 * rcp) + (i + rcp);
        // outer particle start_idx
        int outer_start_idx = cell_start[idx_];
        int outer_end_idx = cell_start[idx_ + 1];

        for (int kk = -1; kk <= 1; kk++) {
            for (int jj = -1; jj <= 1; jj++) {
                for (int ii = -1; ii <= 1; ii++) {
                    int nidx_ = (k + kk + rcp) * (npgrid + 2 * rcp) * (npgrid + 2 * rcp) +
                                (j + jj + rcp) * (npgrid + 2 * rcp) + (i + ii + rcp);
                    int inner_start_idx = cell_start[nidx_];
                    int inner_end_idx = cell_start[nidx_ + 1];
                    for (int outer_idx = outer_start_idx; outer_idx < outer_end_idx; outer_idx++) {
                        float x_outer = xf_sorted[0][outer_idx];
                        float y_outer = xf_sorted[1][outer_idx];
                        float z_outer = xf_sorted[2][outer_idx];

                        for (int inner_idx = inner_start_idx; inner_idx < inner_end_idx;
                             inner_idx++) {
                            float dx = xf_sorted[0][inner_idx] - x_outer;
                            float dy = xf_sorted[1][inner_idx] - y_outer;
                            float dz = xf_sorted[2][inner_idx] - z_outer;

                            float r2 = dx * dx + dy * dy + dz * dz;
                            if (r2 > 0.0f && r2 < pp_range * pp_range) {
                                float rmag = sqrtf(r2);
                                float fpp = F_ra(rmag, app) - F_ra(rmag, pp_range);
                                af[0][outer_idx] += fpp * dx / rmag;
                                af[1][outer_idx] += fpp * dy / rmag;
                                af[2][outer_idx] += fpp * dz / rmag;
                            }
                        }
                    }
                }
            }
        }
    }

    *f2max_t_ = 0.0f;
    vmax_t_[0] = 0.0f;
    vmax_t_[1] = 0.0f;
    vmax_t_[2] = 0.0f;
    int ncount = 0;

    // update vf and store back vp
    for (int k = 0; k < npgrid; k++) {
        for (int j = 0; j < npgrid; j++) {
            for (int i = 0; i < npgrid; i++) {
                int idx_ = (k + rcp) * (npgrid + 2 * rcp) * (npgrid + 2 * rcp) +
                           (j + rcp) * (npgrid + 2 * rcp) + (i + rcp);

                int start_idx = cell_start[idx_];
                int end_idx = cell_start[idx_ + 1];

                for (int idx = start_idx; idx < end_idx; idx++) {
                    ncount += 1;
                    af[0][idx] = af[0][idx] * mass_p_cdm * mass_p_cdm;
                    af[1][idx] = af[1][idx] * mass_p_cdm * mass_p_cdm;
                    af[2][idx] = af[2][idx] * mass_p_cdm * mass_p_cdm;
                    *f2max_t_ = fmaxf(*f2max_t_, af[0][idx] * af[0][idx] + af[1][idx] * af[1][idx] +
                                                     af[2][idx] * af[2][idx]);
                    vp[ip_local[values[idx]]][0] += af[0][idx] * a_mid * dt / 6 / pi;
                    vp[ip_local[values[idx]]][1] += af[1][idx] * a_mid * dt / 6 / pi;
                    vp[ip_local[values[idx]]][2] += af[2][idx] * a_mid * dt / 6 / pi;
                    vmax_t_[0] = fmaxf(vmax_t_[0], fabsf(vp[ip_local[values[idx]]][0]));
                    vmax_t_[1] = fmaxf(vmax_t_[1], fabsf(vp[ip_local[values[idx]]][1]));
                    vmax_t_[2] = fmaxf(vmax_t_[2], fabsf(vp[ip_local[values[idx]]][2]));
                }
            }
        }
    }

    free(cell_start);
    free(_xf_sorted);
    free(_vf_sorted);
    free(keys);
    free(values);

    free(_xf);
    free(_vf);
    free(_af);
    free(ip_local);

    clock_gettime(CLOCK_MONOTONIC, &t2);

    ptotal[itile] = ncount;
    ttotal[itile] = (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec) * 1.0e-9;

    return;
}

void c_pp_force_kernel(int isort[ns3], int ires[ns3], int ixyz3[ns3][6], float apm3[7],
                                  int ratio_sf[7],
                                  int* _rhoc,  // [nnt, nnt, nnt, nt+2*ncb, nt+2*ncb, nt+2*ncb]
                                  int64_t* _idx_b_r,  // [nnt, nnt, nnt, nt+2*ncb, nt+2*ncb]
                                  int16_t* _xp,       // [nplocal, 3]
                                  float* _vp,         // [nplocal, 3]
                                  float mass_p_cdm, float a_mid, float dt, float* f2max,
                                  float vmax[3], double* ptotal, double* ttotal) {
    float _f2max = 0.0f;
#pragma omp parallel for default(shared) num_threads(8) schedule(dynamic) reduction(max : _f2max, vmax[ : 3])
    for (int it = 0; it < ns3; it++) {
        int itile = isort[it] - 1;
        int iapm = ires[itile];

        int tile1[3] = {ixyz3[itile][3], ixyz3[itile][4], ixyz3[itile][5]};
        int utile_shift[3] = {ixyz3[itile][0] - 1, ixyz3[itile][1] - 1, ixyz3[itile][2] - 1};
        int nc1[3] = {(ixyz3[itile][0] - 1) * ntt + 1, (ixyz3[itile][1] - 1) * ntt + 1,
                      (ixyz3[itile][2] - 1) * ntt + 1};
        int nc2[3] = {ixyz3[itile][0] * ntt, ixyz3[itile][1] * ntt, ixyz3[itile][2] * ntt};

        float f2max_team;
        float vmax_team[3];

        c_pp_force_kernel_subtile(itile, iapm, tile1, nc1, nc2, utile_shift, ratio_sf, apm3, _rhoc,
                                  _idx_b_r, _xp, _vp, mass_p_cdm, a_mid, dt, &f2max_team, vmax_team,
                                  ptotal, ttotal);

        _f2max = fmaxf(_f2max, f2max_team);
        vmax[0] = fmaxf(vmax[0], vmax_team[0]);
        vmax[1] = fmaxf(vmax[1], vmax_team[1]);
        vmax[2] = fmaxf(vmax[2], vmax_team[2]);
    }
    *f2max = _f2max;
    return;
}