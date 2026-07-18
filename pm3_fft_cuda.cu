#include "cuda_parameters_generated.h"
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include <cufft.h>
#include <iostream>
#include <mutex>
#include <vector>

extern "C" void c_select_gpu_for_image();
extern "C" float *c_gpu_velocity_acquire(int, int);
extern "C" void c_gpu_velocity_wait_default_stream();
#define CK(x)                                                                  \
  do {                                                                         \
    cudaError_t e = (x);                                                       \
    if (e != cudaSuccess) {                                                    \
      fprintf(stderr, "PM3 CUDA: %s\n", cudaGetErrorString(e));                \
      abort();                                                                 \
    }                                                                          \
  } while (0)
#define FK(x)                                                                  \
  do {                                                                         \
    cufftResult e = (x);                                                       \
    if (e != CUFFT_SUCCESS) {                                                  \
      fprintf(stderr, "PM3 cuFFT: %d\n", (int)e);                              \
      abort();                                                                 \
    }                                                                          \
  } while (0)
namespace {
constexpr int T = 256, ng = CUBE_NG, nnt = CUBE_NNT, nns = CUBE_NNS,
              rcs = CUBE_RATIO_CS;
constexpr int nt = ng / rcs / nnt, ntt = nt / nns, ncb = CUBE_NGB / rcs;
struct CellDesc {
  uint32_t base, count;
  int16_t cell[3], utile[3];
  uint16_t batch;
  uint8_t active;
};
struct Job {
  int it, level, tile[3], utile[3], nc1[3], nc2[3];
};
struct FFTPlan {
  int n = 0, nb = 0;
  cufftHandle fw = 0, bw = 0;
  size_t fw_work = 0, bw_work = 0;
  uint64_t last_use = 0;
};
struct Green {
  const float *host = nullptr;
  size_t count = 0;
  float *device = nullptr;
};
struct WS {
  float *rho = nullptr, *force = nullptr, *f2 = nullptr;
  CellDesc *p = nullptr;
  size_t rho_cap = 0, force_cap = 0, f2_cap = 0, p_cap = 0;
  std::vector<FFTPlan> plans;
  std::vector<Green> greens;
  void *fft_work = nullptr;
  size_t fft_work_cap = 0;
  std::vector<std::vector<Job>> host_job_groups;
  std::vector<CellDesc> host_cells;
  std::vector<float> host_f2;
  uint64_t plan_clock = 0;
} w;
std::mutex m;
int16_t *d_xp = nullptr;
size_t d_xp_capacity = 0;
struct Profile {
  bool on = false;
  int batches = 0, particles = 0;
  int fft_jobs = 0, fft_slots = 0;
  int reserve_grows = 0, plan_hits = 0, plan_misses = 0, plan_evictions = 0;
  int green_hits = 0, green_misses = 0;
  size_t reserve_bytes = 0;
  double schedule_ms = 0, pack_ms = 0, ensure_ms = 0, reserve_ms = 0;
  double plan_create_ms = 0, plan_evict_ms = 0, green_upload_ms = 0;
  double velocity_init_ms = 0, gpu_select_ms = 0;
  float h2d_ms = 0, deposit_ms = 0, fft_forward_ms = 0, green_ms = 0,
        fft_backward_ms = 0, gradient_ms = 0, interpolate_ms = 0, reduce_ms = 0,
        clear_ms = 0, xp_h2d_ms = 0;
  cudaEvent_t start = nullptr, stop = nullptr;
};
template <class F> void gpu_time(Profile *p, float &total, F &&fn) {
  if (!p || !p->on) {
    fn();
    return;
  }
  CK(cudaEventRecord(p->start));
  fn();
  CK(cudaEventRecord(p->stop));
  CK(cudaEventSynchronize(p->stop));
  float ms = 0;
  CK(cudaEventElapsedTime(&ms, p->start, p->stop));
  total += ms;
}
void drop() {
  for (const FFTPlan &plan : w.plans) {
    if (plan.fw)
      FK(cufftDestroy(plan.fw));
    if (plan.bw)
      FK(cufftDestroy(plan.bw));
  }
  if (w.rho)
    CK(cudaFree(w.rho));
  for (const Green &green : w.greens)
    if (green.device)
      CK(cudaFree(green.device));
  if (w.force)
    CK(cudaFree(w.force));
  if (w.f2)
    CK(cudaFree(w.f2));
  if (w.p)
    CK(cudaFree(w.p));
  if (w.fft_work)
    CK(cudaFree(w.fft_work));
  w = WS{};
}
template <class T>
void reserve_device(T *&ptr, size_t &capacity, size_t required,
                    Profile *profile = nullptr) {
  if (required <= capacity)
    return;
  const size_t new_capacity =
      std::max(required, capacity + capacity / 2 + static_cast<size_t>(1024));
  const auto start = std::chrono::steady_clock::now();
  T *replacement = nullptr;
  CK(cudaMalloc(&replacement, new_capacity * sizeof(T)));
  if (ptr)
    CK(cudaFree(ptr));
  ptr = replacement;
  capacity = new_capacity;
  if (profile && profile->on) {
    profile->reserve_grows++;
    profile->reserve_bytes += new_capacity * sizeof(T);
    profile->reserve_ms += std::chrono::duration<double, std::milli>(
                               std::chrono::steady_clock::now() - start)
                               .count();
  }
}
void sync_xp(const int16_t *xp, size_t particle_count, Profile *profile) {
  const size_t elements = particle_count * 3;
  if (!elements)
    return;
  reserve_device(d_xp, d_xp_capacity, elements);
  gpu_time(profile, profile->xp_h2d_ms, [&] {
    CK(cudaMemcpy(d_xp, xp, elements * sizeof(int16_t), cudaMemcpyHostToDevice));
  });
}
Green *green_for(const float *host, size_t count, Profile *profile) {
  for (Green &green : w.greens)
    if (green.host == host && green.count == count) {
      if (profile && profile->on)
        profile->green_hits++;
      return &green;
    }
  if (profile && profile->on)
    profile->green_misses++;
  const auto start = std::chrono::steady_clock::now();
  Green green;
  green.host = host;
  green.count = count;
  CK(cudaMalloc(&green.device, count * sizeof(float)));
  CK(cudaMemcpy(green.device, host, count * sizeof(float),
                cudaMemcpyHostToDevice));
  w.greens.push_back(green);
  if (profile && profile->on)
    profile->green_upload_ms += std::chrono::duration<double, std::milli>(
                                    std::chrono::steady_clock::now() - start)
                                    .count();
  return &w.greens.back();
}
void reserve_fft_work(size_t required) {
  if (required <= w.fft_work_cap)
    return;
  // Plans execute serially in the PM3 default stream.  They can therefore
  // share one work area sized for the largest cached plan.
  void *replacement = nullptr;
  CK(cudaMalloc(&replacement, required));
  for (const FFTPlan &plan : w.plans) {
    if (plan.fw_work)
      FK(cufftSetWorkArea(plan.fw, replacement));
    if (plan.bw_work)
      FK(cufftSetWorkArea(plan.bw, replacement));
  }
  if (w.fft_work)
    CK(cudaFree(w.fft_work));
  w.fft_work = replacement;
  w.fft_work_cap = required;
}
FFTPlan *plan_for(int n, int nb, Profile *profile) {
  for (FFTPlan &plan : w.plans)
    if (plan.n == n && plan.nb == nb) {
      plan.last_use = ++w.plan_clock;
      if (profile && profile->on)
        profile->plan_hits++;
      return &plan;
    }
  if (profile && profile->on)
    profile->plan_misses++;
  // A bounded cache keeps cuFFT work areas from growing with every possible
  // occupancy pattern while preserving the active bucket working set.
  constexpr size_t max_cached_plans = 40;
  if (w.plans.size() == max_cached_plans) {
    auto old = std::min_element(
        w.plans.begin(), w.plans.end(),
        [](const FFTPlan &a, const FFTPlan &b) { return a.last_use < b.last_use; });
    const auto evict_start = std::chrono::steady_clock::now();
    FK(cufftDestroy(old->fw));
    FK(cufftDestroy(old->bw));
    if (profile && profile->on) {
      profile->plan_evictions++;
      profile->plan_evict_ms += std::chrono::duration<double, std::milli>(
                                    std::chrono::steady_clock::now() - evict_start)
                                    .count();
    }
    *old = w.plans.back();
    w.plans.pop_back();
  }
  FFTPlan plan;
  plan.n = n;
  plan.nb = nb;
  plan.last_use = ++w.plan_clock;
  const size_t rr = (size_t)(n + 2) * n * n;
  const size_t kk = (size_t)(n / 2 + 1) * n * n;
  int d[3] = {n, n, n}, ie[3] = {n, n, n + 2}, oe[3] = {n, n, n / 2 + 1};
  const auto create_start = std::chrono::steady_clock::now();
  FK(cufftCreate(&plan.fw));
  FK(cufftCreate(&plan.bw));
  FK(cufftSetAutoAllocation(plan.fw, 0));
  FK(cufftSetAutoAllocation(plan.bw, 0));
  FK(cufftMakePlanMany(plan.fw, 3, d, ie, 1, rr, oe, 1, kk, CUFFT_R2C, nb,
                        &plan.fw_work));
  FK(cufftMakePlanMany(plan.bw, 3, d, oe, 1, kk, ie, 1, rr, CUFFT_C2R, nb,
                        &plan.bw_work));
  reserve_fft_work(std::max(plan.fw_work, plan.bw_work));
  if (plan.fw_work)
    FK(cufftSetWorkArea(plan.fw, w.fft_work));
  if (plan.bw_work)
    FK(cufftSetWorkArea(plan.bw, w.fft_work));
  if (profile && profile->on)
    profile->plan_create_ms += std::chrono::duration<double, std::milli>(
                                   std::chrono::steady_clock::now() - create_start)
                                   .count();
  w.plans.push_back(plan);
  return &w.plans.back();
}
FFTPlan *ensure(int n, int np, int nb, size_t count, Profile *profile) {
  const size_t rr = (size_t)(n + 2) * n * n;
  const size_t ff = (size_t)(np + 2) * (np + 2) * (np + 2);
  reserve_device(w.rho, w.rho_cap, (size_t)nb * rr, profile);
  reserve_device(w.force, w.force_cap, (size_t)nb * 3 * ff, profile);
  reserve_device(w.f2, w.f2_cap, (size_t)nb, profile);
  reserve_device(w.p, w.p_cap, count, profile);
  return plan_for(n, nb, profile);
}
__device__ inline void tsc(float x, int &c, float q[3]) {
  c = (int)floorf(x) + 1;
  float f = x - floorf(x);
  q[0] = .5f * (1 - f) * (1 - f);
  q[2] = .5f * f * f;
  q[1] = 1 - q[0] - q[2];
}
__device__ inline void mx(float *p, float v) {
  int *ip = (int *)p, o = *ip, n = __float_as_int(v);
  while (o < n) {
    int a = o;
    o = atomicCAS(ip, a, n);
    if (o == a)
      break;
  }
}
__global__ void dep(const CellDesc *p, int N, const int16_t *xp, float *r,
                    int n, int np, int rcp, int fb, float mass, size_t stride) {
  int q = blockIdx.x;
  if (q >= N)
    return;
  const CellDesc d = p[q];
  int lo = 1 - fb;
  float *R = r + d.batch * stride;
  for (uint32_t l = threadIdx.x; l < d.count; l += blockDim.x) {
    uint32_t ip = d.base + l;
    int x, y, z;
    float a[3], b[3], c[3];
    float px = ((int16_t)(xp[3 * ip] - (int16_t)0x8000) + 32768.5f) / 65536.f;
    float py =
        ((int16_t)(xp[3 * ip + 1] - (int16_t)0x8000) + 32768.5f) / 65536.f;
    float pz =
        ((int16_t)(xp[3 * ip + 2] - (int16_t)0x8000) + 32768.5f) / 65536.f;
    tsc(((d.cell[0] - ncb) + px) * (rcs * rcp) - d.utile[0] * np, x, a);
    tsc(((d.cell[1] - ncb) + py) * (rcs * rcp) - d.utile[1] * np, y, b);
    tsc(((d.cell[2] - ncb) + pz) * (rcs * rcp) - d.utile[2] * np, z, c);
    for (int k = -1; k < 2; k++)
      for (int j = -1; j < 2; j++)
        for (int i = -1; i < 2; i++) {
          int X = x + i - lo, Y = y + j - lo, Z = z + k - lo;
          if ((unsigned)X < (unsigned)n && (unsigned)Y < (unsigned)n &&
              (unsigned)Z < (unsigned)n)
            atomicAdd(R + ((size_t)Z * n + Y) * (n + 2) + X,
                      mass * a[i + 1] * b[j + 1] * c[k + 1]);
        }
  }
}
__global__ void green(cufftComplex *a, const float *g, size_t nk,
                      size_t total) {
  size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
  if (i < total) {
    float v = g[i % nk];
    a[i].x *= v;
    a[i].y *= v;
  }
}
__global__ void grad(const float *r, float *f, float *f2, int n, int np, int fb,
                     int nb, size_t rs, size_t fs) {
  size_t q = (size_t)blockIdx.x * blockDim.x + threadIdx.x,
         tot = (size_t)nb * fs;
  if (q >= tot)
    return;
  int B = q / fs;
  size_t u = q - B * fs;
  int nf = np + 2, x = u % nf, y = (u / nf) % nf, z = u / (nf * nf),
      lo = 1 - fb, s = n + 2;
  const float *R = r + B * rs;
  auto A = [&](int X, int Y, int Z) {
    return R[((size_t)(Z - lo) * n + (Y - lo)) * s + (X - lo)];
  };
  float h = 1.f / ((float)n * n * n),
        a = (A(x - 2, y, z) - 8 * A(x - 1, y, z) + 8 * A(x + 1, y, z) -
             A(x + 2, y, z)) *
            h / 12.f,
        b = (A(x, y - 2, z) - 8 * A(x, y - 1, z) + 8 * A(x, y + 1, z) -
             A(x, y + 2, z)) *
            h / 12.f,
        c = (A(x, y, z - 2) - 8 * A(x, y, z - 1) + 8 * A(x, y, z + 1) -
             A(x, y, z + 2)) *
            h / 12.f;
  float *F = f + B * 3 * fs;
  F[u] = a;
  F[fs + u] = b;
  F[2 * fs + u] = c;
  if (x >= 1 && x <= np && y >= 1 && y <= np && z >= 1 && z <= np)
    mx(f2 + B, a * a + b * b + c * c);
}
__global__ void interp(const CellDesc *p, int N, const int16_t *xp,
                       const float *f, float *dv, int np, float kick, size_t fs) {
  int q = blockIdx.x;
  if (q >= N || !p[q].active)
    return;
  const CellDesc d = p[q];
    int nf = np + 2;
    const float *F = f + d.batch * 3 * fs;
    for (uint32_t l = threadIdx.x; l < d.count; l += blockDim.x) {
    uint32_t ip0 = d.base + l;
    int x, y, z;
    float a[3], b[3], c[3], X = 0, Y = 0, Z = 0;
    float px = ((int16_t)(xp[3 * ip0] - (int16_t)0x8000) + 32768.5f) / 65536.f,
          py = ((int16_t)(xp[3 * ip0 + 1] - (int16_t)0x8000) + 32768.5f) /
               65536.f,
          pz = ((int16_t)(xp[3 * ip0 + 2] - (int16_t)0x8000) + 32768.5f) /
               65536.f;
    tsc(((d.cell[0] - ncb) + px) * (rcs * (np / (ng / (nnt * nns)))) -
            d.utile[0] * np,
        x, a);
    tsc(((d.cell[1] - ncb) + py) * (rcs * (np / (ng / (nnt * nns)))) -
            d.utile[1] * np,
        y, b);
    tsc(((d.cell[2] - ncb) + pz) * (rcs * (np / (ng / (nnt * nns)))) -
            d.utile[2] * np,
        z, c);
    for (int k = -1; k < 2; k++)
      for (int j = -1; j < 2; j++)
        for (int i = -1; i < 2; i++) {
          int xx = x + i, yy = y + j, zz = z + k;
              if ((unsigned)xx < (unsigned)nf && (unsigned)yy < (unsigned)nf &&
                  (unsigned)zz < (unsigned)nf) {
            size_t u = ((size_t)zz * nf + yy) * nf + xx;
            float v = a[i + 1] * b[j + 1] * c[k + 1];
            X += F[u] * v;
            Y += F[fs + u] * v;
            Z += F[2 * fs + u] * v;
          }
        }
    size_t ip = (size_t)ip0 * 3;
    dv[ip] += X * kick;
    dv[ip + 1] += Y * kick;
    dv[ip + 2] += Z * kick;
  }
}
int fft_batch_bucket(int jobs) {
  // Finer-than-power-of-two buckets retain most of the plan reuse while
  // substantially reducing padded FFT work for partially occupied parents.
  constexpr int buckets[] = {1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64};
  for (int bucket : buckets)
    if (jobs <= bucket)
      return bucket;
  std::fprintf(stderr, "PM3 batch size %d exceeds bucket capacity\n", jobs);
  std::abort();
}
void run_batch(const std::vector<Job> &jobs, int *rh, int64_t *ix, int16_t *xp,
               const float *g, int rcp, float mass, float amid,
               float dt, float *d_dv, bool *velocity_waited, float *fo,
               double *pt, Profile *profile) {
  const int job_count = jobs.size();
  const int nb = fft_batch_bucket(job_count);
  int np = ng / (nnt * nns) * rcp, fb = rcs * rcp, n = np + 2 * fb;
  int(*R)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb][nt + 2 * ncb] =
      (int(*)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb][nt + 2 * ncb]) rh;
  int64_t(*I)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb] =
      (int64_t(*)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb]) ix;
  std::vector<CellDesc> &h = w.host_cells;
  h.clear();
  const auto pack_start = std::chrono::steady_clock::now();
  for (int b = 0; b < job_count; b++) {
    const Job &j = jobs[b];
    int tx = j.tile[0] - 1, ty = j.tile[1] - 1, tz = j.tile[2] - 1, a[3], e[3];
    for (int d = 0; d < 3; d++) {
      a[d] = j.nc1[d] + ncb - 1;
      e[d] = j.nc2[d] + ncb - 1;
    }
    const int q_begin = a[0] - 1, q_end = e[0] + 1;
    for (int k = a[2] - 1; k <= e[2] + 1; k++)
      for (int y = a[1] - 1; y <= e[1] + 1; y++) {
        // I points just after the row.  Compute its suffix once, then update
        // it while q advances through contiguous memory.
        int sum = 0;
        for (int q = q_begin; q < nt + 2 * ncb; q++)
          sum += R[tz][ty][tx][k][y][q];
        const bool active_yz = y >= a[1] && y <= e[1] && k >= a[2] && k <= e[2];
        for (int q = q_begin; q <= q_end; q++) {
          const int cnt = R[tz][ty][tx][k][y][q];
          const int64_t base = I[tz][ty][tx][k][y] - sum;
          const bool ac = active_yz && q >= a[0] && q <= e[0];
          if (ac)
            pt[j.it] += cnt;
          if (cnt)
            h.push_back({(uint32_t)base,
                         (uint32_t)cnt,
                         {(int16_t)q, (int16_t)y, (int16_t)k},
                         {(int16_t)j.utile[0], (int16_t)j.utile[1],
                          (int16_t)j.utile[2]},
                         (uint16_t)b,
                         (uint8_t)ac});
          sum -= cnt;
        }
      }
  }
  if (profile && profile->on) {
    profile->pack_ms += std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - pack_start)
                            .count();
    profile->batches++;
    profile->fft_jobs += job_count;
    profile->fft_slots += nb;
    for (const auto &d : h)
      profile->particles += d.count;
  }
  size_t rs = (size_t)(n + 2) * n * n, nk = (size_t)(n / 2 + 1) * n * n,
         fs = (size_t)(np + 2) * (np + 2) * (np + 2);
  const auto ensure_start = std::chrono::steady_clock::now();
  FFTPlan *plan = ensure(n, np, nb, h.size(), profile);
  Green *green_field = green_for(g, nk, profile);
  if (profile && profile->on)
    profile->ensure_ms += std::chrono::duration<double, std::milli>(
                              std::chrono::steady_clock::now() - ensure_start)
                              .count();
  gpu_time(profile, profile->h2d_ms, [&] {
    if (!h.empty())
      CK(cudaMemcpy(w.p, h.data(), h.size() * sizeof(CellDesc),
                    cudaMemcpyHostToDevice));
  });
  gpu_time(profile, profile->clear_ms, [&] {
    CK(cudaMemset(w.rho, 0, nb * rs * sizeof(float)));
    CK(cudaMemset(w.f2, 0, nb * sizeof(float)));
  });
  int bl = h.size();
  gpu_time(profile, profile->deposit_ms, [&] {
    if (bl)
      dep<<<bl, T>>>(w.p, h.size(), d_xp, w.rho, n, np, rcp, fb, mass, rs);
  });
  gpu_time(profile, profile->fft_forward_ms,
           [&] { FK(cufftExecR2C(plan->fw, w.rho, (cufftComplex *)w.rho)); });
  gpu_time(profile, profile->green_ms, [&] {
    green<<<(nb * nk + T - 1) / T, T>>>((cufftComplex *)w.rho,
                                        green_field->device, nk, nb * nk);
  });
  gpu_time(profile, profile->fft_backward_ms,
           [&] { FK(cufftExecC2R(plan->bw, (cufftComplex *)w.rho, w.rho)); });
  gpu_time(profile, profile->gradient_ms, [&] {
    grad<<<(nb * fs + T - 1) / T, T>>>(w.rho, w.force, w.f2, n, np, fb, nb, rs,
                                       fs);
  });
  if (bl && !*velocity_waited) {
    c_gpu_velocity_wait_default_stream();
    *velocity_waited = true;
  }
  gpu_time(profile, profile->interpolate_ms, [&] {
    if (bl)
      interp<<<bl, T>>>(w.p, h.size(), d_xp, w.force, d_dv,
                        np,
                        amid * dt / (6 * 3.141592653589793f), fs);
  });
  CK(cudaGetLastError());
  std::vector<float> &hf = w.host_f2;
  hf.resize(nb);
  gpu_time(profile, profile->reduce_ms, [&] {
    CK(cudaMemcpy(hf.data(), w.f2, nb * sizeof(float), cudaMemcpyDeviceToHost));
  });
  for (int b = 0; b < nb; b++) {
    *fo = fmaxf(*fo, hf[b]);
  }
}
} // namespace
extern "C" void c_pm3_force_kernel(
    const int *isort, const int *ires, const int (*xyz)[6], const int *ratio,
    int *rh, int64_t *ix, int16_t *xp, const float *g2, const float *g4,
    const float *g6, const float *g8, const float *g12, float * /*vp*/, float mass,
    float amid, float dt, float *f2, float vm[3], double *pt) {
  std::lock_guard<std::mutex> l(m);
  Profile profile;
  profile.on = std::getenv("CUBE_PM3_PROFILE") != nullptr;
  const auto gpu_select_start = std::chrono::steady_clock::now();
  c_select_gpu_for_image();
  if (profile.on) {
    profile.gpu_select_ms = std::chrono::duration<double, std::milli>(
                                std::chrono::steady_clock::now() - gpu_select_start)
                                .count();
    CK(cudaEventCreate(&profile.start));
    CK(cudaEventCreate(&profile.stop));
  }
  int64_t(*spine)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb] =
      (int64_t(*)[nnt][nnt][nt + 2 * ncb][nt + 2 * ncb]) ix;
  const int last = nt + 2 * ncb - 1;
  const int64_t n_used = spine[nnt - 1][nnt - 1][nnt - 1][last][last];
  if (n_used < 0) {
    std::fprintf(stderr, "invalid PM3 particle count\n");
    std::abort();
  }
  sync_xp(xp, static_cast<size_t>(n_used), &profile);
  const auto velocity_init_start = std::chrono::steady_clock::now();
  float *d_dv = c_gpu_velocity_acquire(ratio[5], static_cast<int>(n_used));
  if (profile.on)
    profile.velocity_init_ms += std::chrono::duration<double, std::milli>(
                                       std::chrono::steady_clock::now() -
                                       velocity_init_start)
                                       .count();
  bool velocity_waited = false;
  const float *g[7] = {0, 0, g2, g4, g6, g8, g12};
  *f2 = 0;
  vm[0] = vm[1] = vm[2] = 0;
  constexpr int pm3_levels = 5;
  const int group_count = nnt * nnt * nnt * pm3_levels;
  std::vector<std::vector<Job>> &job_groups = w.host_job_groups;
  if ((int)job_groups.size() != group_count)
    job_groups.resize(group_count);
  for (std::vector<Job> &jobs : job_groups)
    jobs.clear();
  const auto schedule_start = std::chrono::steady_clock::now();
  for (int o = 0; o < nnt * nnt * nnt * nns * nns * nns; o++) {
    const int it = isort[o] - 1, lev = ires[it];
    if (lev < 2 || lev > 6)
      continue;
    const int px = xyz[it][3] - 1, py = xyz[it][4] - 1, pz = xyz[it][5] - 1;
    if ((unsigned)px >= (unsigned)nnt || (unsigned)py >= (unsigned)nnt ||
        (unsigned)pz >= (unsigned)nnt)
      continue;
    const int pa = px + nnt * (py + nnt * pz);
    Job z{};
    z.it = it;
    z.level = lev;
    for (int d = 0; d < 3; d++) {
      z.tile[d] = xyz[it][d + 3];
      z.utile[d] = xyz[it][d] - 1;
      z.nc1[d] = z.utile[d] * ntt + 1;
      z.nc2[d] = (z.utile[d] + 1) * ntt;
    }
    job_groups[pa * pm3_levels + lev - 2].push_back(z);
  }
  if (profile.on)
    profile.schedule_ms = std::chrono::duration<double, std::milli>(
                              std::chrono::steady_clock::now() - schedule_start)
                              .count();
  for (int pa = 0; pa < nnt * nnt * nnt; pa++)
    for (int lev = 2; lev <= 6; lev++) {
      const std::vector<Job> &jobs = job_groups[pa * pm3_levels + lev - 2];
      if (!jobs.empty())
        run_batch(jobs, rh, ix, xp, g[lev], ratio[lev - 1], mass, amid, dt,
                  d_dv, &velocity_waited, f2, pt, &profile);
    }
  if (profile.on) {
    std::fprintf(
        stderr,
        "[PM3 profile] batches=%d particles=%d fft_jobs/slots=%d/%d schedule=%.2fms "
        "pack=%.2fms ensure=%.2fms clear=%.2fms xp_h2d=%.2fms "
        "h2d=%.2f dep=%.2f fft_f=%.2f green=%.2f fft_b=%.2f "
        "grad=%.2f interp=%.2f reduce=%.2f\n"
        "[PM3 ensure] plan hit/miss/evict=%d/%d/%d create=%.2fms "
        "evict=%.2fms; reserve grows=%d bytes=%.2fMiB time=%.2fms; "
        "green hit/miss=%d/%d upload=%.2fms; cache plans=%zu greens=%zu "
        "fft_work=%.2fMiB\n"
        "[PM3 setup] gpu_select=%.2fms velocity_delta_init=%.2fms\n",
        profile.batches, profile.particles, profile.fft_jobs, profile.fft_slots,
        profile.schedule_ms, profile.pack_ms,
        profile.ensure_ms, profile.clear_ms, profile.xp_h2d_ms, profile.h2d_ms,
        profile.deposit_ms, profile.fft_forward_ms,
        profile.green_ms, profile.fft_backward_ms, profile.gradient_ms,
        profile.interpolate_ms, profile.reduce_ms, profile.plan_hits,
        profile.plan_misses, profile.plan_evictions, profile.plan_create_ms,
        profile.plan_evict_ms, profile.reserve_grows,
        profile.reserve_bytes / (1024.0 * 1024.0), profile.reserve_ms,
        profile.green_hits, profile.green_misses, profile.green_upload_ms,
        w.plans.size(), w.greens.size(),
        w.fft_work_cap / (1024.0 * 1024.0), profile.gpu_select_ms,
        profile.velocity_init_ms);
    CK(cudaEventDestroy(profile.start));
    CK(cudaEventDestroy(profile.stop));
  }
}
extern "C" void c_pm3_fft_finalize() {
  std::lock_guard<std::mutex> l(m);
  drop();
  if (d_xp)
    CK(cudaFree(d_xp));
  d_xp = nullptr;
  d_xp_capacity = 0;
}
