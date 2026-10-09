#include "pipe_sim/banded_solver.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

namespace pipe_sim {
extern "C" cudaError_t launch_pipe_banded_solve(
    int,double*,const double*,double*,int*);
}

#define CHECK(call) do {                                      \
    cudaError_t err = (call);                                \
    if (err != cudaSuccess) {                                 \
        std::fprintf(stderr,"%s:%d: %s\n",                    \
            __FILE__,__LINE__,cudaGetErrorString(err));       \
        return 1;                                             \
    }                                                         \
} while (0)

int main()
{
    int comparisons = 0;
    int exact = 0;
    int nonfinite = 0;
    int info_mismatches = 0;

    double max_abs = 0.0;
    double max_rel = 0.0;

    int worst_n = 0;
    int worst_scenario = 0;
    int worst_index = 0;

    const char* worst_kind = "none";

    double worst_cpu = 0.0;
    double worst_gpu = 0.0;

    for (int n : {2,3,10,25,100,202}) {
        constexpr int ldab = 6;
        constexpr int kl = 2;
        constexpr int ku = 1;
        constexpr int kv = 3;

        std::vector<double> ab(ldab*n,0.0);
        std::vector<double> b(n);

        const auto set = [&](int i,int j,double value) {
            const int row = kv+i-j;

            if (row >= 0 && row < ldab)
                ab[j*ldab+row] = value;
        };

        double* dab = nullptr;
        double* db = nullptr;
        double* dx = nullptr;
        int* dinfo = nullptr;

        CHECK(cudaMalloc(&dab,ab.size()*sizeof(double)));
        CHECK(cudaMalloc(&db,b.size()*sizeof(double)));
        CHECK(cudaMalloc(&dx,b.size()*sizeof(double)));
        CHECK(cudaMalloc(&dinfo,sizeof(int)));

        for (int scenario = 0; scenario < 4; ++scenario) {
            std::fill(ab.begin(),ab.end(),0.0);

            for (int i = 0; i < n; ++i) {
                set(i,i,5.0+0.01*i);

                if (i+1 < n) {
                    set(i+1,i,-0.5);
                    set(i,i+1,0.25);
                }

                if (i+2 < n)
                    set(i+2,i,0.125);

                b[i] = 1.0+0.1*i;
            }

            if (scenario == 1) {
                set(0,0,0.01);
                set(1,0,2.0);
            }

            if (scenario == 2) {
                set(0,0,0.0);
                set(1,0,0.0);

                if (n > 2)
                    set(2,0,0.0);
            }

            if (scenario == 3) {
                set(n-1,n-1,0.0);

                if (n >= 2)
                    set(n-1,n-2,0.0);

                if (n >= 3)
                    set(n-1,n-3,0.0);
            }

            std::vector<double> cpu_ab = ab;
            std::vector<double> cpu_x(n,0.0);

            std::vector<double> gpu_ab(ab.size());
            std::vector<double> gpu_x(n);

            const int cpu_info = pipe_sim::solve_banded(
                n,cpu_ab,ldab,kl,ku,b,cpu_x);

            CHECK(cudaMemcpy(
                dab,ab.data(),ab.size()*sizeof(double),
                cudaMemcpyHostToDevice));

            CHECK(cudaMemcpy(
                db,b.data(),b.size()*sizeof(double),
                cudaMemcpyHostToDevice));

            CHECK(pipe_sim::launch_pipe_banded_solve(
                n,dab,db,dx,dinfo));

            CHECK(cudaDeviceSynchronize());

            int gpu_info = -999;

            CHECK(cudaMemcpy(
                &gpu_info,dinfo,sizeof(int),
                cudaMemcpyDeviceToHost));

            CHECK(cudaMemcpy(
                gpu_ab.data(),dab,
                ab.size()*sizeof(double),
                cudaMemcpyDeviceToHost));

            CHECK(cudaMemcpy(
                gpu_x.data(),dx,
                b.size()*sizeof(double),
                cudaMemcpyDeviceToHost));

            ++comparisons;

            if (cpu_info == gpu_info) {
                ++exact;
            } else {
                ++info_mismatches;

                std::printf(
                    "INFO FAIL n=%d scenario=%d cpu=%d gpu=%d\n",
                    n,scenario,cpu_info,gpu_info);
            }

            auto compare = [&](const std::vector<double>& cpu,
                               const std::vector<double>& gpu,
                               const char* kind) {
                for (std::size_t j = 0; j < cpu.size(); ++j) {
                    ++comparisons;

                    if (!std::isfinite(cpu[j]) ||
                        !std::isfinite(gpu[j])) {
                        ++nonfinite;
                        continue;
                    }

                    if (std::memcmp(
                        &cpu[j],&gpu[j],sizeof(double)) == 0)
                        ++exact;

                    const double error =
                        std::abs(cpu[j]-gpu[j]);

                    const double scale =
                        std::max(std::abs(cpu[j]),
                                 std::abs(gpu[j]));

                    const double relative =
                        scale > 0.0 ? error/scale : 0.0;

                    max_abs = std::max(max_abs,error);

                    if (relative > max_rel) {
                        max_rel = relative;
                        worst_n = n;
                        worst_scenario = scenario;
                        worst_index = static_cast<int>(j);
                        worst_kind = kind;
                        worst_cpu = cpu[j];
                        worst_gpu = gpu[j];
                    }
                }
            };

            compare(cpu_ab,gpu_ab,"matrix");
            compare(cpu_x,gpu_x,"rhs/solution");
        }

        CHECK(cudaFree(dab));
        CHECK(cudaFree(db));
        CHECK(cudaFree(dx));
        CHECK(cudaFree(dinfo));
    }

    std::printf(
        "Comparisons: %d\n"
        "Bitwise identical: %d\n"
        "Nonfinite: %d\n"
        "INFO mismatches: %d\n",
        comparisons,exact,nonfinite,info_mismatches);

    std::printf(
        "Maximum absolute error: %.17e\n"
        "Maximum relative error: %.17e\n",
        max_abs,max_rel);

    std::printf(
        "Worst relative discrepancy: "
        "n=%d scenario=%d %s index=%d\n"
        "CPU=%.17e GPU=%.17e\n",
        worst_n,worst_scenario,worst_kind,worst_index,
        worst_cpu,worst_gpu);

    // Diagnostic phase: inspect numerical discrepancies separately.
    if (nonfinite || info_mismatches)
        return 1;

    return 0;
}
