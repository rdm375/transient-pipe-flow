#include "pipe_sim/transient.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

namespace pipe_sim {
extern "C" cudaError_t launch_pipe_residual(
    int, const double*, const double*, const double*,
    TransientBoundary, TransientParameters, double*);
}

#define CHECK(call) do {                                      \
    cudaError_t err = (call);                                 \
    if (err != cudaSuccess) {                                 \
        std::fprintf(stderr, "%s:%d: %s\n",                   \
                     __FILE__, __LINE__,                      \
                     cudaGetErrorString(err));                \
        return 1;                                             \
    }                                                         \
} while (0)

int main()
{
    using namespace pipe_sim;

    int comparisons = 0;
    int exact = 0;
    int nonfinite = 0;

    double max_absolute = 0.0;
    double max_relative = 0.0;

    int worst_n = 0;
    int worst_scenario = 0;
    int worst_component = 0;

    double worst_cpu = 0.0;
    double worst_gpu = 0.0;

    for (int n : {2,3,10,25,100}) {
        const double pi = std::acos(-1.0);

        const TransientParameters par{
            100000.0/n, 15.0, 0.65,
            1.0, pi/4.0, 288.15, 0.9,
            500.0, 1.1e-5, 4.5e-5
        };

        const TransientBoundary bc{
            100.0, 100.0, 105.0, 8.0e6
        };

        std::vector<double> po(n+1);
        std::vector<double> mo(n);
        std::vector<double> u(2*n+2);
        std::vector<double> cpu(2*n+2);
        std::vector<double> gpu(2*n+2);

        double *du=nullptr, *dpo=nullptr;
        double *dmo=nullptr, *dr=nullptr;

        CHECK(cudaMalloc(&du,u.size()*sizeof(double)));
        CHECK(cudaMalloc(&dpo,po.size()*sizeof(double)));
        CHECK(cudaMalloc(&dmo,mo.size()*sizeof(double)));
        CHECK(cudaMalloc(&dr,gpu.size()*sizeof(double)));

        for (int scenario=0; scenario<4; ++scenario) {
            for (int i=0; i<=n; ++i)
                po[i] = 8.0e6-1.0e5*i/n;

            for (int i=0; i<n; ++i) {
                const double flow =
                    (scenario==2) ? -100.0 : 100.0;
                mo[i] = flow+0.1*i/n;
            }

            u[0] = bc.inlet_pressure_new+100.0*scenario;
            u[1] = (scenario==2) ? -100.5 : 100.5;

            for (int i=0; i<n; ++i) {
                u[2*i+2] = mo[i]+0.25;
                u[2*i+3] = po[i+1]+200.0*(scenario+1);
            }

            assemble_residual(n,u,po,mo,bc,par,cpu);

            CHECK(cudaMemcpy(
                du,u.data(),u.size()*sizeof(double),
                cudaMemcpyHostToDevice));

            CHECK(cudaMemcpy(
                dpo,po.data(),po.size()*sizeof(double),
                cudaMemcpyHostToDevice));

            CHECK(cudaMemcpy(
                dmo,mo.data(),mo.size()*sizeof(double),
                cudaMemcpyHostToDevice));

            CHECK(launch_pipe_residual(
                n,du,dpo,dmo,bc,par,dr));

            CHECK(cudaDeviceSynchronize());

            CHECK(cudaMemcpy(
                gpu.data(),dr,gpu.size()*sizeof(double),
                cudaMemcpyDeviceToHost));

            for (std::size_t j=0; j<cpu.size(); ++j) {
                ++comparisons;

                if (!std::isfinite(cpu[j]) ||
                    !std::isfinite(gpu[j])) {
                    ++nonfinite;
                    continue;
                }

                if (std::memcmp(
                        &cpu[j],&gpu[j],sizeof(double))==0)
                    ++exact;

                const double error =
                    std::abs(cpu[j]-gpu[j]);

                const double scale =
                    std::max(std::abs(cpu[j]),
                             std::abs(gpu[j]));

                const double relative =
                    scale>0.0 ? error/scale : 0.0;

                if (error>max_absolute)
                    max_absolute=error;

                if (relative>max_relative) {
                    max_relative=relative;
                    worst_n=n;
                    worst_scenario=scenario;
                    worst_component=static_cast<int>(j);
                    worst_cpu=cpu[j];
                    worst_gpu=gpu[j];
                }
            }
        }

        CHECK(cudaFree(du));
        CHECK(cudaFree(dpo));
        CHECK(cudaFree(dmo));
        CHECK(cudaFree(dr));
    }

    std::printf("Comparisons: %d\n",comparisons);
    std::printf("Bitwise identical: %d\n",exact);
    std::printf("Nonfinite: %d\n",nonfinite);
    std::printf("Maximum absolute error: %.17e\n",
                max_absolute);
    std::printf("Maximum relative error: %.17e\n",
                max_relative);

    std::printf(
        "Worst relative discrepancy: "
        "n=%d scenario=%d component=%d\n"
        "CPU=%.17e GPU=%.17e\n",
        worst_n,worst_scenario,worst_component,
        worst_cpu,worst_gpu);

    if (nonfinite != 0)
        return 1;

    // Diagnostic experiment: numerical acceptance follows
    // investigation of any observed discrepancies.
    return 0;
}
