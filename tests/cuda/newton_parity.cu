#include "pipe_sim/newton.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

namespace pipe_sim {
extern "C" cudaError_t launch_pipe_newton(
    int, double*, const double*, const double*,
    TransientBoundary, TransientParameters,
    NewtonOptions, NewtonResult*);
}

#define CHECK(call) do {                                      \
    cudaError_t err=(call);                                   \
    if (err!=cudaSuccess) {                                   \
        std::fprintf(stderr,"%s:%d: %s\n",                    \
            __FILE__,__LINE__,cudaGetErrorString(err));       \
        return 1;                                             \
    }                                                         \
} while (0)

int main()
{
    using namespace pipe_sim;

    int comparisons=0;
    int exact=0;
    int nonfinite=0;
    int info_mismatches=0;
    int iteration_mismatches=0;
    int state_failures=0;
    int failure_state_changes=0;

    double max_abs=0.0;
    double max_rel=0.0;

    int worst_n=0;
    int worst_scenario=0;
    int worst_index=0;
    double worst_cpu=0.0;
    double worst_gpu=0.0;

    for (int n : {2,3,10,25,100}) {
        const double pi=std::acos(-1.0);

        const TransientParameters par{
            100000.0/n,15.0,0.65,1.0,pi/4.0,
            288.15,0.9,500.0,1.1e-5,4.5e-5
        };

        const int nu=2*n+2;

        double* dstate=nullptr;
        double* dpo=nullptr;
        double* dmo=nullptr;
        NewtonResult* dresult=nullptr;

        CHECK(cudaMalloc(&dstate,nu*sizeof(double)));
        CHECK(cudaMalloc(&dpo,(n+1)*sizeof(double)));
        CHECK(cudaMalloc(&dmo,n*sizeof(double)));
        CHECK(cudaMalloc(&dresult,sizeof(NewtonResult)));

        for (int scenario=0; scenario<6; ++scenario) {
            TransientBoundary bc{
                100.0,100.0,105.0,8.0e6
            };

            NewtonOptions opt{
                1.0e-9,1.0e-12,30,false
            };

            std::vector<double> po(n+1);
            std::vector<double> mo(n,100.0);
            std::vector<double> initial(nu);

            const double re=
                100.0*par.diameter/(par.area*par.viscosity);

            const double arg=
                par.roughness/(3.7*par.diameter)
                +5.74/std::pow(re,0.9);

            const double fd=
                0.25/std::pow(std::log10(arg),2.0);

            const double coef=
                fd*par.z*par.gas_constant*par.temperature
                *100.0*100.0
                /(par.diameter*par.area*par.area);

            for (int i=0; i<=n; ++i) {
                po[i]=std::sqrt(
                    bc.inlet_pressure_new*bc.inlet_pressure_new
                    -coef*i*par.dx);
            }

            initial[0]=po[0];
            initial[1]=100.0;

            for (int i=0; i<n; ++i) {
                initial[2*i+2]=100.0;
                initial[2*i+3]=po[i+1];
            }

            if (scenario==1) {
                initial[0]+=100.0;
                initial[1]+=0.5;

                for (int i=0; i<n; ++i) {
                    initial[2*i+2]+=0.25;
                    initial[2*i+3]+=100.0;
                }
            }

            if (scenario==2)
                opt.max_iterations=0;

            if (scenario==3)
                bc.inlet_pressure_new=-1.0;

            if (scenario==4)
                opt.step_tolerance=0.0;

            if (scenario==5)
                po[n/2]=-1.0;

            std::vector<double> cpu=initial;
            std::vector<double> gpu(nu);

            NewtonWorkspace ws(n);
            const NewtonResult cpu_result=
                newton_solve_banded(
                    n,cpu,po,mo,bc,par,opt,ws);

            CHECK(cudaMemcpy(
                dstate,initial.data(),nu*sizeof(double),
                cudaMemcpyHostToDevice));

            CHECK(cudaMemcpy(
                dpo,po.data(),(n+1)*sizeof(double),
                cudaMemcpyHostToDevice));

            CHECK(cudaMemcpy(
                dmo,mo.data(),n*sizeof(double),
                cudaMemcpyHostToDevice));

            CHECK(launch_pipe_newton(
                n,dstate,dpo,dmo,bc,par,opt,dresult));

            CHECK(cudaDeviceSynchronize());

            NewtonResult gpu_result{};

            CHECK(cudaMemcpy(
                &gpu_result,dresult,sizeof(NewtonResult),
                cudaMemcpyDeviceToHost));

            CHECK(cudaMemcpy(
                gpu.data(),dstate,nu*sizeof(double),
                cudaMemcpyDeviceToHost));

            ++comparisons;
            if (cpu_result.info==gpu_result.info) {
                ++exact;
            } else {
                ++info_mismatches;
                std::printf(
                    "INFO FAIL n=%d scenario=%d cpu=%d gpu=%d\n",
                    n,scenario,cpu_result.info,gpu_result.info);
            }

            ++comparisons;
            if (cpu_result.iterations==gpu_result.iterations) {
                ++exact;
            } else {
                ++iteration_mismatches;
                std::printf(
                    "NITER FAIL n=%d scenario=%d cpu=%d gpu=%d\n",
                    n,scenario,
                    cpu_result.iterations,gpu_result.iterations);
            }

            if (cpu_result.info!=0) {
                if (std::memcmp(
                    gpu.data(),initial.data(),
                    nu*sizeof(double))!=0) {
                    ++failure_state_changes;
                    std::printf(
                        "FAILURE STATE CHANGED n=%d scenario=%d\n",
                        n,scenario);
                }
            }

            for (int i=0; i<nu; ++i) {
                ++comparisons;

                if (!std::isfinite(cpu[i]) ||
                    !std::isfinite(gpu[i])) {
                    ++nonfinite;
                    continue;
                }

                if (std::memcmp(
                    &cpu[i],&gpu[i],sizeof(double))==0)
                    ++exact;

                const double error=std::abs(cpu[i]-gpu[i]);
                const double scale=std::max(
                    std::abs(cpu[i]),std::abs(gpu[i]));

                const double relative=
                    scale>0.0 ? error/scale : 0.0;

                max_abs=std::max(max_abs,error);

                if (relative>max_rel) {
                    max_rel=relative;
                    worst_n=n;
                    worst_scenario=scenario;
                    worst_index=i;
                    worst_cpu=cpu[i];
                    worst_gpu=gpu[i];
                }

                if (error>5.0e-14*scale)
                    ++state_failures;
            }

            std::printf(
                "n=%3d scenario=%d CPU=(%d,%d) GPU=(%d,%d)\n",
                n,scenario,
                cpu_result.info,cpu_result.iterations,
                gpu_result.info,gpu_result.iterations);
        }

        CHECK(cudaFree(dstate));
        CHECK(cudaFree(dpo));
        CHECK(cudaFree(dmo));
        CHECK(cudaFree(dresult));
    }

    std::printf(
        "\nComparisons: %d\n"
        "Bitwise identical: %d\n"
        "Nonfinite: %d\n"
        "INFO mismatches: %d\n"
        "Iteration mismatches: %d\n"
        "State tolerance failures: %d\n"
        "Failure-state changes: %d\n",
        comparisons,exact,nonfinite,info_mismatches,
        iteration_mismatches,state_failures,
        failure_state_changes);

    std::printf(
        "Maximum absolute error: %.17e\n"
        "Maximum relative error: %.17e\n",
        max_abs,max_rel);

    std::printf(
        "Worst relative discrepancy: "
        "n=%d scenario=%d index=%d\n"
        "CPU=%.17e GPU=%.17e\n",
        worst_n,worst_scenario,worst_index,
        worst_cpu,worst_gpu);

    if (nonfinite || info_mismatches ||
        iteration_mismatches || state_failures ||
        failure_state_changes) {
        std::puts("FAIL M14.4a Newton parity");
        return 1;
    }

    std::puts("PASS M14.4a Newton parity");
    return 0;
}
