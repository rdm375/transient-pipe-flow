#include "pipe_sim/transient.hpp"
#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

namespace pipe_sim {
extern "C" cudaError_t launch_pipe_jacobian(int, const double*,
                                             TransientParameters, double*);
}
#define CHECK(call) do { auto err=(call); if(err!=cudaSuccess) { \
  std::fprintf(stderr,"%s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(err)); \
  return 1; } } while(0)

int main() {
    using namespace pipe_sim;
    constexpr int ldab=6;
    int comparisons=0, exact=0, nonfinite=0;
    double max_absolute=0.0, max_relative=0.0;
    int worst_n=0, worst_scenario=0, worst_index=0;
    double worst_cpu=0.0, worst_gpu=0.0;
    for (int n : {2,3,10,25,100}) {
        const TransientParameters par{
            100000.0/n,15.0,0.65,1.0,std::acos(-1.0)/4.0,
            288.15,0.9,500.0,1.1e-5,4.5e-5
        };
        const int nu=2*n+2;
        std::vector<double> u(nu), cpu(ldab*nu), gpu(ldab*nu);
        double *du=nullptr, *dab=nullptr;
        CHECK(cudaMalloc(&du,u.size()*sizeof(double)));
        CHECK(cudaMalloc(&dab,gpu.size()*sizeof(double)));
        for (int scenario=0; scenario<4; ++scenario) {
            u[0]=8.0e6+100.0*scenario;
            u[1]=(scenario==2) ? -100.0 : 100.0;
            for (int i=0; i<n; ++i) {
                u[2*i+2]=((scenario==2) ? -100.0 : 100.0)+0.25*(i+1);
                u[2*i+3]=8.0e6-1.0e5*(i+1)/n+200.0*scenario;
            }
            assemble_jacobian_banded(n,u,par,cpu,ldab);
            CHECK(cudaMemcpy(du,u.data(),u.size()*sizeof(double),cudaMemcpyHostToDevice));
            CHECK(launch_pipe_jacobian(n,du,par,dab));
            CHECK(cudaDeviceSynchronize());
            CHECK(cudaMemcpy(gpu.data(),dab,gpu.size()*sizeof(double),cudaMemcpyDeviceToHost));
            for (std::size_t j=0; j<cpu.size(); ++j) {
                ++comparisons;
                if (!std::isfinite(cpu[j]) || !std::isfinite(gpu[j])) {
                    ++nonfinite;
                    continue;
                }
                if (std::memcmp(&cpu[j],&gpu[j],sizeof(double))==0) ++exact;
                const double error=std::abs(cpu[j]-gpu[j]);
                const double scale=std::max(std::abs(cpu[j]),std::abs(gpu[j]));
                const double relative=(scale>0.0) ? error/scale : 0.0;
                max_absolute=std::max(max_absolute,error);
                if (relative>max_relative) {
                    max_relative=relative;
                    worst_n=n; worst_scenario=scenario;
                    worst_index=static_cast<int>(j);
                    worst_cpu=cpu[j]; worst_gpu=gpu[j];
                }
            }
        }
        CHECK(cudaFree(du));
        CHECK(cudaFree(dab));
    }
    std::printf("Comparisons: %d\nBitwise identical: %d\nNonfinite: %d\n",comparisons,exact,nonfinite);
    std::printf("Maximum absolute error: %.17e\n",max_absolute);
    std::printf("Maximum relative error: %.17e\n",max_relative);
    std::printf("Worst relative discrepancy: n=%d scenario=%d storage_index=%d\n",worst_n,worst_scenario,worst_index);
    std::printf("CPU=%.17e GPU=%.17e\n",worst_cpu,worst_gpu);
    return nonfinite ? 1 : 0; // Diagnostic; acceptance requires discrepancy review.
}
