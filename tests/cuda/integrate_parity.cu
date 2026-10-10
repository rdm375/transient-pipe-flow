#include "pipe_sim/integrate.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

namespace pipe_sim {
extern "C" cudaError_t launch_pipe_integrate(
    int, double*, double*, double*,
    DemandRamp, int, TransientParameters, NewtonOptions,
    double*, double*, double*, double*, double*,
    int*, IntegrationResult*);
}

#define CHECK(call) do {                                      \
    const auto err = (call);                                 \
    if (err != cudaSuccess) {                                 \
        std::fprintf(stderr, "%s:%d: %s\n",                   \
            __FILE__, __LINE__, cudaGetErrorString(err));     \
        return 1;                                             \
    }                                                         \
} while (0)

template<class T>
struct DeviceBuffer {
    T* ptr = nullptr;
    std::size_t count;

    explicit DeviceBuffer(std::size_t n) : count(n) {
        cudaError_t err = cudaMalloc(
            reinterpret_cast<void**>(&ptr), n*sizeof(T));
        if (err != cudaSuccess) {
            std::fprintf(stderr, "cudaMalloc: %s\n",
                         cudaGetErrorString(err));
            std::abort();
        }
    }

    ~DeviceBuffer() { cudaFree(ptr); }

    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    void upload(const std::vector<T>& v) {
        cudaMemcpy(ptr, v.data(), count*sizeof(T),
                   cudaMemcpyHostToDevice);
    }

    void download(std::vector<T>& v) const {
        cudaMemcpy(v.data(), ptr, count*sizeof(T),
                   cudaMemcpyDeviceToHost);
    }
};

struct Comparison {
    long long count = 0;
    long long exact = 0;
    long long failures = 0;
    double max_relative = 0.0;
    double max_absolute = 0.0;

    char worst_relative_name[128] = {};
    char worst_absolute_name[128] = {};

    int worst_relative_index = -1;
    int worst_absolute_index = -1;

    double worst_relative_cpu = 0.0;
    double worst_relative_gpu = 0.0;

    double worst_absolute_cpu = 0.0;
    double worst_absolute_gpu = 0.0;

    void check(const char* name, int k,
               double cpu, double gpu,
               double relative_tolerance = 5.0e-12,
               double absolute_tolerance = 1.0e-9)
    {
        ++count;

        if (!std::isfinite(cpu) || !std::isfinite(gpu)) {
            ++failures;
            std::printf("FAIL nonfinite %s[%d]\n", name, k);
            return;
        }

        if (std::memcmp(&cpu, &gpu, sizeof(double)) == 0)
            ++exact;

        const double error = std::abs(cpu-gpu);
        const double scale =
            std::max(std::abs(cpu), std::abs(gpu));
        const double relative =
            scale > 0.0 ? error/scale : 0.0;

        if (error > max_absolute) {
            max_absolute = error;
            std::snprintf(
                worst_absolute_name,
                sizeof(worst_absolute_name),
                "%s",name);
            worst_absolute_index = k;
            worst_absolute_cpu = cpu;
            worst_absolute_gpu = gpu;
        }

        if (relative > max_relative) {
            max_relative = relative;
            std::snprintf(
                worst_relative_name,
                sizeof(worst_relative_name),
                "%s",name);
            worst_relative_index = k;
            worst_relative_cpu = cpu;
            worst_relative_gpu = gpu;
        }

        const double tolerance =
            std::max(absolute_tolerance,
                     relative_tolerance*scale);

        if (error > tolerance) {
            ++failures;
            if (failures <= 20)
                std::printf(
                    "FAIL %s[%d] CPU=%.17e GPU=%.17e "
                    "error=%.3e tolerance=%.3e\n",
                    name,k,cpu,gpu,error,tolerance);
        }
    }
};

int run_case(int n, int steps, Comparison& cmp)
{
    using namespace pipe_sim;

    const double pi = std::acos(-1.0);
    const double diameter = 1.0;
    const double area = pi*diameter*diameter/4.0;

    const TransientParameters par{
        100000.0/n, 15.0, 0.65,
        diameter, area, 288.15, 0.9,
        500.0, 1.1e-5, 4.5e-5
    };

    const DemandRamp demand{
        8.0e6, 100.0, 105.0, 600.0
    };

    const NewtonOptions opt{
        1.0e-9, 1.0e-12, 30, false
    };

    std::vector<double> pressure(n+1);
    std::vector<double> flow(n,100.0);

    double inlet_flow = 100.0;

    const double reynolds =
        inlet_flow*diameter/(area*par.viscosity);

    const double friction_arg =
        par.roughness/(3.7*diameter) +
        5.74/std::pow(reynolds,0.9);

    const double friction =
        0.25/std::pow(std::log10(friction_arg),2.0);

    const double coefficient =
        friction*par.z*par.gas_constant*
        par.temperature*inlet_flow*inlet_flow/
        (diameter*area*area);

    for (int i=0; i<=n; ++i)
        pressure[i] = std::sqrt(
            demand.inlet_pressure*demand.inlet_pressure -
            coefficient*i*par.dx);

    auto gpu_pressure = pressure;
    auto gpu_flow = flow;
    double gpu_inlet = inlet_flow;

    const int hsize = steps+1;

    std::vector<double> ct(hsize),ci(hsize),co(hsize);
    std::vector<double> cp(hsize),cl(hsize);
    std::vector<int> cn(hsize);

    IntegrationHistory hist{
        ct,ci,co,cp,cl,cn
    };

    IntegrationWorkspace ws(n);

    const IntegrationResult cpu = integrate_transient(
        n,pressure,flow,inlet_flow,demand,
        steps,par,opt,hist,ws);

    DeviceBuffer<double> dp(n+1),df(n),di(1);
    DeviceBuffer<double> dt(hsize),dhi(hsize);
    DeviceBuffer<double> dho(hsize),dhp(hsize);
    DeviceBuffer<double> dhl(hsize);
    DeviceBuffer<int> dhn(hsize);
    DeviceBuffer<IntegrationResult> dr(1);

    dp.upload(gpu_pressure);
    df.upload(gpu_flow);

    CHECK(cudaMemcpy(di.ptr,&gpu_inlet,sizeof(double),
                     cudaMemcpyHostToDevice));

    cudaEvent_t start_event, stop_event;

    CHECK(cudaEventCreate(&start_event));
    CHECK(cudaEventCreate(&stop_event));

    CHECK(cudaEventRecord(start_event));

    CHECK(launch_pipe_integrate(
        n,dp.ptr,df.ptr,di.ptr,demand,steps,par,opt,
        dt.ptr,dhi.ptr,dho.ptr,dhp.ptr,dhl.ptr,
        dhn.ptr,dr.ptr));

    CHECK(cudaEventRecord(stop_event));
    CHECK(cudaEventSynchronize(stop_event));

    float gpu_elapsed_ms = 0.0f;

    CHECK(cudaEventElapsedTime(
        &gpu_elapsed_ms,start_event,stop_event));

    CHECK(cudaEventDestroy(start_event));
    CHECK(cudaEventDestroy(stop_event));

    if (steps == 8000) {
        std::printf(
            "GPU_BENCHMARK n=%d steps=%d "
            "kernel_ms=%.6f ns_per_step=%.3f\n",
            n,steps,
            static_cast<double>(gpu_elapsed_ms),
            static_cast<double>(gpu_elapsed_ms)*1.0e6/steps);
    }

    IntegrationResult gpu{};

    CHECK(cudaMemcpy(&gpu,dr.ptr,sizeof(gpu),
                     cudaMemcpyDeviceToHost));

    dp.download(gpu_pressure);
    df.download(gpu_flow);

    CHECK(cudaMemcpy(&gpu_inlet,di.ptr,sizeof(double),
                     cudaMemcpyDeviceToHost));

    std::vector<double> gt(hsize),gi(hsize),go(hsize);
    std::vector<double> gp(hsize),gl(hsize);
    std::vector<int> gn(hsize);

    dt.download(gt);
    dhi.download(gi);
    dho.download(go);
    dhp.download(gp);
    dhl.download(gl);
    dhn.download(gn);

    std::printf(
        "n=%3d steps=%5d CPU=(%d,%d) GPU=(%d,%d)\n",
        n,steps,cpu.info,cpu.completed_steps,
        gpu.info,gpu.completed_steps);

    if (cpu.info != gpu.info ||
        cpu.completed_steps != gpu.completed_steps) {
        ++cmp.failures;
        return 1;
    }

    for (int i=0; i<=n; ++i)
        cmp.check("pressure",i,pressure[i],gpu_pressure[i]);

    for (int i=0; i<n; ++i)
        cmp.check("flow",i,flow[i],gpu_flow[i]);

    cmp.check("inlet_flow",0,inlet_flow,gpu_inlet);

    cmp.check("max_mass_defect",0,
              cpu.max_mass_defect,gpu.max_mass_defect,
              0.0,1.0e-9);

    for (int k=0; k<=cpu.completed_steps; ++k) {
        cmp.check("time",k,ct[k],gt[k]);
        cmp.check("history_inlet",k,ci[k],gi[k]);
        cmp.check("history_outlet",k,co[k],go[k]);
        cmp.check("history_pressure",k,cp[k],gp[k]);
        cmp.check("history_linepack",k,cl[k],gl[k]);

        ++cmp.count;
        if (cn[k] == gn[k])
            ++cmp.exact;
        else {
            ++cmp.failures;
            std::printf(
                "FAIL iterations[%d] CPU=%d GPU=%d\n",
                k,cn[k],gn[k]);
        }
    }

    return 0;
}

int main()
{
    Comparison cmp;

    for (int n : {2,3,10,25,100}) {
        for (int steps : {0,1,5,20,80}) {
            if (run_case(n,steps,cmp) != 0)
                return 1;
        }
    }

    // Full M13 reference trajectory.
    if (run_case(100,8000,cmp) != 0)
        return 1;

    std::printf(
        "\nComparisons: %lld\n"
        "Bitwise identical: %lld\n"
        "Failures: %lld\n"
        "Maximum absolute error: %.17e\n"
        "Maximum relative error: %.17e\n",
        cmp.count,cmp.exact,cmp.failures,
        cmp.max_absolute,cmp.max_relative);

    std::printf(
        "Worst absolute: %s[%d] "
        "CPU=%.17e GPU=%.17e\n",
        cmp.worst_absolute_name,
        cmp.worst_absolute_index,
        cmp.worst_absolute_cpu,
        cmp.worst_absolute_gpu);

    std::printf(
        "Worst relative: %s[%d] "
        "CPU=%.17e GPU=%.17e\n",
        cmp.worst_relative_name,
        cmp.worst_relative_index,
        cmp.worst_relative_cpu,
        cmp.worst_relative_gpu);

    if (cmp.failures != 0)
        return 1;

    std::puts("PASS M14.5b persistent CUDA trajectory parity");
    return 0;
}
