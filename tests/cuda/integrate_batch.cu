#include "pipe_sim/integrate.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace pipe_sim {
extern "C" cudaError_t launch_pipe_integrate_batch(
    int, int, double*, double*, double*,
    DemandRamp, int, TransientParameters, NewtonOptions,
    double*, double*, double*, double*, double*,
    int*, IntegrationResult*);
}

#define CHECK(call) do {                                  \
    cudaError_t err = (call);                             \
    if (err != cudaSuccess) {                             \
        std::fprintf(stderr, "%s:%d: %s\n",               \
            __FILE__, __LINE__, cudaGetErrorString(err)); \
        return 1;                                         \
    }                                                     \
} while (0)

template<class T>
struct DeviceBuffer {
    T* ptr = nullptr;

    explicit DeviceBuffer(std::size_t count) {
        cudaError_t err = cudaMalloc(
            reinterpret_cast<void**>(&ptr),
            count*sizeof(T));

        if (err != cudaSuccess) {
            std::fprintf(stderr,
                "cudaMalloc: %s\n",
                cudaGetErrorString(err));
            std::abort();
        }
    }

    ~DeviceBuffer() {
        if (ptr)
            cudaFree(ptr);
    }

    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
};

int run_batch(int batch, int steps)
{
    using namespace pipe_sim;

    constexpr int n = 100;

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

    std::vector<double> initial_pressure(n+1);
    std::vector<double> initial_flow(n,100.0);

    const double inlet = 100.0;

    const double reynolds =
        inlet*diameter/(area*par.viscosity);

    const double friction_arg =
        par.roughness/(3.7*diameter) +
        5.74/std::pow(reynolds,0.9);

    const double friction =
        0.25/std::pow(std::log10(friction_arg),2.0);

    const double coefficient =
        friction*par.z*par.gas_constant*
        par.temperature*inlet*inlet/
        (diameter*area*area);

    for (int i=0; i<=n; ++i) {
        initial_pressure[i] = std::sqrt(
            demand.inlet_pressure*demand.inlet_pressure -
            coefficient*i*par.dx);
    }

    // CPU reference: one independent simulation.
    auto cpu_pressure = initial_pressure;
    auto cpu_flow = initial_flow;
    double cpu_inlet = inlet;

    const int stride = steps+1;

    std::vector<double> ct(stride),ci(stride),co(stride);
    std::vector<double> cp(stride),cl(stride);
    std::vector<int> cn(stride);

    IntegrationHistory history{
        ct,ci,co,cp,cl,cn
    };

    IntegrationWorkspace workspace(n);

    const auto cpu_start =
        std::chrono::steady_clock::now();

    const IntegrationResult cpu = integrate_transient(
        n,cpu_pressure,cpu_flow,cpu_inlet,
        demand,steps,par,opt,history,workspace);

    const auto cpu_stop =
        std::chrono::steady_clock::now();

    const double cpu_ms =
        std::chrono::duration<double,std::milli>(
            cpu_stop-cpu_start).count();

    std::printf(
        "CPU_REFERENCE steps=%d elapsed_ms=%.6f "
        "simulations_per_second=%.3f\n",
        steps,cpu_ms,1000.0/cpu_ms);

    if (cpu.info != 0 || cpu.completed_steps != steps) {
        std::fprintf(stderr,
            "CPU reference failed: info=%d steps=%d\n",
            cpu.info,cpu.completed_steps);
        return 1;
    }

    // Structure-of-arrays within each independent pipe.
    std::vector<double> pressure(batch*(n+1));
    std::vector<double> flow(batch*n);
    std::vector<double> inlet_flow(batch,inlet);

    for (int b=0; b<batch; ++b) {
        std::copy(
            initial_pressure.begin(),
            initial_pressure.end(),
            pressure.begin()+b*(n+1));

        std::copy(
            initial_flow.begin(),
            initial_flow.end(),
            flow.begin()+b*n);
    }

    const std::size_t hcount =
        static_cast<std::size_t>(batch)*stride;

    DeviceBuffer<double> dp(pressure.size());
    DeviceBuffer<double> df(flow.size());
    DeviceBuffer<double> di(inlet_flow.size());

    DeviceBuffer<double> dt(hcount);
    DeviceBuffer<double> dhi(hcount);
    DeviceBuffer<double> dho(hcount);
    DeviceBuffer<double> dhp(hcount);
    DeviceBuffer<double> dhl(hcount);

    DeviceBuffer<int> dhn(hcount);
    DeviceBuffer<IntegrationResult> dr(batch);

    CHECK(cudaMemcpy(
        dp.ptr,pressure.data(),
        pressure.size()*sizeof(double),
        cudaMemcpyHostToDevice));

    CHECK(cudaMemcpy(
        df.ptr,flow.data(),
        flow.size()*sizeof(double),
        cudaMemcpyHostToDevice));

    CHECK(cudaMemcpy(
        di.ptr,inlet_flow.data(),
        inlet_flow.size()*sizeof(double),
        cudaMemcpyHostToDevice));

    cudaEvent_t start,stop;

    CHECK(cudaEventCreate(&start));
    CHECK(cudaEventCreate(&stop));

    CHECK(cudaEventRecord(start));

    CHECK(launch_pipe_integrate_batch(
        n,batch,
        dp.ptr,df.ptr,di.ptr,
        demand,steps,par,opt,
        dt.ptr,dhi.ptr,dho.ptr,dhp.ptr,dhl.ptr,
        dhn.ptr,dr.ptr));

    CHECK(cudaEventRecord(stop));
    CHECK(cudaEventSynchronize(stop));

    float elapsed_ms = 0.0f;

    CHECK(cudaEventElapsedTime(
        &elapsed_ms,start,stop));

    CHECK(cudaEventDestroy(start));
    CHECK(cudaEventDestroy(stop));

    std::vector<IntegrationResult> results(batch);

    CHECK(cudaMemcpy(
        results.data(),dr.ptr,
        batch*sizeof(IntegrationResult),
        cudaMemcpyDeviceToHost));

    CHECK(cudaMemcpy(
        pressure.data(),dp.ptr,
        pressure.size()*sizeof(double),
        cudaMemcpyDeviceToHost));

    CHECK(cudaMemcpy(
        flow.data(),df.ptr,
        flow.size()*sizeof(double),
        cudaMemcpyDeviceToHost));

    CHECK(cudaMemcpy(
        inlet_flow.data(),di.ptr,
        inlet_flow.size()*sizeof(double),
        cudaMemcpyDeviceToHost));

    std::vector<double> gpu_time(hcount);
    std::vector<double> gpu_in(hcount);
    std::vector<double> gpu_out(hcount);
    std::vector<double> gpu_pressure_history(hcount);
    std::vector<double> gpu_linepack(hcount);
    std::vector<int> gpu_iterations(hcount);

    CHECK(cudaMemcpy(
        gpu_time.data(),dt.ptr,
        hcount*sizeof(double),
        cudaMemcpyDeviceToHost));

    CHECK(cudaMemcpy(
        gpu_in.data(),dhi.ptr,
        hcount*sizeof(double),
        cudaMemcpyDeviceToHost));

    CHECK(cudaMemcpy(
        gpu_out.data(),dho.ptr,
        hcount*sizeof(double),
        cudaMemcpyDeviceToHost));

    CHECK(cudaMemcpy(
        gpu_pressure_history.data(),dhp.ptr,
        hcount*sizeof(double),
        cudaMemcpyDeviceToHost));

    CHECK(cudaMemcpy(
        gpu_linepack.data(),dhl.ptr,
        hcount*sizeof(double),
        cudaMemcpyDeviceToHost));

    CHECK(cudaMemcpy(
        gpu_iterations.data(),dhn.ptr,
        hcount*sizeof(int),
        cudaMemcpyDeviceToHost));

    long long comparisons = 0;
    long long failures = 0;
    double max_abs = 0.0;

    auto compare = [&](double expected, double actual) {
        ++comparisons;

        if (!std::isfinite(actual)) {
            ++failures;
            return;
        }

        const double error = std::abs(expected-actual);

        max_abs = std::max(max_abs,error);

        const double tolerance = std::max(
            1.0e-9,
            5.0e-12*std::max(
                std::abs(expected),std::abs(actual)));

        if (error > tolerance)
            ++failures;
    };

    for (int b=0; b<batch; ++b) {
        const auto& r = results[b];

        if (r.info != cpu.info ||
            r.completed_steps != cpu.completed_steps)
            ++failures;

        compare(cpu.max_mass_defect,r.max_mass_defect);
        compare(cpu_inlet,inlet_flow[b]);

        for (int i=0; i<=n; ++i)
            compare(
                cpu_pressure[i],
                pressure[b*(n+1)+i]);

        for (int i=0; i<n; ++i)
            compare(
                cpu_flow[i],
                flow[b*n+i]);

        for (int k=0; k<stride; ++k) {
            const int j=b*stride+k;

            compare(ct[k],gpu_time[j]);
            compare(ci[k],gpu_in[j]);
            compare(co[k],gpu_out[j]);
            compare(cp[k],gpu_pressure_history[j]);
            compare(cl[k],gpu_linepack[j]);

            ++comparisons;

            if (cn[k] != gpu_iterations[j])
                ++failures;
        }
    }

    const double seconds =
        static_cast<double>(elapsed_ms)*1.0e-3;

    const double throughput =
        seconds > 0.0 ? batch/seconds : 0.0;

    std::printf(
        "BATCH batch=%3d steps=%5d "
        "kernel_ms=%12.4f "
        "simulations_per_second=%10.4f "
        "comparisons=%lld failures=%lld "
        "max_abs=%.3e\n",
        batch,steps,
        static_cast<double>(elapsed_ms),
        throughput,
        comparisons,failures,max_abs);

    return failures == 0 ? 0 : 1;
}

int main(int argc, char** argv)
{
    int steps = 80;

    if (argc > 1)
        steps = std::atoi(argv[1]);

    if (steps < 0) {
        std::fprintf(stderr,"Invalid steps\n");
        return 1;
    }

    for (int batch : {512,1024,2048}) {
        if (run_batch(batch,steps) != 0)
            return 1;
    }

    std::puts("PASS M14.6b batch parity and throughput");

    return 0;
}
