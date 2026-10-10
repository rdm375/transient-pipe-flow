#include "pipe_sim/integrate.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

extern "C" {

void integrate_transient_(
    const int* n,
    double* po,
    double* mo,
    double* mino,
    const double* pin,
    const double* mdot0,
    const double* mdot1,
    const double* tramp,
    const double* dt,
    const int* nsteps,
    const double* theta,
    const double* dx,
    const double* d,
    const double* a,
    const double* t,
    const double* z,
    const double* rs,
    const double* mu,
    const double* eps,
    const double* rtol,
    const double* stol,
    const int* maxit,
    double* time,
    double* hmin,
    double* hout,
    double* hpout,
    double* hline,
    int* hnit,
    double* hmaxbal,
    int* info);

}

namespace {

int comparisons = 0;
int failures = 0;
double maximum_relative_error = 0.0;

void compare(double cpp, double f77, const std::string& label)
{
    ++comparisons;

    const double scale = std::max(std::abs(cpp), std::abs(f77));
    const double error = std::abs(cpp - f77);
    const double relative = scale > 0.0 ? error / scale : 0.0;

    maximum_relative_error =
        std::max(maximum_relative_error, relative);

    if (!std::isfinite(cpp) || !std::isfinite(f77) ||
        error > 5.0e-14 * scale) {
        ++failures;
        if (failures <= 30) {
            std::cerr << std::setprecision(17)
                      << "FAIL " << label
                      << " cpp=" << cpp
                      << " f77=" << f77 << '\n';
        }
    }
}

void compare_int(int cpp, int f77, const std::string& label)
{
    ++comparisons;

    if (cpp != f77) {
        ++failures;
        if (failures <= 30) {
            std::cerr << "FAIL " << label
                      << " cpp=" << cpp
                      << " f77=" << f77 << '\n';
        }
    }
}

struct History {
    std::vector<double> time, min, mout, pout, line;
    std::vector<int> nit;

    explicit History(int steps)
        : time(steps+1, -999.0),
          min(steps+1, -999.0),
          mout(steps+1, -999.0),
          pout(steps+1, -999.0),
          line(steps+1, -999.0),
          nit(steps+1, -999)
    {
    }

    pipe_sim::IntegrationHistory spans()
    {
        return {time, min, mout, pout, line, nit};
    }
};

void run_case(int n, int steps, int scenario)
{
    using namespace pipe_sim;

    const double pi = std::acos(-1.0);

    TransientParameters par{
        100000.0/n,
        15.0,
        0.65,
        1.0,
        pi/4.0,
        288.15,
        0.9,
        500.0,
        1.1e-5,
        4.5e-5
    };

    DemandRamp demand{8.0e6, 100.0, 105.0, 600.0};
    NewtonOptions options{1.0e-9, 1.0e-12, 30, false};

    if (scenario == 1)
        demand.final_outlet_flow = 95.0;
    if (scenario == 2)
        demand.ramp_duration = 0.0;
    if (scenario == 3)
        demand.inlet_pressure = -1.0;
    if (scenario == 4)
        options.max_iterations = 0;
    if (scenario == 5)
        par.theta = 1.5;

    std::vector<double> pressure(n+1);
    std::vector<double> flow(n, 100.0);

    const double re =
        100.0 * par.diameter / (par.area * par.viscosity);

    const double arg =
        par.roughness / (3.7 * par.diameter) +
        5.74 / std::pow(re, 0.9);

    const double fd =
        0.25 / std::pow(std::log10(arg), 2.0);

    const double coef =
        fd * par.z * par.gas_constant * par.temperature *
        10000.0 / (par.diameter * par.area * par.area);

    for (int i = 0; i <= n; ++i) {
        pressure[i] = std::sqrt(
            8.0e6 * 8.0e6 - coef*i*par.dx);
    }

    std::vector<double> pressure_f77 = pressure;
    std::vector<double> flow_f77 = flow;

    double inlet = 100.0;
    double inlet_f77 = inlet;

    History cpp(steps);
    History f77(steps);

    IntegrationWorkspace workspace(n);

    const IntegrationResult result = integrate_transient(
        n, pressure, flow, inlet,
        demand, steps, par, options,
        cpp.spans(), workspace);

    double maxbal_f77 = -999.0;
    int info_f77 = -999;

    integrate_transient_(
        &n,
        pressure_f77.data(),
        flow_f77.data(),
        &inlet_f77,
        &demand.inlet_pressure,
        &demand.initial_outlet_flow,
        &demand.final_outlet_flow,
        &demand.ramp_duration,
        &par.dt,
        &steps,
        &par.theta,
        &par.dx,
        &par.diameter,
        &par.area,
        &par.temperature,
        &par.z,
        &par.gas_constant,
        &par.viscosity,
        &par.roughness,
        &options.residual_tolerance,
        &options.step_tolerance,
        &options.max_iterations,
        f77.time.data(),
        f77.min.data(),
        f77.mout.data(),
        f77.pout.data(),
        f77.line.data(),
        f77.nit.data(),
        &maxbal_f77,
        &info_f77);

    const std::string label =
        " n=" + std::to_string(n) +
        " steps=" + std::to_string(steps) +
        " scenario=" + std::to_string(scenario);

    compare_int(result.info, info_f77, "INFO" + label);

    for (int i = 0; i <= n; ++i) {
        compare(pressure[i], pressure_f77[i],
                "final pressure " + std::to_string(i) + label);
    }

    for (int i = 0; i < n; ++i) {
        compare(flow[i], flow_f77[i],
                "final flow " + std::to_string(i) + label);
    }

    compare(inlet, inlet_f77, "final inlet" + label);

    // Fortran does not initialize history on invalid input.
    // On a failed timestep, only completed records are defined.
    if (info_f77 != 3) {
        compare(result.max_mass_defect, maxbal_f77,
                "maximum mass defect" + label);

        const int completed = result.completed_steps;

        for (int k = 0; k <= completed; ++k) {
            const std::string at =
                " k=" + std::to_string(k) + label;

            compare(cpp.time[k], f77.time[k], "time" + at);
            compare(cpp.min[k], f77.min[k], "inlet" + at);
            compare(cpp.mout[k], f77.mout[k], "outlet" + at);
            compare(cpp.pout[k], f77.pout[k], "pressure" + at);
            compare(cpp.line[k], f77.line[k], "linepack" + at);
            compare_int(cpp.nit[k], f77.nit[k], "iterations" + at);
        }

        if (info_f77 == 0)
            compare_int(completed, steps, "completed steps" + label);
    }
}

} // namespace

int main()
{
    for (int n : {2, 3, 10, 25, 100}) {
        for (int steps : {0, 1, 10, 80}) {
            for (int scenario = 0; scenario < 6; ++scenario) {
                run_case(n, steps, scenario);
            }
        }
    }

    // M12.5a: M11 long-duration reference workload.
    run_case(100, 8000, 0);

    std::cout << std::scientific
              << std::setprecision(17)
              << "M12.4b integration comparisons: "
              << comparisons << '\n'
              << "Maximum relative error: "
              << maximum_relative_error << '\n';

    if (failures != 0) {
        std::cerr << "FAIL M12.4b: "
                  << failures << " discrepancies\n";
        return 1;
    }

    std::cout << "PASS M12.4b integration parity\n";
    return 0;
}
