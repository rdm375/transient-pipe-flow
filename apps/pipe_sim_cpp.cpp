#include "pipe_sim/integrate.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
    using namespace pipe_sim;

    int n = 100;
    int steps = 8000;

    if (argc > 1)
        n = std::stoi(argv[1]);

    if (argc > 2)
        steps = std::stoi(argv[2]);

    if (n < 2 || n > 100 || steps < 0) {
        std::cerr << "Usage: pipe_sim_cpp [cells 2..100] [steps >= 0]\n";
        return 2;
    }

    const double pi = std::acos(-1.0);

    const double diameter = 1.0;
    const double area = pi * diameter * diameter / 4.0;

    TransientParameters par{
        100000.0 / n,
        15.0,
        0.65,
        diameter,
        area,
        288.15,
        0.9,
        500.0,
        1.1e-5,
        4.5e-5
    };

    DemandRamp demand{
        8.0e6,
        100.0,
        105.0,
        600.0
    };

    NewtonOptions options{
        1.0e-9,
        1.0e-12,
        30,
        false
    };

    std::vector<double> pressure(n+1);
    std::vector<double> flow(n, 100.0);

    double inlet_flow = 100.0;

    const double reynolds =
        inlet_flow * diameter / (area * par.viscosity);

    const double friction_arg =
        par.roughness / (3.7 * diameter) +
        5.74 / std::pow(reynolds, 0.9);

    const double friction =
        0.25 / std::pow(std::log10(friction_arg), 2.0);

    const double coefficient =
        friction * par.z * par.gas_constant *
        par.temperature * inlet_flow * inlet_flow /
        (diameter * area * area);

    for (int i = 0; i <= n; ++i) {
        pressure[i] = std::sqrt(
            demand.inlet_pressure * demand.inlet_pressure -
            coefficient * i * par.dx);
    }

    std::vector<double> time(steps+1);
    std::vector<double> hmin(steps+1);
    std::vector<double> hout(steps+1);
    std::vector<double> hpout(steps+1);
    std::vector<double> hline(steps+1);
    std::vector<int> hnit(steps+1);

    IntegrationHistory history{
        time, hmin, hout, hpout, hline, hnit
    };

    IntegrationWorkspace workspace(n);

    const auto start = std::chrono::steady_clock::now();

    const IntegrationResult result = integrate_transient(
        n,
        pressure,
        flow,
        inlet_flow,
        demand,
        steps,
        par,
        options,
        history,
        workspace);

    const auto stop = std::chrono::steady_clock::now();

    const double elapsed =
        std::chrono::duration<double>(stop-start).count();

    if (result.info != 0) {
        std::cerr << "Integration failed: INFO="
                  << result.info
                  << " completed_steps="
                  << result.completed_steps << '\n';
        return 1;
    }

    long long total_newton_iterations = 0;

    for (int k = 1; k <= steps; ++k)
        total_newton_iterations += hnit[k];

    std::cout << std::setprecision(17);

    std::cout
        << "backend=cpp20\n"
        << "cells=" << n << '\n'
        << "steps=" << steps << '\n'
        << "duration_s=" << steps*par.dt << '\n'
        << "elapsed_s=" << elapsed << '\n'
        << "ns_per_step="
        << (steps > 0 ? elapsed*1.0e9/steps : 0.0) << '\n'
        << "final_outlet_pressure_pa=" << pressure[n] << '\n'
        << "final_inlet_flow_kg_s=" << inlet_flow << '\n'
        << "final_outlet_flow_kg_s=" << hout[steps] << '\n'
        << "final_linepack_kg=" << hline[steps] << '\n'
        << "max_mass_defect_kg_s="
        << result.max_mass_defect << '\n'
        << "total_newton_iterations="
        << total_newton_iterations << '\n';

    return 0;
}
