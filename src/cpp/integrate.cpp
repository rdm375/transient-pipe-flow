#include "pipe_sim/integrate.hpp"
#include "pipe_sim/physics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pipe_sim {

IntegrationWorkspace::IntegrationWorkspace(int n)
    : newton(n),
      state(2*n + 2),
      new_pressure(n + 1),
      new_flow(n)
{
}

double demand_ramp(
    double time,
    double initial_flow,
    double final_flow,
    double ramp_duration)
{
    if (ramp_duration <= 0.0 || time >= ramp_duration)
        return final_flow;

    if (time <= 0.0)
        return initial_flow;

    const double fraction = time / ramp_duration;

    return initial_flow +
           fraction * (final_flow - initial_flow);
}

double linepack_cz(
    int n,
    std::span<const double> pressure,
    const TransientParameters& par)
{
    if (n < 2 || n > 100 ||
        pressure.size() != static_cast<std::size_t>(n+1)) {
        throw std::invalid_argument("Invalid linepack dimensions");
    }

    double mass = 0.0;

    for (int i = 0; i <= n; ++i) {
        double weight = 1.0;

        if (i == 0 || i == n)
            weight = 0.5;

        mass = mass +
            weight * density_cz(
                pressure[i],
                par.temperature,
                par.z,
                par.gas_constant);
    }

    mass = par.area * par.dx * mass;

    return mass;
}

IntegrationResult integrate_transient(
    int n,
    std::span<double> pressure,
    std::span<double> flow,
    double& inlet_flow,
    const DemandRamp& demand,
    int steps,
    const TransientParameters& par,
    const NewtonOptions& opt,
    const IntegrationHistory& hist,
    IntegrationWorkspace& ws)
{
    if (n < 2 || n > 100 || steps < 0)
        return {3, 0, 0.0};

    if (pressure.size() != static_cast<std::size_t>(n+1) ||
        flow.size() != static_cast<std::size_t>(n)) {
        throw std::invalid_argument("Invalid integration state dimensions");
    }

    const auto required = static_cast<std::size_t>(steps) + 1;

    if (hist.time.size() < required ||
        hist.inlet_flow.size() < required ||
        hist.outlet_flow.size() < required ||
        hist.outlet_pressure.size() < required ||
        hist.linepack.size() < required ||
        hist.newton_iterations.size() < required) {
        throw std::invalid_argument("Insufficient integration history capacity");
    }

    if (ws.state.size() != static_cast<std::size_t>(2*n+2) ||
        ws.new_pressure.size() != static_cast<std::size_t>(n+1) ||
        ws.new_flow.size() != static_cast<std::size_t>(n)) {
        throw std::invalid_argument("Integration workspace dimensions mismatch");
    }

    // Reproduce PIPE_VALID and the additional integrator checks.
    if (!std::isfinite(inlet_flow) ||
        !std::isfinite(par.dx) || par.dx <= 0.0 ||
        !std::isfinite(par.dt) || par.dt <= 0.0 ||
        !std::isfinite(par.theta) ||
        par.theta < 0.0 || par.theta > 1.0 ||
        !std::isfinite(par.diameter) || par.diameter <= 0.0 ||
        !std::isfinite(par.area) || par.area <= 0.0 ||
        !std::isfinite(par.temperature) || par.temperature <= 0.0 ||
        !std::isfinite(par.z) || par.z <= 0.0 ||
        !std::isfinite(par.gas_constant) || par.gas_constant <= 0.0 ||
        !std::isfinite(par.viscosity) || par.viscosity <= 0.0 ||
        !std::isfinite(par.roughness) || par.roughness < 0.0 ||
        !std::isfinite(opt.residual_tolerance) ||
        opt.residual_tolerance <= 0.0 ||
        !std::isfinite(opt.step_tolerance) ||
        opt.step_tolerance < 0.0 ||
        opt.max_iterations < 0 ||
        !std::isfinite(demand.inlet_pressure) ||
        demand.inlet_pressure <= 0.0 ||
        !std::isfinite(demand.initial_outlet_flow) ||
        !std::isfinite(demand.final_outlet_flow) ||
        !std::isfinite(demand.ramp_duration) ||
        demand.ramp_duration < 0.0) {
        return {3, 0, 0.0};
    }

    for (double value : pressure) {
        if (!std::isfinite(value) || value <= 0.0)
            return {3, 0, 0.0};
    }

    for (double value : flow) {
        if (!std::isfinite(value))
            return {3, 0, 0.0};
    }

    double old_outlet_flow = demand.initial_outlet_flow;

    hist.time[0] = 0.0;
    hist.inlet_flow[0] = inlet_flow;
    hist.outlet_flow[0] = old_outlet_flow;
    hist.outlet_pressure[0] = pressure[n];
    hist.linepack[0] = linepack_cz(n, pressure, par);
    hist.newton_iterations[0] = 0;

    double max_mass_defect = 0.0;

    for (int k = 1; k <= steps; ++k) {
        const double new_time =
            static_cast<double>(k) * par.dt;

        const double new_outlet_flow = demand_ramp(
            new_time,
            demand.initial_outlet_flow,
            demand.final_outlet_flow,
            demand.ramp_duration);

        const TransientBoundary bc{
            inlet_flow,
            old_outlet_flow,
            new_outlet_flow,
            demand.inlet_pressure
        };

        const NewtonResult step = transient_step(
            n,
            pressure,
            flow,
            bc,
            par,
            opt,
            ws.state,
            ws.newton);

        if (step.info != 0)
            return {step.info, k-1, max_mass_defect};

        ws.new_pressure[0] = ws.state[0];

        for (int i = 0; i < n; ++i) {
            ws.new_flow[i] = ws.state[2*i+2];
            ws.new_pressure[i+1] = ws.state[2*i+3];
        }

        hist.time[k] = new_time;
        hist.inlet_flow[k] = ws.state[1];
        hist.outlet_flow[k] = new_outlet_flow;
        hist.outlet_pressure[k] = ws.new_pressure[n];

        hist.linepack[k] =
            linepack_cz(n, ws.new_pressure, par);

        hist.newton_iterations[k] = step.iterations;

        // Preserve the Fortran arithmetic grouping.
        const double balance_error = std::abs(
            (hist.linepack[k] - hist.linepack[k-1]) / par.dt
            - (par.theta *
                (hist.inlet_flow[k] - hist.outlet_flow[k])
              + (1.0 - par.theta) *
                (hist.inlet_flow[k-1] - hist.outlet_flow[k-1])));

        max_mass_defect =
            std::max(max_mass_defect, balance_error);

        std::copy(
            ws.new_pressure.begin(),
            ws.new_pressure.end(),
            pressure.begin());

        std::copy(
            ws.new_flow.begin(),
            ws.new_flow.end(),
            flow.begin());

        inlet_flow = hist.inlet_flow[k];
        old_outlet_flow = new_outlet_flow;
    }

    return {0, steps, max_mass_defect};
}

} // namespace pipe_sim
