#pragma once

#include "pipe_sim/transient_step.hpp"

#include <span>
#include <vector>

namespace pipe_sim {

struct DemandRamp {
    double inlet_pressure;
    double initial_outlet_flow;
    double final_outlet_flow;
    double ramp_duration;
};

struct IntegrationResult {
    int info = 0;
    int completed_steps = 0;
    double max_mass_defect = 0.0;
};

struct IntegrationHistory {
    std::span<double> time;
    std::span<double> inlet_flow;
    std::span<double> outlet_flow;
    std::span<double> outlet_pressure;
    std::span<double> linepack;
    std::span<int> newton_iterations;
};

struct IntegrationWorkspace {
    NewtonWorkspace newton;
    std::vector<double> state;
    std::vector<double> new_pressure;
    std::vector<double> new_flow;

    explicit IntegrationWorkspace(int n);
};

double demand_ramp(
    double time,
    double initial_flow,
    double final_flow,
    double ramp_duration);

double linepack_cz(
    int n,
    std::span<const double> pressure,
    const TransientParameters& parameters);

IntegrationResult integrate_transient(
    int n,
    std::span<double> pressure,
    std::span<double> flow,
    double& inlet_flow,
    const DemandRamp& demand,
    int steps,
    const TransientParameters& parameters,
    const NewtonOptions& options,
    const IntegrationHistory& history,
    IntegrationWorkspace& workspace);

} // namespace pipe_sim
