#pragma once

#include "pipe_sim/transient.hpp"

#include <span>
#include <vector>

namespace pipe_sim {

struct NewtonOptions {
    double residual_tolerance;
    double step_tolerance;
    int max_iterations;
    bool verbose = false;
};

struct NewtonResult {
    int info = 0;
    int iterations = 0;
};

struct NewtonWorkspace {
    std::vector<double> residual;
    std::vector<double> trial_residual;
    std::vector<double> jacobian;
    std::vector<double> rhs;
    std::vector<double> correction;
    std::vector<double> trial_state;
    std::vector<double> work_state;
    std::vector<double> old_friction;
    std::vector<double> old_friction_factor;
    std::vector<double> initial_guess;

    explicit NewtonWorkspace(int n);
};

NewtonResult newton_solve_banded(
    int n,
    std::span<double> state,
    std::span<const double> old_pressure,
    std::span<const double> old_flow,
    const TransientBoundary& boundary,
    const TransientParameters& parameters,
    const NewtonOptions& options,
    NewtonWorkspace& workspace);

} // namespace pipe_sim
