#pragma once

#include "pipe_sim/newton.hpp"

#include <span>

namespace pipe_sim {

NewtonResult transient_step(
    int n,
    std::span<const double> old_pressure,
    std::span<const double> old_flow,
    const TransientBoundary& boundary,
    const TransientParameters& parameters,
    const NewtonOptions& options,
    std::span<double> new_state,
    NewtonWorkspace& workspace);

} // namespace pipe_sim
