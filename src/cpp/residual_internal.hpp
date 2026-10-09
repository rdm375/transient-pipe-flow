#pragma once

#include "pipe_sim/transient.hpp"

#include <span>

namespace pipe_sim {

void assemble_residual_cached(
    int n,
    std::span<const double> state,
    std::span<const double> old_pressure,
    std::span<const double> old_flow,
    const TransientBoundary& boundary,
    const TransientParameters& parameters,
    std::span<const double> old_friction,
    std::span<double> residual);

} // namespace pipe_sim
