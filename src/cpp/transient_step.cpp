#include "pipe_sim/transient_step.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pipe_sim {

NewtonResult transient_step(
    int n,
    std::span<const double> old_pressure,
    std::span<const double> old_flow,
    const TransientBoundary& bc,
    const TransientParameters& par,
    const NewtonOptions& opt,
    std::span<double> new_state,
    NewtonWorkspace& ws)
{
    if (n < 2 || n > 100)
        return {3, 0};

    const int nu = 2*n + 2;

    if (old_pressure.size() != static_cast<std::size_t>(n+1) ||
        old_flow.size() != static_cast<std::size_t>(n) ||
        new_state.size() != static_cast<std::size_t>(nu)) {
        throw std::invalid_argument("Invalid transient step dimensions");
    }

    // Match the Fortran entry-point validation.
    if (!std::isfinite(bc.inlet_flow_old) ||
        !std::isfinite(bc.outlet_flow_old) ||
        !std::isfinite(bc.outlet_flow_new) ||
        !std::isfinite(bc.inlet_pressure_new) ||
        bc.inlet_pressure_new <= 0.0) {
        return {3, 0};
    }

    // Newton validates the remaining physical parameters and old state.
    // Its caller-visible state is unchanged on failure.
    //
    // Use a separate temporary state to preserve the Fortran contract:
    // NEW_STATE is modified only when INFO == 0.
    //
    // This vector is allocated once per call in the initial port.
    // We will eliminate that allocation during M13 optimization.
    std::vector<double> guess(nu);

    guess[0] = bc.inlet_pressure_new;
    guess[1] = bc.inlet_flow_old;

    for (int i = 0; i < n; ++i) {
        guess[2*i+2] = old_flow[i];
        guess[2*i+3] = old_pressure[i+1];
    }

    const NewtonResult result = newton_solve_banded(
        n, guess, old_pressure, old_flow, bc, par, opt, ws);

    if (result.info == 0) {
        std::copy(
            guess.begin(),
            guess.end(),
            new_state.begin());
    }

    return result;
}

} // namespace pipe_sim
