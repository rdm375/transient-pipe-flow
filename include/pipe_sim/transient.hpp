#pragma once

#include <span>

namespace pipe_sim {

struct TransientParameters {
    double dx;
    double dt;
    double theta;
    double diameter;
    double area;
    double temperature;
    double z;
    double gas_constant;
    double viscosity;
    double roughness;
};

struct TransientBoundary {
    double inlet_flow_old;
    double outlet_flow_old;
    double outlet_flow_new;
    double inlet_pressure_new;
};

// N intervals:
//   old_pressure: N+1
//   old_flow:     N
//   state:        2*N+2
//   residual:     2*N+2
void assemble_residual(
    int n,
    std::span<const double> state,
    std::span<const double> old_pressure,
    std::span<const double> old_flow,
    const TransientBoundary& boundary,
    const TransientParameters& parameters,
    std::span<double> residual);

// LAPACK-style column-major band storage.
// KL=2, KU=1, KV=KL+KU=3.
// Storage index for zero-based row i, column j:
//     ab[j*ldab + (KV+i-j)]
// Requires ldab >= 6 and ab.size() == ldab*(2*n+2).
void assemble_jacobian_banded(
    int n,
    std::span<const double> state,
    const TransientParameters& parameters,
    std::span<double> ab,
    int ldab);

} // namespace pipe_sim
