#include "pipe_sim/newton.hpp"
#include "pipe_sim/banded_solver.hpp"
#include "pipe_sim/physics.hpp"
#include "residual_internal.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace pipe_sim {

NewtonWorkspace::NewtonWorkspace(int n)
{
    if (n < 2 || n > 100)
        throw std::invalid_argument("Newton workspace requires 2 <= n <= 100");

    const int nu = 2*n + 2;

    residual.resize(nu);
    trial_residual.resize(nu);
    jacobian.resize(6*nu);
    rhs.resize(nu);
    correction.resize(nu);
    trial_state.resize(nu);
    work_state.resize(nu);
    old_friction.resize(n);
    initial_guess.resize(nu);
}

NewtonResult newton_solve_banded(
    int n,
    std::span<double> state,
    std::span<const double> old_pressure,
    std::span<const double> old_flow,
    const TransientBoundary& bc,
    const TransientParameters& par,
    const NewtonOptions& opt,
    NewtonWorkspace& ws)
{
    if (n < 2 || n > 100)
        return {3,0};

    const int nu = 2*n + 2;

    if (state.size() != static_cast<std::size_t>(nu) ||
        old_pressure.size() != static_cast<std::size_t>(n+1) ||
        old_flow.size() != static_cast<std::size_t>(n)) {
        throw std::invalid_argument("Invalid Newton input dimensions");
    }

    if (ws.residual.size() != static_cast<std::size_t>(nu) ||
        ws.trial_residual.size() != static_cast<std::size_t>(nu) ||
        ws.jacobian.size() != static_cast<std::size_t>(6*nu) ||
        ws.rhs.size() != static_cast<std::size_t>(nu) ||
        ws.correction.size() != static_cast<std::size_t>(nu) ||
        ws.trial_state.size() != static_cast<std::size_t>(nu) ||
        ws.work_state.size() != static_cast<std::size_t>(nu) ||
        ws.old_friction.size() != static_cast<std::size_t>(n)) {
        throw std::invalid_argument("Newton workspace dimensions mismatch");
    }

    // Provisional validation. The full PIPE_VALID translation
    // will be added after its Fortran implementation is inspected.
    if (!std::isfinite(par.dx) || par.dx <= 0.0 ||
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
        !std::isfinite(opt.step_tolerance) ||
        opt.residual_tolerance <= 0.0 ||
        opt.step_tolerance < 0.0 ||
        opt.max_iterations < 0) {
        return {3,0};
    }

    if (!std::isfinite(bc.inlet_flow_old) ||
        !std::isfinite(bc.outlet_flow_old) ||
        !std::isfinite(bc.outlet_flow_new) ||
        !std::isfinite(bc.inlet_pressure_new) ||
        bc.inlet_pressure_new <= 0.0) {
        return {3,0};
    }

    for (double value : state) {
        if (!std::isfinite(value))
            return {3,0};
    }

    for (double value : old_pressure) {
        if (!std::isfinite(value) || value <= 0.0)
            return {3,0};
    }

    for (double value : old_flow) {
        if (!std::isfinite(value))
            return {3,0};
    }

    // Old-time friction is invariant during Newton iteration.
    // Evaluate it once per solve, using the original arithmetic.
    const double d = par.diameter;
    const double a = par.area;
    const double t = par.temperature;
    const double z = par.z;
    const double rs = par.gas_constant;
    const double mu = par.viscosity;
    const double eps = par.roughness;

    for (int i = 0; i < n; ++i) {
        const double rhol0 = density_cz(old_pressure[i],t,z,rs);
        const double rhor0 = density_cz(old_pressure[i+1],t,z,rs);
        const double rhof0 = 0.5*(rhol0+rhor0);

        const double re0 = reynolds_mass(old_flow[i],d,a,mu);
        const double fd0 = friction_sj(re0,eps,d);

        ws.old_friction[i] =
            friction_source(fd0,d,a,rhof0,old_flow[i]);
    }

    std::copy(state.begin(),state.end(),ws.work_state.begin());

    assemble_residual_cached(
        n,ws.work_state,old_pressure,old_flow,
        bc,par,ws.old_friction,ws.residual);

    constexpr int kl = 2;
    constexpr int ku = 1;
    constexpr int ldab = 6;

    for (int iter = 0; iter < opt.max_iterations; ++iter) {

        double rnorm = 0.0;

        for (int i = 0; i < nu; ++i)
            rnorm = std::max(rnorm,std::abs(ws.residual[i]));

        if (rnorm <= opt.residual_tolerance) {
            std::copy(
                ws.work_state.begin(),
                ws.work_state.end(),
                state.begin());

            return {0,iter};
        }

        assemble_jacobian_banded(
            n,ws.work_state,par,ws.jacobian,ldab);

        for (int i = 0; i < nu; ++i)
            ws.rhs[i] = -ws.residual[i];

        const int linfo = solve_banded(
            nu,ws.jacobian,ldab,kl,ku,
            ws.rhs,ws.correction);

        if (linfo != 0)
            return {2,iter+1};

        double snorm = 0.0;
        double uscale = 1.0;

        for (int i = 0; i < nu; ++i) {
            snorm = std::max(snorm,std::abs(ws.correction[i]));
            uscale = std::max(uscale,std::abs(ws.work_state[i]));
        }

        double lambda = 1.0;
        double rnew = 0.0;
        bool accepted = false;

        for (int ls = 0; ls < 20; ++ls) {

            for (int i = 0; i < nu; ++i) {
                ws.trial_state[i] =
                    ws.work_state[i] + lambda*ws.correction[i];
            }

            assemble_residual_cached(
                n,ws.trial_state,old_pressure,old_flow,
                bc,par,ws.old_friction,ws.trial_residual);

            rnew = 0.0;

            for (int i = 0; i < nu; ++i)
                rnew = std::max(
                    rnew,std::abs(ws.trial_residual[i]));

            if (rnew < rnorm) {
                accepted = true;
                break;
            }

            lambda = 0.5*lambda;
        }

        if (!accepted)
            return {4,iter+1};

        if (opt.verbose) {
            std::cout << std::setw(4) << iter
                      << " " << std::scientific
                      << std::setprecision(4)
                      << rnorm
                      << " " << snorm
                      << " " << lambda << '\n';
        }

        std::copy(
            ws.trial_state.begin(),
            ws.trial_state.end(),
            ws.work_state.begin());

        std::copy(
            ws.trial_residual.begin(),
            ws.trial_residual.end(),
            ws.residual.begin());

        if (lambda*snorm <= opt.step_tolerance*uscale &&
            rnew <= 10.0*opt.residual_tolerance) {

            std::copy(
                ws.work_state.begin(),
                ws.work_state.end(),
                state.begin());

            return {0,iter+1};
        }
    }

    return {1,opt.max_iterations};
}

} // namespace pipe_sim
