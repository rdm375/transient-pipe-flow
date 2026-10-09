#include "pipe_sim/newton.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

extern "C" {

void newton_solve_banded_(
    const int* n,
    double* u,
    const double* po,
    const double* mo,
    const double* mino,
    const double* mouto,
    const double* moutn,
    const double* pinn,
    const double* dx,
    const double* dt,
    const double* theta,
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
    const int* verbose,
    int* info,
    int* niter);

}

namespace {

int comparisons = 0;
int failures = 0;
double max_relative_error = 0.0;

void compare(
    double cpp,
    double f77,
    const std::string& label)
{
    ++comparisons;

    const double scale =
        std::max(std::abs(cpp), std::abs(f77));

    const double error =
        std::abs(cpp-f77);

    const double relative =
        scale > 0.0 ? error/scale : 0.0;

    max_relative_error =
        std::max(max_relative_error, relative);

    if (!std::isfinite(cpp) ||
        !std::isfinite(f77) ||
        error > 5.0e-14*scale) {

        ++failures;

        std::cerr << std::setprecision(17)
                  << "FAIL " << label
                  << " cpp=" << cpp
                  << " f77=" << f77 << '\n';
    }
}

void run_case(int n, int scenario)
{
    using namespace pipe_sim;

    const double pi = std::acos(-1.0);

    const TransientParameters par{
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

    TransientBoundary bc{
        100.0,
        100.0,
        105.0,
        8.0e6
    };

    NewtonOptions opt{
        1.0e-9,
        1.0e-12,
        30,
        false
    };

    std::vector<double> po(n+1);
    std::vector<double> mo(n);

    const int nu = 2*n+2;

    std::vector<double> cpp(nu);
    std::vector<double> f77(nu);

    // Steady-state pressure profile.
    const double re =
        100.0*par.diameter/(par.area*par.viscosity);

    const double arg =
        par.roughness/(3.7*par.diameter)
        + 5.74/std::pow(re,0.9);

    const double fd =
        0.25/std::pow(std::log10(arg),2.0);

    const double coef =
        fd*par.z*par.gas_constant*par.temperature
        *100.0*100.0
        /(par.diameter*par.area*par.area);

    for (int i = 0; i <= n; ++i) {
        po[i] = std::sqrt(
            bc.inlet_pressure_new*bc.inlet_pressure_new
            - coef*i*par.dx);
    }

    std::fill(mo.begin(),mo.end(),100.0);

    cpp[0] = po[0];
    cpp[1] = 100.0;

    for (int i = 0; i < n; ++i) {
        cpp[2*i+2] = 100.0;
        cpp[2*i+3] = po[i+1];
    }

    if (scenario == 1) {
        // Perturbed initial state.
        cpp[0] += 100.0;
        cpp[1] += 0.5;

        for (int i = 0; i < n; ++i) {
            cpp[2*i+2] += 0.25;
            cpp[2*i+3] += 100.0;
        }
    }

    if (scenario == 2) {
        // Deliberately insufficient iteration budget.
        opt.max_iterations = 0;
    }

    if (scenario == 3) {
        // Invalid boundary pressure.
        bc.inlet_pressure_new = -1.0;
    }

    if (scenario == 4) {
        // Zero step tolerance is valid in Fortran.
        opt.step_tolerance = 0.0;
    }

    if (scenario == 5) {
        // Invalid old pressure.
        po[n/2] = -1.0;
    }

    f77 = cpp;

    NewtonWorkspace ws(n);

    const NewtonResult result = newton_solve_banded(
        n,cpp,po,mo,bc,par,opt,ws);

    int info_f77 = -999;
    int niter_f77 = -999;
    const int verbose = 0;

    newton_solve_banded_(
        &n,
        f77.data(),
        po.data(),
        mo.data(),
        &bc.inlet_flow_old,
        &bc.outlet_flow_old,
        &bc.outlet_flow_new,
        &bc.inlet_pressure_new,
        &par.dx,
        &par.dt,
        &par.theta,
        &par.diameter,
        &par.area,
        &par.temperature,
        &par.z,
        &par.gas_constant,
        &par.viscosity,
        &par.roughness,
        &opt.residual_tolerance,
        &opt.step_tolerance,
        &opt.max_iterations,
        &verbose,
        &info_f77,
        &niter_f77);

    ++comparisons;

    if (result.info != info_f77) {
        ++failures;

        std::cerr << "FAIL INFO"
                  << " n=" << n
                  << " scenario=" << scenario
                  << " cpp=" << result.info
                  << " f77=" << info_f77 << '\n';
    }

    ++comparisons;

    if (result.iterations != niter_f77) {
        ++failures;

        std::cerr << "FAIL NITER"
                  << " n=" << n
                  << " scenario=" << scenario
                  << " cpp=" << result.iterations
                  << " f77=" << niter_f77 << '\n';
    }

    for (int i = 0; i < nu; ++i) {
        compare(
            cpp[i],
            f77[i],
            "state n="+std::to_string(n)
            +" scenario="+std::to_string(scenario)
            +" index="+std::to_string(i));
    }
}

} // namespace

int main()
{
    for (int n : {2,3,10,25,100}) {
        for (int scenario = 0; scenario < 6; ++scenario)
            run_case(n,scenario);
    }

    std::cout << std::scientific
              << std::setprecision(17)
              << "M12.3b Newton comparisons: "
              << comparisons << '\n'
              << "Maximum relative error: "
              << max_relative_error << '\n';

    if (failures != 0) {
        std::cerr << "FAIL M12.3b: "
                  << failures << " discrepancies\n";
        return 1;
    }

    std::cout << "PASS M12.3b Newton parity\n";
    return 0;
}
