#include "pipe_sim/physics.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

extern "C" {
double eos_cz_rho_(const double*, const double*,
                   const double*, const double*);
double eos_cz_drhodp_(const double*, const double*,
                      const double*);
double reynolds_mass_(const double*, const double*,
                      const double*, const double*);
double dreynolds_dmdot_(const double*, const double*,
                        const double*, const double*);
double friction_sj_(const double*, const double*, const double*);
double dfriction_sj_dre_(const double*, const double*, const double*);
void friction_sj_value_derivative_(
    const double*, const double*, const double*,
    double*, double*);
double friction_source_(const double*, const double*,
                        const double*, const double*,
                        const double*);
double dfriction_source_dmdot_(
    const double*, const double*, const double*,
    const double*, const double*, const double*);
}

namespace {

int failures = 0;
int comparisons = 0;
double max_relative_error = 0.0;

void compare(const std::string& name, double cpp, double f77)
{
    ++comparisons;

    if (!std::isfinite(cpp) || !std::isfinite(f77)) {
        std::cerr << "FAIL " << name << ": nonfinite result\n";
        ++failures;
        return;
    }

    const double scale = std::max(std::abs(cpp), std::abs(f77));
    const double error = std::abs(cpp - f77);
    const double relative = scale > 0.0 ? error / scale : 0.0;

    max_relative_error = std::max(max_relative_error, relative);

    // Tight numerical parity, without requiring identical libm results.
    if (error > 5.0e-14 * scale) {
        std::cerr << std::setprecision(17)
                  << "FAIL " << name
                  << " cpp=" << cpp
                  << " f77=" << f77
                  << " relative=" << relative << '\n';
        ++failures;
    }
}

} // namespace

int main()
{
    using namespace pipe_sim;

    const double temperature = 288.15;
    const double z = 0.9;
    const double rs = 500.0;
    const double diameter = 1.0;
    const double area = std::acos(-1.0) * diameter * diameter / 4.0;
    const double viscosity = 1.1e-5;

    for (double pressure : {1.0e5, 1.0e6, 8.0e6, 2.0e7}) {
        compare("density",
                density_cz(pressure, temperature, z, rs),
                eos_cz_rho_(&pressure, &temperature, &z, &rs));
    }

    compare("density derivative",
            density_derivative_cz(temperature, z, rs),
            eos_cz_drhodp_(&temperature, &z, &rs));

    for (double mdot : {-200.0, -100.0, -1.0, -0.0,
                         0.0, 1.0, 100.0, 200.0}) {
        compare("Reynolds",
                reynolds_mass(mdot, diameter, area, viscosity),
                reynolds_mass_(&mdot, &diameter, &area, &viscosity));

        compare("Reynolds derivative",
                reynolds_derivative_mdot(mdot, diameter, area, viscosity),
                dreynolds_dmdot_(&mdot, &diameter, &area, &viscosity));
    }

    for (double reynolds : {1.0e3, 1.0e4, 1.0e5,
                             1.0e6, 1.0e7, 1.0e9}) {
        for (double roughness : {0.0, 1.0e-6, 4.5e-5, 1.0e-3}) {
            compare("friction",
                    friction_sj(reynolds, roughness, diameter),
                    friction_sj_(&reynolds, &roughness, &diameter));

            compare("friction derivative",
                    friction_derivative_re(reynolds, roughness, diameter),
                    dfriction_sj_dre_(&reynolds, &roughness, &diameter));

            double fd_f77 = 0.0;
            double dfdre_f77 = 0.0;

            friction_sj_value_derivative_(
                &reynolds, &roughness, &diameter,
                &fd_f77, &dfdre_f77);

            const auto result =
                friction_sj_value_derivative(
                    reynolds, roughness, diameter);

            compare("combined friction", result.value, fd_f77);
            compare("combined derivative",
                    result.derivative_re, dfdre_f77);

            const double density_face = 60.0;
            const double mdot = 100.0;
            const double dfdm = 1.0e-5;

            compare("friction source",
                    friction_source(result.value, diameter, area,
                                    density_face, mdot),
                    friction_source_(&fd_f77, &diameter, &area,
                                     &density_face, &mdot));

            compare("friction source derivative",
                    friction_source_derivative_mdot(
                        result.value, dfdm, diameter, area,
                        density_face, mdot),
                    dfriction_source_dmdot_(
                        &fd_f77, &dfdm, &diameter, &area,
                        &density_face, &mdot));
        }
    }

    std::cout << std::scientific << std::setprecision(17)
              << "M12.1 physics parity comparisons: "
              << comparisons << '\n'
              << "Maximum relative error: "
              << max_relative_error << '\n';

    if (failures != 0) {
        std::cerr << "FAIL M12.1: "
                  << failures << " discrepancies\n";
        return 1;
    }

    std::cout << "PASS M12.1 physics parity\n";
    return 0;
}
