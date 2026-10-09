#include "pipe_sim/transient.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

extern "C" {
void assemble_residual_(
    const int*, const double*, const double*, const double*,
    const double*, const double*, const double*, const double*,
    const double*, const double*, const double*, const double*,
    const double*, const double*, const double*, const double*,
    const double*, const double*, double*);
}

int main()
{
    using namespace pipe_sim;

    int comparisons = 0;
    int failures = 0;
    double max_relative_error = 0.0;

    for (int n : {2, 3, 10, 25, 100}) {
        const double pi = std::acos(-1.0);

        const TransientParameters par{
            100000.0/n, 15.0, 0.65,
            1.0, pi/4.0, 288.15, 0.9,
            500.0, 1.1e-5, 4.5e-5
        };

        const TransientBoundary bc{
            100.0, 100.0, 105.0, 8.0e6
        };

        std::vector<double> po(n+1);
        std::vector<double> mo(n);
        std::vector<double> u(2*n+2);
        std::vector<double> r_cpp(2*n+2);
        std::vector<double> r_f77(2*n+2);

        for (int scenario = 0; scenario < 4; ++scenario) {
            for (int i = 0; i <= n; ++i) {
                po[i] = 8.0e6 - 1.0e5*i/n;
            }

            for (int i = 0; i < n; ++i) {
                const double flow =
                    (scenario == 2) ? -100.0 : 100.0;

                mo[i] = flow + 0.1*i/n;
            }

            u[0] = bc.inlet_pressure_new
                 + 100.0*scenario;
            u[1] = (scenario == 2) ? -100.5 : 100.5;

            for (int i = 0; i < n; ++i) {
                u[2*i+2] = mo[i] + 0.25;
                u[2*i+3] = po[i+1]
                         + 200.0*(scenario+1);
            }

            assemble_residual(
                n, u, po, mo, bc, par, r_cpp);

            assemble_residual_(
                &n, u.data(), po.data(), mo.data(),
                &bc.inlet_flow_old,
                &bc.outlet_flow_old,
                &bc.outlet_flow_new,
                &bc.inlet_pressure_new,
                &par.dx, &par.dt, &par.theta,
                &par.diameter, &par.area,
                &par.temperature, &par.z,
                &par.gas_constant, &par.viscosity,
                &par.roughness, r_f77.data());

            for (std::size_t j = 0; j < r_cpp.size(); ++j) {
                ++comparisons;

                const double scale =
                    std::max(std::abs(r_cpp[j]),
                             std::abs(r_f77[j]));

                const double error =
                    std::abs(r_cpp[j]-r_f77[j]);

                const double relative =
                    scale > 0.0 ? error/scale : 0.0;

                max_relative_error =
                    std::max(max_relative_error, relative);

                if (!std::isfinite(r_cpp[j]) ||
                    !std::isfinite(r_f77[j]) ||
                    error > 5.0e-14*scale) {
                    ++failures;

                    std::cerr << std::setprecision(17)
                              << "FAIL n=" << n
                              << " scenario=" << scenario
                              << " component=" << j
                              << " cpp=" << r_cpp[j]
                              << " f77=" << r_f77[j]
                              << '\n';
                }
            }
        }
    }

    std::cout << std::scientific << std::setprecision(17)
              << "M12.2a residual comparisons: "
              << comparisons << '\n'
              << "Maximum relative error: "
              << max_relative_error << '\n';

    if (failures != 0) {
        std::cerr << "FAIL M12.2a: "
                  << failures << " discrepancies\n";
        return 1;
    }

    std::cout << "PASS M12.2a residual parity\n";
    return 0;
}
