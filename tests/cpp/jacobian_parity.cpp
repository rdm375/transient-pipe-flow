#include "pipe_sim/transient.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

extern "C" {
void assemble_jacobian_banded_(
    const int* n,
    const double* u,
    const double* dx,
    const double* dt,
    const double* theta,
    const double* diameter,
    const double* area,
    const double* temperature,
    const double* z,
    const double* rs,
    const double* viscosity,
    const double* roughness,
    double* ab,
    const int* ldab);
}

int main()
{
    using namespace pipe_sim;

    int comparisons = 0;
    int failures = 0;
    double max_relative_error = 0.0;

    constexpr int ldab = 6;

    for (int n : {2, 3, 10, 25, 100}) {
        const TransientParameters par{
            100000.0/n,
            15.0,
            0.65,
            1.0,
            std::acos(-1.0)/4.0,
            288.15,
            0.9,
            500.0,
            1.1e-5,
            4.5e-5
        };

        const int nu = 2*n+2;

        std::vector<double> u(nu);
        std::vector<double> cpp(ldab*nu);
        std::vector<double> f77(ldab*nu);

        for (int scenario = 0; scenario < 4; ++scenario) {
            u[0] = 8.0e6 + 100.0*scenario;
            u[1] = (scenario == 2) ? -100.0 : 100.0;

            for (int i = 0; i < n; ++i) {
                u[2*i+2] =
                    ((scenario == 2) ? -100.0 : 100.0)
                    + 0.25*(i+1);

                u[2*i+3] =
                    8.0e6
                    - 1.0e5*(i+1)/n
                    + 200.0*scenario;
            }

            assemble_jacobian_banded(
                n,u,par,cpp,ldab);

            assemble_jacobian_banded_(
                &n,u.data(),
                &par.dx,&par.dt,&par.theta,
                &par.diameter,&par.area,
                &par.temperature,&par.z,
                &par.gas_constant,&par.viscosity,
                &par.roughness,
                f77.data(),&ldab);

            for (std::size_t j = 0; j < cpp.size(); ++j) {
                ++comparisons;

                const double scale =
                    std::max(std::abs(cpp[j]),std::abs(f77[j]));

                const double error =
                    std::abs(cpp[j]-f77[j]);

                const double relative =
                    scale > 0.0 ? error/scale : 0.0;

                max_relative_error =
                    std::max(max_relative_error,relative);

                if (!std::isfinite(cpp[j]) ||
                    !std::isfinite(f77[j]) ||
                    error > 5.0e-14*scale) {
                    ++failures;

                    std::cerr << std::setprecision(17)
                              << "FAIL n=" << n
                              << " scenario=" << scenario
                              << " storage_index=" << j
                              << " cpp=" << cpp[j]
                              << " f77=" << f77[j]
                              << '\n';
                }
            }
        }
    }

    std::cout << std::scientific << std::setprecision(17)
              << "M12.2b Jacobian comparisons: "
              << comparisons << '\n'
              << "Maximum relative error: "
              << max_relative_error << '\n';

    if (failures != 0) {
        std::cerr << "FAIL M12.2b: "
                  << failures << " discrepancies\n";
        return 1;
    }

    std::cout << "PASS M12.2b Jacobian parity\n";
    return 0;
}
