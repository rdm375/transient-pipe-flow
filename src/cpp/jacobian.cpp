#include "pipe_sim/transient.hpp"
#include "pipe_sim/physics.hpp"

#include <algorithm>
#include <stdexcept>

namespace pipe_sim {

void assemble_jacobian_banded(
    int n,
    std::span<const double> u,
    const TransientParameters& par,
    std::span<double> ab,
    int ldab)
{
    const int nu = 2*n + 2;
    constexpr int kv = 3;

    if (n < 2 || ldab < 6 ||
        u.size() != static_cast<std::size_t>(nu) ||
        ab.size() != static_cast<std::size_t>(ldab*nu)) {
        throw std::invalid_argument("Invalid Jacobian dimensions");
    }

    std::fill(ab.begin(), ab.end(), 0.0);

    const double dx = par.dx;
    const double dt = par.dt;
    const double theta = par.theta;
    const double d = par.diameter;
    const double a = par.area;
    const double t = par.temperature;
    const double z = par.z;
    const double rs = par.gas_constant;
    const double mu = par.viscosity;
    const double eps = par.roughness;

    const double crho = density_derivative_cz(t,z,rs);

    // Set coefficient at zero-based equation row and state column.
    const auto set = [&](int row, int col, double value) {
        const int band_row = kv + row - col;

        if (band_row < 0 || band_row >= ldab)
            throw std::logic_error("Jacobian coefficient outside band");

        ab[col*ldab + band_row] = value;
    };

    // Inlet pressure BC and inlet half-cell continuity.
    set(0,0,1.0);
    set(1,0,a*dx*crho/(2.0*dt));
    set(1,1,-theta);
    set(1,2,theta);

    for (int i = 0; i < n; ++i) {
        const int im = 2*i + 2;
        const int il = (i == 0) ? 0 : 2*i + 1;
        const int ir = 2*i + 3;
        const int row = 2*i + 2;

        const double rhol = density_cz(u[il],t,z,rs);
        const double rhor = density_cz(u[ir],t,z,rs);
        const double rhof = 0.5*(rhol+rhor);

        const double re = reynolds_mass(u[im],d,a,mu);

        const auto friction =
            friction_sj_value_derivative(re,eps,d);

        const double fd = friction.value;
        const double dfdre = friction.derivative_re;

        const double dredm =
            reynolds_derivative_mdot(u[im],d,a,mu);

        const double dfdm = dfdre*dredm;

        const double sf =
            friction_source(fd,d,a,rhof,u[im]);

        const double dsdm =
            friction_source_derivative_mdot(
                fd,dfdm,d,a,rhof,u[im]);

        const double dsdp =
            -sf*crho/(2.0*rhof);

        set(row,il,theta*(-1.0/dx+dsdp));
        set(row,im,1.0/(a*dt)+theta*dsdm);
        set(row,ir,theta*(1.0/dx+dsdp));

        if (i < n-1) {
            set(row+1,ir,a*dx*crho/dt);
            set(row+1,im,-theta);
            set(row+1,im+2,theta);
        } else {
            set(row+1,ir,a*dx*crho/(2.0*dt));
            set(row+1,im,-theta);
        }
    }
}

} // namespace pipe_sim
