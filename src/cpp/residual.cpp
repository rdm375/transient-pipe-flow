#include "pipe_sim/transient.hpp"
#include "pipe_sim/physics.hpp"

#include <stdexcept>

namespace pipe_sim {

void assemble_residual(
    int n,
    std::span<const double> u,
    std::span<const double> po,
    std::span<const double> mo,
    const TransientBoundary& bc,
    const TransientParameters& par,
    std::span<double> r)
{
    if (n < 2 ||
        u.size() != static_cast<std::size_t>(2*n+2) ||
        po.size() != static_cast<std::size_t>(n+1) ||
        mo.size() != static_cast<std::size_t>(n) ||
        r.size() != static_cast<std::size_t>(2*n+2)) {
        throw std::invalid_argument("Invalid transient array dimensions");
    }

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

    const double cold = 1.0 - theta;
    const double cnew = theta;

    const auto rho = [&](double pressure) {
        return density_cz(pressure, t, z, rs);
    };

    // Inlet pressure boundary condition.
    r[0] = u[0] - bc.inlet_pressure_new;

    // Inlet half control volume.
    r[1] = a*dx/(2.0*dt)
         * (rho(u[0]) - rho(po[0]))
         + cnew*(u[2]-u[1])
         + cold*(mo[0]-bc.inlet_flow_old);

    for (int i = 0; i < n; ++i) {
        const int im = 2*i + 2;
        const int il = (i == 0) ? 0 : 2*i + 1;
        const int ir = 2*i + 3;
        const int row = 2*i + 2;

        // New-time friction source.
        const double rhol = rho(u[il]);
        const double rhor = rho(u[ir]);
        const double rhof = 0.5*(rhol+rhor);

        const double re = reynolds_mass(u[im],d,a,mu);
        const double fd = friction_sj(re,eps,d);
        const double sf = friction_source(fd,d,a,rhof,u[im]);

        // Old-time friction source.
        const double rhol0 = rho(po[i]);
        const double rhor0 = rho(po[i+1]);
        const double rhof0 = 0.5*(rhol0+rhor0);

        const double re0 = reynolds_mass(mo[i],d,a,mu);
        const double fd0 = friction_sj(re0,eps,d);
        const double sf0 = friction_source(fd0,d,a,rhof0,mo[i]);

        // Reduced momentum equation.
        r[row] = (u[im]-mo[i])/(a*dt)
               + cnew*((u[ir]-u[il])/dx+sf)
               + cold*((po[i+1]-po[i])/dx+sf0);

        // Continuity to the right of this face.
        if (i < n-1) {
            r[row+1] = a*dx/dt
                     * (rho(u[ir])-rho(po[i+1]))
                     + cnew*(u[im+2]-u[im])
                     + cold*(mo[i+1]-mo[i]);
        } else {
            // Outlet half control volume.
            r[row+1] = a*dx/(2.0*dt)
                     * (rho(u[ir])-rho(po[n]))
                     + cnew*(bc.outlet_flow_new-u[im])
                     + cold*(bc.outlet_flow_old-mo[i]);
        }
    }
}

} // namespace pipe_sim
