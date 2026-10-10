#include "pipe_sim/transient.hpp"

#include <cuda_runtime.h>
#include <cmath>

namespace pipe_sim {

__global__ void residual_kernel(
    int n,
    const double* u,
    const double* po,
    const double* mo,
    TransientBoundary bc,
    TransientParameters par,
    double* r)
{
    const int i = threadIdx.x;

    const double dx = par.dx;
    const double dt = par.dt;
    const double theta = par.theta;
    const double d = par.diameter;
    const double a = par.area;

    const double cold = 1.0 - theta;
    const double cnew = theta;

    if (i == 0) {
        r[0] = u[0] - bc.inlet_pressure_new;

        const double rho_new =
            u[0] / (par.z * par.gas_constant * par.temperature);

        const double rho_old =
            po[0] / (par.z * par.gas_constant * par.temperature);

        r[1] = a*dx/(2.0*dt)
             * (rho_new-rho_old)
             + cnew*(u[2]-u[1])
             + cold*(mo[0]-bc.inlet_flow_old);
    }

    if (i >= n)
        return;

    const int im = 2*i + 2;
    const int il = (i == 0) ? 0 : 2*i + 1;
    const int ir = 2*i + 3;
    const int row = 2*i + 2;

    const double eos =
        par.z * par.gas_constant * par.temperature;

    const double rhol = u[il] / eos;
    const double rhor = u[ir] / eos;
    const double rhof = 0.5*(rhol+rhor);

    const double re =
        fabs(u[im])*d/(a*par.viscosity);

    const double arg =
        par.roughness/(3.7*d)
        + 5.74/pow(re,0.9);

    const double larg = log10(arg);
    const double fd = 0.25/(larg*larg);

    const double sf =
        fd*u[im]*fabs(u[im])
        /(2.0*d*a*a*rhof);

    const double rhol0 = po[i] / eos;
    const double rhor0 = po[i+1] / eos;
    const double rhof0 = 0.5*(rhol0+rhor0);

    const double re0 =
        fabs(mo[i])*d/(a*par.viscosity);

    const double arg0 =
        par.roughness/(3.7*d)
        + 5.74/pow(re0,0.9);

    const double larg0 = log10(arg0);
    const double fd0 = 0.25/(larg0*larg0);

    const double sf0 =
        fd0*mo[i]*fabs(mo[i])
        /(2.0*d*a*a*rhof0);

    r[row] = (u[im]-mo[i])/(a*dt)
           + cnew*((u[ir]-u[il])/dx+sf)
           + cold*((po[i+1]-po[i])/dx+sf0);

    if (i < n-1) {
        r[row+1] = a*dx/dt
                 * (u[ir]/eos-po[i+1]/eos)
                 + cnew*(u[im+2]-u[im])
                 + cold*(mo[i+1]-mo[i]);
    } else {
        r[row+1] = a*dx/(2.0*dt)
                 * (u[ir]/eos-po[n]/eos)
                 + cnew*(bc.outlet_flow_new-u[im])
                 + cold*(bc.outlet_flow_old-mo[i]);
    }
}

extern "C" cudaError_t launch_pipe_residual(
    int n,
    const double* u,
    const double* po,
    const double* mo,
    TransientBoundary bc,
    TransientParameters par,
    double* r)
{
    if (n < 2 || n > 100 ||
        !u || !po || !mo || !r)
        return cudaErrorInvalidValue;

    residual_kernel<<<1,128>>>(
        n,u,po,mo,bc,par,r);

    return cudaGetLastError();
}

} // namespace pipe_sim
