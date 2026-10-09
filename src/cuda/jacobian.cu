#include "pipe_sim/transient.hpp"
#include <cuda_runtime.h>
#include <cmath>

namespace pipe_sim {
__global__ void jacobian_kernel(int n, const double* u,
                                TransientParameters par, double* ab) {
    constexpr int ldab = 6;
    constexpr int kv = 3;
    const int tid = threadIdx.x;
    const int nu = 2*n + 2;
    for (int j = tid; j < ldab*nu; j += blockDim.x) ab[j] = 0.0;
    __syncthreads();

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
    const double crho = 1.0/(z*rs*t);

    if (tid == 0) {
        ab[0*ldab + kv + 0-0] = 1.0;
        ab[0*ldab + kv + 1-0] = a*dx*crho/(2.0*dt);
        ab[1*ldab + kv + 1-1] = -theta;
        ab[2*ldab + kv + 1-2] = theta;
    }
    if (tid >= n) return;
    const int i = tid;
    const int im = 2*i+2;
    const int il = (i==0) ? 0 : 2*i+1;
    const int ir = 2*i+3;
    const int row = 2*i+2;
    const double rhol = u[il]/(z*rs*t);
    const double rhor = u[ir]/(z*rs*t);
    const double rhof = 0.5*(rhol+rhor);
    const double re = fabs(u[im])*d/(a*mu);
    const double re09 = pow(re,0.9);
    const double arg = eps/(3.7*d)+5.74/re09;
    const double larg = log10(arg);
    const double fd = 0.25/(larg*larg);
    const double ln10 = log(10.0);
    const double dxdre = -0.9*5.74/(re*re09);
    const double dfdre = -0.5*dxdre/(arg*ln10*larg*larg*larg);
    const double dredm = copysign(d/(a*mu),u[im]);
    const double dfdm = dfdre*dredm;
    const double sf = fd*u[im]*fabs(u[im])/(2.0*d*a*a*rhof);
    const double c = 1.0/(2.0*d*a*a*rhof);
    const double dsdm = c*(dfdm*u[im]*fabs(u[im])+2.0*fd*fabs(u[im]));
    const double dsdp = -sf*crho/(2.0*rhof);

    ab[il*ldab + kv + row-il] = theta*(-1.0/dx+dsdp);
    ab[im*ldab + kv + row-im] = 1.0/(a*dt)+theta*dsdm;
    ab[ir*ldab + kv + row-ir] = theta*(1.0/dx+dsdp);
    if (i < n-1) {
        ab[ir*ldab + kv + row+1-ir] = a*dx*crho/dt;
        ab[im*ldab + kv + row+1-im] = -theta;
        ab[(im+2)*ldab + kv + row+1-(im+2)] = theta;
    } else {
        ab[ir*ldab + kv + row+1-ir] = a*dx*crho/(2.0*dt);
        ab[im*ldab + kv + row+1-im] = -theta;
    }
}

extern "C" cudaError_t launch_pipe_jacobian(int n, const double* u,
                                            TransientParameters par,
                                            double* ab) {
    if (n < 2 || n > 100 || !u || !ab) return cudaErrorInvalidValue;
    jacobian_kernel<<<1,128>>>(n,u,par,ab);
    return cudaGetLastError();
}
} // namespace pipe_sim
