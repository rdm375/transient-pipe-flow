#include "pipe_sim/newton.hpp"
#include "pipe_sim/integrate.hpp"

#include <cuda_runtime.h>
#include <cmath>

namespace pipe_sim {
namespace {

constexpr int MAX_N = 100;
constexpr int MAX_NU = 2*MAX_N+2;
constexpr int LDAB = 6;
constexpr int KV = 3;

__device__ double density(double p, const TransientParameters& par)
{
    return p/(par.z*par.gas_constant*par.temperature);
}

__device__ double friction(double m, const TransientParameters& par)
{
    const double re =
        fabs(m)*par.diameter/(par.area*par.viscosity);

    const double arg =
        par.roughness/(3.7*par.diameter)+5.74/pow(re,0.9);

    const double larg = log10(arg);
    return 0.25/(larg*larg);
}

__device__ double friction_source(
    double fd, double m, double rhof,
    const TransientParameters& par)
{
    return fd*m*fabs(m)/
        (2.0*par.diameter*par.area*par.area*rhof);
}

__device__ void residual(
    int n, const double* u, const double* po, const double* mo,
    TransientBoundary bc, TransientParameters par,
    const double* old_sf, const double* old_fd,
    bool initial, double* r)
{
    const double dx = par.dx;
    const double dt = par.dt;
    const double a = par.area;
    const double theta = par.theta;
    const double cold = 1.0-theta;
    const double cnew = theta;

    r[0] = u[0]-bc.inlet_pressure_new;

    r[1] = a*dx/(2.0*dt)*
        (density(u[0],par)-density(po[0],par))
        +cnew*(u[2]-u[1])
        +cold*(mo[0]-bc.inlet_flow_old);

    for (int i=0; i<n; ++i) {
        const int im=2*i+2;
        const int il=(i==0)?0:2*i+1;
        const int ir=2*i+3;
        const int row=2*i+2;

        const double rhol=density(u[il],par);
        const double rhor=density(u[ir],par);
        const double rhof=0.5*(rhol+rhor);

        const double fd =
            (initial && u[im]==mo[i])
            ? old_fd[i] : friction(u[im],par);

        const double sf=friction_source(fd,u[im],rhof,par);

        r[row]=(u[im]-mo[i])/(a*dt)
            +cnew*((u[ir]-u[il])/dx+sf)
            +cold*((po[i+1]-po[i])/dx+old_sf[i]);

        if (i<n-1) {
            r[row+1]=a*dx/dt*
                (density(u[ir],par)-density(po[i+1],par))
                +cnew*(u[im+2]-u[im])
                +cold*(mo[i+1]-mo[i]);
        } else {
            r[row+1]=a*dx/(2.0*dt)*
                (density(u[ir],par)-density(po[n],par))
                +cnew*(bc.outlet_flow_new-u[im])
                +cold*(bc.outlet_flow_old-mo[i]);
        }
    }
}

__device__ void jacobian(
    int n, const double* u, TransientParameters par, double* ab)
{
    const int nu=2*n+2;
    for (int j=0; j<LDAB*nu; ++j) ab[j]=0.0;

    const double dx=par.dx;
    const double dt=par.dt;
    const double theta=par.theta;
    const double d=par.diameter;
    const double a=par.area;
    const double mu=par.viscosity;
    const double eps=par.roughness;
    const double crho=
        1.0/(par.z*par.gas_constant*par.temperature);

    ab[0*LDAB+KV+0-0]=1.0;
    ab[0*LDAB+KV+1-0]=a*dx*crho/(2.0*dt);
    ab[1*LDAB+KV+1-1]=-theta;
    ab[2*LDAB+KV+1-2]=theta;

    for (int i=0; i<n; ++i) {
        const int im=2*i+2;
        const int il=(i==0)?0:2*i+1;
        const int ir=2*i+3;
        const int row=2*i+2;

        const double rhol=density(u[il],par);
        const double rhor=density(u[ir],par);
        const double rhof=0.5*(rhol+rhor);
        const double re=fabs(u[im])*d/(a*mu);
        const double re09=pow(re,0.9);
        const double arg=eps/(3.7*d)+5.74/re09;
        const double larg=log10(arg);
        const double fd=0.25/(larg*larg);
        const double ln10=log(10.0);
        const double dxdre=-0.9*5.74/(re*re09);
        const double dfdre=
            -0.5*dxdre/(arg*ln10*larg*larg*larg);
        const double dredm=copysign(d/(a*mu),u[im]);
        const double dfdm=dfdre*dredm;
        const double sf=friction_source(fd,u[im],rhof,par);
        const double c=1.0/(2.0*d*a*a*rhof);
        const double dsdm=c*
            (dfdm*u[im]*fabs(u[im])+2.0*fd*fabs(u[im]));
        const double dsdp=-sf*crho/(2.0*rhof);

        ab[il*LDAB+KV+row-il]=theta*(-1.0/dx+dsdp);
        ab[im*LDAB+KV+row-im]=1.0/(a*dt)+theta*dsdm;
        ab[ir*LDAB+KV+row-ir]=theta*(1.0/dx+dsdp);

        if (i<n-1) {
            ab[ir*LDAB+KV+row+1-ir]=a*dx*crho/dt;
            ab[im*LDAB+KV+row+1-im]=-theta;
            ab[(im+2)*LDAB+KV+row+1-(im+2)]=theta;
        } else {
            ab[ir*LDAB+KV+row+1-ir]=a*dx*crho/(2.0*dt);
            ab[im*LDAB+KV+row+1-im]=-theta;
        }
    }
}

__device__ int banded_solve(
    int n, double* ab, const double* b, double* x)
{
    constexpr int KL=2;
    constexpr double pivtol=
        100.0*2.220446049250313080847263336181640625e-16;

    for (int i=0; i<n; ++i) x[i]=b[i];

    for (int k=0; k<n-1; ++k) {
        const int imax=(k+KL<n)?k+KL:n-1;
        const int jmax=(k+KV<n)?k+KV:n-1;

        int ip=k;
        double amax=fabs(ab[k*LDAB+KV]);

        for (int i=k+1; i<=imax; ++i) {
            const double candidate=
                fabs(ab[k*LDAB+KV+i-k]);
            if (candidate>amax) {
                amax=candidate;
                ip=i;
            }
        }

        if (amax<=pivtol) return k+1;

        if (ip!=k) {
            for (int j=k; j<=jmax; ++j) {
                const int ik=j*LDAB+KV+k-j;
                const int ii=j*LDAB+KV+ip-j;
                const double tmp=ab[ik];
                ab[ik]=ab[ii];
                ab[ii]=tmp;
            }
            const double tmp=x[k];
            x[k]=x[ip];
            x[ip]=tmp;
        }

        for (int i=k+1; i<=imax; ++i) {
            const int ik=k*LDAB+KV+i-k;
            const double factor=ab[ik]/ab[k*LDAB+KV];
            ab[ik]=0.0;

            for (int j=k+1; j<=jmax; ++j) {
                const int ij=j*LDAB+KV+i-j;
                const int kj=j*LDAB+KV+k-j;
                ab[ij]=ab[ij]-factor*ab[kj];
            }
            x[i]=x[i]-factor*x[k];
        }
    }

    if (fabs(ab[(n-1)*LDAB+KV])<=pivtol) return n;

    for (int i=n-1; i>=0; --i) {
        double sum=x[i];
        const int jmax=(i+KV<n)?i+KV:n-1;
        for (int j=i+1; j<=jmax; ++j)
            sum-=ab[j*LDAB+KV+i-j]*x[j];
        x[i]=sum/ab[i*LDAB+KV];
    }
    return 0;
}

__device__ bool valid(
    int n, const double* state, const double* po,
    const double* mo, TransientBoundary bc,
    TransientParameters par, NewtonOptions opt)
{
    if (n<2 || n>MAX_N) return false;

    if (!isfinite(par.dx) || par.dx<=0.0 ||
        !isfinite(par.dt) || par.dt<=0.0 ||
        !isfinite(par.theta) || par.theta<0.0 ||
        par.theta>1.0 ||
        !isfinite(par.diameter) || par.diameter<=0.0 ||
        !isfinite(par.area) || par.area<=0.0 ||
        !isfinite(par.temperature) || par.temperature<=0.0 ||
        !isfinite(par.z) || par.z<=0.0 ||
        !isfinite(par.gas_constant) || par.gas_constant<=0.0 ||
        !isfinite(par.viscosity) || par.viscosity<=0.0 ||
        !isfinite(par.roughness) || par.roughness<0.0 ||
        !isfinite(opt.residual_tolerance) ||
        !isfinite(opt.step_tolerance) ||
        opt.residual_tolerance<=0.0 ||
        opt.step_tolerance<0.0 ||
        opt.max_iterations<0) return false;

    if (!isfinite(bc.inlet_flow_old) ||
        !isfinite(bc.outlet_flow_old) ||
        !isfinite(bc.outlet_flow_new) ||
        !isfinite(bc.inlet_pressure_new) ||
        bc.inlet_pressure_new<=0.0) return false;

    for (int i=0; i<2*n+2; ++i)
        if (!isfinite(state[i])) return false;

    for (int i=0; i<n+1; ++i)
        if (!isfinite(po[i]) || po[i]<=0.0) return false;

    for (int i=0; i<n; ++i)
        if (!isfinite(mo[i])) return false;

    return true;
}

__device__ NewtonResult newton_device(
    int n, double* state, const double* po, const double* mo,
    TransientBoundary bc, TransientParameters par,
    NewtonOptions opt)
{
    if (!valid(n,state,po,mo,bc,par,opt)) {
        return {3,0};
    }

    const int nu=2*n+2;

    // Device-local workspace; allocated once per kernel invocation.
    double r[MAX_NU];
    double trial_r[MAX_NU];
    double ab[LDAB*MAX_NU];
    double rhs[MAX_NU];
    double correction[MAX_NU];
    double trial[MAX_NU];
    double work[MAX_NU];
    double old_sf[MAX_N];
    double old_fd[MAX_N];

    for (int i=0; i<n; ++i) {
        const double rhol=density(po[i],par);
        const double rhor=density(po[i+1],par);
        const double rhof=0.5*(rhol+rhor);
        const double fd=friction(mo[i],par);
        old_fd[i]=fd;
        old_sf[i]=friction_source(fd,mo[i],rhof,par);
    }

    for (int i=0; i<nu; ++i) work[i]=state[i];

    residual(n,work,po,mo,bc,par,old_sf,old_fd,true,r);

    for (int iter=0; iter<opt.max_iterations; ++iter) {
        double rnorm=0.0;
        for (int i=0; i<nu; ++i)
            rnorm=fmax(rnorm,fabs(r[i]));

        if (rnorm<=opt.residual_tolerance) {
            for (int i=0; i<nu; ++i) state[i]=work[i];
            return {0,iter};
        }

        jacobian(n,work,par,ab);

        for (int i=0; i<nu; ++i) rhs[i]=-r[i];

        const int linfo=banded_solve(nu,ab,rhs,correction);
        if (linfo!=0) {
            return {2,iter+1};
        }

        double snorm=0.0;
        double uscale=1.0;
        for (int i=0; i<nu; ++i) {
            snorm=fmax(snorm,fabs(correction[i]));
            uscale=fmax(uscale,fabs(work[i]));
        }

        double lambda=1.0;
        double rnew=0.0;
        bool accepted=false;

        for (int ls=0; ls<20; ++ls) {
            for (int i=0; i<nu; ++i)
                trial[i]=work[i]+lambda*correction[i];

            residual(n,trial,po,mo,bc,par,
                     old_sf,old_fd,false,trial_r);

            rnew=0.0;
            for (int i=0; i<nu; ++i)
                rnew=fmax(rnew,fabs(trial_r[i]));

            if (rnew<rnorm) {
                accepted=true;
                break;
            }
            lambda=0.5*lambda;
        }

        if (!accepted) {
            return {4,iter+1};
        }

        for (int i=0; i<nu; ++i) {
            work[i]=trial[i];
            r[i]=trial_r[i];
        }

        if (lambda*snorm<=opt.step_tolerance*uscale &&
            rnew<=10.0*opt.residual_tolerance) {
            for (int i=0; i<nu; ++i) state[i]=work[i];
            return {0,iter+1};
        }
    }

    return {1,opt.max_iterations};
}

__global__ void newton_kernel(
    int n, double* state, const double* po, const double* mo,
    TransientBoundary bc, TransientParameters par,
    NewtonOptions opt, NewtonResult* result)
{
    if (blockIdx.x!=0 || threadIdx.x!=0) return;
    *result = newton_device(n,state,po,mo,bc,par,opt);
}

} // namespace

extern "C" cudaError_t launch_pipe_newton(
    int n, double* state, const double* old_pressure,
    const double* old_flow, TransientBoundary bc,
    TransientParameters par, NewtonOptions opt,
    NewtonResult* result)
{
    if (n<2 || n>100 || !state || !old_pressure ||
        !old_flow || !result)
        return cudaErrorInvalidValue;

    newton_kernel<<<1,1>>>(
        n,state,old_pressure,old_flow,bc,par,opt,result);

    return cudaGetLastError();
}


namespace {

__device__ double gpu_demand_ramp(
    double time, DemandRamp demand)
{
    if (demand.ramp_duration <= 0.0 ||
        time >= demand.ramp_duration)
        return demand.final_outlet_flow;

    if (time <= 0.0)
        return demand.initial_outlet_flow;

    const double fraction = time / demand.ramp_duration;

    return demand.initial_outlet_flow +
        fraction * (demand.final_outlet_flow -
                    demand.initial_outlet_flow);
}

__device__ double gpu_linepack(
    int n, const double* pressure,
    TransientParameters par)
{
    double mass = 0.0;

    for (int i = 0; i <= n; ++i) {
        double weight = 1.0;

        if (i == 0 || i == n)
            weight = 0.5;

        mass = mass +
            weight * density(pressure[i], par);
    }

    return par.area * par.dx * mass;
}

__global__ void persistent_integrate_kernel(
    int n,
    double* pressure,
    double* flow,
    double* inlet_flow,
    DemandRamp demand,
    int steps,
    TransientParameters par,
    NewtonOptions opt,
    double* time,
    double* history_inlet,
    double* history_outlet,
    double* history_pressure,
    double* history_linepack,
    int* history_iterations,
    IntegrationResult* result)
{
    if (blockIdx.x != 0 || threadIdx.x != 0)
        return;

    if (n < 2 || n > 100 || steps < 0) {
        *result = {3, 0, 0.0};
        return;
    }

    if (!isfinite(*inlet_flow) ||
        !isfinite(demand.inlet_pressure) ||
        demand.inlet_pressure <= 0.0 ||
        !isfinite(demand.initial_outlet_flow) ||
        !isfinite(demand.final_outlet_flow) ||
        !isfinite(demand.ramp_duration) ||
        demand.ramp_duration < 0.0) {
        *result = {3, 0, 0.0};
        return;
    }

    // Reuse Newton's complete physical and numerical validation.
    // A constructed state is sufficient for the initial validation.
    double state[MAX_NU];

    state[0] = demand.inlet_pressure;
    state[1] = *inlet_flow;

    for (int i = 0; i < n; ++i) {
        state[2*i+2] = flow[i];
        state[2*i+3] = pressure[i+1];
    }

    const TransientBoundary initial_bc{
        *inlet_flow,
        demand.initial_outlet_flow,
        demand.initial_outlet_flow,
        demand.inlet_pressure
    };

    if (!valid(n,state,pressure,flow,
               initial_bc,par,opt)) {
        *result = {3, 0, 0.0};
        return;
    }

    double old_outlet_flow =
        demand.initial_outlet_flow;

    time[0] = 0.0;
    history_inlet[0] = *inlet_flow;
    history_outlet[0] = old_outlet_flow;
    history_pressure[0] = pressure[n];
    history_linepack[0] =
        gpu_linepack(n,pressure,par);
    history_iterations[0] = 0;

    double max_mass_defect = 0.0;

    for (int k = 1; k <= steps; ++k) {
        const double new_time =
            static_cast<double>(k) * par.dt;

        const double new_outlet_flow =
            gpu_demand_ramp(new_time,demand);

        const TransientBoundary bc{
            *inlet_flow,
            old_outlet_flow,
            new_outlet_flow,
            demand.inlet_pressure
        };

        // Match the CPU transient_step initial guess.
        state[0] = bc.inlet_pressure_new;
        state[1] = bc.inlet_flow_old;

        for (int i = 0; i < n; ++i) {
            state[2*i+2] = flow[i];
            state[2*i+3] = pressure[i+1];
        }

        const NewtonResult step = newton_device(
            n,state,pressure,flow,bc,par,opt);

        if (step.info != 0) {
            *result = {
                step.info,k-1,max_mass_defect
            };
            return;
        }

        // Newton succeeded: accept the new timestep.
        pressure[0] = state[0];

        for (int i = 0; i < n; ++i) {
            flow[i] = state[2*i+2];
            pressure[i+1] = state[2*i+3];
        }

        const double new_inlet_flow = state[1];

        time[k] = new_time;
        history_inlet[k] = new_inlet_flow;
        history_outlet[k] = new_outlet_flow;
        history_pressure[k] = pressure[n];

        history_linepack[k] =
            gpu_linepack(n,pressure,par);

        history_iterations[k] = step.iterations;

        const double balance_error = fabs(
            (history_linepack[k] -
             history_linepack[k-1]) / par.dt
            - (par.theta *
                (history_inlet[k] -
                 history_outlet[k])
              + (1.0-par.theta) *
                (history_inlet[k-1] -
                 history_outlet[k-1])));

        max_mass_defect =
            fmax(max_mass_defect,balance_error);

        *inlet_flow = new_inlet_flow;
        old_outlet_flow = new_outlet_flow;
    }

    *result = {0,steps,max_mass_defect};
}

} // namespace

extern "C" cudaError_t launch_pipe_integrate(
    int n,
    double* pressure,
    double* flow,
    double* inlet_flow,
    DemandRamp demand,
    int steps,
    TransientParameters par,
    NewtonOptions opt,
    double* time,
    double* history_inlet,
    double* history_outlet,
    double* history_pressure,
    double* history_linepack,
    int* history_iterations,
    IntegrationResult* result)
{
    if (n < 2 || n > 100 || steps < 0 ||
        !pressure || !flow || !inlet_flow ||
        !time || !history_inlet || !history_outlet ||
        !history_pressure || !history_linepack ||
        !history_iterations || !result)
        return cudaErrorInvalidValue;

    persistent_integrate_kernel<<<1,1>>>(
        n,pressure,flow,inlet_flow,demand,steps,
        par,opt,time,history_inlet,history_outlet,
        history_pressure,history_linepack,
        history_iterations,result);

    return cudaGetLastError();
}

} // namespace pipe_sim
