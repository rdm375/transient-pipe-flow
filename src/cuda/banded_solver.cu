#include <cuda_runtime.h>
#include <cmath>

namespace pipe_sim {

__global__ void banded_solve_kernel(
    int n, double* ab, const double* b, double* x, int* info)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) return;

    constexpr int KL = 2;
    constexpr int KV = 3;
    constexpr int LDAB = 6;
    constexpr double pivtol =
        100.0 * 2.220446049250313080847263336181640625e-16;

    auto get = [&](int i, int j) -> double {
        return ab[j*LDAB + KV + i-j];
    };

    auto set = [&](int i, int j, double value) {
        ab[j*LDAB + KV + i-j] = value;
    };

    for (int i = 0; i < n; ++i)
        x[i] = b[i];

    for (int k = 0; k < n-1; ++k) {
        const int imax = (k+KL < n) ? k+KL : n-1;
        const int jmax = (k+KV < n) ? k+KV : n-1;

        int ip = k;
        double amax = fabs(get(k,k));

        for (int i = k+1; i <= imax; ++i) {
            const double candidate = fabs(get(i,k));

            if (candidate > amax) {
                amax = candidate;
                ip = i;
            }
        }

        if (amax <= pivtol) {
            *info = k+1;
            return;
        }

        if (ip != k) {
            for (int j = k; j <= jmax; ++j) {
                const double tmp = get(k,j);
                set(k,j,get(ip,j));
                set(ip,j,tmp);
            }

            const double tmp = x[k];
            x[k] = x[ip];
            x[ip] = tmp;
        }

        for (int i = k+1; i <= imax; ++i) {
            const double factor = get(i,k)/get(k,k);

            set(i,k,0.0);

            for (int j = k+1; j <= jmax; ++j) {
                const double updated =
                    get(i,j)-factor*get(k,j);

                set(i,j,updated);
            }

            x[i] = x[i]-factor*x[k];
        }
    }

    if (fabs(get(n-1,n-1)) <= pivtol) {
        *info = n;
        return;
    }

    for (int i = n-1; i >= 0; --i) {
        double sum = x[i];

        const int jmax = (i+KV < n) ? i+KV : n-1;

        for (int j = i+1; j <= jmax; ++j)
            sum -= get(i,j)*x[j];

        x[i] = sum/get(i,i);
    }

    *info = 0;
}

extern "C" cudaError_t launch_pipe_banded_solve(
    int n, double* ab, const double* b, double* x, int* info)
{
    if (n < 1 || n > 202 ||
        !ab || !b || !x || !info)
        return cudaErrorInvalidValue;

    banded_solve_kernel<<<1,1>>>(n,ab,b,x,info);

    return cudaGetLastError();
}

} // namespace pipe_sim
