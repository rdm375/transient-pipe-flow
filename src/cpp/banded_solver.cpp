#include "pipe_sim/banded_solver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace pipe_sim {

template<int KL, int KU, int LDAB>
int solve_banded_fixed(
    int n,
    std::span<double> ab,
    std::span<const double> b,
    std::span<double> x)
{
    constexpr int KV = KL + KU;

    const double pivtol =
        100.0 * std::numeric_limits<double>::epsilon();

    const auto get = [&](int i, int j) -> double {
        return ab[j*LDAB + KV + i - j];
    };

    const auto set = [&](int i, int j, double value) {
        ab[j*LDAB + KV + i - j] = value;
    };

    std::copy(b.begin(), b.end(), x.begin());

    for (int k = 0; k < n-1; ++k) {
        const int imax = std::min(n-1, k+KL);
        const int jmax = std::min(n-1, k+KV);

        int ip = k;
        double amax = std::abs(get(k,k));

        for (int i = k+1; i <= imax; ++i) {
            const double tmp = std::abs(get(i,k));

            if (tmp > amax) {
                amax = tmp;
                ip = i;
            }
        }

        if (amax <= pivtol)
            return k+1;

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
                const double tmp =
                    get(i,j) - factor*get(k,j);

                set(i,j,tmp);
            }

            x[i] = x[i] - factor*x[k];
        }
    }

    if (std::abs(get(n-1,n-1)) <= pivtol)
        return n;

    for (int i = n-1; i >= 0; --i) {
        double sum = x[i];

        for (int j = i+1; j <= std::min(n-1,i+KV); ++j)
            sum = sum - get(i,j)*x[j];

        x[i] = sum/get(i,i);
    }

    return 0;
}

int solve_banded(
    int n,
    std::span<double> ab,
    int ldab,
    int kl,
    int ku,
    std::span<const double> b,
    std::span<double> x)
{
    // Preserve the Fortran argument-validation order.
    if (n < 1 || kl < 0 || ku < 0)
        return -1;

    if (ldab < 2*kl + ku + 1)
        return -2;

    if (ab.size() != static_cast<std::size_t>(ldab*n) ||
        b.size() != static_cast<std::size_t>(n) ||
        x.size() != static_cast<std::size_t>(n)) {
        throw std::invalid_argument("Invalid banded solver array dimensions");
    }

    if (kl == 2 && ku == 1 && ldab == 6)
        return solve_banded_fixed<2,1,6>(n,ab,b,x);

    const int kv = kl + ku;
    const double pivtol =
        100.0 * std::numeric_limits<double>::epsilon();

    const auto get = [&](int i, int j) -> double {
        const int row = kv + i - j;

        if (row < 0 || row >= ldab)
            return 0.0;

        return ab[j*ldab + row];
    };

    const auto set = [&](int i, int j, double value) {
        const int row = kv + i - j;

        if (row >= 0 && row < ldab)
            ab[j*ldab + row] = value;
    };

    std::copy(b.begin(), b.end(), x.begin());

    for (int k = 0; k < n-1; ++k) {
        const int imax = std::min(n-1, k+kl);
        const int jmax = std::min(n-1, k+kv);

        // Pivot search.
        int ip = k;
        double amax = std::abs(get(k,k));

        for (int i = k+1; i <= imax; ++i) {
            const double tmp = std::abs(get(i,k));

            if (tmp > amax) {
                amax = tmp;
                ip = i;
            }
        }

        if (amax <= pivtol)
            return k+1;

        // Swap active row segments.
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

        // Eliminate entries below the pivot.
        for (int i = k+1; i <= imax; ++i) {
            const double factor = get(i,k)/get(k,k);

            set(i,k,0.0);

            for (int j = k+1; j <= jmax; ++j) {
                const double tmp =
                    get(i,j) - factor*get(k,j);

                set(i,j,tmp);
            }

            x[i] = x[i] - factor*x[k];
        }
    }

    if (std::abs(get(n-1,n-1)) <= pivtol)
        return n;

    // Back substitution.
    for (int i = n-1; i >= 0; --i) {
        double sum = x[i];

        for (int j = i+1; j <= std::min(n-1,i+kv); ++j) {
            sum = sum - get(i,j)*x[j];
        }

        x[i] = sum/get(i,i);
    }

    return 0;
}

} // namespace pipe_sim
