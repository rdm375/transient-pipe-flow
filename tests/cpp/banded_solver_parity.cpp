#include "pipe_sim/banded_solver.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

extern "C" {
void solve_banded_(
    const int* n,
    double* ab,
    const int* ldab,
    const int* kl,
    const int* ku,
    const double* b,
    double* x,
    int* info);
}

namespace {

int comparisons = 0;
int failures = 0;
double max_relative_error = 0.0;

void compare(
    double cpp,
    double f77,
    const std::string& label)
{
    ++comparisons;

    const double scale =
        std::max(std::abs(cpp),std::abs(f77));

    const double error = std::abs(cpp-f77);

    const double relative =
        scale > 0.0 ? error/scale : 0.0;

    max_relative_error =
        std::max(max_relative_error,relative);

    if (!std::isfinite(cpp) ||
        !std::isfinite(f77) ||
        error > 5.0e-14*scale) {

        ++failures;

        std::cerr << std::setprecision(17)
                  << "FAIL " << label
                  << " cpp=" << cpp
                  << " f77=" << f77 << '\n';
    }
}

void run_case(int n, int scenario)
{
    constexpr int kl = 2;
    constexpr int ku = 1;
    constexpr int ldab = 6;
    constexpr int kv = kl+ku;

    std::vector<double> ab(ldab*n,0.0);
    std::vector<double> b(n);
    std::vector<double> x_cpp(n,0.0);
    std::vector<double> x_f77(n,0.0);

    const auto set = [&](int i,int j,double value) {
        const int row = kv+i-j;

        if (row >= 0 && row < ldab)
            ab[j*ldab+row] = value;
    };

    // Deterministic banded matrices.
    for (int i = 0; i < n; ++i) {
        set(i,i,5.0+0.01*i);

        if (i+1 < n) {
            set(i+1,i,-0.5);
            set(i,i+1,0.25);
        }

        if (i+2 < n)
            set(i+2,i,0.125);

        b[i] = 1.0+0.1*i;
    }

    if (scenario == 1) {
        // Force a row interchange at the first pivot.
        set(0,0,0.01);
        set(1,0,2.0);
    }

    if (scenario == 2) {
        // Singular first column.
        set(0,0,0.0);
        set(1,0,0.0);

        if (n > 2)
            set(2,0,0.0);
    }

    if (scenario == 3) {
        // Singular final pivot.
        // Zero the final matrix row's accessible entries.
        set(n-1,n-1,0.0);

        if (n >= 2)
            set(n-1,n-2,0.0);

        if (n >= 3)
            set(n-1,n-3,0.0);
    }

    std::vector<double> cpp = ab;
    std::vector<double> f77 = ab;

    int info_f77 = -999;

    const int info_cpp = pipe_sim::solve_banded(
        n,cpp,ldab,kl,ku,b,x_cpp);

    solve_banded_(
        &n,f77.data(),&ldab,&kl,&ku,
        b.data(),x_f77.data(),&info_f77);

    ++comparisons;

    if (info_cpp != info_f77) {
        ++failures;

        std::cerr << "FAIL INFO"
                  << " n=" << n
                  << " scenario=" << scenario
                  << " cpp=" << info_cpp
                  << " f77=" << info_f77 << '\n';
    }

    for (std::size_t j = 0; j < cpp.size(); ++j) {
        compare(
            cpp[j],f77[j],
            "matrix n="+std::to_string(n)
            +" scenario="+std::to_string(scenario)
            +" index="+std::to_string(j));
    }

    // On successful solves, compare the final solution.
    // On pivot failure, compare the partially transformed RHS.
    if (info_cpp == info_f77) {
        for (int i = 0; i < n; ++i) {
            compare(
                x_cpp[i],x_f77[i],
                "solution n="+std::to_string(n)
                +" scenario="+std::to_string(scenario)
                +" index="+std::to_string(i));
        }
    }
}

void run_invalid_arguments()
{
    constexpr int n = 3;
    constexpr int ldab = 6;
    constexpr int kl_value = 2;
    constexpr int ku_value = 1;

    std::vector<double> ab(ldab*n,0.0);
    std::vector<double> b(n,1.0);
    std::vector<double> x_cpp(n,0.0);
    std::vector<double> x_f77(n,0.0);

    for (int scenario = 0; scenario < 2; ++scenario) {
        const int test_n = scenario == 0 ? 0 : n;
        const int test_ldab = scenario == 0 ? ldab : 5;

        int info_f77 = -999;

        const int info_cpp = pipe_sim::solve_banded(
            test_n,ab,test_ldab,2,1,b,x_cpp);

        solve_banded_(
            &test_n,ab.data(),&test_ldab,
            // The remaining integer arguments are KL and KU.
            // Both remain valid in these tests.
            &kl_value,&ku_value,
            b.data(),x_f77.data(),&info_f77);

        ++comparisons;

        if (info_cpp != info_f77) {
            ++failures;

            std::cerr << "FAIL invalid argument"
                      << " scenario=" << scenario
                      << " cpp=" << info_cpp
                      << " f77=" << info_f77 << '\n';
        }
    }
}

} // namespace

int main()
{
    for (int n : {2,3,10,25,100}) {
        for (int scenario = 0; scenario < 4; ++scenario)
            run_case(n,scenario);
    }

    run_invalid_arguments();

    std::cout << std::scientific << std::setprecision(17)
              << "M12.3a solver comparisons: "
              << comparisons << '\n'
              << "Maximum relative error: "
              << max_relative_error << '\n';

    if (failures != 0) {
        std::cerr << "FAIL M12.3a: "
                  << failures << " discrepancies\n";
        return 1;
    }

    std::cout << "PASS M12.3a banded solver parity\n";
    return 0;
}
