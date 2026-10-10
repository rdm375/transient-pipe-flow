#pragma once

#include <span>

namespace pipe_sim {

// Solve A*x=b using the existing Fortran-compatible band algorithm.
//
// Matrix storage is column-major:
//   ab[j*ldab + (kl+ku+i-j)] = A(i,j)
//
// ab and x are overwritten; b is unchanged.
//
// Return value:
//    0: success
//   >0: failed pivot (1-based index)
//   -1: invalid n, kl, or ku
//   -2: invalid ldab
//
// Array dimensions are checked separately using exceptions.
int solve_banded(
    int n,
    std::span<double> ab,
    int ldab,
    int kl,
    int ku,
    std::span<const double> b,
    std::span<double> x);

} // namespace pipe_sim
