module pipe_solver_api
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private

  integer, parameter, public :: PIPE_SUCCESS = 0
  integer, parameter, public :: PIPE_INVALID_ARGUMENT = 3

  public :: pipe_step_checked

contains

  subroutine pipe_step_checked(n, po, mo, mino, mouto, moutn, pinn, &
       dx, dt, theta, d, a, temp, z, rs, mu, eps, &
       rtol, stol, maxit, verbose, u, info, niter)

    integer, intent(in) :: n, maxit, verbose

    real(real64), intent(in) :: po(:), mo(:)
    real(real64), intent(in) :: mino, mouto, moutn, pinn
    real(real64), intent(in) :: dx, dt, theta
    real(real64), intent(in) :: d, a, temp, z, rs, mu, eps
    real(real64), intent(in) :: rtol, stol

    real(real64), intent(inout) :: u(:)

    integer, intent(out) :: info, niter

    info = PIPE_INVALID_ARGUMENT
    niter = 0

    if (n < 2 .or. n > 100) return

    if (size(po) < n+1) return
    if (size(mo) < n) return
    if (size(u) < 2*n+2) return

    call transient_step(n, po, mo, mino, mouto, moutn, pinn, &
         dx, dt, theta, d, a, temp, z, rs, mu, eps, &
         rtol, stol, maxit, verbose, u, info, niter)

  end subroutine pipe_step_checked

end module pipe_solver_api
