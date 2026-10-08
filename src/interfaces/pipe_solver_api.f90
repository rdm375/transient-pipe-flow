module pipe_solver_api
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private

  integer, parameter, public :: PIPE_SUCCESS = 0
  integer, parameter, public :: PIPE_INVALID_ARGUMENT = 3

  public :: pipe_step_checked
  public :: pipe_integrate_checked
  public :: pipe_adaptive_checked

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


  subroutine pipe_integrate_checked(n, po, mo, mino, pinn, &
       mdot0, mdot1, tramp, dt, nsteps, theta, dx, d, a, &
       temp, z, rs, mu, eps, rtol, stol, maxit, &
       time, hmin, hout, hpout, hline, hnit, hmaxbal, info)

    integer, intent(in) :: n, nsteps, maxit

    real(real64), intent(inout) :: po(:), mo(:), mino

    real(real64), intent(in) :: pinn, mdot0, mdot1, tramp
    real(real64), intent(in) :: dt, theta, dx, d, a
    real(real64), intent(in) :: temp, z, rs, mu, eps
    real(real64), intent(in) :: rtol, stol

    real(real64), intent(inout) :: time(:), hmin(:), hout(:)
    real(real64), intent(inout) :: hpout(:), hline(:)
    integer, intent(inout) :: hnit(:)

    real(real64), intent(inout) :: hmaxbal
    integer, intent(out) :: info

    info = PIPE_INVALID_ARGUMENT

    if (n < 2 .or. n > 100) return
    if (nsteps < 0) return

    if (size(po) < n+1) return
    if (size(mo) < n) return

    if (size(time) < nsteps+1) return
    if (size(hmin) < nsteps+1) return
    if (size(hout) < nsteps+1) return
    if (size(hpout) < nsteps+1) return
    if (size(hline) < nsteps+1) return
    if (size(hnit) < nsteps+1) return

    call integrate_transient(n, po, mo, mino, pinn, &
         mdot0, mdot1, tramp, dt, nsteps, theta, dx, &
         d, a, temp, z, rs, mu, eps, rtol, stol, maxit, &
         time, hmin, hout, hpout, hline, hnit, hmaxbal, info)

  end subroutine pipe_integrate_checked


  subroutine pipe_adaptive_checked(n, po, mo, mino, &
       pinn, mdot0, mdot1, tramp, duration, theta, dx, &
       d, a, temp, z, rs, mu, eps, rtol, stol, maxit, &
       dtinit, dtmin, dtmax, etol, cflmax, bfrac, &
       maxrej, maxrec, time, dtrec, hmin, hout, &
       hpout, hline, hnit, heta, hbal, hcfl, &
       nacc, nrej, ntotal, nnewton, info)

    integer, intent(in) :: n, maxit, maxrej, maxrec

    real(real64), intent(inout) :: po(:), mo(:), mino

    real(real64), intent(in) :: pinn, mdot0, mdot1, tramp
    real(real64), intent(in) :: duration, theta, dx
    real(real64), intent(in) :: d, a, temp, z, rs, mu, eps
    real(real64), intent(in) :: rtol, stol
    real(real64), intent(in) :: dtinit, dtmin, dtmax
    real(real64), intent(in) :: etol, cflmax, bfrac

    real(real64), intent(inout) :: time(:), dtrec(:)
    real(real64), intent(inout) :: hmin(:), hout(:)
    real(real64), intent(inout) :: hpout(:), hline(:)
    real(real64), intent(inout) :: heta(:), hbal(:), hcfl(:)

    integer, intent(inout) :: hnit(:)

    integer, intent(out) :: nacc, nrej, ntotal, nnewton
    integer, intent(out) :: info

    info = PIPE_INVALID_ARGUMENT

    ! Legacy validation does not define counters on early failure.
    ! Initialize them here for a deterministic checked interface.
    nacc = 0
    nrej = 0
    ntotal = 0
    nnewton = 0

    if (n < 2 .or. n > 100) return
    if (maxrec < 1) return

    if (size(po) < n+1) return
    if (size(mo) < n) return

    if (size(time) < maxrec+1) return
    if (size(dtrec) < maxrec+1) return
    if (size(hmin) < maxrec+1) return
    if (size(hout) < maxrec+1) return
    if (size(hpout) < maxrec+1) return
    if (size(hline) < maxrec+1) return
    if (size(hnit) < maxrec+1) return
    if (size(heta) < maxrec+1) return
    if (size(hbal) < maxrec+1) return
    if (size(hcfl) < maxrec+1) return

    call integrate_transient_adaptive(n, po, mo, mino, &
         pinn, mdot0, mdot1, tramp, duration, theta, &
         dx, d, a, temp, z, rs, mu, eps, rtol, stol, &
         maxit, dtinit, dtmin, dtmax, etol, cflmax, &
         bfrac, maxrej, maxrec, time, dtrec, hmin, &
         hout, hpout, hline, hnit, heta, hbal, hcfl, &
         nacc, nrej, ntotal, nnewton, info)

  end subroutine pipe_adaptive_checked

end module pipe_solver_api
