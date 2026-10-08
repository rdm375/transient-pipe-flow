program m9e2_interface_check
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use pipe_solver_api
  implicit none

  integer :: info, niter, failures
  real(real64) :: po(5), mo(4), u(10)
  real(real64) :: sentinel
  real(real64) :: u_legacy(10), u_checked(10)
  integer :: info_legacy, info_checked
  integer :: niter_legacy, niter_checked

  failures = 0
  sentinel = -1234567.0_real64

  po = 8.0e6_real64
  mo = 100.0_real64
  u = sentinel

  ! N=4 requires PO(5), MO(4), U(10).
  ! First, deliberately supply insufficient output capacity.

  call pipe_step_checked(4, po, mo, 100.0_real64, &
       100.0_real64, 100.0_real64, 8.0e6_real64, &
       1000.0_real64, 60.0_real64, 0.65_real64, &
       1.0_real64, 1.0_real64, 288.15_real64, &
       0.9_real64, 500.0_real64, 1.1e-5_real64, &
       4.5e-5_real64, 1.0e-9_real64, 1.0e-12_real64, &
       30, 0, u(1:9), info, niter)

  if (info /= PIPE_INVALID_ARGUMENT) failures = failures+1
  if (niter /= 0) failures = failures+1
  if (any(transfer(u, [0_int64], size(u)) /= transfer(sentinel, 0_int64))) failures = failures+1

  ! Insufficient pressure capacity.

  call pipe_step_checked(4, po(1:4), mo, 100.0_real64, &
       100.0_real64, 100.0_real64, 8.0e6_real64, &
       1000.0_real64, 60.0_real64, 0.65_real64, &
       1.0_real64, 1.0_real64, 288.15_real64, &
       0.9_real64, 500.0_real64, 1.1e-5_real64, &
       4.5e-5_real64, 1.0e-9_real64, 1.0e-12_real64, &
       30, 0, u, info, niter)

  if (info /= PIPE_INVALID_ARGUMENT) failures = failures+1
  if (any(transfer(u, [0_int64], size(u)) /= transfer(sentinel, 0_int64))) failures = failures+1

  ! Insufficient internal-flow capacity.

  call pipe_step_checked(4, po, mo(1:3), 100.0_real64, &
       100.0_real64, 100.0_real64, 8.0e6_real64, &
       1000.0_real64, 60.0_real64, 0.65_real64, &
       1.0_real64, 1.0_real64, 288.15_real64, &
       0.9_real64, 500.0_real64, 1.1e-5_real64, &
       4.5e-5_real64, 1.0e-9_real64, 1.0e-12_real64, &
       30, 0, u, info, niter)

  if (info /= PIPE_INVALID_ARGUMENT) failures = failures+1
  if (any(transfer(u, [0_int64], size(u)) /= transfer(sentinel, 0_int64))) failures = failures+1


  ! Successful-call parity against the original solver.
  !
  ! Use a non-equilibrium state: increase outlet demand while
  ! retaining the same inlet pressure.

  u_legacy = sentinel
  u_checked = sentinel

  call transient_step(4, po, mo, 100.0_real64, &
       100.0_real64, 105.0_real64, 8.0e6_real64, &
       1000.0_real64, 60.0_real64, 0.65_real64, &
       1.0_real64, 1.0_real64, 288.15_real64, &
       0.9_real64, 500.0_real64, 1.1e-5_real64, &
       4.5e-5_real64, 1.0e-9_real64, 1.0e-12_real64, &
       30, 0, u_legacy, info_legacy, niter_legacy)

  call pipe_step_checked(4, po, mo, 100.0_real64, &
       100.0_real64, 105.0_real64, 8.0e6_real64, &
       1000.0_real64, 60.0_real64, 0.65_real64, &
       1.0_real64, 1.0_real64, 288.15_real64, &
       0.9_real64, 500.0_real64, 1.1e-5_real64, &
       4.5e-5_real64, 1.0e-9_real64, 1.0e-12_real64, &
       30, 0, u_checked, info_checked, niter_checked)

  if (info_legacy /= PIPE_SUCCESS) failures = failures+1
  if (info_checked /= info_legacy) failures = failures+1
  if (niter_checked /= niter_legacy) failures = failures+1

  if (any(transfer(u_checked, [0_int64], size(u_checked)) /= &
          transfer(u_legacy, [0_int64], size(u_legacy)))) then
     failures = failures+1
  endif

  if (info_legacy == PIPE_SUCCESS .and. &
      info_checked == PIPE_SUCCESS) then
     write(*,'(a,i0)') 'M9E.2 parity Newton iterations: ', &
          niter_checked
     write(*,'(a,es14.6)') 'M9E.2 parity outlet pressure: ', &
          u_checked(10)
  endif

  if (failures /= 0) then
     write(*,*) 'FAIL M9E.2 interface capacity checks:',failures
     stop 1
  endif

  write(*,'(a)') 'PASS M9E.2 interface capacity checks'

end program
