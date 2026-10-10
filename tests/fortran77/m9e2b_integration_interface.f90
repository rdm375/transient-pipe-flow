program m9e2b_integration_interface
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use pipe_solver_api
  implicit none

  integer, parameter :: n = 4, ns = 4
  real(real64), parameter :: sentinel = -1234567.0_real64

  real(real64) :: po(n+1), mo(n), mino
  real(real64) :: time(ns+1), hmin(ns+1), hout(ns+1)
  real(real64) :: hpout(ns+1), hline(ns+1), hmaxbal
  integer :: hnit(ns+1)

  real(real64) :: pol(n+1), mol(n), minl
  real(real64) :: tl(ns+1), hil(ns+1), hol(ns+1)
  real(real64) :: hpl(ns+1), hll(ns+1), maxl
  integer :: nil(ns+1)

  integer :: info, infol, failures, k

  failures = 0

  ! Test each undersized history independently.
  do k = 1, 6
     call initialize()

     select case(k)
     case(1)
        call checked(time(1:ns), hmin, hout, hpout, hline, hnit)
     case(2)
        call checked(time, hmin(1:ns), hout, hpout, hline, hnit)
     case(3)
        call checked(time, hmin, hout(1:ns), hpout, hline, hnit)
     case(4)
        call checked(time, hmin, hout, hpout(1:ns), hline, hnit)
     case(5)
        call checked(time, hmin, hout, hpout, hline(1:ns), hnit)
     case(6)
        call checked(time, hmin, hout, hpout, hline, hnit(1:ns))
     end select

     if (info /= PIPE_INVALID_ARGUMENT) failures = failures+1
     call verify_unchanged()
  enddo

  ! Insufficient pressure capacity.
  call initialize()
  call pipe_integrate_checked(n, po(1:n), mo, mino, &
       8.0e6_real64, 100.0_real64, 105.0_real64, &
       120.0_real64, 60.0_real64, ns, 0.65_real64, &
       1000.0_real64, 1.0_real64, 1.0_real64, &
       288.15_real64, 0.9_real64, 500.0_real64, &
       1.1e-5_real64, 4.5e-5_real64, 1.0e-9_real64, &
       1.0e-12_real64, 30, &
       time, hmin, hout, hpout, hline, hnit, hmaxbal, info)

  if (info /= PIPE_INVALID_ARGUMENT) failures = failures+1
  call verify_unchanged()

  ! Insufficient internal mass-flow capacity.
  call initialize()
  call pipe_integrate_checked(n, po, mo(1:n-1), mino, &
       8.0e6_real64, 100.0_real64, 105.0_real64, &
       120.0_real64, 60.0_real64, ns, 0.65_real64, &
       1000.0_real64, 1.0_real64, 1.0_real64, &
       288.15_real64, 0.9_real64, 500.0_real64, &
       1.1e-5_real64, 4.5e-5_real64, 1.0e-9_real64, &
       1.0e-12_real64, 30, &
       time, hmin, hout, hpout, hline, hnit, hmaxbal, info)

  if (info /= PIPE_INVALID_ARGUMENT) failures = failures+1
  call verify_unchanged()

  ! Valid integration: compare complete trajectories bitwise.
  call initialize()

  pol = po
  mol = mo
  minl = mino

  tl = time
  hil = hmin
  hol = hout
  hpl = hpout
  hll = hline
  nil = hnit
  maxl = hmaxbal

  call integrate_transient(n, pol, mol, minl, &
       8.0e6_real64, 100.0_real64, 105.0_real64, &
       120.0_real64, 60.0_real64, ns, 0.65_real64, &
       1000.0_real64, 1.0_real64, 1.0_real64, &
       288.15_real64, 0.9_real64, 500.0_real64, &
       1.1e-5_real64, 4.5e-5_real64, 1.0e-9_real64, &
       1.0e-12_real64, 30, &
       tl, hil, hol, hpl, hll, nil, maxl, infol)

  call checked(time, hmin, hout, hpout, hline, hnit)

  if (infol /= PIPE_SUCCESS) failures = failures+1
  if (info /= infol) failures = failures+1

  if (.not. same_real(po, pol)) failures = failures+1
  if (.not. same_real(mo, mol)) failures = failures+1
  if (.not. same_real([mino], [minl])) failures = failures+1

  if (.not. same_real(time, tl)) failures = failures+1
  if (.not. same_real(hmin, hil)) failures = failures+1
  if (.not. same_real(hout, hol)) failures = failures+1
  if (.not. same_real(hpout, hpl)) failures = failures+1
  if (.not. same_real(hline, hll)) failures = failures+1
  if (any(hnit /= nil)) failures = failures+1
  if (.not. same_real([hmaxbal], [maxl])) failures = failures+1

  if (infol == PIPE_SUCCESS .and. info == PIPE_SUCCESS) then
     write(*,'(a,es14.6)') &
          'M9E.2B final outlet pressure: ', po(n+1)
     write(*,'(a,es14.6)') &
          'M9E.2B maximum mass-balance defect: ', hmaxbal
  endif

  if (failures /= 0) then
     write(*,'(a,i0)') 'FAIL M9E.2B integration interface: ', failures
     stop 1
  endif

  write(*,'(a)') 'PASS M9E.2B integration interface'

contains

  subroutine initialize()
    po = 8.0e6_real64
    mo = 100.0_real64
    mino = 100.0_real64

    time = sentinel
    hmin = sentinel
    hout = sentinel
    hpout = sentinel
    hline = sentinel
    hnit = -999
    hmaxbal = sentinel
  end subroutine initialize

  subroutine checked(t, hi, ho, hp, hl, ni)
    real(real64), intent(inout) :: t(:), hi(:), ho(:)
    real(real64), intent(inout) :: hp(:), hl(:)
    integer, intent(inout) :: ni(:)

    call pipe_integrate_checked(n, po, mo, mino, &
         8.0e6_real64, 100.0_real64, 105.0_real64, &
         120.0_real64, 60.0_real64, ns, 0.65_real64, &
         1000.0_real64, 1.0_real64, 1.0_real64, &
         288.15_real64, 0.9_real64, 500.0_real64, &
         1.1e-5_real64, 4.5e-5_real64, 1.0e-9_real64, &
         1.0e-12_real64, 30, &
         t, hi, ho, hp, hl, ni, hmaxbal, info)
  end subroutine checked

  subroutine verify_unchanged()
    if (.not. same_real(po, &
         spread(8.0e6_real64, 1, n+1))) failures = failures+1

    if (.not. same_real(mo, &
         spread(100.0_real64, 1, n))) failures = failures+1

    if (.not. same_real([mino], &
         [100.0_real64])) failures = failures+1

    if (.not. same_real(time, &
         spread(sentinel, 1, ns+1))) failures = failures+1

    if (.not. same_real(hmin, &
         spread(sentinel, 1, ns+1))) failures = failures+1

    if (.not. same_real(hout, &
         spread(sentinel, 1, ns+1))) failures = failures+1

    if (.not. same_real(hpout, &
         spread(sentinel, 1, ns+1))) failures = failures+1

    if (.not. same_real(hline, &
         spread(sentinel, 1, ns+1))) failures = failures+1

    if (any(hnit /= -999)) failures = failures+1

    if (.not. same_real([hmaxbal], &
         [sentinel])) failures = failures+1
  end subroutine verify_unchanged

  logical function same_real(x, y)
    real(real64), intent(in) :: x(:), y(:)

    same_real = .false.
    if (size(x) /= size(y)) return

    same_real = all(transfer(x, [0_int64], size(x)) == &
                    transfer(y, [0_int64], size(y)))
  end function same_real

end program m9e2b_integration_interface
