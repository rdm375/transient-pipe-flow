module m10_profile
  use, intrinsic :: iso_fortran_env, only: real64, error_unit
  implicit none
  private
  public :: profile_start, profile_stop, profile_report
  integer, save :: counts(4) = 0
  real(real64), save :: elapsed(4) = 0.0_real64
  real(real64), save :: starts(4) = 0.0_real64
contains
  subroutine profile_start(k)
    integer, intent(in) :: k
    counts(k) = counts(k) + 1
    call cpu_time(starts(k))
  end subroutine
  subroutine profile_stop(k)
    integer, intent(in) :: k
    real(real64) :: finish
    call cpu_time(finish)
    elapsed(k) = elapsed(k) + finish - starts(k)
  end subroutine
  subroutine profile_report(total)
    real(real64), intent(in) :: total
    character(len=20), parameter :: labels(4) = [character(len=20) :: &
      'residual', 'jacobian', 'jacobian_copy', 'linear_solve']
    integer :: k
    write(error_unit,'(a)') 'M10_PROFILE_BEGIN'
    do k=1,4
      write(error_unit,'(a,",",i0,",",es24.16)') &
        trim(labels(k)), counts(k), elapsed(k)
    end do
    write(error_unit,'(a,",",i0,",",es24.16)') &
      'other', 0, total - sum(elapsed)
    write(error_unit,'(a,",",i0,",",es24.16)') &
      'total', 0, total
    write(error_unit,'(a)') 'M10_PROFILE_END'
  end subroutine
end module
