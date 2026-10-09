program m13_banded_solver_parity
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: kl = 2, ku = 1, ldab = 6
    integer :: n, scenario
    integer :: comparisons, failures
    real(dp) :: max_error

    comparisons = 0
    failures = 0
    max_error = 0.0_dp

    do n = 2, 100
        do scenario = 0, 5
            call run_case(n, scenario)
        end do
    end do

    print '(a,i0)', 'Comparisons: ', comparisons
    print '(a,es24.16)', 'Maximum relative error: ', max_error

    if (failures /= 0) then
        print '(a,i0)', 'FAIL discrepancies: ', failures
        stop 1
    end if

    print '(a)', 'PASS M13.6a Fortran banded solver parity'

contains

    subroutine set_band(ab, i, j, value)
        real(dp), intent(inout) :: ab(:,:)
        integer, intent(in) :: i, j
        real(dp), intent(in) :: value
        integer :: row

        row = kl + ku + 1 + i - j

        if (row >= 1 .and. row <= ldab) then
            ab(row,j) = value
        end if
    end subroutine

    subroutine compare(a, b, n, scenario, label)
        real(dp), intent(in) :: a, b
        integer, intent(in) :: n, scenario
        character(*), intent(in) :: label
        real(dp) :: scale, error, relative

        comparisons = comparisons + 1

        if (.not. ieee_is_finite(a) .or. &
            .not. ieee_is_finite(b)) then
            failures = failures + 1
            print *, 'FAIL nonfinite ', label, n, scenario
            return
        end if

        scale = max(abs(a), abs(b))
        error = abs(a-b)

        if (scale > 0.0_dp) then
            relative = error/scale
        else
            relative = 0.0_dp
        end if

        max_error = max(max_error, relative)

        if (a /= b) then
            failures = failures + 1
            if (failures <= 20) then
                print *, 'FAIL ', label, n, scenario, a, b
            end if
        end if
    end subroutine

    subroutine run_case(n, scenario)
        integer, intent(in) :: n, scenario

        real(dp) :: ab(ldab,n), original(ldab,n)
        real(dp) :: b(n), x_generic(n), x_special(n)
        real(dp) :: ab_special(ldab,n)
        integer :: i, j, info_generic, info_special

        original = 0.0_dp

        do i = 1, n
            call set_band(original,i,i,5.0_dp+0.01_dp*i)

            if (i+1 <= n) then
                call set_band(original,i+1,i,-0.5_dp)
                call set_band(original,i,i+1,0.25_dp)
            end if

            if (i+2 <= n) then
                call set_band(original,i+2,i,0.125_dp)
            end if

            b(i) = 1.0_dp+0.1_dp*i
        end do

        select case (scenario)
        case (1)
            call set_band(original,1,1,0.01_dp)
            call set_band(original,2,1,2.0_dp)

        case (2)
            call set_band(original,1,1,0.0_dp)
            call set_band(original,2,1,0.0_dp)
            if (n > 2) then
                call set_band(original,3,1,0.0_dp)
            end if

        case (3)
            call set_band(original,n,n,0.0_dp)
            if (n >= 2) then
                call set_band(original,n,n-1,0.0_dp)
            end if
            if (n >= 3) then
                call set_band(original,n,n-2,0.0_dp)
            end if

        case (4)
            j = max(1,n/2)
            call set_band(original,j,j,0.001_dp)
            if (j < n) then
                call set_band(original,j+1,j,3.0_dp)
            end if

        case (5)
            call set_band(original,1,1,1.0e-12_dp)
            call set_band(original,2,1,1.0e-11_dp)
            if (n > 2) then
                call set_band(original,3,1,1.0e-10_dp)
            end if
        end select

        ab = original
        ab_special = original
        x_generic = -999.0_dp
        x_special = -999.0_dp

        call solve_banded(n,ab,ldab,kl,ku,b,x_generic,info_generic)
        call solve_banded_21(n,ab_special,b,x_special,info_special)

        comparisons = comparisons + 1

        if (info_generic /= info_special) then
            failures = failures + 1
            print *, 'FAIL INFO ', n, scenario, &
                     info_generic, info_special
        end if

        do j = 1, n
            do i = 1, ldab
                call compare(ab(i,j),ab_special(i,j), &
                             n,scenario,'matrix')
            end do
        end do

        do i = 1, n
            call compare(x_generic(i),x_special(i), &
                         n,scenario,'solution')
        end do
    end subroutine

end program
