program m13_banded_solver_benchmark
    use iso_fortran_env, only: real64
    implicit none

    integer, parameter :: n = 202
    integer, parameter :: ldab = 6, kl = 2, ku = 1
    integer, parameter :: repeats = 100000

    real(real64) :: original(ldab,n), ab(ldab,n)
    real(real64) :: b(n), x(n), checksum
    real(real64) :: t0, t1
    integer :: i, j, k, method, info, round_number
    integer :: order(2)

    original = 0.0_real64

    do i = 1, n
        call set_band(i,i,5.0_real64+0.01_real64*i)

        if (i < n) then
            call set_band(i+1,i,-0.5_real64)
            call set_band(i,i+1,0.25_real64)
        end if

        if (i+2 <= n) then
            call set_band(i+2,i,0.125_real64)
        end if

        b(i) = 1.0_real64+0.1_real64*i
    end do

    do round_number = 1, 10
        if (mod(round_number,2) == 0) then
            order = [2,1]
        else
            order = [1,2]
        end if

        do j = 1, 2
            method = order(j)
            checksum = 0.0_real64

            call cpu_time(t0)

            do k = 1, repeats
                ab = original

                if (method == 1) then
                    call solve_banded( &
                        n,ab,ldab,kl,ku,b,x,info)
                else
                    call solve_banded_21( &
                        n,ab,b,x,info)
                end if

                if (info /= 0) then
                    print *, 'FAIL solver:', method, info
                    stop 1
                end if

                checksum = checksum + x(mod(k-1,n)+1)
            end do

            call cpu_time(t1)

            if (method == 1) then
                write(*,'(a,i0,a,es14.6,a,es20.12)') &
                    'round=',round_number, &
                    ' generic_seconds=',t1-t0, &
                    ' checksum=',checksum
            else
                write(*,'(a,i0,a,es14.6,a,es20.12)') &
                    'round=',round_number, &
                    ' special_seconds=',t1-t0, &
                    ' checksum=',checksum
            end if
        end do
    end do

contains

    subroutine set_band(i,j,value)
        integer, intent(in) :: i,j
        real(real64), intent(in) :: value
        integer :: row

        row = kl+ku+1+i-j

        if (row >= 1 .and. row <= ldab) then
            original(row,j) = value
        end if
    end subroutine

end program
