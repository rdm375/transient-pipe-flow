program m10_jacobian_numeric_audit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: max_n = 100
    integer, parameter :: max_u = 2*max_n + 2
    integer, parameter :: meshes(6) = [2, 4, 10, 20, 80, 100]
    integer :: k, n, nu, direction, i, j, interval
    integer :: lower, upper, observed, expected
    real(8) :: u(max_u), jac(max_u,max_u)
    real(8) :: dx, dt, theta, diameter, area, temp, z, rs, mu, eps
    real(8) :: flow, value
    logical :: structural
    integer :: failures

    failures = 0

    diameter = 0.8d0
    area = acos(-1.0d0)*diameter**2/4.0d0
    temp = 288.15d0
    z = 0.9d0
    rs = 500.0d0
    mu = 1.1d-5
    eps = 1.0d-5
    dt = 15.0d0
    theta = 0.65d0

    print '(a)', 'M10.3B NUMERICAL JACOBIAN BANDWIDTH AUDIT'

    do k = 1, size(meshes)
        n = meshes(k)
        nu = 2*n + 2
        dx = 100000.0d0 / real(n,8)

        do direction = -1, 1, 2
            flow = real(direction,8)*100.0d0

            u = 0.0d0
            do i = 1, nu
                if (i == 1 .or. (mod(i,2) == 0 .and. i >= 4)) then
                    u(i) = 8.0d6 - 1000.0d0*real(i,8)
                else
                    u(i) = flow
                endif
            enddo

            jac = 0.0d0

            call assemble_jacobian(n,u,dx,dt,theta,diameter,area, &
                 temp,z,rs,mu,eps,jac,max_u)

            lower = 0
            upper = 0
            observed = 0
            expected = 0

            do j = 1, nu
                do i = 1, nu
                    structural = .false.

                    if (i == 1 .and. j == 1) structural = .true.
                    if (i == 2 .and. j >= 1 .and. j <= 3) structural = .true.

                    do interval = 0, n-1
                        block
                            integer :: im, il, ir, row

                            im = 2*interval + 3
                            il = 1
                            if (interval > 0) il = 2*interval + 2
                            ir = 2*interval + 4
                            row = 2*interval + 3

                            if (i == row) then
                                if (j == il .or. j == im .or. j == ir) structural = .true.
                            endif

                            if (i == row+1) then
                                if (j == ir .or. j == im) structural = .true.
                                if (interval < n-1 .and. j == im+2) structural = .true.
                            endif
                        end block
                    enddo

                    if (structural) expected = expected + 1

                    value = jac(i,j)

                    if (.not. ieee_is_finite(value)) then
                        print *, 'FAIL nonfinite:', n, direction, i, j
                        failures = failures + 1
                    endif

                    if (value /= 0.0d0) then
                        observed = observed + 1
                        lower = max(lower,i-j)
                        upper = max(upper,j-i)

                        if (.not. structural) then
                            print *, 'FAIL unexpected nonzero:', n, direction, i, j, value
                            failures = failures + 1
                        endif
                    endif
                enddo
            enddo

            print '(a,i3,a,i4,a,i2,a,i4,a,i4,a,i2,a,i2)', &
                'N=',n,' unknowns=',nu,' direction=',direction, &
                ' structural=',expected,' observed=',observed, &
                ' KL=',lower,' KU=',upper

            if (lower /= 2 .or. upper /= 1) then
                print *, 'FAIL bandwidth:', n, direction
                failures = failures + 1
            endif

            if (expected /= 6*n + 3) then
                print *, 'FAIL structural count:', n, expected
                failures = failures + 1
            endif
        enddo
    enddo

    if (failures /= 0) then
        print *, 'FAIL M10.3B:', failures, 'errors'
        stop 1
    endif

    print '(a)', 'PASS M10.3B: numerical Jacobian bandwidth'
end program
