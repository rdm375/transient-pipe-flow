program m10_direct_banded_assembly
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: max_n=100, max_u=202
    integer, parameter :: kl=2, ku=1, ldab=6
    integer, parameter :: meshes(6)=[2,4,10,20,80,100]

    integer :: n, nu, k, direction, i, j
    integer :: failures, cases, row
    real(8) :: u(max_u), dense(max_u,max_u)
    real(8) :: band(ldab,max_u)
    real(8) :: dx, dt, theta, d, a, t, z, rs, mu, eps
    real(8) :: expected, actual, difference, worst

    failures=0
    cases=0
    worst=0.0d0

    d=0.8d0
    a=acos(-1.0d0)*d*d/4.0d0
    t=288.15d0
    z=0.9d0
    rs=500.0d0
    mu=1.1d-5
    eps=1.0d-5
    dt=15.0d0
    theta=0.65d0

    print '(a)', 'M10.4 DIRECT BANDED ASSEMBLY PARITY'

    do k=1,size(meshes)
        n=meshes(k)
        nu=2*n+2
        dx=100000.0d0/real(n,8)

        do direction=-1,1,2
            cases=cases+1

            u=0.0d0
            u(1)=8.0d6

            do i=2,nu
                if (mod(i,2)==0 .and. i>=4) then
                    u(i)=8.0d6-1000.0d0*real(i,8)
                else
                    u(i)=real(direction,8)*100.0d0
                endif
            enddo

            dense=0.0d0
            band=0.0d0

            call assemble_jacobian(n,u,dx,dt,theta, &
                d,a,t,z,rs,mu,eps,dense,max_u)

            call assemble_jacobian_banded(n,u,dx,dt,theta, &
                d,a,t,z,rs,mu,eps,band,ldab)

            if (.not.all(ieee_is_finite(band))) then
                print *, 'FAIL nonfinite band:',n,direction
                failures=failures+1
                cycle
            endif

            do j=1,nu
                do i=1,nu
                    expected=dense(i,j)
                    actual=0.0d0

                    row=kl+ku+1+i-j

                    if (row>=1 .and. row<=ldab) then
                        actual=band(row,j)
                    endif

                    difference=abs(expected-actual)
                    worst=max(worst,difference)

                    if (difference>0.0d0) then
                        print *, 'FAIL:',n,direction,i,j, &
                            expected,actual
                        failures=failures+1
                    endif
                enddo
            enddo
        enddo

        print '(a,i4,a)', 'N=',n,' complete'
    enddo

    print '(a,i5)', 'Configurations: ',cases
    print '(a,es14.6)', 'Worst absolute difference: ',worst
    print '(a,i5)', 'Failures: ',failures

    if (failures/=0) stop 1

    print '(a)', 'PASS M10.4: direct banded assembly parity'
end program
