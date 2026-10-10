program m10_physical_jacobian_parity
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: max_n=100
    integer, parameter :: max_u=2*max_n+2
    integer, parameter :: kl=2, ku=1
    integer, parameter :: ldab=2*kl+ku+1
    integer, parameter :: meshes(6)=[2,4,10,20,80,100]

    integer :: n, nu, mesh, direction, state, weight
    integer :: i, j, info_d, info_b
    integer :: failures, tests, successes
    real(8) :: u(max_u), jac(max_u,max_u)
    real(8) :: ac(max_u,max_u), ab(ldab,max_u)
    real(8) :: rhs(max_u), xtrue(max_u)
    real(8) :: xd(max_u), xb(max_u)
    real(8) :: dx, dt, theta, diameter, area
    real(8) :: temp, z, rs, mu, eps
    real(8) :: pressure, flow, pressure_drop
    real(8) :: parity, known_error, backward_error
    real(8) :: scale, row_sum, matrix_norm
    real(8) :: residual, denominator
    real(8) :: worst_parity, worst_known, worst_backward

    failures=0
    tests=0
    successes=0

    worst_parity=0.0d0
    worst_known=0.0d0
    worst_backward=0.0d0

    diameter=0.8d0
    area=acos(-1.0d0)*diameter**2/4.0d0

    temp=288.15d0
    z=0.9d0
    rs=500.0d0
    mu=1.1d-5
    eps=1.0d-5

    print '(a)', 'M10.3E PHYSICAL JACOBIAN PARITY'
    print '(a)', 'Production EOS + friction + analytic Jacobian'

    do mesh=1,size(meshes)
        n=meshes(mesh)
        nu=2*n+2
        dx=100000.0d0/real(n,8)

        do direction=-1,1,2
            do state=1,3
                do weight=1,2
                    tests=tests+1

                    select case(weight)
                    case(1)
                        theta=0.60d0
                    case(2)
                        theta=0.70d0
                    end select

                    select case(state)
                    case(1)
                        pressure=5.0d6
                        flow=25.0d0
                        dt=5.0d0
                        pressure_drop=2.0d4
                    case(2)
                        pressure=8.0d6
                        flow=100.0d0
                        dt=15.0d0
                        pressure_drop=5.0d4
                    case(3)
                        pressure=12.0d6
                        flow=250.0d0
                        dt=60.0d0
                        pressure_drop=1.0d5
                    end select

                    flow=flow*real(direction,8)

                    u=0.0d0
                    jac=0.0d0
                    ac=0.0d0
                    ab=0.0d0
                    rhs=0.0d0
                    xtrue=0.0d0
                    xd=0.0d0
                    xb=0.0d0

                    ! Unknown ordering:
                    ! p(0), m(0), m(1/2), p(1), ...
                    !
                    ! In the assembled system, index 1 is inlet
                    ! pressure, index 2 is inlet mass flow,
                    ! odd indices >=3 are interior face flows,
                    ! and even indices >=4 are pressures.

                    u(1)=pressure

                    do i=2,nu
                        if (mod(i,2)==0 .and. i>=4) then
                            u(i)=pressure-pressure_drop * &
                                real(i-2,8)/real(2*n,8)
                        else
                            u(i)=flow*(1.0d0+ &
                                0.05d0*sin(real(i,8)))
                        endif
                    enddo

                    call assemble_jacobian(n,u,dx,dt,theta, &
                        diameter,area,temp,z,rs,mu,eps,jac,max_u)

                    if (.not.all(ieee_is_finite(jac(1:nu,1:nu)))) then
                        print *, 'FAIL nonfinite Jacobian:', &
                            n,direction,state,weight
                        failures=failures+1
                        cycle
                    endif

                    ! Known correction: bounded, nontrivial,
                    ! with components at every unknown.

                    do i=1,nu
                        xtrue(i)=sin(0.37d0*real(i,8)) + &
                            0.2d0*cos(0.13d0*real(i,8))
                    enddo

                    rhs(1:nu)=matmul(jac(1:nu,1:nu),xtrue(1:nu))

                    ac(1:nu,1:nu)=jac(1:nu,1:nu)

                    ! Pack the verified KL=2, KU=1 band.

                    do j=1,nu
                        do i=max(1,j-ku),min(nu,j+kl)
                            ab(kl+ku+1+i-j,j)=jac(i,j)
                        enddo
                    enddo

                    call solve_dense(nu,ac,max_u,rhs,xd,info_d)
                    call solve_banded(nu,ab,ldab,kl,ku, &
                        rhs,xb,info_b)

                    if (info_d/=info_b) then
                        print *, 'FAIL solver status:', &
                            n,direction,state,weight,info_d,info_b
                        failures=failures+1
                        cycle
                    endif

                    if (info_d/=0) then
                        print *, 'FAIL unexpected pivot failure:', &
                            n,direction,state,weight,info_d
                        failures=failures+1
                        cycle
                    endif

                    if (.not.all(ieee_is_finite(xb(1:nu)))) then
                        print *, 'FAIL nonfinite solution:', &
                            n,direction,state,weight
                        failures=failures+1
                        cycle
                    endif

                    successes=successes+1

                    scale=max(1.0d0,maxval(abs(xd(1:nu))))

                    parity=maxval(abs(xb(1:nu)-xd(1:nu)))/scale

                    scale=max(1.0d0,maxval(abs(xtrue(1:nu))))

                    known_error=maxval( &
                        abs(xb(1:nu)-xtrue(1:nu)))/scale

                    matrix_norm=0.0d0

                    do i=1,nu
                        row_sum=sum(abs(jac(i,1:nu)))
                        matrix_norm=max(matrix_norm,row_sum)
                    enddo

                    residual=maxval(abs( &
                        matmul(jac(1:nu,1:nu),xb(1:nu)) &
                        -rhs(1:nu)))

                    denominator=matrix_norm*maxval(abs(xb(1:nu))) &
                        +maxval(abs(rhs(1:nu)))

                    backward_error=residual/max(denominator, &
                        tiny(1.0d0))

                    worst_parity=max(worst_parity,parity)
                    worst_known=max(worst_known,known_error)
                    worst_backward=max(worst_backward,backward_error)

                    if (parity>1.0d-9 .or. &
                        known_error>1.0d-9 .or. &
                        backward_error>1.0d-11) then

                        print '(a,4i5,3(a,es13.5))', &
                            'FAIL case=',n,direction,state,weight, &
                            ' parity=',parity, &
                            ' known=',known_error, &
                            ' backward=',backward_error

                        failures=failures+1
                    endif
                enddo
            enddo
        enddo

        print '(a,i4,a,i4)', &
            'N=',n,' configurations completed=',12
    enddo

    print *
    print '(a,i6)', 'Configurations: ',tests
    print '(a,i6)', 'Successful solves: ',successes
    print '(a,es14.6)', 'Worst dense/banded difference: ', &
        worst_parity
    print '(a,es14.6)', 'Worst known-solution error: ', &
        worst_known
    print '(a,es14.6)', 'Worst relative backward error: ', &
        worst_backward
    print '(a,i6)', 'Failures: ',failures

    if (failures/=0) then
        print '(a)', 'FAIL M10.3E'
        stop 1
    endif

    print '(a)', 'PASS M10.3E: physical Jacobian parity'
end program
