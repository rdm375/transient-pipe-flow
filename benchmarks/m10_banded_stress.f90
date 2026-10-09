program m10_banded_stress
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: max_n=202
    integer, parameter :: kl=2, ku=1
    integer, parameter :: ldab=2*kl+ku+1
    integer, parameter :: sizes(7)=[6,10,22,42,82,162,202]
    integer :: n, s, trial, family, i, j, info_d, info_b
    integer :: failures, tests, seed_size
    integer, allocatable :: seed(:)
    real(8) :: a(max_n,max_n), ac(max_n,max_n)
    real(8) :: ab(ldab,max_n), b(max_n)
    real(8) :: xd(max_n), xb(max_n), xtrue(max_n)
    real(8) :: r, err, res_d, res_b, scale, rhs_scale
    real(8) :: diff, worst_diff, worst_residual
    logical :: finite

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed=17391
    call random_seed(put=seed)
    deallocate(seed)

    failures=0
    tests=0
    worst_diff=0.0d0
    worst_residual=0.0d0

    print '(a)', 'M10.3D BANDED SOLVER STRESS TEST'

    do s=1,size(sizes)
        n=sizes(s)

        do family=1,5
            do trial=1,20
                tests=tests+1

                a=0.0d0
                ac=0.0d0
                ab=0.0d0
                b=0.0d0

                ! Generate matrices with KL=2 and KU=1.
                do j=1,n
                    do i=max(1,j-ku),min(n,j+kl)
                        call random_number(r)
                        a(i,j)=2.0d0*r-1.0d0
                    enddo
                enddo

                select case(family)
                case(1)
                    ! Strongly diagonally dominant.
                    do i=1,n
                        a(i,i)=5.0d0
                    enddo

                case(2)
                    ! Force pivots throughout the interior.
                    do i=1,n-2
                        if (mod(i,3)==1) then
                            a(i,i)=1.0d-10
                            a(i+2,i)=3.0d0
                        endif
                    enddo

                case(3)
                    ! Alternating strong subdiagonal entries.
                    do i=1,n-1
                        if (mod(i,2)==1) then
                            a(i,i)=1.0d-8
                            a(i+1,i)=4.0d0
                        endif
                    enddo

                case(4)
                    ! Unstructured random banded matrix.
                    ! Add diagonal stabilization to avoid
                    ! accidental near-singularity.
                    do i=1,n
                        a(i,i)=a(i,i)+1.5d0
                    enddo

                case(5)
                    ! Guaranteed singular: zero final row.
                    do j=max(1,n-kl),n
                        a(n,j)=0.0d0
                    enddo
                end select

                do i=1,n
                    xtrue(i)=sin(real(i+trial,8))
                enddo

                b(1:n)=matmul(a(1:n,1:n),xtrue(1:n))
                ac(1:n,1:n)=a(1:n,1:n)

                do j=1,n
                    do i=max(1,j-ku),min(n,j+kl)
                        ab(kl+ku+1+i-j,j)=a(i,j)
                    enddo
                enddo

                call solve_dense(n,ac,max_n,b,xd,info_d)
                call solve_banded(n,ab,ldab,kl,ku,b,xb,info_b)

                if (family==5) then
                    if (info_d==0 .or. info_b==0) then
                        print *, 'FAIL singular status:', &
                            n,family,trial,info_d,info_b
                        failures=failures+1
                    endif
                    cycle
                endif

                if (info_d/=info_b) then
                    print *, 'FAIL status:', &
                        n,family,trial,info_d,info_b
                    failures=failures+1
                    cycle
                endif

                if (info_d/=0) cycle

                finite=all(ieee_is_finite(xb(1:n)))
                if (.not.finite) then
                    print *, 'FAIL nonfinite:',n,family,trial
                    failures=failures+1
                    cycle
                endif

                scale=max(1.0d0,maxval(abs(xd(1:n))))
                rhs_scale=max(1.0d0,maxval(abs(b(1:n))))

                diff=maxval(abs(xd(1:n)-xb(1:n)))/scale

                res_d=maxval(abs( &
                    matmul(a(1:n,1:n),xd(1:n))-b(1:n)))/rhs_scale

                res_b=maxval(abs( &
                    matmul(a(1:n,1:n),xb(1:n))-b(1:n)))/rhs_scale

                worst_diff=max(worst_diff,diff)
                worst_residual=max(worst_residual,res_b)

                if (diff>1.0d-10 .or. res_b>1.0d-10) then
                    print '(a,3i5,3(a,es13.5))', &
                        'FAIL n/family/trial=',n,family,trial, &
                        ' difference=',diff, &
                        ' dense_res=',res_d, &
                        ' band_res=',res_b
                    failures=failures+1
                endif
            enddo
        enddo

        print '(a,i4,a,i5)', 'N=',n,' tests completed=',100
    enddo

    print '(a,i6)', 'Total tests: ',tests
    print '(a,es14.6)', 'Worst relative difference: ',worst_diff
    print '(a,es14.6)', 'Worst banded residual: ',worst_residual
    print '(a,i6)', 'Failures: ',failures

    if (failures/=0) then
        print '(a)', 'FAIL M10.3D'
        stop 1
    endif

    print '(a)', 'PASS M10.3D: banded solver stress tests'
end program
