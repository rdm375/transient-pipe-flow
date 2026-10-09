program m10_banded_parity
    implicit none

    integer, parameter :: max_n=202
    integer, parameter :: kl=2, ku=1
    integer, parameter :: ldab=2*kl+ku+1
    integer :: sizes(7) = [6,10,22,42,82,162,202]
    integer :: n, case_id, s, i, j, info_d, info_b
    integer :: failures
    real(8) :: a(max_n,max_n), ac(max_n,max_n)
    real(8) :: ab(ldab,max_n)
    real(8) :: b(max_n), xd(max_n), xb(max_n)
    real(8) :: err, scale, residual, tmp

    failures=0

    print '(a)', 'M10.3C BANDED/DENSE LINEAR SOLVER PARITY'

    do s=1,size(sizes)
        n=sizes(s)

        do case_id=1,3
            a=0.0d0
            ac=0.0d0
            ab=0.0d0
            b=0.0d0

            do j=1,n
                do i=max(1,j-ku),min(n,j+kl)
                    if (i==j) then
                        a(i,j)=8.0d0
                    else
                        a(i,j)=0.15d0*sin(real(3*i+7*j,8))
                    endif
                enddo
            enddo

            if (case_id==2) then
                ! Force a pivot in the first column.
                a(1,1)=1.0d-8
                a(2,1)=2.0d0
            endif

            if (case_id==3) then
                ! Force a pivot two rows below the diagonal.
                a(1,1)=1.0d-8
                a(2,1)=1.0d0
                a(3,1)=3.0d0
            endif

            do i=1,n
                b(i)=sin(real(i,8))
            enddo

            ac(1:n,1:n)=a(1:n,1:n)

            do j=1,n
                do i=max(1,j-ku),min(n,j+kl)
                    ab(kl+ku+1+i-j,j)=a(i,j)
                enddo
            enddo

            call solve_dense(n,ac,max_n,b,xd,info_d)
            call solve_banded(n,ab,ldab,kl,ku,b,xb,info_b)

            if (info_d/=info_b) then
                print *, 'FAIL status',n,case_id,info_d,info_b
                failures=failures+1
                cycle
            endif

            if (info_d/=0) then
                print *, 'FAIL unexpected singularity',n,case_id
                failures=failures+1
                cycle
            endif

            scale=max(1.0d0,maxval(abs(xd(1:n))))
            err=maxval(abs(xd(1:n)-xb(1:n)))/scale

            residual=0.0d0
            do i=1,n
                tmp=-b(i)
                do j=max(1,i-kl),min(n,i+ku)
                    tmp=tmp+a(i,j)*xb(j)
                enddo
                residual=max(residual,abs(tmp))
            enddo

            print '(a,i4,a,i2,a,es12.4,a,es12.4)', &
                'n=',n,' case=',case_id, &
                ' rel_difference=',err,' residual=',residual

            if (err>1.0d-11 .or. residual>1.0d-10) then
                print *, 'FAIL numerical parity',n,case_id
                failures=failures+1
            endif
        enddo
    enddo

    if (failures/=0) then
        print *, 'FAIL M10.3C:',failures,'cases'
        stop 1
    endif

    print '(a)', 'PASS M10.3C: banded/dense solver parity'
end program
