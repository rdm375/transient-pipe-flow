program m9d5_halfcell_regression
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter :: maxn=100
  integer, parameter :: ncases=2
  integer, parameter :: meshes(ncases)=[4,20]
  real(8), parameter :: weights(3)=[0.5d0,0.65d0,1.d0]
  real(8), parameter :: residual_tol=1.d-10
  real(8), parameter :: jac_tol=1.d-12

  integer :: n,nu,ic,it,idir,i,k,failures
  real(8) :: u(2*maxn+2),po(0:maxn),mo(0:maxn-1)
  real(8) :: r(2*maxn+2),jac(2*maxn+2,2*maxn+2)
  real(8) :: pi,d,a,temp,z,rs,mu,eps,dx,dt,theta
  real(8) :: pin,mino,mouto,moutn,crho,storage
  real(8) :: expected_in,expected_out,expected_bc
  real(8) :: expected_j(2*maxn+2)
  real(8) :: err_in,err_out,err_bc,err_j
  real(8) :: max_rerr,max_jerr,flow_sign

  pi=4.d0*atan(1.d0)
  d=1.d0
  a=pi*d*d/4.d0
  temp=288.15d0
  z=0.9d0
  rs=500.d0
  mu=1.1d-5
  eps=4.5d-5
  dt=37.d0
  pin=8.d6
  crho=1.d0/(z*rs*temp)

  failures=0
  max_rerr=0.d0
  max_jerr=0.d0

  write(*,'(a)') 'M9D.5 BOUNDARY HALF-CELL REGRESSION'

  do ic=1,ncases
     n=meshes(ic)
     nu=2*n+2
     dx=100000.d0/dble(n)
     storage=a*dx*crho/(2.d0*dt)

     do it=1,size(weights)
        theta=weights(it)

        do idir=1,2
           if (idir==1) then
              flow_sign=1.d0
           else
              flow_sign=-1.d0
           endif

           ! Independent non-equilibrium old and new states.
           do i=0,n
              po(i)=pin-1000.d0*dble(i)
              if (i==0) then
                 u(1)=pin+75.d0
              else
                 u(2*i+2)=pin-1000.d0*dble(i)+ &
                      50.d0*dble(i+1)
              endif
           enddo

           do i=0,n-1
              mo(i)=flow_sign*(100.d0+0.2d0*dble(i))
              u(2*i+3)=flow_sign*(102.d0+0.3d0*dble(i))
           enddo

           mino=flow_sign*99.d0
           mouto=flow_sign*101.d0
           moutn=flow_sign*106.d0
           u(2)=flow_sign*103.d0

           call assemble_residual(n,u,po,mo,mino, &
                mouto,moutn,pin,dx,dt,theta, &
                d,a,temp,z,rs,mu,eps,r)

           call assemble_jacobian(n,u,dx,dt,theta, &
                d,a,temp,z,rs,mu,eps,jac,2*maxn+2)

           ! Independently evaluated inlet boundary condition.
           expected_bc=u(1)-pin

           ! Independently evaluated inlet half-cell continuity.
           expected_in=storage*(u(1)-po(0)) &
                +theta*(u(3)-u(2)) &
                +(1.d0-theta)*(mo(0)-mino)

           ! Independently evaluated outlet half-cell continuity.
           expected_out=storage*(u(2*n+2)-po(n)) &
                +theta*(moutn-u(2*n+1)) &
                +(1.d0-theta)*(mouto-mo(n-1))

           err_bc=abs(r(1)-expected_bc)
           err_in=abs(r(2)-expected_in)
           err_out=abs(r(nu)-expected_out)

           max_rerr=max(max_rerr,err_bc,err_in,err_out)

           if (max(err_bc,err_in,err_out)>residual_tol) then
              write(*,*) 'FAIL boundary residual',n,it,idir
              failures=failures+1
           endif

           ! Check the complete inlet pressure boundary row.
           expected_j(1:nu)=0.d0
           expected_j(1)=1.d0

           do k=1,nu
              err_j=abs(jac(1,k)-expected_j(k))
              max_jerr=max(max_jerr,err_j)
              if (err_j>jac_tol) failures=failures+1
           enddo

           ! Check the complete inlet continuity row.
           expected_j(1:nu)=0.d0
           expected_j(1)=storage
           expected_j(2)=-theta
           expected_j(3)=theta

           do k=1,nu
              err_j=abs(jac(2,k)-expected_j(k))
              max_jerr=max(max_jerr,err_j)
              if (err_j>jac_tol) failures=failures+1
           enddo

           ! Check the complete outlet continuity row.
           expected_j(1:nu)=0.d0
           expected_j(2*n+1)=-theta
           expected_j(2*n+2)=storage

           do k=1,nu
              err_j=abs(jac(nu,k)-expected_j(k))
              max_jerr=max(max_jerr,err_j)
              if (err_j>jac_tol) failures=failures+1
           enddo

           if (.not.all(ieee_is_finite(r(1:nu)))) then
              write(*,*) 'FAIL nonfinite residual'
              failures=failures+1
           endif

           if (.not.all(ieee_is_finite(jac(1:nu,1:nu)))) then
              write(*,*) 'FAIL nonfinite Jacobian'
              failures=failures+1
           endif
        enddo
     enddo
  enddo

  write(*,'(a,es12.4)') &
       'maximum boundary residual discrepancy = ',max_rerr
  write(*,'(a,es12.4)') &
       'maximum boundary Jacobian discrepancy = ',max_jerr

  if (failures/=0) then
     write(*,*) 'FAIL M9D.5 failures=',failures
     stop 1
  endif

  write(*,'(a)') 'PASS M9D.5 boundary half-cell regression'
end program
