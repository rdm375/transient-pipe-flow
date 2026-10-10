program m9d3_spatial_regression
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter :: nlevels=3
  integer, parameter :: meshes(nlevels)=[20,40,80]
  integer, parameter :: ntimes=2
  integer, parameter :: nsteps(ntimes)=[640,1280]
  integer, parameter :: nsmax=1280

  real(8), parameter :: dt=0.9375d0
  real(8), parameter :: theta=0.5d0

  real(8), parameter :: qmin(3)=[1.75d0,1.70d0,1.75d0]
  real(8), parameter :: qmax(3)=[2.25d0,2.30d0,2.25d0]

  real(8), parameter :: bal_tol=1.d-7
  real(8), parameter :: diff_floor(3)=[1.d-6,1.d-9,1.d-6]

  integer :: j,k,i,n,ns,info
  integer :: nit(0:nsmax)

  real(8) :: po(0:100),mo(0:99)
  real(8) :: thist(0:nsmax),hin(0:nsmax),hout(0:nsmax)
  real(8) :: hpout(0:nsmax),hline(0:nsmax)

  real(8) :: val(3,nlevels)
  real(8) :: diff(3,nlevels-1)
  real(8) :: order(3)

  real(8) :: pi,d,a,temp,z,rs,mu,eps
  real(8) :: dx,pin,mino
  real(8) :: re,fd,coef,bal

  real(8) :: reynolds_mass,friction_sj

  character(len=16), parameter :: names(3) = [ &
       character(len=16) :: &
       'outlet pressure', &
       'inlet flow', &
       'linepack' ]

  pi=4.d0*atan(1.d0)

  d=1.d0
  a=pi*d*d/4.d0
  temp=288.15d0
  z=0.9d0
  rs=500.d0
  mu=1.1d-5
  eps=4.5d-5

  pin=8.d6

  re=reynolds_mass(100.d0,d,a,mu)
  fd=friction_sj(re,eps,d)

  coef=fd*z*rs*temp*100.d0**2/(d*a*a)

  write(*,'(a)') 'M9D.3 SPATIAL CONVERGENCE REGRESSION'
  write(*,'(a,f10.5,a,f6.3)') &
       'dt=',dt,' theta=',theta

  do k=1,ntimes

     ns=nsteps(k)

     write(*,*)
     write(*,'(a,f10.3,a,i6)') &
          'Observation time=',dt*dble(ns),' steps=',ns

     do j=1,nlevels

        n=meshes(j)
        dx=100000.d0/dble(n)

        do i=0,n
           po(i)=sqrt(pin**2-coef*dble(i)*dx)
        enddo

        mo(0:n-1)=100.d0
        mino=100.d0

        call integrate_transient(n,po,mo,mino,pin, &
             100.d0,105.d0,600.d0,dt,ns,theta, &
             dx,d,a,temp,z,rs,mu,eps, &
             1.d-9,1.d-12,30, &
             thist,hin,hout,hpout,hline, &
             nit,bal,info)

        if (info/=0) then
           write(*,*) 'FAIL solver mesh=',n, &
                ' time=',dt*dble(ns),' info=',info
           stop 1
        endif

        val(1,j)=hpout(ns)
        val(2,j)=hin(ns)
        val(3,j)=hline(ns)

        if (.not.all(ieee_is_finite(val(:,j)))) then
           write(*,*) 'FAIL nonfinite result mesh=',n
           stop 1
        endif

        if (.not.ieee_is_finite(bal)) then
           write(*,*) 'FAIL nonfinite mass balance mesh=',n
           stop 1
        endif

        if (bal>bal_tol .or. bal<0.d0) then
           write(*,*) 'FAIL mass balance mesh=',n, &
                ' maxbal=',bal,' tolerance=',bal_tol
           stop 1
        endif

        write(*,'(a,i4,a,3(1x,es18.9),a,es12.4)') &
             'N=',n,' values:',val(:,j),' maxbal=',bal

     enddo

     do j=1,nlevels-1
        diff(:,j)=abs(val(:,j)-val(:,j+1))
     enddo

     write(*,'(a)') 'Spatial convergence:'

     do i=1,3

        if (.not.all(ieee_is_finite(diff(i,:)))) then
           write(*,*) 'FAIL nonfinite differences: ',trim(names(i))
           stop 1
        endif

        if (any(diff(i,:)<=diff_floor(i))) then
           write(*,*) 'FAIL insignificant mesh differences: ', &
                trim(names(i))
           write(*,*) 'differences=',diff(i,:)
           write(*,*) 'floor=',diff_floor(i)
           stop 1
        endif

        order(i)=log(diff(i,1)/diff(i,2))/log(2.d0)

        if (.not.ieee_is_finite(order(i))) then
           write(*,*) 'FAIL nonfinite order: ',trim(names(i))
           stop 1
        endif

        write(*,'(2x,a16,a,f9.5,a,2(1x,es12.4))') &
             names(i),' order=',order(i), &
             ' differences=',diff(i,:)

        if (order(i)<qmin(i) .or. order(i)>qmax(i)) then
           write(*,*) 'FAIL spatial order: ',trim(names(i))
           write(*,*) 'expected range=',qmin(i),qmax(i)
           stop 1
        endif

     enddo

     write(*,'(a,f8.1,a)') &
          'PASS spatial convergence at t=',dt*dble(ns),' s'

  enddo

  write(*,*)
  write(*,'(a)') 'PASS M9D.3 spatial convergence regression'

end program m9d3_spatial_regression
