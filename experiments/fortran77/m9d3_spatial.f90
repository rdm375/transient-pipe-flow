program m9d3_spatial
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter :: nlevels=4, ns=320
  integer, parameter :: meshes(nlevels)=[10,20,40,80]
  integer :: j,i,n,info
  integer :: nit(0:ns)

  real(8) :: po(0:100),mo(0:99)
  real(8) :: thist(0:ns),hin(0:ns),hout(0:ns)
  real(8) :: hpout(0:ns),hline(0:ns)
  real(8) :: val(3,nlevels)

  real(8) :: pi,d,a,temp,z,rs,mu,eps
  real(8) :: dx,pin,dt,theta,mino
  real(8) :: re,fd,coef,bal
  real(8) :: diff(3,nlevels-1),rate

  real(8) :: reynolds_mass,friction_sj

  pi=4.d0*atan(1.d0)
  d=1.d0
  a=pi*d*d/4.d0
  temp=288.15d0
  z=0.9d0
  rs=500.d0
  mu=1.1d-5
  eps=4.5d-5

  pin=8.d6
  theta=0.5d0
  dt=300.d0/dble(ns)

  re=reynolds_mass(100.d0,d,a,mu)
  fd=friction_sj(re,eps,d)
  coef=fd*z*rs*temp*100.d0**2/(d*a*a)

  write(*,'(a)') 'M9D.3 SPATIAL CONVERGENCE'
  write(*,'(a,f8.5,a,f5.2)') &
       'dt=',dt,' theta=',theta

  do j=1,nlevels
     n=meshes(j)
     dx=100000.d0/dble(n)

     do i=0,n
        po(i)=sqrt(pin**2-coef*dble(i)*dx)
     enddo

     mo(0:n-1)=100.d0
     mino=100.d0

     call integrate_transient(n,po,mo,mino,pin,100.d0,105.d0, &
          600.d0,dt,ns,theta,dx,d,a,temp,z,rs,mu,eps, &
          1.d-9,1.d-12,30,thist,hin,hout,hpout,hline, &
          nit,bal,info)

     if (info/=0) then
        write(*,*) 'FAIL mesh=',n,' info=',info
        stop 1
     endif

     val(1,j)=hpout(ns)
     val(2,j)=hin(ns)
     val(3,j)=hline(ns)

     if (.not.all(ieee_is_finite(val(:,j)))) then
        write(*,*) 'FAIL nonfinite result mesh=',n
        stop 1
     endif

     write(*,'(a,i4,a,3(1x,es18.9),a,es12.4)') &
          'N=',n,' values:',val(:,j),' maxbal=',bal
  enddo

  write(*,'(a)') 'Successive mesh differences:'
  write(*,'(a)') &
       'N       outlet_p(Pa)       inlet_m(kg/s)      linepack(kg)'

  do j=1,nlevels-1
     diff(:,j)=abs(val(:,j)-val(:,j+1))
     write(*,'(i4,3(1x,es18.9))') &
          meshes(j),diff(:,j)
  enddo

  write(*,'(a)') 'Three-mesh self-convergence orders:'

  do j=1,nlevels-2
     write(*,'(a,i4,a)',advance='no') &
          ' N=',meshes(j),' orders:'

     do i=1,3
        if (diff(i,j)>0.d0 .and. diff(i,j+1)>0.d0) then
           rate=log(diff(i,j)/diff(i,j+1))/log(2.d0)
           write(*,'(1x,f9.4)',advance='no') rate
        else
           write(*,'(1x,a9)',advance='no') 'n/a'
        endif
     enddo

     write(*,*)
  enddo

  write(*,'(a)') 'M9D.3 characterization complete'
end program
