program m9d4_inventory_regression
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter :: nm=4, ns=360
  integer, parameter :: meshes(nm)=[10,20,40,80]
  real(8), parameter :: dt=60.d0, theta=0.65d0
  real(8), parameter :: mdot0=100.d0, mdot1=105.d0
  real(8), parameter :: tramp=600.d0
  real(8), parameter :: step_tol=1.d-7
  real(8), parameter :: inventory_tol=1.d-5

  integer :: j,k,i,n,info
  integer :: nit(0:ns)

  real(8) :: po(0:100),mo(0:99)
  real(8) :: time(0:ns),hmin(0:ns),hout(0:ns)
  real(8) :: hpout(0:ns),hline(0:ns)

  real(8) :: pi,d,a,temp,z,rs,mu,eps
  real(8) :: pin,dx,mino,re,fd,coef
  real(8) :: expected,imbalance,weighted
  real(8) :: inventory,step_defect,inventory_defect
  real(8) :: max_step,max_inventory,max_physical
  real(8) :: reported_maxbal
  real(8) :: reynolds_mass,friction_sj,demand_ramp

  pi=4.d0*atan(1.d0)
  d=1.d0
  a=pi*d*d/4.d0
  temp=288.15d0
  z=0.9d0
  rs=500.d0
  mu=1.1d-5
  eps=4.5d-5
  pin=8.d6

  re=reynolds_mass(mdot0,d,a,mu)
  fd=friction_sj(re,eps,d)
  coef=fd*z*rs*temp*mdot0**2/(d*a*a)

  write(*,'(a)') 'M9D.4 BOUNDARY AND INVENTORY REGRESSION'

  do j=1,nm
     n=meshes(j)
     dx=100000.d0/dble(n)

     do i=0,n
        po(i)=sqrt(pin**2-coef*dble(i)*dx)
     enddo

     mo(0:n-1)=mdot0
     mino=mdot0

     call integrate_transient(n,po,mo,mino,pin, &
          mdot0,mdot1,tramp,dt,ns,theta, &
          dx,d,a,temp,z,rs,mu,eps, &
          1.d-11,1.d-12,30, &
          time,hmin,hout,hpout,hline, &
          nit,reported_maxbal,info)

     if (info/=0) then
        write(*,*) 'FAIL integration N=',n,' info=',info
        stop 1
     endif

     if (.not.all(ieee_is_finite(hmin))) stop 1
     if (.not.all(ieee_is_finite(hout))) stop 1
     if (.not.all(ieee_is_finite(hline))) stop 1
     if (.not.all(ieee_is_finite(hpout))) stop 1
     if (.not.all(ieee_is_finite(time))) stop 1

     if (abs(po(0)-pin)>1.d-8) then
        write(*,*) 'FAIL inlet pressure N=',n
        stop 1
     endif

     inventory=hline(0)
     max_step=0.d0
     max_inventory=0.d0
     max_physical=0.d0

     if (abs(time(0))>1.d-12) then
        write(*,*) 'FAIL initial time N=',n
        stop 1
     endif

     if (abs(hout(0)-mdot0)>1.d-12) then
        write(*,*) 'FAIL initial outlet demand N=',n
        stop 1
     endif

     max_physical=abs(hmin(0)-hout(0))

     do k=1,ns
        expected=demand_ramp(time(k),mdot0,mdot1,tramp)

        if (abs(hout(k)-expected)>1.d-12) then
           write(*,*) 'FAIL outlet demand N=',n,' step=',k
           stop 1
        endif

        if (abs(time(k)-dt*dble(k))>1.d-9) then
           write(*,*) 'FAIL time history N=',n,' step=',k
           stop 1
        endif

        imbalance=hmin(k)-hout(k)
        max_physical=max(max_physical,abs(imbalance))

        weighted=theta*(hmin(k)-hout(k)) &
             +(1.d0-theta)*(hmin(k-1)-hout(k-1))

        step_defect=(hline(k)-hline(k-1))/dt-weighted
        max_step=max(max_step,abs(step_defect))

        inventory=inventory+dt*weighted
        inventory_defect=inventory-hline(k)
        max_inventory=max(max_inventory,abs(inventory_defect))

        if (abs(step_defect)>step_tol) then
           write(*,*) 'FAIL step conservation N=',n,' k=',k
           write(*,*) 'defect=',step_defect
           stop 1
        endif

        if (abs(inventory_defect)>inventory_tol) then
           write(*,*) 'FAIL inventory reconstruction N=',n,' k=',k
           write(*,*) 'defect=',inventory_defect
           stop 1
        endif
     enddo

     if (.not.ieee_is_finite(reported_maxbal)) then
        write(*,*) 'FAIL nonfinite reported balance N=',n
        stop 1
     endif

     if (abs(reported_maxbal-max_step)>1.d-10) then
        write(*,*) 'FAIL reported maximum balance N=',n
        stop 1
     endif

     if (max_physical<1.d-3) then
        write(*,*) 'FAIL insufficient physical imbalance N=',n
        stop 1
     endif

     if (hline(ns)>=hline(0)) then
        write(*,*) 'FAIL linepack did not decrease N=',n
        stop 1
     endif

     write(*,'(a,i4,a,es12.4,a,es12.4,a,es12.4)') &
          'N=',n, &
          ' max step defect=',max_step, &
          ' max inventory defect=',max_inventory, &
          ' max physical imbalance=',max_physical
  enddo

  write(*,'(a)') 'PASS M9D.4 boundary and inventory regression'
end program
