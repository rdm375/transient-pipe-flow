program m9d2_temporal
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter :: n=20, nlevels=7
  integer :: j,k,i,ns,info,nit(0:640),theta_idx
  integer :: failures
  real(8) :: q,lower,upper,dcoarse,dfine

  real(8) :: po(0:n),mo(0:n-1)
  real(8) :: t_hist(0:640),hin(0:640)
  real(8) :: hout(0:640),hpout(0:640),hline(0:640)

  real(8) :: pi,d,a,temp,z,rs,mu,eps,dx,pin,dt,theta
  real(8) :: re,fd,coef,mino,bal,rate
  real(8) :: reynolds_mass,friction_sj

  real(8) :: val(3,nlevels,3)
  real(8) :: err(3,nlevels-1)
  real(8) :: theta_values(3)

  theta_values=[0.5d0,0.65d0,1.d0]
  failures=0

  pi=4.d0*atan(1.d0)
  d=1.d0
  a=pi*d*d/4.d0
  temp=288.15d0
  z=0.9d0
  rs=500.d0
  mu=1.1d-5
  eps=4.5d-5
  dx=100000.d0/dble(n)
  pin=8.d6

  re=reynolds_mass(100.d0,d,a,mu)
  fd=friction_sj(re,eps,d)
  coef=fd*z*rs*temp*100.d0**2/(d*a*a)

  do theta_idx=1,3
     theta=theta_values(theta_idx)

     do j=1,nlevels
        ns=10*2**(j-1)
        dt=300.d0/dble(ns)

        do i=0,n
           po(i)=sqrt(pin**2-coef*dble(i)*dx)
        enddo

        mo=100.d0
        mino=100.d0

        call integrate_transient(n,po,mo,mino,pin,100.d0,105.d0, &
             600.d0,dt,ns,theta,dx,d,a,temp,z,rs,mu,eps, &
             1.d-9,1.d-12,30,t_hist,hin,hout,hpout,hline, &
             nit,bal,info)

        if (info /= 0) then
           write(*,*) 'FAIL theta=',theta,' dt=',dt, &
                      ' info=',info
           stop 1
        endif

        val(1,j,theta_idx)=hpout(ns)
        val(2,j,theta_idx)=hin(ns)
        val(3,j,theta_idx)=hline(ns)
     enddo

     do j=1,nlevels
        do k=1,3
           if (.not.ieee_is_finite(val(k,j,theta_idx))) then
              write(*,*) 'FAIL nonfinite observable:', &
                   theta_idx,j,k
              stop 1
           endif
        enddo
     enddo

     write(*,'(a,f5.2)') 'THETA ',theta
     write(*,'(a)') &
          'dt(s)   outlet_p_error(Pa)  inlet_m_error(kg/s)  linepack_error(kg)'

     do j=1,nlevels-1
        err(:,j)=abs(val(:,j,theta_idx)-val(:,nlevels,theta_idx))

        write(*,'(f8.5,3(1x,es18.9))') &
             300.d0/dble(10*2**(j-1)),err(:,j)
     enddo

     write(*,'(a)') &
          'Three-resolution self-convergence orders:'

     do j=1,nlevels-2
        write(*,'(a,f8.5,a)',advance='no') &
             ' dt=',300.d0/dble(10*2**(j-1)),' orders:'

        do k=1,3
           if (abs(val(k,j+1,theta_idx)- &
                   val(k,j+2,theta_idx)) > 0.d0 .and. &
               abs(val(k,j,theta_idx)- &
                   val(k,j+1,theta_idx)) > 0.d0) then

              rate=log(abs(val(k,j,theta_idx)- &
                           val(k,j+1,theta_idx))/ &
                       abs(val(k,j+1,theta_idx)- &
                           val(k,j+2,theta_idx)))/log(2.d0)

              write(*,'(1x,f8.4)',advance='no') rate
           else
              write(*,'(1x,a8)',advance='no') 'n/a'
           endif
        enddo

        write(*,*)
     enddo

     write(*,'(a)') &
          'Observed orders (coarse-to-fine; reference finest):'

     do j=1,nlevels-2
        write(*,'(a,f8.5,a)',advance='no') &
             ' dt=',300.d0/dble(10*2**(j-1)),' orders:'

        do k=1,3
           if (err(k,j+1)>0.d0 .and. err(k,j)>0.d0) then
              rate=log(err(k,j)/err(k,j+1))/log(2.d0)
              write(*,'(1x,f8.4)',advance='no') rate
           else
              write(*,'(1x,a8)',advance='no') 'n/a'
           endif
        enddo

        write(*,*)
     enddo
     ! Three-resolution regression at the finest available scales.
     !
     ! j=5 compares dt=1.875, 0.9375 and 0.46875 seconds.
     ! k=1: outlet pressure
     ! k=2: inlet mass flow
     ! k=3: linepack

     j=nlevels-2

     do k=1,3
        dcoarse=abs(val(k,j,theta_idx)-val(k,j+1,theta_idx))
        dfine=abs(val(k,j+1,theta_idx)-val(k,j+2,theta_idx))

        if (.not.ieee_is_finite(dcoarse) .or. &
            .not.ieee_is_finite(dfine)) then
           write(*,*) 'FAIL nonfinite difference:',theta_idx,k
           failures=failures+1
           cycle
        endif

        if (dcoarse<=0.d0 .or. dfine<=0.d0) then
           write(*,*) 'FAIL degenerate difference:',theta_idx,k
           failures=failures+1
           cycle
        endif

        q=log(dcoarse/dfine)/log(2.d0)

        if (theta_idx==1) then
           if (k==2) then
              lower=1.50d0
              upper=2.30d0
           else
              lower=1.75d0
              upper=2.25d0
           endif
        else
           if (k==2) then
              if (theta_idx==2) then
                 lower=0.70d0
              else
                 lower=0.65d0
              endif
              upper=1.15d0
           else
              lower=0.85d0
              upper=1.15d0
           endif
        endif

        if (.not.ieee_is_finite(q)) then
           write(*,*) 'FAIL nonfinite order:',theta_idx,k
           failures=failures+1
        else if (q<lower .or. q>upper) then
           write(*,'(a,f5.2,a,i1,a,f9.5)') &
                'FAIL theta=',theta, &
                ' observable=',k,' order=',q
           failures=failures+1
        endif
     enddo

  enddo

  if (failures/=0) then
     write(*,*) 'M9D.2 FAIL:',failures,' convergence checks'
     stop 1
  endif

  write(*,'(a)') 'M9D.2 PASS: temporal convergence'
end program
