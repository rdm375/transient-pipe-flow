program m10_cpu_baseline
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer :: n, ns, info, i, argc, total_newton
  integer, allocatable :: nit(:)
  real(real64), allocatable :: po(:), mo(:)
  real(real64), allocatable :: th(:), hi(:), ho(:), hp(:), hl(:)
  real(real64) :: dt, duration, dx, pin, mino
  real(real64) :: d, area, temp, z, rs, mu, eps
  real(real64) :: re, fd, coef, balance
  real(real64) :: start_cpu, end_cpu, elapsed
  real(real64) :: reynolds_mass, friction_sj
  character(len=64) :: arg

  argc = command_argument_count()
  if (argc /= 3) then
     write(*,'(a)') 'Usage: m10_cpu_baseline N DT DURATION'
     stop 2
  endif

  call get_command_argument(1,arg)
  read(arg,*,iostat=info) n
  if (info /= 0) stop 2

  call get_command_argument(2,arg)
  read(arg,*,iostat=info) dt
  if (info /= 0) stop 2

  call get_command_argument(3,arg)
  read(arg,*,iostat=info) duration
  if (info /= 0) stop 2

  if (n < 2 .or. n > 100) stop 2
  if (.not.ieee_is_finite(dt)) stop 2
  if (.not.ieee_is_finite(duration)) stop 2
  if (dt <= 0.0_real64 .or. duration <= 0.0_real64) stop 2

  ns = nint(duration/dt)
  if (ns < 1) stop 2
  if (abs(real(ns,real64)*dt-duration) > &
       1.0e-10_real64*duration) stop 2

  allocate(po(0:n),mo(0:n-1))
  allocate(th(0:ns),hi(0:ns),ho(0:ns))
  allocate(hp(0:ns),hl(0:ns),nit(0:ns))

  d = 1.0_real64
  area = acos(-1.0_real64)*d*d/4.0_real64
  temp = 288.15_real64
  z = 0.9_real64
  rs = 500.0_real64
  mu = 1.1e-5_real64
  eps = 4.5e-5_real64

  pin = 8.0e6_real64
  mino = 100.0_real64
  dx = 100000.0_real64/real(n,real64)

  re = reynolds_mass(mino,d,area,mu)
  fd = friction_sj(re,eps,d)
  coef = fd*z*rs*temp*mino**2/(d*area**2)

  do i=0,n
     po(i) = sqrt(pin**2-coef*real(i,real64)*dx)
  enddo
  mo = mino

  call cpu_time(start_cpu)

  call integrate_transient(n,po,mo,mino,pin, &
       100.0_real64,105.0_real64,600.0_real64, &
       dt,ns,0.65_real64,dx,d,area,temp,z,rs,mu,eps, &
       1.0e-9_real64,1.0e-12_real64,30, &
       th,hi,ho,hp,hl,nit,balance,info)

  call cpu_time(end_cpu)
  elapsed = end_cpu-start_cpu

  if (info /= 0) then
     write(*,'(a,i0)') 'SOLVER_FAILURE info=',info
     stop 1
  endif

  if (.not.ieee_is_finite(balance)) stop 1
  if (.not.ieee_is_finite(hp(ns))) stop 1
  if (.not.ieee_is_finite(hl(ns))) stop 1

  total_newton = sum(nit(1:ns))

  write(*,'(a)') &
       'n,dt,duration,steps,cpu_seconds,newton_total,' // &
       'newton_per_step,outlet_pressure,linepack,max_mass_defect'

  write(*,'(i0,",",es24.16,",",es24.16,",",i0,",",' // &
       'es24.16,",",i0,",",es24.16,",",es24.16,","' // &
       ',es24.16,",",es24.16)') &
       n,dt,duration,ns,elapsed,total_newton, &
       real(total_newton,real64)/real(ns,real64), &
       hp(ns),hl(ns),balance

end program m10_cpu_baseline
