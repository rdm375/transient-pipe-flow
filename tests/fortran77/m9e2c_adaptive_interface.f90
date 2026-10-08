program m9e2c_adaptive_interface
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use pipe_solver_api
  implicit none

  integer, parameter :: n=20, maxrec=4000
  real(real64), parameter :: sentinel=-1234567.0_real64
  integer, parameter :: intsentinel=-999

  real(real64) :: po(n+1), mo(n), mino
  real(real64) :: pol(n+1), mol(n), minl
  real(real64) :: time(maxrec+1), dtrec(maxrec+1)
  real(real64) :: hmin(maxrec+1), hout(maxrec+1)
  real(real64) :: hpout(maxrec+1), hline(maxrec+1)
  real(real64) :: heta(maxrec+1), hbal(maxrec+1)
  real(real64) :: hcfl(maxrec+1)

  real(real64) :: tl(maxrec+1), dl(maxrec+1)
  real(real64) :: hil(maxrec+1), hol(maxrec+1)
  real(real64) :: hpl(maxrec+1), hll(maxrec+1)
  real(real64) :: hel(maxrec+1), hbl(maxrec+1)
  real(real64) :: hcl(maxrec+1)

  integer :: hnit(maxrec+1), nil(maxrec+1)
  integer :: nacc,nrej,ntotal,nnewton,info
  integer :: al,rl,totl,nwl,infol
  integer :: failures,k,j
  real(real64) :: pi,d,a,temp,z,rs,mu,eps,dx
  real(real64) :: reynolds,fd,coef,x,cum,defect
  real(real64), external :: reynolds_mass, friction_sj

  failures=0

  pi=4.0_real64*atan(1.0_real64)
  d=1.0_real64
  a=pi*d*d/4.0_real64
  temp=288.15_real64
  z=0.9_real64
  rs=500.0_real64
  mu=1.1e-5_real64
  eps=4.5e-5_real64
  dx=100000.0_real64/real(n,real64)

  reynolds=reynolds_mass(100.0_real64,d,a,mu)
  fd=friction_sj(reynolds,eps,d)
  coef=fd*z*rs*temp*100.0_real64**2/(d*a*a)

  ! Each of the ten history arrays is checked independently.
  do k=1,10
     call initialize()

     select case(k)
     case(1)
        call checked(time(:maxrec),dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl,maxrec)
     case(2)
        call checked(time,dtrec(:maxrec),hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl,maxrec)
     case(3)
        call checked(time,dtrec,hmin(:maxrec),hout,hpout, &
             hline,hnit,heta,hbal,hcfl,maxrec)
     case(4)
        call checked(time,dtrec,hmin,hout(:maxrec),hpout, &
             hline,hnit,heta,hbal,hcfl,maxrec)
     case(5)
        call checked(time,dtrec,hmin,hout,hpout(:maxrec), &
             hline,hnit,heta,hbal,hcfl,maxrec)
     case(6)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline(:maxrec),hnit,heta,hbal,hcfl,maxrec)
     case(7)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit(:maxrec),heta,hbal,hcfl,maxrec)
     case(8)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta(:maxrec),hbal,hcfl,maxrec)
     case(9)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal(:maxrec),hcfl,maxrec)
     case(10)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl(:maxrec),maxrec)
     end select

     call verify_invalid()
  enddo

  ! Insufficient pressure capacity.
  call initialize()
  call pipe_adaptive_checked(n,po(:n),mo,mino, &
       8.0e6_real64,100.0_real64,105.0_real64, &
       600.0_real64,3600.0_real64,0.65_real64,dx, &
       d,a,temp,z,rs,mu,eps,1.0e-11_real64, &
       1.0e-12_real64,30,7.5_real64,1.0e-5_real64, &
       120.0_real64,0.01_real64,0.0_real64, &
       0.125_real64,20,maxrec, &
       time,dtrec,hmin,hout,hpout,hline,hnit,heta,hbal,hcfl, &
       nacc,nrej,ntotal,nnewton,info)
  call verify_invalid()

  ! Insufficient mass-flow capacity.
  call initialize()
  call pipe_adaptive_checked(n,po,mo(:n-1),mino, &
       8.0e6_real64,100.0_real64,105.0_real64, &
       600.0_real64,3600.0_real64,0.65_real64,dx, &
       d,a,temp,z,rs,mu,eps,1.0e-11_real64, &
       1.0e-12_real64,30,7.5_real64,1.0e-5_real64, &
       120.0_real64,0.01_real64,0.0_real64, &
       0.125_real64,20,maxrec, &
       time,dtrec,hmin,hout,hpout,hline,hnit,heta,hbal,hcfl, &
       nacc,nrej,ntotal,nnewton,info)
  call verify_invalid()

  ! Invalid MAXREC.
  call initialize()
  call checked(time,dtrec,hmin,hout,hpout, &
       hline,hnit,heta,hbal,hcfl,0)
  call verify_invalid()

  ! Successful integration: compare legacy and checked interfaces.
  call initialize()

  pol=po
  mol=mo
  minl=mino

  tl=sentinel
  dl=sentinel
  hil=sentinel
  hol=sentinel
  hpl=sentinel
  hll=sentinel
  hel=sentinel
  hbl=sentinel
  hcl=sentinel
  nil=intsentinel

  call integrate_transient_adaptive(n,pol,mol,minl, &
       8.0e6_real64,100.0_real64,105.0_real64, &
       600.0_real64,3600.0_real64,0.65_real64,dx, &
       d,a,temp,z,rs,mu,eps,1.0e-11_real64, &
       1.0e-12_real64,30,7.5_real64,1.0e-5_real64, &
       120.0_real64,0.01_real64,0.0_real64, &
       0.125_real64,20,maxrec, &
       tl,dl,hil,hol,hpl,hll,nil,hel,hbl,hcl, &
       al,rl,totl,nwl,infol)

  call checked(time,dtrec,hmin,hout,hpout, &
       hline,hnit,heta,hbal,hcfl,maxrec)

  if (infol /= PIPE_SUCCESS) failures=failures+1
  if (info /= infol) failures=failures+1

  if (nacc /= al) failures=failures+1
  if (nrej /= rl) failures=failures+1
  if (ntotal /= totl) failures=failures+1
  if (nnewton /= nwl) failures=failures+1

  if (.not.same_real(po,pol)) failures=failures+1
  if (.not.same_real(mo,mol)) failures=failures+1
  if (.not.same_real([mino],[minl])) failures=failures+1

  if (infol == PIPE_SUCCESS .and. info == PIPE_SUCCESS) then
     j=nacc+1

     if (nacc < 30 .or. nacc > 1000) failures=failures+1
     if (ntotal /= nacc+nrej) failures=failures+1
     if (abs(time(j)-3600.0_real64)>1.0e-8_real64) &
          failures=failures+1

     if (.not.same_real(time(:j),tl(:j))) failures=failures+1
     if (.not.same_real(dtrec(:j),dl(:j))) failures=failures+1
     if (.not.same_real(hmin(:j),hil(:j))) failures=failures+1
     if (.not.same_real(hout(:j),hol(:j))) failures=failures+1
     if (.not.same_real(hpout(:j),hpl(:j))) failures=failures+1
     if (.not.same_real(hline(:j),hll(:j))) failures=failures+1
     if (.not.same_real(heta(:j),hel(:j))) failures=failures+1
     if (.not.same_real(hbal(:j),hbl(:j))) failures=failures+1
     if (.not.same_real(hcfl(:j),hcl(:j))) failures=failures+1
     if (any(hnit(:j)/=nil(:j))) failures=failures+1

     ! Independently reconstruct the inventory change.
     cum=0.0_real64
     defect=0.0_real64

     do k=2,j
        x=dtrec(k)*(0.65_real64*(hmin(k)-hout(k))+ &
             0.35_real64*(hmin(k-1)-hout(k-1)))

        cum=cum+x
        defect=defect+hbal(k)

        if (abs(hline(k)-hline(k-1)-x)>1.0e-7_real64) &
             failures=failures+1
        if (abs(hbal(k))>1.0e-7_real64) failures=failures+1
        if (abs(time(k)-time(k-1)-dtrec(k))>1.0e-9_real64) &
             failures=failures+1
        if (dtrec(k)<=0.0_real64) failures=failures+1
        if (heta(k)>1.00000001_real64) failures=failures+1
     enddo

     if (abs(hline(j)-hline(1)-cum)>1.0e-5_real64) &
          failures=failures+1
     if (abs(defect)>1.0e-5_real64) failures=failures+1
     if (abs(cum)<100.0_real64) failures=failures+1

     write(*,'(a,i0)') 'M9E.2C accepted steps: ',nacc
     write(*,'(a,i0)') 'M9E.2C rejected steps: ',nrej
     write(*,'(a,i0)') 'M9E.2C Newton iterations: ',nnewton
     write(*,'(a,es14.6)') &
          'M9E.2C final outlet pressure: ',po(n+1)
     write(*,'(a,es14.6)') &
          'M9E.2C inventory change: ',hline(j)-hline(1)
  endif

  if (failures /= 0) then
     write(*,'(a,i0)') &
          'FAIL M9E.2C adaptive interface: ',failures
     stop 1
  endif

  print '(a)', 'PASS M9E.2C adaptive interface'

contains

  subroutine initialize()
    integer :: i

    mino=100.0_real64

    do i=1,n+1
       x=real(i-1,real64)*dx
       po(i)=sqrt(8.0e6_real64**2-coef*x)
    enddo

    mo=100.0_real64

    time=sentinel
    dtrec=sentinel
    hmin=sentinel
    hout=sentinel
    hpout=sentinel
    hline=sentinel
    heta=sentinel
    hbal=sentinel
    hcfl=sentinel
    hnit=intsentinel

    nacc=-1
    nrej=-1
    ntotal=-1
    nnewton=-1
    info=-1
  end subroutine initialize

  subroutine checked(t,dt,hi,ho,hp,hl,ni,he,hb,hc,capacity)
    real(real64), intent(inout) :: t(:),dt(:),hi(:),ho(:)
    real(real64), intent(inout) :: hp(:),hl(:),he(:),hb(:),hc(:)
    integer, intent(inout) :: ni(:)
    integer, intent(in) :: capacity

    call pipe_adaptive_checked(n,po,mo,mino, &
         8.0e6_real64,100.0_real64,105.0_real64, &
         600.0_real64,3600.0_real64,0.65_real64,dx, &
         d,a,temp,z,rs,mu,eps,1.0e-11_real64, &
         1.0e-12_real64,30,7.5_real64,1.0e-5_real64, &
         120.0_real64,0.01_real64,0.0_real64, &
         0.125_real64,20,capacity, &
         t,dt,hi,ho,hp,hl,ni,he,hb,hc, &
         nacc,nrej,ntotal,nnewton,info)
  end subroutine checked

  subroutine verify_invalid()
    real(real64) :: expected(n+1)
    integer :: i

    if (info /= PIPE_INVALID_ARGUMENT) failures=failures+1
    if (nacc /= 0) failures=failures+1
    if (nrej /= 0) failures=failures+1
    if (ntotal /= 0) failures=failures+1
    if (nnewton /= 0) failures=failures+1

    do i=1,n+1
       expected(i)=sqrt(8.0e6_real64**2- &
            coef*real(i-1,real64)*dx)
    enddo

    if (.not.same_real(po,expected)) failures=failures+1
    if (.not.same_real(mo, &
         spread(100.0_real64,1,n))) failures=failures+1
    if (.not.same_real([mino],[100.0_real64])) failures=failures+1

    if (.not.same_real(time, &
         spread(sentinel,1,maxrec+1))) failures=failures+1
    if (.not.same_real(dtrec, &
         spread(sentinel,1,maxrec+1))) failures=failures+1
    if (.not.same_real(hmin, &
         spread(sentinel,1,maxrec+1))) failures=failures+1
    if (.not.same_real(hout, &
         spread(sentinel,1,maxrec+1))) failures=failures+1
    if (.not.same_real(hpout, &
         spread(sentinel,1,maxrec+1))) failures=failures+1
    if (.not.same_real(hline, &
         spread(sentinel,1,maxrec+1))) failures=failures+1
    if (.not.same_real(heta, &
         spread(sentinel,1,maxrec+1))) failures=failures+1
    if (.not.same_real(hbal, &
         spread(sentinel,1,maxrec+1))) failures=failures+1
    if (.not.same_real(hcfl, &
         spread(sentinel,1,maxrec+1))) failures=failures+1

    if (any(hnit/=intsentinel)) failures=failures+1
  end subroutine verify_invalid

  logical function same_real(lhs,rhs)
    real(real64), intent(in) :: lhs(:),rhs(:)

    same_real=.false.
    if (size(lhs)/=size(rhs)) return

    same_real=all(transfer(lhs,[0_int64],size(lhs)) == &
                  transfer(rhs,[0_int64],size(rhs)))
  end function same_real

end program m9e2c_adaptive_interface
