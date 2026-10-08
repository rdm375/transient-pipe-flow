program m9e2d_schedule_interface
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use pipe_solver_api
  implicit none

  integer, parameter :: n=20, maxrec=4000
  real(real64), parameter :: sentinel=-1234567.0_real64
  integer, parameter :: isentinel=-999

  real(real64) :: po(n+1),mo(n),mino
  real(real64) :: pol(n+1),mol(n),minl
  real(real64) :: ts(4),ps(4),qs(4),atol(4),rtoly(4)
  real(real64) :: time(maxrec+1),dtrec(maxrec+1)
  real(real64) :: hmin(maxrec+1),hout(maxrec+1)
  real(real64) :: hpout(maxrec+1),hline(maxrec+1)
  real(real64) :: heta(maxrec+1),hbal(maxrec+1),hcfl(maxrec+1)
  real(real64) :: tl(maxrec+1),dl(maxrec+1)
  real(real64) :: hil(maxrec+1),hol(maxrec+1)
  real(real64) :: hpl(maxrec+1),hll(maxrec+1)
  real(real64) :: hel(maxrec+1),hbl(maxrec+1),hcl(maxrec+1)
  integer :: hnit(maxrec+1),nil(maxrec+1)
  integer :: nacc,nrej,ntotal,nnewton,info
  integer :: al,rl,totl,nwl,infol
  integer :: failures,k,j,scenario,ns,mr
  real(real64) :: pi,d,a,temp,z,rs,mu,eps,dx
  real(real64) :: re,fd,coef,x,cum,defect
  real(real64), external :: reynolds_mass,friction_sj

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

  re=reynolds_mass(100.0_real64,d,a,mu)
  fd=friction_sj(re,eps,d)
  coef=fd*z*rs*temp*100.0_real64**2/(d*a*a)

  call initialize(1)

  ! Test each of the ten output histories independently.
  do k=1,10
     call initialize(1)
     select case(k)
     case(1)
        call checked(time(:maxrec),dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl)
     case(2)
        call checked(time,dtrec(:maxrec),hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl)
     case(3)
        call checked(time,dtrec,hmin(:maxrec),hout,hpout, &
             hline,hnit,heta,hbal,hcfl)
     case(4)
        call checked(time,dtrec,hmin,hout(:maxrec),hpout, &
             hline,hnit,heta,hbal,hcfl)
     case(5)
        call checked(time,dtrec,hmin,hout,hpout(:maxrec), &
             hline,hnit,heta,hbal,hcfl)
     case(6)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline(:maxrec),hnit,heta,hbal,hcfl)
     case(7)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit(:maxrec),heta,hbal,hcfl)
     case(8)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta(:maxrec),hbal,hcfl)
     case(9)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal(:maxrec),hcfl)
     case(10)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl(:maxrec))
     end select
     call verify_invalid()
  enddo

  ! Independently undersize all five schedule/scaling arrays.
  do k=11,15
     call initialize(1)
     select case(k)
     case(11)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl,tsin=ts(:3))
     case(12)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl,psin=ps(:3))
     case(13)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl,qsin=qs(:3))
     case(14)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl,atin=atol(:3))
     case(15)
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl,rtin=rtoly(:3))
     end select
     call verify_invalid()
  enddo

  ! Independently undersize pressure and mass-flow states.
  do k=16,17
     call initialize(1)
     if (k==16) then
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl,poin=po(:n))
     else
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl,moin=mo(:n-1))
     endif
     call verify_invalid()
  enddo

  ! Invalid schedule count and record count.
  do k=18,19
     call initialize(1)
     if (k==18) then
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl,nsin=1)
     else
        call checked(time,dtrec,hmin,hout,hpout, &
             hline,hnit,heta,hbal,hcfl,mrin=0)
     endif
     call verify_invalid()
  enddo

  ! Two successful cases: constant and varying inlet pressure.
  do scenario=1,2
     call initialize(scenario)
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
     nil=isentinel
     call integrate_transient_adaptive_schedule(n,pol,mol,minl, &
          4,ts,ps,qs,3600.0_real64,0.65_real64,dx, &
          d,a,temp,z,rs,mu,eps,1.0e-11_real64, &
          1.0e-12_real64,30,7.5_real64,1.0e-5_real64, &
          120.0_real64,0.01_real64,0.0_real64, &
          0.125_real64,20,maxrec,atol,rtoly, &
          tl,dl,hil,hol,hpl,hll,nil,hel,hbl,hcl, &
          al,rl,totl,nwl,infol)

     call checked(time,dtrec,hmin,hout,hpout, &
          hline,hnit,heta,hbal,hcfl)

     call require(infol==PIPE_SUCCESS,'legacy success')
     call require(info==infol,'status parity')
     call require(nacc==al,'accepted parity')
     call require(nrej==rl,'rejected parity')
     call require(ntotal==totl,'attempt parity')
     call require(nnewton==nwl,'Newton parity')
     call require(ntotal==nacc+nrej,'attempt accounting')
     call require(same(po,pol),'pressure state parity')
     call require(same(mo,mol),'flow state parity')
     call require(same([mino],[minl]),'inlet flow parity')

     if (info/=PIPE_SUCCESS) cycle
     j=nacc+1

     call require(nacc>=30.and.nacc<=1000,'accepted range')
     call require(abs(time(j)-3600.0_real64)<1.0e-8_real64, &
          'terminal time')

     call require(same(time(:j),tl(:j)),'time parity')
     call require(same(dtrec(:j),dl(:j)),'dt parity')
     call require(same(hmin(:j),hil(:j)),'inlet history parity')
     call require(same(hout(:j),hol(:j)),'outlet history parity')
     call require(same(hpout(:j),hpl(:j)),'pressure history parity')
     call require(same(hline(:j),hll(:j)),'linepack parity')
     call require(same(heta(:j),hel(:j)),'indicator parity')
     call require(same(hbal(:j),hbl(:j)),'balance parity')
     call require(same(hcfl(:j),hcl(:j)),'CFL parity')
     call require(all(hnit(:j)==nil(:j)),'iteration history parity')

     cum=0.0_real64
     defect=0.0_real64
     do k=2,j
        x=dtrec(k)*(0.65_real64*(hmin(k)-hout(k))+ &
             0.35_real64*(hmin(k-1)-hout(k-1)))
        cum=cum+x
        defect=defect+hbal(k)

        call require(abs(hline(k)-hline(k-1)-x)<1.0e-7_real64, &
             'local mass conservation')
        call require(abs(hbal(k))<1.0e-7_real64, &
             'stored balance defect')
        call require(abs(time(k)-time(k-1)-dtrec(k))< &
             1.0e-9_real64,'time increment')
        call require(dtrec(k)>0.0_real64,'positive timestep')
        call require(heta(k)<=1.00000001_real64, &
             'accepted indicator')

        if (time(k-1)<600.0_real64-1.0e-8_real64) &
             call require(time(k)<=600.0_real64+1.0e-8_real64, &
             'first schedule knot')
        if (time(k-1)<1800.0_real64-1.0e-8_real64) &
             call require(time(k)<=1800.0_real64+1.0e-8_real64, &
             'second schedule knot')
     enddo

     call require(abs(hline(j)-hline(1)-cum)<1.0e-5_real64, &
          'integrated mass conservation')
     call require(abs(defect)<1.0e-5_real64, &
          'integrated numerical defect')
     call require(abs(cum)>100.0_real64, &
          'nontrivial inventory change')

     print '(a,i0,a,i0,a,i0,a,i0)', &
          'M9E.2D scenario=',scenario, &
          ' accepted=',nacc,' rejected=',nrej, &
          ' Newton=',nnewton
     print '(a,es14.6)', &
          'M9E.2D outlet pressure: ',po(n+1)
     print '(a,es14.6)', &
          'M9E.2D linepack change: ',hline(j)-hline(1)
  enddo

  if (failures/=0) then
     print '(a,i0)', 'FAIL M9E.2D checks: ',failures
     stop 1
  endif
  print '(a)', 'PASS M9E.2D schedule interface'

contains

  subroutine initialize(which)
    integer, intent(in) :: which
    integer :: i

    mino=100.0_real64
    do i=1,n+1
       x=real(i-1,real64)*dx
       po(i)=sqrt(8.0e6_real64**2-coef*x)
    enddo
    mo=100.0_real64

    ts=[0.0_real64,600.0_real64,1800.0_real64,3600.0_real64]
    ps=8.0e6_real64
    if (which==2) ps(3:4)=8.001e6_real64
    qs=[100.0_real64,105.0_real64,105.0_real64,105.0_real64]
    atol=[8.0e6_real64,100.0_real64,5000.0_real64,5.0_real64]
    rtoly=0.0_real64

    time=sentinel
    dtrec=sentinel
    hmin=sentinel
    hout=sentinel
    hpout=sentinel
    hline=sentinel
    heta=sentinel
    hbal=sentinel
    hcfl=sentinel
    hnit=isentinel
    nacc=-1
    nrej=-1
    ntotal=-1
    nnewton=-1
    info=-1
  end subroutine initialize

  subroutine checked(t,dt,hi,ho,hp,hl,ni,he,hb,hc, &
       tsin,psin,qsin,atin,rtin,poin,moin,nsin,mrin)
    real(real64), intent(inout) :: t(:),dt(:),hi(:),ho(:)
    real(real64), intent(inout) :: hp(:),hl(:),he(:),hb(:),hc(:)
    integer, intent(inout) :: ni(:)
    real(real64), optional, intent(in) :: tsin(:),psin(:),qsin(:)
    real(real64), optional, intent(in) :: atin(:),rtin(:)
    real(real64), optional, intent(inout) :: poin(:),moin(:)
    integer, optional, intent(in) :: nsin,mrin


    ! The default arguments must alias the actual state and inputs.
    ! Optional undersized arrays are passed directly in each branch.
    ns=4
    mr=maxrec
    if (present(nsin)) ns=nsin
    if (present(mrin)) mr=mrin

    if (present(poin)) then
       call invoke(t,dt,hi,ho,hp,hl,ni,he,hb,hc,poin,mo,ts,ps,qs,atol,rtoly)
    else if (present(moin)) then
       call invoke(t,dt,hi,ho,hp,hl,ni,he,hb,hc,po,moin,ts,ps,qs,atol,rtoly)
    else if (present(tsin)) then
       call invoke(t,dt,hi,ho,hp,hl,ni,he,hb,hc,po,mo,tsin,ps,qs,atol,rtoly)
    else if (present(psin)) then
       call invoke(t,dt,hi,ho,hp,hl,ni,he,hb,hc,po,mo,ts,psin,qs,atol,rtoly)
    else if (present(qsin)) then
       call invoke(t,dt,hi,ho,hp,hl,ni,he,hb,hc,po,mo,ts,ps,qsin,atol,rtoly)
    else if (present(atin)) then
       call invoke(t,dt,hi,ho,hp,hl,ni,he,hb,hc,po,mo,ts,ps,qs,atin,rtoly)
    else if (present(rtin)) then
       call invoke(t,dt,hi,ho,hp,hl,ni,he,hb,hc,po,mo,ts,ps,qs,atol,rtin)
    else
       call invoke(t,dt,hi,ho,hp,hl,ni,he,hb,hc,po,mo,ts,ps,qs,atol,rtoly)
    endif

  end subroutine checked

  subroutine invoke(t,dt,hi,ho,hp,hl,ni,he,hb,hc, &
       p,m,tsv,psv,qsv,av,rv)
      real(real64), intent(inout) :: t(:),dt(:),hi(:),ho(:)
      real(real64), intent(inout) :: hp(:),hl(:),he(:),hb(:),hc(:)
      integer, intent(inout) :: ni(:)
      real(real64), intent(inout) :: p(:),m(:)
      real(real64), intent(in) :: tsv(:),psv(:),qsv(:),av(:),rv(:)

      call pipe_adaptive_schedule_checked(n,p,m,mino, &
           ns,tsv,psv,qsv,3600.0_real64,0.65_real64,dx, &
           d,a,temp,z,rs,mu,eps,1.0e-11_real64, &
           1.0e-12_real64,30,7.5_real64,1.0e-5_real64, &
           120.0_real64,0.01_real64,0.0_real64, &
           0.125_real64,20,mr,av,rv, &
           t,dt,hi,ho,hp,hl,ni,he,hb,hc, &
           nacc,nrej,ntotal,nnewton,info)
  end subroutine invoke

  subroutine verify_invalid()
    real(real64) :: expected(n+1)
    integer :: i

    call require(info==PIPE_INVALID_ARGUMENT,'invalid status')
    call require(nacc==0.and.nrej==0.and. &
         ntotal==0.and.nnewton==0,'invalid counters')

    do i=1,n+1
       expected(i)=sqrt(8.0e6_real64**2- &
            coef*real(i-1,real64)*dx)
    enddo

    call require(same(po,expected),'unchanged pressure')
    call require(same(mo,spread(100.0_real64,1,n)),'unchanged flow')
    call require(same([mino],[100.0_real64]),'unchanged inlet flow')
    call require(same(time,spread(sentinel,1,maxrec+1)),'unchanged time')
    call require(same(dtrec,spread(sentinel,1,maxrec+1)),'unchanged dt')
    call require(same(hmin,spread(sentinel,1,maxrec+1)),'unchanged inlet history')
    call require(same(hout,spread(sentinel,1,maxrec+1)),'unchanged outlet history')
    call require(same(hpout,spread(sentinel,1,maxrec+1)),'unchanged pressure history')
    call require(same(hline,spread(sentinel,1,maxrec+1)),'unchanged linepack')
    call require(same(heta,spread(sentinel,1,maxrec+1)),'unchanged indicator')
    call require(same(hbal,spread(sentinel,1,maxrec+1)),'unchanged balance')
    call require(same(hcfl,spread(sentinel,1,maxrec+1)),'unchanged CFL')
    call require(all(hnit==isentinel),'unchanged iterations')
  end subroutine verify_invalid

  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not.condition) then
       failures=failures+1
       print '(a,a)', 'FAIL: ',message
    endif
  end subroutine require

  logical function same(lhs,rhs)
    real(real64), intent(in) :: lhs(:),rhs(:)
    same=.false.
    if (size(lhs)/=size(rhs)) return
    same=all(transfer(lhs,[0_int64],size(lhs))== &
             transfer(rhs,[0_int64],size(rhs)))
  end function same

end program m9e2d_schedule_interface
