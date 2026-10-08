program m9e3_failure_interface
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use pipe_solver_api
  implicit none
  integer, parameter :: n=20, mr=1
  real(real64) :: po(n+1),mo(n),mino,base(n+1)
  real(real64) :: tm(mr+1),dt(mr+1),hi(mr+1),ho(mr+1)
  real(real64) :: hp(mr+1),hl(mr+1),he(mr+1),hb(mr+1),hc(mr+1)
  integer :: hn(mr+1),na,nr,nt,nn,info,i
  real(real64) :: pi,d,a,temp,z,rs,mu,eps,dx,re,fd,coef
  real(real64) :: ts(2),ps(2),qs(2),at(4),ry(4)
  real(real64), external :: reynolds_mass,friction_sj

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
  do i=1,n+1
     base(i)=sqrt(8.0e6_real64**2-coef*real(i-1,real64)*dx)
  enddo
  ts=[0.0_real64,3600.0_real64]
  ps=8.0e6_real64
  qs=[100.0_real64,105.0_real64]
  at=[8.0e6_real64,100.0_real64,5000.0_real64,5.0_real64]
  ry=0.0_real64

  ! Capacity failure after one accepted timestep: ramp wrapper.
  call reset()
  call ramp(30,1)
  call capacity_checks('ramp')

  ! Same progress contract with piecewise-linear schedules.
  call reset()
  call schedule(30,1)
  call capacity_checks('schedule')

  ! A failed first solve must not commit any physical state.
  call reset()
  call ramp(0,0)
  call require(info==6,'ramp first-step status')
  call require(na==0.and.nr==1.and.nt==1,'ramp first-step accounting')
  call require(same(po,base),'ramp initial pressure preserved')
  call require(same(mo,spread(100.0_real64,1,n)).and.same([mino],[100.0_real64]), &
       'ramp initial mass flow preserved')

  call reset()
  call schedule(0,0)
  call require(info==6,'schedule first-step status')
  call require(na==0.and.nr==1.and.nt==1,'schedule first-step accounting')
  call require(same(po,base),'schedule initial pressure preserved')
  call require(same(mo,spread(100.0_real64,1,n)).and.same([mino],[100.0_real64]), &
       'schedule initial mass flow preserved')

  ! Checked capacity rejection must leave all state and histories untouched.
  call reset()
  call ramp(30,1,undersized=.true.)
  call require(info==PIPE_INVALID_ARGUMENT,'invalid capacity status')
  call require(na==0.and.nr==0.and.nt==0.and.nn==0, &
       'invalid capacity counters')
  call require(same(po,base),'invalid capacity pressure preserved')
  call require(same(tm,spread(-12345.0_real64,1,size(tm))),'invalid capacity histories preserved')
  call require(all(hn==-999),'invalid capacity iteration history preserved')

  print '(a)', 'PASS M9E.3 checked failure-progress interface'
contains
  subroutine reset()
    po=base
    mo=100.0_real64
    mino=100.0_real64
    tm=-12345.0_real64
    dt=-12345.0_real64
    hi=-12345.0_real64
    ho=-12345.0_real64
    hp=-12345.0_real64
    hl=-12345.0_real64
    he=-12345.0_real64
    hb=-12345.0_real64
    hc=-12345.0_real64
    hn=-999
    na=-1; nr=-1; nt=-1; nn=-1; info=-1
  end subroutine

  subroutine ramp(iter,rej,undersized)
    integer, intent(in) :: iter,rej
    logical, optional, intent(in) :: undersized
    if (present(undersized)) then
       if (undersized) then
          call invoke_ramp(tm(:mr),iter,rej)
          return
       endif
    endif
    call invoke_ramp(tm,iter,rej)
  end subroutine

  subroutine invoke_ramp(tvec,iter,rej)
    real(real64), intent(inout) :: tvec(:)
    integer, intent(in) :: iter,rej
    call pipe_adaptive_checked(n,po,mo,mino, &
         8.0e6_real64,100.0_real64,105.0_real64, &
         600.0_real64,3600.0_real64,0.65_real64,dx, &
         d,a,temp,z,rs,mu,eps,1.0e-11_real64, &
         1.0e-12_real64,iter,7.5_real64,1.0e-5_real64, &
         120.0_real64,0.01_real64,0.0_real64, &
         0.125_real64,rej,mr, &
         tvec,dt,hi,ho,hp,hl,hn,he,hb,hc,na,nr,nt,nn,info)
  end subroutine

  subroutine schedule(iter,rej)
    integer, intent(in) :: iter,rej
    call pipe_adaptive_schedule_checked(n,po,mo,mino, &
         2,ts,ps,qs,3600.0_real64,0.65_real64,dx, &
         d,a,temp,z,rs,mu,eps,1.0e-11_real64, &
         1.0e-12_real64,iter,7.5_real64,1.0e-5_real64, &
         120.0_real64,0.01_real64,0.0_real64, &
         0.125_real64,rej,mr,at,ry, &
         tm,dt,hi,ho,hp,hl,hn,he,hb,hc,na,nr,nt,nn,info)
  end subroutine

  subroutine capacity_checks(label)
    character(*), intent(in) :: label
    call require(info==7,label//' capacity status')
    call require(na==1,label//' accepted count')
    call require(nt==na+nr,label//' attempt accounting')
    call require(nn>=0,label//' Newton accounting')
    call require(same(tm(1:1),[0.0_real64]).and.tm(2)>0.0_real64, &
         label//' valid accepted times')
    call require(dt(2)>0.0_real64,label//' accepted dt')
    call require(same([mino],[hi(2)]),label//' committed inlet flow')
    call require(same([po(n+1)],[hp(2)]),label//' committed outlet pressure')
    call require(hn(1)==0,label//' initial iteration history')
  end subroutine

  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not.condition) then
       print '(a,a)', 'FAIL M9E.3: ',message
       stop 1
    endif
  end subroutine

  logical function same(lhs,rhs)
    real(real64), intent(in) :: lhs(:),rhs(:)
    same=all(transfer(lhs,[0_int64],size(lhs)) == &
         transfer(rhs,[0_int64],size(rhs)))
  end function
end program
