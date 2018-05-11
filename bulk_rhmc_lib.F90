module remez
  use params
  implicit none
  save

  real(dp) :: anum2(0:ndiag), aden2(ndiag)
  real(dp) :: bnum2(0:ndiag), bden2(ndiag)
  real(dp) :: anum4(0:ndiag), aden4(ndiag)
  real(dp) :: bnum4(0:ndiag), bden4(ndiag)
end module remez

module remezg
  use params
  implicit none
  save

  real(dp) :: anum2g(0:ndiagg), aden2g(ndiagg)
  real(dp) :: bnum2g(0:ndiagg), bden2g(ndiagg)
  real(dp) :: anum4g(0:ndiagg), aden4g(ndiagg)
  real(dp) :: bnum4g(0:ndiagg), bden4g(ndiagg)
end module remezg

module trial
  use params
  implicit none
  save

  complex(dp) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
  real :: theta(ksizex_l, ksizey_l, ksizet_l, 3)
  real :: pp(ksizex_l, ksizey_l, ksizet_l, 3)
end module trial

module gauge
  use params
  implicit none
  save
  
  real :: theta(ksizex_l, ksizey_l, ksizet_l, 3)
end module gauge

module vector
  use params
  implicit none
  save

  complex(dp) :: X(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
end module vector

module gforce
  use params
  implicit none
  save

  real :: dSdpi(ksizex_l, ksizey_l, ksizet_l, 3)
end module 

  

module param
  implicit none
  save

  real :: ancg, ancgh, ancgf, ancgpf
  real :: ancgpv, ancghpv, ancgfpv, ancgpfpv
end module param

module dum1
  use params
  implicit none
  save

  complex(dp) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
  real :: ps(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 2)
end module dum1

module phizero
  use params
  implicit none
  save

  complex(dp) :: Phi0(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4, 25)
end module phizero

module dirac
  use params
  implicit none
  save

  complex(dp) :: gamval(6,4)
  integer :: gamin(6,4)
end module dirac

module qmrherm_scratch
  use params
  implicit none

!   complex :: vtild(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
!   complex :: q(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
!   complex :: pm1(kthird, ksizex_l, ksizey_l, ksizet_l, 4,ndiagq)
!   complex :: qm1(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
!   complex :: p(kthird, ksizex_l, ksizey_l, ksizet_l, 4,ndiagq)
!   complex :: x3(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
!   complex :: R(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
!   complex x1(kthird, ksizex_l, ksizey_l, ksizet_l, 4,ndiagq)
!   complex :: x2(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: vtild(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
  complex(dp) :: q(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
  complex(dp) :: pm1(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
  complex(dp) :: qm1(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: p(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
  complex(dp) :: x3(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
  complex(dp) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
  complex(dp) :: x1(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
  complex(dp) :: x2(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
end module qmrherm_scratch

module dwf3d_lib
  use params
  implicit none
  save

  ! Random numbers
  real(dp) :: seed

! Useful constants
  real, parameter :: One = 1.0
  real, parameter :: tpi = 2.0*acos(-1.0)
contains

  subroutine dwf3d_main
    use random
    use remez
    use remezg
    use trial, ut=>u, thetat=>theta
    use gauge
    use vector, X1=>X
    use gforce
    use param
    use dum1
    use comms
!*******************************************************************
!    Rational Hybrid Monte Carlo algorithm for bulk Thirring Model with Domain Wall
!         fermions
!
!    requires operation of QMR on complex vectors to determine
!    (Mdagger M)**-1  Phi 
!
!    requires input partial fraction coefficients from Remez algorithm
!    of Clarke & Kennedy
!
!    the "third" direction is actually indexed 4 in the code - sorry!
!      (should have used C I know)
!
!           { 1 - hermition mass term psibar psi
!    imass= { 3 - antih.    mass term i psibar gamma_3 psi
!           { 5 - antih.    mass term i psibar gamma_5 psi
!     linear combinations of above require code modification
!
!    code exploits fact that gamma_3 is diagonal in Dirac basis used - speeds up 
!    evolution slightly, and measurement by factor of two. 
!
!    Fermion expectation values are measured using a noisy estimator.
!    on the Wilson matrix, which has dimension 4*kvol*kthird
!    inversions done using congrad, and matrix multiplies with dslash, dslashd
!
!    Pauli-Villars determinant defined using a hermitian mass m_h=One
!
!    trajectory length is random with mean dt*iterl
!    The code runs for a fixed number iter2 of trajectories.
!
!    Phi: pseudofermion field 
!    am: bare fermion mass 
!    actiona: running average of total action
!
!                                               SJH February 2017
!*******************************************************************
    real, parameter :: respbp=0.000001, rescgg=0.000001
    real, parameter :: rescga=0.000000001
    real, parameter :: rescgm=0.000000001
    integer, parameter :: itermax=1000
!     complex :: Phi(kthird,0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4, Nf)
!     complex qq,qbqb
!     complex u
!     complex a,b
    complex(dp) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)!
! Currently unused
!    complex(dp) :: qq,qbqb
!    complex(dp) :: u
!    complex(dp) :: a,b
    real(dp) :: H0,H1,S0,S1,dH,dS,hg,hp
    real :: action, paction, gaction
    real :: vel2, x, ytest, atraj
    real :: dt, am, y, traj, proby
    real :: actiona, vel2a, pbp, pbpa, yav, yyav
    real :: ancgm, ancgma
    integer :: imass, iter, iterl, iter2, i, ia, idirac, ithird
    integer :: naccp, ipbp, itot, isweep, itercg, mu
!*******************************************************************
!     variables to keep track of MPI requests
!*******************************************************************
#ifdef MPI
    type(MPI_Request) :: reqs_ps(12)
#else
    integer :: reqs_ps
#endif
!
!*******************************************************************
!     input
!*******************************************************************
    complex(dp), parameter :: zi=(0.0,1.0)
    ibound=-1
#ifdef MPI
    call init_MPI
#endif
!*******************************************************************
!     end of input
!*******************************************************************
!*******************************************************************
!     check qmrherm is going to be OK
!*******************************************************************
    if (ndiagg.gt.ndiag) then
       print *, 'The qmrherm_scratch module currently requires ndiag be greater than ndiagg.'
       print *, 'Please adjust it and recompile.'
       call exit(1)
    endif
    if (ip_global .eq. 0) then
       open(unit=7,file='output',status='unknown')
       open(unit=98,file='control',status='unknown')
    end if
    open(unit=25,file='midout',status='old')
    open(unit=36,file='remez2',status='old')
    open(unit=37,file='remez4',status='old')
    open(unit=38,file='remez2g',status='old')
    open(unit=39,file='remez4g',status='old')
    if(iread.eq.1) then
       call sread
    endif
    read(25,*) dt,beta,am3,am,imass,iterl,iter2
    close(25)
! set a new seed by hand...
    if(iseed.ne.0)then
       seed=4139764973254.0
    endif
    if (ip_global .eq. 0) then
       write(7,*) 'seed: ', seed
    end if
    call init_random(seed)
!     write(6,*) 'ran: ', y,idum
!*******************************************************************
!     initialization
!     istart.lt.0 : start from tape
!     istart=0    : ordered start
!     istart=1    : random start
!*******************************************************************
    call init(istart)
!  read in Remez coefficients
    read(36,*) anum2(0)
    read(37,*) anum4(0)
    read(38,*) anum2g(0)
    read(39,*) anum4g(0)
    do i=1,ndiag
       read(36,*) anum2(i),aden2(i)
       read(37,*) anum4(i),aden4(i)
    enddo
    do i=1,ndiagg
       read(38,*) anum2g(i),aden2g(i)
       read(39,*) anum4g(i),aden4g(i)
    enddo
    read(36,*) bnum2(0)
    read(37,*) bnum4(0)
    read(38,*) bnum2g(0)
    read(39,*) bnum4g(0)
    do i=1,ndiag
       read(36,*) bnum2(i),bden2(i)
       read(37,*) bnum4(i),bden4(i)
    enddo
    do i=1,ndiagg
       read(38,*) bnum2g(i),bden2g(i)
       read(39,*) bnum4g(i),bden4g(i)
    enddo
!*******************************************************************
!     print heading
!*******************************************************************
    traj=iterl*dt
    proby=1.0/float(iterl)
    if (ip_global .eq. 0) then
!       write(6, 9001)ksize,ksizet,kthird,Nf,dt,traj,ndiag,ndiagg, &
!      & iter2,beta,am3,am,imass
       write(7, 9001)ksize,ksizet,kthird,Nf,dt,traj,ndiag,ndiagg, &
            & iter2,beta,am3,am,imass
9001   format(' ksize=',i3,' ksizet=',i3,/ &
          ,' kthird=',i3,/ &
          ,' Nf =',i3,/ &
          ,' time step: dt=',f6.4,' trajectory length=',f9.6,/ &
          ,' Remez ndiag: action =',i3,' guidance=',i3,/ &
          ,' # trajectories=',i6,' beta=',f9.6,/ &
          ,' am3=',f6.4,' am=',f6.4/ &
          ,' imass=',i2)
#ifdef MPI
       write(7, 9002) NP_X, NP_Y, NP_T, ksizex_l, ksizey_l, ksizet_l
9002   format(" NP_X=", i3, " NP_Y=", i3, " NP_T=", i3,/ &
            & " ksizex_l=", i3, " ksizey_l=", i3, " ksizet_l=", i3)
#endif
!     write(6,9004) rescgg,rescga,respbp
       write(7,9004) rescgg,rescga,respbp
9004   format(' Stopping residuals: guidance: ',e11.4,' acceptance: ', &
            &     e11.4,' estimator: ',e11.4)
!     write(6,9044) rescgm
       write(7,9044) rescgm
9044   format(' Stopping residuals: meson: ',e11.4)
       call rranget(seed, 1, 1, 1)
! c     write(6,*) 'seed: ', seed
       write(7,*) 'seed: ', seed
    end if
!*******************************************************************
!       initialize for averages
!*******************************************************************
    actiona=0.0
    vel2a=0.0
    pbpa=0.0
    ancg=0.0
    ancgh=0.0
    ancgf=0.0
    ancgpf=0.0
    ancgpv=0.0
    ancghpv=0.0
    ancgfpv=0.0
    ancgpfpv=0.0
    ancgma=0.0
    yav=0.0
    yyav=0.0 
    naccp=0
    ipbp=0
    itot=0
!*******************************************************************
!     start of classical evolution
!*******************************************************************
    do isweep=1,iter2
! uncomment line below to go straight to measurement
!     goto 666
!*******************************************************************
!     initialise trial fields
!*******************************************************************
       thetat = theta
!
       call coef(ut,thetat)
!*******************************************************************
!  Pseudofermion fields: Phi = {MdaggerM(1)}^-1/4 * {MdaggerM(m)}^1/4 * R, where
!   R is gaussian
!*******************************************************************
       do ia=1,Nf
!
          do idirac=1,4
             do ithird=1,kthird
                call gauss0(ps, reqs_ps)
#ifdef MPI
                call complete_halo_update(reqs_ps)
#endif
                R(ithird,:,:,:,idirac) = cmplx(ps(:,:,:,1), ps(:,:,:,2))
             enddo
          enddo
!
!  For now Phi = {MdaggerM}^0.25 * R
!
          call qmrherm(R,rescga,itercg,am,imass,anum4,aden4,ndiag,0,isweep,0)
          ancgpf=ancgpf+float(itercg)
!
          R = X1
!
          call qmrherm(R,rescga,itercg,One,1,bnum4,bden4,ndiag,0,isweep,0)
          ancgpfpv=ancgpfpv+float(itercg)
!
          Phi = X1
!
       enddo
!*******************************************************************
!     heatbath for p 
!*******************************************************************
!  for some occult reason this write statement is needed to ensure compatibility with earlier versions
!     write(6,*) idum
!     write(98,*) idum
       do mu=1,3
          call gaussp(ps, reqs_ps)
#ifdef MPI
          call complete_halo_update(reqs_ps)
#endif
          pp(:,:,:,mu) = ps(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1)
       enddo
!     write(6,*) idum
!*******************************************************************
!  call to Hamiltonian
!     
       call hamilton(Phi, H0, hg, hp, S0, rescga, isweep, 0, am, imass)
       if(isweep.eq.1) then
          action = real(S0) / kvol
          gaction = real(hg) / kvol
          paction = real(hp) / kvol
       endif
!     goto 501
!*******************************************************************
!      half-step forward for p
!*******************************************************************
       call force(Phi,rescgg,am,imass,isweep,0)
       pp = pp - 0.5 * dt * dSdpi
!*******************************************************************
!     start of main loop for classical time evolution
!*******************************************************************
       do iter=1,itermax
!
!  step (i) st(t+dt)=st(t)+p(t+dt/2)*dt;
!
          thetat = thetat + dt * pp
!
!  step (ii)  p(t+3dt/2)=p(t+dt/2)-dSds(t+dt)*dt (1/2 step on last iteration)
!
          call coef(ut,thetat)
          call force(Phi,rescgg,am,imass,isweep,iter)
!
! test for end of random trajectory
! 
          if (ip_global .eq. 0) then
             ytest = rano(yran, idum, 1, 1, 1)
          end if
#ifdef MPI
          call MPI_Bcast(ytest, 1, MPI_Real, 0, comm)
#endif
          if(ytest.lt.proby)then
             pp = pp - 0.5 * dt * dSdpi
             itot = itot + iter 
             goto 501
          else
             pp = pp - dt * dSdpi
          endif
! 
       end do
!**********************************************************************
!  Monte Carlo step: accept new fields with probability=
!              min(1,exp(H0-H1))
!**********************************************************************
501    continue 
       call hamilton(Phi, H1, hg, hp, S1, rescga, isweep, -1, am, imass)
       dH = H0 - H1
       dS = S0 - S1
       if (ip_global .eq. 0) then
          write(98,*) dH,dS
       end if
       y = exp(real(dH))      
       yav = yav + y 
       yyav = yyav + y*y 
!
       if(dH.lt.0.0)then
          x = rano(yran, idum, 1, 1, 1)
#ifdef MPI
          call MPI_Bcast(x, 1, MPI_Real, 0, comm)
#endif
          if(x.gt.y)goto 600
       endif
!
!     step accepted: set s=st
!
       theta = thetat
       naccp = naccp+1
       action = real(S1) / kvol
       gaction = real(hg) /kvol
       paction = real(hp) /kvol
600    continue
       if (ip_global .eq. 0) then
          write(11,*) isweep,gaction,paction
       end if
       actiona=actiona+action 
       vel2 = sum(pp * pp)
#ifdef MPI
       call MPI_AllReduce(MPI_In_Place, vel2, 1, MPI_Real, MPI_Sum, comm)
#endif
       vel2 = vel2 / (3 * kvol)
       vel2a = vel2a + vel2
!
!     uncomment to disable measurements
!     goto 601
!666    continue
       if((isweep/iprint)*iprint.eq.isweep)then
          thetat = theta
          call coef(ut,thetat)
          call measure(pbp,respbp,ancgm,am,imass)
!        call meson(rescgm,itercg,ancgm,am,imass)
          pbpa=pbpa+pbp
          ancgma=ancgma+ancgm
          ipbp=ipbp+1
!        write(11,*) pbp
!        write(6,*) isweep,':  ',pbp,ancgm
       endif
!
       if((isweep/icheck)*icheck.eq.isweep)then
          call rranget(seed, 1, 1, 1)
          if(iwrite.eq.1) then
             call swrite
          endif
          flush(100)
          flush(200)
!     flush(302)
!     flush(400)
!     flush(500)
!     flush(501)
          if(imass.ne.1)then
!        flush(401)
!        flush(402)
!        flush(403)
          endif
!     write(7,9023) seed
       endif
!
    end do
!*******************************************************************
!     end of main loop
!*******************************************************************
    actiona=actiona/iter2 
    vel2a=vel2a/iter2 
    pbpa=pbpa/ipbp
    ancg=ancg/(Nf*itot)
    ancgh=ancgh/(2*Nf*iter2)
    ancgpf=ancgpf/(Nf*iter2)
    ancgpv=ancgpv/(Nf*itot)
    ancgf=ancgf/(Nf*itot)
    ancgfpv=ancgfpv/(Nf*itot)
    ancghpv=ancghpv/(2*Nf*iter2)
    ancgpfpv=ancgpfpv/(iter2*Nf)
    ancgma=ancgma/ipbp
    yav=yav/iter2
    yyav=yyav/iter2-yav*yav 
    yyav=sqrt(yyav/(iter2-1)) 
    atraj=dt*itot/iter2 
!*******************************************************************
!     print global averages
!*******************************************************************
    if (ip_global .eq. 0) then
!     write(6, 9022) iter2,naccp,atraj,yav,yyav,ancg,ancgpv,ancgh,ancghpv,ancgf,
!    & ancgfpv,ancgpf,ancgpfpv,pbpa,vel2a,actiona
       write(7, 9022) iter2,naccp,atraj,yav,yyav, &
            & ancg,ancgpv,ancgh,ancghpv,ancgf,ancgfpv,ancgpf,ancgpfpv, &
            & pbpa,ancgma,vel2a,actiona
9022   format(' averages for last ',i6,' trajectories',/  &
            & 1x,' # of acceptances: ',i6,' average trajectory length= ',f8.3/ &
            & 1x,' <exp-dH>=',e11.4,' +/-',e10.3/ &
            & 1x,' av. # QMR itr.'/ &
            & 1x,'     guidance: DWF  ',f9.3,'; PV  ',f9.3/ &
            & 1x,'   acceptance: DWF  ',f9.3,'; PV  ',f9.3/ &
            & 1x,'        force: DWF  ',f9.3,'; PV  ',f9.3/ &
            & 1x,'pseudofermion: DWF  ',f9.3,'; PV  ',f9.3/ &
            & 1x,' psibarpsi=',e11.3/ &
            & 1x,' av. # QMR itr.',f9.3// &
            & 1x,' mean square velocity=',e10.3,'; action per site=',e10.3//)
       write(7, 9024)
       write(7, 9024)
9024   format(1x)
!
       close(11)
    end if
!
    if(iwrite.eq.1) then
       call rranget(seed, 1, 1, 1)
       call swrite
       write(7,*) 'seed: ', idum
    endif
!
  end subroutine dwf3d_main
!******************************************************************
!   calculate dSds for gauge fields at each intermediate time
!******************************************************************
  subroutine force(Phi,res1,am,imass,isweep,iter)
    use remezg
    use trial
    use vector, X1=>X
    use gforce
    use param
    use comms

    complex(dp), intent(in) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4, Nf)
    real, intent(in) :: res1, am
    integer, intent(in) :: imass, isweep, iter
!     complex Phi(kferm,Nf),X2(kferm)
!     complex X1,u
    complex(dp) :: X2(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer :: ia, itercg
!
!     write(6,111)
!111 format(' Hi from force')
!
    dSdpi = 0.0
!
! uncomment this line to quench the fermions!
!     return
! pseudofermion action is
!   Phi^dagger {MdaggerM(1)}^1/4 {MdaggerM(m)})^-1/2 {MdaggerM(1)}^1/4 Phi
!
    do ia=1,Nf
!
       X2 = Phi(:, :, :, :, :, ia)

       call qmrherm(X2, res1, itercg, One, 1, anum4g, aden4g, ndiagg, 1, &
            & isweep, iter)
       ancgpv=ancgpv+float(itercg)

       X2 = X1
!
       call qmrherm(X2, res1, itercg, am, imass, bnum2g, bden2g, ndiagg, 0, &
            & isweep, iter)
       ancg=ancg+float(itercg)
!     write(111,*) itercg
       X2 = X1
!
!  evaluates -X2dagger * d/dpi[{MdaggerM(m)}^1/2] * X2
       call qmrherm(X2, res1, itercg, am, imass, anum2g, aden2g, ndiagg, 2, &
            & isweep, iter)
       ancgf=ancgf+float(itercg)

!     write(113,*) itercg
!  evaluates +2Re{Phidagger * d/dpi[{MdaggerM(1)}^1/4] * X2}
       call qmrherm(X2, res1, itercg, One, 1, anum4g, aden4g, ndiagg, 3, &
            & isweep, iter)
       ancgfpv=ancgfpv+float(itercg)
!
    enddo
!
    if (ibound.eq.-1 .and. ip_t.eq.(np_t - 1)) then
       dSdpi(:, :, ksizet_l, 3) = -dSdpi(:, :, ksizet_l, 3)
    endif
!
    dSdpi = dSdpi + beta * Nf * theta
!
    return
  end subroutine force
!******************************************************************
!   Evaluation of Hamiltonian function
!******************************************************************
  subroutine hamilton(Phi, h, hg, hp, s, res2, isweep, iflag, am, imass)
    use remez
    use trial, only: theta, pp
    use vector, X1=>X
    use dum1
    use param
    use comms
    complex(dp), intent(in) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4, Nf)
    real(dp), intent(out) :: h, hg, hp, s
    real, intent(in) :: res2, am
    integer, intent(in) :: isweep, iflag, imass
!     complex, intent(in) :: Phi(kthird, ksizex_l, ksizey_l, ksizet_l, 4, Nf)
!     complex X1,R
    real(dp) :: hf
    integer :: itercg, ia
!     write(6,111)
!111 format(' Hi from hamilton')
!
    hf=0.0
!
    hp = 0.5 * sum(pp ** 2)
#ifdef MPI
    call MPI_AllReduce(MPI_In_Place, hp, 1, MPI_Double_Precision, MPI_Sum, comm)
#endif

    hg = 0.5 * Nf * beta * sum(theta ** 2)
#ifdef MPI
    call MPI_AllReduce(MPI_In_Place, hg, 1, MPI_Double_Precision, MPI_Sum, comm)
#endif
    h = hg + hp
! 
! uncomment these lines to quench the fermions!
!     write(6,*) isweep,':  hg', hg,'   hp', hp,'   h',h
!     return
!         
!  pseudofermion action is
!   Phi^dagger {MdaggerM(1)}^1/4 {MdaggerM(m)})^-1/2 {MdaggerM(1)}^1/4 Phi
!
    do ia = 1,Nf
!
       R = Phi(:, :, :, :, :, ia)

       call qmrherm(R, res2, itercg, One, 1, anum4, aden4, ndiag, 0, isweep, iflag)
       ancghpv=ancghpv+float(itercg)
!
       R = X1
!
       call qmrherm(R, res2, itercg, am, imass, bnum2, bden2, ndiag, 0, isweep, iflag)
       ancgh=ancgh+float(itercg)
!
       hf = hf + sum(real(conjg(R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)) &
       &        * X1(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)))
!
    enddo
#ifdef MPI
! hf is built up from zero during the loop so only needs to be summed across
! all partitions at this point
    call MPI_AllReduce(MPI_In_Place, hf, 1, MPI_Double_Precision, MPI_Sum, comm)
#endif
!
    h = hg + hp + hf
!     write(6,*) isweep,':  hg', hg,'   hp', hp,'   hf', hf, '   h',h
    s = hg + hf
!
    return
  end subroutine hamilton

!******************************************************************
!    multisolver matrix inversion via Lanczos technique
!  eg. Golub & van Loan "Matrix Computations" 9.3.1
!       solves (MdaggerM+diag)*x=Phi for ndiag different values of diag
!   iflag=0: simply evaluates X = {MdaggerM}^p * Phi
!   can be interchanged with congrad for p=-1
!   iflag=1: in addition updates Phi0 register needed for PV force term
!   iflag=2: evaluates DWF force term
!   iflag=3: evaluates PV force term
!*****************************************************************m
  subroutine qmrherm(Phi, res, itercg, am, imass, anum, aden, ndiagq, iflag, isweep, &
       & iter)
    use trial, only: u
    use vector
    use gforce
    use phizero
    use qmrherm_scratch
    use comms
    complex(dp), intent(in) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass, ndiagq, iflag, isweep, iter
    real(dp), intent(in) :: anum(0:ndiagq), aden(ndiagq)
    real, intent(in) :: res, am
    integer, intent(out) :: itercg
!
    real :: alphatild
    real(dp) :: coeff
!      
    real(dp) :: alpha(ndiagq)
    real(dp) :: amu(ndiagq), d(ndiagq), dm1(ndiagq)
    real(dp) :: rho(ndiagq), rhom1(ndiagq)
    real(dp) :: betaq, betaq0, phimod
    real :: resid, rhomax, arelax
    integer :: niter, idiag
#ifdef MPI
    type(MPI_Request), dimension(12) :: reqs_X2, reqs_vtild, reqs_Phi0, reqs_R, reqs_x
#else
    integer :: reqs_X2, reqs_vtild, reqs_Phi0, reqs_R, reqs_x
#endif
!
!     write(6,111)
!111 format(' Hi from qmrherm')
!
    resid=sqrt(kthird*ksize*ksize*ksizet*4*res*res)
!     write(6,*) iflag, resid
    itercg=0
!
!   initialise r=Phi
!
    R = Phi
    qm1 = cmplx(0.0, 0.0)
    x = anum(0) * Phi

    betaq = sum(abs(R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)) ** 2)
#ifdef MPI
    call MPI_AllReduce(MPI_In_Place, betaq, 1, MPI_Double_Precision, MPI_Sum, comm)
#endif
    betaq = sqrt(betaq)
    phimod=betaq
!     write(6,*) '|| Phi || = ', phimod
!
    do niter=1,max_qmr_iters
       itercg=itercg+1
!
!  Lanczos steps
!
       q = R / betaq

       call dslash(vtild,q,u,am,imass)
#ifdef MPI
! No way to hide communications here unfortunately
       call start_halo_update_5(4, vtild, 1, reqs_vtild)
       call complete_halo_update(reqs_vtild)
#else
       call update_halo_5(4, vtild)
#endif

       call dslashd(x3,vtild,u,am,imass)
!
       alphatild = sum(real(conjg(q(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)) & 
       &                * x3(:,1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)))
#ifdef MPI
       call MPI_AllReduce(MPI_In_Place, alphatild, 1, MPI_Real, MPI_Sum, comm)
#endif
!
       R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = &
            & x3(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) &
            & - alphatild * q(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) &
            & - betaq * qm1
#ifdef MPI
! R will be needed at the start of the next iteration to compute q
! so start updating the boundary
       call start_halo_update_5(4, R, 2, reqs_R)
#else
       call update_halo_5(4, R)
#endif
       qm1 = q(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
!
       betaq0 = betaq
       betaq = sum(abs(R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l,:)) ** 2)
#ifdef MPI
       call MPI_AllReduce(MPI_In_Place, betaq, 1, MPI_Double_Precision, MPI_Sum, comm)
#endif
       betaq = sqrt(betaq)
!
       alpha = alphatild + aden
!
       if(niter.eq.1)then
          d = alpha
          rho = betaq0 / alpha
          rhom1 = rho
          do idiag = 1, ndiagq
             p(:, :, :, :, :, idiag) = q(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
             x1(:, :, :, :, :, idiag) = rho(idiag) * p(:, :, :, :, :, idiag)
          enddo
          pm1 = p
       else
          amu = betaq0 / d
          dm1 = d
          d = alpha - betaq0 * amu
          rho = -amu * dm1 * rhom1 / d
          do idiag = 1, ndiagq
             p(:, :, :, :, :, idiag) = q(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) &
                  & - amu(idiag) * pm1(:, :, :, :, :, idiag)
          enddo
          pm1 = p
!     Convergence criterion (a bit ad hoc for now...)
          rhomax = real(maxval(abs(phimod * rho)))
          rhom1 = rho
          do idiag = 1, ndiagq
             x1(:, :, :, :, :, idiag) = &
                  & x1(:, :, :, :, :, idiag) &
                  & + rho(idiag) * p(:, :, :, :, :, idiag)
          enddo

!     check to see whether the residual is acceptable for all ndiagq....
!     criterion is a bit ad hoc -- relaxing by a factor arelax improves code
!     stability and leads to quicker convergence
          arelax=2.0
          if(rhomax .lt. arelax * resid) then
!     if(rhomax.lt.resid) then
!     call testinv(Phi,resmax,itercg,am,imass,x1,aden,ndiagq)
!     convergence based on || residual || not working well in single precision...
!     if(resmax.lt.resid) goto 8
             exit
          endif
       endif
!     
!     end of loop over iter
#ifdef MPI
! R will be needed at the start of the next iteration to compute q
! so start updating the bounddary
       call complete_halo_update(reqs_R)
#endif
    enddo
    if (niter .gt. max_qmr_iters .and. ip_global .eq. 0) then
       write(7,*) 'QMRniterc!, isweep,iter,iflag,imass,anum,ndiagq = ' &
       &        ,isweep, iter, iflag, imass, anum(0), ndiagq
    end if
!     
    if(iflag.lt.2)then
!     Now evaluate solution x=(MdaggerM)^p * Phi
       do idiag=1,ndiagq
          x(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = &
               & x(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) &
               & + anum(idiag) * x1(:, :, :, :, :, idiag)
       enddo
#ifdef MPI
! x is a saved module variable, so must be updated to avoid polluting the parent function
! could this in principle be moved outside the function so we don't do it unnecessarily?
! but in that case we wouldn't be able to hide the communications
       call start_halo_update_5(4, x, 3, reqs_x)
#else
       call update_halo_5(4, x)
#endif
!     
!  update phi0 block if required...
       if(iflag.eq.1) then
          Phi0(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :, 1:ndiagq) = &
               & X1(:, :, :, :, :, 1:ndiagq)
#ifdef MPI
! No way to hide communications here unfortunately
! In principle this could be better interleaved with the x update
! but that would add extra branching, and this section is messy enough already
          call start_halo_update_6(4, ndiagq, Phi0, 4, reqs_Phi0)
          call complete_halo_update(reqs_Phi0)
#else
          call update_halo_6(4, ndiagq, Phi0)
#endif
       endif
#ifdef MPI
       call complete_halo_update(reqs_x)
#endif
!     
    else
!
       do idiag=1, ndiagq
!
!  X2 = M*X1
          R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = X1(:, :, :, :, :, idiag)
#ifdef MPI
! No way to hide communications here unfortunately
          call start_halo_update_5(4, R, 5, reqs_R)
          call complete_halo_update(reqs_R)
#else
          call update_halo_5(4, R)
#endif

! Communication of X2 generated here can be hidden if iflag isn't 2, while R is updated
          call dslash(X2, R, u, am, imass)
#ifdef MPI
          call start_halo_update_5(4, X2, 6, reqs_X2)
#else
          call update_halo_5(4, X2)
#endif
!
          if(iflag.eq.2)then
             coeff=anum(idiag)
#ifdef MPI
             call complete_halo_update(reqs_X2)
#endif
             call derivs(R, X2, coeff, 0)
          else
             coeff=-anum(idiag)
             R = Phi0(:, :, :, :, :, idiag)
#ifdef MPI
             call complete_halo_update(reqs_X2)
#endif
             call derivs(R, X2, coeff, 0)
!
! Communication of X2 generated here can be hidden while R is updated
             call dslash(X2, R, u, am, imass)
#ifdef MPI
             call start_halo_update_5(4, X2, 7, reqs_X2)
#else
             call update_halo_5(4, X2)
#endif
!
             R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = x1(: ,:, :, :, :, idiag)
#ifdef MPI
             call start_halo_update_5(4, R, 8, reqs_R)
             call complete_halo_update(reqs_X2)
             call complete_halo_update(reqs_R)
#else
             call update_halo_5(4, R)
#endif
             call derivs(X2, R, coeff, 1)
          endif
!
       enddo
    endif
!
!
    return
  end subroutine qmrherm
!**********************************************************************
!  iflag = 0 : evaluates Rdagger*(Mdagger)'*X2
!  iflag = 1 : evaluates Rdagger*(M)'*X2
!**********************************************************************
  subroutine derivs(R,X2,anum,iflag)
    use gforce
    use dirac
!      complex, intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
!      complex, intent(in) :: X2(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: X2(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real(dp), intent(in) :: anum
    integer, intent(in) :: iflag

!      complex(dp) :: tzi
    real :: tzi_real
    integer :: ix, iy, it, ixup, iyup, itup, idirac, mu
    integer :: igork1
!
!     write(6,111)
!111 format(' Hi from derivs')

!     dSdpi=dSdpi-Re(Rdagger *(d(Mdagger)dp)* X2)
!     Cf. Montvay & Muenster (7.215)
!      tzi=cmplx(0.0,2*anum)
    tzi_real = 2 * real(anum)
!     factor of 2 picks up second term in M&M (7.215)
!
    do mu = 1,3
       ixup = kdelta(1, mu)
       iyup = kdelta(2, mu)
       itup = kdelta(3, mu)

       do idirac=1,4
          do it = 1,ksizet_l
             do iy = 1,ksizey_l
                do ix = 1,ksizex_l
                   dSdpi(ix,iy,it,mu) = &
                   &     dSdpi(ix,iy,it,mu) + tzi_real * real(akappa) * sum(aimag( &
                   &       conjg(R(:,ix,iy,it,idirac)) * &
                   &         X2(:,ix+ixup,iy+iyup,it+itup,idirac)) &
                   &     - aimag(conjg(R(:,ix+ixup,iy+iyup,it+itup,idirac)) * &
                   &         X2(:,ix,iy,it,idirac)))
                enddo
             enddo
          enddo
!
          igork1=gamin(mu,idirac)
          if(iflag.eq.0)then
             do it = 1,ksizet_l
                do iy = 1,ksizey_l
                   do ix = 1,ksizex_l
                      dSdpi(ix,iy,it,mu) = &
                      &     dSdpi(ix,iy,it,mu)+ tzi_real * sum(aimag(gamval(mu,idirac)* &
                      &(conjg(R(:,ix,iy,it,idirac)) * &
                      &        X2(:, ix+ixup,iy+iyup,it+itup,igork1) &
                      &+conjg(R(:,ix+ixup,iy+iyup,it+itup,idirac))* &
                      &             X2(:,ix,iy,it,igork1))))
                   enddo
                enddo
             enddo
          else
             do it = 1,ksizet_l
                do iy = 1,ksizey_l
                   do ix = 1,ksizex_l
                      dSdpi(ix,iy,it,mu) = &
                      &     dSdpi(ix,iy,it,mu)- tzi_real * sum(aimag(gamval(mu,idirac)* &
                      &(conjg(R(:,ix,iy,it,idirac)) * &
                      &        X2(:,ix+ixup,iy+iyup,it+itup,igork1) &
                      &+conjg(R(:,ix+ixup,iy+iyup,it+itup,idirac)) * &
                      &             X2(:,ix,iy,it,igork1))))
                   enddo
                enddo
             enddo
          endif
!
       enddo
    enddo
!
    return
  end subroutine derivs
! c******************************************************************
! c   Calculates residual for testing purposes....
! c   needs to run with double precision vectors to be useful.....
! c******************************************************************
!       subroutine testinv(Phi,resmax,itercg,am,imass,x,aden,ndiag)
!       parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
!       parameter(kferm=4*kthird*kvol)
!       common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
!       common/para/bbb,am3,ibound
!       complex Phi(kferm)
!       complex(dp) vtild(kferm)
!       complex(dp) x3(kferm)
!       complex(dp) x(kferm,ndiag),x1(kferm),x2(kferm)
!       complex(dp) u
! c     complex vtild(kferm)
! c     complex x3(kferm)
! c     complex x(kferm,ndiag),x1(kferm),x2(kferm)
! c      complex u
!       real(dp) residual
!       real(dp) aden(ndiag)
! c
!       write(6,111)
! 111   format(' Hi from testinv')
! c
!       resmax=0.0
! c
!       do idiag=1,ndiag
!       residual=0.0
!       do i=1,kferm
!       x3(i)=x(i,idiag)
!       enddo
!       call dslash(x2,x3,u,am,imass)
!       call dslashd(x1,x2,u,am,imass)
!       do i=1,kferm
!       vtild(i)=x1(i)+aden(idiag)*x3(i)-Phi(i)
!       residual=residual+conjg(vtild(i))*vtild(i)
!       enddo
! c     residual=sqrt(residual)
!       if(residual.gt.resmax) resmax=residual
! c
!       write(6,*) idiag, 'itercg = ',itercg, ' residual = ',residual
!       enddo
! c
!       resmax=sqrt(resmax)
! c
!       return
!       end
!******************************************************************
!    matrix inversion via conjugate gradient algorithm
!       solves (Mdagger)Mx=Phi, 
!           NB. no even/odd partitioning
!******************************************************************
  subroutine congrad(Phi,res,itercg,am,imass)
    use trial, only: u
    use vector
    use comms
    complex(dp), intent(in) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
!     complex, intent(in) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real, intent(in) :: res, am
    integer, intent(out) :: itercg
    integer, intent(in) :: imass

    integer, parameter :: niterc=kthird*kvol
!     complex x1(kferm),x2(kferm),p(kferm),r(kferm)
    complex(dp) :: x1(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: x2(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: p(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: r(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real :: resid
    real :: betacg, betacgn, betacgd, alpha, alphan, alphad
    integer :: nx
#ifdef MPI
    type(MPI_Request), dimension(12) :: reqs_x1, reqs_r
#else
    integer :: reqs_x1, reqs_r
#endif
!     write(6,111)
!111 format(' Hi from congrad')
!
    resid = 4 * ksize * ksize * ksizet * kthird * res * res
    itercg = 0
    alphan = 0.0
!   initialise p=x, r=Phi(na)
    p = x
    r = Phi
!
    betacgd=1.0
    alpha=1.0
!
    do nx=1,niterc
       itercg=itercg+1
!
!  x1=Mp
       call dslash(x1,p,u,am,imass)
#ifdef MPI
       call start_halo_update_5(4, x1, 8, reqs_x1)
#endif
!
       if(nx.ne.1)then
!
!   alpha=(r,r)/(p,(Mdagger)Mp)
!   Don't need x1's halo at this point
          alphad = real(sum(abs(x1(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)) ** 2))
#ifdef MPI
          call MPI_AllReduce(MPI_In_Place, alphad, 1, MPI_Real, MPI_Sum, comm)
#endif       
          alpha = alphan / alphad
!     
!   x=x+alpha*p
          x = x + alpha * p
       end if

!   Now we need x1's halo, so ensure communication is finished
#ifdef MPI
       call complete_halo_update(reqs_x1)
#else
       call update_halo_5(4, x1)
#endif
!     
!   x2=(Mdagger)x1=(Mdagger)Mp
       call dslashd(x2, x1, u, am, imass)
!
!   r=r-alpha*(Mdagger)Mp
!   Use x2 with wrong halo since we can't hide communications here
       r = r - alpha * x2

!   Now update halo for r instead since we can hide communication during the summation
!   x2 is discarded so we no longer care about its halo
#ifdef MPI
       call start_halo_update_5(4, r, 10, reqs_r)
#endif

!   betacg=(r_k+1,r_k+1)/(r_k,r_k)
       betacgn = real(sum(abs(r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)) ** 2))
#ifdef MPI
       call MPI_AllReduce(MPI_In_Place, betacgn, 1, MPI_Real, MPI_Sum, comm)
#endif       
       betacg = betacgn / betacgd
       betacgd = betacgn
       alphan = betacgn
!
       if(nx.eq.1) betacg=0.0
!
!   p=r+betacg*p
!   Now the correct value of r is needed to avoid having to communicate p as well
#ifdef MPI
       call complete_halo_update(reqs_r)
#else
       call update_halo_5(4, r)
#endif
       p = r + betacg * p
       if(betacgn.lt.resid) exit
    end do
!     write(6,1000)
    if (nx.gt.niterc .and. ip_global.eq.0) then
       write(7,1000)
1000   format(' # iterations of congrad exceeds niterc')
    end if
!    print *, nx
!     write(6,*) itercg
    return
  end subroutine congrad

!*****************************************************************
!   Calculate fermion expectation values via a noisy estimator
!   -matrix inversion via conjugate gradient algorithm
!       solves Mx=x1
!     (Numerical Recipes section 2.10 pp.70-73)   
!*******************************************************************
  subroutine measure(psibarpsi, res, aviter, am, imass)
    use trial, only: u
    use vector, xi=>x
    use comms
    real, intent(out) :: psibarpsi, aviter
    real, intent(in) :: res, am
    integer, intent(in) :: imass
    integer, parameter :: knoise = 10
!     complex :: x(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
!     complex :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
!     complex :: psibarpsi1,psibarpsi2
    complex(dp) :: x(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: psibarpsi1,psibarpsi2
    real(dp) :: cnum(0:1),cden(1)
    real :: ps(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 2)
    real :: pt(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 2)
    real(dp) :: pbp(knoise)
    integer :: idsource, idsource2, inoise
    integer :: iter, itercg
    real :: susclsing
#ifdef MPI
    type(MPI_Request), dimension(12) :: reqs_ps, reqs_pt, reqs_Phi
#else
    integer :: reqs_ps, reqs_pt, reqs_Phi
#endif
!     write(6,*) 'hi from measure'
!
    iter=0
!     pbp=0.0
    cnum(0)=0.0
    cnum(1)=1.0
    cden(1)=0.0
!
    do inoise=1,knoise
!
!     set up noise
       call gauss0(ps, reqs_ps)
       psibarpsi1=(0.0,0.0)
       call gauss0(pt, reqs_pt)
       psibarpsi2=(0.0,0.0)
!
       do idsource=1,2
!
!  source on domain wall at ithird=1
!   
          x = cmplx(0.0, 0.0)
#ifdef MPI
!  We started a halo update in the gauss0 call;
!  Now we need it to be complete if it isn't already
          if (idsource .eq. 1) then
             call complete_halo_update(reqs_ps)
          end if
#endif
          if(imass.ne.5)then
             x(:, :, :, idsource) = cmplx(ps(:,:,:,1), ps(:,:,:,2))
          else
             x(:, :, :, idsource+2) = cmplx(ps(:,:,:,1), ps(:,:,:,2))
          endif
!
          xi = cmplx(0.0, 0.0)
          if(imass.ne.5)then
             xi(1, :, :, :, :) = x
          else
!     xi = 0.5(1+gamma_4)*gamma_5*eta on DW at ithird=1
             xi(1, :, :, :, 1:2) = -x(:, :, :, 3:4)
!     xi = 0.5(1+gamma_4)*eta on DW at ithird=1
          endif
!
! Phi= Mdagger*xi
!
          call dslashd(Phi, xi, u, am, imass)
#ifdef MPI
! No way to hide communications here unfortunately
          call start_halo_update_5(4, Phi, 11, reqs_Phi)
          call complete_halo_update(reqs_Phi)
#else
          call update_halo_5(4, Phi)
#endif
!     call qmrherm(Phi,res,itercg,am,imass,cnum,cden,1,0)
          call congrad(Phi, res, itercg, am, imass)
          iter = iter + itercg
!
          if(imass.ne.5)then
!     pbp1 = x^dagger (0.5(1+gamma_4)) xi(kthird)
             psibarpsi1=psibarpsi1 &
             &           + sum(conjg(x(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource)) * &
             &               xi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource))
          else
!     pbp1 = x^dagger (0.5(1-gamma_4)) xi(1)
             psibarpsi1=psibarpsi1 &
             &           + sum(conjg(x(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource+2)) &
             &                 * xi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource+2))
          endif
!
!
!  source on domain wall at ithird=kthird
          idsource2=idsource+2
!
          x = cmplx(0.0, 0.0)
#ifdef MPI
!  Again, if this isn't finished by now we have to wait for it
          if (idsource .eq. 1) then
             call complete_halo_update(reqs_pt)
          end if
#endif
          if(imass.ne.5)then
             x(:, :, :, idsource2) = cmplx(pt(:,:,:,1), pt(:,:,:,2))
          else
             x(:, :, :, idsource2 - 2) = cmplx(pt(:,:,:,1), pt(:,:,:,2))
          endif
!
          xi = cmplx(0.0, 0.0)
          if(imass.ne.5)then
!   xi = 0.5(1-gamma_4)*eta on DW at ithird=kthird
             xi(kthird, :, :, :, :) = x
          else
!   xi = 0.5(1-gamma_4)*gamma_5*eta on DW at ithird=kthird
             xi(kthird, :, :, :, 3:4) = -x(:, :, :, 1:2)
          endif
!     
! Phi= Mdagger*xi
!
          call dslashd(Phi,xi,u,am,imass)
#ifdef MPI
! No way to hide communications here unfortunately
          call start_halo_update_5(4, Phi, 12, reqs_Phi)
          call complete_halo_update(reqs_Phi)
#else
          call update_halo_5(4, Phi)
#endif
!
! xi= (M)**-1 * Phi
!
!     call qmrherm(Phi,res,itercg,am,imass,cnum,cden,1,0)
          call congrad(Phi,res,itercg,am,imass)
          iter = iter + itercg
!     
          if(imass.ne.5)then
! pbp2= - x^dagger (0.5(1-gamma_4)) xi(1)
             psibarpsi2=psibarpsi2 &
             &           +sum(conjg(x(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource2)) &
             &                * xi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource2))
          else
! pbp2= - x^dagger (0.5(1-gamma_4)) xi(kthird)
             psibarpsi2=psibarpsi2 &
             &           +sum(conjg(x(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource)) &
             &           * xi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idsource))
          endif
!
!  end trace on Dirac indices....
       enddo
#ifdef MPI
!  The sums psibarpsi1 and psibarpsi2 are initialised to 0 outside the loop
!  So can be summed up per process, then collected here at the end
       call MPI_AllReduce(MPI_In_Place, psibarpsi1, 1, MPI_Double_Complex, &
            & MPI_Sum, comm)
       call MPI_AllReduce(MPI_In_Place, psibarpsi2, 1, MPI_Double_Complex, &
            & MPI_Sum, comm)
#endif
!
       if(imass.eq.1)then
          psibarpsi1 = psibarpsi1 / kvol
          psibarpsi2 = psibarpsi2 / kvol
          pbp(inoise) = psibarpsi1 + psibarpsi2
       elseif(imass.eq.3)then
          psibarpsi1 = cmplx(0.0,-1.0) * psibarpsi1 / kvol
          psibarpsi2 = cmplx(0.0,+1.0) * psibarpsi2 / kvol
          pbp(inoise) = psibarpsi1 + psibarpsi2
       elseif(imass.eq.5)then
          psibarpsi1 = cmplx(0.0,-1.0) * psibarpsi1 / kvol
          psibarpsi2 = cmplx(0.0,-1.0) * psibarpsi2 / kvol
          pbp(inoise) = psibarpsi1 + psibarpsi2
       endif
!        write(6,*) real(psibarpsi1),aimag(psibarpsi1), real(psibarpsi2),aimag(psibarpsi2)
       if (ip_global .eq. 0) then
          write(100,*) real(psibarpsi1), aimag(psibarpsi1), &
               & real(psibarpsi2), aimag(psibarpsi2)
       end if
!
! end loop on noise
    enddo
!
    susclsing=0.0
!
    psibarpsi = real(sum(pbp))
    do inoise=1,knoise
       susclsing = susclsing + real(sum(pbp(inoise) * pbp(inoise+1:knoise)))
    enddo
    psibarpsi = psibarpsi / knoise
    susclsing = 2 * kvol * susclsing / (knoise * (knoise-1))
    if (ip_global .eq. 0) then
       write(200,*) psibarpsi, susclsing
    end if
    aviter = float(iter) / (4*knoise)
    return
  end subroutine measure

!******************************************************************
!   Calculate meson correlators using point sources on domain walls
!   -matrix inversion via conjugate gradient algorithm
!       solves Mx=x1
!     (Numerical Recipes section 2.10 pp.70-73)   
!*******************************************************************
  subroutine meson(res,itercg,aviter,am,imass)
    use random
    use vector, xi=>x
    use dirac
    use trial
    
    real, intent(in) :: res, am
    integer, intent(out) :: itercg
    real, intent(out) :: aviter
    integer, intent(in) :: imass
    !     complex x(kvol,4),x0(kvol,4),Phi(kthird,kvol,4)
    !     complex xi,gamval
    !     complex prop00(kvol,3:4,1:2),prop0L(kvol,3:4,3:4)
    complex(dp) :: x(ksizex_l, ksizey_l, ksizet_l, 4),
    complex(dp) :: x0(ksizex_l, ksizey_l, ksizet_l, 4)
    complex(dp) :: Phi(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
    complex(dp) :: prop00(ksizex_l, ksizey_l, ksizet_l, 3:4, 1:2)
    complex(dp) :: prop0L(ksizex_l, ksizey_l, ksizet_l, 3:4, 3:4)
!    complex :: prop00n(ksizex_l, ksizey_l, ksizet_l, 3:4, 1:2)
!    complex :: prop0Ln(ksizex_l, ksizey_l, ksizet_l, 3:4, 3:4)
    real :: cpm(0:ksizet-1), cmm(0:ksizet-1)
!    complex :: cpmn(0:ksizet-1),cmmn(0:ksizet-1)
!    complex :: cferm1(0:ksizet-1), cferm2(0:ksizet-1)
    complex(dp) :: cferm1(0:ksizet-1), cferm2(0:ksizet-1)
    real :: ps(ksizex_l, ksizey_l, ksizet_l, 2)
    real :: c, chim, chip
    real :: ix, iy, it
    integer :: ip_xxx, ip_yyy, ip_ttt, ixxx_l, iyyy_l, ittt_l
    !     write(6,*) 'hi from meson'
    !      
    nsource=5
    nsmear=10
    c=0.25
    iter=0
!
    cpm = (0.0,0.0)
    cmm = (0.0,0.0)
!    cpmn = (0.0,0.0)
!    cmmn = (0.0,0.0)
    cferm1 = (0.0,0.0)
    cferm2 = (0.0,0.0)
    !
    !      susceptibility
    chim = 0.0
    chip = 0.0
    !
    do ksource=1,nsource
       !
       !   random location for +m source
       ip_xxx = int(np_x * rano(yran, idum, 1, 1, 1))
       ip_yyy = int(np_y * rano(yran, idum, 1, 1, 1))
       ip_ttt = int(np_t * rano(yran, idum, 1, 1, 1))
       ixxx_l = int(ksizex_l * rano(yran, idum, 1, 1, 1)) + 1
       iyyy_l = int(ksizey_l * rano(yran, idum, 1, 1, 1)) + 1
       ittt_l = int(ksizet_l * rano(yran, idum, 1, 1, 1)) + 1
!!!! COMMUNICATE
!       write(6,*) ixxx,iyyy,ittt, isource
       !
       do idsource=3,4
          !  source on domain wall at ithird=1
          xi = (0.0, 0.0)
          x = (0.0,0.0)
          !  wall source
          !  if (ip_t .eq. np_t) then
          !    x(:, :, ksizet_l, :) = cmplx(1.0,0.0) / ksize2
          !  end if
          !  point source at fixed site, spin...
          if (ip_xxx .eq. ip_x .and. ip_yyy .eq. ip_y .and. ip_ttt .eq. ip_t) then
             x(ixxx_l, iyyy_l, ittt_l, idsource) = cmplx(1.0,0.0)
          end if
          !
          ! now smear it.....
          !
          do ismear=1,nsmear
             call dslash2d(x0, x, u)
!!!! comms
             x = (1.0-c) * x + c * x0
          enddo
          !
          !
          !   xi = x  on DW at ithird=1
          !
          xi(1, :, :, :, :) = x
          !
          ! Phi= Mdagger*xi
          !
          call dslashd(Phi, xi, u, am, imass)
!!!! comms
          !
          !  preconditioning (no,really)
          call dslashd(xi, Phi, u, am, imass)
!!!! comms
          !  
          ! xi= (MdaggerM)**-1 * Phi 
          !
          ! call congrad(Phi,res,itercg,am,imass) 
          iter=iter+itercg
          !
          prop00(:, :, :, idsource, 1:2) = xi(1, :, :, :, 1:2)
          prop0L(:, :, :, idsource, 3:4) = xi(kthird, :, :, :, 3:4)
          !
          ! if(imass.ne.1)then
          !  now evaluate with sign of mass reversed (not needed for hermitian mass term)
          ! am=-am
          !  source on domain wall at ithird=1
          ! xi = (0.0,0.0)
          ! if (ip_t .eq. np_t) then
          !   xi(1, :, :, :, :) = x
          ! end if
          !
          ! Phi= Mdagger*xi
          !
          ! call dslashd(Phi,xi,u,am,imass)
!!!!COMMS
          !
          ! call dslashd(xi,Phi,u,am,imass)
          !  
!!!! COMMS
          ! xi= (MdaggerM)**-1 * Phi 
          !
          ! call congrad(Phi,res,itercg,am,imass) 
          ! iter=iter+itercg
          !
          ! prop00n(:, :, :, idsource, 1:2) = xi(1, :, :, :, 1:2)
          ! prop0Ln(:, :, :, idsource, 3:4) = xi(kthird, :, :, :, 3:4)
          !
          ! am=-am
          !
          !  end loop on source Dirac index....
       enddo
       !
       !  Now tie up the ends....
       !
       !  First C+-
       !
       !  now evaluate the trace (exploiting projection)
!!!! SOME MESSY COMMS HERE
       do it=0,ksizet-1
          itt=mod((ittt+it-1),ksizet)+1
          cpm(it) = cpm(it) + &
               & sum(prop00(:, :, itt, 3:4, 1:2) * conjg(prop00(:, :, itt, 3:4, 1:2)))
       enddo
       !
       ! if(imass.ne.1)then
       !     do it=0,ksizet-1
       !         itt=mod((ittt+it-1),ksizet)+1
       !         cpmn(it)=cpmn(it) &
       !        & + sum(prop00(:, :, itt, 3:4, 1:2)*conjg(prop00n(:, :, itt, 3:4, 1:2)))
       !     enddo
       ! endif
       !
       !  next C--
       !  now evaluate the trace exploiting projection
       do it=0,ksizet-1
          itt=mod((ittt+it-1),ksizet)+1
          cmm(it) = cmm(it) + &
               & sum(prop0L(:, :, itt, 3:4, 3:4) * conjg(prop0L(:, :, itt, 3:4, 3:4)))
       enddo
          !
          !     if(imass.ne.1)then
          !     do id1=3,4
          !     do id2=3,4
          !     do it=0,ksizet-1
          !     itt=mod((ittt+it-1),ksizet)+1
          !     ioff=(itt-1)*ksize2
          !     do i=1,ksize2
          !     cmmn(it)=cmmn(it) &
          !    &   +prop0L(i+ioff,id1,id2)*conjg(prop0Ln(i+ioff,id1,id2))
          !     enddo
          !     enddo
          !     enddo
          !     enddo
          !     endif
          !
          !    now the fermion propagator
          !  = tr{ P_-*Psi(0,1)Psibar(x,Ls) + gamma_0*P_-*Psi(0,1)Psibar(x,1) }
          do idd=3,4
             do it=0,ksizet-1
                itt=mod((ittt+it-1),ksizet)+1
                ! correct for apbc
                if(itt.ge.ittt)then
                   isign=1
                else
                   isign=ibound
                endif
                !
                ioff=(itt-1)*ksize2
                do i=1,ksize2
                   cferm1(it)=cferm1(it) + &
                        & + isign * akappa * prop0L(i+ioff,idd,idd)
                   cferm2(it) = cferm2(it) + &
                        & isign * gamval(3,idd) * prop00(i+ioff,idd,gamin(3,idd))
                enddo
             enddo
          enddo
          !
          !  finish loop over sources
       enddo
       !
       do it=0,ksizet-1
          cpm(it)=cpm(it)/nsource
          cmm(it)=cmm(it)/nsource
          !  Cf. (54) of 1507.07717
          chim=chim+2*(cpm(it)+cmm(it))
       enddo
       !     if(imass.ne.1)then
       !       if(imass.eq.3)then
       !         do it=0,ksizet-1
       !           cpmn(it)=cpmn(it)/nsource
       !           cmmn(it)=cmmn(it)/nsource
       !  Cf. (54),(61) of 1507.07717
       !           chip=chip-2*(cpmn(it)-cmmn(it))
       !         enddo
       !       else
       !         do it=0,ksizet-1
       !           cpmn(it)=cpmn(it)/nsource
       !           cmmn(it)=cmmn(it)/nsource
       !  Cf. (64),(65) of 1507.07717
       !           chip=chip-2*(cpm(it)-cmm(it))
       !         enddo
       !       endif
       !     endif
       !
       do it=0,ksizet-1
          write(302,*) it, cpm(it), cmm(it)
          write(500,*) it, real(cferm1(it)), aimag(cferm1(it))
          write(501,*) it, real(cferm2(it)), aimag(cferm2(it))
       enddo
       !     write(6,*) chim
       write(400,*) chim
       !     if(imass.ne.1)then
       !     do it=0,ksizet-1
       !     write(402,*) it, real(cpmn(it)), real(cmmn(it))
       !     write(403,*) it, aimag(cpmn(it)), aimag(cmmn(it))
       !     enddo
       !     write(401,*) chip
       !     endif
       !
       !     if(imass.eq.1)then
       aviter=float(iter)/(2*nsource)
       !     else
       !     aviter=float(iter)/(4*nsource)
       !     endif
       !
       return
    end do
  end subroutine meson
!******************************************************************
!
  subroutine sread
    use random
    use gauge
#ifdef MPI
    use comms
    type(MPI_File) :: mpi_fh
    type(MPI_Status) :: status
    
    call MPI_File_Open(comm, 'con', MPI_Mode_Rdonly, &
         & MPI_Info_Null, mpi_fh)
! Get the configuration
    call MPI_File_Set_View(mpi_fh, 0_8, MPI_Real, mpiio_type, "native", &
         & MPI_Info_Null)
    call MPI_File_Read_All(mpi_fh, theta, 3 * ksizex_l * ksizey_l * ksizet_l, &
         & MPI_Real, status)
    call MPI_File_Close(mpi_fh)
! Get the seed
    if (ip_global.eq.0) then
       open(unit=10, file='con', status='old', form='unformatted', access='stream')
       call fseek(10, 3 * ksize * ksize * ksizet * 4 + 4, 0)
       read (10) seed
       close(10)
    end if
    call MPI_Bcast(seed, 1, MPI_Double_Precision, 0, comm)
#else
    open(unit=10, file='con', status='old', form='unformatted')
    read (10) theta, seed
    close(10)
#endif
    return
  end subroutine sread
!
  subroutine swrite
    use random
    use gauge
#ifdef MPI
    use comms
    type(MPI_File) :: mpi_fh
    type(MPI_Status) :: status
    
! Write theta
    call MPI_File_Open(comm, 'con', MPI_Mode_Wronly + MPI_Mode_Create, &
         & MPI_Info_Null, mpi_fh)
    call MPI_File_Set_View(mpi_fh, 0_8, MPI_Real, mpiio_type, "native", &
         & MPI_Info_Null)
    call MPI_File_Write_All(mpi_fh, theta, 3 * ksizex_l * ksizey_l * ksizet_l, &
         & MPI_Real, status)
    call MPI_File_Close(mpi_fh)

! Write seed in serial
    if (ip_global.eq.0) then
       open(unit=31, file='con', status='old', form='unformatted', access='stream')
       call fseek(31, 3 * ksize * ksize * ksizet * 4 + 4, 0)
! Manually compute the effective record length to be compatible with serial Fortran
       write (31), seed, 3 * ksize * ksize * ksizet * 4 + 8
       close(31)
    end if
#else
    open(unit=31, file='con', status='unknown', form='unformatted')
    write (31) theta, seed
    close(31)
#endif
    return
  end subroutine swrite
!
  subroutine init(nc)
    use random
    use gauge
    use dirac
!*******************************************************************
!     sets initial values
!     nc=0 cold start
!     nc=1 hot start
!     nc<0 no initialization
!*******************************************************************
    integer, intent(in) :: nc
!     complex one,zi
    complex(dp), parameter :: one=(1.0, 0.0), zi=(0.0, 1.0)
    integer :: ix, iy, it, mu
    real :: g
!
!*******************************************************************
!  calculate constants
!*******************************************************************
!      call addrc
!*******************************************************************
!    setup Dirac algebra
!*******************************************************************
!
!     gamma_1
!
    gamval(1,1)=-zi
    gamval(1,2)=-zi
    gamval(1,3)= zi
    gamval(1,4)= zi
!
    gamin(1,1)=4
    gamin(1,2)=3
    gamin(1,3)=2
    gamin(1,4)=1
!
!     gamma_2
!
    gamval(2,1)=-one
    gamval(2,2)= one
    gamval(2,3)= one
    gamval(2,4)=-one
!
    gamin(2,1)=4
    gamin(2,2)=3
    gamin(2,3)=2
    gamin(2,4)=1
!
!     gamma_3
!
    gamval(3,1)=-zi
    gamval(3,2)= zi
    gamval(3,3)= zi
    gamval(3,4)=-zi
!
    gamin(3,1)=3
    gamin(3,2)=4
    gamin(3,3)=1
    gamin(3,4)=2
!
!     gamma_4
!
    gamval(4,1)= one
    gamval(4,2)= one
    gamval(4,3)= -one
    gamval(4,4)= -one
!
    gamin(4,1)=1
    gamin(4,2)=2
    gamin(4,3)=3
    gamin(4,4)=4
!
!     gamma_5 = gamma_1 * gamma_2 * gamma_3 * gamma_4
!
    gamval(5,1)=-one
    gamval(5,2)=-one
    gamval(5,3)=-one
    gamval(5,4)=-one
!
    gamin(5,1)=3
    gamin(5,2)=4
    gamin(5,3)=1
    gamin(5,4)=2
!
!     gamma_4 * gamma_5 (called gamma_3 gamma_5 in notes)
    gamval(6,1)=-one
    gamval(6,2)=-one
    gamval(6,3)= one
    gamval(6,4)= one
!
    gamin(6,1)=3
    gamin(6,2)=4
    gamin(6,3)=1
    gamin(6,4)=2
!
!
    gamval = gamval * akappa
!
    if(nc.lt.0) return
!
!     initialize gauge fields
!
    if(nc .eq. 1)goto 40
!     (else cold start)
    theta = 0.0
    return
!
40  continue
    g=0.05
    do mu=1,3
       do it = 1, ksizet_l
          do iy = 1, ksizey_l
             do ix = 1, ksizex_l
!               theta(ix, iy, it, mu) = 2.0 * g * rranf(ix, iy, it) - 1.0
                theta(ix, iy, it, mu) = 2.0 * g * rano(yran,idum, ix, iy, it) - 1.0
             enddo
          enddo
       enddo
    enddo
    return
  end subroutine init
!******************************************************************
!   calculate compact links from non-compact links
!******************************************************************
  subroutine coef(u, theta)
    use comms
!
    complex(dp), intent(out) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    real, intent(in) :: theta(ksizex_l, ksizey_l, ksizet_l, 3)
!
!     u(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = exp(cmplx(0.0, theta))
    u(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = (1.0 + cmplx(0.0, theta))
!
!  anti-p.b.c. in timelike direction
    if(ibound.eq.-1 .and. ip_t .eq. np_t)then
       u(:, :, ksizet_l, 3) = -u(:, :, ksizet_l, 3)
    end if
!      
!!!!    call complete_halo_update_4(3, u)
    return
  end subroutine coef
!**********************************************************************
! calculate vector of gaussian random numbers with unit variance
! to refresh momenta
!   Numerical Recipes pp.203
!**********************************************************************
  subroutine gaussp(ps, reqs)
    use random
    use comms
    real, intent(out) :: ps(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 2)
#ifdef MPI
    type(MPI_Request), intent(out) :: reqs(12)
#else
    integer, intent(out), optional :: reqs
#endif
    integer ix, iy, it
    real :: theta
!     write(6,1)
!1   format(' Hi from gaussp')
    do it = 1, ksizet_l
       do iy = 1, ksizey_l
          do ix = 1, ksizex_l
             ps(ix, iy, it, 2) = sqrt(-2.0 * log(rano(yran, idum, ix, iy, it)))
          end do
       end do
    end do
    do it = 1, ksizet_l
       do iy = 1, ksizey_l
          do ix = 1, ksizex_l
             theta = tpi * rano(yran, idum, ix, iy, it)
             ps(ix, iy, it, 1) = ps(ix, iy, it, 2) * sin(theta)
             ps(ix, iy, it, 2) = ps(ix, iy, it, 2) * cos(theta)
          end do
       end do
    end do
#ifdef MPI
    call start_halo_update_4_real(2, ps, 13, reqs)
#else
    call update_halo_4_real(2, ps)
#endif
    return
  end subroutine gaussp
!**********************************************************************
! calculate vector of gaussian random numbers with unit variance
! to generate pseudofermion fields R
!   Numerical Recipes pp.203
!**********************************************************************
  subroutine gauss0(ps, reqs)
    use random
    use comms
    real, intent(out) :: ps(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 2)
#ifdef MPI
    type(MPI_Request), intent(out) :: reqs(12)
#else
    integer, intent(out), optional :: reqs
#endif
    integer :: ix, iy, it
    real :: theta
!     write(6,1)
!1   format(' Hi from gauss0')
    do it = 1, ksizet_l
       do iy = 1, ksizey_l
          do ix = 1, ksizex_l
             ps(ix, iy, it, 2) = sqrt(-log(rano(yran, idum, ix, iy, it)))
          end do
       end do
    end do
    do it = 1, ksizet_l
       do iy = 1, ksizey_l
          do ix = 1, ksizex_l
             theta = tpi * rano(yran, idum, ix, iy, it)
             ps(ix, iy, it, 1) = ps(ix, iy, it, 2) * sin(theta)
             ps(ix, iy, it, 2) = ps(ix, iy, it, 2) * cos(theta)
          end do
       end do
    end do
#ifdef MPI
    call start_halo_update_4_real(2, ps, 14, reqs)
#else
    call update_halo_4_real(2, ps)
#endif
    return
  end subroutine gauss0

!***********************************************************************
  pure subroutine dslash(Phi,R,u,am,imass)
    use dirac
!
!     calculates Phi = M*R
!
!     complex, intent(in) :: u(0:ksize+1,0:ksize+1,0:ksizet+1,3)
!     complex, intent(in) :: Phi(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
!     complex, intent(in) :: R(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
!     complex :: zkappa
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: zkappa
    real :: diag
    integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork
!     write(6,*) 'hi from dslash'
!
!     diagonal term
    diag=(3.0-am3)+1.0
    Phi=diag*R
!      
!     Wilson term (hermitian) and Dirac term (antihermitian)
    do mu=1,3
       ixup = kdelta(1, mu)
       iyup = kdelta(2, mu)
       itup = kdelta(3, mu)

       do idirac=1,4
          igork=gamin(mu,idirac)
          do it = 1,ksizet_l
             do iy = 1,ksizey_l
                do ix = 1,ksizex_l
                   Phi(:,ix,iy,it,idirac)=Phi(:,ix,iy,it,idirac) &
! Wilson term (hermitian)
                   &    -akappa*(u(ix,iy,it,mu) &
                   &              * R(:, ix+ixup, iy+iyup, it+itup, idirac) &
                   &             + conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                   &              * R(:, ix-ixup, iy-iyup, it-itup, idirac)) &
! Dirac term (antihermitian)
                   &     + gamval(mu,idirac) * &
                   &       (u(ix,iy,it,mu) &
                   &         * R(:, ix+ixup, iy+iyup, it+itup, igork) &
                   &        - conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                   &         * R(:, ix-ixup, iy-iyup, it-itup, igork))
                enddo
             enddo
          enddo
       enddo
    enddo
!
!  s-like term exploiting projection
    Phi(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
         & = Phi(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
         & - R(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
    Phi(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
         & = Phi(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
         & - R(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
!
!  Mass term (couples the two walls unless imass=5)
    if (imass.eq.1) then
       zkappa=cmplx(am,0.0)
       Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
            & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
            & + zkappa * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
       Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
            & Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) + &
            & zkappa * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
    elseif (imass.eq.3) then
       zkappa=cmplx(0.0,-am)
       Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
            & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
            & - zkappa * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
       Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
            & Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
            & + zkappa * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
    elseif (imass.eq.5) then
       zkappa=cmplx(0.0,-am)
!         do idirac=3,4
!         igork=gamin(5,idirac)
       Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
            & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
            & - zkappa * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
!        Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idirac) = &
!            & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idirac) &
!            & + 2 * zkappa * gamval(5,idirac) * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, igork)
!         enddo
!         do idirac=1,2
!         igork=gamin(5,idirac)
       Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
            & Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
            & - zkappa * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
!        Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idirac) = &
!            & Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, idirac) 
!            & + 2 * zkappa * gamval(5,idirac) * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, igork)
!         enddo
    endif
!
    return
  end subroutine dslash
!***********************************************************************
  pure subroutine dslashd(Phi,R,u,am,imass)
    use dirac
!
!     calculates Phi = Mdagger*R
!
!     complex, intent(in) ::  u(0:ksize+1,0:ksize+1,0:ksizet+1,3)
!     complex, intent(out) :: Phi(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
!     complex, intent(in) :: R(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
!     complex :: zkappa
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: zkappa
    real :: diag
    integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork
!     write(6,*) 'hi from dslashd'
!
!     diagonal term (hermitian)
    diag=(3.0-am3)+1.0
    Phi = diag * R
!
!     Wilson term (hermitian) and Dirac term (antihermitian)
    do mu=1,3
       ixup = kdelta(1, mu)
       iyup = kdelta(2, mu)
       itup = kdelta(3, mu)

       do idirac=1,4
          igork=gamin(mu,idirac)
          do it = 1,ksizet_l
             do iy = 1,ksizey_l
                do ix = 1,ksizex_l
                   Phi(:,ix,iy,it,idirac)=Phi(:,ix,iy,it,idirac) &
! Wilson term (hermitian)
                   &    - akappa * (u(ix,iy,it,mu) &
                   &              * R(:, ix+ixup, iy+iyup, it+itup, idirac) &
                   &             + conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                   &              * R(:, ix-ixup, iy-iyup, it-itup, idirac)) &
! Dirac term (antihermitian)
                   &    - gamval(mu,idirac) * &
                   &       (u(ix,iy,it,mu) &
                   &         * R(:, ix+ixup, iy+iyup, it+itup, igork) &
                   &        - conjg(u(ix-ixup, iy-iyup, it-itup, mu)) &
                   &         * R(:, ix-ixup, iy-iyup, it-itup, igork))
                enddo
             enddo
          enddo
       enddo
    enddo
!
!  s-like term exploiting projection
    Phi(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
         & = Phi(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
         & - R(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
    Phi(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
         & = Phi(2:kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
         & - R(1:kthird-1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
!
!  Mass term (couples the two walls unless imass=5) 
    if(imass.eq.1)then
       zkappa=cmplx(am,0.0)
       Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
            & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
            & + zkappa * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
       Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
            & Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
            & + zkappa * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
    elseif(imass.eq.3)then
       zkappa = cmplx(0.0,am)
       Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
            & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
            & + zkappa * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
       Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
            & Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
            & - zkappa * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
    elseif(imass.eq.5)then
       zkappa = cmplx(0.0,am)
       Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) = &
            & Phi(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2) &
            & - zkappa * R(kthird, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4)
       Phi(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) = &
            & Phi(1,1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 3:4) &
            & - zkappa * R(1, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, 1:2)
    endif
!      call complete_halo_update_5(4, Phi)
!
    return
  end subroutine dslashd

!***********************************************************************
  pure subroutine dslash2d(Phi,R,u)
    use dirac
!
!     calculates Phi = M*R
!
!     complex, intent(in) :: u(ksizex_l, ksizey_l, ksizet_l, 3)
!     complex, intent(out) :: Phi(ksizex_l,ksizey_l,ksizet_l, 4)
!     complex, intent(in) :: R(ksizex_l,ksizey_l,ksizet_l,4)
    complex(dp), intent(in) ::  u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(out) :: Phi(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer :: ix, iy, it, idirac, mu, ixup, iyup, igork
    real :: diag

!     write(6,*) 'hi from dslash2d'
!
!     diagonal term
    diag=2.0
    Phi = diag * R

!     Wilson and Dirac terms
    do mu=1,2
       ixup=kdelta(1,mu)
       iyup=kdelta(2,mu)
!
       do idirac=1,4
          igork=gamin(mu,idirac)
          do it=1,ksizet_l
             do iy=1,ksizey_l
                do ix=1,ksizex_l
                   Phi(ix,iy,it,idirac) = &
! Wilson term
                   &    Phi(ix,iy,it,idirac) &
                   &    - akappa * (u(ix,iy,it,mu) * R(ix+ixup, iy+iyup, it, idirac) &
                   &             + conjg(u(ix-ixup, iy-iyup, it, mu)) &
                   &              * R(ix-ixup, iy-iyup, it, idirac)) &
! Dirac term
                   &     + gamval(mu,idirac) * &
                   &      (u(ix,iy,it,mu)*R(ix+ixup, iy+iyup, it, igork) &
                   &       - conjg(u(ix-ixup, iy-iyup, it,mu)) &
                   &        * R(ix-ixup, iy-iyup, it, igork))
                enddo
             enddo
          enddo
       enddo
    enddo
!    call complete_halo_update_4(4, Phi)
!
    return
  end subroutine dslash2d
!

!***********************************************************************
!   A Kronecker delta function
!   Useful for calculating coordinate offsets
!***********************************************************************
  pure integer function kdelta(nu, mu)
    integer, intent(in) :: nu
    integer, intent(in) :: mu

    kdelta=merge(1,0,nu==mu)
  end function kdelta

end module dwf3d_lib
