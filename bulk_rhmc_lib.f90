module params
  ! Type definitions
  integer, parameter :: k4b=selected_int_kind(9)
  integer, parameter :: dp=kind(1.d0)

  ! Lattice parameters
  integer, parameter :: ksize=12, ksizet=12, kthird=24
  integer, parameter :: kvol=ksize*ksize*ksizet
  integer, parameter :: ndiag=25, ndiagg=12
  integer, parameter :: Nf=1
  real, parameter :: akappa = 0.5
end module params

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

module dwf3d_lib
  use params
  implicit none

  ! Random numbers
  real(dp) :: seed

! Useful constants
  real, parameter :: One = 1.0

contains

  subroutine dwf3d_main
    use random
    use remez
    use remezg
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
    common/gauge/theta(ksize, ksize, ksizet, 3)
    common/trial/ut(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3), &
         &     thetat(ksize, ksize, ksizet, 3), &
         &     pp(ksize, ksize, ksizet, 3)
    common /para/beta,am3,ibound
    common/gforce/dSdpi(ksize, ksize, ksizet, 3)
!       common /neighb/id(kvol,3),iu(kvol,3)
    common/param/ancg,ancgh,ancgf,ancgpf
    common/parampv/ancgpv,ancghpv,ancgfpv,ancgpfpv
    common/trans/tpi 
    common/dum1/ R(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4), &
         &     ps(0:ksize+1, 0:ksize+1, 0:ksizet+1, 2)
    common/vector/X1(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
!     complex :: Phi(kthird,0:ksize+1, 0:ksize+1, 0:ksizet+1, 4, Nf)
!     complex R,qq,qbqb
!     complex u,ut,X1
!     complex a,b
    complex(dp) :: Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)!
    complex(dp) :: R,qq,qbqb
    complex(dp) :: u,ut,X1
    complex(dp) :: a,b
    real(dp) :: H0,H1,S0,S1,dH,dS,hg,hp
    real :: action, paction, gaction
    real :: vel2, x, ytest, atraj
    real :: thetat, theta, dSdpi
    real :: tpi, dt, beta, am3, am, y, traj, proby
    real :: actiona, vel2a, pbp, pbpa, yav, yyav
    real :: ps, pp
    real :: ancg, ancgh, ancgf, ancgpf, ancgm, &
         &     ancgpv, ancgfpv, ancghpv, ancgpfpv, ancgma
    integer :: imass, iter, iterl, iter2, i, ia, idirac, ithird
    integer :: naccp, ipbp, itot, isweep, itercg, mu
!
    integer :: ibound
!*******************************************************************
!     input
!*******************************************************************
    integer, parameter :: istart=-1
    integer, parameter :: iread=1
    integer, parameter :: iwrite=0
    integer, parameter :: iprint=5
    integer, parameter :: iseed=1
    integer, parameter :: icheck=100
    complex(dp), parameter :: zi=(0.0,1.0)
!      integer(k4b) :: idum
    ibound=-1
    tpi=2.0*acos(-1.0)
!*******************************************************************
!     end of input
!*******************************************************************
    open(unit=7,file='output',status='unknown')
    open(unit=25,file='midout',status='unknown')
    open(unit=98,file='control',status='unknown')
    open(unit=36,file='remez2',status='unknown')
    open(unit=37,file='remez4',status='unknown')
    open(unit=38,file='remez2g',status='unknown')
    open(unit=39,file='remez4g',status='unknown')
    if(iread.eq.1) then
       call sread
    endif
    read(25,*) dt,beta,am3,am,imass,iterl,iter2
    close(25)
! set a new seed by hand...
!      if(iseed.ne.0)then
!         idum=-413976497
!      endif
    if(iseed.ne.0)then
       seed=4139764973254.0
    endif
    write(7,*) 'seed: ', seed !idum
    call rranset(seed)
    idum=-1
    y=rano(yran,idum)
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
!     write(6, 9001)ksize,ksizet,kthird,Nf,dt,traj,ndiag,ndiagg, &
!    & iter2,beta,am3,am,imass
    write(7, 9001)ksize,ksizet,kthird,Nf,dt,traj,ndiag,ndiagg, &
         & iter2,beta,am3,am,imass
9001 format(' ksize=',i3,' ksizet=',i3,/ &
         ,' kthird=',i3,/ &
         ,' Nf =',i3,/ &
         ,' time step: dt=',f6.4,' trajectory length=',f9.6,/ &
         ,' Remez ndiag: action =',i3,' guidance=',i3,/ &
         ,' # trajectories=',i6,' beta=',f9.6,/ &
         ,' am3=',f6.4,' am=',f6.4/ &
         ,' imass=',i2)
!     write(6,9004) rescgg,rescga,respbp
    write(7,9004) rescgg,rescga,respbp
9004 format(' Stopping residuals: guidance: ',e11.4,' acceptance: ', &
         &     e11.4,' estimator: ',e11.4)
!     write(6,9044) rescgm
    write(7,9044) rescgm
9044 format(' Stopping residuals: meson: ',e11.4)
    call rranget(seed)
! c     write(6,*) 'seed: ', seed
    write(7,*) 'seed: ', seed
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
                call gauss0(ps)
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
          call gaussp(ps)
          pp(:,:,:,mu) = ps(1:ksize, 1:ksize, 1:ksize, 1)
       enddo
!     write(6,*) idum
!*******************************************************************
!  call to Hamiltonian
!     
       call hamilton(Phi, H0, hg, hp, S0, rescga, isweep, 0, am, imass)
       if(isweep.eq.1) then
          action = S0 / kvol
          gaction = hg / kvol
          paction = hp / kvol
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
          ytest=rano(yran,idum)
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
       write(98,*) dH,dS
       y = exp(dH)      
       yav = yav + y 
       yyav = yyav + y*y 
       if(dH.lt.0.0)then
          x=rano(yran,idum)
          if(x.gt.y)goto 600
       endif
!
!     step accepted: set s=st
!
       theta = thetat
       naccp = naccp+1
       action = S1/kvol
       gaction = hg/kvol
       paction = hp/kvol
600    continue
       write(11,*) isweep,gaction,paction
       actiona=actiona+action 
       vel2 = sum(pp * pp) / (3 * kvol)
       vel2a = vel2a + vel2
!
!     uncomment to disable measurements
!     goto 601
666    if((isweep/iprint)*iprint.eq.isweep)then
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
          call rranget(seed)
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
!     write(6, 9022) iter2,naccp,atraj,yav,yyav,ancg,ancgpv,ancgh,ancghpv,ancgf,
!    & ancgfpv,ancgpf,ancgpfpv,pbpa,vel2a,actiona
    write(7, 9022) iter2,naccp,atraj,yav,yyav, &
         & ancg,ancgpv,ancgh,ancghpv,ancgf,ancgfpv,ancgpf,ancgpfpv, &
         & pbpa,ancgma,vel2a,actiona
9022 format(' averages for last ',i6,' trajectories',/  &
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
9024 format(1x)
!
    close(11)
!
    if(iwrite.eq.1) then
       call rranget(seed)
       call swrite
       write(7,*) 'seed: ', idum
    endif
!
    stop
  end subroutine dwf3d_main
!******************************************************************
!   calculate dSds for gauge fields at each intermediate time
!******************************************************************
  subroutine force(Phi,res1,am,imass,isweep,iter)
    use remezg
    complex(dp), intent(in) :: Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4, Nf)
    real, intent(in) :: res1, am
    integer, intent(in) :: imass, isweep, iter
    common/phi0/Phi0(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4, 25)
    common/trial/u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3), &
         &     theta(ksize, ksize, ksizet, 3), &
         &     pp(ksize, ksize, ksizet, 3)
    common/para/beta,am3,ibound
    common/param/ancg,ancgh,ancgf,ancgpf
    common/parampv/ancgpv,ancghpv,ancgfpv,ancgpfpv
    common/gforce/dSdpi(ksize, ksize, ksizet, 3)
    common/vector/X1(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
!     complex Phi(kferm,Nf),X2(kferm)
!     complex X1,u,Phi0
    complex(dp) :: X2(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp) :: X1,u,Phi0
    real :: ancg, ancgh, ancgf, ancgpf
    real :: ancgpv, ancghpv, ancgfpv, ancgpfpv
    real :: theta, pp, beta, am3, dSdpi
    integer :: ibound, ia, itercg
!
!     write(6,111)
111 format(' Hi from force')
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

       call qmrherm(X2,res1,itercg,One,1,anum4g,aden4g,ndiagg,1,isweep,iter)
       ancgpv=ancgpv+float(itercg)

       X2 = X1
!
       call qmrherm(X2,res1,itercg,am,imass,bnum2g,bden2g,ndiagg,0,isweep,iter)
       ancg=ancg+float(itercg)
!     write(111,*) itercg
       X2 = X1
!
!  evaluates -X2dagger * d/dpi[{MdaggerM(m)}^1/2] * X2
       call qmrherm(X2,res1,itercg,am,imass,anum2g,aden2g,ndiagg,2,isweep,iter)
       ancgf=ancgf+float(itercg)

!     write(113,*) itercg
!  evaluates +2Re{Phidagger * d/dpi[{MdaggerM(1)}^1/4] * X2}
       call qmrherm(X2,res1,itercg,One,1,anum4g,aden4g,ndiagg,3,isweep,iter)
       ancgfpv=ancgfpv+float(itercg)
!
    enddo
!
    if(ibound.eq.-1)then
       dSdpi(:, :, ksizet, 3) = -dSdpi(:, :, ksizet, 3)
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
    complex(dp), intent(in) :: Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4, Nf)
    real(dp), intent(out) :: h, hg, hp, s
    real, intent(in) :: res2, am
    integer, intent(in) :: isweep, iflag, imass
    common/trial/u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3), &
    &     theta(ksize, ksize, ksizet, 3), &
    &     pp(ksize, ksize, ksizet, 3)
    common/param/ancg,ancgh,ancgf,ancgpf
    common/parampv/ancgpv,ancghpv,ancgfpv,ancgpfpv
    common /para/beta,am3,ibound
    common/vector/X1(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    common/dum1/R(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4), &
    &     ps(0:ksize+1, 0:ksize+1, 0:ksizet+1, 2)
!     complex, intent(in) :: Phi(kthird, ksize, ksize, ksizet, 4, Nf)
!     complex X1,R
!     complex u
    complex(dp) :: X1,R
    complex(dp) :: u
    real(dp) :: hf
    real :: ancg, ancgh, ancgf, ancgpf, ancgma
    real :: ancgpv, ancgfpv, ancghpv, ancgpfpv
    real :: pp, ps, beta, theta, am3
    integer :: itercg, ia, ibound
!     write(6,111)
111 format(' Hi from hamilton')
!
    hf=0.0
!
    hp = 0.5 * sum(pp ** 2)

!      print *, pp, theta

    hg = 0.5 * Nf * beta * sum(theta ** 2)
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
       hf = hf + sum(real(conjg(R(:, 1:ksize, 1:ksize, 1:ksizet, :)) &
       &        * X1(:, 1:ksize, 1:ksize, 1:ksizet, :)))
!
    enddo
!
    h = hg + hp + hf
!     write(6,*) isweep,':  hg', hg,'   hp', hp,'   hf', hf,
!    &   '   h',h
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
  subroutine qmrherm(Phi, res, itercg, am, imass, anum, aden, ndiagq, iflag, isweep, iter)
    implicit none
    complex(dp), intent(in) :: Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    integer, intent(in) :: imass, ndiagq, iflag, isweep, iter
    real(dp), intent(in) :: anum(0:ndiagq), aden(ndiagq)
    real, intent(in) :: res, am
    integer, intent(out) :: itercg
!
    integer, parameter :: niterc=7500
    common/trial/u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3), &
    &             theta(kthird, ksize, ksize, ksizet, 3), &
    &             pp(kthird, ksize, ksize, ksizet, 3)
    complex(dp) :: u
    real :: theta, pp
    common/para/bbb,am3,ibound
    real :: bbb, am3
    integer :: ibound
    common/vector/x(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp) :: x
    common/gforce/dSdpi(ksize, ksize, ksizet,3)
    real :: dSdpi
    common/phi0/Phi0(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4, 25)
    complex(dp) :: Phi0
!     complex, intent(in) :: Phi(kthird, ksize, ksize, ksizet, 4)
!     complex :: x,u,Phi0
    real :: alphatild
    real(dp) :: coeff
!      
!     complex :: vtild(kthird, ksize, ksize, ksizet, 4)
!     complex :: q(kthird, ksize, ksize, ksizet, 4)
!     complex :: pm1(kthird, ksize, ksize, ksizet, 4,ndiagq)
!     complex :: qm1(kthird, ksize, ksize, ksizet, 4)
!     complex :: p(kthird, ksize, ksize, ksizet, 4,ndiagq)
!     complex :: x3(kthird, ksize, ksize, ksizet, 4)
!     complex :: R(kthird, ksize, ksize, ksizet, 4)
!     complex x1(kthird, ksize, ksize, ksizet, 4,ndiagq)
!     complex :: x2(kthird, ksize, ksize, ksizet, 4)
!     real alpha(ndiagq),beta
!     real amu(ndiagq),d(ndiagq),dm1(ndiagq),rho(ndiagq),rhom1(ndiagq)
    complex(dp) :: vtild(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp) :: q(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp) :: pm1(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4, ndiagq)
    complex(dp) :: qm1(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp) :: p(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1,4, ndiagq)
    complex(dp) :: x3(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp) :: R(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp) :: x1(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4, ndiagq)
    complex(dp) :: x2(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    real(dp) :: alpha(ndiagq), beta, beta0, phimod
    real(dp) :: amu(ndiagq), d(ndiagq), dm1(ndiagq)
    real(dp) :: rho(ndiagq), rhom1(ndiagq)
!      
    real :: resid, rhomax, arelax
    integer :: niter, idiag, ix, iy, it
!
!     write(6,111)
111 format(' Hi from qmrherm')
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

    beta = sqrt(sum(abs(R(:,1:ksize,1:ksize,1:ksizet,:)) ** 2))
    phimod=beta
!     write(6,*) '|| Phi || = ', phimod
!
    do niter=1,niterc
       itercg=itercg+1
!
!  Lanczos steps
!
       q = R / beta

       call dslash(vtild,q,u,am,imass)
       call update_halo_5(4, vtild)
       call dslashd(x3,vtild,u,am,imass)
       call update_halo_5(4, x3)
!
       alphatild = sum(real(conjg(q(:,1:ksize,1:ksize,1:ksizet,:)) & 
       &                * x3(:,1:ksize,1:ksize,1:ksizet,:)))
!
       R = x3 - alphatild * q - beta * qm1
       qm1 = q
!
       beta0=beta
       beta = sqrt(sum(abs(R(:,1:ksize,1:ksize,1:ksizet,:)) ** 2))
!
       alpha = alphatild + aden
!
       if(niter.eq.1)then
          d = alpha
          rho = beta0 / alpha
          rhom1 = rho
          do idiag = 1, ndiagq
             p(:, :, :, :, :, idiag) = q
             pm1(:, :, :, :, :, idiag) = q
             x1(:, :, :, :, :, idiag) = rho(idiag) * q
          enddo
       else
          amu = beta0 / d
          dm1 = d
          d = alpha - beta0 * amu
          rho = -amu * dm1 * rhom1 / d
          do idiag = 1, ndiagq
             p(:, :, :, :, :, idiag) = q - amu(idiag) * pm1(:, :, :, :, :, idiag)
          enddo
          pm1 = p
!     Convergence criterion (a bit ad hoc for now...)
          rhomax = maxval(abs(phimod * rho))
          rhom1 = rho
          do idiag = 1, ndiagq
             x1(:, :, :, :, :, idiag) = x1(:, :, :, :, :, idiag) &
             &           + rho(idiag) * p(:, :, :, :, :, idiag)
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
    enddo
    if (niter .gt. niterc) then
       write(7,*) 'QMRniterc!, isweep,iter,iflag,imass,anum,ndiagq = ' &
       &        ,isweep, iter, iflag, imass, anum(0), ndiagq
    end if
!     
    if(iflag.lt.2)then
!     Now evaluate solution x=(MdaggerM)^p * Phi
       do idiag=1,ndiagq
          x = x + anum(idiag) * x1(:, :, :, :, :, idiag)
       enddo
!     
!  update phi0 block if required...
       if(iflag.eq.1) then
          Phi0(:, :, :, :, :, 1:ndiagq) = X1(:, :, :, :, :, 1:ndiagq)
       endif
!     
    else
!
       do idiag=1, ndiagq
!
!  X2 = M*X1
          R = X1(:, :, :, :, :, idiag)
          call dslash(X2, R, u, am, imass)
          call update_halo_5(4, X2)
!
          if(iflag.eq.2)then
             coeff=anum(idiag)
             call derivs(R, X2, coeff, 0)
          else
             coeff=-anum(idiag)
             R = Phi0(:, :, :, :, :, idiag)
             call derivs(R, X2, coeff, 0)
!
             call dslash(X2, R, u, am, imass)
             call update_halo_5(4, X2)
!
             R = x1(: ,:, :, :, :, idiag)
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
    implicit none
!      complex, intent(in) :: R(kthird, 0:ksize+1, 0:ksize+1, 
!    &                                  0:ksizet+1, 4)
!      complex, intent(in) :: X2(kthird, 0:ksize+1, 0:ksize+1, 
!    &                                   0:ksizet+1, 4)

    complex(dp), intent(in) :: R(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp), intent(in) :: X2(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    real(dp), intent(in) :: anum
    integer, intent(in) :: iflag


    common/dirac/gamval(6,4),gamin(6,4)
    common/gforce/dSdpi(ksize,ksize,ksizet,3)

!      complex :: gamval
    complex(dp) :: gamval
    integer :: gamin
    real :: dSdpi

!      complex(dp) :: tzi
    real(dp) :: tzi_real
    integer :: ix, iy, it, ixup, iyup, itup, idirac, ithird, mu
    integer :: igork1
!
!     write(6,111)
111 format(' Hi from derivs')

!     dSdpi=dSdpi-Re(Rdagger *(d(Mdagger)dp)* X2)
!     Cf. Montvay & Muenster (7.215)
!      tzi=cmplx(0.0,2*anum)
    tzi_real = 2 * anum
!     factor of 2 picks up second term in M&M (7.215)
!
    do mu = 1,3
       ixup = kdelta(1, mu)
       iyup = kdelta(2, mu)
       itup = kdelta(3, mu)

       do idirac=1,4
!      do ithird=1,kthird
!
          do it = 1,ksizet
             do iy = 1,ksize
                do ix = 1,ksize
                   dSdpi(ix,iy,it,mu) = &
                   &     dSdpi(ix,iy,it,mu) + tzi_real * akappa * sum(dimag( &
                   &       conjg(R(:,ix,iy,it,idirac)) * &
                   &         X2(:,ix+ixup,iy+iyup,it+itup,idirac)) &
                   &     - dimag(conjg(R(:,ix+ixup,iy+iyup,it+itup,idirac)) * &
                   &         X2(:,ix,iy,it,idirac)))
                enddo
             enddo
          enddo
!
          igork1=gamin(mu,idirac)
          if(iflag.eq.0)then
             do it = 1,ksizet
                do iy = 1,ksize
                   do ix = 1,ksize
                      dSdpi(ix,iy,it,mu) = &
                      &     dSdpi(ix,iy,it,mu)+ tzi_real * sum(dimag(gamval(mu,idirac)* &
                      &(conjg(R(:,ix,iy,it,idirac)) * &
                      &        X2(:, ix+ixup,iy+iyup,it+itup,igork1) &
                      &+conjg(R(:,ix+ixup,iy+iyup,it+itup,idirac))* &
                      &             X2(:,ix,iy,it,igork1))))
                   enddo
                enddo
             enddo
          else
             do it = 1,ksizet
                do iy = 1,ksize
                   do ix = 1,ksize
                      dSdpi(ix,iy,it,mu) = &
                      &     dSdpi(ix,iy,it,mu)- tzi_real * sum(dimag(gamval(mu,idirac)* &
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
!      enddo
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
    complex(dp), intent(in) :: Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
!     complex, intent(in) :: Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    real, intent(in) :: res, am
    integer, intent(out) :: itercg
    integer, intent(in) :: imass

    integer, parameter :: niterc=kthird*kvol
    common/trial/u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3), &
    &     thetas(ksize, ksize, ksizet, 3), &
    &     pp(ksize, ksize, ksizet, 3)
    real :: thetas, pp
    common/para/bbb,am3,ibound
    real :: bbb, am3
    integer :: ibound
    common/vector/x(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
!     complex x,u
!     complex x1(kferm),x2(kferm),p(kferm),r(kferm)
    complex(dp) x,u
    complex(dp) :: x1(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp) :: x2(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp) :: p(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp) :: r(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    real :: resid
    real :: beta, betan, betad, alpha, alphan, alphad
    integer :: nx
!     write(6,111)
111 format(' Hi from congrad')
!
    resid = 4 * ksize * ksize * ksizet * kthird * res * res
    itercg = 0
    alphan = 0.0
!
    do nx=1,niterc
       itercg=itercg+1
       if(nx.gt.1) goto 51
!
!   initialise p=x, r=Phi(na)
       p = x
       r = Phi
!
       betad=1.0
       alpha=1.0
51     alphad=0.0
!
!  x1=Mp
       call dslash(x1,p,u,am,imass)
       call update_halo_5(4, x1)
!
       if(nx.ne.1)then
!
!   alpha=(r,r)/(p,(Mdagger)Mp)
          alphad = sum(abs(x1(:, 1:ksize, 1:ksize, 1:ksize, :)) ** 2)
          alpha = alphan / alphad
!     
!   x=x+alpha*p
          x = x + alpha * p
       end if
!     
!   x2=(Mdagger)x1=(Mdagger)Mp
       call dslashd(x2, x1, u, am, imass)
       call update_halo_5(4, x2)
!
!   r=r-alpha*(Mdagger)Mp
       r = r - alpha * x2
!
!   beta=(r_k+1,r_k+1)/(r_k,r_k)
       betan = sum(abs(r(:, 1:ksize, 1:ksize, 1:ksizet, :)) ** 2)
       beta = betan / betad
       betad = betan
       alphan = betan
!
       if(nx.eq.1) beta=0.0
!
!   p=r+beta*p
       p = r + beta * p
       if(betan.lt.resid) exit
    end do
!     write(6,1000)
    if (nx.gt.niterc) then
       write(7,1000)
1000   format(' # iterations of congrad exceeds niterc')
    end if
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
    real, intent(out) :: psibarpsi, aviter
    real, intent(in) :: res, am
    integer, intent(in) :: imass
    integer, parameter :: knoise = 10
    common/trial/u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3), &
    &     thetat(ksize, ksize, ksizet, 3), &
    &     pp(ksize, ksize, ksizet, 3)
    real :: thetat, pp
    common/para/beta,am3,ibound
    real :: beta, am3
    integer :: ibound
    common/dirac/gamval(6,4),gamin(6,4)
    common/vector/xi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
!      common/ranseed/idum
!      integer(k4b) :: idum
!     complex :: x(0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
!     complex :: Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
!     complex :: xi,gamval
!     complex :: psibarpsi1,psibarpsi2
!     complex :: u
    complex(dp) :: x(0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp) :: Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp) :: xi,gamval
    complex(dp) :: psibarpsi1,psibarpsi2
    complex(dp) :: u
    real(dp) :: cnum(0:1),cden(1)
    real :: ps(0:ksize+1, 0:ksize+1, 0:ksizet+1, 2)
    real :: pt(0:ksize+1, 0:ksize+1, 0:ksizet+1, 2)
    real(dp) :: pbp(knoise)
    integer :: gamin
    integer :: idsource, idsource2, idirac, inoise, jnoise, ithird
    integer :: iter, itercg
    real :: susclsing
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
       call gauss0(ps)
       psibarpsi1=(0.0,0.0)
       call gauss0(pt)
       psibarpsi2=(0.0,0.0)
!
       do idsource=1,2
!
!  source on domain wall at ithird=1
!   
          x = cmplx(0.0, 0.0)
          if(imass.ne.5)then
             x(:, :, :, idsource) = cmplx(ps(:,:,:,1), ps(:,:,:,2))
          else
             x(:, :, :, idsource) = cmplx(ps(:,:,:,1), ps(:,:,:,2))
          endif
!
          xi = cmplx(0.0, 0.0)
          if(imass.ne.5)then
             xi(1, :, :, :, :) = x
          else
!     xi = 0.5(1+gamma_4)*gamma_5*eta on DW at ithird=1
             do idirac=1,2
                xi(1, :, :, :, idirac) = -x(:, :, :, idirac+2)
             enddo
!     xi = 0.5(1+gamma_4)*eta on DW at ithird=1
          endif
!
! Phi= Mdagger*xi
!
          call dslashd(Phi, xi, u, am, imass)
          call update_halo_5(4, Phi)
!     call qmrherm(Phi,res,itercg,am,imass,cnum,cden,1,0)
          call congrad(Phi, res, itercg, am, imass)
          iter = iter + itercg
!
          if(imass.ne.5)then
!     pbp1 = x^dagger (0.5(1+gamma_4)) xi(kthird)
             psibarpsi1=psibarpsi1 &
             &           + sum(conjg(x(1:ksize, 1:ksize, 1:ksizet, idsource)) * &
             &               xi(kthird, 1:ksize, 1:ksize, 1:ksizet, idsource))
          else
!     pbp1 = x^dagger (0.5(1-gamma_4)) xi(1)
             psibarpsi1=psibarpsi1 &
             &           + sum(conjg(x(1:ksize, 1:ksize, 1:ksize, idsource+2)) &
             &                 * xi(1, 1:ksize, 1:ksize, 1:ksize, idsource+2))
          endif
!
!
!  source on domain wall at ithird=kthird
          idsource2=idsource+2
!
          x = cmplx(0.0, 0.0)
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
          call update_halo_5(4, Phi)
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
             &           +sum(conjg(x(1:ksize, 1:ksize, 1:ksizet, idsource2)) &
             &                * xi(1, 1:ksize, 1:ksize, 1:ksizet, idsource2))
          else
! pbp2= - x^dagger (0.5(1-gamma_4)) xi(kthird)
             psibarpsi2=psibarpsi2 &
             &           +sum(conjg(x(1:ksize, 1:ksize, 1:ksizet, idsource)) &
             &           * xi(kthird, 1:ksize, 1:ksize, 1:ksizet, idsource))
          endif
!
!  end trace on Dirac indices....
       enddo
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
!        write(6,*) real(psibarpsi1),aimag(psibarpsi1),
!    &       real(psibarpsi2),aimag(psibarpsi2)
       write(100,*) real(psibarpsi1), aimag(psibarpsi1), real(psibarpsi2), aimag(psibarpsi2)
!
! end loop on noise
    enddo
!
    psibarpsi=0.0
    susclsing=0.0
!
    psibarpsi = sum(pbp)
    do inoise=1,knoise
       susclsing = susclsing + &
       &        sum(pbp(inoise) * pbp(inoise+1:knoise))
    enddo
    psibarpsi = psibarpsi / knoise
    susclsing = 2 * kvol * susclsing / (knoise * (knoise-1))
    write(200,*) psibarpsi, susclsing
    aviter = float(iter) / (4*knoise)
    return
  end subroutine measure
! c******************************************************************
! c   Calculate meson correlators using point sources on domain walls
! c   -matrix inversion via conjugate gradient algorithm
! c       solves Mx=x1
! c     (Numerical Recipes section 2.10 pp.70-73)   
! c*******************************************************************
!       subroutine meson(res,itercg,aviter,am,imass)
!       use random
!       parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
!       parameter(ksize2=ksize*ksize)
!       parameter(akappa=0.5)
!       common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
!       common/para/beta,am3,ibound
!       common/dirac/gamval(6,4),gamin(6,4)
!       common /neighb/id(kvol,3),iu(kvol,3)
!       common/vector/xi(kthird,kvol,4)
!       common/ranseed/idum
!       common/v/v(97)
! c     complex x(kvol,4),x0(kvol,4),Phi(kthird,kvol,4)
! c     complex xi,gamval
! c     complex prop00(kvol,3:4,1:2),prop0L(kvol,3:4,3:4)
!       complex(dp) x(kvol,4),x0(kvol,4),Phi(kthird,kvol,4)
!       complex(dp) xi,gamval
!       complex(dp) prop00(kvol,3:4,1:2),prop0L(kvol,3:4,3:4)
! c     complex prop00n(kvol,3:4,1:2),prop0Ln(kvol,3:4,3:4)
!       real cpm(0:ksizet-1),cmm(0:ksizet-1)
! c     complex cpmn(0:ksizet-1),cmmn(0:ksizet-1)
! c     complex cferm1(0:ksizet-1), cferm2(0:ksizet-1)
! c     complex u
!       complex(dp) cferm1(0:ksizet-1), cferm2(0:ksizet-1)
!       complex(dp) u
!       real ps(kvol,2)
!       real ran
!       integer gamin
! c     write(6,*) 'hi from meson'
! c      
!       nsource=5
!       nsmear=10
!       c=0.25
!       iter=0
! c
!       do it=0,ksizet-1
!       cpm(it)=(0.0,0.0)
!       cmm(it)=(0.0,0.0)
! c     cpmn(it)=(0.0,0.0)
! c     cmmn(it)=(0.0,0.0)
!       cferm1(it)=(0.0,0.0)
!       cferm2(it)=(0.0,0.0)
!       enddo
! c
! c      susceptibility
!       chim=0.0
!       chip=0.0
! c
!       do ksource=1,nsource
! c
! c   random location for +m source
!       ixxx=int(ksize*rano(yran,idum))+1
!       iyyy=int(ksize*rano(yran,idum))+1
!       ittt=int(ksizet*rano(yran,idum))+1
!       isource=ixxx+ksize*((iyyy-1)+ksize*(ittt-1))
! c     write(6,*) ixxx,iyyy,ittt, isource
! c
! c  loop over Dirac index of source...   
!       do idsource=3,4
! c
! c  source on domain wall at ithird=1
!       do 300 ithird=1,kthird
!       if(ithird.eq.1)then
!       do idirac=1,4
!       do i=1,kvol
!       x(i,idirac)=(0.0,0.0)
!       enddo
! c  wall source
! c     ioff=(ittt-1)*ksize2
! c     do i=1,ksize2
! c     x(i+ioff,idirac)=(1.0,0.0)/ksize2
! c     enddo
! c
!       enddo
! c  point source at fixed site, spin...
!       x(isource,idsource)=(1.0,0.0)
! c
! c now smear it.....
! c
!       do ismear=1,nsmear
!       call dslash2d(x0,x,u)
!       do idirac=1,4
!       do i=1,kvol
!       x(i,idirac)=(1.0-c)*x(i,idirac)+c*x0(i,idirac)
!       enddo
!       enddo
!       enddo
! c
!       else
!       do idirac=1,4
!       do i=1,kvol
!       xi(ithird,i,idirac)=(0.0,0.0)
!       enddo
!       enddo
!       endif
! 300   continue
! c
! c   xi = x  on DW at ithird=1
! c
!       do idirac=1,4
!       do i=1,kvol
!       xi(1,i,idirac)=x(i,idirac)
!       enddo
!       enddo
! c
! c Phi= Mdagger*xi
! c
!       call dslashd(Phi,xi,u,am,imass)
! c
! c  preconditioning (no,really)
!       call dslashd(xi,Phi,u,am,imass)
! c  
! c xi= (MdaggerM)**-1 * Phi 
! c
! c     call congrad(Phi,res,itercg,am,imass) 
!       iter=iter+itercg
! c
!       do idsink=1,2
!       idsink2=idsink+2
!       do i=1,kvol
!       prop00(i,idsource,idsink)=xi(1,i,idsink)
!       prop0L(i,idsource,idsink2)=xi(kthird,i,idsink2)
!       enddo
!       enddo
! c
! c     if(imass.ne.1)then
! c  now evaluate with sign of mass reversed (not needed for hermitian mass term)
! c     am=-am
! c  source on domain wall at ithird=1
! c     do 400 ithird=1,kthird
! c     if(ithird.eq.1)then
! c     do idirac=1,4
! c     do i=1,kvol
! c     xi(1,i,idirac)=x(i,idirac)
! c     enddo
! c     enddo
! c     else
! c     do idirac=1,4
! c     do i=1,kvol
! c     xi(ithird,i,idirac)=(0.0,0.0)
! c     enddo
! c     enddo
! c     endif
! 400   continue
! c
! c Phi= Mdagger*xi
! c
! c     call dslashd(Phi,xi,u,am,imass)
! c
! c     call dslashd(xi,Phi,u,am,imass)
! c  
! c xi= (MdaggerM)**-1 * Phi 
! c
! c     call congrad(Phi,res,itercg,am,imass) 
! c     iter=iter+itercg
! c
! c     do idsink=1,2
! c     idsink2=idsink+2
! c     do i=1,kvol
! c     prop00n(i,idsource,idsink)=xi(1,i,idsink)
! c     prop0Ln(i,idsource,idsink2)=xi(kthird,i,idsink2)
! c     enddo
! c     enddo
! c
! c     am=-am
! c     endif
! c
! c  end loop on source Dirac index....
!       enddo
! c
! c  Now tie up the ends....
! c
! c  First C+-
! c
! c  now evaluate the trace (exploiting projection)
!       do id1=3,4
!       do id2=1,2
!       do it=0,ksizet-1
!       itt=mod((ittt+it-1),ksizet)+1
!       ioff=(itt-1)*ksize2
!       do i=1,ksize2
!       cpm(it)=cpm(it)
!      &        +prop00(i+ioff,id1,id2)*conjg(prop00(i+ioff,id1,id2))
!       enddo
!       enddo
!       enddo
!       enddo
! c
! c     if(imass.ne.1)then
! c     do id1=3,4
! c     do id2=1,2
! c     do it=0,ksizet-1
! c     itt=mod((ittt+it-1),ksizet)+1
! c     ioff=(itt-1)*ksize2
! c     do i=1,ksize2
! c     cpmn(it)=cpmn(it)
! c    &        +prop00(i+ioff,id1,id2)*conjg(prop00n(i+ioff,id1,id2))
! c     enddo
! c     enddo
! c     enddo
! c     enddo
! c     endif
! c
! c  next C--
! c  now evaluate the trace exploiting projection
!       do id1=3,4
!       do id2=3,4
!       do it=0,ksizet-1
!       itt=mod((ittt+it-1),ksizet)+1
!       ioff=(itt-1)*ksize2
!       do i=1,ksize2
!       cmm(it)=cmm(it)
!      &   +prop0L(i+ioff,id1,id2)*conjg(prop0L(i+ioff,id1,id2))
!       enddo
!       enddo
!       enddo
!       enddo
! c
! c     if(imass.ne.1)then
! c     do id1=3,4
! c     do id2=3,4
! c     do it=0,ksizet-1
! c     itt=mod((ittt+it-1),ksizet)+1
! c     ioff=(itt-1)*ksize2
! c     do i=1,ksize2
! c     cmmn(it)=cmmn(it)
! c    &   +prop0L(i+ioff,id1,id2)*conjg(prop0Ln(i+ioff,id1,id2))
! c     enddo
! c     enddo
! c     enddo
! c     enddo
! c     endif
! c
! c    now the fermion propagator
! c  = tr{ P_-*Psi(0,1)Psibar(x,Ls) + gamma_0*P_-*Psi(0,1)Psibar(x,1) }
!       do idd=3,4
!       do it=0,ksizet-1
!       itt=mod((ittt+it-1),ksizet)+1
! c correct for apbc
!       if(itt.ge.ittt)then
!         isign=1
!       else
!         isign=ibound
!       endif
! c
!       ioff=(itt-1)*ksize2
!       do i=1,ksize2
!       cferm1(it)=cferm1(it)
!      & +isign*akappa*prop0L(i+ioff,idd,idd)
!       cferm2(it)=cferm2(it)
!      & +isign*gamval(3,idd)*prop00(i+ioff,idd,gamin(3,idd))
!       enddo
!       enddo
!       enddo
! c
! c  finish loop over sources
!       enddo
! c
!       do it=0,ksizet-1
!       cpm(it)=cpm(it)/nsource
!       cmm(it)=cmm(it)/nsource
! c  Cf. (54) of 1507.07717
!       chim=chim+2*(cpm(it)+cmm(it))
!       enddo
! c     if(imass.ne.1)then
! c       if(imass.eq.3)then
! c         do it=0,ksizet-1
! c           cpmn(it)=cpmn(it)/nsource
! c           cmmn(it)=cmmn(it)/nsource
! c  Cf. (54),(61) of 1507.07717
! c           chip=chip-2*(cpmn(it)-cmmn(it))
! c         enddo
! c       else
! c         do it=0,ksizet-1
! c           cpmn(it)=cpmn(it)/nsource
! c           cmmn(it)=cmmn(it)/nsource
! c  Cf. (64),(65) of 1507.07717
! c           chip=chip-2*(cpm(it)-cmm(it))
! c         enddo
! c       endif
! c     endif
! c
!       do it=0,ksizet-1
!       write(302,*) it, cpm(it), cmm(it)
!       write(500,*) it, real(cferm1(it)), aimag(cferm1(it))
!       write(501,*) it, real(cferm2(it)), aimag(cferm2(it))
!       enddo
! c     write(6,*) chim
!       write(400,*) chim
! c     if(imass.ne.1)then
! c     do it=0,ksizet-1
! c     write(402,*) it, real(cpmn(it)), real(cmmn(it))
! c     write(403,*) it, aimag(cpmn(it)), aimag(cmmn(it))
! c     enddo
! c     write(401,*) chip
! c     endif
! c
! c     if(imass.eq.1)then
!       aviter=float(iter)/(2*nsource)
! c     else
! c     aviter=float(iter)/(4*nsource)
! c     endif
! c
!       return
!       end
!*******************************************************************
!
  subroutine sread
    use random
    implicit none
    common/gauge/ theta(ksize, ksize, ksizet, 3)
!     common/ranseed/ idum
    real :: theta
!      integer(k4b) :: idum
    open(unit=10, file='con', status='unknown', form='unformatted')
    read (10) theta, seed
!      idum = -idum
    close(10)
    return
  end subroutine sread
!
  subroutine swrite
    use random
    implicit none
    common/gauge/ theta(ksize, ksize, ksizet, 3)
!      common/ranseed/ idum
    real :: theta
!      integer(k4b) :: idum
    open(unit=31, file='con', status='unknown', form='unformatted')
    write (31) theta, seed
    close(31)
    return
  end subroutine swrite
!
  subroutine init(nc)
    use random
!*******************************************************************
!     sets initial values
!     nc=0 cold start
!     nc=1 hot start
!     nc<0 no initialization
!*******************************************************************
    implicit none
    integer, intent(in) :: nc
    common/gauge/theta(ksize, ksize, ksizet, 3), seed
    common/dirac/gamval(6,4),gamin(6,4)
!      common/ranseed/idum
!     complex gamval,one,zi
    complex(dp) :: gamval,one,zi
    real :: theta
    real(dp) :: seed
    integer :: gamin
!      integer(k4b) :: idum
    integer :: ix, iy, it, mu
    real :: g
!
!
    one=(1.0,0.0)
    zi=(0.0,1.0)
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
       do it = 1, ksizet
          do iy = 1, ksize
             do ix = 1, ksize
                theta(ix, iy, it, mu) = 2.0 * g * rranf() - 1.0
!              theta(ix, iy, it, mu) = 2.0 * g * rano(yran,idum) - 1.0
             enddo
          enddo
       enddo
    enddo
    return
  end subroutine init
!******************************************************************
!   calculate compact links from non-compact links
!******************************************************************
  pure subroutine coef(u,theta)
    implicit none
    common/para/beta,am3,ibound
!
    complex(dp), intent(inout) :: u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3)
    real, intent(in) :: theta(ksize, ksize, ksizet, 3)
    real :: beta, am3
    integer :: ibound
    integer :: ix, iy, it, mu
!
!     u(1:ksize, 1:ksize, 1:ksizet, :) = exp(cmplx(0.0, theta))
    u(1:ksize, 1:ksize, 1:ksizet, :) = (1.0 + cmplx(0.0, theta))
!
!  anti-p.b.c. in timelike direction
    if(ibound.eq.-1)then
       u(:, :, ksizet, 3) = -u(:, :, ksizet, 3)
    end if
!      
    call update_halo_4(3, u)
    return
  end subroutine coef
!**********************************************************************
! calculate vector of gaussian random numbers with unit variance
! to refresh momenta
!   Numerical Recipes pp.203
!**********************************************************************
  subroutine gaussp(ps)
    use random
    implicit none
    common/trans/tpi 
!      common/ranseed/idum
    real, intent(out) :: ps(0:ksize+1, 0:ksize+1, 0:ksizet+1, 2)
    real :: tpi
    integer ix, iy, it
!      integer(k4b) :: idum
    real :: theta
!     write(6,1)
1   format(' Hi from gaussp')
    do it = 1, ksizet
       do iy = 1, ksize
          do ix = 1, ksize
             ps(ix, iy, it, 2) = sqrt(-2.0 * log(rano(yran,idum)))
          end do
       end do
    end do
    do it = 1, ksizet
       do iy = 1, ksize
          do ix = 1, ksize
             theta = tpi * rano(yran,idum)
             ps(ix, iy, it, 1) = ps(ix, iy, it, 2) * sin(theta)
             ps(ix, iy, it, 2) = ps(ix, iy, it, 2) * cos(theta)
          end do
       end do
    end do
    call update_halo_4_real(2, ps)

    return
  end subroutine gaussp
!**********************************************************************
! calculate vector of gaussian random numbers with unit variance
! to generate pseudofermion fields R
!   Numerical Recipes pp.203
!**********************************************************************
  subroutine gauss0(ps)
    use random
    implicit none
    common/trans/tpi 
!      common/ranseed/idum
    real :: tpi
!      integer(k4b) :: idum
    real, intent(out) :: ps(0:ksize+1, 0:ksize+1, 0:ksizet+1, 2)
    integer :: ix, iy, it
    real :: theta
!     write(6,1)
1   format(' Hi from gauss0')
    do it = 1, ksizet
       do iy = 1, ksize
          do ix = 1, ksize
             ps(ix, iy, it, 2) = sqrt(-log(rano(yran,idum)))
          end do
       end do
    end do
    do it = 1, ksizet
       do iy = 1, ksize
          do ix = 1, ksize
             theta = tpi * rano(yran,idum)
             ps(ix, iy, it, 1) = ps(ix, iy, it, 2) * sin(theta)
             ps(ix, iy, it, 2) = ps(ix, iy, it, 2) * cos(theta)
          end do
       end do
    end do
    call update_halo_4_real(2, ps)
    return
  end subroutine gauss0

!*****************************************
!  Random number generator Numerical recipes B7
!
  real function ran(idum)
    implicit none

    integer, parameter :: k4b=selected_int_kind(9)
    integer(k4b), intent(inout) :: idum
    integer(k4b), parameter :: IA=16807, IM=2147483647, IQ=127773, IR=2836
    real, save :: am
    integer(k4b), save :: ix=-1, iy=-1, k
    if (idum <= 0 .or. iy < 0) then
       am = nearest(1.0,-1.0) / IM
       iy = ior(ieor(888889999, abs(idum)), 1)
       ix = ieor(777755555, abs(idum))
       idum = abs(idum) + 1
    end if
    ix = ieor(ix, ishft(ix, 13))
    ix = ieor(ix, ishft(ix, -17))
    ix = ieor(ix, ishft(ix, 5))
    k = iy / IQ
    iy = IA * (iy - k * IQ) - IR * k
    if (iy < 0) iy = iy + IM
    ran = am * ior(iand(IM, ieor(ix, iy)), 1)
  end function ran
!***********************************************************************
  pure subroutine dslash(Phi,R,u,am,imass)
!
!     calculates Phi = M*R
!
    implicit none
    common/para/beta,am3,ibound
    common/dirac/gamval(6,4),gamin(6,4)
!     complex, intent(in) :: u(0:ksize+1,0:ksize+1,0:ksizet+1,3)
!     complex, intent(in) :: Phi(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
!     complex, intent(in) :: R(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
!     complex gamval
!     complex zkappa
    complex(dp), intent(in) :: u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: gamval
    complex(dp) :: zkappa
    integer :: gamin
    real :: beta, am3, diag
    integer :: ibound
    integer :: ixup, iyup, itup, ix, iy, it, ithird, idirac, mu, igork
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
          do it = 1,ksizet
             do iy = 1,ksize
                do ix = 1,ksize
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
    Phi(1:kthird-1,:,:,:,3:4) = Phi(1:kthird-1,:,:,:,3:4) - R(2:kthird,:,:,:,3:4)
    Phi(2:kthird,:,:,:,1:2) = Phi(2:kthird,:,:,:,1:2) - R(1:kthird-1,:,:,:,1:2)
!
!  Mass term (couples the two walls unless imass=5)
    if (imass.eq.1) then
       zkappa=cmplx(am,0.0)
       Phi(kthird, :, :, :, 3:4) = Phi(kthird, :, :, :, 3:4) + zkappa * R(1, :, :, :, 3:4)
       Phi(1, :, :, :, 1:2) = Phi(1, :, :, :, 1:2) + zkappa * R(kthird, :, :, :, 1:2)
    elseif (imass.eq.3) then
       zkappa=cmplx(0.0,-am)
       Phi(kthird,:, :, :, 3:4) = Phi(kthird, :, :, :, 3:4) - zkappa * R(1, :, :, :, 3:4)
       Phi(1, :, :, :, 1:2) = Phi(1, :, :, :, 1:2) + zkappa * R(kthird, :, :, :, 1:2)
    elseif (imass.eq.5) then
       zkappa=cmplx(0.0,-am)
!         do idirac=3,4
!         igork=gamin(5,idirac)
       Phi(kthird, :, :, :, 3:4) = Phi(kthird, :, :, :, 3:4) &
            &     - zkappa * R(kthird, :, :, :, 1:2)
!        Phi(kthird,:,:,:,idirac)=Phi(kthird,:,:,:,idirac) &
!           &     + 2 * zkappa*gamval(5,idirac)*R(kthird,:,:,:,igork)
!         enddo
!         do idirac=1,2
!         igork=gamin(5,idirac)
       Phi(1, :, :, :, 1:2) = Phi(1, :, :, :, 1:2) - zkappa * R(1, :, :, :, 3:4)
!        Phi(1,:,:,:,idirac)=Phi(1,:,:,:,idirac)+2*zkappa*gamval(5,idirac)*R(1,:,:,:,igork)
!         enddo
    endif
!
    return
  end subroutine dslash
!***********************************************************************
  pure subroutine dslashd(Phi,R,u,am,imass)
!
!     calculates Phi = Mdagger*R
!
    implicit none
    common/para/beta,am3,ibound
    common/dirac/gamval(6,4),gamin(6,4)
!     complex, intent(in) ::  u(0:ksize+1,0:ksize+1,0:ksizet+1,3)
!     complex, intent(out) :: Phi(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
!     complex, intent(in) R(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
!     complex gamval
!     complex zkappa
    complex(dp), intent(in) :: u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3)
    complex(dp), intent(out) :: Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp), intent(in) :: R(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: gamval
    complex(dp) :: zkappa
    integer :: gamin
    real :: beta, am3, diag
    integer :: ibound
    integer :: ixup, iyup, itup, ix, iy, it, ithird, idirac, mu, igork
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
          do it = 1,ksizet
             do iy = 1,ksize
                do ix = 1,ksize
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
    Phi(1:kthird-1,:,:,:,1:2) = Phi(1:kthird-1,:,:,:,1:2) - R(2:kthird,:,:,:,1:2)
    Phi(2:kthird,:,:,:,3:4) = Phi(2:kthird,:,:,:,3:4) - R(1:kthird-1,:,:,:,3:4)
!
!  Mass term (couples the two walls unless imass=5) 
    if(imass.eq.1)then
       zkappa=cmplx(am,0.0)
       Phi(kthird,:,:,:,1:2) = Phi(kthird,:,:,:,1:2) + zkappa * R(1,:,:,:,1:2)
       Phi(1,:,:,:,3:4) = Phi(1,:,:,:,3:4) + zkappa * R(kthird,:,:,:,3:4)
    elseif(imass.eq.3)then
       zkappa = cmplx(0.0,am)
       Phi(kthird,:,:,:,1:2) = Phi(kthird,:,:,:,1:2) + zkappa * R(1,:,:,:,1:2)
       Phi(1,:,:,:,3:4) = Phi(1,:,:,:,3:4) - zkappa * R(kthird,:,:,:,3:4)
    elseif(imass.eq.5)then
       zkappa = cmplx(0.0,am)
       Phi(kthird,:,:,:,1:2) = Phi(kthird,:,:,:,1:2) - zkappa * R(kthird,:,:,:,3:4)
       Phi(1,:,:,:,3:4) = Phi(1,:,:,:,3:4) - zkappa * R(1,:,:,:,1:2)
    endif
!      call update_halo_5(4, Phi)
!
    return
  end subroutine dslashd

!***********************************************************************
  pure subroutine dslash2d(Phi,R,u)
!
!     calculates Phi = M*R
!
    implicit none
!      integer :: kdelta
    common/dirac/gamval(6,4),gamin(6,4)
!     complex u(ksize,ksize,ksizet,3)
!     complex Phi(ksize,ksize,ksizet,4),R(ksize,ksize,ksizet,4)
!     complex gamval
    complex(dp), intent(in) ::  u(0:ksize+1,0:ksize+1,0:ksizet+1,3)
    complex(dp), intent(out) :: Phi(0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
    complex(dp), intent(in) :: R(0:ksize+1,0:ksize+1,0:ksizet+1,4)
    complex(dp) :: gamval
    integer :: gamin
    integer :: ix, iy, it, idirac, mu, ixup, iyup, itup, igork
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
          do it=1,ksizet
             do iy=1,ksize
                do ix=1,ksize
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
    call update_halo_4(4, Phi)
!
    return
  end subroutine dslash2d
!
!***********************************************************************
!   Update boundary terms
!***********************************************************************
  pure subroutine update_halo_4(size4, Array)
!     
    implicit none
!
    integer, intent(in) :: size4
    complex(dp), intent(inout) :: Array(0:ksize+1, 0:ksize+1, 0:ksizet+1, size4)
!
    Array(0,:,:,:) = Array(ksize,:,:,:)
    Array(ksize+1,:,:,:) = Array(1,:,:,:)
    Array(:,0,:,:) = Array(:,ksize,:,:)
    Array(:,ksize+1,:,:) = Array(:,1,:,:)
    Array(:,:,0,:) = Array(:,:,ksize,:)
    Array(:,:,ksize+1,:) = Array(:,:,1,:)
!      
    return
!      
  end subroutine update_halo_4
!***********************************************************************
  pure subroutine update_halo_4_real(size4, Array)
!     
    implicit none
!
    integer, intent(in) :: size4
    real, intent(inout) :: Array(0:ksize+1, 0:ksize+1, 0:ksizet+1, size4)
!
    Array(0,:,:,:) = Array(ksize,:,:,:)
    Array(ksize+1,:,:,:) = Array(1,:,:,:)
    Array(:,0,:,:) = Array(:,ksize,:,:)
    Array(:,ksize+1,:,:) = Array(:,1,:,:)
    Array(:,:,0,:) = Array(:,:,ksize,:)
    Array(:,:,ksize+1,:) = Array(:,:,1,:)
!      
    return
!      
  end subroutine update_halo_4_real
!***********************************************************************
  pure subroutine update_halo_5(size5, Array)
!     
    implicit none
!
    integer, intent(in) :: size5
    complex(dp), intent(inout) :: Array(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, size5)
!
    Array(:,0,:,:,:) = Array(:,ksize,:,:,:)
    Array(:,ksize+1,:,:,:) = Array(:,1,:,:,:)
    Array(:,:,0,:,:) = Array(:,:,ksize,:,:)
    Array(:,:,ksize+1,:,:) = Array(:,:,1,:,:)
    Array(:,:,:,0,:) = Array(:,:,:,ksize,:)
    Array(:,:,:,ksize+1,:) = Array(:,:,:,1,:)
!      
    return
!      
  end subroutine update_halo_5
!***********************************************************************
  pure subroutine update_halo_6(size5, size6, Array)
!     
    implicit none
!
    integer, intent(in) :: size5, size6
    complex(dp), intent(inout) :: Array(kthird, 0:ksize+1, 0:ksize+1, &
         &                              0:ksizet+1, size5, size6)
!
    Array(:,0,:,:,:,:) = Array(:,ksize,:,:,:,:)
    Array(:,ksize+1,:,:,:,:) = Array(:,1,:,:,:,:)
    Array(:,:,0,:,:,:) = Array(:,:,ksize,:,:,:)
    Array(:,:,ksize+1,:,:,:) = Array(:,:,1,:,:,:)
    Array(:,:,:,0,:,:) = Array(:,:,:,ksize,:,:)
    Array(:,:,:,ksize+1,:,:) = Array(:,:,:,1,:,:)
!      
    return
!      
  end subroutine update_halo_6

!***********************************************************************
!   A Kronecker delta function
!   Useful for calculating coordinate offsets
!***********************************************************************
  pure integer function kdelta(nu, mu)
    implicit none
    integer, intent(in) :: nu
    integer, intent(in) :: mu

    kdelta=merge(1,0,nu==mu)
  end function kdelta

end module dwf3d_lib


