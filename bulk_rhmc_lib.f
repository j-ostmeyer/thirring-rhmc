      module random
! Random numbers
      real :: yran
      integer :: idum
      real :: v(97)

      DOUBLE PRECISION, PRIVATE :: DS(2),    DM(2)
      DOUBLE PRECISION, PRIVATE :: DX24,     DX48
      DATA      DS     /  1665 1885.D0, 286 8876.D0  /
      DATA      DM     /  1518 4245.D0, 265 1554.D0  /
      DATA      DX24   /  1677 7216.D0  /
      DATA      DX48   /  281 4749 7671 0656.D0  /

      contains

c*****************************************
c  Random number generator Numerical recipes 7.1
c
      real function rano(y,i)
      real, intent(out) :: y
      integer, intent(inout) :: i
      integer :: j
      real :: dum
c     
      if(i.lt.0)then
         i=1
         do j=1,97
            dum=rranf()
         enddo
         do j=1,97
            v(j)=rranf()
         enddo
         y=rranf()
      endif
c     
      j=1+int(97.0*y)
      if(j.gt.97) j=97
      if(j.lt.1) j=1
c     write(6,*) j,y
c     write(6,*) 'problems with rano'
c     stop
c     endif
      y=v(j)
      rano=y
      v(j)=rranf()
      return
      end function
C========================================================================
C     
      SUBROUTINE RRANGET(LSEED)
      DOUBLE PRECISION, INTENT(OUT) :: LSEED
      LSEED  =  G900GT()
      RETURN
      END SUBROUTINE

      SUBROUTINE RRANSET(LSEED)
      DOUBLE PRECISION, INTENT(IN) :: LSEED
      DOUBLE PRECISION DUMMY
      DUMMY  =  G900ST(LSEED)
      RETURN
      END SUBROUTINE

      

      REAL FUNCTION RRANF()
      RRANF = SNGL(DRANF())
      RETURN
      END FUNCTION

      DOUBLE PRECISION FUNCTION DRANF()
      DOUBLE PRECISION    DL,       DC,       DU,       DR
      DL  =  DS(1) * DM(1)
      DC  =  DINT(DL/DX24)
      DL  =  DL - DC*DX24
      DU  =  DS(1)*DM(2) + DS(2)*DM(1) + DC
      DS(2)  =  DU - DINT(DU/DX24)*DX24
      DS(1)  =  DL
      DR     =  (DS(2)*DX24 + DS(1)) / DX48
      DRANF  =  DR
      RETURN
      END FUNCTION

      DOUBLE PRECISION FUNCTION G900GT()
      G900GT  =  DS(2)*DX24 + DS(1)
      RETURN
      END FUNCTION

      DOUBLE PRECISION FUNCTION G900ST(DSEED)
      DOUBLE PRECISION, INTENT(IN) :: DSEED
      DS(2)  =  DINT(DSEED/DX24)
      DS(1)  =  DSEED - DS(2)*DX24
      G900ST =  DS(1)
      RETURN
      END FUNCTION
c***********************************************************************
      
      end module random     


      module dwf3d_lib
      implicit none

! Type definitions
      integer, parameter :: k4b=selected_int_kind(9)
      integer, parameter :: dp=kind(1.d0)

! Lattice parameters
      integer, parameter :: ksize=12, ksizet=12, kthird=24
      integer, parameter :: kvol=ksize*ksize*ksizet
      integer, parameter :: Nf=1
      real, parameter :: akappa = 0.5

! Random numbers
      real(dp) :: seed

! Useful constants
      real, parameter :: One = 1.0

      contains

      subroutine dwf3d_main
      use random
c*******************************************************************
c    Rational Hybrid Monte Carlo algorithm for bulk Thirring Model with Domain Wall
c         fermions
c
c    requires operation of QMR on complex vectors to determine
c    (Mdagger M)**-1  Phi 
c
c    requires input partial fraction coefficients from Remez algorithm
c    of Clarke & Kennedy
c
c    the "third" direction is actually indexed 4 in the code - sorry!
c      (should have used C I know)
c
c           { 1 - hermition mass term psibar psi
c    imass= { 3 - antih.    mass term i psibar gamma_3 psi
c           { 5 - antih.    mass term i psibar gamma_5 psi
c     linear combinations of above require code modification
c
c    code exploits fact that gamma_3 is diagonal in Dirac basis used - speeds up 
c    evolution slightly, and measurement by factor of two. 
c
c    Fermion expectation values are measured using a noisy estimator.
c    on the Wilson matrix, which has dimension 4*kvol*kthird
c    inversions done using congrad, and matrix multiplies with dslash, dslashd
c
c    Pauli-Villars determinant defined using a hermitian mass m_h=One
c
c    trajectory length is random with mean dt*iterl
c    The code runs for a fixed number iter2 of trajectories.
c
c    Phi: pseudofermion field 
c    am: bare fermion mass 
c    actiona: running average of total action
c
c                                               SJH February 2017
c*******************************************************************
      real, parameter :: respbp=0.000001, rescgg=0.000001
      real, parameter :: rescga=0.000000001
      real, parameter :: rescgm=0.000000001
      integer, parameter :: itermax=1000
      integer, parameter :: ndiag=25, ndiagg=12
      common/gauge/theta(ksize, ksize, ksizet, 3)
      common/trial/ut(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3),
     &     thetat(ksize, ksize, ksizet, 3),
     &     pp(ksize, ksize, ksizet, 3)
      common /para/beta,am3,ibound
      common/remez2/anum2(0:ndiag),aden2(ndiag),
     &              bnum2(0:ndiag),bden2(ndiag)
      common/remez4/anum4(0:ndiag),aden4(ndiag),
     &              bnum4(0:ndiag),bden4(ndiag)
      common/remez2g/anum2g(0:ndiagg),aden2g(ndiagg),
     &              bnum2g(0:ndiagg),bden2g(ndiagg)
      common/remez4g/anum4g(0:ndiagg),aden4g(ndiagg),
     &              bnum4g(0:ndiagg),bden4g(ndiagg)
      common/gforce/dSdpi(ksize, ksize, ksizet, 3)
C       common /neighb/id(kvol,3),iu(kvol,3)
      common/param/ancg,ancgh,ancgf,ancgpf
      common/parampv/ancgpv,ancghpv,ancgfpv,ancgpfpv
      common/trans/tpi 
      common/dum1/ R(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4),
     &     ps(0:ksize+1, 0:ksize+1, 0:ksizet+1, 2)
      common/vector/X1(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
!      common/ranseed/idum
C       common/v/v(97)
c     complex :: Phi(kthird,0:ksize+1, 0:ksize+1, 0:ksizet+1, 4, Nf)!,
!     &     X0(kthird,kvol,4)
c     complex R,qq,qbqb
c     complex u,ut,X1
c     complex a,b
      complex(dp) :: Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)!,
!     &     X0(kthird,kvol,4)
      complex(dp) :: R,qq,qbqb
      complex(dp) :: u,ut,X1
      complex(dp) :: a,b
      real(dp) :: H0,H1,S0,S1,dH,dS,hg,hp
      real(dp) :: anum2,aden2,bnum2,bden2
      real(dp) :: anum4,aden4,bnum4,bden4
      real(dp) :: anum2g,aden2g,bnum2g,bden2g
      real(dp) :: anum4g,aden4g,bnum4g,bden4g
      real :: action, paction, gaction
      real :: vel2, x, ytest, atraj
      real :: thetat, theta, dSdpi
      real :: tpi, dt, beta, am3, am, y, traj, proby
      real :: actiona, vel2a, pbp, pbpa, yav, yyav
      real :: ps, pp
      real :: ancg, ancgh, ancgf, ancgpf, ancgm,
     &     ancgpv, ancgfpv, ancghpv, ancgpfpv, ancgma
      integer :: imass, iter, iterl, iter2, i, ia, idirac, ithird
      integer :: naccp, ipbp, itot, isweep, itercg, mu
c
      integer :: ibound
c*******************************************************************
c     input
c*******************************************************************
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
c*******************************************************************
c     end of input
c*******************************************************************
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
c set a new seed by hand...
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
      print *,y
c     write(6,*) 'ran: ', y,idum
c*******************************************************************
c     initialization
c     istart.lt.0 : start from tape
c     istart=0    : ordered start
c     istart=1    : random start
c*******************************************************************
      call init(istart)
c  read in Remez coefficients
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
c*******************************************************************
c     print heading
c*******************************************************************
      traj=iterl*dt
      proby=1.0/float(iterl)
c     write(6, 9001)ksize,ksizet,kthird,Nf,dt,traj,ndiag,ndiagg,
c    & iter2,beta,am3,am,imass
      write(7, 9001)ksize,ksizet,kthird,Nf,dt,traj,ndiag,ndiagg,
     & iter2,beta,am3,am,imass
9001  format(' ksize=',i3,' ksizet=',i3,/
     1 ,' kthird=',i3,/
     1 ,' Nf =',i3,/
     1 ,' time step: dt=',f6.4,' trajectory length=',f9.6,/
     1 ,' Remez ndiag: action =',i3,' guidance=',i3,/
     1 ,' # trajectories=',i6,' beta=',f9.6,/
     1 ,' am3=',f6.4,' am=',f6.4/
     1 ,' imass=',i2)
c     write(6,9004) rescgg,rescga,respbp
      write(7,9004) rescgg,rescga,respbp
9004  format(' Stopping residuals: guidance: ',e11.4,' acceptance: ',
     &     e11.4,' estimator: ',e11.4)
c     write(6,9044) rescgm
      write(7,9044) rescgm
9044  format(' Stopping residuals: meson: ',e11.4)
       call rranget(seed)
! c     write(6,*) 'seed: ', seed
       write(7,*) 'seed: ', seed
c*******************************************************************
c       initialize for averages
c*******************************************************************
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
c*******************************************************************
c     start of classical evolution
c*******************************************************************
      do 601 isweep=1,iter2
c uncomment line below to go straight to measurement
c     goto 666
c*******************************************************************
c     initialise trial fields
c*******************************************************************
C       do 2007 mu=1,3
C       do 2007 i=1,kvol
       thetat = theta
C 2007  continue
c
      call coef(ut,thetat)
c*******************************************************************
c  Pseudofermion fields: Phi = {MdaggerM(1)}^-1/4 * {MdaggerM(m)}^1/4 * R, where
c   R is gaussian
c*******************************************************************
      do ia=1,Nf
c
      do idirac=1,4
         do ithird=1,kthird
            call gauss0(ps)
            R(ithird,:,:,:,idirac) = cmplx(ps(:,:,:,1), ps(:,:,:,2))
         enddo
      enddo
c
c  For now Phi = {MdaggerM}^0.25 * R
c
      call qmrherm(R,rescga,itercg,am,imass,anum4,aden4,ndiag,
     &             0,isweep,0)
      ancgpf=ancgpf+float(itercg)
c
      R = X1
c
      call qmrherm(R,rescga,itercg,One,1,bnum4,bden4,ndiag,
     &             0,isweep,0)
      ancgpfpv=ancgpfpv+float(itercg)
c
      Phi = X1
c
      enddo
c*******************************************************************
c     heatbath for p 
c*******************************************************************
c  for some occult reason this write statement is needed to ensure compatibility with earlier versions
c     write(6,*) idum
c     write(98,*) idum
      do mu=1,3
         call gaussp(ps)
         pp(:,:,:,mu) = ps(1:ksize, 1:ksize, 1:ksize, 1)
      enddo
c     write(6,*) idum
c*******************************************************************
c  call to Hamiltonian
c     
      call hamilton(Phi,
     &      H0,hg,hp,S0,rescga,isweep,0,am,imass)
      if(isweep.eq.1) then
         action = S0 / kvol
         gaction = hg / kvol
         paction = hp / kvol
      endif 
c     goto 501
c*******************************************************************
c      half-step forward for p
c*******************************************************************
      call force(Phi,rescgg,am,imass,isweep,0)
      pp = pp - 0.5 * dt * dSdpi
c*******************************************************************
c     start of main loop for classical time evolution
c*******************************************************************
      do 500 iter=1,itermax
c
c  step (i) st(t+dt)=st(t)+p(t+dt/2)*dt;
c
         thetat = thetat + dt * pp
c
c  step (ii)  p(t+3dt/2)=p(t+dt/2)-dSds(t+dt)*dt (1/2 step on last iteration)
c
         call coef(ut,thetat)
         call force(Phi,rescgg,am,imass,isweep,iter)
c
c test for end of random trajectory
c 
         ytest=rano(yran,idum)
         if(ytest.lt.proby)then
            pp = pp - 0.5 * dt * dSdpi
            itot = itot + iter 
            goto 501
         else
            pp = pp - dt * dSdpi
         endif
c 
500   continue
c**********************************************************************
c  Monte Carlo step: accept new fields with probability=
c              min(1,exp(H0-H1))
c**********************************************************************
501   continue 
      call hamilton(Phi,
     &       H1,hg,hp,S1,rescga,isweep,-1,am,imass)
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
c
c     step accepted: set s=st
c
      theta = thetat
      naccp = naccp+1
      action = S1/kvol
      gaction = hg/kvol
      paction = hp/kvol
600   continue
      write(11,*) isweep,gaction,paction
      actiona=actiona+action 
      vel2 = sum(pp * pp) / (3 * kvol)
      vel2a = vel2a + vel2
c
c     uncomment to disable measurements
c     goto 601
666   if((isweep/iprint)*iprint.eq.isweep)then
         thetat = theta
         call coef(ut,thetat)
         call measure(pbp,respbp,ancgm,am,imass)
c        call meson(rescgm,itercg,ancgm,am,imass)
         pbpa=pbpa+pbp
         ancgma=ancgma+ancgm
         ipbp=ipbp+1
c        write(11,*) pbp
c        write(6,*) isweep,':  ',pbp,ancgm
      endif
c
      if((isweep/icheck)*icheck.eq.isweep)then
      call rranget(seed)
      if(iwrite.eq.1) then
         call swrite
      endif
      flush(100)
      flush(200)
c     flush(302)
c     flush(400)
c     flush(500)
c     flush(501)
      if(imass.ne.1)then
c        flush(401)
c        flush(402)
c        flush(403)
      endif
c     write(7,9023) seed
      endif
c
601   continue
c*******************************************************************
c     end of main loop
c*******************************************************************
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
c*******************************************************************
c     print global averages
c*******************************************************************
c     write(6, 9022) iter2,naccp,atraj,yav,yyav,ancg,ancgpv,ancgh,ancghpv,ancgf,
c    & ancgfpv,ancgpf,ancgpfpv,pbpa,vel2a,actiona
      write(7, 9022) iter2,naccp,atraj,yav,yyav,
     & ancg,ancgpv,ancgh,ancghpv,ancgf,ancgfpv,ancgpf,ancgpfpv,
     & pbpa,ancgma,vel2a,actiona
9022  format(' averages for last ',i6,' trajectories',/ 
     & 1x,' # of acceptances: ',i6,' average trajectory length= ',f8.3/
     & 1x,' <exp-dH>=',e11.4,' +/-',e10.3/
     2 1x,' av. # QMR itr.'/
     & 1x,'     guidance: DWF  ',f9.3,'; PV  ',f9.3/
     & 1x,'   acceptance: DWF  ',f9.3,'; PV  ',f9.3/
     & 1x,'        force: DWF  ',f9.3,'; PV  ',f9.3/
     & 1x,'pseudofermion: DWF  ',f9.3,'; PV  ',f9.3/
     1 1x,' psibarpsi=',e11.3/
     2 1x,' av. # QMR itr.',f9.3//
     & 1x,' mean square velocity=',e10.3,'; action per site=',e10.3//)
      write(7, 9024)
      write(7, 9024)
9024  format(1x)
c
      close(11)
c
      if(iwrite.eq.1) then
      call rranget(seed)
      call swrite
      write(7,*) 'seed: ', idum
      endif
c
      stop
      end subroutine
c******************************************************************
c   calculate dSds for gauge fields at each intermediate time
c******************************************************************
      subroutine force(Phi,res1,am,imass,isweep,iter)
      complex(dp), intent(in) :: 
     &     Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4, Nf)
      real, intent(in) :: res1, am
      integer, intent(in) :: imass, isweep, iter
      integer, parameter :: ndiagg=12, ndiag=ndiagg
      common/remez2g/anum2(0:ndiag),aden2(ndiag),
     &              bnum2(0:ndiag),bden2(ndiag)
      common/remez4g/anum4(0:ndiag),aden4(ndiag),
     &              bnum4(0:ndiag),bden4(ndiag)
      common/phi0/Phi0(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4, 25)
      common/trial/u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3),
     &     theta(ksize, ksize, ksizet, 3),
     &     pp(ksize, ksize, ksizet, 3)
      common/para/beta,am3,ibound
      common/param/ancg,ancgh,ancgf,ancgpf
      common/parampv/ancgpv,ancghpv,ancgfpv,ancgpfpv
      common/gforce/dSdpi(ksize, ksize, ksizet, 3)
      common/vector/X1(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
c     complex Phi(kferm,Nf),X2(kferm)
c     complex X1,u,Phi0
      complex(dp) :: X2(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex(dp) :: X1,u,Phi0
      real(dp) :: anum2,aden2,bnum2,bden2
      real(dp) :: anum4,aden4,bnum4,bden4
      real :: ancg, ancgh, ancgf, ancgpf
      real :: ancgpv, ancghpv, ancgfpv, ancgpfpv
      real :: theta, pp, beta, am3, dSdpi
      integer :: ibound, ia, itercg
c
c     write(6,111)
111   format(' Hi from force')
c
      dSdpi = 0.0
c
c uncomment this line to quench the fermions!
c     return
c pseudofermion action is
c   Phi^dagger {MdaggerM(1)}^1/4 {MdaggerM(m)})^-1/2 {MdaggerM(1)}^1/4 Phi
c
      do ia=1,Nf
c
         X2 = Phi(:, :, :, :, :, ia)

         call qmrherm(X2,res1,itercg,One,1,anum4,aden4,ndiag,
     &        1,isweep,iter)
         ancgpv=ancgpv+float(itercg)

         X2 = X1
c
         call qmrherm(X2,res1,itercg,am,imass,bnum2,bden2,ndiag,
     &        0,isweep,iter)
         ancg=ancg+float(itercg)
c     write(111,*) itercg
         X2 = X1
c
c  evaluates -X2dagger * d/dpi[{MdaggerM(m)}^1/2] * X2
         call qmrherm(X2,res1,itercg,am,imass,anum2,aden2,ndiag,
     &        2,isweep,iter)
         ancgf=ancgf+float(itercg)

c     write(113,*) itercg
c  evaluates +2Re{Phidagger * d/dpi[{MdaggerM(1)}^1/4] * X2}
         call qmrherm(X2,res1,itercg,One,1,anum4,aden4,ndiag,
     &        3,isweep,iter)
         ancgfpv=ancgfpv+float(itercg)
c
      enddo
c
      if(ibound.eq.-1)then
         dSdpi(:, :, ksizet, 3) = -dSdpi(:, :, ksizet, 3)
      endif
c
      dSdpi = dSdpi + beta * Nf * theta
c
      return
      end subroutine
c******************************************************************
c   Evaluation of Hamiltonian function
c******************************************************************
      subroutine hamilton(Phi,
     &       h,hg,hp,s,res2,isweep,iflag,am,imass)
      complex(dp), intent(in) :: Phi(kthird, 0:ksize+1, 0:ksize+1, 
     &                              0:ksizet+1, 4, Nf)
      real(dp), intent(out) :: h, hg, hp, s
      real, intent(in) :: res2, am
      integer, intent(in) :: isweep, iflag, imass
      integer, parameter :: ndiag=25
      common/trial/u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3),
     &     theta(ksize, ksize, ksizet, 3),
     &     pp(ksize, ksize, ksizet, 3)
      common/remez2/anum2(0:ndiag),aden2(ndiag),
     &              bnum2(0:ndiag),bden2(ndiag)
      common/remez4/anum4(0:ndiag),aden4(ndiag),
     &              bnum4(0:ndiag),bden4(ndiag)
      common/param/ancg,ancgh,ancgf,ancgpf
      common/parampv/ancgpv,ancghpv,ancgfpv,ancgpfpv
      common /para/beta,am3,ibound
      common/vector/X1(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      common/dum1/R(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4),
     &     ps(0:ksize+1, 0:ksize+1, 0:ksizet+1, 2)
c     complex, intent(in) :: Phi(kthird, ksize, ksize, ksizet, 4, Nf)
c     complex X1,R
c     complex u
      complex(dp) :: X1,R
      complex(dp) :: u
      real(dp) :: hf
      real(dp) :: anum2,aden2,bnum2,bden2
      real(dp) :: anum4,aden4,bnum4,bden4
      real :: ancg, ancgh, ancgf, ancgpf, 
     &     ancgpv, ancgfpv, ancghpv, ancgpfpv, ancgma
      real :: pp, ps, beta, theta, am3
      integer :: itercg, ia, ibound
c     write(6,111)
111   format(' Hi from hamilton')
c
      hf=0.0
c
      hp = 0.5 * sum(pp ** 2)

!      print *, pp, theta

      hg = 0.5 * Nf * beta * sum(theta ** 2)
      h = hg + hp
c 
c uncomment these lines to quench the fermions!
c     write(6,*) isweep,':  hg', hg,'   hp', hp,'   h',h
c     return
c         
c  pseudofermion action is
c   Phi^dagger {MdaggerM(1)}^1/4 {MdaggerM(m)})^-1/2 {MdaggerM(1)}^1/4 Phi
c
      do ia = 1,Nf
c
         R = Phi(:, :, :, :, :, ia)

         call qmrherm(R,res2,itercg,One,1,anum4,aden4,ndiag,
     &        0,isweep,iflag)
         ancghpv=ancghpv+float(itercg)
c
         R = X1
c
         call qmrherm(R,res2,itercg,am,imass,bnum2,bden2,ndiag,
     &        0,isweep,iflag)
         ancgh=ancgh+float(itercg)
c
         hf = hf + sum(real(conjg(R(:, 1:ksize, 1:ksize, 1:ksizet, :))
     &        * X1(:, 1:ksize, 1:ksize, 1:ksizet, :)))
c
      enddo
c
      h = hg + hp + hf
c     write(6,*) isweep,':  hg', hg,'   hp', hp,'   hf', hf,
c    &   '   h',h
      s = hg + hf
c
      return
      end subroutine

c******************************************************************
c    multisolver matrix inversion via Lanczos technique
c  eg. Golub & van Loan "Matrix Computations" 9.3.1
c       solves (MdaggerM+diag)*x=Phi for ndiag different values of diag
c   iflag=0: simply evaluates X = {MdaggerM}^p * Phi
c   can be interchanged with congrad for p=-1
c   iflag=1: in addition updates Phi0 register needed for PV force term
c   iflag=2: evaluates DWF force term
c   iflag=3: evaluates PV force term
c*****************************************************************m
      subroutine qmrherm(Phi,res,itercg,am,imass,anum,aden,ndiag,
     &                   iflag,isweep,iter)
      implicit none
      complex(dp), intent(in) :: Phi(kthird, 0:ksize+1, 
     &                              0:ksize+1, 0:ksizet+1, 4)
      integer, intent(in) :: imass, ndiag, iflag, isweep, iter
      real(dp), intent(in) :: anum(0:ndiag), aden(ndiag)
      real, intent(in) :: res, am
      integer, intent(out) :: itercg
c
      integer, parameter :: niterc=7500
      common/trial/u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3),
     &             theta(kthird, ksize, ksize, ksizet, 3),
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
c     complex, intent(in) :: Phi(kthird, ksize, ksize, ksizet, 4)
c     complex :: x,u,Phi0
      real :: alphatild
      real(dp) :: coeff
c      
c     complex :: vtild(kthird, ksize, ksize, ksizet, 4)
c     complex :: q(kthird, ksize, ksize, ksizet, 4)
c     complex :: pm1(kthird, ksize, ksize, ksizet, 4,ndiag)
c     complex :: qm1(kthird, ksize, ksize, ksizet, 4)
c     complex :: p(kthird, ksize, ksize, ksizet, 4,ndiag)
c     complex :: x3(kthird, ksize, ksize, ksizet, 4)
c     complex :: R(kthird, ksize, ksize, ksizet, 4)
c     complex x1(kthird, ksize, ksize, ksizet, 4,ndiag)
c     complex :: x2(kthird, ksize, ksize, ksizet, 4)
c     real alpha(ndiag),beta
c     real amu(ndiag),d(ndiag),dm1(ndiag),rho(ndiag),rhom1(ndiag)
      complex(dp) :: vtild(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex(dp) :: q(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex(dp) :: pm1(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4, 
     &                  ndiag)
      complex(dp) :: qm1(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex(dp) :: p(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1,4,ndiag)
      complex(dp) :: x3(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex(dp) :: R(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex(dp) :: x1(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4, 
     &                 ndiag)
      complex(dp) :: x2(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      real(dp) :: alpha(ndiag), beta, beta0, phimod
      real(dp) :: amu(ndiag), d(ndiag), dm1(ndiag)
      real(dp) :: rho(ndiag), rhom1(ndiag)
c      
      real :: resid, rhomax, arelax
      integer :: niter, idiag, ix, iy, it
c
c     write(6,111)
111   format(' Hi from qmrherm')
c
      resid=sqrt(kthird*ksize*ksize*ksizet*4*res*res)
c     write(6,*) iflag, resid
      itercg=0
c
c   initialise r=Phi
c
      R = Phi
      qm1 = cmplx(0.0, 0.0)
      x = anum(0) * Phi

      beta = sqrt(sum(abs(R(:,1:ksize,1:ksize,1:ksizet,:)) ** 2))
      phimod=beta
c     write(6,*) '|| Phi || = ', phimod
c
      do niter=1,niterc
      itercg=itercg+1
c
c  Lanczos steps
c
      q = R / beta
      
      call dslash(vtild,q,u,am,imass)
      call update_halo_5(4, vtild)
      call dslashd(x3,vtild,u,am,imass)
      call update_halo_5(4, x3)
c
      alphatild = sum(real(conjg(q(:,1:ksize,1:ksize,1:ksizet,:)) 
     &                * x3(:,1:ksize,1:ksize,1:ksizet,:)))
c
      R = x3 - alphatild * q - beta * qm1
      qm1 = q
c
      beta0=beta
      beta = sqrt(sum(abs(R(:,1:ksize,1:ksize,1:ksizet,:)) ** 2))
c
      alpha = alphatild + aden
c
      if(niter.eq.1)then
         d = alpha
         rho = beta0 / alpha
         rhom1 = rho
         do idiag = 1, ndiag
            p(:, :, :, :, :, idiag) = q
            pm1(:, :, :, :, :, idiag) = q
            x1(:, :, :, :, :, idiag) = rho(idiag) * q
         enddo
      else
         amu = beta0 / d
         dm1 = d
         d = alpha - beta0 * amu
         rho = -amu * dm1 * rhom1 / d
         do idiag = 1, ndiag
            p(:, :, :, :, :, idiag) =
     &           q - amu(idiag) * pm1(:, :, :, :, :, idiag)
         enddo
         pm1 = p
c     Convergence criterion (a bit ad hoc for now...)
         rhomax = maxval(abs(phimod * rho))
         rhom1 = rho
         do idiag = 1, ndiag
            x1(:, :, :, :, :, idiag) = x1(:, :, :, :, :, idiag) 
     &           + rho(idiag) * p(:, :, :, :, :, idiag)
         enddo
         
c     check to see whether the residual is acceptable for all ndiag....
c     criterion is a bit ad hoc -- relaxing by a factor arelax improves code
c     stability and leads to quicker convergence
         arelax=2.0
         if(rhomax .lt. arelax * resid) then
c     if(rhomax.lt.resid) then
c     call testinv(Phi,resmax,itercg,am,imass,x1,aden,ndiag)
c     convergence based on || residual || not working well in single precision...
c     if(resmax.lt.resid) goto 8
            exit
         endif
      endif
c     
c     end of loop over iter
      enddo
      if (niter .gt. niterc) then
         write(7,*) 'QMRniterc!, isweep,iter,iflag,imass,anum,ndiag = '
     &        ,isweep, iter, iflag, imass, anum(0), ndiag
      end if
c     
      if(iflag.lt.2)then
c     Now evaluate solution x=(MdaggerM)^p * Phi
         do idiag=1,ndiag
            x = x + anum(idiag) * x1(:, :, :, :, :, idiag)
         enddo
c     
c  update phi0 block if required...
         if(iflag.eq.1) then
            Phi0(:, :, :, :, :, 1:ndiag) = X1(:, :, :, :, :, 1:ndiag)
         endif
c     
      else
c
      do idiag=1, ndiag
c
c  X2 = M*X1
        R = X1(:, :, :, :, :, idiag)
        call dslash(X2, R, u, am, imass)
        call update_halo_5(4, X2)
c
        if(iflag.eq.2)then
          coeff=anum(idiag)
          call derivs(R, X2, coeff, 0)
        else
          coeff=-anum(idiag)
          R = Phi0(:, :, :, :, :, idiag)
          call derivs(R, X2, coeff, 0)
c
          call dslash(X2, R, u, am, imass)
          call update_halo_5(4, X2)
c
          R = x1(: ,:, :, :, :, idiag)
          call derivs(X2, R, coeff, 1)
        endif
c
      enddo
      endif
c
c
      return
      end subroutine
c**********************************************************************
c  iflag = 0 : evaluates Rdagger*(Mdagger)'*X2
c  iflag = 1 : evaluates Rdagger*(M)'*X2
c**********************************************************************
      subroutine derivs(R,X2,anum,iflag)
      implicit none
c      complex, intent(in) :: R(kthird, 0:ksize+1, 0:ksize+1, 
c    &                                  0:ksizet+1, 4)
c      complex, intent(in) :: X2(kthird, 0:ksize+1, 0:ksize+1, 
c    &                                   0:ksizet+1, 4)

      complex(dp), intent(in) :: R(kthird, 0:ksize+1, 0:ksize+1, 
     &                            0:ksizet+1, 4)
      complex(dp), intent(in) :: X2(kthird, 0:ksize+1, 0:ksize+1, 
     &                            0:ksizet+1, 4)
      real(dp), intent(in) :: anum
      integer, intent(in) :: iflag

      
      common/dirac/gamval(6,4),gamin(6,4)
      common/gforce/dSdpi(ksize,ksize,ksizet,3)

c      complex :: gamval
      complex(dp) :: gamval
      integer :: gamin
      real :: dSdpi

c      complex(dp) :: tzi
      real(dp) :: tzi_real
      integer :: ix, iy, it, ixup, iyup, itup, idirac, ithird, mu
      integer :: igork1
c
c     write(6,111)
111   format(' Hi from derivs')

c     dSdpi=dSdpi-Re(Rdagger *(d(Mdagger)dp)* X2)
c     Cf. Montvay & Muenster (7.215)
c      tzi=cmplx(0.0,2*anum)
      tzi_real = 2 * anum
c     factor of 2 picks up second term in M&M (7.215)
c
      do mu = 1,3
      ixup = kdelta(1, mu)
      iyup = kdelta(2, mu)
      itup = kdelta(3, mu)

      do idirac=1,4
c      do ithird=1,kthird
c
      do it = 1,ksizet
      do iy = 1,ksize
      do ix = 1,ksize
      dSdpi(ix,iy,it,mu)=
     &     dSdpi(ix,iy,it,mu) + tzi_real * akappa * sum(dimag(
     &conjg(R(:,ix,iy,it,idirac))*
     & X2(:,ix+ixup,iy+iyup,it+itup,idirac))
     &-dimag(conjg(R(:,ix+ixup,iy+iyup,it+itup,idirac))*
     &  X2(:,ix,iy,it,idirac)))
      enddo
      enddo
      enddo
c
      igork1=gamin(mu,idirac)
      if(iflag.eq.0)then
      do it = 1,ksizet
      do iy = 1,ksize
      do ix = 1,ksize
      dSdpi(ix,iy,it,mu)=
     &     dSdpi(ix,iy,it,mu)+ tzi_real * sum(dimag(gamval(mu,idirac)*
     &(conjg(R(:,ix,iy,it,idirac))*
     &        X2(:, ix+ixup,iy+iyup,it+itup,igork1)
     &+conjg(R(:,ix+ixup,iy+iyup,it+itup,idirac))*
     &             X2(:,ix,iy,it,igork1))))
      enddo
      enddo
      enddo
      else
      do it = 1,ksizet
      do iy = 1,ksize
      do ix = 1,ksize
      dSdpi(ix,iy,it,mu)=
     &     dSdpi(ix,iy,it,mu)- tzi_real * sum(dimag(gamval(mu,idirac)*
     &(conjg(R(:,ix,iy,it,idirac))*
     &        X2(:,ix+ixup,iy+iyup,it+itup,igork1)
     &+conjg(R(:,ix+ixup,iy+iyup,it+itup,idirac))*
     &             X2(:,ix,iy,it,igork1))))
      enddo
      enddo
      enddo
      endif
c
      enddo
c      enddo
      enddo
c
      return
      end subroutine
C c******************************************************************
C c   Calculates residual for testing purposes....
C c   needs to run with double precision vectors to be useful.....
C c******************************************************************
C       subroutine testinv(Phi,resmax,itercg,am,imass,x,aden,ndiag)
C       parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
C       parameter(kferm=4*kthird*kvol)
C       common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
C       common/para/bbb,am3,ibound
C       complex Phi(kferm)
C       complex(dp) vtild(kferm)
C       complex(dp) x3(kferm)
C       complex(dp) x(kferm,ndiag),x1(kferm),x2(kferm)
C       complex(dp) u
C c     complex vtild(kferm)
C c     complex x3(kferm)
C c     complex x(kferm,ndiag),x1(kferm),x2(kferm)
C c      complex u
C       real(dp) residual
C       real(dp) aden(ndiag)
C c
C       write(6,111)
C 111   format(' Hi from testinv')
C c
C       resmax=0.0
C c
C       do idiag=1,ndiag
C       residual=0.0
C       do i=1,kferm
C       x3(i)=x(i,idiag)
C       enddo
C       call dslash(x2,x3,u,am,imass)
C       call dslashd(x1,x2,u,am,imass)
C       do i=1,kferm
C       vtild(i)=x1(i)+aden(idiag)*x3(i)-Phi(i)
C       residual=residual+conjg(vtild(i))*vtild(i)
C       enddo
C c     residual=sqrt(residual)
C       if(residual.gt.resmax) resmax=residual
C c
C       write(6,*) idiag, 'itercg = ',itercg, ' residual = ',residual
C       enddo
C c
C       resmax=sqrt(resmax)
C c
C       return
C       end
c******************************************************************
c    matrix inversion via conjugate gradient algorithm
c       solves (Mdagger)Mx=Phi, 
c           NB. no even/odd partitioning
c******************************************************************
      subroutine congrad(Phi,res,itercg,am,imass)
      complex(dp), intent(in) :: 
     &     Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
c     complex, intent(in) ::
c    &     Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      real, intent(in) :: res, am
      integer, intent(out) :: itercg
      integer, intent(in) :: imass

      integer, parameter :: niterc=kthird*kvol
      common/trial/u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3),
     &     thetas(ksize, ksize, ksizet, 3),
     &     pp(ksize, ksize, ksizet, 3)
      real :: thetas, pp
      common/para/bbb,am3,ibound
      real :: bbb, am3
      integer :: ibound
      common/vector/x(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
c     complex x,u
c     complex x1(kferm),x2(kferm),p(kferm),r(kferm)
      complex(dp) x,u
      complex(dp) :: x1(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex(dp) :: x2(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex(dp) :: p(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex(dp) :: r(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      real :: resid
      real :: beta, betan, betad, alpha, alphan, alphad
      integer :: nx
c     write(6,111)
111   format(' Hi from congrad')
c
      resid = 4 * ksize * ksize * ksizet * kthird * res * res
      itercg = 0
      alphan = 0.0
c
      do nx=1,niterc
         itercg=itercg+1
         if(nx.gt.1) goto 51
c
c   initialise p=x, r=Phi(na)
         p = x
         r = Phi
c
         betad=1.0
         alpha=1.0
 51      alphad=0.0
c
c  x1=Mp
         call dslash(x1,p,u,am,imass)
         call update_halo_5(4, x1)
c
         if(nx.ne.1)then
c
c   alpha=(r,r)/(p,(Mdagger)Mp)
            alphad = sum(abs(x1(:, 1:ksize, 1:ksize, 1:ksize, :)) ** 2)
            alpha = alphan / alphad
c     
c   x=x+alpha*p
            x = x + alpha * p
         end if
c     
c   x2=(Mdagger)x1=(Mdagger)Mp
         call dslashd(x2, x1, u, am, imass)
         call update_halo_5(4, x2)
c
c   r=r-alpha*(Mdagger)Mp
         r = r - alpha * x2
c
c   beta=(r_k+1,r_k+1)/(r_k,r_k)
         betan = sum(abs(r(:, 1:ksize, 1:ksize, 1:ksizet, :)) ** 2)
         beta = betan / betad
         betad = betan
         alphan = betan
c
         if(nx.eq.1) beta=0.0
c
c   p=r+beta*p
         p = r + beta * p
         if(betan.lt.resid) exit
      end do
c     write(6,1000)
      if (nx.gt.niterc) then
         write(7,1000)
 1000    format(' # iterations of congrad exceeds niterc')
      end if
c     write(6,*) itercg
      return
      end subroutine

c*****************************************************************
c   Calculate fermion expectation values via a noisy estimator
c   -matrix inversion via conjugate gradient algorithm
c       solves Mx=x1
c     (Numerical Recipes section 2.10 pp.70-73)   
c*******************************************************************
      subroutine measure(psibarpsi, res, aviter, am, imass)
      real, intent(out) :: psibarpsi, aviter
      real, intent(in) :: res, am
      integer, intent(in) :: imass
      integer, parameter :: knoise = 10
      common/trial/u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3),
     &     thetat(ksize, ksize, ksizet, 3),
     &     pp(ksize, ksize, ksizet, 3)
      real :: thetat, pp
      common/para/beta,am3,ibound
      real :: beta, am3
      integer :: ibound
      common/dirac/gamval(6,4),gamin(6,4)
      common/vector/xi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
!      common/ranseed/idum
!      integer(k4b) :: idum
c     complex :: x(0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
c     complex :: Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
c     complex :: xi,gamval
c     complex :: psibarpsi1,psibarpsi2
c     complex :: u
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
c     write(6,*) 'hi from measure'
c
      iter=0
c     pbp=0.0
      cnum(0)=0.0
      cnum(1)=1.0
      cden(1)=0.0
c
      do inoise=1,knoise
c
c     set up noise
      call gauss0(ps)
      psibarpsi1=(0.0,0.0)
      call gauss0(pt)
      psibarpsi2=(0.0,0.0)
c
      do idsource=1,2
c
c  source on domain wall at ithird=1
c   
         x = cmplx(0.0, 0.0)
         if(imass.ne.5)then
            x(:, :, :, idsource) = cmplx(ps(:,:,:,1), ps(:,:,:,2))
         else
            x(:, :, :, idsource) = cmplx(ps(:,:,:,1), ps(:,:,:,2))
         endif
c
         xi = cmplx(0.0, 0.0)
         if(imass.ne.5)then
            xi(1, :, :, :, :) = x
         else
c     xi = 0.5(1+gamma_4)*gamma_5*eta on DW at ithird=1
            do idirac=1,2
               xi(1, :, :, :, idirac) = -x(:, :, :, idirac+2)
            enddo
c     xi = 0.5(1+gamma_4)*eta on DW at ithird=1
         endif
c
c Phi= Mdagger*xi
c
         call dslashd(Phi, xi, u, am, imass)
         call update_halo_5(4, Phi)
c     call qmrherm(Phi,res,itercg,am,imass,cnum,cden,1,0)
         call congrad(Phi, res, itercg, am, imass)
         iter = iter + itercg
c
         if(imass.ne.5)then
c     pbp1 = x^dagger (0.5(1+gamma_4)) xi(kthird)
            psibarpsi1=psibarpsi1
     &           + sum(conjg(x(1:ksize, 1:ksize, 1:ksizet, idsource)) * 
     &               xi(kthird, 1:ksize, 1:ksize, 1:ksizet, idsource))
         else
c     pbp1 = x^dagger (0.5(1-gamma_4)) xi(1)
            psibarpsi1=psibarpsi1
     &           + sum(conjg(x(1:ksize, 1:ksize, 1:ksize, idsource+2))
     &                 * xi(1, 1:ksize, 1:ksize, 1:ksize, idsource+2))
         endif
c
c
c  source on domain wall at ithird=kthird
         idsource2=idsource+2
c
         x = cmplx(0.0, 0.0)
         if(imass.ne.5)then
            x(:, :, :, idsource2) = cmplx(pt(:,:,:,1), pt(:,:,:,2))
         else
            x(:, :, :, idsource2 - 2) = cmplx(pt(:,:,:,1), pt(:,:,:,2))
         endif
c
         xi = cmplx(0.0, 0.0)
         if(imass.ne.5)then
c   xi = 0.5(1-gamma_4)*eta on DW at ithird=kthird
            xi(kthird, :, :, :, :) = x
         else
c   xi = 0.5(1-gamma_4)*gamma_5*eta on DW at ithird=kthird
            xi(kthird, :, :, :, 3:4) = -x(:, :, :, 1:2)
         endif
c     
c Phi= Mdagger*xi
c
         call dslashd(Phi,xi,u,am,imass)
         call update_halo_5(4, Phi)
c
c xi= (M)**-1 * Phi
c
c     call qmrherm(Phi,res,itercg,am,imass,cnum,cden,1,0)
         call congrad(Phi,res,itercg,am,imass)
         iter = iter + itercg
c     
         if(imass.ne.5)then
c pbp2= - x^dagger (0.5(1-gamma_4)) xi(1)
            psibarpsi2=psibarpsi2
     &           +sum(conjg(x(1:ksize, 1:ksize, 1:ksizet, idsource2))
     &                * xi(1, 1:ksize, 1:ksize, 1:ksizet, idsource2))
         else
c pbp2= - x^dagger (0.5(1-gamma_4)) xi(kthird)
            psibarpsi2=psibarpsi2
     &           +sum(conjg(x(1:ksize, 1:ksize, 1:ksizet, idsource))
     &           * xi(kthird, 1:ksize, 1:ksize, 1:ksizet, idsource))
         endif
c
c  end trace on Dirac indices....
      enddo
c
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
c        write(6,*) real(psibarpsi1),aimag(psibarpsi1),
c    &       real(psibarpsi2),aimag(psibarpsi2)
         write(100,*) real(psibarpsi1),aimag(psibarpsi1),
     &       real(psibarpsi2),aimag(psibarpsi2)
c
c end loop on noise
      enddo
c
      psibarpsi=0.0
      susclsing=0.0
c
      psibarpsi = sum(pbp)
      do inoise=1,knoise
         susclsing = susclsing + 
     &        sum(pbp(inoise) * pbp(inoise+1:knoise))
      enddo
      psibarpsi = psibarpsi / knoise
      susclsing = 2 * kvol * susclsing / (knoise * (knoise-1))
      write(200,*) psibarpsi, susclsing
      aviter = float(iter) / (4*knoise)
      return
      end subroutine
C c******************************************************************
C c   Calculate meson correlators using point sources on domain walls
C c   -matrix inversion via conjugate gradient algorithm
C c       solves Mx=x1
C c     (Numerical Recipes section 2.10 pp.70-73)   
C c*******************************************************************
C       subroutine meson(res,itercg,aviter,am,imass)
C       use random
C       parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
C       parameter(ksize2=ksize*ksize)
C       parameter(akappa=0.5)
C       common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
C       common/para/beta,am3,ibound
C       common/dirac/gamval(6,4),gamin(6,4)
C       common /neighb/id(kvol,3),iu(kvol,3)
C       common/vector/xi(kthird,kvol,4)
C       common/ranseed/idum
C       common/v/v(97)
C c     complex x(kvol,4),x0(kvol,4),Phi(kthird,kvol,4)
C c     complex xi,gamval
C c     complex prop00(kvol,3:4,1:2),prop0L(kvol,3:4,3:4)
C       complex(dp) x(kvol,4),x0(kvol,4),Phi(kthird,kvol,4)
C       complex(dp) xi,gamval
C       complex(dp) prop00(kvol,3:4,1:2),prop0L(kvol,3:4,3:4)
C c     complex prop00n(kvol,3:4,1:2),prop0Ln(kvol,3:4,3:4)
C       real cpm(0:ksizet-1),cmm(0:ksizet-1)
C c     complex cpmn(0:ksizet-1),cmmn(0:ksizet-1)
C c     complex cferm1(0:ksizet-1), cferm2(0:ksizet-1)
C c     complex u
C       complex(dp) cferm1(0:ksizet-1), cferm2(0:ksizet-1)
C       complex(dp) u
C       real ps(kvol,2)
C       real ran
C       integer gamin
C c     write(6,*) 'hi from meson'
C c      
C       nsource=5
C       nsmear=10
C       c=0.25
C       iter=0
C c
C       do it=0,ksizet-1
C       cpm(it)=(0.0,0.0)
C       cmm(it)=(0.0,0.0)
C c     cpmn(it)=(0.0,0.0)
C c     cmmn(it)=(0.0,0.0)
C       cferm1(it)=(0.0,0.0)
C       cferm2(it)=(0.0,0.0)
C       enddo
C c
C c      susceptibility
C       chim=0.0
C       chip=0.0
C c
C       do ksource=1,nsource
C c
C c   random location for +m source
C       ixxx=int(ksize*rano(yran,idum))+1
C       iyyy=int(ksize*rano(yran,idum))+1
C       ittt=int(ksizet*rano(yran,idum))+1
C       isource=ixxx+ksize*((iyyy-1)+ksize*(ittt-1))
C c     write(6,*) ixxx,iyyy,ittt, isource
C c
C c  loop over Dirac index of source...   
C       do idsource=3,4
C c
C c  source on domain wall at ithird=1
C       do 300 ithird=1,kthird
C       if(ithird.eq.1)then
C       do idirac=1,4
C       do i=1,kvol
C       x(i,idirac)=(0.0,0.0)
C       enddo
C c  wall source
C c     ioff=(ittt-1)*ksize2
C c     do i=1,ksize2
C c     x(i+ioff,idirac)=(1.0,0.0)/ksize2
C c     enddo
C c
C       enddo
C c  point source at fixed site, spin...
C       x(isource,idsource)=(1.0,0.0)
C c
C c now smear it.....
C c
C       do ismear=1,nsmear
C       call dslash2d(x0,x,u)
C       do idirac=1,4
C       do i=1,kvol
C       x(i,idirac)=(1.0-c)*x(i,idirac)+c*x0(i,idirac)
C       enddo
C       enddo
C       enddo
C c
C       else
C       do idirac=1,4
C       do i=1,kvol
C       xi(ithird,i,idirac)=(0.0,0.0)
C       enddo
C       enddo
C       endif
C 300   continue
C c
C c   xi = x  on DW at ithird=1
C c
C       do idirac=1,4
C       do i=1,kvol
C       xi(1,i,idirac)=x(i,idirac)
C       enddo
C       enddo
C c
C c Phi= Mdagger*xi
C c
C       call dslashd(Phi,xi,u,am,imass)
C c
C c  preconditioning (no,really)
C       call dslashd(xi,Phi,u,am,imass)
C c  
C c xi= (MdaggerM)**-1 * Phi 
C c
C c     call congrad(Phi,res,itercg,am,imass) 
C       iter=iter+itercg
C c
C       do idsink=1,2
C       idsink2=idsink+2
C       do i=1,kvol
C       prop00(i,idsource,idsink)=xi(1,i,idsink)
C       prop0L(i,idsource,idsink2)=xi(kthird,i,idsink2)
C       enddo
C       enddo
C c
C c     if(imass.ne.1)then
C c  now evaluate with sign of mass reversed (not needed for hermitian mass term)
C c     am=-am
C c  source on domain wall at ithird=1
C c     do 400 ithird=1,kthird
C c     if(ithird.eq.1)then
C c     do idirac=1,4
C c     do i=1,kvol
C c     xi(1,i,idirac)=x(i,idirac)
C c     enddo
C c     enddo
C c     else
C c     do idirac=1,4
C c     do i=1,kvol
C c     xi(ithird,i,idirac)=(0.0,0.0)
C c     enddo
C c     enddo
C c     endif
C 400   continue
C c
C c Phi= Mdagger*xi
C c
C c     call dslashd(Phi,xi,u,am,imass)
C c
C c     call dslashd(xi,Phi,u,am,imass)
C c  
C c xi= (MdaggerM)**-1 * Phi 
C c
C c     call congrad(Phi,res,itercg,am,imass) 
C c     iter=iter+itercg
C c
C c     do idsink=1,2
C c     idsink2=idsink+2
C c     do i=1,kvol
C c     prop00n(i,idsource,idsink)=xi(1,i,idsink)
C c     prop0Ln(i,idsource,idsink2)=xi(kthird,i,idsink2)
C c     enddo
C c     enddo
C c
C c     am=-am
C c     endif
C c
C c  end loop on source Dirac index....
C       enddo
C c
C c  Now tie up the ends....
C c
C c  First C+-
C c
C c  now evaluate the trace (exploiting projection)
C       do id1=3,4
C       do id2=1,2
C       do it=0,ksizet-1
C       itt=mod((ittt+it-1),ksizet)+1
C       ioff=(itt-1)*ksize2
C       do i=1,ksize2
C       cpm(it)=cpm(it)
C      &        +prop00(i+ioff,id1,id2)*conjg(prop00(i+ioff,id1,id2))
C       enddo
C       enddo
C       enddo
C       enddo
C c
C c     if(imass.ne.1)then
C c     do id1=3,4
C c     do id2=1,2
C c     do it=0,ksizet-1
C c     itt=mod((ittt+it-1),ksizet)+1
C c     ioff=(itt-1)*ksize2
C c     do i=1,ksize2
C c     cpmn(it)=cpmn(it)
C c    &        +prop00(i+ioff,id1,id2)*conjg(prop00n(i+ioff,id1,id2))
C c     enddo
C c     enddo
C c     enddo
C c     enddo
C c     endif
C c
C c  next C--
C c  now evaluate the trace exploiting projection
C       do id1=3,4
C       do id2=3,4
C       do it=0,ksizet-1
C       itt=mod((ittt+it-1),ksizet)+1
C       ioff=(itt-1)*ksize2
C       do i=1,ksize2
C       cmm(it)=cmm(it)
C      &   +prop0L(i+ioff,id1,id2)*conjg(prop0L(i+ioff,id1,id2))
C       enddo
C       enddo
C       enddo
C       enddo
C c
C c     if(imass.ne.1)then
C c     do id1=3,4
C c     do id2=3,4
C c     do it=0,ksizet-1
C c     itt=mod((ittt+it-1),ksizet)+1
C c     ioff=(itt-1)*ksize2
C c     do i=1,ksize2
C c     cmmn(it)=cmmn(it)
C c    &   +prop0L(i+ioff,id1,id2)*conjg(prop0Ln(i+ioff,id1,id2))
C c     enddo
C c     enddo
C c     enddo
C c     enddo
C c     endif
C c
C c    now the fermion propagator
C c  = tr{ P_-*Psi(0,1)Psibar(x,Ls) + gamma_0*P_-*Psi(0,1)Psibar(x,1) }
C       do idd=3,4
C       do it=0,ksizet-1
C       itt=mod((ittt+it-1),ksizet)+1
C c correct for apbc
C       if(itt.ge.ittt)then
C         isign=1
C       else
C         isign=ibound
C       endif
C c
C       ioff=(itt-1)*ksize2
C       do i=1,ksize2
C       cferm1(it)=cferm1(it)
C      & +isign*akappa*prop0L(i+ioff,idd,idd)
C       cferm2(it)=cferm2(it)
C      & +isign*gamval(3,idd)*prop00(i+ioff,idd,gamin(3,idd))
C       enddo
C       enddo
C       enddo
C c
C c  finish loop over sources
C       enddo
C c
C       do it=0,ksizet-1
C       cpm(it)=cpm(it)/nsource
C       cmm(it)=cmm(it)/nsource
C c  Cf. (54) of 1507.07717
C       chim=chim+2*(cpm(it)+cmm(it))
C       enddo
C c     if(imass.ne.1)then
C c       if(imass.eq.3)then
C c         do it=0,ksizet-1
C c           cpmn(it)=cpmn(it)/nsource
C c           cmmn(it)=cmmn(it)/nsource
C c  Cf. (54),(61) of 1507.07717
C c           chip=chip-2*(cpmn(it)-cmmn(it))
C c         enddo
C c       else
C c         do it=0,ksizet-1
C c           cpmn(it)=cpmn(it)/nsource
C c           cmmn(it)=cmmn(it)/nsource
C c  Cf. (64),(65) of 1507.07717
C c           chip=chip-2*(cpm(it)-cmm(it))
C c         enddo
C c       endif
C c     endif
C c
C       do it=0,ksizet-1
C       write(302,*) it, cpm(it), cmm(it)
C       write(500,*) it, real(cferm1(it)), aimag(cferm1(it))
C       write(501,*) it, real(cferm2(it)), aimag(cferm2(it))
C       enddo
C c     write(6,*) chim
C       write(400,*) chim
C c     if(imass.ne.1)then
C c     do it=0,ksizet-1
C c     write(402,*) it, real(cpmn(it)), real(cmmn(it))
C c     write(403,*) it, aimag(cpmn(it)), aimag(cmmn(it))
C c     enddo
C c     write(401,*) chip
C c     endif
C c
C c     if(imass.eq.1)then
C       aviter=float(iter)/(2*nsource)
C c     else
C c     aviter=float(iter)/(4*nsource)
C c     endif
C c
C       return
C       end
c*******************************************************************
c
      subroutine sread
      use random
      implicit none
      common/gauge/ theta(ksize, ksize, ksizet, 3)
!     common/ranseed/ idum
      real :: theta
!      integer(k4b) :: idum
      open(unit=10,file='con',
     1     status='unknown',form='unformatted')
      read (10) theta, seed
!      idum = -idum
      close(10)
      return
      end subroutine
c
      subroutine swrite
      use random
      implicit none
      common/gauge/ theta(ksize, ksize, ksizet, 3)
!      common/ranseed/ idum
      real :: theta
!      integer(k4b) :: idum
      open(unit=31,file='con',
     1     status='unknown',form='unformatted')
      write (31) theta, seed
      close(31)
      return
      end subroutine
c
      subroutine init(nc)
      use random
c*******************************************************************
c     sets initial values
c     nc=0 cold start
c     nc=1 hot start
c     nc<0 no initialization
c*******************************************************************
      implicit none
      integer, intent(in) :: nc
      common/gauge/theta(ksize, ksize, ksizet, 3), seed
      common/dirac/gamval(6,4),gamin(6,4)
!      common/ranseed/idum
c     complex gamval,one,zi
      complex(dp) :: gamval,one,zi
      real :: theta
      real(dp) :: seed
      integer :: gamin
!      integer(k4b) :: idum
      integer :: ix, iy, it, mu
      real :: g
c
c
      one=(1.0,0.0)
      zi=(0.0,1.0)
c*******************************************************************
c  calculate constants
c*******************************************************************
C      call addrc
c*******************************************************************
c    setup Dirac algebra
c*******************************************************************
c
c     gamma_1
c
      gamval(1,1)=-zi
      gamval(1,2)=-zi
      gamval(1,3)= zi
      gamval(1,4)= zi
c
      gamin(1,1)=4
      gamin(1,2)=3
      gamin(1,3)=2
      gamin(1,4)=1
c
c     gamma_2
c
      gamval(2,1)=-one
      gamval(2,2)= one
      gamval(2,3)= one
      gamval(2,4)=-one
c
      gamin(2,1)=4
      gamin(2,2)=3
      gamin(2,3)=2
      gamin(2,4)=1
c
c     gamma_3
c
      gamval(3,1)=-zi
      gamval(3,2)= zi
      gamval(3,3)= zi
      gamval(3,4)=-zi
c
      gamin(3,1)=3
      gamin(3,2)=4
      gamin(3,3)=1
      gamin(3,4)=2
c
c     gamma_4
c
      gamval(4,1)= one
      gamval(4,2)= one
      gamval(4,3)= -one
      gamval(4,4)= -one
c
      gamin(4,1)=1
      gamin(4,2)=2
      gamin(4,3)=3
      gamin(4,4)=4
c
c     gamma_5 = gamma_1 * gamma_2 * gamma_3 * gamma_4
c
      gamval(5,1)=-one
      gamval(5,2)=-one
      gamval(5,3)=-one
      gamval(5,4)=-one
c
      gamin(5,1)=3
      gamin(5,2)=4
      gamin(5,3)=1
      gamin(5,4)=2
c
c     gamma_4 * gamma_5 (called gamma_3 gamma_5 in notes)
      gamval(6,1)=-one
      gamval(6,2)=-one
      gamval(6,3)= one
      gamval(6,4)= one
c
      gamin(6,1)=3
      gamin(6,2)=4
      gamin(6,3)=1
      gamin(6,4)=2
c
c
      gamval = gamval * akappa
c
      if(nc.lt.0) return
c
c     initialize gauge fields
c
      if(nc .eq. 1)goto 40
c     (else cold start)
      theta = 0.0
      return
c
40    continue
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
      end subroutine
c******************************************************************
c   calculate compact links from non-compact links
c******************************************************************
      pure subroutine coef(u,theta)
      implicit none
      common/para/beta,am3,ibound
c
      complex(dp), intent(inout) :: u(0:ksize+1, 0:ksize+1, 
     &                               0:ksizet+1, 3)
      real, intent(in) :: theta(ksize, ksize, ksizet, 3)
      real :: beta, am3
      integer :: ibound
      integer :: ix, iy, it, mu
c
c     u(1:ksize, 1:ksize, 1:ksizet, :) = exp(cmplx(0.0, theta))
      u(1:ksize, 1:ksize, 1:ksizet, :) = (1.0 + cmplx(0.0, theta))
c
c  anti-p.b.c. in timelike direction
      if(ibound.eq.-1)then
        u(:, :, ksizet, 3) = -u(:, :, ksizet, 3)
      end if
c      
      call update_halo_4(3, u)
      return
      end subroutine
c**********************************************************************
c calculate vector of gaussian random numbers with unit variance
c to refresh momenta
c   Numerical Recipes pp.203
c**********************************************************************
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
c     write(6,1)
1     format(' Hi from gaussp')
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
      end subroutine
c**********************************************************************
c calculate vector of gaussian random numbers with unit variance
c to generate pseudofermion fields R
c   Numerical Recipes pp.203
c**********************************************************************
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
c     write(6,1)
1     format(' Hi from gauss0')
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
      end subroutine

c*****************************************
c  Random number generator Numerical recipes B7
c
      real function ran(idum)
      implicit none

      integer, parameter :: k4b=selected_int_kind(9)
      integer(k4b), intent(inout) :: idum
      integer(k4b), parameter :: IA=16807, IM=2147483647,
     &                           IQ=127773, IR=2836
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
c***********************************************************************
      pure subroutine dslash(Phi,R,u,am,imass)
c
c     calculates Phi = M*R
c
      implicit none
      common/para/beta,am3,ibound
      common/dirac/gamval(6,4),gamin(6,4)
c     complex, intent(in) :: u(0:ksize+1,0:ksize+1,0:ksizet+1,3)
c     complex, intent(in) :: Phi(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
c     complex, intent(in) :: R(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
c     complex gamval
c     complex zkappa
      complex(dp), intent(in) :: u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3)
      complex(dp), intent(out) :: Phi(kthird, 0:ksize+1,
     &                                 0:ksize+1, 0:ksizet+1, 4)
      complex(dp), intent(in) :: R(kthird, 0:ksize+1, 0:ksize+1,
     &                            0:ksizet+1, 4)
      integer, intent(in) :: imass
      real, intent(in) :: am
      complex(dp) :: gamval
      complex(dp) :: zkappa
      integer :: gamin
      real :: beta, am3, diag
      integer :: ibound
      integer :: ixup, iyup, itup, ix, iy, it, ithird, idirac, mu, igork
c     write(6,*) 'hi from dslash'
c
c     diagonal term
      diag=(3.0-am3)+1.0
      Phi=diag*R
c      
c     Wilson term (hermitian) and Dirac term (antihermitian)
      do mu=1,3
      ixup = kdelta(1, mu)
      iyup = kdelta(2, mu)
      itup = kdelta(3, mu)

      do idirac=1,4
      igork=gamin(mu,idirac)
      do it = 1,ksizet
      do iy = 1,ksize
      do ix = 1,ksize
      Phi(:,ix,iy,it,idirac)=Phi(:,ix,iy,it,idirac)
c Wilson term (hermitian)
     &    -akappa*(u(ix,iy,it,mu)
     &              *R(:, ix+ixup, iy+iyup, it+itup, idirac)
     &             +conjg(u(ix-ixup, iy-iyup, it-itup, mu))
     &              *R(:, ix-ixup, iy-iyup, it-itup, idirac))
c Dirac term (antihermitian)
     &     +gamval(mu,idirac)*
     &       (u(ix,iy,it,mu)
     &         *R(:, ix+ixup, iy+iyup, it+itup, igork)
     &        -conjg(u(ix-ixup, iy-iyup, it-itup, mu))
     &         *R(:, ix-ixup, iy-iyup, it-itup, igork))
      enddo
      enddo
      enddo
      enddo
      enddo
c
c  s-like term exploiting projection
      Phi(1:kthird-1,:,:,:,3:4)=Phi(1:kthird-1,:,:,:,3:4)
     &   -R(2:kthird,:,:,:,3:4)
      Phi(2:kthird,:,:,:,1:2)=Phi(2:kthird,:,:,:,1:2)
     &    -R(1:kthird-1,:,:,:,1:2)
c
c  Mass term (couples the two walls unless imass=5)
      if (imass.eq.1) then
         zkappa=cmplx(am,0.0)
         Phi(kthird, :, :, :, 3:4) = Phi(kthird, :, :, :, 3:4)
     &                               + zkappa * R(1, :, :, :, 3:4)
         Phi(1, :, :, :, 1:2) = Phi(1, :, :, :, 1:2)
     &                          + zkappa * R(kthird, :, :, :, 1:2)
      elseif (imass.eq.3) then
         zkappa=cmplx(0.0,-am)
         Phi(kthird,:, :, :, 3:4) = Phi(kthird, :, :, :, 3:4)
     &                              - zkappa * R(1, :, :, :, 3:4)
         Phi(1, :, :, :, 1:2) = Phi(1, :, :, :, 1:2)
     &                          + zkappa * R(kthird, :, :, :, 1:2)
      elseif (imass.eq.5) then
         zkappa=cmplx(0.0,-am)
c         do idirac=3,4
c         igork=gamin(5,idirac)
         Phi(kthird, :, :, :, 3:4) = Phi(kthird, :, :, :, 3:4)
     &                               - zkappa * R(kthird, :, :, :, 1:2)
c        Phi(kthird,:,:,:,idirac)=Phi(kthird,:,:,:,idirac)
c    &           +2*zkappa*gamval(5,idirac)*R(kthird,:,:,:,igork)
c         enddo
c         do idirac=1,2
c         igork=gamin(5,idirac)
         Phi(1, :, :, :, 1:2) = Phi(1, :, :, :, 1:2)
     &                          - zkappa * R(1, :, :, :, 3:4)
c        Phi(1,:,:,:,idirac)=Phi(1,:,:,:,idirac)
c    &        +2*zkappa*gamval(5,idirac)*R(1,:,:,:,igork)
c         enddo
      endif
c
      return
      end subroutine
c***********************************************************************
      pure subroutine dslashd(Phi,R,u,am,imass)
c
c     calculates Phi = Mdagger*R
c
      implicit none
      common/para/beta,am3,ibound
      common/dirac/gamval(6,4),gamin(6,4)
c     complex, intent(in) ::  u(0:ksize+1,0:ksize+1,0:ksizet+1,3)
c     complex, intent(out) :: Phi(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
c     complex, intent(in) R(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
c     complex gamval
c     complex zkappa
      complex(dp), intent(in) :: u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3)
      complex(dp), intent(out) :: Phi(kthird, 0:ksize+1,
     &                                 0:ksize+1, 0:ksizet+1, 4)
      complex(dp), intent(in) :: R(kthird, 0:ksize+1, 0:ksize+1,
     &                            0:ksizet+1, 4)
      integer, intent(in) :: imass
      real, intent(in) :: am
      complex(dp) :: gamval
      complex(dp) :: zkappa
      integer :: gamin
      real :: beta, am3, diag
      integer :: ibound
      integer :: ixup, iyup, itup, ix, iy, it, ithird, idirac, mu, igork
c     write(6,*) 'hi from dslashd'
c
c     diagonal term (hermitian)
      diag=(3.0-am3)+1.0
      Phi = diag * R
c
c     Wilson term (hermitian) and Dirac term (antihermitian)
      do mu=1,3
      ixup = kdelta(1, mu)
      iyup = kdelta(2, mu)
      itup = kdelta(3, mu)

      do idirac=1,4
      igork=gamin(mu,idirac)
      do it = 1,ksizet
      do iy = 1,ksize
      do ix = 1,ksize
      Phi(:,ix,iy,it,idirac)=Phi(:,ix,iy,it,idirac)
c Wilson term (hermitian)
     &    -akappa*(u(ix,iy,it,mu)
     &              *R(:, ix+ixup, iy+iyup, it+itup, idirac)
     &             +conjg(u(ix-ixup, iy-iyup, it-itup, mu))
     &              *R(:, ix-ixup, iy-iyup, it-itup, idirac))
c Dirac term (antihermitian)
     &     -gamval(mu,idirac)*
     &       (u(ix,iy,it,mu)
     &         *R(:, ix+ixup, iy+iyup, it+itup, igork)
     &        -conjg(u(ix-ixup, iy-iyup, it-itup, mu))
     &         *R(:, ix-ixup, iy-iyup, it-itup, igork))
      enddo
      enddo
      enddo
      enddo
      enddo
c
c  s-like term exploiting projection
      Phi(1:kthird-1,:,:,:,1:2)=Phi(1:kthird-1,:,:,:,1:2)
     &   -R(2:kthird,:,:,:,1:2)
      Phi(2:kthird,:,:,:,3:4)=Phi(2:kthird,:,:,:,3:4)
     &   -R(1:kthird-1,:,:,:,3:4)
c
c  Mass term (couples the two walls unless imass=5) 
      if(imass.eq.1)then
         zkappa=cmplx(am,0.0)
         Phi(kthird,:,:,:,1:2)=Phi(kthird,:,:,:,1:2)
     &                         +zkappa*R(1,:,:,:,1:2)
         Phi(1,:,:,:,3:4)=Phi(1,:,:,:,3:4)+zkappa*R(kthird,:,:,:,3:4)
      elseif(imass.eq.3)then
         zkappa=cmplx(0.0,am)
         Phi(kthird,:,:,:,1:2)=Phi(kthird,:,:,:,1:2)
     &                         +zkappa*R(1,:,:,:,1:2)
         Phi(1,:,:,:,3:4)=Phi(1,:,:,:,3:4)
     &                       -zkappa*R(kthird,:,:,:,3:4)
      elseif(imass.eq.5)then
         zkappa=cmplx(0.0,am)
         Phi(kthird,:,:,:,1:2)=Phi(kthird,:,:,:,1:2)
     &                        -zkappa*R(kthird,:,:,:,3:4)
         Phi(1,:,:,:,3:4)=Phi(1,:,:,:,3:4)-zkappa*R(1,:,:,:,1:2)
      endif
!      call update_halo_5(4, Phi)
c
      return
      end subroutine

c***********************************************************************
      pure subroutine dslash2d(Phi,R,u)
c
c     calculates Phi = M*R
c
      implicit none
C      integer :: kdelta
      common/dirac/gamval(6,4),gamin(6,4)
c     complex u(ksize,ksize,ksizet,3)
c     complex Phi(ksize,ksize,ksizet,4),R(ksize,ksize,ksizet,4)
c     complex gamval
      complex(dp), intent(in) ::  u(0:ksize+1,0:ksize+1,0:ksizet+1,3)
      complex(dp), intent(out) :: 
     &     Phi(0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex(dp), intent(in) :: R(0:ksize+1,0:ksize+1,0:ksizet+1,4)
      complex(dp) :: gamval
      integer :: gamin
      integer :: ix, iy, it, idirac, mu, ixup, iyup, itup, igork
      real :: diag
      
c     write(6,*) 'hi from dslash2d'
c
c     diagonal term
      diag=2.0
      Phi = diag * R
      
c     Wilson and Dirac terms
      do mu=1,2
      ixup=kdelta(1,mu)
      iyup=kdelta(2,mu)
c
      do idirac=1,4
      igork=gamin(mu,idirac)
      do it=1,ksizet
      do iy=1,ksize
      do ix=1,ksize
      Phi(ix,iy,it,idirac)=
c Wilson term
     &    Phi(ix,iy,it,idirac)
     &    -akappa*(u(ix,iy,it,mu)*R(ix+ixup, iy+iyup, it, idirac)
     &             +conjg(u(ix-ixup, iy-iyup, it, mu))
     &              *R(ix-ixup, iy-iyup, it, idirac))
c Dirac term
     &     +gamval(mu,idirac)*
     &      (u(ix,iy,it,mu)*R(ix+ixup, iy+iyup, it, igork)
     &       -conjg(u(ix-ixup, iy-iyup, it,mu))
     &        *R(ix-ixup, iy-iyup, it, igork))
      enddo
      enddo
      enddo
      enddo
      enddo
      call update_halo_4(4, Phi)
c
      return
      end subroutine
c
c***********************************************************************
c   Update boundary terms
c***********************************************************************
      pure subroutine update_halo_4(size4, Array)
c     
      implicit none
c
      integer, intent(in) :: size4
      complex(dp), intent(inout) :: Array(0:ksize+1, 0:ksize+1,
     &                                   0:ksizet+1, size4)
c
      Array(0,:,:,:) = Array(ksize,:,:,:)
      Array(ksize+1,:,:,:) = Array(1,:,:,:)
      Array(:,0,:,:) = Array(:,ksize,:,:)
      Array(:,ksize+1,:,:) = Array(:,1,:,:)
      Array(:,:,0,:) = Array(:,:,ksize,:)
      Array(:,:,ksize+1,:) = Array(:,:,1,:)
c      
      return
c      
      end subroutine
c***********************************************************************
      pure subroutine update_halo_4_real(size4, Array)
c     
      implicit none
c
      integer, intent(in) :: size4
      real, intent(inout) :: Array(0:ksize+1, 0:ksize+1,
     &                             0:ksizet+1, size4)
c
      Array(0,:,:,:) = Array(ksize,:,:,:)
      Array(ksize+1,:,:,:) = Array(1,:,:,:)
      Array(:,0,:,:) = Array(:,ksize,:,:)
      Array(:,ksize+1,:,:) = Array(:,1,:,:)
      Array(:,:,0,:) = Array(:,:,ksize,:)
      Array(:,:,ksize+1,:) = Array(:,:,1,:)
c      
      return
c      
      end subroutine
c***********************************************************************
      pure subroutine update_halo_5(size5, Array)
c     
      implicit none
c
      integer, intent(in) :: size5
      complex(dp), intent(inout) :: Array(kthird, 0:ksize+1, 0:ksize+1,
     &                                   0:ksizet+1, size5)
c
      Array(:,0,:,:,:) = Array(:,ksize,:,:,:)
      Array(:,ksize+1,:,:,:) = Array(:,1,:,:,:)
      Array(:,:,0,:,:) = Array(:,:,ksize,:,:)
      Array(:,:,ksize+1,:,:) = Array(:,:,1,:,:)
      Array(:,:,:,0,:) = Array(:,:,:,ksize,:)
      Array(:,:,:,ksize+1,:) = Array(:,:,:,1,:)
c      
      return
c      
      end subroutine
c***********************************************************************
      pure subroutine update_halo_6(size5, size6, Array)
c     
      implicit none
c
      integer, intent(in) :: size5, size6
      complex(dp), intent(inout) :: Array(kthird, 0:ksize+1, 0:ksize+1,
     &                                   0:ksizet+1, size5, size6)
c
      Array(:,0,:,:,:,:) = Array(:,ksize,:,:,:,:)
      Array(:,ksize+1,:,:,:,:) = Array(:,1,:,:,:,:)
      Array(:,:,0,:,:,:) = Array(:,:,ksize,:,:,:)
      Array(:,:,ksize+1,:,:,:) = Array(:,:,1,:,:,:)
      Array(:,:,:,0,:,:) = Array(:,:,:,ksize,:,:)
      Array(:,:,:,ksize+1,:,:) = Array(:,:,:,1,:,:)
c      
      return
c      
      end subroutine

c***********************************************************************
c   A Kronecker delta function
c   Useful for calculating coordinate offsets
c***********************************************************************
      pure integer function kdelta(nu, mu)
        implicit none
        integer, intent(in) :: nu
        integer, intent(in) :: mu

        kdelta=merge(1,0,nu==mu)
      end function

      end module dwf3d_lib


