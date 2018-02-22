      module purefunctions
      interface
        pure function kdelta(nu,mu)
          integer, intent(in) :: mu
          integer, intent(in) :: nu
          integer :: kdelta
        end function
      end interface
      interface
        pure subroutine update_halo_4(size4, Array)
          integer, parameter :: ksize=12,ksizet=12
          integer, intent(in) :: size4
          complex*16, intent(inout) :: Array(0:ksize+1, 0:ksize+1,
     &                                       0:ksizet+1, size4)
        end subroutine
      end interface
      interface
        pure subroutine update_halo_4_real(size4, Array)
          integer, parameter :: ksize=12,ksizet=12
          integer, intent(in) :: size4
          real, intent(inout) :: Array(0:ksize+1, 0:ksize+1,
     &                                       0:ksizet+1, size4)
        end subroutine
      end interface

      interface
        pure subroutine update_halo_5(size5, Array)
c     
          integer, parameter :: ksize=12,ksizet=12,kthird=24
c
          integer, intent(in) :: size5
          complex*16, intent(inout) :: Array(kthird,0:ksize+1,0:ksize+1,
     &                                       0:ksizet+1, size5)
        end subroutine
      end interface
      end module

C       subroutine dwf3d_main
C c*******************************************************************
C c    Rational Hybrid Monte Carlo algorithm for bulk Thirring Model with Domain Wall
C c         fermions
C c
C c    requires operation of QMR on complex vectors to determine
C c    (Mdagger M)**-1  Phi 
C c
C c    requires input partial fraction coefficients from Remez algorithm
C c    of Clarke & Kennedy
C c
C c    the "third" direction is actually indexed 4 in the code - sorry!
C c      (should have used C I know)
C c
C c           { 1 - hermition mass term psibar psi
C c    imass= { 3 - antih.    mass term i psibar gamma_3 psi
C c           { 5 - antih.    mass term i psibar gamma_5 psi
C c     linear combinations of above require code modification
C c
C c    code exploits fact that gamma_3 is diagonal in Dirac basis used - speeds up 
C c    evolution slightly, and measurement by factor of two. 
C c
C c    Fermion expectation values are measured using a noisy estimator.
C c    on the Wilson matrix, which has dimension 4*kvol*kthird
C c    inversions done using congrad, and matrix multiplies with dslash, dslashd
C c
C c    Pauli-Villars determinant defined using a hermitian mass m_h=One
C c
C c    trajectory length is random with mean dt*iterl
C c    The code runs for a fixed number iter2 of trajectories.
C c
C c    Phi: pseudofermion field 
C c    am: bare fermion mass 
C c    actiona: running average of total action
C c
C c                                               SJH February 2017
C c*******************************************************************
C       parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
C       parameter(respbp=0.000001,rescgg=0.000001,rescga=0.000000001)
C       parameter(rescgm=0.000000001)
C       parameter(itermax=1000)
C       parameter(One=1.0)
C       parameter(ndiag=25,ndiagg=12,Nf=1)
C       common/gauge/theta(kvol,3),seed
C       common/trial/ut(kvol,3),thetat(kvol,3),pp(kvol,3)
C       common /para/beta,am3,ibound
C       common/remez2/anum2(0:ndiag),aden2(ndiag),
C      &              bnum2(0:ndiag),bden2(ndiag)
C       common/remez4/anum4(0:ndiag),aden4(ndiag),
C      &              bnum4(0:ndiag),bden4(ndiag)
C       common/remez2g/anum2g(0:ndiagg),aden2g(ndiagg),
C      &              bnum2g(0:ndiagg),bden2g(ndiagg)
C       common/remez4g/anum4g(0:ndiagg),aden4g(ndiagg),
C      &              bnum4g(0:ndiagg),bden4g(ndiagg)
C       common/gforce/dSdpi(kvol,3)
C       common /neighb/id(kvol,3),iu(kvol,3)
C       common/param/ancg,ancgh,ancgf,ancgpf
C       common/parampv/ancgpv,ancghpv,ancgfpv,ancgpfpv
C       common/trans/tpi 
C       common/dum1/ R(kthird,kvol,4),ps(kvol,2)
C       common/vector/X1(kthird,kvol,4)
C       common/ranseed/idum
C       common/v/v(97)
C c     complex Phi(kthird,kvol,4,Nf),X0(kthird,kvol,4)
C c     complex R,zi,qq,qbqb
C c     complex u,ut,X1
C c     complex a,b
C       complex*16 Phi(kthird,kvol,4,Nf),X0(kthird,kvol,4)
C       complex*16 R,zi,qq,qbqb
C       complex*16 u,ut,X1
C       complex*16 a,b
C       real ran
C       real*8 H0,H1,S0,S1,dH,dS,hg,hp
C       real*8 seed
C       real*8 anum2,aden2,bnum2,bden2
C       real*8 anum4,aden4,bnum4,bden4
C       real*8 anum2g,aden2g,bnum2g,bden2g
C       real*8 anum4g,aden4g,bnum4g,bden4g
C c*******************************************************************
C c     input
C c*******************************************************************
C       ibound=-1
C       istart=-1
C       iread=1
C       iwrite=0
C       iprint=5
C       iseed=0
C       icheck=100
C       zi=(0.0,1.0)
C       tpi=2.0*acos(-1.0)
C c*******************************************************************
C c     end of input
C c*******************************************************************
C       open(unit=7,file='output',status='unknown')
C       open(unit=25,file='midout',status='unknown')
C       open(unit=98,file='control',status='unknown')
C       open(unit=36,file='remez2',status='unknown')
C       open(unit=37,file='remez4',status='unknown')
C       open(unit=38,file='remez2g',status='unknown')
C       open(unit=39,file='remez4g',status='unknown')
C       if(iread.eq.1) then
C       call sread
C       endif
C       read(25,*) dt,beta,am3,am,imass,iterl,iter2
C       close(25)
C c set a new seed by hand...
C       if(iseed.ne.0)then
C       seed=4139764973254.0
C       endif
C       write(7,*) 'seed: ', seed
C       call rranset(seed)
C       idum=-1
C       y=ran(idum)
C c     write(6,*) 'ran: ', y,idum
C c*******************************************************************
C c     initialization
C c     istart.lt.0 : start from tape
C c     istart=0    : ordered start
C c     istart=1    : random start
C c*******************************************************************
C       call init(istart)
C c  read in Remez coefficients
C       read(36,*) anum2(0)
C       read(37,*) anum4(0)
C       read(38,*) anum2g(0)
C       read(39,*) anum4g(0)
C       do i=1,ndiag
C       read(36,*) anum2(i),aden2(i)
C       read(37,*) anum4(i),aden4(i)
C       enddo
C       do i=1,ndiagg
C       read(38,*) anum2g(i),aden2g(i)
C       read(39,*) anum4g(i),aden4g(i)
C       enddo
C       read(36,*) bnum2(0)
C       read(37,*) bnum4(0)
C       read(38,*) bnum2g(0)
C       read(39,*) bnum4g(0)
C       do i=1,ndiag
C       read(36,*) bnum2(i),bden2(i)
C       read(37,*) bnum4(i),bden4(i)
C       enddo
C       do i=1,ndiagg
C       read(38,*) bnum2g(i),bden2g(i)
C       read(39,*) bnum4g(i),bden4g(i)
C       enddo
C c*******************************************************************
C c     print heading
C c*******************************************************************
C       traj=iterl*dt
C       proby=1.0/float(iterl)
C c     write(6, 9001)ksize,ksizet,kthird,Nf,dt,traj,ndiag,ndiagg,
C c    & iter2,beta,am3,am,imass
C       write(7, 9001)ksize,ksizet,kthird,Nf,dt,traj,ndiag,ndiagg,
C      & iter2,beta,am3,am,imass
C 9001  format(' ksize=',i3,' ksizet=',i3,/
C      1 ,' kthird=',i3,/
C      1 ,' Nf =',i3,/
C      1 ,' time step: dt=',f6.4,' trajectory length=',f9.6,/
C      1 ,' Remez ndiag: action =',i3' guidance=',i3,/
C      1 ,' # trajectories=',i6,' beta=',f9.6,/
C      1 ,' am3=',f6.4,' am=',f6.4/
C      1 ,' imass=',i2)
C c     write(6,9004) rescgg,rescga,respbp
C       write(7,9004) rescgg,rescga,respbp
C 9004  format(' Stopping residuals: guidance: ',e11.4,' acceptance: ',
C      &     e11.4,' estimator: ',e11.4)
C c     write(6,9044) rescgm
C       write(7,9044) rescgm
C 9044  format(' Stopping residuals: meson: ',e11.4)
C       call rranget(seed)
C c     write(6,*) 'seed: ', seed
C       write(7,*) 'seed: ', seed
C c*******************************************************************
C c       initialize for averages
C c*******************************************************************
C       actiona=0.0
C       vel2a=0.0
C       pbpa=0.0
C       ancg=0.0
C       ancgh=0.0
C       ancgf=0.0
C       ancgpf=0.0
C       ancgpv=0.0
C       ancghpv=0.0
C       ancgfpv=0.0
C       ancgpfpv=0.0
C       ancgma=0.0
C       yav=0.0
C       yyav=0.0 
C       naccp=0
C       ipbp=0
C       itot=0
C c*******************************************************************
C c     start of classical evolution
C c*******************************************************************
C       do 601 isweep=1,iter2
C c uncomment line below to go straight to measurement
C c     goto 666
C c*******************************************************************
C c     initialise trial fields
C c*******************************************************************
C       do 2007 mu=1,3
C       do 2007 i=1,kvol
C       thetat(i,mu)=theta(i,mu)
C 2007  continue
C c
C       call coef(ut,thetat)
C c*******************************************************************
C c  Pseudofermion fields: Phi = {MdaggerM(1)}^-1/4 * {MdaggerM(m)}^1/4 * R, where
C c   R is gaussian
C c*******************************************************************
C       do ia=1,Nf
C c
C       do idirac=1,4
C       do ithird=1,kthird
C       call gauss0(ps)
C       do i=1,kvol
C       R(ithird,i,idirac)=cmplx(ps(i,1),ps(i,2))
C       enddo
C       enddo
C       enddo
C c
C c  For now Phi = {MdaggerM}^0.25 * R
C c
C       call qmrherm(R,rescga,itercg,am,imass,anum4,aden4,ndiag,
C      &             0,isweep,0)
C       ancgpf=ancgpf+float(itercg)
C       do idirac=1,4
C       do i=1,kvol
C       do ithird=1,kthird
C       R(ithird,i,idirac)=X1(ithird,i,idirac)
C       enddo
C       enddo
C       enddo
C       call qmrherm(R,rescga,itercg,One,1,bnum4,bden4,ndiag,
C      &             0,isweep,0)
C       ancgpfpv=ancgpfpv+float(itercg)
C c
C       do idirac=1,4
C       do i=1,kvol
C       do ithird=1,kthird
C       Phi(ithird,i,idirac,ia)=X1(ithird,i,idirac)
C       enddo
C       enddo
C       enddo
C c
C       enddo
C c*******************************************************************
C c     heatbath for p 
C c*******************************************************************
C c  for some occult reason this write statement is needed to ensure compatibility with earlier versions
C c     write(6,*) idum
C c     write(98,*) idum
C       do mu=1,3
C       call gaussp(ps)
C       do i=1,kvol
C       pp(i,mu)=ps(i,1)
C       enddo
C       enddo
C c     write(6,*) idum
C c*******************************************************************
C c  call to Hamiltonian
C c     
C       call hamilton(Phi,
C      &      H0,hg,hp,S0,rescga,isweep,0,am,imass)
C       if(isweep.eq.1) then
C       action=S0/kvol
C       gaction=hg/kvol
C       paction=hp/kvol
C       endif 
C c     goto 501
C c*******************************************************************
C c      half-step forward for p
C c*******************************************************************
C       call force(Phi,rescgg,am,imass,isweep,0)
C       d=dt*0.5
C       do 2004 mu=1,3
C       do 2004 i=1,kvol
C       pp(i,mu)=pp(i,mu)-dSdpi(i,mu)*d
C 2004  continue
C c*******************************************************************
C c     start of main loop for classical time evolution
C c*******************************************************************
C       do 500 iter=1,itermax
C c
C c  step (i) st(t+dt)=st(t)+p(t+dt/2)*dt;
C c
C       d=dt
C       do 2001 mu=1,3
C       do 2001 i=1,kvol
C       thetat(i,mu)=thetat(i,mu)+d*pp(i,mu)
C 2001  continue
C c
C c  step (ii)  p(t+3dt/2)=p(t+dt/2)-dSds(t+dt)*dt (1/2 step on last iteration)
C c
C       call coef(ut,thetat)
C       call force(Phi,rescgg,am,imass,isweep,iter)
C c
C c test for end of random trajectory
C c 
C       ytest=ran(idum)
C       if(ytest.lt.proby)then
C       d=dt*0.5
C       do 2005 mu=1,3
C       do 2005 i=1,kvol
C       pp(i,mu)=pp(i,mu)-d*dSdpi(i,mu)
C 2005  continue
C       itot=itot+iter 
C       goto 501
C       else
C       d=dt
C       do 3005 mu=1,3
C       do 3005 i=1,kvol
C       pp(i,mu)=pp(i,mu)-d*dSdpi(i,mu)
C 3005  continue
C       endif
C c 
C 500   continue
C c**********************************************************************
C c  Monte Carlo step: accept new fields with probability=
C c              min(1,exp(H0-H1))
C c**********************************************************************
C 501   continue 
C       call hamilton(Phi,
C      &       H1,hg,hp,S1,rescga,isweep,-1,am,imass)
C       dH=H0-H1
C       dS=S0-S1
C       write(98,*) dH,dS
C       y=exp(dH)      
C       yav=yav+y 
C       yyav=yyav+y*y 
C       if(dH.lt.0.0)then
C       x=ran(idum)
C       if(x.gt.y)goto 600
C       endif
C c
C c     step accepted: set s=st
C c
C       do 2006 mu=1,3
C       do 2006 i=1,kvol
C       theta(i,mu)=thetat(i,mu)
C 2006  continue      
C       naccp=naccp+1
C       action=S1/kvol
C       gaction=hg/kvol
C       paction=hp/kvol
C 600   continue
C       write(11,*) isweep,gaction,paction
C       actiona=actiona+action 
C       vel2=0.0
C       do 457 mu=1,3
C       do 457 i=1,kvol
C       vel2=vel2+pp(i,mu)*pp(i,mu)
C 457   continue
C       vel2=vel2/(3*kvol)
C       vel2a=vel2a+vel2
C c
C c     uncomment to disable measurements
C c     goto 601
C 666   if((isweep/iprint)*iprint.eq.isweep)then
C       do 2066 mu=1,3
C       do 2066 i=1,kvol
C       thetat(i,mu)=theta(i,mu)
C 2066  continue      
C       call coef(ut,thetat)
C       call measure(pbp,respbp,ancgm,am,imass)
C c     call meson(rescgm,itercg,ancgm,am,imass)
C       pbpa=pbpa+pbp
C       ancgma=ancgma+ancgm
C       ipbp=ipbp+1
C c     write(11,*) pbp
C c     write(6,*) isweep,':  ',pbp,ancgm
C       endif
C c
C       if((isweep/icheck)*icheck.eq.isweep)then
C       call rranget(seed)
C       if(iwrite.eq.1) then
C       call swrite
C       endif
C       flush(100)
C       flush(200)
C c     flush(302)
C c     flush(400)
C c     flush(500)
C c     flush(501)
C       if(imass.ne.1)then
C c     flush(401)
C c     flush(402)
C c     flush(403)
C       endif
C c     write(7,9023) seed
C       endif
C c
C 601   continue
C c*******************************************************************
C c     end of main loop
C c*******************************************************************
C       actiona=actiona/iter2 
C       vel2a=vel2a/iter2 
C       pbpa=pbpa/ipbp
C       ancg=ancg/(Nf*itot)
C       ancgh=ancgh/(2*Nf*iter2)
C       ancgpf=ancgpf/(Nf*iter2)
C       ancgpv=ancgpv/(Nf*itot)
C       ancgf=ancgf/(Nf*itot)
C       ancgfpv=ancgfpv/(Nf*itot)
C       ancghpv=ancghpv/(2*Nf*iter2)
C       ancgpfpv=ancgpfpv/(iter2*Nf)
C       ancgma=ancgma/ipbp
C       yav=yav/iter2
C       yyav=yyav/iter2-yav*yav 
C       yyav=sqrt(yyav/(iter2-1)) 
C       atraj=dt*itot/iter2 
C c*******************************************************************
C c     print global averages
C c*******************************************************************
C c     write(6, 9022) iter2,naccp,atraj,yav,yyav,ancg,ancgpv,ancgh,ancghpv,ancgf,
C c    & ancgfpv,ancgpf,ancgpfpv,pbpa,vel2a,actiona
C       write(7, 9022) iter2,naccp,atraj,yav,yyav,
C      & ancg,ancgpv,ancgh,ancghpv,ancgf,ancgfpv,ancgpf,ancgpfpv,
C      & pbpa,ancgma,vel2a,actiona
C 9022  format(' averages for last ',i6,' trajectories',/ 
C      & 1x,' # of acceptances: ',i6,' average trajectory length= ',f8.3/
C      & 1x,' <exp-dH>=',e11.4,' +/-',e10.3/
C      2 1x,' av. # QMR itr.'/
C      & 1x,'     guidance: DWF  ',f9.3,'; PV  ',f9.3/
C      & 1x,'   acceptance: DWF  ',f9.3,'; PV  ',f9.3/
C      & 1x,'        force: DWF  ',f9.3,'; PV  ',f9.3/
C      & 1x,'pseudofermion: DWF  ',f9.3,'; PV  ',f9.3/
C      1 1x,' psibarpsi=',e11.3/
C      2 1x,' av. # QMR itr.',f9.3//
C      & 1x,' mean square velocity=',e10.3,'; action per site=',e10.3//)
C       write(7, 9024)
C       write(7, 9024)
C 9024  format(1x)
C c
C       close(11)
C c
C       if(iwrite.eq.1) then
C       call rranget(seed)
C       call swrite
C       write(7,*) 'seed: ', seed
C       endif
C c
C       stop
C       end
C c******************************************************************
C c   calculate dSds for gauge fields at each intermediate time
C c******************************************************************
C       subroutine force(Phi,res1,am,imass,isweep,iter)
C       parameter(ksize=12,ksizet=12,kthird=24,kvol=ksize*ksize*ksizet)
C       parameter(kvol2=ksize*ksize)
C       parameter(kferm=4*kvol*kthird)
C       parameter(ndiagg=12,ndiag=ndiagg)
C       parameter(One=1.0)
C       parameter(Nf=1)
C       common/remez2g/anum2(0:ndiag),aden2(ndiag),
C      &              bnum2(0:ndiag),bden2(ndiag)
C       common/remez4g/anum4(0:ndiag),aden4(ndiag),
C      &              bnum4(0:ndiag),bden4(ndiag)
C       common/phi0/Phi0(kferm,25)
C       common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
C       common/para/beta,am3,ibound
C       common/param/ancg,ancgh,ancgf,ancgpf
C       common/parampv/ancgpv,ancghpv,ancgfpv,ancgpfpv
C       common/gforce/dSdpi(kvol,3)
C       common/vector/X1(kferm)
C c     complex Phi(kferm,Nf),X2(kferm)
C c     complex X1,u,Phi0
C       complex*16 Phi(kferm,Nf),X2(kferm)
C       complex*16 X1,u,Phi0
C       real*8 anum2,aden2,bnum2,bden2
C       real*8 anum4,aden4,bnum4,bden4
C c
C c     write(6,111)
C 111   format(' Hi from force')
C c
C       do mu=1,3
C       do i=1,kvol
C          dSdpi(i,mu)=0.0
C       enddo
C       enddo
C c
C c uncomment this line to quench the fermions!
C c     return
C c pseudofermion action is
C c   Phi^dagger {MdaggerM(1)}^1/4 {MdaggerM(m)})^-1/2 {MdaggerM(1)}^1/4 Phi
C c
C       do ia=1,Nf
C c
C       do i=1,kferm
C       X2(i)=Phi(i,ia)
C       enddo
C       call qmrherm(X2,res1,itercg,One,1,anum4,aden4,ndiag,
C      &             1,isweep,iter)
C       ancgpv=ancgpv+float(itercg)
C       do i=1,kferm
C       X2(i)=X1(i)
C       enddo
C c
C       call qmrherm(X2,res1,itercg,am,imass,bnum2,bden2,ndiag,
C      &             0,isweep,iter)
C       ancg=ancg+float(itercg)
C c     write(111,*) itercg
C       do i=1,kferm
C       X2(i)=X1(i)
C       enddo
C c
C c  evaluates -X2dagger * d/dpi[{MdaggerM(m)}^1/2] * X2
C       call qmrherm(X2,res1,itercg,am,imass,anum2,aden2,ndiag,
C      &             2,isweep,iter)
C       ancgf=ancgf+float(itercg)
C c     write(113,*) itercg
C c  evaluates +2Re{Phidagger * d/dpi[{MdaggerM(1)}^1/4] * X2}
C       call qmrherm(X2,res1,itercg,One,1,anum4,aden4,ndiag,
C      &             3,isweep,iter)
C       ancgfpv=ancgfpv+float(itercg)
C c
C       enddo
C c
C       if(ibound.eq.-1)then
C       ioffset=(ksizet-1)*kvol2
C       do i=ioffset+1,kvol
C       dSdpi(i,3)=-dSdpi(i,3)
C       enddo
C       endif
C c
C       b=beta*Nf
C       do mu=1,3
C       do i=1,kvol
C          dSdpi(i,mu)=dSdpi(i,mu)+b*theta(i,mu)
C       enddo
C       enddo
C c
C       return
C       end
C c******************************************************************
C c   Evaluation of Hamiltonian function
C c******************************************************************
C       subroutine hamilton(Phi,
C      &       h,hg,hp,s,res2,isweep,iflag,am,imass)
C       parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
C       parameter(kmom=3*kvol,kferm=4*kvol*kthird)
C       parameter(ndiag=25)
C       parameter(One=1.0)
C       parameter(Nf=1)
C       common/trial/u(kvol,3),theta(kmom),pp(kmom)
C       common/remez2/anum2(0:ndiag),aden2(ndiag),
C      &              bnum2(0:ndiag),bden2(ndiag)
C       common/remez4/anum4(0:ndiag),aden4(ndiag),
C      &              bnum4(0:ndiag),bden4(ndiag)
C       common/param/ancg,ancgh,ancgf,ancgpf
C       common/parampv/ancgpv,ancghpv,ancgfpv,ancgpfpv
C       common /para/beta,am3,ibound
C       common/vector/ X1(kferm)      
C       common/dum1/R(kferm),ps(kvol,2)
C c     complex Phi(kferm,Nf)
C c     complex X1,R
C c     complex u
C       complex*16 Phi(kferm,Nf)
C       complex*16 X1,R
C       complex*16 u
C       real*8 hp,hg,hf,h,s
C       real*8 anum2,aden2,bnum2,bden2
C       real*8 anum4,aden4,bnum4,bden4
C c     write(6,111)
C 111   format(' Hi from hamilton')
C c
C       hp=0.0
C       hg=0.0
C       hf=0.0
C c
C       do 22 imom=1,kmom
C       hp=hp+pp(imom)*pp(imom)
C 22    continue
C c
C       hp=hp*0.5
C c
C       do imom=1,kmom
C       hg=hg+theta(imom)*theta(imom)
C       enddo
C       hg=0.5*Nf*beta*hg
C       h=hg+hp
C c 
C c uncomment these lines to quench the fermions!
C c     write(6,*) isweep,':  hg', hg,'   hp', hp,'   h',h
C c     return
C c         
C c  pseudofermion action is
C c   Phi^dagger {MdaggerM(1)}^1/4 {MdaggerM(m)})^-1/2 {MdaggerM(1)}^1/4 Phi
C c
C       do ia=1,Nf
C c
C       do i=1,kferm
C       R(i)=Phi(i,ia)
C       enddo
C       call qmrherm(R,res2,itercg,One,1,anum4,aden4,ndiag,
C      &             0,isweep,iflag)
C       ancghpv=ancghpv+float(itercg)
C       do i=1,kferm
C       R(i)=X1(i)
C       enddo
C c
C       call qmrherm(R,res2,itercg,am,imass,bnum2,bden2,ndiag,
C      &             0,isweep,iflag)
C       ancgh=ancgh+float(itercg)
C c
C       do 4 iferm=1,kferm
C       hf=hf+conjg(R(iferm))*X1(iferm)
C 4     continue
C c
C       enddo
C c
C       h=hg+hp+hf
C c     write(6,*) isweep,':  hg', hg,'   hp', hp,'   hf', hf,
C c    &   '   h',h
C       s=hg+hf
C c
C       return
C       end              
C c******************************************************************
C c    multisolver matrix inversion via Lanczos technique
C c  eg. Golub & van Loan "Matrix Computations" 9.3.1
C c       solves (MdaggerM+diag)*x=Phi for ndiag different values of diag
C c   iflag=0: simply evaluates X = {MdaggerM}^p * Phi
C c   can be interchanged with congrad for p=-1
C c   iflag=1: in addition updates Phi0 register needed for PV force term
C c   iflag=2: evaluates DWF force term
C c   iflag=3: evaluates PV force term
C c*****************************************************************m
C       subroutine qmrherm(Phi,res,itercg,am,imass,anum,aden,ndiag,
C      &                   iflag,isweep,iter)
C       parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
C       parameter(kferm=4*kthird*kvol)
C c     parameter(niterc=10*kferm)
C       parameter(niterc=7500)
C       common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
C       common/para/bbb,am3,ibound
C       common/vector/x(kferm)
C       common/gforce/dSdpi(kvol,3)
C       common/phi0/Phi0(kferm,25)
C c     complex Phi(kferm)
C c     complex x,u,Phi0
C       complex*16 Phi(kferm)
C       complex*16 x,u,Phi0
C       real*8 anum(0:ndiag),aden(ndiag),coeff
C c     complex vtild(kferm),q(kferm),pm1(kferm,ndiag)
C c     complex qm1(kferm),p(kferm,ndiag),x3(kferm),R(kferm)
C c     complex x1(kferm,ndiag),x2(kferm)
C c     real alpha(ndiag),beta
C c     real amu(ndiag),d(ndiag),dm1(ndiag),rho(ndiag),rhom1(ndiag)
C       complex*16 vtild(kferm),q(kferm),pm1(kferm,ndiag)
C       complex*16 qm1(kferm),p(kferm,ndiag),x3(kferm),R(kferm)
C       complex*16 x1(kferm,ndiag),x2(kferm)
C       real*8 alpha(ndiag),beta
C       real*8 amu(ndiag),d(ndiag),dm1(ndiag),rho(ndiag),rhom1(ndiag)
C c
C c     write(6,111)
C 111   format(' Hi from qmrherm')
C c
C       resid=sqrt(kferm*res*res)
C c     write(6,*) iflag, resid
C       itercg=0
C c
C c   initialise r=Phi
C c
C       do i=1,kferm
C          R(i)=Phi(i)
C          qm1(i)=(0.0,0.0)
C          x(i)=anum(0)*Phi(i)
C       enddo
C       beta=0.0
C       do i=1,kferm
C          beta=beta+conjg(r(i))*r(i)
C       enddo
C       beta=sqrt(beta)
C       phimod=beta
C c     write(6,*) '|| Phi || = ', phimod
C c
C       do niter=1,niterc
C       itercg=itercg+1
C c
C c  Lanczos steps
C c
C       do i=1,kferm
C          q(i)=R(i)/beta
C       enddo
C c
C       call dslash(vtild,q,u,am,imass)
C       call dslashd(x3,vtild,u,am,imass)
C c
C       alphatild=0.0
C       do i=1,kferm
C          alphatild=alphatild+conjg(q(i))*x3(i)
C       enddo
C c
C       do i=1,kferm
C          R(i)=x3(i)-alphatild*q(i)-beta*qm1(i)
C          qm1(i)=q(i)
C       enddo
C c
C       beta0=beta
C       beta=0.0
C       do i=1,kferm
C          beta=beta+conjg(R(i))*R(i)
C       enddo
C       beta=sqrt(beta)
C c
C       do idiag=1,ndiag
C       alpha(idiag)=alphatild+aden(idiag)
C       enddo
C c
C       if(niter.eq.1)then
C            do idiag=1,ndiag
C              d(idiag)=alpha(idiag)
C              do i=1,kferm
C               p(i,idiag)=q(i)
C               pm1(i,idiag)=p(i,idiag)
C              enddo
C              rho(idiag)=beta0/alpha(idiag)
C              rhom1(idiag)=rho(idiag)
C              do i=1,kferm
C               x1(i,idiag)=rho(idiag)*q(i)
C              enddo
C            enddo
C       else
C            do idiag=1,ndiag
C              amu(idiag)=beta0/d(idiag)
C              dm1(idiag)=d(idiag)
C              d(idiag)=alpha(idiag)-beta0*amu(idiag)
C              do i=1,kferm
C               p(i,idiag)=q(i)-amu(idiag)*pm1(i,idiag)
C               pm1(i,idiag)=p(i,idiag)
C              enddo
C              rho(idiag)=-amu(idiag)*dm1(idiag)*rhom1(idiag)/d(idiag)
C c  Convergence criterion (a bit ad hoc for now...)
C              if(idiag.eq.1)then
C                rhomax=abs(phimod*rho(idiag))
C              else
C                if(abs(phimod*rho(idiag)).gt.rhomax)
C      &             rhomax=abs(phimod*rho(idiag))
C              endif
C              rhom1(idiag)=rho(idiag)
C              do i=1,kferm
C                x1(i,idiag)=x1(i,idiag)+rho(idiag)*p(i,idiag)
C              enddo
C            enddo
C c  check to see whether the residual is acceptable for all ndiag....
C c  criterion is a bit ad hoc -- relaxing by a factor arelax improves code
C c  stability and leads to quicker convergence
C            arelax=2.0
C            if(rhomax.lt.arelax*resid) then
C c          if(rhomax.lt.resid) then
C c          call testinv(Phi,resmax,itercg,am,imass,x1,aden,ndiag)
C c  convergence based on || residual || not working well in single precision...
C c          if(resmax.lt.resid) goto 8
C              goto 8
C            endif
C       endif
C c
C c   end of loop over iter
C       enddo
C       write(7,*) 'QMRniterc!, isweep,iter,iflag,imass,anum,ndiag = '
C      &        ,isweep,iter,iflag,imass,anum(0),ndiag
C 8     continue
C c
C       if(iflag.lt.2)then
C c  Now evaluate solution x=(MdaggerM)^p * Phi
C       do idiag=1,ndiag
C       do i=1,kferm
C       x(i)=x(i)+anum(idiag)*x1(i,idiag)
C       enddo
C       enddo
C c
C c  update phi0 block if required...
C       if(iflag.eq.1) then
C       do idiag=1,ndiag
C       do i=1,kferm
C       Phi0(i,idiag)=x1(i,idiag)
C       enddo
C       enddo
C       endif
C c
C       else
C c
C       do idiag=1,ndiag
C c
C c  X2 = M*X1
C       do i=1,kferm
C       R(i)=x1(i,idiag)
C       enddo
C       call dslash(X2,R,u,am,imass)
C c
C       if(iflag.eq.2)then
C           coeff=anum(idiag)
C           call derivs(R,X2,coeff,0)
C       else
C           coeff=-anum(idiag)
C           do i=1,kferm
C           R(i)=Phi0(i,idiag)
C           enddo
C           call derivs(R,X2,coeff,0)
C c
C           call dslash(X2,R,u,am,imass)
C           do i=1,kferm
C           R(i)=x1(i,idiag)
C           enddo
C           call derivs(X2,R,coeff,1)
C       endif
C c
C       enddo
C       endif
C c
C c
C       return
C       end
c**********************************************************************
c  iflag = 0 : evaluates Rdagger*(Mdagger)'*X2
c  iflag = 1 : evaluates Rdagger*(M)'*X2
c**********************************************************************
      subroutine derivs(R,X2,anum,iflag)
      use purefunctions
      implicit none
c      complex, intent(in) :: R(kthird, ksize, ksize, ksizet, 4)
c      complex, intent(in) :: X2(kthird, ksize, ksize, ksizet, 4)
      integer, parameter :: ksize=12, ksizet=12, kthird=24
      real, parameter :: akappa = 0.5

      complex*16, intent(in) :: R(kthird, 0:ksize+1, 0:ksize+1, 
     &                            0:ksizet+1, 4)
      complex*16, intent(in) :: X2(kthird, 0:ksize+1, 0:ksize+1, 
     &                            0:ksizet+1, 4)
      real*8, intent(in) :: anum
      integer, intent(in) :: iflag

      
      common/dirac/gamval(6,4),gamin(6,4)
      common/gforce/dSdpi(ksize,ksize,ksizet,3)

c      complex :: gamval
      complex*16 :: gamval
      integer :: gamin
      real :: dSdpi

c      complex*16 :: tzi
      real*8 :: tzi_real
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
      end
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
C       complex*16 vtild(kferm)
C       complex*16 x3(kferm)
C       complex*16 x(kferm,ndiag),x1(kferm),x2(kferm)
C       complex*16 u
C c     complex vtild(kferm)
C c     complex x3(kferm)
C c     complex x(kferm,ndiag),x1(kferm),x2(kferm)
C c      complex u
C       real*8 residual
C       real*8 aden(ndiag)
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
C c******************************************************************
C c    matrix inversion via conjugate gradient algorithm
C c       solves (Mdagger)Mx=Phi, 
C c           NB. no even/odd partitioning
C c******************************************************************
C       subroutine congrad(Phi,res,itercg,am,imass)
C       parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
C       parameter(kferm=4*kthird*kvol)
C       parameter(niterc=kthird*kvol)
C       common/trial/u(kvol,3),thetas(kvol,3),pp(kvol,3)
C       common/para/bbb,am3,ibound
C       common/vector/x(kferm)
C c     complex Phi(kferm)
C c     complex x,u
C c     complex x1(kferm),x2(kferm),p(kferm),r(kferm)
C       complex*16 Phi(kferm)
C       complex*16 x,u
C       complex*16 x1(kferm),x2(kferm),p(kferm),r(kferm)
C c     write(6,111)
C 111   format(' Hi from congrad')
C c
C       resid=kferm*res*res
C       itercg=0
C c
C       do 1 nx=1,niterc
C       itercg=itercg+1
C       if(nx.gt.1) goto 51
C c
C c   initialise p=x, r=Phi(na)
C c
C       do 2 i=1,kferm
C       p(i)=x(i)
C       r(i)=Phi(i)
C 2     continue
C       betad=1.0
C       alpha=1.0
C 51    alphad=0.0
C c
C c  x1=Mp
C c
C       call dslash(x1,p,u,am,imass)
C c
C       if(nx.eq.1) goto 201
C c
C c   alpha=(r,r)/(p,(Mdagger)Mp)
C c 
C       alphad=0.0
C       do 31 i=1,kferm
C       alphad=alphad+conjg(x1(i))*x1(i)
C 31    continue
C       alpha=alphan/alphad
C c      
C c   x=x+alpha*p
C c
C       do 4 i=1,kferm
C       x(i)=x(i)+alpha*p(i)
C 4     continue
C 201   continue
C c     
C c   x2=(Mdagger)x1=(Mdagger)Mp
C c
C       call dslashd(x2,x1,u,am,imass)
C c
C c   r=r-alpha*(Mdagger)Mp
C c
C       do 6 i=1,kferm
C       r(i)=r(i)-alpha*x2(i)
C 6     continue
C c
C c   beta=(r_k+1,r_k+1)/(r_k,r_k)
C c
C       betan=0.0 
C       do 61 i=1,kferm
C       betan=betan+conjg(r(i))*r(i) 
C 61    continue 
C       beta=betan/betad
C       betad=betan
C       alphan=betan
C c
C       if(nx.eq.1) beta=0.0
C c
C c   p=r+beta*p
C c
C       do 7 i=1,kferm
C       p(i)=r(i)+beta*p(i)
C 7     continue
C       if(betan.lt.resid) goto 8
C 1     continue
C c     write(6,1000)
C       write(7,1000)
C 1000  format(' # iterations of congrad exceeds niterc')
C 8     continue
C c     write(6,*) itercg
C       return
C       end      
C c*****************************************************************
C c   Calculate fermion expectation values via a noisy estimator
C c   -matrix inversion via conjugate gradient algorithm
C c       solves Mx=x1
C c     (Numerical Recipes section 2.10 pp.70-73)   
C c*******************************************************************
C       subroutine measure(psibarpsi,res,aviter,am,imass)
C       parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
C       parameter(akappa=0.5)
C       parameter(knoise=10)
C       common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
C       common/para/beta,am3,ibound
C       common/dirac/gamval(6,4),gamin(6,4)
C       common /neighb/id(kvol,3),iu(kvol,3)
C       common/vector/xi(kthird,kvol,4)
C       common/ranseed/idum
C       common/v/v(97)
C c     complex x(kvol,4), Phi(kthird,kvol,4)
C c     complex xi,gamval
C c     complex psibarpsi1,psibarpsi2
C c     complex u
C       complex*16 x(kvol,4), Phi(kthird,kvol,4)
C       complex*16 xi,gamval
C       complex*16 psibarpsi1,psibarpsi2
C       complex*16 u
C       real*8 cnum(0:1),cden(1)
C       real ps(kvol,2),pt(kvol,2)
C       real*8 pbp(knoise)
C       integer gamin
C c     write(6,*) 'hi from measure'
C c
C       iter=0
C c     pbp=0.0
C       cnum(0)=0.0
C       cnum(1)=1.0
C       cden(1)=0.0
C c
C       do inoise=1,knoise
C c
C c     set up noise
C       call gauss0(ps)
C       psibarpsi1=(0.0,0.0)
C       call gauss0(pt)
C       psibarpsi2=(0.0,0.0)
C c
C       do idsource=1,2
C c
C c  source on domain wall at ithird=1
C c   
C       if(imass.ne.5)then
C          do idirac=1,4
C            if(idirac.eq.idsource)then
C              do i=1,kvol
C                x(i,idirac)=cmplx(ps(i,1),ps(i,2))
C              enddo
C            else
C              do i=1,kvol
C                x(i,idirac)=(0.0,0.0)
C              enddo
C            endif
C          enddo
C       else
C          do idirac=1,4
C            if(idirac.eq.idsource+2)then
C              do i=1,kvol
C                x(i,idirac)=cmplx(ps(i,1),ps(i,2))
C              enddo
C            else
C              do i=1,kvol
C                x(i,idirac)=(0.0,0.0)
C              enddo
C            endif
C          enddo
C       endif
C c
C c
C       if(imass.ne.5)then
C          do ithird=1,kthird
C          if(ithird.eq.1)then
C             do idirac=1,4
C             do i=1,kvol
C                xi(ithird,i,idirac)=x(i,idirac)
C             enddo
C             enddo
C          else
C             do idirac=1,4
C             do i=1,kvol
C                xi(ithird,i,idirac)=(0.0,0.0)
C             enddo
C             enddo
C          endif
C          enddo
C       else
C c   xi = 0.5(1+gamma_4)*gamma_5*eta on DW at ithird=1
C          do ithird=1,kthird
C          if(ithird.eq.1)then
C             do idirac=1,2
C             do i=1,kvol
C                xi(ithird,i,idirac)=-x(i,idirac+2)
C                xi(ithird,i,idirac+2)=(0.0,0.0)
C             enddo
C             enddo
C          else
C c   xi = 0.5(1+gamma_4)*eta on DW at ithird=1
C             do idirac=1,4
C             do i=1,kvol
C                xi(ithird,i,idirac)=(0.0,0.0)
C             enddo
C             enddo
C          endif
C          enddo
C       endif
C c
C c Phi= Mdagger*xi
C c
C       call dslashd(Phi,xi,u,am,imass)
C c     call qmrherm(Phi,res,itercg,am,imass,cnum,cden,1,0)
C       call congrad(Phi,res,itercg,am,imass)
C       iter=iter+itercg
C c
C       if(imass.ne.5)then
C c  pbp1 = x^dagger (0.5(1+gamma_4)) xi(kthird)
C          do i=1,kvol
C          psibarpsi1=psibarpsi1
C      &    +conjg(x(i,idsource))*xi(kthird,i,idsource)
C          enddo
C       else
C c  pbp1 = x^dagger (0.5(1-gamma_4)) xi(1)
C          do i=1,kvol
C          psibarpsi1=psibarpsi1
C      &    +conjg(x(i,idsource+2))*xi(1,i,idsource+2)
C          enddo
C       endif
C c
C c
C c  source on domain wall at ithird=kthird
C       idsource2=idsource+2
C c
C       if(imass.ne.5)then
C          do idirac=1,4
C          do i=1,kvol
C          if(idirac.eq.idsource2)then
C              x(i,idirac)=cmplx(pt(i,1),pt(i,2))
C          else
C              x(i,idirac)=(0.0,0.0)
C          endif
C          enddo
C          enddo
C       else
C          do idirac=1,4
C          do i=1,kvol
C          if(idirac.eq.idsource2-2)then
C              x(i,idirac)=cmplx(pt(i,1),pt(i,2))
C          else
C              x(i,idirac)=(0.0,0.0)
C          endif
C          enddo
C          enddo
C       endif
C c
C       if(imass.ne.5)then
C c   xi = 0.5(1-gamma_4)*eta on DW at ithird=kthird
C          do ithird=1,kthird
C          if(ithird.eq.kthird)then
C             do idirac=1,4
C             do i=1,kvol
C                xi(ithird,i,idirac)=x(i,idirac)
C             enddo
C             enddo
C          else
C             do idirac=1,4
C             do i=1,kvol
C                xi(ithird,i,idirac)=(0.0,0.0)
C             enddo
C             enddo
C          endif
C          enddo
C       else
C c   xi = 0.5(1-gamma_4)*gamma_5*eta on DW at ithird=kthird
C          do ithird=1,kthird
C          if(ithird.eq.kthird)then
C             do idirac=1,2
C             do i=1,kvol
C                xi(ithird,i,idirac+2)=-x(i,idirac)
C                xi(ithird,i,idirac)=(0.0,0.0)
C             enddo
C             enddo
C          else
C             do idirac=1,4
C             do i=1,kvol
C                xi(ithird,i,idirac)=(0.0,0.0)
C             enddo
C             enddo
C          endif
C          enddo
C       endif
C c
C c Phi= Mdagger*xi
C c
C       call dslashd(Phi,xi,u,am,imass)
C c
C c xi= (M)**-1 * Phi
C c
C c     call qmrherm(Phi,res,itercg,am,imass,cnum,cden,1,0)
C       call congrad(Phi,res,itercg,am,imass)
C       iter=iter+itercg
C c
C       if(imass.ne.5)then
C c pbp2= - x^dagger (0.5(1-gamma_4)) xi(1)
C          do i=1,kvol
C          psibarpsi2=psibarpsi2
C      &     +conjg(x(i,idsource2))*xi(1,i,idsource2)
C          enddo
C       else
C c pbp2= - x^dagger (0.5(1-gamma_4)) xi(kthird)
C          do i=1,kvol
C          psibarpsi2=psibarpsi2
C      &     +conjg(x(i,idsource))*xi(kthird,i,idsource)
C          enddo
C       endif
C c
C c  end trace on Dirac indices....
C       enddo
C c
C       if(imass.eq.1)then
C          psibarpsi1=psibarpsi1/kvol
C          psibarpsi2=psibarpsi2/kvol
C          pbp(inoise)=psibarpsi1+psibarpsi2
C       elseif(imass.eq.3)then
C          psibarpsi1=(0.0,-1.0)*psibarpsi1/kvol
C          psibarpsi2=(0.0,+1.0)*psibarpsi2/kvol
C          pbp(inoise)=psibarpsi1+psibarpsi2
C       elseif(imass.eq.5)then
C          psibarpsi1=(0.0,-1.0)*psibarpsi1/kvol
C          psibarpsi2=(0.0,-1.0)*psibarpsi2/kvol
C          pbp(inoise)=psibarpsi1+psibarpsi2
C       endif
C c        write(6,*) real(psibarpsi1),aimag(psibarpsi1),
C c    &       real(psibarpsi2),aimag(psibarpsi2)
C          write(100,*) real(psibarpsi1),aimag(psibarpsi1),
C      &       real(psibarpsi2),aimag(psibarpsi2)
C c
C c end loop on noise
C       enddo
C c
C       psibarpsi=0.0
C       susclsing=0.0
C c
C       do inoise=1,knoise
C       psibarpsi=psibarpsi+pbp(inoise)
C          do jnoise=knoise,inoise+1,-1
C          susclsing=susclsing+pbp(inoise)*pbp(jnoise)
C          enddo
C       enddo
C       psibarpsi=psibarpsi/knoise
C       susclsing=2*kvol*susclsing/(knoise*(knoise-1))
C          write(200,*) psibarpsi,susclsing
C       aviter=float(iter)/(4*knoise)
C       return
C       end
C c******************************************************************
C c   Calculate meson correlators using point sources on domain walls
C c   -matrix inversion via conjugate gradient algorithm
C c       solves Mx=x1
C c     (Numerical Recipes section 2.10 pp.70-73)   
C c*******************************************************************
C       subroutine meson(res,itercg,aviter,am,imass)
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
C       complex*16 x(kvol,4),x0(kvol,4),Phi(kthird,kvol,4)
C       complex*16 xi,gamval
C       complex*16 prop00(kvol,3:4,1:2),prop0L(kvol,3:4,3:4)
C c     complex prop00n(kvol,3:4,1:2),prop0Ln(kvol,3:4,3:4)
C       real cpm(0:ksizet-1),cmm(0:ksizet-1)
C c     complex cpmn(0:ksizet-1),cmmn(0:ksizet-1)
C c     complex cferm1(0:ksizet-1), cferm2(0:ksizet-1)
C c     complex u
C       complex*16 cferm1(0:ksizet-1), cferm2(0:ksizet-1)
C       complex*16 u
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
C       ixxx=int(ksize*ran(idum))+1
C       iyyy=int(ksize*ran(idum))+1
C       ittt=int(ksizet*ran(idum))+1
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
      implicit none
      integer, parameter :: ksize=12, ksizet=12
      common/gauge/ theta(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3), seed
      real :: theta
      real*8 :: seed
      open(unit=10,file='con',
     1     status='unknown',form='unformatted')
      read (10) theta(1:ksize, 1:ksize, 1:ksizet, :), seed
      close(10)
      call update_halo_4_real(3, theta)
      return
      end
c
      subroutine swrite
      implicit none
      integer, parameter :: ksize=12, ksizet=12
      common/gauge/ theta(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3), seed
      real :: theta
      real*8 seed
      open(unit=31,file='con',
     1     status='unknown',form='unformatted')
      write (31) theta(1:ksize, 1:ksize, 1:ksizet, :), seed
      close(31)
      return
      end
c
      subroutine init(nc)
c*******************************************************************
c     sets initial values
c     nc=0 cold start
c     nc=1 hot start
c     nc<0 no initialization
c*******************************************************************
      implicit none
      integer, parameter :: ksize=12, ksizet=12, kthird=24
      real, parameter :: akappa=0.5
      integer, intent(in) :: nc
      common/gauge/theta(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3), seed
      common/dirac/gamval(6,4),gamin(6,4)
      common/ranseed/idum
c     complex gamval,one,zi
      complex*16 :: gamval,one,zi
      real :: ran
      real :: theta
      real*8 :: seed
      integer :: gamin
      integer :: idum
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
c             theta(ix, iy, it, mu) = 2.0 * g * rranf() - 1.0
              theta(ix, iy, it, mu) = 2.0 * g * ran(idum) - 1.0
            enddo
          enddo
        enddo
      enddo
      call update_halo_4(3, theta)
      return
      end
c******************************************************************
c   calculate compact links from non-compact links
c******************************************************************
      pure subroutine coef(u,theta)
      use purefunctions
      implicit none
      integer, parameter :: ksize=12, ksizet=12
      common/para/beta,am3,ibound
c     complex u(kvol,3)
      complex*16, intent(inout) :: u(0:ksize+1, 0:ksize+1, 
     &                               0:ksizet+1, 3)
      real, intent(in) :: theta(ksize, ksize, ksizet, 3)
      real :: beta, am3
      integer :: ibound
      integer :: ix, iy, it, mu
c
      do mu=1, 3
        do it = 1, ksizet
          do iy = 1, ksize
            do ix = 1, ksize
c             u(ix, iy, it, mu) = exp(cmplx(0.0, theta(ix, iy, it, mu)))
              u(ix, iy, it, mu) = (1.0 + 
     &                             cmplx(0.0, theta(ix, iy, it, mu)))
            enddo
          enddo
        enddo
      enddo
c
c  anti-p.b.c. in timelike direction
      if(ibound.eq.-1)then
        u(:, :, ksizet, 3) = -u(:, :, ksizet, 3)
      end if
c      
      call update_halo_4(3, u)
      call update_halo_4_real(3, theta)
      return
      end
c**********************************************************************
c calculate vector of gaussian random numbers with unit variance
c to refresh momenta
c   Numerical Recipes pp.203
c**********************************************************************
      pure subroutine gaussp(ps)
      use purefunctions
      implicit none
      integer, parameter :: ksize=12, ksizet=12
      common/trans/tpi 
      common/ranseed/idum
      real :: ran
      real, intent(out) :: ps(0:ksize+1, 0:ksize+1, 0:ksizet+1, 2)
      real :: tpi
      integer ix, iy, it
      integer :: idum
      real :: theta
c     write(6,1)
1     format(' Hi from gaussp')
      do it = 1, ksizet
        do iy = 1, ksize
          do ix = 1, ksize
            ps(ix, iy, it, 2) = sqrt(-2.0 * log(ran(idum)))
          end do
        end do
      end do
      do it = 1, ksizet
        do iy = 1, ksize
          do ix = 1, ksize
            theta = tpi * ran(idum)
            ps(ix, iy, it, 1) = ps(ix, iy, it, 2) * sin(theta)
            ps(ix, iy, it, 2) = ps(ix, iy, it, 2) * cos(theta)
          end do
        end do
      end do
      call update_halo_4_real(2, ps)

      return
      end      
c**********************************************************************
c calculate vector of gaussian random numbers with unit variance
c to generate pseudofermion fields R
c   Numerical Recipes pp.203
c**********************************************************************
      pure subroutine gauss0(ps)
      use purefunctions
      implicit none
      integer, parameter :: ksize=12, ksizet=12
      common/trans/tpi 
      common/ranseed/idum
      real :: ran
      real :: tpi
      integer :: idum
      real, intent(out) :: ps(0:ksize+1, 0:ksize+1, 0:ksizet+1, 2)
      integer :: ix, iy, it
      real :: theta
c     write(6,1)
1     format(' Hi from gauss0')
      do it = 1, ksizet
        do iy = 1, ksize
          do ix = 1, ksize
            ps(ix, iy, it, 2) = sqrt(-log(ran(idum)))
          end do
        end do
      end do
      do it = 1, ksizet
        do iy = 1, ksize
          do ix = 1, ksize
            theta = tpi * ran(idum)
            ps(ix, iy, it, 1) = ps(ix, iy, it, 2) * sin(theta)
            ps(ix, iy, it, 2) = ps(ix, iy, it, 2) * cos(theta)
          end do
        end do
      end do
      call update_halo_4_real(2, ps)
1001  continue 
      return
      end      

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
      use purefunctions
      implicit none
      integer, parameter :: ksize=12,ksizet=12,kthird=24
      real, parameter :: akappa=0.5
      common/para/beta,am3,ibound
      common/dirac/gamval(6,4),gamin(6,4)
c     complex, intent(in) :: u(0:ksize+1,0:ksize+1,0:ksizet+1,3)
c     complex, intent(in) :: Phi(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
c     complex, intent(in) :: R(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
c     complex gamval
c     complex zkappa
      complex*16, intent(in) :: u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3)
      complex*16, intent(inout) :: Phi(kthird, 0:ksize+1,
     &                                 0:ksize+1, 0:ksizet+1, 4)
      complex*16, intent(in) :: R(kthird, 0:ksize+1, 0:ksize+1,
     &                            0:ksizet+1, 4)
      integer, intent(in) :: imass
      real, intent(in) :: am
      complex*16 :: gamval
      complex*16 :: zkappa
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
      end
c***********************************************************************
      pure subroutine dslashd(Phi,R,u,am,imass)
c
c     calculates Phi = Mdagger*R
c
      use purefunctions
      implicit none
      integer, parameter :: ksize=12,ksizet=12,kthird=24
      real, parameter :: akappa=0.5
      common/para/beta,am3,ibound
      common/dirac/gamval(6,4),gamin(6,4)
c     complex, intent(in) ::  u(0:ksize+1,0:ksize+1,0:ksizet+1,3)
c     complex, intent(inout) :: Phi(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
c     complex, intent(in) R(kthird,0:ksize+1,0:ksize+1,0:ksizet+1,4)
c     complex gamval
c     complex zkappa
      complex*16, intent(in) :: u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3)
      complex*16, intent(inout) :: Phi(kthird, 0:ksize+1,
     &                                 0:ksize+1, 0:ksizet+1, 4)
      complex*16, intent(in) :: R(kthird, 0:ksize+1, 0:ksize+1,
     &                            0:ksizet+1, 4)
      integer, intent(in) :: imass
      real, intent(in) :: am
      complex*16 :: gamval
      complex*16 :: zkappa
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
      end

c***********************************************************************
      pure subroutine dslash2d(Phi,R,u)
      use purefunctions
c
c     calculates Phi = M*R
c
      implicit none
C      integer :: kdelta
      integer, parameter :: ksize=12, ksizet=12
      real, parameter :: akappa=0.5
      common/dirac/gamval(6,4),gamin(6,4)
c     complex u(ksize,ksize,ksizet,3)
c     complex Phi(ksize,ksize,ksizet,4),R(ksize,ksize,ksizet,4)
c     complex gamval
      complex*16, intent(in) ::  u(0:ksize+1,0:ksize+1,0:ksizet+1,3)
      complex*16, intent(inout) :: Phi(0:ksize+1,0:ksize+1,0:ksizet+1,4)
      complex*16, intent(in) :: R(0:ksize+1,0:ksize+1,0:ksizet+1,4)
      complex*16 :: gamval
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
      end
c
c***********************************************************************
c   Update boundary terms
c***********************************************************************
      pure subroutine update_halo_4(size4, Array)
c     
      implicit none
      integer, parameter :: ksize=12,ksizet=12
c
      integer, intent(in) :: size4
      complex*16, intent(inout) :: Array(0:ksize+1, 0:ksize+1,
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
      integer, parameter :: ksize=12,ksizet=12
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
      integer, parameter :: ksize=12,ksizet=12,kthird=24
c
      integer, intent(in) :: size5
      complex*16, intent(inout) :: Array(kthird, 0:ksize+1, 0:ksize+1,
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
c   A Kronecker delta function
c   Useful for calculating coordinate offsets
c***********************************************************************
      pure integer function kdelta(nu, mu)
        implicit none
        integer, intent(in) :: nu
        integer, intent(in) :: mu

        kdelta=merge(1,0,nu==mu)
      end

