      subroutine qmrherm(Phi,res,itercg,am,imass,anum,aden,ndiag,
     &                   iflag,isweep,iter)
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
      parameter(kferm=4*kthird*kvol)
c     parameter(niterc=10*kferm)
      parameter(niterc=7500)
      common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
      common/para/bbb,am3,ibound
      common/vector/x(kferm)
      common/gforce/dSdpi(kvol,3)
      common/phi0/Phi0(kferm,25)
c     complex Phi(kferm)
c     complex x,u,Phi0
      complex*16 Phi(kferm)
      complex*16 x,u,Phi0
      real*8 anum(0:ndiag),aden(ndiag),coeff
c     complex vtild(kferm),q(kferm),pm1(kferm,ndiag)
c     complex qm1(kferm),p(kferm,ndiag),x3(kferm),R(kferm)
c     complex x1(kferm,ndiag),x2(kferm)
c     real alpha(ndiag),beta
c     real amu(ndiag),d(ndiag),dm1(ndiag),rho(ndiag),rhom1(ndiag)
      complex*16 vtild(kferm),q(kferm),pm1(kferm,ndiag)
      complex*16 qm1(kferm),p(kferm,ndiag),x3(kferm),R(kferm)
      complex*16 x1(kferm,ndiag),x2(kferm)
      real*8 alpha(ndiag),beta
      real*8 amu(ndiag),d(ndiag),dm1(ndiag),rho(ndiag),rhom1(ndiag)
c
c     write(6,111)
111   format(' Hi from qmrherm')
c
      resid=sqrt(kferm*res*res)
c     write(6,*) iflag, resid
      itercg=0
c
c   initialise r=Phi
c
      do i=1,kferm
         R(i)=Phi(i)
         qm1(i)=(0.0,0.0)
         x(i)=anum(0)*Phi(i)
      enddo
      beta=0.0
      do i=1,kferm
         beta=beta+conjg(r(i))*r(i)
      enddo
      beta=sqrt(beta)
      phimod=beta
c     write(6,*) '|| Phi || = ', phimod
c
      do niter=1,niterc
      itercg=itercg+1
c
c  Lanczos steps
c
      do i=1,kferm
         q(i)=R(i)/beta
      enddo
c
      call dslash(vtild,q,u,am,imass)
      call dslashd(x3,vtild,u,am,imass)
c
      alphatild=0.0
      do i=1,kferm
         alphatild=alphatild+conjg(q(i))*x3(i)
      enddo
c
      do i=1,kferm
         R(i)=x3(i)-alphatild*q(i)-beta*qm1(i)
         qm1(i)=q(i)
      enddo
c
      beta0=beta
      beta=0.0
      do i=1,kferm
         beta=beta+conjg(R(i))*R(i)
      enddo
      beta=sqrt(beta)
c
      do idiag=1,ndiag
      alpha(idiag)=alphatild+aden(idiag)
      enddo
c
      if(niter.eq.1)then
           do idiag=1,ndiag
             d(idiag)=alpha(idiag)
             do i=1,kferm
              p(i,idiag)=q(i)
              pm1(i,idiag)=p(i,idiag)
             enddo
             rho(idiag)=beta0/alpha(idiag)
             rhom1(idiag)=rho(idiag)
             do i=1,kferm
              x1(i,idiag)=rho(idiag)*q(i)
             enddo
           enddo
      else
           do idiag=1,ndiag
             amu(idiag)=beta0/d(idiag)
             dm1(idiag)=d(idiag)
             d(idiag)=alpha(idiag)-beta0*amu(idiag)
             do i=1,kferm
              p(i,idiag)=q(i)-amu(idiag)*pm1(i,idiag)
              pm1(i,idiag)=p(i,idiag)
             enddo
             rho(idiag)=-amu(idiag)*dm1(idiag)*rhom1(idiag)/d(idiag)
c  Convergence criterion (a bit ad hoc for now...)
             if(idiag.eq.1)then
               rhomax=abs(phimod*rho(idiag))
             else
               if(abs(phimod*rho(idiag)).gt.rhomax)
     &             rhomax=abs(phimod*rho(idiag))
             endif
             rhom1(idiag)=rho(idiag)
             do i=1,kferm
               x1(i,idiag)=x1(i,idiag)+rho(idiag)*p(i,idiag)
             enddo
           enddo
c  check to see whether the residual is acceptable for all ndiag....
c  criterion is a bit ad hoc -- relaxing by a factor arelax improves code
c  stability and leads to quicker convergence
           arelax=2.0
           if(rhomax.lt.arelax*resid) then
c          if(rhomax.lt.resid) then
c          call testinv(Phi,resmax,itercg,am,imass,x1,aden,ndiag)
c  convergence based on || residual || not working well in single precision...
c          if(resmax.lt.resid) goto 8
             goto 8
           endif
      endif
c
c   end of loop over iter
      enddo
      write(7,*) 'QMRniterc!, isweep,iter,iflag,imass,anum,ndiag = '
     &        ,isweep,iter,iflag,imass,anum(0),ndiag
8     continue
c
      if(iflag.lt.2)then
c  Now evaluate solution x=(MdaggerM)^p * Phi
      do idiag=1,ndiag
      do i=1,kferm
      x(i)=x(i)+anum(idiag)*x1(i,idiag)
      enddo
      enddo
c
c  update phi0 block if required...
      if(iflag.eq.1) then
      do idiag=1,ndiag
      do i=1,kferm
      Phi0(i,idiag)=x1(i,idiag)
      enddo
      enddo
      endif
c
      else
c
      do idiag=1,ndiag
c
c  X2 = M*X1
      do i=1,kferm
      R(i)=x1(i,idiag)
      enddo
      call dslash(X2,R,u,am,imass)
c
      if(iflag.eq.2)then
          coeff=anum(idiag)
          call derivs(R,X2,coeff,0)
      else
          coeff=-anum(idiag)
          do i=1,kferm
          R(i)=Phi0(i,idiag)
          enddo
          call derivs(R,X2,coeff,0)
c
          call dslash(X2,R,u,am,imass)
          do i=1,kferm
          R(i)=x1(i,idiag)
          enddo
          call derivs(X2,R,coeff,1)
      endif
c
      enddo
      endif
c
c
      return
      end
