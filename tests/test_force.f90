program test_force
      implicit none
! function to test
      external :: force

! general parameters
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksize*ksize*ksizet)
      parameter(kvol2=ksize*ksize)
      parameter(kferm=4*kvol*kthird)
      parameter(ndiagg=12,ndiag=ndiagg)
      parameter(One=1.0)
      parameter(Nf=1)
      integer :: ksize, ksizet, kthird, kvol, kvol2, kferm, ndiagg, ndiag, Nf
      real :: One

! common blocks to function
      common/remez2g/anum2(0:ndiag),aden2(ndiag),
      &              bnum2(0:ndiag),bden2(ndiag)
      common/remez4g/anum4(0:ndiag),aden4(ndiag),
      &              bnum4(0:ndiag),bden4(ndiag)
      common/phi0/Phi0(kferm,25)
      common/trial/u(kvol,3),theta(kvol,3),pp(kvol,3)
      common/para/beta,am3,ibound
      common/param/ancg,ancgh,ancgf,ancgpf
      common/parampv/ancgpv,ancghpv,ancgfpv,ancgpfpv
      common/gforce/dSdpi(kvol,3)
      common/vector/X1(kferm)
      real*8 anum2,aden2,bnum2,bden2
      real*8 anum4,aden4,bnum4,bden4
      real*8 anum2g,aden2g,bnum2g,bden2g
      real*8 anum4g,aden4g,bnum4g,bden4g
      complex*16 Phi0, u, X1
      real theta, pp, beta, am3, ibound
      real ancg, ancgh, ancgh, ancgpf
      real ancgpv, ancghpv, ancgfpv, ancgpfpv
      real dSdpi

! initialise function parameters
      complex*16 Phi(kferm,Nf)

      integer :: i, j
      do i = 1,kferm
         do j = 1,Nf
            idx = (i-1) * Nf + j
            Phi(i, j) = ((-1) ** (idx / 2), (-1) ** ((idx + 1) / 2))
         enddo
      enddo
      
      real*8 rescgg = 0.000001
      real*8 am = 0.05
      integer imass = 3
      integer isweep = 1
      integer iter = 1
      
! initialise common variables
      beta = 0.4
      am3 = 1.0
      ibound = -1
      ancg=0.0
      ancgh=0.0
      ancgf=0.0
      ancgpf=0.0
      ancgpv=0.0
      ancghpv=0.0
      ancgfpv=0.0
      ancgpfpv=0.0
      ancgma=0.0

      do i = 0,ndiag
         anum2(i) = 2 ** (-i)
         anum4(i) = 3 ** (-i)
         bnum2(i) = 5 ** (-i)
         bnum4(i) = 7 ** (-i)
      enddo
      do i = 1,ndiag
         aden2(i) = 11 ** (-i)
         aden4(i) = 13 ** (-i)
         bden2(i) = 17 ** (-i)
         bden4(i) = 19 ** (-i)
      enddo
      
      do i = 1,kferm
         do j = 1,25
            idx = (i - 1) * 25 + j
            Phi0(i,j) = 1.1 * ((-1) ** (idx / 2), (-1) ** ((idx + 1) / 2))
         enddo
         X1(i) = ((-1) ** (i / 2), (-1) ** ((i + 1) / 2))
      enddo

      do i = 1,kvol
         do j = 1,3
            idx = (i - 1) * 3 + j
            u(i, j) = 1.3 * ((-1) ** (idx / 2), (-1) ** ((idx + 1) / 2))
            theta(i, j) = (-1) ** idx
            pp(i, j) = (-1) ** idx * 1.1
            dSdpi(i, j) = (-1) ** idx * 1.3
         enddo
      enddo

      call force(Phi, res1, am, imass, isweep, iter)
      
      do i = 1,kferm,(kferm/10)
         print *,'Phi(', i, ',', 1, ') = ', Phi(i, 1)
      enddo
end program
