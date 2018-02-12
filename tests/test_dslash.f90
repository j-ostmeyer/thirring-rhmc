program test_dslash
      implicit none
! function to test
      external :: dslash

! supporting functions
      external :: init

! general parameters
      integer :: ksize, ksizet, kthird, kvol
      complex, parameter :: iunit = cmplx(0, 1)
      parameter(ksize=12,ksizet=12,kthird=24,kvol=ksizet*ksize*ksize)
      real*8, parameter :: tau = 8 * atan(1.0_8)

! common blocks to function
      common/para/beta,am3,ibound
      common/dirac/gamval(6,4),gamin(6,4)
      common /neighb/id(kvol,3),iu(kvol,3)
      real :: beta, am3
      integer :: ibound, istart
      complex*16 :: gamval
      integer :: gamin
      integer :: iu, id

! initialise function parameters
      complex*16 Phi(kthird,kvol,4), R(kthird,kvol,4)
      complex*16 u(kvol, 3)
      real*8 :: am = 0.05
      integer :: imass = 3

      integer :: i, j, k, l, idx
      real*8 :: repart, impart
      integer :: idxmax = 4 * kvol * kthird
      do k = 1,4
         do j = 1,kvol
            do i = 1,kthird
               idx = (k-1) * kvol * kthird + (j-1) * kthird + i
               Phi(i, j, k) = 1.1 * exp(iunit * idx * tau / idxmax)
               R(i, j, k) = 1.3 * exp(iunit * idx * tau / idxmax)
            enddo
         enddo
      enddo
      do j = 1,3
         do i = 1,kvol
            idx = (k-1) * kvol * kthird + (j-1) * kthird + i
            u(i, j) = exp(iunit * idx * tau / idxmax)
         enddo
      enddo
      
! initialise common variables
      beta = 0.4
      am3 = 1.0
      ibound = -1
      istart = -1
      call init(istart)
      
! call function
      call dslash(Phi, R, u, am, imass)
      
      do i = 1,10
         j = 1 + i * (kthird - 1) / 10
         k = 1 + i * (kvol - 1) / 10
         l = 1 + i * (4 - 1) / 10
         print *,'Phi(', j, ',', k, ',', l, ') = ', Phi(j, k, l)
      enddo
end program
