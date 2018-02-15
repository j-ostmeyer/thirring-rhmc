program test_dslashd
      implicit none
! function to test
      external :: dslashd

! supporting functions
      external :: init
      external :: update_halo_4
      external :: update_halo_5

! general parameters
      logical :: generate = .false.
      integer :: timing_loops = 1000
      integer, parameter :: ksize=12, ksizet=12, kthird=24
      complex, parameter :: iunit = cmplx(0, 1)
      real*8, parameter :: tau = 8 * atan(1.0_8)
      complex*16 :: acc_sum = 0.
      real*8 :: acc_max = 0.

! common blocks to function
      common/para/beta,am3,ibound
      common/dirac/gamval(6,4),gamin(6,4)
      real :: beta, am3
      integer :: ibound, istart
      complex*16 :: gamval
      integer :: gamin
      integer :: iu, id

! initialise function parameters
      complex*16 u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3)
      complex*16 Phi(kthird,0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex*16 Phiref(kthird,ksize, ksize, ksizet, 4)
      complex*16 R(kthird,0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex*16 diff(kthird,ksize, ksize, ksizet, 4)

      real*8, parameter :: am = 0.05
      integer, parameter :: imass = 3

      integer :: i, j, l, ix, iy, it, ithird
      integer, parameter :: idxmax = 4 * ksize * ksize * ksizet * kthird
      integer :: idx = 0
      do j = 1,4
         do it = 1,ksizet
            do iy = 1,ksize
               do ix = 1,ksize
                  do ithird = 1,kthird
                     idx = idx + 1
                     Phi(ithird, ix, iy, it, j) = 1.1 * exp(iunit * idx * tau / idxmax)
                     R(ithird, ix, iy, it, j) = 1.3 * exp(iunit * idx * tau / idxmax)
                  enddo
               enddo
            enddo
         enddo
      enddo
      do j = 1,3
         idx = 4 * ksize * ksize * ksizet * kthird + (j-1) * kthird
         do it = 1,ksizet
            do iy = 1,ksize
               do ix = 1,ksize
                  idx = idx + 1
                  u(ix, iy, it, j) = exp(iunit * idx * tau / idxmax)
               enddo
            enddo
         enddo
      enddo
!      call update_halo_5(4, Phi)
      call update_halo_5(4, R)
      call update_halo_4(3, u)

! initialise common variables
      beta = 0.4
      am3 = 1.0
      ibound = -1
      istart = -1
      call init(istart)
! call function
      do i = 1,timing_loops
         call dslashd(Phi, R, u, am, imass)
      end do
! check output
!      do i = 1,10
!         j = 1 + i * (kvol - 1) / 10
!         l = 1 + i * (4 - 1) / 10
!         print *,'Phi(', j, ',', l, ') = ', Phi(j, l)
!      enddo

      open(3, file='test_dslashd.dat', form="unformatted", access="sequential")
      if (generate) then
         write(3) Phi(:,1:ksize,1:ksize,1:ksizet,:)
      else
         read(3) Phiref
         
         diff = Phi(:,1:ksize,1:ksize,1:ksizet,:) - Phiref
         print *, 'sum delta = ', sum(diff)
         print *, 'max delta = ', maxval(abs(diff))
      end if
end program
