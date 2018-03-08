program test_qmrherm
      implicit none
! function to test
      external :: qmrherm

! supporting functions
      external :: init
      external :: update_halo_4
      external :: update_halo_5
      external :: update_halo_6

! general parameters
      logical :: generate = .false.
      integer :: timing_loops = 100
      integer, parameter :: ksize=12, ksizet=12, kthird=24
      integer, parameter :: ndiag = 25
      complex, parameter :: iunit = cmplx(0, 1)
      real*8, parameter :: tau = 8 * atan(1.0_8)
      complex*16 :: acc_sum = 0.
      real*8 :: acc_max = 0.

! common blocks to function
      common/trial/u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3), &
           theta(kthird, ksize, ksize, ksizet, 3), &
           pp(kthird, ksize, ksize, ksizet, 3)
      common/para/beta,am3,ibound
      common/vector/x(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      common/phi0/Phi0(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4, 25)
      common/dirac/gamval(6,4),gamin(6,4)
      common/gforce/dSdpi(ksize, ksize, ksizet, 3)
      real :: beta, am3
      integer :: ibound, istart
      complex*16 :: gamval, x, Phi0
      integer :: gamin
      integer :: iu, id
      real :: dSdpi
      real :: theta, pp

! initialise function parameters
      complex*16 :: u
      complex*16 Phi(kthird,0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex*16 X2(kthird,0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex*16 :: Phi0_ref(kthird, ksize, ksize, ksizet, 4, 25)
      complex*16 :: x_ref(kthird, ksize, ksize, ksizet, 4)
      complex*16 :: delta_x(kthird, ksize, ksize, ksizet, 4)
      complex*16 :: ratio_x(kthird, ksize, ksize, ksizet, 4)
      complex*16 :: delta_Phi0(kthird, ksize, ksize, ksizet, 4, 25)
      complex*16 R(kthird,0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      real :: diff(ksize, ksize, ksizet, 3)

      integer :: imass, iflag, isweep, iter
      real*8 :: anum(0:ndiag), aden(ndiag), coeff
      real :: res, am
      integer :: itercg
      
      integer :: i, j, l, ix, iy, it, ithird
      integer, parameter :: idxmax = 4 * ksize * ksize * ksizet * kthird
      integer :: idx = 0

      res = 0.1
      am = 0.05
      imass = 3
      iflag = 0
      isweep = 1
      iter = 0

      anum(0) = 0.5
      do i = 1, ndiag
         anum(i) = exp(iunit * i * tau / ndiag)
         aden(i) = exp(-iunit * 0.5 * i * tau / ndiag)
      enddo
      do j = 1,4
         do it = 1,ksizet
            do iy = 1,ksize
               do ix = 1,ksize
                  do ithird = 1,kthird
                     idx = idx + 1
                     Phi(ithird, ix, iy, it, j) = 1.1 * exp(iunit * idx * tau / idxmax)
                     R(ithird, ix, iy, it, j) = 1.3 * exp(iunit * idx * tau / idxmax)
                     do l = 1, 25
                        Phi0_ref(ithird, ix, iy, it, j, l) = 1.7 * exp(1.0) * exp(iunit * idx * tau / idxmax) + l
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      idx = 0
      do j = 1,3
         do it = 1,ksizet
            do iy = 1,ksize
               do ix = 1,ksize
                  idx = idx + 1
                  u(ix, iy, it, j) = exp(iunit * idx * tau / idxmax)
                  dSdpi(ix, iy, it, j) = tau * exp(iunit * idx * tau / idxmax)
               enddo
            enddo
         enddo
      enddo
      call update_halo_6(4, 25, Phi0)
      call update_halo_5(4, Phi)
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
         Phi0(:, 1:ksize, 1:ksize, 1:ksizet,:, :) = Phi0_ref
         call qmrherm(Phi,res,itercg,am,imass,anum,aden,ndiag,iflag,isweep,iter)
!         print *, "CG iterations:", itercg
      end do
! check output
      open(3, file='test_qmrherm_x.dat', form="unformatted", access="sequential")
      open(4, file='test_qmrherm_Phi0.dat', form="unformatted", access="sequential")
      if (generate) then
         write(3) x(:,1:ksize,1:ksize,1:ksizet,:)
         write(4) Phi0(:,1:ksize,1:ksize,1:ksizet,:,:)
      else
         read(3) x_ref
         read(4) Phi0_ref
         delta_x = x(:,1:ksize,1:ksize,1:ksizet,:) - x_ref
         delta_Phi0 = Phi0(:,1:ksize,1:ksize,1:ksizet,:,:) - Phi0_ref
         ratio_x = x_ref / x(:, 1:ksize, 1:ksize, 1:ksizet, :)
         print *, 'sum delta = ', sum(delta_x), sum(delta_Phi0)
         print *, 'max delta = ', maxval(abs(delta_x)), maxval(abs(delta_Phi0))
         print *, 'sum ratio = ', sum(ratio_x), 'max = ', maxval(abs(ratio_x))
!         do j=1,3
!            do it=1,ksizet
!               do iy=1,ksize
!                  do ix=1,ksize
!                     print *, ix, iy, it, j, diff(ix, iy, it, j)
!                  enddo
!               enddo
!            enddo
!         enddo
      
   end if
end program
