program test_measure
      use dwf3d_lib
      implicit none

! general parameters
      logical :: generate = .false.
      integer :: timing_loops = 1
      integer, parameter :: ndiag = 12
      complex, parameter :: iunit = cmplx(0, 1)
      real*8, parameter :: tau = 8 * atan(1.0_8)
      complex(dp) :: acc_sum = 0.
      real*8 :: acc_max = 0.

! common blocks to function
      common/trial/u(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3), &
           theta(ksize, ksize, ksizet, 3), &
           pp(ksize, ksize, ksizet, 3)
      common/para/beta,am3,ibound
      common/vector/x(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      common/dirac/gamval(6,4),gamin(6,4)
      common/gforce/dSdpi(ksize, ksize, ksizet, 3)
      common/remez2g/anum2(0:ndiag),aden2(ndiag),bnum2(0:ndiag),bden2(ndiag)
      common/remez4g/anum4(0:ndiag),aden4(ndiag),bnum4(0:ndiag),bden4(ndiag)
      common/param/ancg,ancgh,ancgf,ancgpf
      common/parampv/ancgpv,ancghpv,ancgfpv,ancgpfpv
      common/ranseed/idum
      integer :: idum = -1
      real :: beta, am3
      integer :: ibound, istart
      complex(dp) :: gamval, x
      integer :: gamin
      integer :: iu, id
      real :: dSdpi
      real :: dSdpi_ref(ksize, ksize, ksizet, 3)
      real :: theta, pp
      real :: ancg,ancgh,ancgf,ancgpf
      real :: ancgpv,ancghpv,ancgfpv,ancgpfpv
      real*8 :: anum2, aden2, bnum2, bden2, anum4, aden4, bnum4, bden4
      complex(dp) :: u

! initialise function parameters
      complex(dp) :: Phi(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      complex(dp) :: Phi0(kthird, 0:ksize+1, 0:ksize+1, 0:ksizet+1, 4, 25)
      complex(dp) :: Phi0_ref(kthird, ksize, ksize, ksizet, 4, 25)
      complex(dp) :: x_ref(kthird, ksize, ksize, ksizet, 4)
      complex(dp) :: delta(kthird, ksize, ksize, ksizet, 4)
      complex(dp) :: R(kthird,0:ksize+1, 0:ksize+1, 0:ksizet+1, 4)
      real :: diff(ksize, ksize, ksizet, 3)
      real psibarpsi, aviter
      integer :: imass, iflag, isweep, iter
      real(dp) :: anum(0:ndiag), aden(ndiag), coeff
      real :: res, am
      real(dp) :: h, hg, hp, s
      integer :: itercg
      
      integer :: i, j, l, ix, iy, it, ithird
      integer, parameter :: idxmax = 4 * ksize * ksize * ksizet * kthird
      integer :: idx = 0

      h = 0
      hg = 0
      hp = 0
      s = 0
      res = ran(idum)
      res = 0.1
      am = 0.05
      imass = 3
      iflag = 0
      isweep = 1
      iter = 0

      anum2(0) = 0.5
      anum4(0) = 0.51
      bnum2(0) = 0.49
      bnum4(0) = 0.53
      do i = 1, ndiag
         anum2(i) = 0.4 * exp(iunit * i * tau / ndiag)
         aden2(i) = 0.4 * exp(-iunit * 0.5 * i * tau / ndiag)
         anum4(i) = 0.41 * exp(iunit * i * tau / ndiag)
         aden4(i) = 0.41 * exp(-iunit * 0.5 * i * tau / ndiag)
      enddo
      do j = 1,4
         do it = 1,ksizet
            do iy = 1,ksize
               do ix = 1,ksize
                  do ithird = 1,kthird
                     idx = idx + 1
                     Phi(ithird, ix, iy, it, j) = 1.1 * exp(iunit * idx * tau / idxmax)
                     R(ithird, ix, iy, it, j) = 1.3 * exp(iunit * idx * tau / idxmax)
                     X(ithird, ix, iy, it, j) = 0.5 * exp(1.0) * exp(iunit * idx * tau / idxmax)
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
                  theta(ix, iy, it, j) = 1.9 * exp(iunit * idx * tau / idxmax)
                  pp(ix, iy, it, j) = -1.1 * exp(iunit * idx * tau / idxmax)
                  dSdpi(ix, iy, it, j) = tau * exp(iunit * idx * tau / idxmax)
               enddo
            enddo
         enddo
      enddo
      call update_halo_5(4, Phi)
      call update_halo_5(4, R)
      call update_halo_5(4, X)
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
         h = 0
         hg = 0
         hp = 0
         s = 0
         call update_halo_6(4, 25, Phi0)
         call measure(psibarpsi, res, aviter, am, imass)
      end do
      print *, psibarpsi, aviter
! check output
!      open(3, file='test_congrad.dat', form="unformatted", access="sequential")
!      if (generate) then
!         write(3) x(:, 1:ksize, 1:ksize, 1:ksizet, :)
!      else
!         read(3) x_ref
!         delta = x_ref - x(:, 1:ksize, 1:ksize, 1:ksizet, :)
!         print *, 'sum delta = ', sum(delta)
!         print *, 'max delta = ', maxval(abs(delta))
!   end if
end program
