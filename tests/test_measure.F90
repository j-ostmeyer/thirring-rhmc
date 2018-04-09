program test_measure
      use dwf3d_lib
      use trial, only: u
      use vector
      use comms
      use random
      implicit none

! general parameters
      integer :: timing_loops = 1
      complex, parameter :: iunit = cmplx(0, 1)
      real*8, parameter :: tau = 8 * atan(1.0_8)
      complex(dp) :: acc_sum = 0.
      real*8 :: acc_max = 0.
      integer :: istart

! initialise function parameters
      real psibarpsi, aviter
      integer :: imass, iflag, isweep, iter
      real :: res, am
      real(dp) :: h, hg, hp, s
      integer :: itercg
      
      integer :: i, j, l, ix, iy, it, ithird
      integer, parameter :: idxmax = 4 * ksize * ksize * ksizet * kthird
      integer :: idx = 0
#ifdef MPI
      integer, dimension(12) :: reqs_x, reqs_u

      call init_MPI
#endif
      seed = 4139764973254.0
      call init_random(seed)
      res = 0.1
      am = 0.05
      imass = 3
      iflag = 0
      isweep = 1
      iter = 0
      am3 = 1.0

      do j = 1,4
         do it = 1,ksizet_l
            do iy = 1,ksizey_l
               do ix = 1,ksizex_l
                  do ithird = 1,kthird
                     idx = ithird + (ip_x * ksizex_l + ix - 1) * kthird &
                          & + (ip_y * ksizey_l + iy - 1) * kthird * ksize &
                          & + (ip_t * ksizet_l + it - 1) * kthird * ksize * ksize &
                          & + (j - 1) * kthird * ksize * ksize * ksizet
                     X(ithird, ix, iy, it, j) = 0.5 * exp(1.0) * exp(iunit*idx*tau / idxmax)
                  enddo
               enddo
            enddo
         enddo
      enddo
#ifdef MPI
      call start_halo_update_5(4, x, 0, reqs_x)
#endif
      do j = 1,3
         do it = 1,ksizet_l
            do iy = 1,ksizey_l
               do ix = 1,ksizex_l
                  idx = ip_x * ksizex_l + ix &
                       & + (ip_y * ksizey_l + iy - 1) * ksize &
                       & + (ip_t * ksizet_l + it - 1) * ksize * ksize &
                       & + (j - 1) * ksize * ksize * ksizet
                  u(ix, iy, it, j) = exp(iunit * idx * tau / idxmax)
               enddo
            enddo
         enddo
      enddo
#ifdef MPI
      call start_halo_update_4(4, u, 0, reqs_u)
      call complete_halo_update(reqs_x)
      call complete_halo_update(reqs_u)
#else
      call update_halo_5(4, X)
      call update_halo_4(3, u)
#endif

! initialise common variables
      istart = -1
      call init(istart)
! call function
      do i = 1,timing_loops
         call measure(psibarpsi, res, aviter, am, imass)
      end do
      if (ip_global .eq. 0) then
         print *, psibarpsi, aviter
      end if
#ifdef MPI
      call MPI_Finalize(ierr)
#endif
end program
