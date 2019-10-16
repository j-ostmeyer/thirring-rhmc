#include "test_utils.fh"
program test_measure
  use dwf3d_lib
  use trial, only: u
  use vector
  use comms
  use comms4
  use comms5
  use gammamatrices
  use random
  use measure_module
  use test_utils
  implicit none

  ! general parameters
  integer :: timing_loops = 1
  complex, parameter :: iunit = cmplx(0, 1)
  real*8, parameter :: tau = 8*atan(1.0_8)

  ! initialise function parameters
  real psibarpsi, aviter
  integer :: imass, iflag, isweep, iter
  real :: res, am

  integer :: i, j, ix, iy, it, ithird
  integer, parameter :: idxmax = 4*ksize*ksize*ksizet*kthird
  integer :: idx = 0
#ifdef MPI
  integer, dimension(12) :: reqs_x, reqs_u
  integer :: ierr
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

  do j = 1, 4
    do it = 1, ksizet_l
      do iy = 1, ksizey_l
        do ix = 1, ksizex_l
          do ithird = 1, kthird
            idx = ithird + (ip_x*ksizex_l + ix - 1)*kthird &
              & + (ip_y*ksizey_l + iy - 1)*kthird*ksize &
              & + (ip_t*ksizet_l + it - 1)*kthird*ksize*ksize &
              & + (j - 1)*kthird*ksize*ksize*ksizet
            X(ithird, ix, iy, it, j) = 0.5*exp(1.0)*exp(iunit*idx*tau/idxmax)
          enddo
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call start_halo_update_5(4, x, 0, reqs_x)
#endif
  do j = 1, 3
    do it = 1, ksizet_l
      do iy = 1, ksizey_l
        do ix = 1, ksizex_l
          idx = ip_x*ksizex_l + ix &
            & + (ip_y*ksizey_l + iy - 1)*ksize &
            & + (ip_t*ksizet_l + it - 1)*ksize*ksize &
            & + (j - 1)*ksize*ksize*ksizet
          u(ix, iy, it, j) = exp(iunit*idx*tau/idxmax)
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call start_halo_update_4(3, u, 0, reqs_u) ! to check.
  call complete_halo_update(reqs_x)
  call complete_halo_update(reqs_u)
#else
  call update_halo_5(4, X)
  call update_halo_4(3, u)
#endif

  ! initialise common variables

  call init_gammas()
  ! call function
  do i = 1, timing_loops
    call measure(psibarpsi, res, aviter, am, imass)
  end do
#ifdef SITE_RANDOM
  ! differing random numbers will throw off stochastic estimates like these
  check_float_equality(psibarpsi, 2.504295e-4, 0.001, 'psibarpsi', 'test_measure')
#else
  if (ip_global .eq. 0) then
    write (6, *) "This test is not supposed to work if SITE_RANDOM is not defined"
  endif
  check_float_equality(psibarpsi, 2.504295e-4, 0.001, 'psibarpsi', 'test_measure')
#endif
  check_equality(aviter, 5, 'aviter', 'test_measure')
#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program
