#include "test_utils.fh"
program test_measure
  ! use dwf3d_lib
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

  ! NOTICE! When dwf3d is done remove the following and uncomment the use dwf3d_lib
  real(dp) :: seed

  ! general parameters
  integer :: timing_loops = 1
  complex, parameter :: iunit = cmplx(0, 1)
  real*8, parameter :: tau = 8*atan(1.0_8)

  ! initialise function parameters
  complex(dp) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  real psibarpsi, aviter
  integer :: imass, iflag, isweep, iter
  real :: res, am

  integer :: i, j, ix, iy, it, ithird
  integer, parameter :: idxmax = 4*ksize*ksize*ksizet*kthird
  integer :: idx = 0
#ifdef MPI
  integer, dimension(16) :: reqs_x, reqs_Phi
  integer, dimension(12) :: reqs_u
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

  call generate_starting_state(Phi, X, u, reqs_Phi, reqs_x)
 
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
