#include "test_utils.fh"
program test_measure
  use gdbhook
  use dwf3d_lib
  use trial, only: u
  use vector
  use dirac
  use comms
  use comms4
  use comms5
  use gammamatrices
  use measure_module
  use test_utils
  use gdbhook
  implicit none

  ! general parameters
  integer :: timing_loops = 1

  ! initialise function parameters
  complex(dp) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

  ! initialise function parameters
  real :: aviter
  integer :: iflag = 0

  integer :: i
  real :: res, am
  integer :: imass, itercg
#ifdef MPI
  integer, dimension(16) :: reqs_Phi
  integer :: ierr

  call init_MPI
  call gdb_wait()
#endif
  seed = 4139764973254.0
  idum = -1
  call rranset(seed, 1, 1, 1)
  seed = rano(yran, idum, 1, 1, 1)

  res = 0.1
  am = 0.05
  imass = 3
  iflag = 0

  call generate_starting_state_Phi(Phi, u, reqs_Phi)

  ! initialise common variables
  seed = rano(yran, idum, 1, 1, 1)
  seed = rano(yran, idum, 1, 1, 1)

  ! call function
  do i = 1, timing_loops
    x = (0.D0, 0.D0)
    call meson(res, itercg, aviter, am, imass)
  end do

#ifdef SITE_RANDOM
  ! differing random numbers will throw off stochastic estimates like these
  check_float_equality(psibarpsi, 2.504295e-4, 0.001, 'psibarpsi', 'test_meson')
#else
  if (ip_global .eq. 0) then
    write (6, *) "This test is not supposed to work if SITE_RANDOM is not defined"
  endif
  check_float_equality(psibarpsi, 2.504295e-4, 0.001, 'psibarpsi', 'test_meson')
#endif

  check_equality(aviter, 5, 'aviter', 'test_meson')

#ifdef MPI
    endif! if(ip_global .eq. 0) then
    call MPI_Finalize(ierr)
#endif
  end if
  close (3)
end program
