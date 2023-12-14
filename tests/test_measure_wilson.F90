#include "test_utils.fh"
program test_measure_wilson
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
  integer :: i, imass, timing_loops = 1

  ! initialise function parameters
  complex(dp) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  real psibarpsi, aviter
  integer :: iflag, isweep, iter
  real :: res, am

  integer, parameter :: idxmax = 4*ksize*ksize*ksizet*kthird
  integer :: idx = 0
#ifdef MPI
  integer, dimension(16) :: reqs_x, reqs_Phi
  integer :: ierr
  call init_MPI
#endif
  seed = 4139764973254.0
  call init_random(seed)
  res = 0.1
  am = 0.05
  iflag = 0
  isweep = 1
  iter = 0
  am3 = 1.0
  imass = 3

#ifdef MPI
  call MPI_Barrier(comm, ierr)
#endif
  ! Phi is passed here simply to reduce the complexity of generate_starting_state_Phi_and_X
  call generate_starting_state_Phi_and_X(Phi, reqs_Phi, u, X, reqs_x)

  ! call function
  do i = 1, timing_loops
    call measure_wilson(psibarpsi, res, aviter, am, imass)
  end do

#ifndef SITE_RANDOM
  if (ip_global .eq. 0) then
    write (6, *) "This test is not supposed to work if SITE_RANDOM is not defined"
  endif
#endif
  check_float_equality(psibarpsi, -5.4226685E-02, 0.001, 'psibarpsi', 'test_measure_wilson')
  check_float_equality(aviter, 10.9, 0.001, 'aviter', 'test_measure_wilson')

#ifdef MPI
  call MPI_Finalize(ierr)
#endif

end program
