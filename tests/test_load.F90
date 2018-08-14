#include "test_utils.fh"
program test_load
  use dwf3d_lib
  use gauge
  use comms
  implicit none

  real :: sumtheta, maxtheta, mintheta

  ! initialise MPI
#ifdef MPI
  integer :: ierr
  call init_MPI
#endif

  ! call function
  call sread

  ! check output
  sumtheta = sum(theta)
  maxtheta = maxval(theta)
  mintheta = minval(theta)
#ifdef MPI
  call MPI_AllReduce(MPI_IN_PLACE, sumtheta, 1, MPI_REAL, MPI_SUM, comm,ierr)
  call MPI_AllReduce(MPI_IN_PLACE, maxtheta, 1, MPI_REAL, MPI_MAX, comm,ierr)
  call MPI_AllReduce(MPI_IN_PLACE, mintheta, 1, MPI_REAL, MPI_MIN, comm,ierr)
#endif
  check_float_equality(sumtheta, -185.5681, 0.001, 'sum', 'test_load')
  check_float_equality(maxtheta, 5.015248, 0.001, 'max', 'test_load')
  check_float_equality(mintheta, -5.267587, 0.001, 'min', 'test_load')
  check_float_equality(seed, 196829593468928., 0.001, 'seed', 'test_load')
#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program
