! Checking that the nonlocal parts of dirac split work as intended
! NO MPI really used here
! phi=1 here
#include "test_utils.fh"
program test_dirac_split
  use mpi
  use comms
  use dirac_split
  implicit none
  integer :: lc
  integer :: ierr

#ifdef MPI
  call init_MPI
#endif
  call get_dslash_work_ordering(dslash_work_ordering,.false.)
  call get_dslash_work_ordering(dslashd_work_ordering,.true.)
  do lc=1,27*7
    write(20+ip_global,"(4I4)") dslash_work_ordering(:,lc)
  enddo
  do lc=1,27*7
    write(60+ip_global,"(4I4)") dslashd_work_ordering(:,lc)
  enddo
 
#ifdef MPI
  call MPI_Finalize(ierr)
#endif

end program


