! Checking that the nonlocal parts of dirac split work as intended
! NO MPI really used here
! phi=1 here
#include "test_utils.fh"
program test_dirac_split
  use mpi_f08
  use comms
  use dirac_split
  implicit none
  integer :: lc

#ifdef MPI
  call init_MPI
#endif
  if(ip_global.eq.0)then
    call get_dslash_work_ordering(dslash_work_ordering)
    do lc=1,27*7
      print*,dslash_work_ordering(:,lc)
    enddo
  endif
 
#ifdef MPI
  call MPI_Finalize(ierr)
#endif

end program


