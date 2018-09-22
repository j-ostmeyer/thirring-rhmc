! Checking that the MPI_datatype sizes match the sizes of the partitions
#include "test_utils.fh"
program test_dirac_split
  use mpi_f08
  use params
  use partitioning
  use comms_partitioning
  use comms
  use dirac_split
  implicit none
  integer :: ierr

  !BUFfer To Send
  complex(dp) :: bufts(kthird,0:ksizex_l+1,0:ksizey_l+1,0:ksizet_l+1,4)
  !BUFfer To Recv
  complex(dp) :: buftr(kthird,0:ksizex_l+1,0:ksizey_l+1,0:ksizet_l+1,4)


#ifdef MPI
  call init_MPI
#endif
  call init_partitioning()
  call init_dirac_border_types(dirac_border_dts,border_partitions_list)
  call init_dirac_halo_types(dirac_halo_dts,halo_partitions_list)

  ! TODO


 #ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program
