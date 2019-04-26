! Checking that the MPI_datatype sizes match the sizes of the partitions
#include "test_utils.fh"
program test_comms_partitioning
  use mpi
  use params
  use partitioning
  use comms_partitioning
  use comms
  implicit none
  integer :: ierr

  integer :: total_border_datasize, tempbds
  integer :: exp_total_border_datasize
  integer :: total_halo_datasize, temphds
  integer :: exp_total_halo_datasize

  ! 4 : dirac index
  ! kthird : obvious
  ! 16 : size if double complex in bytes
  integer, parameter :: size_factor=4*kthird*16
  integer :: ih,ib

  integer :: tempes ! TEMPorary Expected Size
  integer :: partition_volume




#ifdef MPI
  call init_MPI
#endif
  call init_partitioning()


  call MPI_Barrier(comm,ierr)
  call MPI_Barrier(comm,ierr)
  call init_dirac_border_types(dirac_border_dts,border_partitions_list)

  ! check border
  ! border volume
  exp_total_border_datasize = ksizet_l*ksizey_l*ksizex_l - &
    & (ksizet_l-2)*(ksizey_l-2)*(ksizex_l-2)
  ! mult by size-per-site
  exp_total_border_datasize = exp_total_border_datasize * size_factor

  total_border_datasize = 0
  do ib=1,26
    call MPI_type_size(dirac_border_dts(ib),tempbds,ierr)
    total_border_datasize = tempbds + total_border_datasize
    tempes = partition_volume(border_partitions_list(ib)%chunk)*size_factor
    if((tempes.ne.tempbds).and.(ip_global.eq.0))then
      print*,"Border Partition wrong datasize", ib
      print*,tempes,tempbds
    endif
  enddo

  if((total_border_datasize.ne.exp_total_border_datasize).and.(ip_global.eq.0))then
    print*,"Total border datasize not correct"
    print*,total_border_datasize,exp_total_border_datasize
  endif

  call MPI_Barrier(comm,ierr)
  call MPI_Barrier(comm,ierr)
  call init_dirac_halo_types(dirac_halo_dts,halo_partitions_list)


  
  exp_total_halo_datasize = 2*ksizex_l*ksizey_l+ &
                         &  2*ksizey_l*ksizet_l+ &
                         &  2*ksizet_l*ksizex_l
  
  exp_total_halo_datasize = exp_total_halo_datasize * size_factor
  total_halo_datasize = 0
  do ih=1,54
    call MPI_type_size(dirac_halo_dts(ih),temphds,ierr)
    total_halo_datasize = temphds + total_halo_datasize
    tempes = partition_volume(halo_partitions_list(ih)%chunk)*size_factor
    if((tempes.ne.temphds).and.(ip_global.eq.0))then
      print*,"Halo Partition wrong datasize", ih
      print*,tempes,temphds
    endif
 
  enddo

  if((total_halo_datasize.ne.exp_total_halo_datasize).and.(ip_global.eq.0))then
    print*,"Total halo datasize not correct"
    print*,total_halo_datasize,exp_total_halo_datasize
  endif

#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program

function partition_volume(p) result(vol)
  integer, intent(in) :: p(2,3)
  integer             :: vol

  vol = p(2,1)-p(1,1)+1
  vol = vol*(p(2,2)-p(1,2)+1)
  vol = vol*(p(2,3)-p(1,3)+1)
end function
