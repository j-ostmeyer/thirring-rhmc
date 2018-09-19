#include "test_utils.fh"
program test_dslash
  use partitioning
  use comms
  implicit none

#ifdef MPI
  call init_MPI
#endif
  call dochecks()
#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program

  subroutine dochecks()
    integer:: united_border_partitions(2,3,2,3)
    integer :: id1,id2
    integer :: is1,is2
    integer :: vol,tempvol,expected_vol
    integer :: p1(2,3), p2(2,3)
    logical :: check

    get_united_border_partitions(united_border_partitions)
    vol = 0
    do id1=1,3
      do is1=1,2
        p1 = united_border_partitions(:,:,is1,id1)
        partition_check(check,p1)
        if(.not.check) then
          print *, "partition is bad:", id1,is1
        endif
        partition_volume(tempvol,p)
        vol = vol + tempvol
      enddo
    enddo

    expected_vol = ksizet_l*ksizey_l*ksizex_l
    if(vol.ne.expected_vol) then
      print *, "Total volume incorrect:",vol,expected_vol
    endif

    ! intersection check
    do id1=1,3
      do id2=1,3
        do is1=1,2
          do is2=1,2
            p1 = united_border_partitions(:,:,is1,id1)
            p2 = united_border_partitions(:,:,is2,id2)
            partition_intersect(check,p1,p2)
            if(.not.check) then
              print *, "partitions intersect", id1,is1,id2,is2
            endif
          enddo
        enddo
      enddo
    enddo







  end subroutine


