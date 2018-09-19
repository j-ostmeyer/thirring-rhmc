#include "test_utils.fh"
program test_locvol_partitioning
  use mpi_f08
  use params
  use locvol_partitioning
  use comms
  implicit none
  integer :: all_border_partitions(2,3,-1:1,-1:1,-1:1)
  integer :: id1,id2
  integer :: is1,is2
  integer :: ipx1,ipy1,ipt1
  integer :: ipx2,ipy2,ipt2
  integer :: vol,tempvol,expected_vol
  integer :: p1(2,3), p2(2,3)
  logical :: check
  integer :: ierr

 
#ifdef MPI
  call init_MPI
#endif
  call init_partitions()
  vol = 0
  do id1=1,3
    do is1=1,2
      p1 = united_border_partitions(:,:,is1,id1)
      call partition_check(check,p1)
      if(.not.check) then
        print *, "partition is bad:", id1,is1
      endif
      call partition_volume(tempvol,p1)
      vol = vol + tempvol
    enddo
  enddo

  expected_vol = ksizet_l*ksizey_l*ksizex_l - (ksizet_l-2)*(ksizey_l-2)*(ksizex_l-2)
  if(vol.ne.expected_vol) then
    print *, "Total border volume incorrect:",vol,expected_vol
  endif

  call partition_volume(tempvol,bulk_partition)
  expected_vol = (ksizet_l-2)*(ksizey_l-2)*(ksizex_l-2)
  if(tempvol.ne.expected_vol) then
    print *, "bulk volume incorrect:",vol,expected_vol
  endif

#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
  call get_all_partitions(all_border_partitions)
  do ipx1=-1,1
    do ipy1=-1,1
      do ipt1=-1,1
        do ipx2=-1,1
          do ipy2=-1,1
            do ipt2=-1,1
              if((ipx1.ne.ipx2).and.(ipy1.ne.ipy2).and.(ipt1.ne.ipt2)) then
                p1 = all_border_partitions(:,:,ipx1,ipy1,ipt1)
                p2 = all_border_partitions(:,:,ipx2,ipy2,ipt2)
                call partition_intersect(check,p1,p2)
                if (ip_global.eq.0) then
                  if((check)) then
                    write(6,"(A10,6I4)"), "problem",ipx1,ipy1,ipt1,ipx2,ipy2,ipt2
                    print *, p1
                    print *, p2
                  endif
                endif
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  ! intersection check
  do id1=1,3
    do id2=1,3
      do is1=1,2
        do is2=1,2
          if((id1.ne.id2).and.(is1.ne.is2)) then
            p1 = united_border_partitions(:,:,is1,id1)
            p2 = united_border_partitions(:,:,is2,id2)
            call partition_intersect(check,p1,p2)
            if (ip_global.eq.0) then
              if((check)) then
                print *, "partitions intersect", id1,is1,id2,is2
                print *, p1
                print *, p2
              endif
            endif
          endif
        enddo
      enddo
    enddo
  enddo

  do id1=1,3
    do is1=1,2
      p1 = united_border_partitions(:,:,is1,id1)
      call partition_intersect(check, bulk_partition,p1)
      if (ip_global.eq.0) then
        if((check)) then
          print *, "partitions intersect", id1,is1
          print *, p1
          print *, bulk_partition
        endif
      endif
    enddo
  enddo

 

#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program

subroutine partition_check(check,p)
  integer, intent(in) :: p(2,3)
  logical, intent(out) :: check
  integer :: idir

  check = .true.
  do idir=1,3
    check = check.and.(p(2,idir).ge.p(1,idir))
  enddo
  if(.not.check) then
    print *, p
  endif
end subroutine

subroutine partition_volume(vol,p)
  integer, intent(in) :: p(2,3)
  integer, intent(out) :: vol

  vol = p(2,1)-p(1,1)+1
  vol = vol*(p(2,2)-p(1,2)+1)
  vol = vol*(p(2,3)-p(1,3)+1)
end subroutine

subroutine partition_intersect(check,p1,p2)
  integer, intent(in) :: p1(2,3),p2(2,3)
  logical, intent(out) :: check
  logical :: acheck1, acheck2
  integer :: idir

  check = .true.
  do idir=1,3
    acheck1 = (p1(1,idir).le.p2(2,idir)).and.(p1(2,idir).ge.p2(1,idir))
    acheck2 = (p1(2,idir).ge.p2(1,idir)).and.(p1(2,idir).le.p2(2,idir))
    check = (acheck1.or.acheck2).and.check
  enddo

end subroutine

