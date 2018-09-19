#include "test_utils.fh"
program test_partitioning
  use mpi_f08
  use params
  use partitioning
  use comms
  implicit none
  integer :: all_partitions(2,3,-2:2,-2:2,-2:2)
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
  call init_partitions_and_neighs()

  do ipx1=-2,2
    do ipy1=-2,2
      do ipt1=-2,2
        do ipx2=-2,2
          do ipy2=-2,2
            do ipt2=-2,2
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

  vol = 0
  do ipx1=-2,2
    do ipy1=-2,2
      do ipt1=-2,2
        p1 = all_border_partitions(:,:,ipx1,ipy1,ipt1)
        partition_volume(tempvol,p1)
        vol = vol + tempvol
      enddo
    enddo
  enddo
  expected_vol = (ksizet_l+2)*(ksizey_l+2)*(ksizex_l+2)
  if(expected_vol.ne.vol) then
    if(ip_global.ne.0)then
      print*, "Volume is not correct."
    endif
  endif


 

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

