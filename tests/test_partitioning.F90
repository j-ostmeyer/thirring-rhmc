! Checking that partitions do not intersect
! but that the sum of the volumes is correct
#include "test_utils.fh"
program test_partitioning
  use mpi
  use params
  use partitioning
  use comms
  implicit none
  integer :: ipx1,ipy1,ipt1
  integer :: ipx2,ipy2,ipt2
  integer :: vol,expected_vol
  integer :: p1(2,3), p2(2,3)
  integer :: ierr
  logical :: partition_intersect
  integer :: partition_volume

 
#ifdef MPI
  call init_MPI
#endif
  call init_partitioning()

  do ipx1=-2,2
    do ipy1=-2,2
      do ipt1=-2,2
        do ipx2=-2,2
          do ipy2=-2,2
            do ipt2=-2,2
              if((ipx1.ne.ipx2).and.(ipy1.ne.ipy2).and.(ipt1.ne.ipt2)) then
                p1 = all_partitions(:,:,ipx1,ipy1,ipt1)
                p2 = all_partitions(:,:,ipx2,ipy2,ipt2)
                if (ip_global.eq.0) then
                  if(partition_intersect(p1,p2)) then
                    write(6,"(A25,6I4)") "Intersection not null:",ipx1,ipy1,ipt1,ipx2,ipy2,ipt2
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

  ! local volume + halo
  vol = 0
  do ipx1=-2,2
    do ipy1=-2,2
      do ipt1=-2,2
        p1 = all_partitions(:,:,ipx1,ipy1,ipt1)
        vol = vol + partition_volume(p1)
      enddo
    enddo
  enddo
  expected_vol = (ksizet_l+2)*(ksizey_l+2)*(ksizex_l+2)
  if(expected_vol.ne.vol) then
    if(ip_global.ne.0)then
      print*, "Local+Halo Volume is not correct."
    endif
  endif

  ! local volume
  vol = 0
  do ipx1=-1,1
    do ipy1=-1,1
      do ipt1=-1,1
        p1 = all_partitions(:,:,ipx1,ipy1,ipt1)
        vol = vol + partition_volume(p1)
      enddo
    enddo
  enddo
  expected_vol = ksizet_l*ksizey_l*ksizex_l
  if(expected_vol.ne.vol) then
    if(ip_global.ne.0)then
      print*, "Local Volume is not correct."
    endif
  endif

  
  expected_vol = (ksizet_l-2)*(ksizey_l-2)*(ksizex_l-2)
  if(expected_vol.ne.partition_volume(all_partitions(:,:,0,0,0))) then
    if(ip_global.ne.0)then
      print*, "Bulk Volume is not correct."
    endif
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

function partition_intersect(p1,p2) result(check)
  integer, intent(in) :: p1(2,3),p2(2,3)
  logical             :: check
  logical :: acheck1, acheck2
  integer :: idir

  check = .true.
  do idir=1,3
    acheck1 = (p1(1,idir).le.p2(2,idir)).and.(p1(2,idir).ge.p2(1,idir))
    acheck2 = (p1(2,idir).ge.p2(1,idir)).and.(p1(2,idir).le.p2(2,idir))
    check = (acheck1.or.acheck2).and.check
  enddo

end function

