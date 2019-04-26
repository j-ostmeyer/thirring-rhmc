! Testing that each rank is referenced the correct number of times
! in the global structures in the partitioning module
#include "test_utils.fh"
program test_partitioning
  use mpi
  use params
  use partitioning
  use comms
  implicit none
  integer :: ierr

  integer,allocatable :: mpi_other_rank_count_cube(:)
  integer,allocatable :: mpi_other_rank_count_list(:)
  integer,allocatable :: expected_mpi_other_rank_count(:)

  integer :: ipxs(2,3)
  integer :: idir
  integer :: ipx,ipy,ipt
  integer :: ihp

  integer :: nn,nns(3)
  integer :: ipart
  integer :: discr

  integer :: irank


#ifdef MPI
  call init_MPI
#endif
  call init_partitioning()

  allocate(mpi_other_rank_count_cube(0:np_global-1))
  allocate(mpi_other_rank_count_list(0:np_global-1))
  allocate(expected_mpi_other_rank_count(0:np_global-1))

  expected_mpi_other_rank_count = 0
  do idir=1,3
    call MPI_Cart_Shift(comm, idir-1, 1, ipxs(1,idir), ipxs(2,idir),ierr)
    expected_mpi_other_rank_count(ipxs(1,idir)) = 9 + expected_mpi_other_rank_count(ipxs(1,idir))
    expected_mpi_other_rank_count(ipxs(2,idir)) = 9 + expected_mpi_other_rank_count(ipxs(2,idir))
  enddo

  ! checking the cube version - border
  mpi_other_rank_count_cube = 0
  do ipx=-1,1
    do ipy=-1,1
      do ipt=-1,1
        if((ipx**2+ipy**2+ipt**2).ne.0) then
          nn = border_partitions_cube(ipx,ipy,ipt)%nn
          nns = border_partitions_cube(ipx,ipy,ipt)%nns
          do idir=1,nn
            mpi_other_rank_count_cube(nns(idir)) = mpi_other_rank_count_cube(nns(idir))+1
          enddo
        endif
      enddo
    enddo
  enddo

  do irank=0,np_global-1
    call MPI_Barrier(comm,ierr)
    if(irank.eq.ip_global) then
      if(sum((mpi_other_rank_count_cube-expected_mpi_other_rank_count)**2).ne.0)then
        print *, "Issue with mpi_other_rank_count_cube - border - rank",ip_global
        print *,mpi_other_rank_count_cube
        print *,expected_mpi_other_rank_count
      endif 
    endif
  enddo


  ! checking the list version - border
  mpi_other_rank_count_list = 0
  do ipart=1,26
    nn = border_partitions_list(ipart)%nn
    nns = border_partitions_list(ipart)%nns
    do idir=1,nn
      mpi_other_rank_count_list(nns(idir)) = mpi_other_rank_count_list(nns(idir))+1
    enddo
  enddo

  do irank=0,np_global-1
    call MPI_Barrier(comm,ierr)
    if(irank.eq.ip_global) then
      if(sum((mpi_other_rank_count_list-expected_mpi_other_rank_count)**2).ne.0)then
        print *, "Issue with mpi_other_rank_count_list - border - rank", ip_global
        print *,mpi_other_rank_count_list
        print *,expected_mpi_other_rank_count
      endif 
    endif
  enddo

  ! checking the cube version - halo
  mpi_other_rank_count_cube = 0
  do ipx=-2,2
    do ipy=-2,2
      do ipt=-2,2
        discr = ipx**2+ipy**2+ipt**2
        if((discr.ge.4).and.(discr.le.6)) then
          nn = halo_partitions_cube(ipx,ipy,ipt)%nn
          mpi_other_rank_count_cube(nn) = mpi_other_rank_count_cube(nn)+1
        endif
      enddo
    enddo
  enddo


  do irank=0,np_global-1
    call MPI_Barrier(comm,ierr)
    if(irank.eq.ip_global) then
      if(sum((mpi_other_rank_count_cube-expected_mpi_other_rank_count)**2).ne.0)then
        print *, "Issue with mpi_other_rank_count_cube - halo - rank", ip_global
        print *,mpi_other_rank_count_cube
        print *,expected_mpi_other_rank_count
      endif 
    endif 
  enddo


  ! checking the list version - halo
  mpi_other_rank_count_list = 0
  do ipart=1,54
    nn = halo_partitions_list(ipart)%nn
    mpi_other_rank_count_list(nn) = mpi_other_rank_count_list(nn)+1
  enddo

  do irank=0,np_global-1
    call MPI_Barrier(comm,ierr)
    if(irank.eq.ip_global) then
      if(sum((mpi_other_rank_count_list-expected_mpi_other_rank_count)**2).ne.0)then
        print *, "Issue with mpi_other_rank_count_list - halo - rank", ip_global
        print *,mpi_other_rank_count_list
        print *,expected_mpi_other_rank_count
      endif 
    endif 
  enddo


  ! checking ahpsr member in border partitions.

  ! checking the cube version - border
  mpi_other_rank_count_cube = 0
  do ipx=-1,1
    do ipy=-1,1
      do ipt=-1,1
        do idir=-3,3
            ihp = border_partitions_cube(ipx,ipy,ipt)%ahpsr(idir)
            if(ihp.ne.0) then
              nn = halo_partitions_list(ihp)%nn
              mpi_other_rank_count_cube(nn) = mpi_other_rank_count_cube(nn)+1
            endif
        enddo
      enddo
    enddo
  enddo

  do irank=0,np_global-1
    call MPI_Barrier(comm,ierr)
    if(irank.eq.ip_global) then
      if(sum((mpi_other_rank_count_cube-expected_mpi_other_rank_count)**2).ne.0)then
        print *, "Issue with mpi_other_rank_count_cube - border - rank",ip_global
        print *,mpi_other_rank_count_cube
        print *,expected_mpi_other_rank_count
      endif 
    endif
  enddo

  ! checking the list version - border
  mpi_other_rank_count_list = 0
  do ipart=1,26
    do idir=-3,3
      ihp = border_partitions_list(ipart)%ahpsr(idir)
      if(ihp.ne.0) then
        nn = halo_partitions_list(ihp)%nn
        mpi_other_rank_count_list(nn) = mpi_other_rank_count_list(nn)+1
      endif
    enddo
  enddo

  do irank=0,np_global-1
    call MPI_Barrier(comm,ierr)
    if(irank.eq.ip_global) then
      if(sum((mpi_other_rank_count_list-expected_mpi_other_rank_count)**2).ne.0)then
        print *, "Issue with mpi_other_rank_count_list - border - rank (ahpsr)", ip_global
        print *,mpi_other_rank_count_list
        print *,expected_mpi_other_rank_count
      endif 
    endif
  enddo

 
  deallocate(mpi_other_rank_count_cube)
  deallocate(mpi_other_rank_count_list)
  deallocate(expected_mpi_other_rank_count)



#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program

