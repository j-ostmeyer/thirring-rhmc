! Testing communications with the partitions on 2 dummy arrays.
#include "test_utils.fh"
program test_partitioning
  use mpi_f08
  use params
  use partitioning
  use comms_partitioning
  use comms
  implicit none
  !BUFfer To Send
  complex(dp) :: bufts(kthird,0:ksizex_l+1,0:ksizey_l+1,0:ksizet_l+1,4)
  !BUFfer To Recv
  complex(dp) :: buftr(kthird,0:ksizex_l+1,0:ksizey_l+1,0:ksizet_l+1,4)


  type(MPI_Request) :: rreqs(54)! Recv REQuestS
  type(MPI_Request) :: sreqs(54)! Send REQuestS

  integer :: ix0,iy0,it0
  integer :: ix,iy,it

  integer :: ierr

#ifdef MPI
  call init_MPI
#endif
  call init_partitioning()
  call init_dirac_hb_types()

  call get_dirac_sendreqs(sreqs,bufts)
  call get_dirac_recvreqs(rreqs,buftr)

  ix0 = ip_x*ksizex_l
  iy0 = ip_y*ksizey_l
  it0 = ip_t*ksizet_l


  ! test x-dir
  ! filling send buffer
  do ix=0,ksizex_l+1
    bufts(:,ix,:,:,:) = ix0+ix
  enddo
  ! taking care of the local part
  buftr(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
    & bufts(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)

  ! taking care of the halos
  MPI_StartAll(54,rreqs,ierr)
  MPI_StartAll(54,sreqs,ierr)

  MPI_Waitall(54,rreqs,MPI_Statuses_Ignore,ierr)

  if(.not.all(buftr.eq.bufts))then
    print *,"X dir:problem on rank",ip_x,ip_y,ip_t
  endif

  ! test y-dir
  ! filling send buffer
  do iy=0,ksizey_l+1
    bufts(:,:,iy,:,:) = iy0+iy
  enddo
  ! taking care of the local part
  buftr(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
    & bufts(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)

  ! taking care of the halos
  MPI_StartAll(54,rreqs,ierr)
  MPI_StartAll(54,sreqs,ierr)

  MPI_Waitall(54,rreqs,MPI_Statuses_Ignore,ierr)

  if(.not.all(buftr.eq.bufts))then
    print *,"Y dir:problem on rank",ip_x,ip_y,ip_t
  endif

  ! test t-dir
  ! filling send buffer
  do it=0,ksizet_l+1
    bufts(:,:,:,it,:) = it0+it
  enddo
  ! taking care of the local part
  buftr(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
    & bufts(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)

  ! taking care of the halos
  MPI_StartAll(54,rreqs,ierr)
  MPI_StartAll(54,sreqs,ierr)

  MPI_Waitall(54,rreqs,MPI_Statuses_Ignore,ierr)

  if(.not.all(buftr.eq.bufts))then
    print *,"T dir:problem on rank",ip_x,ip_y,ip_t
  endif

#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program

