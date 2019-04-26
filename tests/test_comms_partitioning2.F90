! Testing communications with the partitions on 2 dummy arrays.
#include "test_utils.fh"
program test_partitioning
  use params
  use partitioning
  use comms_partitioning
  use comms
  implicit none
  !BUFfer To Send
  complex(dp) :: bufts(kthird,0:ksizex_l+1,0:ksizey_l+1,0:ksizet_l+1,4)
  !BUFfer To Recv
  complex(dp) :: buftr(kthird,0:ksizex_l+1,0:ksizey_l+1,0:ksizet_l+1,4)

  integer :: rreqs(54)! Recv REQuestS
  integer :: sreqs(54)! Send REQuestS

  integer :: ix0,iy0,it0
  integer :: ix,iy,it

  integer :: ierr
  logical :: check_cube_eq

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
    bufts(:,ix,:,:,:) = modulo(ix0+ix,KSIZE)
  enddo
  ! taking care of the local part
  buftr(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
    & bufts(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)

  ! taking care of the halos
  call MPI_StartAll(54,rreqs,ierr)
  call MPI_StartAll(54,sreqs,ierr)

  call MPI_Waitall(54,rreqs,MPI_Statuses_Ignore,ierr)
  call MPI_Waitall(54,sreqs,MPI_Statuses_Ignore,ierr)

  if(.not.check_cube_eq(buftr,bufts))then
    print *,"X dir:problem on rank",ip_x,ip_y,ip_t,sum((buftr-bufts)**2)
    print*,int(real(buftr(1,:,1,1,1))),ip_x
  endif


  ! test y-dir
  ! filling send buffer
  do iy=0,ksizey_l+1
    bufts(:,:,iy,:,:) = modulo(iy0+iy,KSIZE)
  enddo
  ! taking care of the local part
  buftr(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
    & bufts(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)

  ! taking care of the halos
  call MPI_StartAll(54,rreqs,ierr)
  call MPI_StartAll(54,sreqs,ierr)

  call MPI_Waitall(54,rreqs,MPI_Statuses_Ignore,ierr)
  call MPI_Waitall(54,sreqs,MPI_Statuses_Ignore,ierr)

  if(.not.check_cube_eq(buftr,bufts))then
    print *,"Y dir:problem on rank",ip_x,ip_y,ip_t,sum((buftr-bufts)**2)
  endif

  ! test t-dir
  ! filling send buffer
  do it=0,ksizet_l+1
    bufts(:,:,:,it,:) = modulo(it0+it,KSIZE)
  enddo
  ! taking care of the local part
  buftr(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
    & bufts(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)

  ! taking care of the halos
  call MPI_StartAll(54,rreqs,ierr)
  call MPI_StartAll(54,sreqs,ierr)

  call MPI_Waitall(54,rreqs,MPI_Statuses_Ignore,ierr)
  call MPI_Waitall(54,sreqs,MPI_Statuses_Ignore,ierr)

  if(.not.check_cube_eq(buftr,bufts))then
    print *,"T dir:problem on rank",ip_x,ip_y,ip_t,sum((buftr-bufts)**2)
  endif

#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program

function check_cube_eq(buf1,buf2) result(check)
  use params
  complex(dp),intent(in) :: buf1(kthird,0:ksizex_l+1,0:ksizey_l+1,0:ksizet_l+1,4)
  complex(dp),intent(in) :: buf2(kthird,0:ksizex_l+1,0:ksizey_l+1,0:ksizet_l+1,4)

  logical :: check

  check = all(buf2(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:).eq.&
    & buf1(:,1:ksizex_l,1:ksizey_l,1:ksizet_l,:))

  ! down x
  check = check.and.all(buf2(:,0,1:ksizey_l,1:ksizet_l,:).eq.&
    & buf1(:,0,1:ksizey_l,1:ksizet_l,:))
  !up x
  check = check.and.all(buf2(:,1+ksizex_l,1:ksizey_l,1:ksizet_l,:).eq.&
    & buf1(:,1+ksizex_l,1:ksizey_l,1:ksizet_l,:))

  !down y
  check = check.and.all(buf2(:,1:ksizex_l,0,1:ksizet_l,:).eq.&
    & buf1(:,1:ksizex_l,0,1:ksizet_l,:))
  !up y
  check = check.and.all(buf2(:,1:ksizex_l,1+ksizey_l,1:ksizet_l,:).eq.&
    & buf1(:,1:ksizex_l,1+ksizey_l,1:ksizet_l,:))

  !down t
  check = check.and.all(buf2(:,1:ksizex_l,1:ksizey_l,0,:).eq.&
    & buf1(:,1:ksizex_l,1:ksizey_l,0,:))
  !up t
  check = check.and.all(buf2(:,1:ksizex_l,1:ksizey_l,1+ksizet_l,:).eq.&
    & buf1(:,1:ksizex_l,1:ksizey_l,1+ksizet_l,:))

end function
