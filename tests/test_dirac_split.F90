! Checking that the nonlocal parts of dirac split work as intended
! NO MPI really used here
#include "test_utils.fh"
program test_dirac_split
  use dwf3d_lib, only : init
  use mpi_f08
  use params
  use partitioning
  use dirac_split
  use dirac
  implicit none
  integer :: ierr
  integer :: ix,iy,it
  integer :: mu,iv
  integer :: ibpx,ibpy,ibpt
  integer :: tchunk(2,3)
  integer :: tchunk_s(2,3) ! tchunk Shifted
  integer :: slicedir

  integer :: idirac,igork


  !BUFfer To 
  complex(dp), target :: phi(kthird,0:ksizex_l+1,0:ksizey_l+1,0:ksizet_l+1,4)
  complex(dp), pointer :: phimod(:,:,:,:,:)
  !BUFfer To Recv
  complex(dp):: r(kthird,0:ksizex_l+1,0:ksizey_l+1,0:ksizet_l+1,4)
  complex(dp),allocatable :: res(:,:,:,:,:)
  complex(dp),allocatable :: tres(:,:,:,:,:)
  ! links
  complex(dp) :: u(0:ksizex_l+1,0:ksizey_l+1,0:ksizet_l+1,3)

#ifdef MPI
  call init_MPI
#endif
  call init(0)
 
  call init_partitioning()

  ! check u=1
  u  = 1 ! so it's just a translation
  do slicedir=1,1!3
    if(slicedir.eq.1)then
      do ix=0,ksizex_l+1
        r(:,ix,:,:,:) = ix
      enddo
    elseif(slicedir.eq.2)then
      do iy=0,ksizey_l+1
        r(:,:,iy,:,:) = iy
      enddo
    elseif(slicedir.eq.3)then
      do it=0,ksizet_l+1
        r(:,:,:,it,:) = it
      enddo
    endif
    do ibpx=-1,-1!1
      do ibpy=-1,-1!1
        do ibpt=-1,-1!1
          tchunk = border_partitions_cube(ibpx,ibpy,ibpt)%chunk
          if(ip_global.eq.0)then
            print*,tchunk
          endif
          allocate(res(kthird,&
            tchunk(1,1):tchunk(2,1),&
            tchunk(1,2):tchunk(2,2),&
            tchunk(1,3):tchunk(2,3), 4))
          allocate(tres(kthird,&
            tchunk(1,1):tchunk(2,1),&
            tchunk(1,2):tchunk(2,2),&
            tchunk(1,3):tchunk(2,3), 4))

          do mu=1,2!3
            do iv=-1,1,2
              ! shifted tchunk
              tchunk_s = tchunk
              tchunk_s(:,mu) = tchunk(:,mu) + iv
              if(ip_global.eq.0)then
                print*,"tcs",tchunk_s
              endif

              ! the modified portion of phi
              phimod => phi(:, tchunk(1,1):tchunk(2,1), &
                tchunk(1,2):tchunk(2,2), &
                tchunk(1,3):tchunk(2,3),:)

              ! DSLASH
              ! check with init = .true.
              phi= 1

              call dslash_split_nonlocal(phi,r,u,tchunk,mu,iv,.true.)

              ! computing expected result
              res = r(:, tchunk_s(1,1):tchunk_s(2,1), &
                tchunk_s(1,2):tchunk_s(2,2), &
                tchunk_s(1,3):tchunk_s(2,3),:)
              !if(ip_global.eq.0)then
              !  print*,res(1,:,:,:,:)
              !  print*,"----------------"
              !  print*,"tchunk",tchunk
              !  print*,"tchunk_s",tchunk_s
              !  print*,'r011',r(1,0,1,1,:)
              !  print*,'r211',r(1,2,1,1,:)
              !endif


              do idirac=1,4
                igork=gamin(mu,idirac)
                tres(:,:,:,:,idirac)=gamval(mu,idirac)*res(:,:,:,:,igork)
              enddo
              res = -akappa*res + tres

              if(.not.all(phimod.eq.res).and.(ip_global.eq.0))then
                write(6,"(A25,6I4)") "issue 1 with partition,dir,v",&
                &                slicedir,ibpx,ibpy,ibpt,mu,iv
                print*,"exp",res(1,:,:,:,:)
                print*,"phi",phimod(1,:,:,:,:)
              endif

              phimod = 1

              if(.not.all(phi.eq.1).and.(ip_global.eq.0))then
                 write(6,"(A25,6I4)") "issue 2 with partition,dir,v",&
                &                slicedir,ibpx,ibpy,ibpt,mu,iv
                print*,"exp",res(1,:,:,:,:)
                print*,"phi",phimod(1,:,:,:,:)
            endif

              phi= 1
              ! check with init = .false.
              call dslash_split_nonlocal(phi,r,u,tchunk,mu,iv,.false.)

              ! computing expected result
              res = res + 1 ! same as before, but added to the old value of phi

              if(.not.all(phimod.eq.res).and.(ip_global.eq.0))then
                 write(6,"(A25,6I4)") "issue 3 with partition,dir,v",&
                &                slicedir,ibpx,ibpy,ibpt,mu,iv
                print*,"exp",res(1,:,:,:,:)
                print*,"phi",phimod(1,:,:,:,:)
              endif

              phimod = 1

              if(.not.all(phi.eq.1).and.(ip_global.eq.0))then
                 write(6,"(A25,6I4)") "issue 4 with partition,dir,v",&
                &                slicedir,ibpx,ibpy,ibpt,mu,iv
                print*,"exp",res(1,:,:,:,:)
                print*,"phi",phimod(1,:,:,:,:)
              endif


              ! DSLASHD
              ! check with init = .true.
              phi= 1
              call dslashd_split_nonlocal(phi,r,u,tchunk,mu,iv,.true.)

              ! computing expected result
              res = r(:,tchunk_s(1,1):tchunk_s(2,1), &
                tchunk_s(1,2):tchunk_s(2,2), &
                tchunk_s(1,3):tchunk_s(2,3),:)
              do idirac=1,4
                igork=gamin(mu,idirac)
                tres(:,:,:,:,idirac)=-gamval(mu,idirac)*res(:,:,:,:,igork)
              enddo
              res = -akappa*res + tres


              if(.not.all(phimod.eq.res).and.(ip_global.eq.0))then
                 write(6,"(A25,6I4)") "issue 5 with partition,dir,v",&
                &                slicedir,ibpx,ibpy,ibpt,mu,iv
                print*,"exp",res(1,:,:,:,:)
                print*,"phi",phimod(1,:,:,:,:)
               endif
 
              phimod = 1

              if(.not.all(phi.eq.1).and.(ip_global.eq.0))then
                 write(6,"(A25,6I4)") "issue 6 with partition,dir,v",&
                &                slicedir,ibpx,ibpy,ibpt,mu,iv
                print*,"exp",res(1,:,:,:,:)
                print*,"phi",phimod(1,:,:,:,:)
               endif
 
              phi= 1
              ! check with init = .false.
              call dslashd_split_nonlocal(phi,r,u,tchunk,mu,iv,.false.)

              ! computing expected result
              res = res + 1 ! same as before, but added to the old value of phi

              if(.not.all(phimod.eq.res).and.(ip_global.eq.0))then
                 write(6,"(A25,6I4)") "issue 7 with partition,dir,v",&
                &                slicedir,ibpx,ibpy,ibpt,mu,iv
                print*,"exp",res(1,:,:,:,:)
                print*,"phi",phimod(1,:,:,:,:)
               endif
 
              phi(:, tchunk(1,1):tchunk(2,1), &
                tchunk(1,2):tchunk(2,2), &
                tchunk(1,3):tchunk(2,3),:) = 1

              if(.not.all(phi.eq.1).and.(ip_global.eq.0))then
                 write(6,"(A25,6I4)") "issue 8 with partition,dir,v",&
                &                slicedir,ibpx,ibpy,ibpt,mu,iv
                print*,"exp",res(1,:,:,:,:)
                print*,"phi",phimod(1,:,:,:,:)
               endif
            enddo
          enddo
          deallocate(tres) 
          deallocate(res) 
        enddo
      enddo
    enddo
  enddo
#ifdef MPI
  call MPI_Finalize(ierr)
#endif

end program


