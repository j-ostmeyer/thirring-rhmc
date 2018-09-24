! Checking that the nonlocal parts of dirac split work as intended
! NO MPI used here
#include "test_utils.fh"
program test_dirac_split
  use params
  use partitioning
  use dirac_split
  use dirac
  implicit none
  integer :: ierr
  integer :: ix,iy,it
  integer :: mu,v
  integer :: ibpx,ibpy,ipbt
  integer :: tchunk(2,3)
  integer :: tchunk_s(2,3) ! tchunk Shifted
  integer :: slicedir

  integer :: idirac,igork


  !BUFfer To 
  complex(dp) :: phi(kthird,0:ksizex_l+1,0:ksizey_l+1,0:ksizet_l+1,4)
  complex(dp), pointer :: phimod(:,:,:,:,:)
  !BUFfer To Recv
  complex(dp) :: rs(kthird,0:ksizex_l+1,0:ksizey_l+1,0:ksizet_l+1,4,3)
  complex(dp),pointer :: r(:,:,:,:,:)
  complex(dp),allocatable :: res(:,:,:,:,:)
  complex(dp),allocatable :: tres(:,:,:,:,:)
  ! links
  complex(dp) :: us(0:ksizex_l+1,0:ksizey_l+1,0:ksizet_l+1,3,3)
  complex(dp),pointer :: u(:,:,:,:)

  
  call init_partitioning()

  ! check u=1
  us  = 1 ! so it's just a translation
  u => us(1)
  rs  = 0
  do ix=1,ksizex_l
    rs(:,ix,:,:,:,1) = ix
  enddo
  do iy=1,ksizey_l
    rs(:,:,iy,:,:,2) = iy
  enddo
  do it=1,ksizet_l
    rs(:,:,:,it,:,3) = it
  enddo

  do slicedir=1,3
    r => rs(:,:,:,:,:,slicedir)
    do ibpx=-1,1
      do ibpy=-1,1
        do ibpt=-1,1
          tchunk = border_partitions_cube(ibpx,ibpy,ipbt)
         allocate(res(kthird,&
            tchunk(1,1):tchunk(2,1),&
            tchunk(1,2):tchunk(2,2),&
            tchunk(1,3):tchunk(2,3), 4))
          allocate(tres(kthird,&
            tchunk(1,1):tchunk(2,1),&
            tchunk(1,2):tchunk(2,2),&
            tchunk(1,3):tchunk(2,3), 4))

          do mu=1,3
            do iv=1,2
              ! shifted tchunk
              tchunk_s = tchunk
              tchunk_s(:,mu) = tchunk(:,mu) + (iv*2-3)
 
              ! DSLASH
              ! check with init = .true.
              phi= 1
              call dslash_split_nonlocal(phi,r,u,tchunk,mu,v,.true.)

              ! computing expected result
              res = r(:, tchunk(1,1):tchunk(2,1), &
                tchunk(1,2):tchunk(2,2), &
                tchunk(1,3):tchunk(2,3),:)
              do idirac=1,4
                igork=gamin(mu,idirac)
                tres(:,:,:,:,idirac)=gamval(mu,idirac)*res(:,:,:,:,igork)
              enddo
              res = akappa*res + tres

              ! the modified portion of phi
              phimod => phi(:, tchunk_s(1,1):tchunk_s(2,1), &
                tchunk_s(1,2):tchunk_s(2,2), &
                tchunk_s(1,3):tchunk_s(2,3),:)

              if(.not.(phimod.eq.res))then
                print*, "d - issue 1 with partition,dir,v", ibpx,ibpy,ipbt,mu,v
              endif

              phimod = 1

              if(.not.all(phi.eq.1))then
                print*, "d - issue 2 with partition,dir,v", ibpx,ibpy,ipbt,mu,v
              endif

              phi= 1
              ! check with init = .false.
              call dslash_split_nonlocal(phi,r,u,tchunk,mu,v,.false.)

              ! computing expected result
              res = res + 1 ! same as before, but added to the old value of phi

              if(.not.(phimod.eq.res))then
                print*, "d - issue 3 with partition,dir,v", ibpx,ibpy,ipbt,mu,v
              endif

              phimod = 1

              if(.not.all(phi.eq.1))then
                print*, "d - issue 4 with partition,dir,v", ibpx,ibpy,ipbt,mu,v
              endif

              ! DSLASHD
              ! check with init = .true.
              phi= 1
              call dslashd_split_nonlocal(phi,r,u,tchunk,mu,v,.true.)

              ! computing expected result
              res = r(:,tchunk(1,1):tchunk(2,1), &
                tchunk(1,2):tchunk(2,2), &
                tchunk(1,3):tchunk(2,3),:)
              do idirac=1,4
                igork=gamin(mu,idirac)
                tres(:,:,:,:,idirac)=-gamval(mu,idirac)*res(:,:,:,:,igork)
              enddo
              res = -akappa*res + tres


              if(.not.(phimod.eq.res))then
                print*, "d - issue 1 with partition,dir,v", ibpx,ibpy,ipbt,mu,v
              endif

              phimod = 1

              if(.not.all(phi.eq.1))then
                print*, "d - issue 2 with partition,dir,v", ibpx,ibpy,ipbt,mu,v
              endif

              phi= 1
              ! check with init = .false.
              call dslashd_split_nonlocal(phi,r,u,tchunk,mu,v,.false.)

              ! computing expected result
              res = res + 1 ! same as before, but added to the old value of phi

              if(.not.(phimod.eq.res))then
                print*, "d - issue 3 with partition,dir,v", ibpx,ibpy,ipbt,mu,v
              endif

              phi(:, tchunk_s(1,1):tchunk_s(2,1), &
                tchunk_s(1,2):tchunk_s(2,2), &
                tchunk_s(1,3):tchunk_s(2,3),:) = 1

              if(.not.all(phi.eq.1))then
                print*, "d - issue 4 with partition,dir,v", ibpx,ibpy,ipbt,mu,v
              endif
            enddo
          enddo
          deallocate(tres) 
          deallocate(res) 
        enddo
      enddo
    enddo
  enddo
end program


subroutine gamult(output,input,mu,)
  use dirac
  use params
  complex(dp), intent(in)  ::  input(kthird,0:ksizex_l+1,0:ksizey_l:1,0:ksizet_l+1,4)
  complex(dp), intent(out) :: output(kthird,0:ksizex_l+1,0:ksizey_l:1,0:ksizet_l+1,4)
  integer, intent(in) :: mu

  integer ::idirak,igork

  do idirac=1,4
    igork=gamin(mu,idirac)
    output(:,:,:,:,idirac)= gamval(mu,idirac) * input(:, :, :, :, igork)
  enddo

  end subroutine

