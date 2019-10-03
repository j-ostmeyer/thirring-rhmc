pure function pid(fip_x, fip_y, fip_t) result(id)
  use params
  implicit none
  integer, intent(in) :: fip_x, fip_y, fip_t
  integer :: id

  id = fip_x + fip_y * np_x + fip_t * np_x * np_y
end function pid

program test_halo_4
  use params
  use comms
  use comms4
  implicit none
  integer :: ix, iy, it, i5, i=0
  integer :: pid
  integer :: passed_basic = 0
  complex(dp) :: test_array(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
#ifdef MPI
  integer :: reqs(12)
  integer :: ierr

  call init_MPI
#endif

  test_array(0,:,:,:) = cmplx(-1,-1)
  test_array(:,0,:,:) = cmplx(-2,-1)
  test_array(:,:,0,:) = cmplx(-3,-1)
  test_array(ksizex_l+1,:,:,:) = cmplx(-1,-2)
  test_array(:,ksizey_l+1,:,:) = cmplx(-2,-2)
  test_array(:,:,ksizet_l+1,:) = cmplx(-3,-2)

  ! Set up local arrays
  do i5=1,4
    do it=1,ksizet_l
      do iy=1,ksizey_l
        do ix=1,ksizex_l
          i = i + 1
          test_array(ix, iy, it, i5) = cmplx(i, pid(ip_x, ip_y, ip_t))
        end do
      end do
    end do
  end do

  ! Communicate
#ifdef MPI
  call start_halo_update_4(4, test_array, 0, reqs)
  call complete_halo_update(reqs)
#else
  call update_halo_4(4, test_array)
#endif

  ! Check output
  if (real(test_array(1,1,1,1)) .ne. real(test_array(ksizex_l+1,1,1,1)) .or. &
    nint(aimag(test_array(ksizex_l+1,1,1,1))) &
    .ne. pid(modulo(ip_x+1, np_x), ip_y, ip_t)) then
    print *, "Negative x update failed on process", ip_x, ip_y, ip_t
  end if
  if (real(test_array(1,1,1,1)) .ne. real(test_array(1,ksizey_l+1,1,1)) .or. &
    nint(aimag(test_array(1,ksizey_l+1,1,1))) &
    .ne. pid(ip_x, modulo(ip_y+1, np_y), ip_t)) then
    print *, "Negative y update failed on process", ip_x, ip_y, ip_t
  end if
  if (real(test_array(1,1,1,1)) .ne. real(test_array(1,1,ksizet_l+1,1)) .or. &
    nint(aimag(test_array(1,1,ksizet_l+1,1))) &
    .ne. pid(ip_x, ip_y, modulo(ip_t+1, np_t))) then
    print *, "Negative t update failed on process", ip_x, ip_y, ip_t
  end if
  if (real(test_array(ksizex_l,1,1,1)) .ne. real(test_array(0,1,1,1)) .or. &
    nint(aimag(test_array(0,1,1,1))) &
    .ne. pid(modulo(ip_x-1, np_x), ip_y, ip_t)) then
    print *, "Positive x update failed on process", ip_x, ip_y, ip_t
  end if
  if (real(test_array(1,ksizey_l,1,1)) .ne. real(test_array(1,0,1,1)) .or. &
    nint(aimag(test_array(1,0,1,1))) &
    .ne. pid(ip_x, modulo(ip_y-1, np_y), ip_t)) then
    print *, "Positive y update failed on process", ip_x, ip_y, ip_t
  endif
  if (real(test_array(1,1,ksizet_l,1)) .ne. real(test_array(1,1,0,1)) .or. &
    nint(aimag(test_array(1,1,0,1))) &
    .ne. pid(ip_x, ip_y,modulo(ip_t-1, np_t))) then
    print *, "Positive t update failed on process", ip_x, ip_y, ip_t
  else
    passed_basic = 1
  end if
  if (passed_basic.eq.1 .and. &
    (real(test_array(1,1,ksizet_l,4)) .ne. real(test_array(1,1,0,4)) .or. &
    nint(aimag(test_array(1,1,0,4))) &
    .ne. pid(ip_x, ip_y, modulo(ip_t-1, np_t)))) then
    print *, "Extent is wrong in size4 direction; process", ip_x, ip_y, ip_t
  end if


#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program test_halo_4
