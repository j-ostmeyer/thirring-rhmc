#include "test_utils.fh"
program test_dslash
  use params
  use mpi
  ! use dwf3d_lib
  use diracWilson
  use comms
  use comms4
  use comms5
  use test_utils

  implicit none

  ! general parameters
  logical :: generate = .false.
  integer :: i, ierr, imass_index, imass, timing_loops = 1
  integer, dimension(3) :: imasses = (/1,3,5/)
  character(len=4) :: imass_char
  character(len=*), parameter :: test_prefix = 'test_dslash_'
  integer, parameter :: kthird_l = kthird ! HACK TODO remove me when migrating to master
  integer, parameter :: ip_third = 0 ! HACK TODO remove me when migrating to master

  ! initialise function parameters
  complex(dp) u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
  complex(dp) Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  integer, dimension(16) :: reqs_Phi

#ifdef MPI
  call init_MPI
#endif 

  do imass_index = 1, size(imasses)
    imass_char = ''
    imass = imasses(imass_index)
    write(imass_char, '(I1)') imass
    if (ip_global == 0) then
      print *, ' imass: ', imass_char
    end if
#ifdef MPI
    call MPI_Barrier(comm, ierr)
#endif
    call generate_starting_state(Phi, R, u, reqs_Phi)

    call run_dslash(Phi, R, u, imass, timing_loops, reqs_Phi)
    if (generate) then
      call generate_data(Phi, test_prefix // trim(imass_char))
    else
      call validate_results(Phi, test_prefix // trim(imass_char))
    end if
  end do

#ifdef MPI
  call MPI_Finalize(ierr)
#endif

contains
  subroutine generate_starting_state(Phi, R, u, reqs_Phi)
    complex(dp) Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    integer, dimension(16), intent(inout) :: reqs_Phi

    complex, parameter :: iunit = cmplx(0, 1)
    real(dp), parameter :: tau = 8*atan(1.0_8)

    integer :: i, j, ix, iy, it, ithird
    integer, parameter :: idxmax = 4*ksize*ksize*ksizet*kthird
    integer :: idx
#ifdef MPI
    integer, dimension(16) :: reqs_R
    integer, dimension(12) :: reqs_u
    integer :: ierr
#endif
    do j = 1, 4
      do it = 1, ksizet_l
        do iy = 1, ksizey_l
          do ix = 1, ksizex_l
            do ithird = 1, kthird_l
              idx = ip_third*kthird_l + ithird &
                    + (ip_x*ksizex_l + ix - 1)*kthird &
                    + (ip_y*ksizey_l + iy - 1)*kthird*ksize &
                    + (ip_t*ksizet_l + it - 1)*kthird*ksize*ksize &
                    + (j - 1)*kthird*ksize*ksize*ksizet

              Phi(ithird, ix, iy, it, j) = 1.1*exp(iunit*idx*tau/idxmax)
              R(ithird, ix, iy, it, j) = 1.3*exp(iunit*idx*tau/idxmax)
            enddo
          enddo
        enddo
      enddo
    enddo
#ifdef MPI
    call start_halo_update_5(4, R, 0, reqs_R)
    call start_halo_update_5(4, Phi, 1, reqs_Phi)
#endif
      do j = 1, 3
        do it = 1, ksizet_l
          do iy = 1, ksizey_l
            do ix = 1, ksizex_l
              idx = ip_x*ksizex_l + ix &
                    + (ip_y*ksizey_l + iy - 1)*ksize &
                    + (ip_t*ksizet_l + it - 1)*ksize*ksize &
                    + (j - 1)*ksize*ksize*ksizet

              u(ix, iy, it, j) = exp(iunit*idx*tau/idxmax)
            enddo
          enddo
        enddo
      enddo
#ifdef MPI
    call start_halo_update_4(3, u, 1, reqs_u)
    call complete_halo_update(reqs_R)
    call complete_halo_update(reqs_Phi)
    call MPI_WaitAll(12, reqs_u, MPI_STATUSES_IGNORE, ierr)
#else
    call update_halo_5(4, R)
    call update_halo_5(4, Phi)
    call update_halo_4(3, u)
#endif

    ! initialise common variables
    beta = 0.4
    am3 = 1.0
    ibound = -1

    call init_gammas()
  end subroutine generate_starting_state

  subroutine run_dslash(Phi, R, u, imass, timing_loops, reqs_Phi)
    complex(dp) Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    integer, intent(in) :: timing_loops
    integer, intent(in) :: imass
    integer, dimension(16), intent(inout) :: reqs_Phi
    
    real, parameter :: am = 0.05

    ! call function
    do i = 1, timing_loops
      call dslash(Phi, R, u, am, imass)
#ifdef MPI
      call start_halo_update_5(4, Phi, 2, reqs_Phi)
      call complete_halo_update(reqs_Phi)
#else
      call update_halo_5(4, Phi)
#endif
    end do
  end subroutine run_dslash

  subroutine generate_data(Phi, test_name)
    complex(dp) Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) PhiToWrite(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
    character(len=*), intent(in) :: test_name

    PhiToWrite(:, :, :, :, :) = Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
    write_file(PhiToWrite, trim(test_name)//'.dat', MPI_Double_Complex)
  end subroutine generate_data

  subroutine validate_results(Phi, test_name)
    complex(dp) Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) Phiref(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
    complex(dp) diff(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
    complex(dp) sum_diff
    real(dp) max_diff
    character(len=*), intent(in) :: test_name

    read_file(Phiref, trim(test_name)//'.dat', MPI_Double_Complex)

    diff = Phi(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) - Phiref
    
    check_max(diff, 1e-11, 'Phi', max_diff, MPI_Double_Precision, trim(test_name))
    check_sum(diff, 1e-11, 'Phi', sum_diff, MPI_Double_Complex, trim(test_name))
  end subroutine validate_results
end program
