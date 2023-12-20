#include "test_utils.fh"
program test_dslashd_reqs_wilson
  use params
  use mpi
  ! use dwf3d_lib
  use dirac
  use comms
  use comms4
  use comms5
  use test_utils

  implicit none

  integer :: ierr, imass_index, imass
  character(len=4) :: imass_char
  character(len=*), parameter :: test_prefix = 'test_dslashd_reqs_wilson_'
  complex(dp) u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
  complex(dp) Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  integer, dimension(16) :: reqs_Phi, reqs_R

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
    call generate_starting_state(Phi, reqs_Phi, u, R, reqs_R)

    call run_dslashd(Phi, R, u, imass, reqs_Phi, reqs_R)
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
  subroutine run_dslashd(Phi, R, u, imass, reqs_Phi, reqs_R)
    complex(dp) Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    integer, intent(in) :: imass
    integer, dimension(16), intent(inout) :: reqs_Phi, reqs_R
    
    real, parameter :: am = 0.05

    ! call function
#ifdef MPI
    call dslashd_wilson(Phi, R, u, am, imass, reqs_R)
    call start_halo_update_5(4, Phi, 2, reqs_Phi)
    call complete_halo_update(reqs_Phi)
#else
    print *, 'Error!! this test requires MPI'
    call exit(-1)
#endif
  end subroutine run_dslashd

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
