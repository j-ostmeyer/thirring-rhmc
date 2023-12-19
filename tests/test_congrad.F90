#include "test_utils.fh"
program test_congrad
  use dwf3d_lib
  use trial, only: u
  use vector
  use dirac
  use comms
  use comms4
  use comms5
  use gammamatrices
  use measure_module
  use test_utils
  implicit none

  ! general parameters
  integer :: i, ierr, iflag, imass, isweep, timing_loops = 1
  complex, parameter :: iunit = cmplx(0, 1)
  real*8, parameter :: tau = 8*atan(1.0_8)

  ! initialise function parameters
  complex(dp) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: xToWrite(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: x_ref(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: diff(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: sum_diff
  real(dp) :: max_diff

  real :: res, am
  integer :: itercg

#ifdef MPI
  integer, dimension(16) :: reqs_X, reqs_Phi

  call init_MPI
#endif

  res = 0.1
  am = 0.05
  imass = 3
  iflag = 0
  isweep = 1

  call generate_starting_state_Phi_and_X(Phi, reqs_Phi, u, X, reqs_X)

  ! call function
  do i = 1, timing_loops
    call congrad(Phi, res, itercg, am, imass)
  end do

  ! check output
  if (generate) then
    xToWrite(:, :, :, :, :) = x(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
    write_file(xToWrite, 'test_congrad.dat', MPI_Double_Complex)
  else
    read_file(x_ref, 'test_congrad.dat', MPI_Double_Complex)

    diff = x_ref - x(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
    
    check_equality(itercg, 27, 'itercg', "test_congrad")
    check_sum(diff, 2, 'x', sum_diff, MPI_Double_Complex, "test_congrad")
    check_max(diff, 5e-2, 'x', max_diff, MPI_Double_Precision, "test_congrad")
  end if

  call MPI_Finalize(ierr)

end program
