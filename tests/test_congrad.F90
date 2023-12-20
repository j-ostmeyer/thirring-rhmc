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

  integer :: ierr, imass, itercg
  real :: res, am
  complex(dp) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) :: xToWrite(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: x_ref(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: diff(kthird_l, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: sum_diff
  real(dp) :: max_diff
  integer, dimension(16) :: reqs_X, reqs_Phi

#ifdef MPI
  call init_MPI
#endif

  ! Setup state
  res = 0.1
  am = 0.05
  imass = 3
  call generate_starting_state_Phi_and_X(Phi, reqs_Phi, u, X, reqs_X)

  ! call function
  call congrad(Phi, res, itercg, am, imass)

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

#ifdef MPI
  call MPI_Finalize(ierr)
#endif

end program
