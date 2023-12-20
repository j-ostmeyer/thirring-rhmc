#include "test_utils.fh"

program test_derivs_shamir
  use dirac
  use gforce
  use comms
  use comms4
  use comms5
  use test_utils
  use derivs_module, only: derivs_shamir
  implicit none

  integer :: ierr, imass
  real(dp), parameter :: tau = 8*atan(1.0_8)
  real(dp), parameter :: anum = tau
  real :: am
  integer, parameter :: iflag = 0
  integer, dimension(16) :: reqs_R, reqs_Phi, reqs_X2
  complex(dp) u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
  complex(dp) Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  complex(dp) X2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  real(dp) :: dSdpi_orig(ksizex_l, ksizey_l, ksizet_l, 3)
  real(dp) :: dSdpi_ref(ksizex_l, ksizey_l, ksizet_l, 3)
  complex(dp) R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  real(dp) :: diff(ksizex_l, ksizey_l, ksizet_l, 3), sum_diff, max_diff

#ifdef MPI
  call init_MPI
#endif

  ! Setup state
  am = 0.05
  imass = 3
  call generate_starting_state(Phi, reqs_Phi, u, R, reqs_R, X2, reqs_X2, dSdpi_orig)
  dSdpi = dSdpi_orig

  ! call function
  call derivs_shamir(R, X2, anum, iflag, am, imass)

  ! check output
  if (generate) then
    print *, "Generating .dat file..."
    write_file(dSdpi, 'test_derivs_shamir.dat', MPI_Double_Precision)
  else
    read_file(dSdpi_ref, 'test_derivs_shamir.dat', MPI_Double_Precision)

    diff = dSdpi - dSdpi_ref

#ifdef MPI
    ! diff will now have duplicates for the same (x,y,t,mu)
    ! due to the parallelization along third dimension
    ! and the MPI_AllReduce operation. For this the sum_diff
    ! will be larger than expected (np_third times larger).
    ! Because of that we divide by np_third
    diff = diff / np_third
#endif

    check_sum(diff, 0.3, 'dSdpi', sum_diff, MPI_Double_Precision, "test_derivs_shamir")
    check_max(diff, 0.01, 'dSdpi', max_diff, MPI_Double_Precision, "test_derivs_shamir")

  end if

#ifdef MPI
  call MPI_Finalize(ierr)
#endif

end program test_derivs_shamir
