module inverter_utils
  use params
  implicit none
  complex(dp) :: xout(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

contains

  !   For the multi-shift inverter
  subroutine dirac_op_shifted(xout, xin, u, am, imass, shift)
    complex(dp), intent(out) :: xout(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: xin(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    real, intent(in) :: am, shift
    integer, intent(in) :: imass

    call dirac_operator(xout, xin, u, am, imass)
    xout = (xout + shift*xin)
  end subroutine

  subroutine dirac_operator(xout, xin, u, am, imass)
    use dirac, only: dslash, dslashd
    use comms
#ifdef MPI
    use comms5, only: start_halo_update_5
#else
    use comms5, only: update_halo_5
#endif
    complex(dp), intent(out) :: xout(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: xin(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    real, intent(in) :: am
    integer, intent(in) :: imass

#ifdef MPI
    integer, dimension(16) :: reqs_xtemp, reqs_xout
#endif
    complex(dp) :: xtemp(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

    call dslash(xtemp, xin, u, am, imass)
#ifdef MPI
    call start_halo_update_5(4, xtemp, 8, reqs_xtemp)
    call complete_halo_update(reqs_xtemp)
#else
    call update_halo_5(4, xtemp)
#endif
    call dslashd(xout, xtemp, u, am, imass)
#ifdef MPI
    call start_halo_update_5(4, xout, 10, reqs_xout)
    call complete_halo_update(reqs_xout)
#else
    call update_halo_5(4, xout)
#endif

  end subroutine dirac_operator

end module
