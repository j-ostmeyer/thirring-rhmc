module gauge
  use params
  implicit none
  save

  real :: theta(ksizex_l, ksizey_l, ksizet_l, 3)

contains
!******************************************************************
!   calculate compact links from non-compact links
!******************************************************************
  subroutine coef(u, theta)
    use comms4, only: start_halo_update_4
    use comms_common, only: ip_t
    use comms, only: complete_halo_update
    use mpi
    implicit none
!
    complex(dp), intent(out) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    real, intent(in) :: theta(ksizex_l, ksizey_l, ksizet_l, 3)
#ifdef MPI
    integer, dimension(12) :: reqs_u
    integer :: ierr
#endif

!     u(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = exp(cmplx(0.0, theta))
    u(1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = (1.0 + cmplx(0.0, theta))
!
!  anti-p.b.c. in timelike direction
    if (ibound .eq. -1 .and. ip_t .eq. (np_t - 1)) then
      u(:, :, ksizet_l, 3) = -u(:, :, ksizet_l, 3)
    end if
!
!!!!    call complete_halo_update_4(3, u)
#ifdef MPI
    call start_halo_update_4(3, u, 3, reqs_u)
    ! call complete_halo_update(reqs_u)
    call MPI_WaitAll(12, reqs_u, MPI_Statuses_Ignore, ierr)
#else
    call update_halo_4(3, u)
#endif
    return
  end subroutine coef

end module gauge
