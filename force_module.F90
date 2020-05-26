module force_module
  implicit none
  real, parameter :: One = 1.0
contains

!******************************************************************
!   calculate dSds for gauge fields at each intermediate time
!******************************************************************
  subroutine force(Phi, res1, am, imass)
    use comms
    use counters, only: ancgpv, ancg, ancgf, ancgfpv
    use gforce, only: dSdpi
    use params, only: kthird_l, ksizex_l, ksizey_l, ksizet_l, dp
    use qmrherm_module, only: qmrherm
    use remezg
    use trial
    implicit none

    complex(dp), intent(in) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4, Nf)
    real, intent(in) :: res1, am
    integer, intent(in) :: imass
!     complex Phi(kferm,Nf),X2(kferm)
!     complex X1,u
    complex(dp) :: X2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: Xresult(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer :: ia, itercg

    dSdpi = 0.0

! uncomment this line to quench the fermions!
!     return

! pseudofermion action is
!   Phi^dagger {MdaggerM(1)}^1/4 {MdaggerM(m)})^-1/2 {MdaggerM(1)}^1/4 Phi
!
    do ia = 1, Nf
!
      X2 = Phi(:, :, :, :, :, ia)

      call qmrherm(X2, Xresult, res1, itercg, One, 1, anum4g, aden4g, ndiagg, 1, spmd)
      ancgpv = ancgpv + float(itercg)

      X2 = Xresult
!
      call qmrherm(X2, Xresult, res1, itercg, am, imass, bnum2g, bden2g, ndiagg, 0, spmd)
      ancg = ancg + float(itercg)
!     write(111,*) itercg
      X2 = Xresult
!
!  evaluates -X2dagger * d/dpi[{MdaggerM(m)}^1/2] * X2
      call qmrherm(X2, Xresult, res1, itercg, am, imass, anum2g, aden2g, ndiagg, 2, spmd)
      ancgf = ancgf + float(itercg)

!     write(113,*) itercg
!  evaluates +2Re{Phidagger * d/dpi[{MdaggerM(1)}^1/4] * X2}
      call qmrherm(X2, Xresult, res1, itercg, One, 1, anum4g, aden4g, ndiagg, 3, spmd)
      ancgfpv = ancgfpv + float(itercg)
!
    enddo
!
    if (ibound .eq. -1 .and. ip_t .eq. (np_t - 1)) then
      dSdpi(:, :, ksizet_l, 3) = -dSdpi(:, :, ksizet_l, 3)
    endif

    dSdpi = dSdpi + beta*Nf*theta
!
    return
  end subroutine force

end module force_module
