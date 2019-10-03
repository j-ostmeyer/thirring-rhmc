module trial
  use params
  implicit none
  save

  complex(dp) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
  real :: theta(ksizex_l, ksizey_l, ksizet_l, 3)
  real :: pp(ksizex_l, ksizey_l, ksizet_l, 3)
end module trial
