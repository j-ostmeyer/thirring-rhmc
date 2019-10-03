module dum1
  use params
  implicit none
  save

  complex(dp) :: R(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
  real :: ps(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 2)
end module dum1
