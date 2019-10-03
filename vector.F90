module vector
  use params
  implicit none
  save

  complex(dp) :: X(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
end module vector
