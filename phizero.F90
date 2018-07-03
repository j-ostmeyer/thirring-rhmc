module phizero
  use params
  implicit none
  save

  complex(dp) :: Phi0(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4,ndiag)
end module phizero

