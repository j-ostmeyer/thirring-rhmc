module remez
  use params
  implicit none
  save

  real(dp) :: anum2(0:ndiag), aden2(ndiag)
  real(dp) :: bnum2(0:ndiag), bden2(ndiag)
  real(dp) :: anum4(0:ndiag), aden4(ndiag)
  real(dp) :: bnum4(0:ndiag), bden4(ndiag)

end module remez

