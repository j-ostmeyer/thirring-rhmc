module remezg
  use params
  implicit none
  save

  real(dp) :: anum2g(0:ndiagg), aden2g(ndiagg)
  real(dp) :: bnum2g(0:ndiagg), bden2g(ndiagg)
  real(dp) :: anum4g(0:ndiagg), aden4g(ndiagg)
  real(dp) :: bnum4g(0:ndiagg), bden4g(ndiagg)
end module remezg
