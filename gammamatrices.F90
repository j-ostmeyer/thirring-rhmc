module gammamatrices
  use params
  implicit none
  complex(dp) :: gamval(6,4)
  integer :: gamin(6,4)
contains 

  pure integer function kdelta(nu, mu)
    integer, intent(in) :: nu
    integer, intent(in) :: mu

    kdelta=merge(1,0,nu==mu)
  end function kdelta

end module gammamatrices
