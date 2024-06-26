module gammamatrices
  use params
  implicit none
  complex(dp) :: gamval(6, 4)
  complex(sp) :: gamvalf(6, 4)
  integer :: gamin(6, 4)
contains

  pure integer function kdelta(nu, mu)
    implicit none
    integer, intent(in) :: nu
    integer, intent(in) :: mu

    kdelta = merge(1, 0, nu == mu)
  end function kdelta

  subroutine init_gammas()
    implicit none
    complex(dp), parameter :: one = (1.0, 0.0), zi = (0.0, 1.0)
    complex(sp), parameter :: onef = (1.0, 0.0), zif = (0.0, 1.0)
!*******************************************************************
!    setup Dirac algebra
!*******************************************************************
!
!     gamma_1
!
    gamval(1, 1) = -zi
    gamval(1, 2) = -zi
    gamval(1, 3) = zi
    gamval(1, 4) = zi
!
    gamvalf(1, 1) = -zif
    gamvalf(1, 2) = -zif
    gamvalf(1, 3) = zif
    gamvalf(1, 4) = zif
!
    gamin(1, 1) = 4
    gamin(1, 2) = 3
    gamin(1, 3) = 2
    gamin(1, 4) = 1
!
!     gamma_2
!
    gamval(2, 1) = -one
    gamval(2, 2) = one
    gamval(2, 3) = one
    gamval(2, 4) = -one
!
    gamvalf(2, 1) = -onef
    gamvalf(2, 2) = onef
    gamvalf(2, 3) = onef
    gamvalf(2, 4) = -onef

!
    gamin(2, 1) = 4
    gamin(2, 2) = 3
    gamin(2, 3) = 2
    gamin(2, 4) = 1
!
!     gamma_3
!
    gamval(3, 1) = -zi
    gamval(3, 2) = zi
    gamval(3, 3) = zi
    gamval(3, 4) = -zi
!
    gamvalf(3, 1) = -zif
    gamvalf(3, 2) = zif
    gamvalf(3, 3) = zif
    gamvalf(3, 4) = -zif
!
    gamin(3, 1) = 3
    gamin(3, 2) = 4
    gamin(3, 3) = 1
    gamin(3, 4) = 2
!
!     gamma_4
!
    gamval(4, 1) = one
    gamval(4, 2) = one
    gamval(4, 3) = -one
    gamval(4, 4) = -one
!
    gamvalf(4, 1) = onef
    gamvalf(4, 2) = onef
    gamvalf(4, 3) = -onef
    gamvalf(4, 4) = -onef
!
    gamin(4, 1) = 1
    gamin(4, 2) = 2
    gamin(4, 3) = 3
    gamin(4, 4) = 4
!
!     gamma_5 = gamma_1 * gamma_2 * gamma_3 * gamma_4
!
    gamval(5, 1) = -one
    gamval(5, 2) = -one
    gamval(5, 3) = -one
    gamval(5, 4) = -one
!
    gamvalf(5, 1) = -onef
    gamvalf(5, 2) = -onef
    gamvalf(5, 3) = -onef
    gamvalf(5, 4) = -onef
!
    gamin(5, 1) = 3
    gamin(5, 2) = 4
    gamin(5, 3) = 1
    gamin(5, 4) = 2
!
!     gamma_4 * gamma_5 (called gamma_3 gamma_5 in notes)
    gamval(6, 1) = -one
    gamval(6, 2) = -one
    gamval(6, 3) = one
    gamval(6, 4) = one
!
    gamvalf(6, 1) = -onef
    gamvalf(6, 2) = -onef
    gamvalf(6, 3) = onef
    gamvalf(6, 4) = onef
!
    gamin(6, 1) = 3
    gamin(6, 2) = 4
    gamin(6, 3) = 1
    gamin(6, 4) = 2
!
!
    gamval = gamval*akappa
    gamvalf = gamvalf*akappaf
!

  end subroutine init_gammas

end module gammamatrices
