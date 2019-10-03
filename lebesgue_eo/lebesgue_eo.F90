module lebesgue_eo
  implicit none
  ! lebesgue factors
  integer, parameter :: flength = 4
  integer :: xf(flength), yf(flength), tf(flength)

contains

  subroutine find_factors(l, factors)
    implicit none
    integer, intent(in) :: l
    integer, intent(out) :: factors(flength)

    integer :: tmpl, tmpfactor, ifactor

    factors = 1
    tmpl = l
    ifactor = 1
    do while (tmpl .gt. 1) .and. (ifactor .le. flength)
      ! finding a factor
      tmpfactor = 2
      do while (mod(tmpl, tmpfactor) .ne. 0)
        tmpfactor = 1 + tmpfactor
      enddo
      tmpl = tmpl/tmpfactor
      factors(ifactor) = tmpfactor
      ifactor = 1 + ifactor
    enddo

  end subroutine find_factors

  subroutine xyt_lebeo(x, y, t, xf, yf, tfh, lx, ly, lth, lebeo)
    integer, intent(in) :: x
    integer, intent(in) :: y
    integer, intent(in) :: t
    integer, intent(in) :: xf(flength)
    integer, intent(in) :: yf(flength)
    integer, intent(in) :: tfh(flength)
    integer, intent(in) :: lx
    integer, intent(in) :: ly
    integer, intent(in) :: lth

    integer, intent(out) :: lebeo

    integer :: par
    integer :: ii
    integer :: tx, ty, tth

    par = mod(x + y + t, 2)

    lebeo = 1 + par*lx*ly*lth

    do ii = 1, flength

    enddo

    end module lebesgue_eo
