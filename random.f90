! This module abstracts out the random number functionality found in the original
! version of the code. Interface is changed to be compatible with the site-based
! version in site_random.f90.

module random
  use comms
  implicit none
  save

  ! Shared externally-visible state, previously held as common blocks
  real :: yran
  integer :: idum
  real :: v(97)

  ! Internal parameters and state of the generator
  DOUBLE PRECISION, PRIVATE :: DS(2)
  DOUBLE PRECISION, PRIVATE, PARAMETER :: DM(2) = (/15184245.D0, 2651554.D0/)
  DOUBLE PRECISION, PRIVATE, PARAMETER :: DX24 = 16777216.D0
  DOUBLE PRECISION, PRIVATE, PARAMETER :: DX48 = 281474976710656.D0
  DATA DS/16651885.D0, 2868876.D0/

contains

  !*****************************************
  !  Random number generator Numerical recipes 7.1
  !
  ! Generate a random number with good entropy from a generator with bad entropy.
  real function rano(y, i, ix, iy, it)
    real, intent(out) :: y
    integer, optional, intent(in) :: ix, iy, it
    integer, intent(inout) :: i
    integer :: j
    real :: dum
    !
    if (i .lt. 0) then
      i = 1
      do j = 1, 97
        dum = rranf()
      enddo
      do j = 1, 97
        v(j) = rranf()
      enddo

      y = rranf()
    endif
    !
    j = 1 + int(97.0*y)
    if (j .gt. 97) j = 97
    if (j .lt. 1) j = 1
    !     write(6,*) j,y
    !     write(6,*) 'problems with rano'
    !     stop
    !     endif
    y = v(j)
    rano = y
    v(j) = rranf()
    return
  end function rano
  !========================================================================
  !
  ! Get the state of the generator
  SUBROUTINE RRANGET(LSEED, ix, iy, it)
    integer, optional, intent(in) :: ix, iy, it
    DOUBLE PRECISION, INTENT(OUT) :: LSEED
    LSEED = G900GT()
    RETURN
  END SUBROUTINE RRANGET

  ! Set the state of the generator
  SUBROUTINE RRANSET(LSEED, ix, iy, it)
    integer, optional, intent(in) :: ix, iy, it
    DOUBLE PRECISION, INTENT(IN) :: LSEED
    DOUBLE PRECISION DUMMY
    DUMMY = G900ST(LSEED)
    RETURN
  END SUBROUTINE RRANSET

  ! Get a random single-precision real number
  REAL FUNCTION RRANF(ix, iy, it)
    integer, optional, intent(in) :: ix, iy, it
    RRANF = SNGL(DRANF())
    RETURN
  END FUNCTION RRANF

  ! Get a random double-precision real number
  ! Has less entropy than rano
  DOUBLE PRECISION FUNCTION DRANF(ix, iy, it)
    integer, optional, intent(in) :: ix, iy, it
    DOUBLE PRECISION DL, DC, DU, DR
    DL = DS(1)*DM(1)
    DC = DINT(DL/DX24)
    DL = DL - DC*DX24
    DU = DS(1)*DM(2) + DS(2)*DM(1) + DC
    DS(2) = DU - DINT(DU/DX24)*DX24
    DS(1) = DL
    DR = (DS(2)*DX24 + DS(1))/DX48
    DRANF = DR
    RETURN
  END FUNCTION DRANF

  ! Actually get the state of the generator
  DOUBLE PRECISION FUNCTION G900GT(ix, iy, it)
    integer, optional, intent(in) :: ix, iy, it
    G900GT = DS(2)*DX24 + DS(1)
    RETURN
  END FUNCTION G900GT

  ! Actually set the state of the generator
  DOUBLE PRECISION FUNCTION G900ST(DSEED, ix, iy, it)
    integer, optional, intent(in) :: ix, iy, it
    DOUBLE PRECISION, INTENT(IN) :: DSEED
    DS(2) = DINT(DSEED/DX24)
    DS(1) = DSEED - DS(2)*DX24
    G900ST = DS(1)
    RETURN
  END FUNCTION G900ST
  !***********************************************************************
  ! Seed the generator and call rano to initialise its state
  subroutine init_random(seed)
    real(dp), intent(in) :: seed
    real :: y

    call rranset(seed + ip_global)
    idum = -1
    y = rano(yran, idum)
  end subroutine init_random

end module random
