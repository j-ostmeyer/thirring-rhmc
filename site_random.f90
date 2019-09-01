! This module matches the original abstracted-out random module in terms of
! its interface, but allows a specific random number chain to be generated
! for each site, allowing identical results to be obtained on disparate
! parallelisations.

module random
  use params
  use comms
  implicit none
  save

  ! Shared externally-visible state, previously held as common blocks
  real :: yran(ksizex_l, ksizey_l, ksizet_l)
  integer :: idum(ksizex_l, ksizey_l, ksizet_l)
  real :: v(97, ksizex_l, ksizey_l, ksizet_l)

  ! Internal parameters and state of the generator
  DOUBLE PRECISION, PRIVATE :: DS(2, ksizex_l, ksizey_l, ksizet_l) = 0.0
  DOUBLE PRECISION, PRIVATE, PARAMETER :: DM(2) = (/15184245.D0, 2651554.D0/)
  DOUBLE PRECISION, PRIVATE, PARAMETER :: DX24 = 16777216.D0
  DOUBLE PRECISION, PRIVATE, PARAMETER :: DX48 = 281474976710656.D0
  !DATA      DS     /  16651885.D0, 2868876.D0  /
  ! RRANSET must be called for each (ix, iy, it) before first call to rano to initialise DS
  ! Provided this is done, this set of DS is not needed

contains

  !*****************************************
  !  Random number generator Numerical recipes 7.1
  !
  ! Generate a random number with good entropy from a generator with bad entropy.
  real function rano(y, i, ix, iy, it)
    real, intent(out) :: y(ksizex_l, ksizey_l, ksizet_l)
    integer, intent(inout) :: i(ksizex_l, ksizey_l, ksizet_l)
    integer, intent(in) :: ix, iy, it
    integer :: j
    real :: dum
    !
    if (i(ix, iy, it) .lt. 0) then
      i(ix, iy, it) = 1
      do j = 1, 97
        dum = rranf(ix, iy, it)
      enddo
      do j = 1, 97
        v(j, ix, iy, it) = rranf(ix, iy, it)
      enddo

      y(ix, iy, it) = rranf(ix, iy, it)
    endif
    !
    j = 1 + int(97.0*y(ix, iy, it))
    if (j .gt. 97) j = 97
    if (j .lt. 1) j = 1
    !     write(6,*) j,y
    !     write(6,*) 'problems with rano'
    !     stop
    !     endif
    y(ix, iy, it) = v(j, ix, iy, it)
    rano = y(ix, iy, it)
    v(j, ix, iy, it) = rranf(ix, iy, it)
    return
  end function rano
  !========================================================================
  !
  ! Calculate an offset for the random seed based on the site location
  function seed_offset(ix, iy, it)
    integer :: seed_offset
    integer, intent(in) :: ix, iy, it
    seed_offset = ip_x*ksizex_l + ix - 1 &
      & + ksize*(ip_y*ksizey_l + iy - 1) &
      & + ksize*ksize*(ip_t*ksizet_l + it - 1)
    return
  end function seed_offset

  ! Get the state of the generator
  SUBROUTINE RRANGET(LSEED, ix, iy, it)
    DOUBLE PRECISION, INTENT(OUT) :: LSEED
    integer, intent(in) :: ix, iy, it
    LSEED = G900GT(ix, iy, it)
    RETURN
  END SUBROUTINE RRANGET

  ! Set the state of the generator
  SUBROUTINE RRANSET(LSEED, ix, iy, it)
    DOUBLE PRECISION, INTENT(IN) :: LSEED
    DOUBLE PRECISION DUMMY
    integer, intent(in) :: ix, iy, it
    DUMMY = G900ST(seed_offset(ix, iy, it) + LSEED, ix, iy, it)
    RETURN
  END SUBROUTINE RRANSET

  ! Get a random single-precision real number
  REAL FUNCTION RRANF(ix, iy, it)
    integer, intent(in) :: ix, iy, it
    RRANF = SNGL(DRANF(ix, iy, it))
    RETURN
  END FUNCTION RRANF

  ! Get a random double-precision real number
  ! Has less entropy than rano
  DOUBLE PRECISION FUNCTION DRANF(ix, iy, it)
    integer, intent(in) :: ix, iy, it
    DOUBLE PRECISION DL, DC, DU, DR
    DL = DS(1, ix, iy, it)*DM(1)
    DC = DINT(DL/DX24)
    DL = DL - DC*DX24
    DU = DS(1, ix, iy, it)*DM(2) + DS(2, ix, iy, it)*DM(1) + DC
    DS(2, ix, iy, it) = DU - DINT(DU/DX24)*DX24
    DS(1, ix, iy, it) = DL
    DR = (DS(2, ix, iy, it)*DX24 + DS(1, ix, iy, it))/DX48
    DRANF = DR
    RETURN
  END FUNCTION DRANF

  ! Actually get the state of the generator
  DOUBLE PRECISION FUNCTION G900GT(ix, iy, it)
    integer, intent(in) :: ix, iy, it
    G900GT = DS(2, ix, iy, it)*DX24 + DS(1, ix, iy, it)
    RETURN
  END FUNCTION G900GT

  ! Actually set the state of the generator
  DOUBLE PRECISION FUNCTION G900ST(DSEED, ix, iy, it)
    DOUBLE PRECISION, INTENT(IN) :: DSEED
    integer, intent(in) :: ix, iy, it
    DS(2, ix, iy, it) = DINT(DSEED/DX24)
    DS(1, ix, iy, it) = DSEED - DS(2, ix, iy, it)*DX24
    G900ST = DS(1, ix, iy, it)
    RETURN
  END FUNCTION G900ST
  !***********************************************************************
  ! Seed the generator and call rano to initialise its state
  subroutine init_random(seed)
    real(dp), intent(in) :: seed
    integer :: ix, iy, it
    real :: y

    do it = 1, ksizet_l
      do iy = 1, ksizey_l
        do ix = 1, ksizex_l
          call rranset(seed, ix, iy, it)
          idum(ix, iy, it) = -1
          y = rano(yran, idum, ix, iy, it)
        end do
      end do
    end do
  end subroutine init_random

end module random
