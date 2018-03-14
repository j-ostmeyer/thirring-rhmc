module random
! Random numbers
  real :: yran
  integer :: idum
  real :: v(97)

  DOUBLE PRECISION, PRIVATE :: DS(2),    DM(2)
  DOUBLE PRECISION, PRIVATE :: DX24,     DX48
  DATA      DS     /  16651885.D0, 2868876.D0  /
  DATA      DM     /  15184245.D0, 2651554.D0  /
  DATA      DX24   /  16777216.D0  /
  DATA      DX48   /  281474976710656.D0  /

contains

!*****************************************
!  Random number generator Numerical recipes 7.1
!
  real function rano(y,i)
    real, intent(out) :: y
    integer, intent(inout) :: i
    integer :: j
    real :: dum
    !     
    if(i.lt.0)then
       i=1
       do j=1,97
          dum=rranf()
       enddo
       do j=1,97
          v(j)=rranf()
       enddo

       y=rranf()
    endif
    !     
    j=1+int(97.0*y)
    if(j.gt.97) j=97
    if(j.lt.1) j=1
    !     write(6,*) j,y
    !     write(6,*) 'problems with rano'
    !     stop
    !     endif
    y=v(j)
    rano=y
    v(j)=rranf()
    return
  end function rano
!========================================================================
!     
  SUBROUTINE RRANGET(LSEED)
    DOUBLE PRECISION, INTENT(OUT) :: LSEED
    LSEED  =  G900GT()
    RETURN
  END SUBROUTINE RRANGET

  SUBROUTINE RRANSET(LSEED)
    DOUBLE PRECISION, INTENT(IN) :: LSEED
    DOUBLE PRECISION DUMMY
    DUMMY  =  G900ST(LSEED)
    RETURN
  END SUBROUTINE RRANSET



  REAL FUNCTION RRANF()
    RRANF = SNGL(DRANF())
    RETURN
  END FUNCTION RRANF

  DOUBLE PRECISION FUNCTION DRANF()
    DOUBLE PRECISION    DL,       DC,       DU,       DR
    DL  =  DS(1) * DM(1)
    DC  =  DINT(DL/DX24)
    DL  =  DL - DC*DX24
    DU  =  DS(1)*DM(2) + DS(2)*DM(1) + DC
    DS(2)  =  DU - DINT(DU/DX24)*DX24
    DS(1)  =  DL
    DR     =  (DS(2)*DX24 + DS(1)) / DX48
    DRANF  =  DR
    RETURN
  END FUNCTION DRANF

  DOUBLE PRECISION FUNCTION G900GT()
    G900GT  =  DS(2)*DX24 + DS(1)
    RETURN
  END FUNCTION G900GT

  DOUBLE PRECISION FUNCTION G900ST(DSEED)
    DOUBLE PRECISION, INTENT(IN) :: DSEED
    DS(2)  =  DINT(DSEED/DX24)
    DS(1)  =  DSEED - DS(2)*DX24
    G900ST =  DS(1)
    RETURN
  END FUNCTION G900ST
!***********************************************************************

end module random
