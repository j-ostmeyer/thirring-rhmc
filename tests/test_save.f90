program test_dslash
      implicit none
! function to test
      external :: swrite

! supporting functions
      external :: sread

! general parameters
      integer, parameter :: ksize=12, ksizet=12

! common blocks to function
      common/gauge/ theta(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3), seed
      real :: theta
      real*8 :: seed

! setup
      call sread

! call function
      call swrite
end program
