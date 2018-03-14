program test_dslash
      use dwf3d_lib
      implicit none

! common blocks to function
      common/gauge/ theta(0:ksize+1, 0:ksize+1, 0:ksizet+1, 3), seed
      real :: theta
      real*8 :: seed

! call function
      call sread

! check output
      print *, 'sum = ', sum(theta)
      print *, 'max = ', maxval(theta)
      print *, 'min = ', minval(theta)
end program
