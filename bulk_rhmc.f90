program dwf3d
  use dwf3d_lib
  implicit none
  call verify_kernel_choice()
  call dwf3d_main
end program

