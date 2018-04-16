program test_save
      use dwf3d_lib
      use gauge
      use comms
      implicit none

! setup
#ifdef MPI
      call init_MPI
#endif
      call sread

! call function
      call swrite

#ifdef MPI
      call MPI_Finalize
#endif
end program
