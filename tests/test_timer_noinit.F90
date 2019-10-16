program test_timer_noinit
  use timer, only: get_time_from_start
#ifdef MPI
  use mpi
  implicit none
  integer :: ierr
#else
  implicit none
#endif
  real :: time

#ifdef MPI
  call MPI_Init(ierr)
#endif
  ! Program must crash here
  time = get_time_from_start()
#ifdef MPI
  ! program should not get here
  call MPI_Finalize(ierr)
#endif

end program test_timer_noinit
