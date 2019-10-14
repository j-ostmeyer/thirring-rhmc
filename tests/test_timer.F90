program test_timer_noinit
  use timer, only: get_time_from_start, initialise
#ifdef MPI
  use mpi
  implicit none
  integer :: ierr
#else
  implicit none
#endif
  real :: time
  integer, parameter :: sleep_time = 5

#ifdef MPI
  call MPI_Init(ierr)
#endif
  call initialise
  call sleep(sleep_time)
  time = get_time_from_start()

  if (abs(time - sleep_time) .ge. 1.0) then
    print *, "error: get_time_from_start did not behave correctly"
    print *, time, ' vs ', sleep_time
  endif

  call MPI_Finalize(ierr)

end program test_timer_noinit
