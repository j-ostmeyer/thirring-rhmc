module timer
  implicit none
  real :: last_time, initial_time
  logical :: initialised = .false.

  private :: last_time, initial_time, initialised

contains

  subroutine initialise()
#ifdef MPI
    use mpi
    use comms_common
#endif
    implicit none
    integer :: count, count_rate, count_max

    integer :: ierr

    call system_clock(count, count_rate, count_max)
    initial_time = real(count)/count_rate

#ifdef MPI
    call MPI_AllReduce(MPI_In_Place, initial_time, 1, MPI_Real, MPI_Min, comm, ierr)
#endif
    last_time = initial_time

    initialised = .true.

  end subroutine initialise

  function get_time() result(time_from_start)
#ifdef MPI
    use mpi
    use comms_common
#endif
    implicit none
    real :: time, time_from_start

    integer :: count, count_rate, count_max
    integer :: ierr

    if (.not. initialised) then
#ifdef MPI
      if (ip_global .eq. 0) then
#endif
        print *, "Timer not initialised!"
#ifdef MPI
      endif
#endif
#ifdef MPI
      call MPI_Abort(comm, 1, ierr)
#endif
      stop
    endif

    call system_clock(count, count_rate, count_max)
    time = real(count)/count_rate
#ifdef MPI
    call MPI_AllReduce(MPI_In_Place, time, 1, MPI_Real, MPI_Max, comm, ierr)
#endif
    do while (time .lt. last_time)
      time = time + real(count_max)/count_rate
    enddo
    last_time = time

    time_from_start = time - initial_time

  end function get_time

end module timer
