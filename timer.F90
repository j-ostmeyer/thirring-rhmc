! The philosophy of this module is the following:
! 1. in case we're using MPI, we can rely on MPI_Wtime and all is
!    easy and straightforward.
! 2. Otherwise, we must rely on system_clock, which gives a value
!    which wraps around and is always smaller than count_max, so we
!    need to compensate for that.
! Note: since this program is mostly run with MPI, case 2 is likely to
!       be less tested.

module timer
  implicit none

  real :: initial_time
  logical :: initialised = .false.
  private :: initial_time, initialised

#ifndef MPI
  real :: last_time
  private :: last_time
#endif

contains

  subroutine initialise()
#ifdef MPI
    use mpi
    use comms_common
    implicit none
    integer :: ierr
#else
    implicit none
    integer :: count, count_rate, count_max
#endif

#ifdef MPI
    initial_time = real(MPI_Wtime())
    call MPI_AllReduce(MPI_In_Place, initial_time, 1, MPI_Real, MPI_Min, comm, ierr)
#else
    call system_clock(count, count_rate, count_max)
    initial_time = real(count)/count_rate
    last_time = initial_time
#endif

    initialised = .true.

  end subroutine initialise

  function get_time_from_start() result(time_from_start)
#ifdef MPI
    use mpi
    use comms_common
    implicit none
    integer :: ierr
#else
    implicit none
    integer :: count, count_rate, count_max
#endif
    real :: time, time_from_start

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

#ifdef MPI
    time = real(MPI_Wtime())
    call MPI_AllReduce(MPI_In_Place, time, 1, MPI_Real, MPI_Max, comm, ierr)
#else
    call system_clock(count, count_rate, count_max)
    time = real(count)/count_rate
    do while (time .lt. last_time)
      time = time + real(count_max)/count_rate
    enddo
    last_time = time
#endif

    time_from_start = time - initial_time

  end function get_time_from_start

end module timer
