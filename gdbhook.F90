module gdbhook

contains

  subroutine gdb_wait()
#ifdef MPI
    use mpi
    implicit none
    logical :: flag
    integer :: rank, ierr
    ! for the use of get_environment_variable
    character :: env_value(100)
    integer :: env_status

    flag = .true.
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    call get_environment_variable(name='GDBHOOK', &
                                  value=env_value, &
                                  status=env_status)
    if (rank .eq. 0 .and. env_status .eq. 0) then
      print *, 'Waiting on master for intervention with GDB...'
      print *, 'Attach to master with gdb -p', getpid()
      do while (flag)
        call sleep(1)
      enddo
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
  end subroutine gdb_wait
end module gdbhook
