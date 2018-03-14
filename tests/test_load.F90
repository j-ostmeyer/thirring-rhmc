program test_load
      use dwf3d_lib
      use gauge
#ifdef MPI
      use mpi
      use mpi_variables
#endif
      implicit none

      real :: sumtheta, maxtheta, mintheta

! initialise MPI
#ifdef MPI
      call init_MPI
#endif

! call function
      call sread

! check output
      sumtheta = sum(theta)
      maxtheta = maxval(theta)
      mintheta = minval(theta)
#ifdef MPI
      call MPI_AllReduce(MPI_IN_PLACE, sumtheta, 1, MPI_REAL, MPI_SUM, comm, ierr)
      call MPI_AllReduce(MPI_IN_PLACE, maxtheta, 1, MPI_REAL, MPI_MAX, comm, ierr)
      call MPI_AllReduce(MPI_IN_PLACE, mintheta, 1, MPI_REAL, MPI_MIN, comm, ierr)
#endif
      if (ip_global .eq. 0) then
         print *, 'sum = ', sumtheta
         print *, 'max = ', maxtheta
         print *, 'min = ', mintheta
      end if
end program
