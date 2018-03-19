module mpi_variables
  implicit none
  save

  integer :: ip_x, ip_y, ip_t, ip_global, np_global
  integer :: comm, ierr
  integer :: mpiio_type
end module mpi_variables
