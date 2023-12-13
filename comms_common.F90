module comms_common
  implicit none
  save
  integer :: ip_x, ip_y, ip_t, ip_global, np_global
  integer :: ip_third = 0 ! helps align measureW with master
  integer :: comm, comm_grp_third
  integer :: mpiio_type
end module comms_common
