module test_utils
  use params
  use random
  use dirac
#ifdef MPI
  use mpi
  use comms
  use comms4
  use comms5
#endif

  implicit none

  ! general parameters
  logical :: generate = .false.
  integer, dimension(3) :: imasses = (/1,3,5/)

contains

subroutine generate_starting_state_Phi(Phi, u, reqs_Phi)
    implicit none
    ! Required inputs
    complex(dp), intent(inout) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(inout) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    integer, dimension(16), intent(inout) :: reqs_Phi

    integer, dimension(16) :: reqs_X, reqs_R
    integer, dimension(12) :: reqs_u
    complex(dp) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: X(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

    call generate_starting_state_Phi_and_R_and_X(Phi, R, u, X, reqs_Phi, reqs_R, reqs_X)
end subroutine generate_starting_state_Phi

subroutine generate_starting_state_Phi_and_R(Phi, R, u, reqs_Phi, reqs_R)
    implicit none
    ! Required inputs
    complex(dp), intent(inout) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(inout) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    complex(dp), intent(inout) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer, dimension(16), intent(inout) :: reqs_Phi, reqs_R

    integer, dimension(16) :: reqs_X
    integer, dimension(12) :: reqs_u
    complex(dp) :: X(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

    call generate_starting_state_Phi_and_R_and_X(Phi, R, u, X, reqs_Phi, reqs_R, reqs_X)
end subroutine generate_starting_state_Phi_and_R

subroutine generate_starting_state_Phi_and_X(Phi, X, u, reqs_Phi, reqs_X)
    implicit none
    ! Required inputs
    complex(dp), intent(inout) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(inout) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    complex(dp), intent(inout) :: X(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer, dimension(16), intent(inout) :: reqs_Phi, reqs_X

    integer, dimension(16) :: reqs_R
    integer, dimension(12) :: reqs_u
    complex(dp) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

    call generate_starting_state_Phi_and_R_and_X(Phi, R, u, X, reqs_Phi, reqs_R, reqs_X)
end subroutine generate_starting_state_Phi_and_X

subroutine generate_starting_state_Phi_and_R_and_X(Phi, R, u, X, reqs_Phi, reqs_R, reqs_X)
    implicit none
    ! Required inputs
    complex(dp), intent(inout) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(inout) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(inout) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    complex(dp), intent(inout) :: X(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer, dimension(16), intent(inout) :: reqs_Phi, reqs_R, reqs_X

    real :: dSdpi_ref(ksizex_l, ksizey_l, ksizet_l, 3)
    complex, parameter :: iunit = cmplx(0, 1)
    real(dp), parameter :: tau = 8*atan(1.0_8)

    integer :: i, j, ix, iy, it, ithird
    integer, parameter :: idxmax = 4*ksize*ksize*ksizet*kthird
    integer :: idx
#ifdef MPI
    integer, dimension(12) :: reqs_u
    integer :: ierr
#endif

    call generate_starting_state_Phi_and_R_and_X_and_dsdpi_ref(Phi, R, u, X, reqs_Phi, reqs_R, reqs_X, dSdpi_ref)

end subroutine generate_starting_state_Phi_and_R_and_X

subroutine generate_starting_state_Phi_and_R_and_X_and_dsdpi_ref(Phi, R, u, X, reqs_Phi, reqs_R, reqs_X, dSdpi_ref)
    implicit none
    ! Required inputs
    complex(dp), intent(inout) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(inout) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(inout) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    complex(dp), intent(inout) :: X(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real, intent(inout) :: dSdpi_ref(ksizex_l, ksizey_l, ksizet_l, 3)
    integer, dimension(16), intent(inout) :: reqs_Phi, reqs_R, reqs_X

    complex, parameter :: iunit = cmplx(0, 1)
    real(dp), parameter :: tau = 8*atan(1.0_8)

    integer :: i, j, ix, iy, it, ithird
    integer, parameter :: idxmax = 4*ksize*ksize*ksizet*kthird
    integer :: idx
#ifdef MPI
    integer, dimension(12) :: reqs_u
    integer :: ierr
#endif

    do j = 1, 4
      do it = 1, ksizet_l
        do iy = 1, ksizey_l
          do ix = 1, ksizex_l
            do ithird = 1, kthird_l
              idx = ip_third*kthird_l + ithird &
                    + (ip_x*ksizex_l + ix - 1)*kthird &
                    + (ip_y*ksizey_l + iy - 1)*kthird*ksize &
                    + (ip_t*ksizet_l + it - 1)*kthird*ksize*ksize &
                    + (j - 1)*kthird*ksize*ksize*ksizet

              Phi(ithird, ix, iy, it, j) = 1.1*exp(iunit*idx*tau/idxmax)
              R(ithird, ix, iy, it, j) = 1.3*exp(iunit*idx*tau/idxmax)
              X(ithird, ix, iy, it, j) = 0.5*exp(1.0)*exp(iunit*idx*tau/idxmax)
            enddo
          enddo
        enddo
      enddo
    enddo
#ifdef MPI
    call start_halo_update_5(4, R, 0, reqs_R)
    call start_halo_update_5(4, Phi, 1, reqs_Phi)
    call start_halo_update_5(4, X, 0, reqs_X)
#endif
      do j = 1, 3
        do it = 1, ksizet_l
          do iy = 1, ksizey_l
            do ix = 1, ksizex_l
              idx = ip_x*ksizex_l + ix &
                    + (ip_y*ksizey_l + iy - 1)*ksize &
                    + (ip_t*ksizet_l + it - 1)*ksize*ksize &
                    + (j - 1)*ksize*ksize*ksizet

              u(ix, iy, it, j) = exp(iunit*idx*tau/idxmax)
              dSdpi_ref(ix, iy, it, j) = real(tau*exp(iunit*idx*tau/idxmax), sp)
            enddo
          enddo
        enddo
      enddo
#ifdef MPI
    call start_halo_update_4(3, u, 1, reqs_u)
    call complete_halo_update(reqs_R)
    call complete_halo_update(reqs_Phi)
    call complete_halo_update(reqs_X)
    call MPI_WaitAll(12, reqs_u, MPI_STATUSES_IGNORE, ierr)
#else
    call update_halo_5(4, R)
    call update_halo_5(4, Phi)
    call update_halo_5(4, X)
    call update_halo_4(3, u)
#endif

    ! initialise common variables
    beta = 0.4
    am3 = 1.0
    ibound = -1

    call init_gammas()
  end subroutine generate_starting_state_Phi_and_R_and_X_and_dsdpi_ref

#ifdef MPI

  subroutine rw_file_mpi(array, array_shape, rank, filename, mpi_dtype, write_out)
    integer, intent(in) :: array_shape(:), rank
    type(*), dimension(..), intent(inout) :: array
    integer, intent(in) :: mpi_dtype
    character(len=*), intent(in) :: filename
    logical, intent(in) :: write_out
    integer, dimension(rank) :: global_size, local_size, start
    integer :: size_index
    integer :: count
    integer :: offset
    integer :: mpi_fh
    integer :: local_mpiio_type
    integer :: status(MPI_STATUS_SIZE)
    integer :: mode
    integer :: ierr


    if (write_out) then
      if (ip_global == 0) then
        print *, "  Creating dat file ", filename
      end if
      mode = MPI_Mode_Wronly + MPI_Mode_Create
    else
      mode = MPI_Mode_Rdonly
    end if

    global_size = array_shape
    local_size = array_shape
    start = 0

    ! Work out which elements of local_size will differ from global_size
    if (rank .eq. 4) then
      offset = 0
    else
      offset = 1
      ! Less changes if we modify the first index here
      global_size(0 + offset) = kthird
      start(0 + offset) = ip_third*kthird_l
    end if

    ! Work out how many elements to transfer
    count = 1
    do size_index = 1, rank
      count = count*local_size(size_index)
    end do

    global_size(1 + offset) = ksize
    global_size(2 + offset) = ksize
    global_size(3 + offset) = ksizet

    start(1 + offset) = ip_x*ksizex_l
    start(2 + offset) = ip_y*ksizey_l
    start(3 + offset) = ip_t*ksizet_l

    call MPI_Type_Create_Subarray(rank, global_size, local_size, start, &
                                  MPI_Order_Fortran, mpi_dtype, local_mpiio_type, &
                                  ierr)
    call MPI_Type_Commit(local_mpiio_type, ierr)
    call MPI_File_Open(comm, filename, mode, MPI_Info_Null, mpi_fh, ierr)
    call MPI_File_Set_View(mpi_fh, 4_8, mpi_dtype, local_mpiio_type, &
                           "native", MPI_Info_Null, ierr)
    if (write_out) then
      call MPI_File_Write_All(mpi_fh, array, count, mpi_dtype, status, ierr)
    else
      call MPI_File_Read_All(mpi_fh, array, count, mpi_dtype, status, ierr)
    end if

    call MPI_File_Close(mpi_fh, ierr)
    call MPI_Type_Free(local_mpiio_type, ierr)

    return
  end subroutine rw_file_mpi
#endif
end module
