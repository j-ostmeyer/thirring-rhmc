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

subroutine generate_starting_state_Phi_and_X(Phi, reqs_Phi, u , X, reqs_X)
    implicit none
    ! Required inputs
    complex(dp), intent(inout) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(inout) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    complex(dp), intent(inout) :: X(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    integer, dimension(16), intent(inout) :: reqs_Phi, reqs_X

    integer, dimension(16) :: reqs_R
    complex(dp) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)

    call generate_starting_state_Phi_and_R_and_X(Phi, reqs_Phi, u, R, reqs_R, X, reqs_X)
end subroutine generate_starting_state_Phi_and_X

subroutine generate_starting_state(Phi, u, reqs_Phi, R, reqs_R, X, reqs_X, dSdpi_ref, Phi0_orig)
    implicit none
    ! Required inputs
    complex(dp), intent(inout) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(inout) :: u(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 3)
    integer, dimension(16), intent(inout) :: reqs_Phi
    ! Optional inputs
    complex(dp), intent(inout), optional :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(inout), optional :: X(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(inout), optional :: Phi0_orig(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4, 25)
    real, intent(inout), optional :: dSdpi_ref(ksizex_l, ksizey_l, ksizet_l, 3)
    integer, dimension(16), intent(inout), optional :: reqs_R, reqs_X

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
              if (present(R)) then
                R(ithird, ix, iy, it, j) = 1.3*exp(iunit*idx*tau/idxmax)
              end if
              if (present(X)) then
                X(ithird, ix, iy, it, j) = 0.5*exp(1.0)*exp(iunit*idx*tau/idxmax)
              end if
              if (present(Phi0_orig)) then
                do l = 1, 25
                  Phi0_orig(ithird, ix, iy, it, j, l) = 1.7*exp(1.0) &
                                                        *exp(iunit*idx*tau/idxmax) + l
                end do
              end if 
            enddo
          enddo
        enddo
      enddo
    enddo
#ifdef MPI
    call start_halo_update_5(4, Phi, 1, reqs_Phi)
    if (present(R)) then
      call start_halo_update_5(4, R, 0, reqs_R)
    end if
    if (present(X)) then
      call start_halo_update_5(4, X, 0, reqs_X)
    end if
    if (present(Phi0_orig)) then
      call start_halo_update_6(4, 25, Phi0_orig, 2, reqs_Phi0)
    end if 
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
              if (present(dSdpi_ref)) then
                dSdpi_ref(ix, iy, it, j) = real(tau*exp(iunit*idx*tau/idxmax), sp)
              end if
            enddo
          enddo
        enddo
      enddo
#ifdef MPI
    call start_halo_update_4(3, u, 1, reqs_u)
    call complete_halo_update(reqs_Phi)
    if (present(R)) then
      call complete_halo_update(reqs_R)
    end if
    if (present(X)) then
      call complete_halo_update(reqs_X)
    end if
    if (present(Phi0_orig)) then
      call complete_halo_update(reqs_Phi0)
    end if 
    call MPI_WaitAll(12, reqs_u, MPI_STATUSES_IGNORE, ierr)
#else
    call update_halo_4(3, u)
    call update_halo_5(4, Phi)
    if (present(R)) then
      call update_halo_5(4, R)
    end if
    if (present(X)) then
      call update_halo_5(4, X)
    end if
    if (present(Phi0_orig)) then
      call update_halo_6(4, 25, Phi0_orig)
    end if 
#endif

    ! initialise common variables
    beta = 0.4
    am3 = 1.0
    ibound = -1

    call init_gammas()
  end subroutine generate_starting_state

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
