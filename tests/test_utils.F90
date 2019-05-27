module test_utils
  use params
  use random
#ifdef MPI
  use comms
#endif
  implicit none

contains

#ifdef MPI

  subroutine gdbwait()
    use comms
    implicit none
    logical,volatile :: wfi ! Wait For Intervention
    integer :: ierr

    wfi = .true.

    if(ip_global.eq.0)then 
      do while(wfi)
        call sleep(1)
      enddo
    endif
    call MPI_Barrier(MPI_Comm_World,ierr)
  end subroutine gdbwait

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
    end if

    ! Work out how many elements to transfer
    count = 1
    do size_index = 1, rank
      count = count * local_size(size_index)
    end do

    global_size(1 + offset) = ksize
    global_size(2 + offset) = ksize
    global_size(3 + offset) = ksizet

    start(1 + offset) = ip_x * ksizex_l
    start(2 + offset) = ip_y * ksizey_l
    start(3 + offset) = ip_t * ksizet_l

    call MPI_Type_Create_Subarray(rank, global_size, local_size, &
      start, MPI_Order_Fortran, mpi_dtype, local_mpiio_type,ierr)
    call MPI_Type_Commit(local_mpiio_type,ierr)   
    call MPI_File_Open(comm, filename, mode, &
      MPI_Info_Null, mpi_fh,ierr)
    call MPI_File_Set_View(mpi_fh, 4_8, mpi_dtype, &
      local_mpiio_type, "native", MPI_Info_Null,ierr)
    if (write_out) then
      call MPI_File_Write_All(mpi_fh, array, count, mpi_dtype, status,ierr)
    else
      call MPI_File_Read_All(mpi_fh, array, count, mpi_dtype, status,ierr)
    end if
    call MPI_File_Close(mpi_fh,ierr)
    call MPI_Type_Free(local_mpiio_type,ierr)
    return
  end subroutine rw_file_mpi
#endif
end module
