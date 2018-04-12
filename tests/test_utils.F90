module test_utils
  use params
  use random
  use gauge
#ifdef MPI
  use comms
#endif
  implicit none
  
#ifdef MPI
  contains

  subroutine rw_file_mpi(array, array_shape, rank, filename, mpi_dtype, write_out)
    integer, intent(in) :: array_shape(:), rank
    type(*), dimension(..), intent(inout) :: array
    type(MPI_Datatype), intent(in) :: mpi_dtype
    character(len=*), intent(in) :: filename
    logical, intent(in) :: write_out
    integer, dimension(rank) :: global_size, local_size, start  
    integer :: size_index
    integer :: count = 1
    integer :: offset = 1
    type(MPI_File) :: mpi_fh
    type(MPI_Datatype) :: local_mpiio_type
    type(MPI_Status) :: status
    integer :: mode

    if (write_out) then
       mode = MPI_Mode_Wronly + MPI_Mode_Create
    else
       mode = MPI_Mode_Rdonly
    end if
    
    global_size = array_shape
    local_size = array_shape
    start = 0

    if (rank .eq. 4) then
       offset = 0
    end if

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
         start, MPI_Order_Fortran, mpi_dtype, local_mpiio_type)
    call MPI_Type_Commit(local_mpiio_type)   
    call MPI_File_Open(comm, filename, mode, &
         MPI_Info_Null, mpi_fh)
    call MPI_File_Set_View(mpi_fh, 4_8, mpi_dtype, &
        local_mpiio_type, "native", MPI_Info_Null)
    if (write_out) then
       call MPI_File_Write_All(mpi_fh, array, count, mpi_dtype, status)
    else
       call MPI_File_Read_All(mpi_fh, array, count, mpi_dtype, status)
    end if
    call MPI_File_Close(mpi_fh)
    return
  end subroutine rw_file_mpi
#endif


end module
