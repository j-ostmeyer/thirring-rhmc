module comms
  use comms_common
  use params
  use mpi
  implicit none
contains
  subroutine init_halo_types()
    use comms5, only: init_halo_types_5
    use comms5_sp, only: init_halo_types_5_sp
    use comms4, only: init_halo_types_4, init_halo_types_4_real
    use comms4_sp, only: init_halo_types_4_sp, init_halo_types_4_real_sp
    use comms6, only: init_halo_types_6

    call init_halo_types_4
    call init_halo_types_4_sp
    call init_halo_types_4_real
    call init_halo_types_4_real_sp
    call init_halo_types_5
    call init_halo_types_5_sp
    call init_halo_types_6
  end subroutine init_halo_types

  subroutine complete_halo_update(reqs)
    integer, intent(inout) :: reqs(12)
    integer :: ierr
    call MPI_Waitall(12, reqs, MPI_Statuses_Ignore, ierr)
  end subroutine complete_halo_update

  !***********************************************************************
  !   Initialise MPI variables
  !***********************************************************************
  subroutine init_MPI()
    integer :: coords(3)
#ifdef WITH_MUST
    integer, parameter :: must_rank = 1
#else
    integer, parameter :: must_rank = 0
#endif
    integer :: ierr

    call MPI_init(ierr)

    ! Check that we have the right number of processes
    call MPI_comm_size(MPI_COMM_WORLD, np_global, ierr)
    call MPI_comm_rank(MPI_COMM_WORLD, ip_global, ierr)
    if (np_global .ne. NP_X*NP_Y*NP_T + must_rank) then
      if(ip_global.eq.0) then
        print *, "MPI dimensionality mismatch: ", NP_X, "*", NP_Y, "*", NP_T, "!=", np_global
      end if
      call MPI_finalize(ierr)
      call exit(2)
    end if
    if(ip_global.eq.0) then
      print *, "Initialising MPI with grid ", NP_X, "*", NP_Y, "*", NP_T
    end if
 
    ! Set up a Cartesian communicator; periodic boundaries, allow reordering
    call MPI_cart_create(MPI_COMM_WORLD, 3, (/NP_X, NP_Y, NP_T/), &
      & (/.true., .true., .true./), .true., comm, ierr)

    ! Know where I am
    call MPI_cart_coords(comm, ip_global, 3, coords, ierr)
    ip_x = coords(1)
    ip_y = coords(2)
    ip_t = coords(3)

    ! Prepare file format for MPI-IO
    call MPI_Type_Create_Subarray(4, &! dimensionality
                                  (/ksize, ksize, ksizet, 3/), &! global volume
                                  (/ksizex_l, ksizey_l, ksizet_l, 3/), &! local volume
                                  (/ip_x*ksizex_l, ip_y*ksizey_l, ip_t*ksizet_l, 0/), &! start location
                                  MPI_Order_Fortran, &! array ordering
                                  MPI_Real, &! datatype to store
                                  mpiio_type, ierr) ! type descriptor for this subarray type
    call MPI_Type_Commit(mpiio_type, ierr)

    ! Prepare all of the halo types
    call init_halo_types
    return
  end subroutine init_MPI
end module comms
