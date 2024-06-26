module comms5
  use comms_common
  use params
  use mpi

  implicit none

  integer :: halo_5_xup_send(4:4), halo_5_xdn_send(4:4)
  integer :: halo_5_xup_recv(4:4), halo_5_xdn_recv(4:4)
  integer :: halo_5_yup_send(4:4), halo_5_ydn_send(4:4)
  integer :: halo_5_yup_recv(4:4), halo_5_ydn_recv(4:4)
  integer :: halo_5_tup_send(4:4), halo_5_tdn_send(4:4)
  integer :: halo_5_tup_recv(4:4), halo_5_tdn_recv(4:4)
  integer :: halo_5_thirdup_send(4:4), halo_5_thirddn_send(4:4)
  integer :: halo_5_thirdup_recv(4:4), halo_5_thirddn_recv(4:4)

contains

  ! Initialise a single halo type.
  ! direction: direction to be communicated; x=0, y=1, t=2, third=-1
  ! position: coordinate in direction of the region
  ! i5, i6: size in extra directions
  ! datatype: datatype to transfer
  ! typetarget: array element to put type into;
  !   should be halo_{4,4_real,5,6}_{x,y,t}{up,dn}_{send,recv}(size{4,5}[, size6])
  subroutine init_single_halo_type_5(direction, position, size5, typetarget)
    implicit none
    integer, intent(in) :: direction, position, size5
    integer, intent(out) :: typetarget
    integer, dimension(5) :: sizes, subsizes, starts
    integer :: ierr

    sizes = (/kthird_l + 2, ksizex_l + 2, ksizey_l + 2, ksizet_l + 2, size5/)
    subsizes = (/kthird_l, ksizex_l, ksizey_l, ksizet_l, size5/)
    subsizes(direction + 2) = 1
    starts = (/1, 1, 1, 1, 0/)
    starts(direction + 2) = position

    call MPI_Type_Create_Subarray(5, sizes, subsizes, starts, MPI_Order_Fortran, &
                                  MPI_Double_Complex, typetarget, ierr)
    call MPI_Type_Commit(typetarget, ierr)
    return
  end subroutine init_single_halo_type_5

  subroutine init_halo_types_5()
    implicit none
    integer, parameter :: nsize5 = 1
    integer :: size5(nsize5) = (/4/)
    integer :: i5

    do i5 = 1, nsize5
      call init_single_halo_type_5(0, 1, size5(i5), halo_5_xdn_send(size5(i5)))
      call init_single_halo_type_5(0, 0, size5(i5), halo_5_xdn_recv(size5(i5)))
      call init_single_halo_type_5(0, ksizex_l, size5(i5), halo_5_xup_send(size5(i5)))
      call init_single_halo_type_5(0, ksizex_l + 1, size5(i5), halo_5_xup_recv(size5(i5)))

      call init_single_halo_type_5(1, 1, size5(i5), halo_5_ydn_send(size5(i5)))
      call init_single_halo_type_5(1, 0, size5(i5), halo_5_ydn_recv(size5(i5)))
      call init_single_halo_type_5(1, ksizey_l, size5(i5), halo_5_yup_send(size5(i5)))
      call init_single_halo_type_5(1, ksizey_l + 1, size5(i5), halo_5_yup_recv(size5(i5)))

      call init_single_halo_type_5(2, 1, size5(i5), halo_5_tdn_send(size5(i5)))
      call init_single_halo_type_5(2, 0, size5(i5), halo_5_tdn_recv(size5(i5)))
      call init_single_halo_type_5(2, ksizet_l, size5(i5), halo_5_tup_send(size5(i5)))
      call init_single_halo_type_5(2, ksizet_l + 1, size5(i5), halo_5_tup_recv(size5(i5)))

      call init_single_halo_type_5(-1, 1, size5(i5), halo_5_thirddn_send(size5(i5)))
      call init_single_halo_type_5(-1, 0, size5(i5), halo_5_thirddn_recv(size5(i5)))
      call init_single_halo_type_5(-1, kthird_l, size5(i5), halo_5_thirdup_send(size5(i5)))
      call init_single_halo_type_5(-1, kthird_l + 1, size5(i5), halo_5_thirdup_recv(size5(i5)))
    end do
  end subroutine init_halo_types_5

  subroutine init_halo_update_5(size5, Array, tag, reqs)
    implicit none
    integer, intent(in) :: size5, tag
    complex(dp), intent(inout) :: Array(0:kthird_l + 1, &
                                        0:ksizex_l + 1, &
                                        0:ksizey_l + 1, &
                                        0:ksizet_l + 1, &
                                        size5)
    integer, intent(out) :: reqs(16)
    integer :: ip_xup, ip_xdn, ip_yup, ip_ydn, ip_tup, ip_tdn, ip_thirdup, ip_thirddn
    integer :: tag_offset
    integer :: ierr

    tag_offset = 8*tag

    ! Start send and receive in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup, ierr)
    call MPI_Send_init(Array, 1, halo_5_xup_send(size5), ip_xup, 0 + tag_offset, &
                       comm, reqs(1), ierr)
    call MPI_Recv_init(Array, 1, halo_5_xdn_recv(size5), ip_xdn, 0 + tag_offset, &
                       comm, reqs(2), ierr)
    call MPI_Send_init(Array, 1, halo_5_xdn_send(size5), ip_xdn, 1 + tag_offset, &
                       comm, reqs(3), ierr)
    call MPI_Recv_init(Array, 1, halo_5_xup_recv(size5), ip_xup, 1 + tag_offset, &
                       comm, reqs(4), ierr)

    ! Start send and receive in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup, ierr)
    call MPI_Send_init(Array, 1, halo_5_yup_send(size5), ip_yup, 2 + tag_offset, &
                       comm, reqs(5), ierr)
    call MPI_Recv_init(Array, 1, halo_5_ydn_recv(size5), ip_ydn, 2 + tag_offset, &
                       comm, reqs(6), ierr)
    call MPI_Send_init(Array, 1, halo_5_ydn_send(size5), ip_ydn, 3 + tag_offset, &
                       comm, reqs(7), ierr)
    call MPI_Recv_init(Array, 1, halo_5_yup_recv(size5), ip_yup, 3 + tag_offset, &
                       comm, reqs(8), ierr)

    ! Start send and receive in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup, ierr)
    call MPI_Send_init(Array, 1, halo_5_tup_send(size5), ip_tup, 4 + tag_offset, &
                       comm, reqs(9), ierr)
    call MPI_Recv_init(Array, 1, halo_5_tdn_recv(size5), ip_tdn, 4 + tag_offset, &
                       comm, reqs(10), ierr)
    call MPI_Send_init(Array, 1, halo_5_tdn_send(size5), ip_tdn, 5 + tag_offset, &
                       comm, reqs(11), ierr)
    call MPI_Recv_init(Array, 1, halo_5_tup_recv(size5), ip_tup, 5 + tag_offset, &
                       comm, reqs(12), ierr)

    ! Start send and receive in "third" direction
    call MPI_Cart_Shift(comm, 3, 1, ip_thirddn, ip_thirdup, ierr)
    call MPI_Send_init(Array, 1, halo_5_thirdup_send(size5), ip_thirdup, 6 + tag_offset, &
                       comm, reqs(13), ierr)
    call MPI_Recv_init(Array, 1, halo_5_thirddn_recv(size5), ip_thirddn, 6 + tag_offset, &
                       comm, reqs(14), ierr)
    call MPI_Send_init(Array, 1, halo_5_thirddn_send(size5), ip_thirddn, 7 + tag_offset, &
                       comm, reqs(15), ierr)
    call MPI_Recv_init(Array, 1, halo_5_thirdup_recv(size5), ip_thirdup, 7 + tag_offset, &
                       comm, reqs(16), ierr)

  end subroutine init_halo_update_5

  subroutine start_halo_update_5(size5, Array, tag, reqs)
    implicit none
    integer, intent(in) :: size5, tag
    complex(dp), intent(inout) :: Array(0:kthird_l + 1, &
                                        0:ksizex_l + 1, &
                                        0:ksizey_l + 1, &
                                        0:ksizet_l + 1, &
                                        size5)
    integer, intent(out) :: reqs(16)
    integer :: ip_xup, ip_xdn, ip_yup, ip_ydn, ip_tup, ip_tdn, ip_thirdup, ip_thirddn
    integer :: tag_offset
    integer :: ierr

    tag_offset = 8*tag

    ! Start send and receive in x direction
    call MPI_Cart_Shift(comm, 0, 1, ip_xdn, ip_xup, ierr)
    call MPI_Isend(Array, 1, halo_5_xup_send(size5), ip_xup, 0 + tag_offset, &
                   comm, reqs(1), ierr)
    call MPI_Irecv(Array, 1, halo_5_xdn_recv(size5), ip_xdn, 0 + tag_offset, &
                   comm, reqs(2), ierr)
    call MPI_Isend(Array, 1, halo_5_xdn_send(size5), ip_xdn, 1 + tag_offset, &
                   comm, reqs(3), ierr)
    call MPI_Irecv(Array, 1, halo_5_xup_recv(size5), ip_xup, 1 + tag_offset, &
                   comm, reqs(4), ierr)

    ! Start send and receive in y direction
    call MPI_Cart_Shift(comm, 1, 1, ip_ydn, ip_yup, ierr)
    call MPI_Isend(Array, 1, halo_5_yup_send(size5), ip_yup, 2 + tag_offset, &
                   comm, reqs(5), ierr)
    call MPI_Irecv(Array, 1, halo_5_ydn_recv(size5), ip_ydn, 2 + tag_offset, &
                   comm, reqs(6), ierr)
    call MPI_Isend(Array, 1, halo_5_ydn_send(size5), ip_ydn, 3 + tag_offset, &
                   comm, reqs(7), ierr)
    call MPI_Irecv(Array, 1, halo_5_yup_recv(size5), ip_yup, 3 + tag_offset, &
                   comm, reqs(8), ierr)

    ! Start send and receive in t direction
    call MPI_Cart_Shift(comm, 2, 1, ip_tdn, ip_tup, ierr)
    call MPI_Isend(Array, 1, halo_5_tup_send(size5), ip_tup, 4 + tag_offset, &
                   comm, reqs(9), ierr)
    call MPI_Irecv(Array, 1, halo_5_tdn_recv(size5), ip_tdn, 4 + tag_offset, &
                   comm, reqs(10), ierr)
    call MPI_Isend(Array, 1, halo_5_tdn_send(size5), ip_tdn, 5 + tag_offset, &
                   comm, reqs(11), ierr)
    call MPI_Irecv(Array, 1, halo_5_tup_recv(size5), ip_tup, 5 + tag_offset, &
                   comm, reqs(12), ierr)

    ! Start send and receive in "third" direction
    call MPI_Cart_Shift(comm, 3, 1, ip_thirddn, ip_thirdup, ierr)
    call MPI_Isend(Array, 1, halo_5_thirdup_send(size5), ip_thirdup, 6 + tag_offset, &
                   comm, reqs(13), ierr)
    call MPI_Irecv(Array, 1, halo_5_thirddn_recv(size5), ip_thirddn, 6 + tag_offset, &
                   comm, reqs(14), ierr)
    call MPI_Isend(Array, 1, halo_5_thirddn_send(size5), ip_thirddn, 7 + tag_offset, &
                   comm, reqs(15), ierr)
    call MPI_Irecv(Array, 1, halo_5_thirdup_recv(size5), ip_thirdup, 7 + tag_offset, &
                   comm, reqs(16), ierr)

  end subroutine start_halo_update_5

end module comms5
