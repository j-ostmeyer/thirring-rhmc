module remez_common_subroutines
  use params
  implicit none

contains
  subroutine read_remez_1(fileid, maxidx, arrnum, arrden)
    implicit none
    integer, intent(in) :: fileid
    integer, intent(in) :: maxidx
    real(dp), intent(out) :: arrnum(0:maxidx), arrden(1:maxidx)

    integer :: i

    read (fileid, *) arrnum(0)
    do i = 1, maxidx
      read (fileid, *) arrnum(i), arrden(i)
    enddo
  end subroutine read_remez_1

  subroutine read_remez_file(filename, maxidx, arrAnum, arrBnum, arrAden, arrBden)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: maxidx
    real(dp), intent(out) :: arrAnum(0:maxidx), arrBnum(0:maxidx)
    real(dp), intent(out) :: arrAden(maxidx), arrBden(maxidx)

    integer, parameter :: fileid = 36

    open (unit=fileid, file=filename, status='old', action='read')

    call read_remez_1(fileid=fileid, maxidx=maxidx, arrnum=arrAnum, arrden=arrAden)
    call read_remez_1(fileid=fileid, maxidx=maxidx, arrnum=arrBnum, arrden=arrBden)

    close (unit=fileid)
  end subroutine read_remez_file

end module remez_common_subroutines
