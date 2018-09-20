#include "test_utils.fh"
program test_partitioning
  use mpi_f08
  use partitioning
  implicit none
  integer :: ierr
  integer :: ips(3)
  integer :: issindex

#ifdef MPI
  call init_MPI
#endif
  call init_partitions_and_neighs()

  issindex = 1
  if(border_cl(0,0,0).ne.0) then
    print *,"Issue ", issindex
  endif
  issindex = issindex + 1

  if(.not.check_hbass((/-2,0,0/),(/1,0,0/)))then 
    print *,"Issue ", issindex
  endif
  issindex = issindex + 1

  if(.not.check_hbass((/-2,1,1/),(/1,1,1/)))then 
  endif
    print *,"Issue ", issindex
  issindex = issindex + 1

  if(.not.check_hbass((/1,-2,1/),(/1,1,1/)))then 
  endif
    print *,"Issue ", issindex
  issindex = issindex + 1

  if(.not.check_hbass((/1,1,-2/),(/1,1,1/)))then 
  endif
    print *,"Issue ", issindex
  issindex = issindex + 1

#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program


function check_hbass(hips,bips) result(check)
  use partitioning
  integer, intent(in) :: hips(3)
  integer, intent(in) :: bips(3)
  logical,            :: check

  integer :: thips(3)
  integer :: ibp,ihp
  integer :: nhp, idhp

  ibp = border_cl(bips(1),bips(2),bips(3))
  nhp = bhass(ibp)%nhp
  check = .false.
  do idhp=1,nhp 
    ihp = bhass(ibp)%hps(idhp)
    thips = halo_cl(:,ihp)
    if(all(thips.eq.hips)) then
      check = .true.
    endif
  enddo
end function


