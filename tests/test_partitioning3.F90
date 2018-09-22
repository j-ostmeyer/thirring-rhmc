! Testing associations between border and halo partitions
! (member ahpss in global structure arrays border_partition_{list,cube}
#include "test_utils.fh"
program test_partitioning
  use partitioning
  use comms
  implicit none
  integer :: ierr
  integer :: issindex
  logical :: check_ahpss
  logical :: check_ahpsr

#ifdef MPI
  call init_MPI
#endif
  call init_partitioning()

  issindex = 1
  if(border_cl(0,0,0).ne.0) then
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpss((/-2,0,0/),(/1,0,0/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpss((/-2,1,1/),(/1,1,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpss((/1,-2,1/),(/1,1,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpss((/1,1,-2/),(/1,1,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpss((/1,-1,-2/),(/1,-1,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpss((/1,2,1/),(/1,-1,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpss((/-2,-1,1/),(/1,-1,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpss((/-2,0,1/),(/1,0,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpss((/1,0,-2/),(/1,0,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1


  if(.not.check_ahpsr((/2,0,0/),(/1,0,0/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpsr((/2,1,1/),(/1,1,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpsr((/1,2,1/),(/1,1,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpsr((/1,1,2/),(/1,1,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpsr((/1,-1,2/),(/1,-1,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpsr((/1,-2,1/),(/1,-1,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpsr((/2,-1,1/),(/1,-1,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpsr((/2,0,1/),(/1,0,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1

  if(.not.check_ahpsr((/1,0,2/),(/1,0,1/)))then 
    if(ip_global.eq.0) then
      print *,"Issue ", issindex
    endif
  endif
  issindex = issindex + 1


#ifdef MPI
  call MPI_Finalize(ierr)
#endif
end program


function check_ahpss(hips,bips) result(check)
  use partitioning
  use comms
  integer, intent(in) :: hips(3)
  integer, intent(in) :: bips(3)
  logical             :: check

  integer :: thips(3)
  integer :: ibp,ihp
  integer :: nhp, idhp

  ibp = border_cl(bips(1),bips(2),bips(3))
  nhp = border_partitions_list(ibp)%nn
  check = .false.
  do idhp=1,nhp ! scanning through associations
    ihp = border_partitions_list(ibp)%ahpss(idhp)
    thips = halo_lc(:,ihp)
    if(all(thips.eq.hips)) then
      check = .true.
    endif
  enddo
end function

function check_ahpsr(hips,bips) result(check)
  use partitioning
  use comms
  integer, intent(in) :: hips(3)
  integer, intent(in) :: bips(3)
  logical             :: check

  integer :: thips(3)
  integer :: ibp,ihp
  integer :: idir

  ibp = border_cl(bips(1),bips(2),bips(3))
  check = .false.
  do idir=-3,3
    ihp = border_partitions_list(ibp)%ahpsr(idir)
    if(ihp.ne.0)then
      thips = halo_lc(:,ihp)
      if(all(thips.eq.hips)) then
        check = .true.
      endif
    endif
  enddo
end function


