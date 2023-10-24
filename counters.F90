module counters
  implicit none
  save

  real :: ancg, ancgh, ancgf, ancgpf
  real :: ancgpv, ancghpv, ancgfpv, ancgpfpv
  real :: action_average, vel2a, pbp, pbp_average, y_average, ysq_average
  real :: ancgm, ancgm_average
  real :: ancgmS, ancgmS_average,pbpS
  integer :: naccp, ipbp, itot
contains

  subroutine init_counters()
    implicit none

    action_average = 0.0
    vel2a = 0.0
    pbp_average = 0.0

    ancg = 0.0
    ancgf = 0.0
    ancgfpv = 0.0
    ancgh = 0.0
    ancghpv = 0.0
    ancgpf = 0.0
    ancgpfpv = 0.0
    ancgpv = 0.0

    ancgm_average = 0.0
    y_average = 0.0
    ysq_average = 0.0
    naccp = 0
    ipbp = 0
    itot = 0

  end subroutine init_counters

  subroutine final_averages(Nf, iter2)
    implicit none
    integer, intent(in) :: Nf, iter2

    action_average = action_average/iter2
    vel2a = vel2a/iter2
    if (ipbp .ne. 0) then
      pbp_average = pbp_average/ipbp
      ancgm_average = ancgm_average/ipbp
    else
      pbp_average = 0
      ancgm_average = 0
    endif
    ancg = ancg/(Nf*itot)
    ancgpv = ancgpv/(Nf*itot)
    ancgf = ancgf/(Nf*itot)
    ancgfpv = ancgfpv/(Nf*itot)

    ancgh = ancgh/(2*Nf*iter2)
    ancgpf = ancgpf/(Nf*iter2)
    ancghpv = ancghpv/(2*Nf*iter2)
    ancgpfpv = ancgpfpv/(Nf*iter2)
    y_average = y_average/iter2
    ysq_average = ysq_average/iter2 - y_average*y_average
    if (iter2 .gt. 1) then
      ysq_average = sqrt(ysq_average/(iter2 - 1))
    else
      ysq_average = 0.0
    endif

  end subroutine final_averages

end module counters
