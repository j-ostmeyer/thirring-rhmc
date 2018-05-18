module inverter_checks
    use params
implicit none
    complex(dp) :: xout(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

contains

!   For the multi-shift inverter
subroutine dirac_op_shifted(xout,xin,am,imass,shift,num)
    use dwf3d_lib
    complex(dp),intent(in) :: xin(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp),intent(out) :: xout(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real, intent(in) :: am,shift,num
    integer, intent(in) :: imass

    call dirac_operator(xout,xin,am,imass)
    xout = (xout + shift*xin)!/num
end subroutine 

subroutine dirac_operator(xout,xin,am,imass)
    use trial, only: u
    use dwf3d_lib
    use comms
    complex(dp),intent(in) :: xin(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp),intent(out) :: xout(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real, intent(in) :: am
    integer, intent(in) :: imass

#ifdef MPI
    type(MPI_Request), dimension(12) :: reqs_xtemp, reqs_xout
#endif
    complex(dp) :: xtemp(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
 
    call dslash(xtemp,xin,u,am,imass)
#ifdef MPI
    call start_halo_update_5(4, xtemp, 8, reqs_xtemp)
    call complete_halo_update(reqs_xtemp)
#else
    call update_halo_5(4, xtemp)
#endif
    call dslashd(xout, xtemp, u, am, imass)
#ifdef MPI
    call start_halo_update_5(4, xout, 10, reqs_xout)
    call complete_halo_update(reqs_xout)
#else
    call update_halo_5(4, xout)
#endif

end subroutine dirac_operator 

subroutine check_diff(diffnorm2,x1,x2)
    use param
    use comms
    complex(dp),intent(in) :: x1(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp),intent(in) :: x2(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real(dp), intent(out) :: diffnorm2
    
    diffnorm2 = sum(abs(x1(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) - &
     &  x2(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))**2)

#ifdef MPI
    call MPI_AllReduce(MPI_In_Place, diffnorm2, 1, MPI_Real, MPI_Sum, comm)
#endif

end subroutine check_diff
end module
