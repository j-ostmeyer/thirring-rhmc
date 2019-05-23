#ifdef SCOREPINST
#include "scorep/SCOREP_User.inc"
#endif

module qmrherm_module
  use params
  implicit none
  complex(dp) :: vtild(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
  complex(dp) :: q(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
  complex(dp) :: pm1(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
  complex(dp) :: qm1(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
  complex(dp) :: p(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
  complex(dp) :: x3(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
  complex(dp) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
  complex(dp) :: x1(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
  complex(dp) :: x2(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

  complex(dp),save :: Phi0(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4,ndiag)
  logical :: printall

contains

  ! From https://arxiv.org/abs/hep-lat/9612014
  ! Krylov space solvers for shifted linear systems, B. Jegerlehner, 1996
  subroutine multishift_solver(u,am,imass,ndiagq,aden,output,input,res,&
    &maxcg,cg_return)
    use dirac
    use params 
    use comms
    use inverter_utils
    ! subroutine parameters
    complex(dp),intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    real, intent(in) :: am
    integer, intent(in) :: imass
    integer, intent(in) :: ndiagq
    real(dp), intent(in) :: aden(ndiagq)
    complex(dp) :: output(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)
    complex(dp) :: input(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real, intent(in) :: res
    integer, intent(in) :: maxcg
    integer, intent(out) :: cg_return

    
    ! temporary variables - large vectors
    complex(dp) ::         r(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) ::         h(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) ::         s(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) ::         p(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: shiftferm(kthird, ksizex_l, ksizey_l, ksizet_l, 4, ndiag)

    ! temporary variables - short vectors
    real(dp) :: zeta_i(ndiag)
    real(dp) :: zeta_ii(ndiag)
    real(dp) :: zeta_iii(ndiag)
    real(dp) :: omegas(ndiag)
    real(dp) :: gammas(ndiag)
    logical :: flags(ndiag) ! flag(i) = .true. : output(i) still has not converged.

    ! temporary variables - scalars
    real(dp) :: alpha,delta,lambda,omega,omega_save,gammag,fact

    ! temporary variables - utilities
    integer :: ishift,maxishift
    real(dp) :: source_norm
#ifdef MPI
    integer, dimension(12) :: reqs_h,reqs_p
    integer :: ierr
    real(dp) :: dp_reduction ! DEBUG
#endif

    r = input ! vector
    p = r     ! vector
    delta = sum(abs(r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l,:))**2)
#ifdef MPI
    call MPI_AllReduce(delta, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm,ierr)
    delta = dp_reduction
#endif

    source_norm = delta
    omega = 1.0d0
    alpha = 0.0d0
    output = 0.0 ! vector
    ! Setting up persistent communication requests
    call init_halo_update_5(4, h, 1, reqs_h)
    call init_halo_update_5(4, p, 2, reqs_p)

    do ishift=1,ndiagq
      flags(ishift) = .true.
      shiftferm(:,:,:,:,:,ishift) = input(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l,:)
      zeta_i(ishift) = 1.0d0
      zeta_ii(ishift) = 1.0d0
      gammas(ishift) = 0.0d0
    enddo
    gammag = 0.0d0
    cg_return = 0
    maxishift = ndiagq

    do while((maxishift .gt. 0).and.(cg_return.lt.maxcg))
      cg_return = cg_return + 1 


      ! DIRAC OPERATOR
      call dirac_operator(s,p,u,am,imass)
!      call dslash(h,p,u,am,imass)
!#ifdef MPI
!      call MPI_Startall(12,reqs_h,ierr)
!      !call complete_halo_update(reqs_h) ! Now this call happens in dslashd
!      call dslashd(s,h,u,am,imass,reqs_h)
!#else
!      call update_halo_5(4, h)
!      call dslashd(s,h,u,am,imass)
!#endif

      alpha = sum(real(conjg(p(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)) & 
        &                * s(:,1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)))

#ifdef MPI
      call MPI_AllReduce(alpha, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm,ierr)
      alpha = dp_reduction
#endif
      
      omega_save = omega
      omega = -delta/alpha

      do ishift=1,maxishift
          zeta_iii(ishift) = (zeta_i(ishift)*zeta_ii(ishift)*omega_save)/ &      
            & ( omega*gammag*(zeta_i(ishift)-zeta_ii(ishift))+ &
            &   zeta_i(ishift)*omega_save*(1.0-aden(ishift)*omega) ) ! TO CHECK definition of aden
          omegas(ishift)=omega*zeta_iii(ishift)/zeta_ii(ishift)   
      enddo

      do ishift=1,maxishift
        output(:,:,:,:,:,ishift) = output(:,:,:,:,:,ishift)-&
          &shiftferm(:,:,:,:,:,ishift)*omegas(ishift) 
      enddo

      r = r + omega*s
      lambda = sum(abs(r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l,:))**2)
#ifdef MPI
      call MPI_AllReduce(lambda, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm,ierr)
      lambda = dp_reduction
#endif
      gammag = lambda/delta

      p = r + gammag*p
      call MPI_Startall(12,reqs_p,ierr)

      gammas(1:ishift)=gammag*zeta_iii(1:ishift)*omegas(1:ishift)/(zeta_ii(1:ishift)*omega)

      do ishift=1,maxishift
        shiftferm(:,:,:,:,:,ishift) = shiftferm(:,:,:,:,:,ishift) * gammas(ishift)+&
         & r(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l,:) * zeta_iii(ishift)
      enddo

!      test : block
!        ! CHECK : Real residual seems to decrease faster than the one computed...
!        complex(dp) :: xout(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
!        complex(dp) :: xin(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
!        complex(dp) :: check(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
!#ifdef MPI
!        integer, dimension(12) :: reqs_xin
!#endif
!
!        real(dp) :: checksum
!
!        if(mod(cg_return,10)==0) then
!          do ishift =1,ndiagq
!            xin(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = output(:,:,:,:,:,ishift)
!#ifdef MPI
!            call start_halo_update_5(4, xin, 10, reqs_xin)
!            call complete_halo_update(reqs_xin)
!#else
!            call update_halo_5(4, xin)
!#endif
!            call dirac_op_shifted(xout,xin,u,am,imass,real(aden(ishift)))
!
!            check = input(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) - &
!              & xout(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
!            checksum = sum(abs(check)**2)
!#ifdef MPI
!            call MPI_AllReduce(checksum, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm,ierr)
!            checksum = dp_reduction
!#endif
!            if(ip_global.eq.0) then
!              print*, ishift,sqrt(delta*zeta_i(ishift)**2),checksum,res,real(aden(ishift))
!            endif
!          enddo
!        endif
!      end block test
 
      maxishift = 0
      do ishift=1,ndiagq
        if(flags(ishift))then
          if(sqrt(delta*zeta_ii(ishift)**2) < res) then
            flags(ishift) = .false.
          else
            maxishift = ishift 
          endif
        endif
      enddo
      zeta_i(1:maxishift)  = zeta_ii(1:maxishift)
      zeta_ii(1:maxishift) = zeta_iii(1:maxishift)

      delta = lambda

      call complete_halo_update(reqs_p)

      !if(ip_global.eq.0) then
      !   print*, "iteration", cg_return, "of", maxcg
      !endif
    enddo


  end subroutine multishift_solver


  subroutine qmrherm(Phi,X, res, itercg, am, imass, anum, aden, ndiagq, iflag, isweep, &
      & iter)
    use params
    use trial, only: u
    use gforce
    use comms
    use dirac
    use derivs_module
    complex(dp), intent(in) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(out) :: X(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass, ndiagq, iflag, isweep, iter
    real(dp), intent(in) :: anum(0:ndiagq), aden(ndiagq)
    real, intent(in) :: res, am
    integer, intent(out) :: itercg
    real(dp) :: coeff
    integer :: niter, idiag
#ifdef MPI
    integer, dimension(12) :: reqs_X2, reqs_Phi0, reqs_R, reqs_x
    integer :: ierr
#endif

    call multishift_solver(u,am,imass,ndiagq,aden,x1,Phi,res,max_qmr_iters,itercg)
 
    if(iflag.lt.2)then
      !     Now evaluate solution x=(MdaggerM)^p * Phi
      do idiag=1,ndiagq
        x(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = &
          & x(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) &
          & + anum(idiag) * x1(:, :, :, :, :, idiag)
      enddo
#ifdef MPI
      ! x is a saved module variable, so must be updated to avoid polluting the parent function
      ! could this in principle be moved outside the function so we don't do it unnecessarily?
      ! but in that case we wouldn't be able to hide the communications
      call start_halo_update_5(4, x, 3, reqs_x)
#else
      call update_halo_5(4, x)
#endif
      !     
      !  update phi0 block if required...
      if(iflag.eq.1) then
        Phi0(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :, 1:ndiagq) = &
          & x1(:, :, :, :, :, 1:ndiagq)
#ifdef MPI
        ! No way to hide communications here unfortunately
        ! In principle this could be better interleaved with the x update
        ! but that would add extra branching, and this section is messy enough already
        call start_halo_update_6(4, ndiagq, Phi0, 4, reqs_Phi0)
        call complete_halo_update(reqs_Phi0)
#else
        call update_halo_6(4, ndiagq, Phi0)
#endif
      endif! if(iflag.eq.1) then
#ifdef MPI
      call complete_halo_update(reqs_x)
#endif
      !     
    else! if(iflag.lt.2)then
      !
      do idiag=1, ndiagq
        !
        !  X2 = M*X1
        R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = x1(:, :, :, :, :, idiag)
#ifdef MPI
        ! No way to hide communications here unfortunately
        call start_halo_update_5(4, R, 5, reqs_R)
        call complete_halo_update(reqs_R)
#else
        call update_halo_5(4, R)
#endif

        ! Communication of X2 generated here can be hidden if iflag isn't 2, while R is updated
        call dslash(X2, R, u, am, imass)
#ifdef MPI
        call start_halo_update_5(4, X2, 6, reqs_X2)
#else
        call update_halo_5(4, X2)
#endif
        !
        if(iflag.eq.2)then
          coeff=anum(idiag)
#ifdef MPI
          call complete_halo_update(reqs_X2)
#endif
          call derivs(R, X2, coeff, 0)
        else! if(iflag.eq.2)then
          coeff=-anum(idiag)
          R = Phi0(:, :, :, :, :, idiag)
#ifdef MPI
          call complete_halo_update(reqs_X2)
#endif
          call derivs(R, X2, coeff, 0)
          !
          ! Communication of X2 generated here can be hidden while R is updated
          call dslash(X2, R, u, am, imass)
#ifdef MPI
          call start_halo_update_5(4, X2, 7, reqs_X2)
#else
          call update_halo_5(4, X2)
#endif
          !
          R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = x1(: ,:, :, :, :, idiag)
#ifdef MPI
          call start_halo_update_5(4, R, 8, reqs_R)
          call complete_halo_update(reqs_X2)
          call complete_halo_update(reqs_R)
#else
          call update_halo_5(4, R)
#endif
          call derivs(X2, R, coeff, 1)
        endif! if(iflag.eq.2)then
      enddo! do idiag=1, ndiagq
    endif !if(iflag.lt.2)then , else

    if (ip_global .eq. 0 .and. printall) then
      print*, "Qmrherm iterations,res:", itercg, res
    endif
    return
  end subroutine qmrherm
  !**********************************************************************
  !  iflag = 0 : evaluates Rdagger*(Mdagger)'*X2
  !  iflag = 1 : evaluates Rdagger*(M)'*X2
  !**********************************************************************
end module qmrherm_module
