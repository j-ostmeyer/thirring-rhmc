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
  !******************************************************************
  !    multisolver matrix inversion via Lanczos technique
  !  eg. Golub & van Loan "Matrix Computations" 9.3.1
  ! http://web.mit.edu/ehliu/Public/sclark/Golub%20G.H.,%20Van%20Loan%20C.F.-%20Matrix%20Computations.pdf
  !       solves (MdaggerM+diag)*x=Phi for ndiag different values of diag
  !   iflag=0: simply evaluates X = {MdaggerM}^p * Phi
  !   can be interchanged with congrad for p=-1
  !   iflag=1: in addition updates Phi0 register needed for PV force term
  !   iflag=2: evaluates DWF force term
  !   iflag=3: evaluates PV force term
  !*****************************************************************m
  subroutine qmrherm(Phi, res, itercg, am, imass, anum, aden, ndiagq, iflag, isweep, &
      & iter)
    use trial, only: u
    use vector
    use gforce
    use comms
    use dirac
    complex(dp), intent(in) :: Phi(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass, ndiagq, iflag, isweep, iter
    real(dp), intent(in) :: anum(0:ndiagq), aden(ndiagq)
    real, intent(in) :: res, am
    integer, intent(out) :: itercg
    !
    real(dp) :: alphatild
    real(dp) :: coeff
    !      
    real(dp) :: alpha(ndiagq)
    real(dp) :: amu(ndiagq), d(ndiagq), dm1(ndiagq)
    real(dp) :: rho(ndiagq), rhom1(ndiagq)
    real(dp) :: betaq, betaq0, phimod
    real :: resid, rhomax, arelax
    integer :: niter, idiag
    logical :: go_on
    #ifdef MPI
    type(MPI_Request), dimension(12) :: reqs_X2, reqs_vtild, reqs_Phi0, reqs_R, reqs_x
    integer :: ierr
    real(dp) :: dp_reduction ! DEBUG
    #endif
    !
    !     write(6,111)
    !111 format(' Hi from qmrherm')
    !
    resid=sqrt(kthird*ksize*ksize*ksizet*4*res*res)

    itercg=0
    !
    !   initialise r=Phi
    !
    R = Phi
    qm1 = cmplx(0.0, 0.0)
    x = anum(0) * Phi

    betaq = sum(abs(R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)) ** 2)
    #ifdef MPI
    !call MPI_AllReduce(MPI_In_Place, betaq, 1, MPI_Double_Precision, MPI_Sum, comm,ierr) ! DEBUG
    call MPI_AllReduce(betaq, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm,ierr) ! DEBUG
    betaq = dp_reduction


    ! Setting up persistent communication requests
    call init_halo_update_5(4, vtild, 1, reqs_vtild)
    call init_halo_update_5(4, R, 2, reqs_R)
    #endif 
    betaq = sqrt(betaq)
    phimod=betaq
    !
    !do niter=1,20
    niter = 0
    go_on = .true.
    do while(niter.lt.max_qmr_iters .and. go_on )
      niter=niter+1
      itercg=itercg+1
      !
      !  Lanczos steps
      !
      q = R / betaq

      call dslash(vtild,q,u,am,imass)
      #ifdef MPI
      ! No way to hide communications here unfortunately
      !call start_halo_update_5(4, vtild, 1, reqs_vtild)
      call MPI_Startall(12,reqs_vtild,ierr)
      !call complete_halo_update(reqs_vtild) ! Now this call happens in dslashd
      call dslashd(x3,vtild,u,am,imass,reqs_vtild)
      #else
      call update_halo_5(4, vtild)
      call dslashd(x3,vtild,u,am,imass)
      #endif

      !
      alphatild = sum(real(conjg(q(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)) & 
        &                * x3(:,1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)))
      #ifdef MPI
      !call MPI_AllReduce(MPI_In_Place, alphatild, 1, MPI_Double_Precision, MPI_Sum, comm,ierr)
      call MPI_AllReduce(alphatild, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm,ierr) ! DEBUG
      alphatild = dp_reduction ! DEBUG
      #endif
      !
      R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = &
        & x3(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) &
        & - alphatild * q(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) &
        & - betaq * qm1
      #ifdef MPI
      ! R will be needed at the start of the next iteration to compute q
      ! so start updating the boundary
      call MPI_Startall(12,reqs_R,ierr)
      !call start_halo_update_5(4, R, 2, reqs_R)
      #else
      call update_halo_5(4, R)
      #endif
      qm1 = q(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
      !
      betaq0 = betaq
      betaq = sum(abs(R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l,:)) ** 2)
      #ifdef MPI
      !call MPI_AllReduce(MPI_In_Place, betaq, 1, MPI_Double_Precision, MPI_Sum, comm,ierr)
      call MPI_AllReduce(betaq, dp_reduction, 1, MPI_Double_Precision, MPI_Sum, comm,ierr) ! DEBUG
      betaq = dp_reduction ! DEBUG
      #endif
      betaq = sqrt(betaq)
      !
      alpha = alphatild + aden
      !
      if(niter.eq.1)then
        d = alpha
        rho = betaq0 / alpha
        rhom1 = rho
        do idiag = 1, ndiagq
          p(:, :, :, :, :, idiag) = q(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :)
          x1(:, :, :, :, :, idiag) = rho(idiag) * p(:, :, :, :, :, idiag)
        enddo
        pm1 = p
      else
        amu = betaq0 / d
        dm1 = d
        d = alpha - betaq0 * amu
        rho = -amu * dm1 * rhom1 / d
        do idiag = 1, ndiagq
          p(:, :, :, :, :, idiag) = q(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) &
            & - amu(idiag) * pm1(:, :, :, :, :, idiag)
        enddo
        pm1 = p
        !     Convergence criterion (a bit ad hoc for now...)
        rhomax = real(maxval(abs(phimod * rho)))
        rhom1 = rho
        do idiag = 1, ndiagq
          x1(:, :, :, :, :, idiag) = &
            & x1(:, :, :, :, :, idiag) &
            & + rho(idiag) * p(:, :, :, :, :, idiag)
        enddo

        !     check to see whether the residual is acceptable for all ndiagq....
        !     criterion is a bit ad hoc -- relaxing by a factor arelax improves code
        !     stability and leads to quicker convergence
        arelax=2.0
        if(rhomax .lt. arelax * resid) then
          !     if(rhomax.lt.resid) then
          !     call testinv(Phi,resmax,itercg,am,imass,x1,aden,ndiagq)
          !     convergence based on || residual || not working well in single precision...
          !     if(resmax.lt.resid) goto 8
          go_on = .false.
        endif
      endif
      !     
      #ifdef MPI
      ! R will be needed at the start of the next iteration to compute q
      ! so start updating the bounddary
      call complete_halo_update(reqs_R)
      #endif
    enddo! do while(niter.lt.max_qmr_iters .and. go_on )



    if (niter.gt.max_qmr_iters) then
      #ifdef MPI
      if (ip_global .eq. 0) then
        #endif
        write(7,*) 'QMRniterc!, niter, isweep,iter,iflag,imass,anum,ndiagq = ', &
        &   niter,isweep, iter, iflag, imass, anum(0), ndiagq
        #ifdef MPI
      end if
      #endif
    endif
    !  
    8   continue
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
          & X1(:, :, :, :, :, 1:ndiagq)
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
        R(:, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :) = X1(:, :, :, :, :, idiag)
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
      write(6,*) "Qmrherm iterations,res:", itercg, res
    endif
    return
  end subroutine qmrherm
  !**********************************************************************
  !  iflag = 0 : evaluates Rdagger*(Mdagger)'*X2
  !  iflag = 1 : evaluates Rdagger*(M)'*X2
  !**********************************************************************
  subroutine derivs(R,X2,anum,iflag)
    use gforce
    use dirac
    !      complex, intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    !      complex, intent(in) :: X2(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: X2(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real(dp), intent(in) :: anum
    integer, intent(in) :: iflag

    !      complex(dp) :: tzi
    real :: tzi_real
    integer :: ix, iy, it, ixup, iyup, itup, idirac, mu
    integer :: igork1
    !
    !     write(6,111)
    !111 format(' Hi from derivs')

    !     dSdpi=dSdpi-Re(Rdagger *(d(Mdagger)dp)* X2)
    !     Cf. Montvay & Muenster (7.215)
    !      tzi=cmplx(0.0,2*anum)
    tzi_real = 2 * real(anum)
    !     factor of 2 picks up second term in M&M (7.215)
    !
    do mu = 1,3
      ixup = kdelta(1, mu)
      iyup = kdelta(2, mu)
      itup = kdelta(3, mu)

      do idirac=1,4
        do it = 1,ksizet_l
          do iy = 1,ksizey_l
            do ix = 1,ksizex_l
              dSdpi(ix,iy,it,mu) = &
                &     dSdpi(ix,iy,it,mu) + tzi_real * real(akappa) * sum(aimag( &
                &       conjg(R(:,ix,iy,it,idirac)) * &
                &         X2(:,ix+ixup,iy+iyup,it+itup,idirac)) &
                &     - aimag(conjg(R(:,ix+ixup,iy+iyup,it+itup,idirac)) * &
                &         X2(:,ix,iy,it,idirac)))
            enddo
          enddo
        enddo
        !
        igork1=gamin(mu,idirac)
        if(iflag.eq.0)then
          do it = 1,ksizet_l
            do iy = 1,ksizey_l
              do ix = 1,ksizex_l
                dSdpi(ix,iy,it,mu) = &
                  &     dSdpi(ix,iy,it,mu)+ tzi_real * sum(aimag(gamval(mu,idirac)* &
                  &(conjg(R(:,ix,iy,it,idirac)) * &
                  &        X2(:, ix+ixup,iy+iyup,it+itup,igork1) &
                  &+conjg(R(:,ix+ixup,iy+iyup,it+itup,idirac))* &
                  &             X2(:,ix,iy,it,igork1))))
              enddo
            enddo
          enddo
        else
          do it = 1,ksizet_l
            do iy = 1,ksizey_l
              do ix = 1,ksizex_l
                dSdpi(ix,iy,it,mu) = &
                  &     dSdpi(ix,iy,it,mu)- tzi_real * sum(aimag(gamval(mu,idirac)* &
                  &(conjg(R(:,ix,iy,it,idirac)) * &
                  &        X2(:,ix+ixup,iy+iyup,it+itup,igork1) &
                  &+conjg(R(:,ix+ixup,iy+iyup,it+itup,idirac)) * &
                  &             X2(:,ix,iy,it,igork1))))
              enddo
            enddo
          enddo
        endif
        !
      enddo
    enddo
    !
    return
  end subroutine derivs

end module qmrherm_module
