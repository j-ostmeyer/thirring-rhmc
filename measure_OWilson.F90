module measure_OWilson ! note the imass=3 mass term has changed sign to be
                       ! consistent with imass=1
  use params
  use gammamatrices
  use dirac ! has Dwilson and DWilsonD
  use comms
  implicit none
  save
  real(dp) :: coeffs(kthird)
  logical,parameter :: VB_OW=.true.

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine odslash(Phi,R,u,am,imass)
    implicit none ! calculates Phi = M*R
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(inout) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Mslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    complex(dp) :: zkappa
    real :: diag
    integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork
    integer il

    ! diagonal
    do il=1,kthird
      Rslice=R(il,:,:,:,:)
      call DWilson(Phislice,Rslice,u,-am3)
      Phi(il,:,:,:,:)=coeffs(il)*Phislice+Rslice
    end do

      !  upper diagonal
      do il=1,kthird-1
        Rslice=cmplx(0.0,0.0)
        Rslice(:,:,:,3:4)=R(il+1,:,:,:,3:4)
        call DWilson(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &    + coeffs(il)*Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)
      end do

      ! lower diagonal
      do il=2,kthird
        Rslice=cmplx(0.0,0.0)
        Rslice(:,:,:,1:2)=R(il-1,:,:,:,1:2)
        call DWilson(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &  +   coeffs(il)*Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)   &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)
      end do
      
      !  Mass term - top right
      Rslice=cmplx(0.0,0.0)
      Rslice(:,:,:,1:2)=R(kthird,:,:,:,1:2)
      call DWilson(Phislice,Rslice,u,-am3)
      if (imass.eq.1) then
        Mslice(:,:,:,:) = -am*(coeffs(1)*Phislice-Rslice)
      elseif (imass.eq.3) then
!DB        Mslice(:,:,:,:) = -cmplx(0.0,am)*(coeffs(1)*Phislice-Rslice)
        Mslice(:,:,:,:) = cmplx(0.0,am)*(coeffs(1)*Phislice-Rslice)
      endif
      Phi(1,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
  &                Phi(1,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)  &
  & +             Mslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)

      !  Mass term - bottom left
      Rslice=cmplx(0.0,0.0)
      Rslice(:,:,:,3:4)=R(1,:,:,:,3:4)
      call DWilson(Phislice,Rslice,u,-am3)
      if (imass.eq.1) then
        Mslice(:,:,:,:)=-am*(coeffs(kthird)*Phislice-Rslice)
      elseif (imass.eq.3) then
!DB        Mslice(:,:,:,:)=cmplx(0.0,am)*(coeffs(kthird)*Phislice-Rslice)
        Mslice(:,:,:,:)=-cmplx(0.0,am)*(coeffs(kthird)*Phislice-Rslice)
      endif
      Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) = &
  &                 Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,:) &
  &            +        Mslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,:)


    !
    return
  end subroutine odslash
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine odslashd(Phi,R,u,am,imass,reqs_R)
    use comms, only : complete_halo_update
    implicit none
    complex(dp), intent(in) :: u(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 3)
    complex(dp), intent(inout) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp), intent(in) :: R(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: Rslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Phislice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: Mslice(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: am
    integer :: ixup, iyup, itup, ix, iy, it, idirac, mu, igork
    integer, dimension(16),intent(inout), optional :: reqs_R
    integer il
    complex(dp) :: zkappa
    real :: diag

    if(present(reqs_r)) then
      call complete_halo_update(reqs_R)
    endif

    do il=1,kthird
      Rslice=R(il,:,:,:,:)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Phi(il,:,:,:,:)=coeffs(il)*Phislice+Rslice
    enddo

      ! upper diagonal
      do il=1,kthird-1
        Rslice=R(il+1,:,:,:,:)
        call DWilsonD(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)   &
 &  +   coeffs(il+1)*Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)   &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)
      enddo

      ! lower diagonal
      do il=2,kthird
        Rslice=R(il-1,:,:,:,:)
        call DWilsonD(Phislice,Rslice,u,-am3)
        Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) = &
 &      Phi(il,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)   &
 &  +   coeffs(il-1)*Phislice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)  &
 &    - Rslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)
      enddo

      ! upper right mass term  
      Rslice=R(kthird,:,:,:,:)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Mslice=coeffs(kthird)*Phislice-Rslice
      if (imass.eq.1) then
        zkappa = -cmplx(am,0.0)
      elseif (imass.eq.3) then
!DB        zkappa = -cmplx(0.0,am)
        zkappa = cmplx(0.0,am)
      endif
      Phi(1,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4) = &
  &                Phi(1,1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)  &
  & +      zkappa*Mslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,3:4)

      ! lower left mass term  
      Rslice=R(1,:,:,:,:)
      call DWilsonD(Phislice,Rslice,u,-am3)
      Mslice=coeffs(1)*Phislice-Rslice
      if (imass.eq.1) then
        zkappa = -cmplx(am,0.0)
      elseif (imass.eq.3) then
!DB        zkappa = cmplx(0.0,am)
        zkappa = -cmplx(0.0,am)
      endif
      Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) = &
  &                 Phi(kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2) &
  &            + zkappa*Mslice(1:ksizex_l,1:ksizey_l,1:ksizet_l,1:2)

    return
  end subroutine odslashd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine prepZolo()
  use zolomodule
  use ratfuncs
  implicit none

  real(dp) lmin,lmax
  type(zolotarev) :: zolo
  type(sgnratfunc) :: srf
  integer i

  if (VB_OW) then ; if (ip_global.eq.0) then ; print *,"prep zolo" ; end if ; end if
  lmin=1d-1
  lmax=5
  call setZolo(lmin,lmax,kthird,zolo)
  call getRoots(zolo) 
  coeffs=1d0/zolo%roots
  if (VB_OW) then ; if (ip_global.eq.0) then ; print *,"zolo coeffs (omega):",coeffs ; end if ; end if

  return
  end subroutine prepZolo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine congradO(Phi, res, itercg, am, imass, iterations)
    use trial, only: u
    use vector
    use comms5, only: init_halo_update_5,start_halo_update_5
    use comms_common, only: comm
    use comms
    use params
    implicit none
    complex(dp), intent(in) :: Phi(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real, intent(in) :: res, am
    integer, intent(out) :: itercg
    integer, intent(in) :: imass
    integer, intent(out), optional :: iterations
    complex(dp) :: x1(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: x2(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: p(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    complex(dp) :: r(0:kthird_l + 1, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 4)
    real :: resid
    real(dp) :: betacg, betacgn, betacgd, alpha, alphan, alphad
    integer :: i
    integer, dimension(16) :: reqs_x1, reqs_r
    integer :: ierr
    logical,parameter :: VB_CG=.false.

    resid = 4*ksize*ksize*ksizet*kthird*res*res
    itercg = 0
    alphan = 0.0

!   initialise
    r = Phi 
    x=r 

    call odslash(x1, r, u, am, imass) 
    call start_halo_update_5(4, x1, 8, reqs_x1) 
    call complete_halo_update(reqs_x1)
    call odslashd(x2, x1, u, am, imass) 
    r=Phi-x2 
    call start_halo_update_5(4, r, 10, reqs_r) 
    call complete_halo_update(reqs_r)
    p=r
    betacgn=sum(conjg(r(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))*r(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))
    call MPI_AllReduce(MPI_In_Place, betacgn, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
    if (betacgn.lt.resid) goto 8
      
    do i=1,niterc
        itercg=itercg+1
        call odslash(x1, p, u, am, imass) 
        call start_halo_update_5(4, x1, 28, reqs_x1) 
!        call complete_halo_update(reqs_x1) ! completed in odslashd
        call odslashd(x2, x1, u, am, imass,reqs_x1)
        alphan=sum(conjg(r(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))*r(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))
        call MPI_AllReduce(MPI_In_Place, alphan, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
        alphad=sum(conjg(p(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))*x2(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))
        call MPI_AllReduce(MPI_In_Place, alphad, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
        alpha=alphan/alphad

        x=x+alpha*p
        r=r-alpha*x2
        call start_halo_update_5(4, r, 30, reqs_r) 
!        call complete_halo_update(reqs_r)
        betacgn=sum(conjg(r(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))*r(1:kthird_l, 1:ksizex_l, 1:ksizey_l, 1:ksizet_l, :))
        call MPI_AllReduce(MPI_In_Place, betacgn, 1, MPI_Double_Precision, MPI_Sum, comm, ierr)
        call complete_halo_update(reqs_r)

        if (VB_CG) then ; if (ip_global.eq.0) then ; print *,betacgn,resid ; end if ; endif

        if (betacgn.lt.resid) goto 8
        betacg=betacgn/alphan
        p=r+betacg*p
    end do
        print *,"WARNING: congrad did not converge"
8     continue

    if (present(iterations)) then
      iterations = i
    endif

    return
  end subroutine congradO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine IDDW(IDR,R,bmass,imass) ! IDR= IDDW.R
    use vector
    use trial, only: u
    use comms5, only: start_halo_update_5
    use comms
    implicit none
    complex(dp), intent(out) :: IDR(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real, intent(in) :: bmass
    integer, intent(in) :: imass
    integer :: itercg,iterations
    integer, dimension(16) :: reqsA
  
    call odslashd(IDR,R,u,bmass,imass) ! IDR = DDW^dag . R
    call start_halo_update_5(4, IDR, 21, reqsA) 
    call complete_halo_update(reqsA)
    call congradO(IDR, respbp, itercg, bmass, imass, iterations) ! X = (DDW)^-1 . (DDW^dag)^-1 . IDR
    IDR=X

  end subroutine IDDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine TEST_IDDW() ! check err=(1-IDDW.DDW).R
    use vector
    use trial, only: u
    use comms5, only: start_halo_update_5
    use comms
    implicit none
    complex(dp) :: DIFF(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
    complex(dp) :: DR(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: R(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: FIN(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real :: bmass,err
    integer :: imass
    integer, dimension(16) :: reqs
  
    bmass=0.01
    imass=1
    R=0
    R(1,:,:,:,:)=1.0

    call IDDW(DR,R,bmass,imass) 
    call odslash(FIN,DR,u,bmass,imass) 

    DIFF=FIN(1:kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)-R(1:kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)
    err=maxval(maxval(maxval(maxval(maxval(abs(DIFF),1),1),1),1),1)
    print *,ip_global,"TEST IDDW:",err
    return
  end subroutine TEST_IDDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine IDDWD(IDR,R,bmass,imass) ! IDR = IDDWdag.R = DDW . (DDW)^-1 . (DDW^dag)^-1 . Phi
    use vector
    use trial, only: u
    use comms5, only: start_halo_update_5
    use comms
    implicit none
!    complex(dp), intent(out) :: IDR(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
!    complex(dp), intent(in) :: R(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(out) :: IDR(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real, intent(in) :: bmass
    integer, intent(in) :: imass
    real :: res
    integer :: itercg,iterations
    integer, dimension(16) :: reqsA
  
    call congradO(R, res, itercg, bmass, imass, iterations) ! X = (DDW)^-1 . (DDW^dag)^-1 . Phi
    call odslash(IDR,X,u,bmass,imass) ! IDR = DDW . X
    call start_halo_update_5(4, IDR, 23, reqsA) 
    call complete_halo_update(reqsA)

  end subroutine IDDWD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine KDDW(DR,R,bmass,imass) ! DR = KDDW.R = Cdag . DDW^-1(1) . DDW(m) . C . R
    use trial, only: u                  
    use comms, only : complete_halo_update
    use comms5, only: start_halo_update_5
    implicit none
    complex(dp), intent(out) :: DR(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(in) :: R(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: TMP(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: bmass
    real :: one=1.0
    integer :: mtype
    integer, dimension(16) :: reqs
    integer l,d,jx,jy,jt

    mtype=imass
    call PermM(R,TMP,.false.) ! TMP = C.R
    call odslash(DR,TMP,u,bmass,mtype) ! DR = DDW . TMP
    call start_halo_update_5(4, DR, 53, reqs) 
    call complete_halo_update(reqs)
    mtype=1
    call IDDW(TMP,DR,one,mtype) ! TMP = IDDW . DR
    call PermM(TMP,DR,.true.) ! DR = Cdag . TMP

  return
  end subroutine KDDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine IKDDW(DR,R,bmass,imass) ! DR = IKDDW.R = Cdag . DDW^-1(m) . DDW(1) . C . R
    use trial, only: u                    
    use comms, only : complete_halo_update
    use comms5, only: start_halo_update_5
    implicit none
    complex(dp),intent(out) :: DR(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp),intent(in) :: R(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: TMP(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: bmass
    real :: one=1.0
    integer :: mtype
    integer, dimension(16) :: reqs

    mtype=1
    call PermM(R,TMP,.false.) ! TMP = C.R
    call odslash(DR,TMP,u,one,mtype) ! DR = DDW . TMP
    call start_halo_update_5(4, DR, 55, reqs) 
    call complete_halo_update(reqs)
    mtype=imass
    call IDDW(TMP,DR,bmass,mtype) ! TMP = IDDW . DR
    call PermM(TMP,DR,.true.) ! DR = Cdag . TMP

  return
  end subroutine IKDDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine TEST_IKDDW() !
    use vector
    use trial, only: u
    use comms5, only: start_halo_update_5
    use comms
    implicit none
    complex(dp) :: DIFF(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
    complex(dp) :: DR(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: R(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: FIN(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real :: bmass,err
    integer :: imass
    integer, dimension(16) :: reqs
  
    bmass=0.01
    imass=3
    R=1
    call KDDW(DR,R,bmass,imass) !
    call start_halo_update_5(4, DR, 21, reqs) 
    call complete_halo_update(reqs)
    call IKDDW(FIN,DR,bmass,imass) !
    
    DIFF=FIN(1:kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)-R(1:kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)
    err=maxval(maxval(maxval(maxval(maxval(abs(DIFF),1),1),1),1),1)
    print *,"Test IKDDW:",err

  end subroutine TEST_IKDDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine KDDWD(DR,R,bmass,imass) ! DR = KDDWdag.R = Cdag . DDWdag(m) . DDWdag^-1(1) . C . R
    use trial, only: u          
    use comms, only : complete_halo_update
    implicit none
    complex(dp),intent(out) :: DR(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp),intent(in) :: R(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: TMP(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer, intent(in) :: imass
    real, intent(in) :: bmass
    real :: one=1.0
    integer :: mtype

    call PermM(R,TMP,.false.) ! TMP = C.R
    mtype=1
    call IDDWD(DR,TMP,one,imass) ! DR = IDDWdag(1) . TMP
    mtype=imass
    call odslashd(TMP,DR,u,bmass,mtype) ! TMP = DDWdag(m) . DR
    call PermM(DR,TMP,.true.) ! DR = Cdag . TMP

  return
  end subroutine KDDWD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine KDDW4(DR,R,bmass,imass)
    implicit none
    complex(dp), intent(in) :: R(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(out) :: DR(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real,intent(in) :: bmass
    integer,intent(in) :: imass
    complex(dp) :: DR5(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: R5(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

    R5=(0.0,0.0)
    R5(1,:,:,:,:)=R
    call KDDW(DR5,R5,bmass,imass)
    DR=DR5(1,:,:,:,:)
      
    return
  end subroutine KDDW4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine IKDDW4(DR,R,bmass,imass)
    implicit none
    complex(dp), intent(in) :: R(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(out) :: DR(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real,intent(in) :: bmass
    integer,intent(in) :: imass
    complex(dp) :: DR5(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: R5(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer ix,iy,it,d

    R5=(0.0,0.0)
    R5(1,:,:,:,:)=R
    call IKDDW(DR5,R5,bmass,imass)
    DR=DR5(1,:,:,:,:)
      
    return
  end subroutine IKDDW4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine TEST_KDDW4() !
    use vector
    use trial, only: u
    use comms5, only: start_halo_update_5
    use comms
    implicit none
    complex(dp) :: DR(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: R(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real :: bmass,err
    integer :: imass
    integer d,jx,jy,jt

    bmass=0.00
    imass=1
    R=1
    print *,"R:"
    call PRINT_ARRAY(R) !
    call KDDW4(DR,R,bmass,imass) 
    print *,"KDDW4:"
    call PRINT_ARRAY(DR) !
    print *,"IKDDW4:"
    call IKDDW4(DR,R,bmass,imass) 
    call PRINT_ARRAY(DR) !
 
    return
  end subroutine TEST_KDDW4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine KDDWD4(DR,R,bmass,imass)
    implicit none
    complex(dp), intent(in) :: R(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp), intent(out) :: DR(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real,intent(in) :: bmass
    integer,intent(in) :: imass
    complex(dp) :: DR5(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: R5(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

    R5=(0.0,0.0)
    R5(1,:,:,:,:)=R
    call KDDWD(DR5,R5,bmass,imass)
    DR=DR5(1,:,:,:,:)
      
    return
  end subroutine KDDWD4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Pplus(R,PR) 
      implicit none
      complex(dp), intent(in) :: R(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp), intent(out) :: PR(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

      PR=0
      PR(:,:,:,1:2)=R(:,:,:,1:2)
      
      return
      end subroutine Pplus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Pminus(R,PR)
      implicit none
      complex(dp), intent(in) :: R(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp), intent(out) :: PR(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)

      PR=0
      PR(:,:,:,3:4)=R(:,:,:,3:4)

      return
      end subroutine Pminus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine PermM(R,PR,DAGGER) ! calculates DR(:,:,l) = Pminus.R(:,:,l)+Pplus.R(:,:,l+1)
      implicit none                 !     if DAGGER=.true. then
                                    ! calculates DR(:,:,l) = Pminus.R(:,:,l)+Pplus.R(:,:,l-1)
      complex(dp), intent(in) :: R(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp), intent(out) :: PR(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: TMP(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      logical DAGGER
      integer l

      do l=1,kthird
        call Pminus(R(l,:,:,:,:),PR(l,:,:,:,:))
      enddo
      do l=1,kthird
        call Pplus(R(l,:,:,:,:),TMP)
        if (.not.DAGGER) then
          if (l.eq.1) then
            PR(kthird,:,:,:,:)=PR(kthird,:,:,:,:)+TMP
          else
            PR(l-1,:,:,:,:)=PR(l-1,:,:,:,:)+TMP
          endif
        else
          if (l.eq.kthird) then
            PR(1,:,:,:,:)=PR(1,:,:,:,:)+TMP
          else
            PR(l+1,:,:,:,:)=PR(l+1,:,:,:,:)+TMP
          endif
        endif
      enddo

      return
      end subroutine PermM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine TEST_PermM() !
    use vector
    use trial, only: u
    use comms5, only: start_halo_update_5
    use comms
    implicit none
    complex(dp) :: DIFF(kthird, ksizex_l, ksizey_l, ksizet_l, 4)
    complex(dp) :: DR(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: R(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    complex(dp) :: FIN(0:kthird_l+1, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    real :: bmass,err
    integer :: imass
    integer, dimension(16) :: reqs
  
    bmass=0.01
    imass=3
    R=1
    call PermM(R,DR,.false.) 
    call PermM(DR,FIN,.true.) 
    
    DIFF=FIN(1:kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)-R(1:kthird,1:ksizex_l,1:ksizey_l,1:ksizet_l,:)
    err=maxval(maxval(maxval(maxval(maxval(abs(DIFF),1),1),1),1),1)
    print *,"PermM Test:",err

  end subroutine TEST_PermM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine measureKDDW4(bmass,imass,pbptot)
      use trial
      implicit none
      real,intent(in) :: bmass
      integer,intent(in) :: imass
      complex(dp),intent(out) :: pbptot
      complex(dp) :: pbp
      integer :: idx,Nnoise

      if (VB_OW) then ; print *,ip_global,"measureKDDW4" ; endif
      Nnoise=10
      pbptot=0.0
      if (VB_OW) then
        if (ip_global.eq.0) then
          open(unit=31,file="cond.dat",status='unknown')
        end if
      end if

      do idx=1,Nnoise
        call evalSingleCondNoise_KDDW4(bmass,imass,pbp)
        pbptot=pbptot+pbp
        if (VB_OW) then
          if (ip_global.eq.0) then
            write(31,*) pbp
          end if
          print *,"pbp comp: ",pbp
        endif
      enddo
      pbptot=pbptot/Nnoise

      if (ip_global.eq.0) then
        print *,"pbp instance: ",pbptot
        if (VB_OW) then
          write(31,*) "average:",pbptot
          close(31)
        end if
      endif

      return
      end subroutine measureKDDW4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine measure_wilsonKDDW(pbp, res, aviter, bmass, imass, isweep)
      use trial
      implicit none
      real,intent(out) :: pbp,aviter
      real,intent(in) :: res,bmass
      integer,intent(in) :: imass
      integer,intent(in),optional :: isweep
      complex(dp) :: cpbp
      integer :: idx
      real(dp) :: vpbp(knoise)
      real :: susclsing
      integer inoise

      if (ip_global.eq.0) then ; print *,"measure_wilsonKDDW" ; end if

      do idx=1,knoise
       if (ip_global.eq.0) then ; print *,ip_global,"measure_wilsonKDDW" ; end if
        call evalSingleCondNoise_KDDW4(bmass,imass,cpbp)
        vpbp(idx)=real(cpbp,dp)
        if (ip_global.eq.0) then ; print *,"pbp noise:",cpbp ; end if
      enddo
      pbp=sum(vpbp)/knoise
      if (ip_global.eq.0) then ; print *,"pbp instance:",pbp ; end if

      call output_OW(pbp,vpbp,isweep)
      aviter = 0 ! not currently calculated - float(iter)/(4*knoise)
      return
      end subroutine measure_wilsonKDDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine output_OW(pbp,vpbp,isweep)
      use trial
      implicit none
      real,intent(in) :: pbp
      real(dp),intent(in) :: vpbp(knoise)
      integer,intent(in),optional :: isweep
      real :: susclsing
      integer inoise
     
      susclsing = 0.0
      do inoise = 1, knoise
        susclsing = susclsing + real(sum(vpbp(inoise)*vpbp(inoise + 1:knoise)))
      enddo
      susclsing = 2*kvol*susclsing/(knoise*(knoise - 1))

      if (ip_global .eq. 0) then
        open (unit=200, file='fort.200', action='write', position='append')
        if (present(isweep)) then
          write (200, '(I5,2E15.7E3)') isweep, pbp, susclsing
        else
          write (200, '(2E15.7E3)') pbp, susclsing
        endif
        close (200)
      end if

      if (ip_global.eq.0) then ; print *,"pbp instance: ",pbp ;endif

      return
      end subroutine output_OW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine evalSingleCondNoise_KDDW4(bmass,imass,pbp)
      use comms
      use comms4
      use gaussian
      implicit none
      real,intent(in) :: bmass
      integer,intent(in) :: imass
      complex(dp),intent(out) :: pbp
      complex(dp) :: eta(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: IDR(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
      complex(dp) :: TMP(1:ksizex_l, 1:ksizey_l, 1:ksizet_l)
      integer idx,ierr
      complex(dp) trcomp,denom
      real :: ps(0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, 2)
      integer, dimension(12) :: reqs_ps
      integer, dimension(12) :: reqs_eta
      real :: one=1.0

      if (VB_OW) then ; if (ip_global.eq.0) then ; print *,"evalSingleCondNoise" ; end if ; end if

      pbp=0.0
      
      do idx=1,4
        eta=0
        call gauss0(ps, reqs_ps)
        call MPI_WaitAll(12, reqs_ps, MPI_Statuses_Ignore, ierr)
        eta(:, :, :, idx) = cmplx(ps(:, :, :, 1), ps(:, :, :, 2))
        call IKDDW4(IDR,eta,bmass,imass)
!        IDR=conjg(eta)*IDR
        IDR=conjg(eta)*(IDR-eta)

        if (imass.eq.1) then
          denom=one-bmass
        elseif (imass.eq.3) then ! assumes gamma3=diag(1 1 -1 -1)
          if (idx.eq.1 .or. (idx.eq.2)) then
            denom=cmplx(-bmass,one)
          else
            denom=cmplx(-bmass,-one)
          endif
        endif
!        TMP=IDR(1:ksizex_l,1:ksizey_l,1:ksizet_l,idx)-one
        TMP=IDR(1:ksizex_l,1:ksizey_l,1:ksizet_l,idx)

        trcomp=sum(TMP)/denom
        pbp=pbp+trcomp
      end do
      pbp=pbp/kvol

      call MPI_AllReduce(MPI_In_Place, pbp, 1, MPI_Double_Complex, &
        & MPI_Sum, comm, ierr)

      return
      end subroutine evalSingleCondNoise_KDDW4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testCondensate()
      use dum1
      use trial,only: u
      use gauge
      use gaussian
      implicit none
      real :: bmass,aviter
      integer :: imass,Ninstance
      complex(dp) :: cpbp,cpbptot
      real :: rpbp,rpbptot
      integer :: idx
      integer :: reqs_ps(12),ierr

      bmass=0.01
      imass=3

      Ninstance=20
      if (ip_global.eq.0) then ; print *,"test Condensate" ; end if
      if (ip_global.eq.0) then ; print *,Ninstance,"instances" ; end if
      cpbptot=0.0
      rpbptot=0.0
      if (ip_global.eq.0) then
        open(unit=31,file="cond.dat",status='unknown')
      end if

      do idx=1,Ninstance

        call gauss0(ps, reqs_ps)
        call MPI_WaitAll(12, reqs_ps, MPI_Statuses_Ignore, ierr)
        theta(:,:,:,1)=ps(1:ksizex_l,1:ksizey_l,1:ksizet_l,1)/3.0       
        theta(:,:,:,2)=ps(1:ksizex_l,1:ksizey_l,1:ksizet_l,2)/3.0   
        call gauss0(ps, reqs_ps)
        call MPI_WaitAll(12, reqs_ps, MPI_Statuses_Ignore, ierr)
        theta(:,:,:,3)=ps(1:ksizex_l,1:ksizey_l,1:ksizet_l,3)/3.0   
        call coef(u,theta)       

        call measureKDDW4(bmass,imass,cpbp)
        call measure_wilsonKDDW(rpbp, respbp, aviter, bmass, imass)
        cpbptot=cpbptot+cpbp
        rpbptot=rpbptot+rpbp
        if (ip_global.eq.0) then
          write(31,*) cpbp,rpbp
          print *,"pbp comp: ",cpbp,rpbp
        end if
      enddo
      cpbptot=cpbptot/Ninstance
      rpbptot=rpbptot/Ninstance

      if (ip_global.eq.0) then
        print *,"pbp average: ",cpbptot,rpbptot
        write(31,*) "average:",cpbptot,rpbptot
        close(31)
      endif

      return
      end subroutine testCondensate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine PRINT_ARRAY(A) !
    use comms_common
    use comms
    implicit none
    complex(dp),intent(in) :: A(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer d,jx,jy,jt
    integer ip,px,py,pt,ierr

    if (ip_global.eq.0) then
      open(unit=31,file="data.dat",status='unknown')
      print *,"T1:"
      write(31,*) "Array:"
      close(31)
    end if

    do d=1,4
      if (ip_global.eq.0) then
        print *,"d:",d
          open(unit=31,file="data.dat",status='unknown',access='append')
        write(31,*) "d:",d
         close(31) 
      endif

    do pt=0,np_t-1
        do jt=1,ksizet_l
      do py=0,np_y-1
          do jy=1,ksizey_l
        do px=0,np_x-1
      
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
          if ((pt.eq.ip_t).and.(py.eq.ip_y).and.(px.eq.ip_x)) then

          open(unit=31,file="data.dat",status='unknown',access='append')

            do jx=1,ksizex_l
              print *,ip_global,A(jx,jy,jt,d)
              write(31,*) ip_global,A(jx,jy,jt,d)
            end do
         
         close(31) 

          end if
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
        end do
       end do
      end do
     end do
    end do

    end do

  end subroutine PRINT_ARRAY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine PRINT_ARRAY5(A) !
    use comms_common
    use comms
    implicit none
    complex(dp),intent(in) :: A(0:kthird_l+1,0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, 4)
    integer d,jx,jy,jt,l
    integer ip,px,py,pt,ierr

   

    if (ip_global.eq.0) then
      open(unit=31,file="data5.dat",status='unknown')
      print *,"T1:"
      write(31,*) "Array:"
      close(31)
    end if

    do l=1,kthird

      if (ip_global.eq.0) then
        print *,"l:",l
          open(unit=31,file="data5.dat",status='unknown',access='append')
        write(31,*) "l:",l
         close(31) 
      endif
    do d=1,4

      if (ip_global.eq.0) then
        print *,"d:",d
          open(unit=31,file="data5.dat",status='unknown',access='append')
        write(31,*) "d:",d
         close(31) 
      endif

    do pt=0,np_t-1
        do jt=1,ksizet_l
      do py=0,np_y-1
          do jy=1,ksizey_l
        do px=0,np_x-1
      
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
          if ((pt.eq.ip_t).and.(py.eq.ip_y).and.(px.eq.ip_x)) then
          open(unit=31,file="data5.dat",status='unknown',access='append')

            do jx=1,ksizex_l
              print *,ip_global,A(l,jx,jy,jt,d)
              write(31,*) ip_global,A(l,jx,jy,jt,d)
            end do
         
         close(31) 

          end if
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
        end do
       end do
      end do
     end do
    end do

    end do
    end do
  end subroutine PRINT_ARRAY5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module measure_OWilson

