!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module polycoeffsmod
      use ellipticmodule
      implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      type poly
!       regular polynomial - ordered c0,c1,c2 ...
        real(prc),allocatable,dimension(:) :: coeffs
      end type poly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      type extpoly
!       extended poly including negative powers 
        integer :: ol,ou
        real(prc),allocatable,dimension(:) :: cfs
      end type 

      type(poly) :: pfcs
      end module polycoeffsmod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module ratfuncs
      use ellipticmodule
      use polycoeffsmod
      implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      type simplepf
!       basic partial fraction
        real(prc) front
        real(prc),allocatable,dimension(:) :: num,denom
      end type simplepf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      type factoredpoly
!       factored polynomial (x-z0)(z-z1)...
        real(prc),allocatable,dimension(:) :: zeros
c        real(prc),allocatable,dimension(:) :: denom
      end type factoredpoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      type gratfunc
!       rational function with regular polys
        type(poly) :: num
        type(poly) :: denom
      end type gratfunc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      type factoredratfunc
!       rational function with factored polys
        type(factoredpoly) :: num
        type(factoredpoly) :: denom
      end type factoredratfunc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      type partfrac
!       partial fraction - with regular poly at front
        type(poly) :: front
        type(factoredpoly) :: denom
        real(prc),allocatable,dimension(:) :: pf
      end type partfrac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      type sgnratfunc
!       rational sgn function in factored and partial fraction forms
        real(prc) :: mult
!       S = m.x. mult_i (num(i)+x^2) / mult_j (denom(j)+x^2)
        type(factoredratfunc) :: frf
!       S = m.x.(front(x^2)+\sum_j  coeffs(j) / (denom(j)+x^2))
        type(partfrac) :: pfrf
!       also include inverse partial fraction
        type(partfrac) :: ipf
!       also include an extended polynomial
!       S=x.\sum_j=-a^b c_j (x'x)^j
        type(extpoly) :: epoly
!       and the basic partial fraction for remez x^-0.5
!       S(=x^-0.5) = front+\sum_i (num(i)/(denom(i)+x)
        type(simplepf) :: spf
      end type sgnratfunc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer,parameter :: VB=0 ! verbose
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine scaleSRatFunc(RF,alpha)
      implicit none
      type(sgnratfunc) :: RF
      real(prc) :: alpha

!     rescale factored form
c      associate(m => RF%mult)
c      associate(cnum=>RF%frf%num%zeros,cdenom => rf%frf%denom%zeros)
c      m=alpha*m
c      cnum=cnum/(alpha*alpha)
c      cdenom=cdenom/(alpha*alpha)
c      end associate

      return
      end subroutine scaleSRatFunc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine printSRatFunc(RF)
      implicit none
      type(sgnratfunc) :: RF

      print *,'Signum Function Rational Approximation:'
      print *,'Factored Form: S = m.x. mult_i (x^2-num(i)) / mult_j (x^2
     &-denom(j)+x^2)'
      print *,'m=',RF%mult
      print *,'num=',RF%frf%num%zeros
      print *,'denom=',RF%frf%denom%zeros
      print *,' '
      print *,'Partial Fraction Form: S = m.x.(front(x^2)+\sum_j  coeffs
     &(j) / (denom(j)+x^2))'
      print *,'mult:',RF%mult
      print *,'front:',RF%pfrf%front%coeffs
      print *,'denom:',RF%pfrf%denom%zeros
      print *,'coeffs:',RF%pfrf%pf

      return
      end subroutine printSRatFunc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function evalSGNFactor(x,rf)
!     assumes form f(x)=mx \prod_i (x^2-num(i))/\prod_j (x^2-denom(j))
      implicit none
      real(prc),intent(in) :: x
      type(sgnratfunc),intent(in) :: rf
      real(prc) f 
      integer nn,nd,j

      nn=size(rf%frf%num%zeros)
      nd=size(rf%frf%denom%zeros)
      f=rf%mult*x
      if (nd.ge.nn) then
        do j=1,nn
          f=f*(x*x-rf%frf%num%zeros(j))
          f=f/(x*x-rf%frf%denom%zeros(j))
        end do
        do j=nn+1,nd
          f=f/(x*x-rf%frf%denom%zeros(j))
        end do
      else
        print *,'nn.gt.nd'
        stop
      end if
      evalSGNFactor=f
      return
      end function evalSGNFactor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function evalISGNFactor(x,rf)
!     assumes form f(x)=mx \prod_i (x^2-num(i))/\prod_j (x^2-denom(j))
      implicit none
      real(prc),intent(in) :: x
      type(sgnratfunc),intent(in) :: rf
      real(prc) f 
      integer nn,nd,j

      nn=size(rf%frf%num%zeros)
      nd=size(rf%frf%denom%zeros)
      f=one/(rf%mult*x)
      if (nd.ge.nn) then
        do j=1,nn
          f=f/(x*x-rf%frf%num%zeros(j))
          f=f*(x*x-rf%frf%denom%zeros(j))
        end do
        do j=nn+1,nd
          f=f*(x*x-rf%frf%denom%zeros(j))
        end do
      else
        print *,'nn.gt.nd'
        stop
      end if
      evalISGNFactor=f
      return
      end function evalISGNFactor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function evalSGNPart(x,rf)
      use pfmodule
      implicit none
      real(prc),intent(in) :: x
      type(sgnratfunc),intent(in) :: rf
      real(prc) f 
      integer np,nd,j

      nd=size(rf%pfrf%denom%zeros)
      f=zero
      if (allocated(rf%pfrf%front%coeffs)) then
        f=evalPoly(x*x,rf%pfrf%front%coeffs)
      end if
      do j=1,nd
        f=f+rf%pfrf%pf(j)/(x*x-rf%pfrf%denom%zeros(j))
      end do
      f=f*rf%mult*x
      evalSGNPart=f
      return
      end function evalSGNPart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function evalISGNPart(x,rf)
      use pfmodule
      implicit none
      real(prc),intent(in) :: x
      type(sgnratfunc),intent(in) :: rf
      real(prc) f 
      integer np,nd,j

      nd=size(rf%ipf%denom%zeros)
      f=zero
      if (allocated(rf%ipf%front%coeffs)) then
        f=evalPoly(x*x,rf%ipf%front%coeffs)
      end if
      do j=1,nd
        f=f+rf%ipf%pf(j)/(x*x-rf%ipf%denom%zeros(j))
      end do
      f=f/(rf%mult*x)
      evalISGNPart=f
      return
      end function evalISGNPart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setHTcoeffs(n,rf)
!      use numbers
      use pfmodule
      implicit none
      integer,intent(in) :: n
      type(sgnratfunc) :: rf
      integer m,j
      integer nn,nd,nf

      if (VB.eq.1) then; print *,'set new HT factored coeffs' ; end if
      call deallocSGNType(rf)

      if (mod(n,2) .eq. 0) then ! even
        m=n/2-1
        allocate(rf%frf%num%zeros(m))
        allocate(rf%frf%denom%zeros(m+1))
        do j=1,m
          rf%frf%num%zeros(j)=-tan(j*pi/n)*tan(j*pi/n)
        end do
        do j=1,m+1
         rf%frf%denom%zeros(j)= 
     &           -tan((j*one-one/2)*pi/n)*tan((j*one-one/2)*pi/n)
        end do
        rf%mult=one*n
      else ! odd
        m=(n-1)/2
        allocate(rf%frf%num%zeros(m))
        allocate(rf%frf%denom%zeros(m))
        do j=1,m
          rf%frf%num%zeros(j)=-tan(j*pi/n)*tan(j*pi/n)
        end do
        do j=1,m
        rf%frf%denom%zeros(j)=
     &           -tan((j*one-one/2)*pi/n)*tan((j*one-one/2)*pi/n)
        end do
        rf%mult=one/n
      end if

      nn=size(rf%frf%num%zeros)
      nd=size(rf%frf%denom%zeros)
      if(VB.eq.1) then ; print *,'nn:',nn,'nd:',nd ; end if

!     make partial fractions
      if(VB.eq.1) then ; print *,'make partial fractions' ; end if
      nf=nn-nd+1
      if(VB.eq.1) then ; print *,nn,nd,nf ; end if
      if (nf.ge.1) then
        allocate(rf%pfrf%front%coeffs(nf))
      end if
      allocate(rf%pfrf%pf(nd))
      allocate(rf%pfrf%denom%zeros(nd))
      call makePartFrac(rf%frf%num%zeros,rf%frf%denom%zeros,
     &                    rf%pfrf%front%coeffs,rf%pfrf%pf)
      rf%pfrf%denom%zeros=rf%frf%denom%zeros
      if(VB.eq.1) then ; print *,rf%pfrf%front%coeffs ; end if
      if(VB.eq.1) then ; print *,rf%pfrf%pf ; end if

      if(VB.eq.1) then ; print *,'make inverse partial fractions';end if
      nf=nd-nn+1
      if (nf.ge.1) then
        allocate(rf%ipf%front%coeffs(nf))
      end if
      allocate(rf%ipf%pf(nn))
      allocate(rf%ipf%denom%zeros(nn))
      call makePartFrac(rf%frf%denom%zeros,rf%frf%num%zeros,
     &                    rf%ipf%front%coeffs,rf%ipf%pf)
      rf%ipf%denom%zeros=rf%frf%num%zeros

      if(VB.eq.1) then ; print *,'HT coeffs set' ; endif
      return
      end subroutine setHTcoeffs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine deallocSGNType(RF)
      implicit none
      type(sgnratfunc) :: RF

!     factored parts
      if (allocated(rf%frf%num%zeros)) then
        deallocate(rf%frf%num%zeros)
        deallocate(rf%frf%denom%zeros)
      end if
!     partial fraction parts, and inverse
      if (allocated(rf%pfrf%front%coeffs)) then
        deallocate(rf%pfrf%front%coeffs)
      end if
      if (allocated(rf%ipf%front%coeffs)) then
        deallocate(rf%ipf%front%coeffs)
      end if
      if (allocated(rf%pfrf%pf)) then
        deallocate(rf%pfrf%pf)
        deallocate(rf%pfrf%denom%zeros)
        deallocate(rf%ipf%pf)
        deallocate(rf%ipf%denom%zeros)
      end if
      if (allocated(rf%epoly%cfs)) then
        deallocate(rf%epoly%cfs)
      endif

      return
      end subroutine deallocSGNType
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setZoloCoeffs(Ns,rf,xmin,xmax)
!      use numbers
      use zolomodule
      use pfmodule
      implicit none
      integer,intent(in) :: Ns
      type(sgnratfunc) :: rf
      real(prc) xmin,xmax
      integer m,j
      integer nn,nd,nf,n,nnum,ndenom
      type(zolotarev) zolo

      if(VBS) then ; print *,'set new Zolo factored coeffs',Ns ; end if
      call setZolo(xmin,xmax,Ns,zolo)

      if (modulo(Ns,2).eq.0) then ! Ns even
        nnum=Ns/2-1
        ndenom=Ns/2
c        n=Ns/2
      else
        nnum=Ns/2
        ndenom=Ns/2
c        n=Ns/2
      endif

      call deallocSGNType(rf)
      allocate(rf%frf%num%zeros(nnum)) ! n-1
      allocate(rf%frf%denom%zeros(ndenom)) ! n
      associate(num => rf%frf%num%zeros,denom => rf%frf%denom%zeros)
      denom=zolo%cdenom*zolo%xmin*zolo%xmin
      num=zolo%cnum(1:nnum)*zolo%xmin*zolo%xmin
      rf%mult=zolo%mult/zolo%xmin
c      do j=1,n-1
c        rf%mult=rf%mult*denom(j)/num(j)
c      end do
c      rf%mult=-rf%mult*denom(n)
      do j=1,nnum
        rf%mult=rf%mult*denom(j)/num(j)
      end do
      do j=nnum+1,ndenom
        rf%mult=-rf%mult*denom(j)
      end do

      nn=size(rf%frf%num%zeros)
      nd=size(rf%frf%denom%zeros)
      if(VBS)then
        print *,'nn:',nn
        print *,rf%frf%num%zeros
        print *,'nd:',nd
        print *,rf%frf%denom%zeros
      endif

!     make partial fractions
      nf=nn-nd+1
      if(VBS)then
        print *,'make partial fractions'
        print *,nn,nd,nf
      endif
      if (nf.ge.1) then
        allocate(rf%pfrf%front%coeffs(nf))
      end if
      allocate(rf%pfrf%pf(nd))
      allocate(rf%pfrf%denom%zeros(nd))
      call makePartFrac(rf%frf%num%zeros,rf%frf%denom%zeros,
     &                    rf%pfrf%front%coeffs,rf%pfrf%pf)
      rf%pfrf%denom%zeros=rf%frf%denom%zeros
      if(VBS)then
        print *,rf%pfrf%front%coeffs
        print *,rf%pfrf%pf
        print *,'make inverse partial fractions'
      endif
      nf=nd-nn+1
      if (nf.ge.1) then
        allocate(rf%ipf%front%coeffs(nf))
      end if
      allocate(rf%ipf%pf(nn))
      allocate(rf%ipf%denom%zeros(nn))
      call makePartFrac(rf%frf%denom%zeros,rf%frf%num%zeros,
     &                    rf%ipf%front%coeffs,rf%ipf%pf)
      rf%ipf%denom%zeros=rf%frf%num%zeros

      end associate

      return
      end subroutine setZoloCoeffs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      subroutine setPolyCoeffs(polyt,func,order,r1,r2)
!      use pfmodule
!      use polyapproxmod
!      use polycoeffsmod
!      implicit none
!      type(poly),intent(out) :: polyt
!      integer,intent(in) :: func,order
!      real(prc),intent(in) :: r1,r2 
! 
!      if (allocated(polyt%coeffs)) then
!        deallocate(polyt%coeffs)
!      endif
!      call getPolyCoeffs(func,order,r1,r2,polyt%coeffs)
!
!      return
!      end subroutine setPolyCoeffs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      subroutine setPolyRFCoeffs(SRF,func,order,r1,r2)
!      use pfmodule
!      use polyapproxmod
!      use polycoeffsmod
!      implicit none
!      type(sgnratfunc),intent(out) :: SRF
!      integer,intent(in) :: func,order
!      real(prc),intent(in) :: r1,r2 
!
!      call deallocSGNType(SRF)
!      call setPolyCoeffs(SRF%pfrf%front,func,order,r1,r2)
!
!      return
!      end subroutine setPolyRFCoeffs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testHTFunctions(Nht,xL,xR)
!      use pacc
!      use numbers
      implicit none
      integer Nht
      real(prc) xL,xR
      real(prc) mx,x,Sf,Sp
      type(sgnratfunc) :: HT,HTs
      integer Nz,j,scale
      real(prc) alpha

      print *,"test Hyperbolic Tanh Functions"
      call setHTcoeffs(Nht,HT)

      call printSRatFunc(HT)

      Nz=100
      mx=(xR**(one/real(Nz/2,prc))-one)
      print *,mx

      x=xL
      open(unit=10,file='HT.dat',form='formatted',status='unknown')
      do j=1,Nz+1
        Sf=evalSGNFactor(x,HT)
        Sp=evalSGNPart(x,HT)
        print *,x,Sf,Sp,one-Sf,one-Sp
        write(10,*) x,Sf,Sp,one-Sf,one-Sp
        x=x*(one+mx)
      end do
      close(10)

      scale=2
      alpha=real(scale,prc);
      call setHTcoeffs(Nht/scale,HTs)
c      alpha=1d0/alpha
      x=xL
      open(unit=10,file='HTs2.dat',form='formatted',status='unknown')
      do j=1,Nz+1
        Sf=evalSGNFactor(alpha*x,HTs)
        Sp=evalSGNPart(alpha*x,HTs)
        print *,x,Sf,Sp,one-Sf,one-Sp
        write(10,*) x,Sf,Sp,one-Sf,one-Sp
        x=x*(one+mx)
      end do
      close(10)

      scale=4
      alpha=real(scale,prc);
      call setHTcoeffs(Nht/scale,HTs)
c      alpha=1d0/alpha
      x=xL
      open(unit=10,file='HTs4.dat',form='formatted',status='unknown')
      do j=1,Nz+1
        Sf=evalSGNFactor(alpha*x,HTs)
        Sp=evalSGNPart(alpha*x,HTs)
        print *,x,Sf,Sp,one-Sf,one-Sp
        write(10,*) x,Sf,Sp,one-Sf,one-Sp
        x=x*(one+mx)
      end do
      close(10)

      scale=8
      alpha=real(scale,prc);
      call setHTcoeffs(Nht/scale,HTs)
c      alpha=1d0/alpha
      x=xL
      open(unit=10,file='HTs8.dat',form='formatted',status='unknown')
      do j=1,Nz+1
        Sf=evalSGNFactor(alpha*x,HTs)
        Sp=evalSGNPart(alpha*x,HTs)
        print *,x,Sf,Sp,one-Sf,one-Sp
        write(10,*) x,Sf,Sp,one-Sf,one-Sp
        x=x*(one+mx)
      end do
      close(10)

      return
      end subroutine testHTFunctions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testZoloFunctions(Nz,lmin,lmax)
!      use pacc
!      use numbers
      use zolomodule
      implicit none
      integer Nz
      real(prc) lmin,lmax
      real(prc) dx,x,Sf,Sp,xdash
      type(sgnratfunc) :: Zolo
      integer Nx,j

      print *,"test Zolotarev Functions",Nz,lmin,lmax
      call setZoloCoeffs(Nz,Zolo,lmin,lmax)
      call printSRatFunc(Zolo)
      Nx=100
c      dx=(lmax-lmin)/(Nx-1)
      x=lmin
      dx=(log(lmax/lmin)-log(one))/(Nx-1)
      open(unit=10,file='Zolo.dat',form='formatted',status='unknown')
      xdash=0
      do j=1,Nx
        Sf=evalSGNFactor(x,Zolo)
        Sp=evalSGNPart(x,Zolo)
        write(10,*) x,Sf,Sp,one-Sf,one-Sp
        print *,x,Sf,Sp,one-Sf,one-Sp
c        x=x+dx
        xdash=xdash+dx
        x=lmin*exp(xdash)
      end do
      close(10)

      return
      end subroutine testZoloFunctions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine plotRationalFunction(RF,xmin,xmax,fname)
!      use pacc
!      use numbers
      implicit none
      type(sgnratfunc) :: RF
      real(prc) xmin,xmax
      character fname*6
      integer Nz,j
      real(prc) dx,x,Sf,ISf,Sp,ISp

      print *,"plot Rational Function"
      Nz=30
      dx=(xmax-xmin)/(Nz-1)
      x=one
      if(fname.ne.'none') then
        open(unit=10,file=fname,status='unknown')
      endif
      do j=1,Nz
        Sf=evalSGNFactor(x,RF)
        Sp=evalSGNPart(x,RF)
        print *,x,Sf,Sp,one-Sf,one-Sp
        if(fname.ne.'none') then
          write(10,*) x,Sf,Sp,one-Sf,one-Sp
        endif
        x=x+dx
      end do
      if(fname.ne.'none') then
        close(10)
      endif

      return
      end subroutine plotRationalFunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module ratfuncs
