      module pfmodule
      use ellipticmodule
!      use pacc
!      use numbers
!     partial fraction from factor of from (x-n0)(x-n1).../(x-d1)(x-d2)...
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function evalfpoly(x,zeros)
!     factor coeffs of poly 
      implicit none
      real(prc) x
      real(prc),dimension(:) :: zeros
      integer nc,j

      nc=size(zeros)
      evalfpoly=one
      do j=1,nc
        evalfpoly=evalfpoly*(x-zeros(j))
      end do
      return
      end function evalfpoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function evalpoly(x,coeffs)
!     standard poly coeffs 
      implicit none
      real(prc) x
      real(prc),dimension(:) :: coeffs
      integer nc,j

      nc=size(coeffs)
      if (nc.eq.1) then
        evalpoly=coeffs(1)
        return
      endif
      evalpoly=x*coeffs(nc)+coeffs(nc-1)
      do j=nc-2,1,-1
        evalpoly=x*evalpoly+coeffs(j)
      end do
      return
      end function evalpoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine expandPoly(zeros,epoly)
!     expand factored poly
      implicit none
      real(prc),dimension(:) :: zeros,epoly
      integer j,nn,nd,i

      nn=size(zeros)
      epoly=0
      epoly(1)=-zeros(1)
      epoly(2)=one
      do j=2,nn
        epoly(j+1)=1
        do i=j,2,-1
          epoly(i)=epoly(i-1)-epoly(i)*zeros(j)
        end do
        epoly(1)=-epoly(1)*zeros(j)
      end do

      return
      end subroutine expandPoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dividePoly(enum,edenom,front,dnum)
      implicit none
      real(prc),dimension(:) :: enum,edenom
      real(prc),dimension(:) :: front
      real(prc),dimension(:) :: dnum
      real(prc),allocatable,dimension(:) :: tnum
      integer j,nn,nd,nf,k
      nn=size(enum)
      nd=size(edenom)
      if (nn.lt.nd) then
        dnum=enum
        front=0
        return
      else
        allocate(tnum(nn))
        tnum=enum
        nf=1+nn-nd
        do j=1,nf
          front(nf+1-j)=tnum(nn+1-j)/edenom(nd)
c          print *,'front',front
          do k=1,nd
        tnum(nn+1-j+1-k)=tnum(nn+1-j+1-k)-front(nf+1-j)*edenom(nd+1-k)
          end do
c          print *,'tnum',tnum
        end do
        dnum(1:nd-1)=tnum(1:nd-1)
        deallocate(tnum)
      end if
      return

      end subroutine dividePoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makePartFrac(num,denom,front,pfcoeffs)
!     assumes denom has distinct real coefficients
!     passes in factored form (coefficients are zeros)
      implicit none
      real(prc),dimension(:) :: num,denom
      real(prc),dimension(:) :: front,pfcoeffs
      real(prc),allocatable,dimension(:) :: enum,edenom,pnum
      real(prc) z
      integer j,nn,nd,nf

      nn=size(num)
      nd=size(denom)
      allocate(enum(nn+1))
      call expandPoly(num,enum)
      allocate(edenom(nd+1))
      call expandPoly(denom,edenom)
      nf=nn-nd+1
      if (nf.ge.1) then
        allocate(pnum(nd))
        call dividePoly(enum,edenom,front,pnum)
        do j=1,nd
          z=denom(j)
          pfcoeffs(j)=evalpoly(z,pnum)/polyDeriv(z,denom)
        end do
        deallocate(pnum)
      else
        do j=1,nd
          z=denom(j)
          pfcoeffs(j)=evalfpoly(z,num)/polyDeriv(z,denom)
        end do
      end if
      deallocate(enum)
      deallocate(edenom)
      
      return
      end subroutine makePartFrac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function polyDeriv(x,zeros)
!     factor coeffs of poly 
      implicit none
      real(prc) x
      real(prc),dimension(:) :: zeros
      integer nc
      integer j,i
      real(prc) st

      nc=size(zeros)
      polyderiv=0
      do j=1,nc
        st=1
        do i=1,nc
          if (i.ne.j) then
            st=st*(x-zeros(i))
          end if
        end do
        polyderiv=polyderiv+st
      end do
      return
      end function polyDeriv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testPartFrac()
      implicit none

      real(prc) num1(2),num2(4)
      real(prc) denom(3)
      real(prc) pfront1(1),pfront2(2)
      real(prc) pf(3)

      num1(1)=4d0
      num1(2)=5d0

      denom(1)=1d0
      denom(2)=2d0
      denom(3)=3d0

      print *,num1
      print *,denom
      call makePartFrac(num1,denom,pfront1,pf)
      print *,pf
      print *,'should be',6,-6,1
      print *,pfront1
      print *,'should be unassigned'

      num2(1:2)=num1(1:2)
      num2(3)=6d0
      num2(4)=7d0

      print *,num2
      print *,denom
      call makePartFrac(num2,denom,pfront2,pf)
      print *,pf
      print *,'should be',180,-120,12
      print *,pfront2
      print *,'should be',-16,1

      return
      end subroutine testPartFrac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testDividePoly()
      implicit none
      real(prc) :: num(5) = (/0d0,1d0,1d0,1d0,1d0/)
      real(prc) :: denom(3) = (/1d0,0d0,1d0/)
      real(prc) :: front(3)
      real(prc) :: dnum(2)

      print *,'num',num
      print *,'denom',denom
      call dividePoly(num,denom,front,dnum)
      print *,'front',front
      print *,'should be:',0,1,1
      print *,'dnum',dnum
      print *,'should be:',0,0

      return
      end subroutine testDividePoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module pfmodule
