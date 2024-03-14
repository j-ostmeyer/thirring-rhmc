!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module zolomodule
      use ellipticmodule
      implicit none
      logical,parameter :: VBS=.false.

      type zolotarev
        integer N
        real(zprc) xmin,xmax
        real(zprc) k
        real(zprc) kdash
        real(zprc) bigK
        real(zprc) bigKdash
        real(zprc),allocatable :: cnum(:)
        real(zprc),allocatable :: cdenom(:)
        real(zprc) M
        real(zprc) C
        real(zprc) lambda
        real(zprc) mult
        real(zprc),allocatable :: roots(:)
        
      end type zolotarev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setZolo(lmin,lmax,N,zolo)
!      use ellipticmodule
      implicit none

      real(prc) lmin,lmax
      integer N
      type(zolotarev) :: zolo
      real(zprc) beta
      real(zprc) secondkind

      beta=lmax/lmin
      zolo%xmin=lmin
      zolo%xmax=lmax
      zolo%N = N
      zolo%k = one/beta
      zolo%kdash = sqrt(1-zolo%k)*sqrt(1+zolo%k)
      call COMELP(zolo%k,zolo%bigK,secondkind) 
      call COMELP(zolo%kdash,zolo%bigKdash,secondkind) 
      if(VBS)then
        print *,"N:",zolo%N
        print *,"k:",zolo%k
        print *,"kdash:",zolo%kdash
        print *,"bigK:",zolo%bigK
        print *,"bigKdash:",zolo%bigKdash
      endif
      allocate(zolo%cnum(N/2))
      allocate(zolo%cdenom(N/2))
      call setInnerCoeffs(zolo)
      zolo%M = calcM(zolo)
      zolo%C = calcC(zolo)
      zolo%lambda = calcLambda(zolo)
      if(VBS)then
        print *,"cnum:",zolo%cnum
        print *,"cdenom:",zolo%cdenom
        print *,"M:",zolo%M
        print *,"C:",zolo%C
        print *,"lambda:",zolo%lambda
      endif
      zolo%mult=2*zolo%lambda/(one+zolo%lambda)/zolo%M

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setInnerCoeffs(zolo)
!      use ellipticmodule
      implicit none
      type(zolotarev) zolo
      integer l
      real(zprc) argC,JargC
      real(zprc) argCdash,JargCdash
      real(zprc) ecn,edn,eph

      do l=0,zolo%N/2-1
        argC = (2*l+2)*zolo%bigKdash/(zolo%N);
        call jelp(argC,zolo%kdash,JargC, ecn, edn, eph )
        zolo%cnum(l+1) = -JargC/(one-JargC)*JargC/(one+JargC)
        if(VBS)then ; print *,"cnum(l):",zolo%cnum(l+1) ; endif
      end do
      do l=0,zolo%N/2-1
        argCdash = (2*l+1)*zolo%bigKdash/(zolo%N);
        call jelp(argCdash,zolo%kdash,JargCdash, ecn, edn, eph )
c        zolo%cdenom(l) = -JargCdash*JargCdash/(1-JargCdash*JargCdash);
      zolo%cdenom(l+1)=-JargCdash/(one-JargCdash)*
     &                                    JargCdash/(1+JargCdash)
c        print *,"argCdash: ",argCdash,"JargCdash: ",JargCdash
        if(VBS)then ; print *,"cdenom(l):",zolo%cdenom(l+1) ; endif
      end do
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function Rsgn(x,zolo) 
!      use numbers
      implicit none
      real(prc) x
      type(zolotarev) zolo
      integer l
      Rsgn=one
      do l=1,zolo%N-1
        Rsgn=Rsgn*(one-x/zolo%cnum(l))/(one-x/zolo%cdenom(l))
      end do
      Rsgn=Rsgn/(one-x/zolo%cdenom(zolo%N));
      return
      end function Rsgn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function R(x,zolo) 
!      use numbers
      implicit none
      real(prc) x
      type(zolotarev) zolo
      integer l
      R=one
      do l=1,zolo%N/2
        R=R*(one-x/zolo%cnum(l))/(one-x/zolo%cdenom(l))
      end do
      return
      end function R
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function calcM(zolo)
!      use numbers
      implicit none
      type(zolotarev) zolo

c      calcM=Rsgn(one,zolo)
      calcM=R(one,zolo)
      return
      end function calcM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function calcC(zolo) 
!      use ellipticmodule
      implicit none
      type(zolotarev) zolo
      real(zprc) bKD
      real(zprc) num
      real(zprc) denom
      real(zprc) esn,ecn,edn,eph
      real(zprc) dsn,dcn,ddn,dph

cc      bKD = zolo%bigKdash/(2*zolo%N)
      bKD = zolo%bigKdash/(zolo%N)
c      num = J.sn(bigK)*J.dnd(bKD);
      call jelp(zolo%bigK,zolo%k, esn, ecn, edn, eph )
      call jelp(bKD,zolo%kdash, dsn, dcn, ddn, dph )
      num = esn*ddn;
c      denom = J.dn2(bigK)*J.snd2(bKD);
      denom = edn*edn*dsn*dsn;
      calcC = num/(1.0-denom)
      return
      end function calcC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function calcLambda(zolo) 
      implicit none
      type(zolotarev) zolo
      real(prc) C2

      C2=zolo%C*zolo%C
c      calcLambda = zolo%M/(zolo%C*Rsgn(C2,zolo));
      calcLambda = zolo%M/(zolo%C*R(C2,zolo));
      return
      end function calcLambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      real(prc) function zsgn(x,zolo) 
c      use numbers
c      implicit none
c      real(prc) x
c      type(zolotarev) zolo
c      
c      x=x/zolo%xmin
c      zsgn=zolo%mult*x*Rsgn(x*x,zolo)
c      return
c      end function zsgn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function zsgn(x,zolo) 
!      use numbers
      implicit none
      real(prc) x
      type(zolotarev) zolo
      real(prc) arg
      
      arg=x/zolo%xmin*x/zolo%xmin
      zsgn=x/zolo%xmin*zolo%mult*R(arg,zolo)
      return
      end function zsgn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getRoots(zolo) 
!      use ellipticmodule
      implicit none
! this uses Chiu's formula for the roots, but uses the lambda calculated
! with Kennedys formulation (ie no theta function)
      type(zolotarev) zolo
      real(zprc) l,M,Malt,lkappa
      real(zprc) argT,argB,arg
      real(zprc) asn,c1,c2,omega
      real(zprc) snda,sndb,sn2,t1,t2,t3
      real(zprc),allocatable :: vs(:)
      integer s,so2,Ns
        
c      print *,"getRoots"
      l=zolo%lambda;
      M=zolo%M  
      Ns=zolo%N
c      print *,l,M,Ns
  
      allocate(zolo%roots(Ns))
c      print *,"lambda",l 
c      print *,"M",M 
  
      Malt=1.0;
      do s=1,Ns/2
c        print *,"s",s
        argT = (2*s-1)*zolo%bigKdash/Ns;
        argB = (2*s)*zolo%bigKdash/Ns;
        call jelp(argT,zolo%kdash, snda, t1, t2, t3 )
        call jelp(argB,zolo%kdash, sndb, t1, t2, t3 )
        Malt = Malt*snda/sndb*snda/sndb;
c       Malt = Malt*zoloCoeffs.J.snd2(argT)/zoloCoeffs.J.snd2(argB);
      end do
c      print *,"altM",Malt
c      print *,"ialtM",1.0/Malt 
      allocate(vs(Ns))
      arg=sqrt( (1.0+3.0*l)/( (1.0+l)**3 ) );
      lkappa=sqrt((1.0+l)*(1.0-l));
c      print *,"arg",arg,"lkappa",lkappa
      asn = arcsn(arg,lkappa,10000);
c      asn = zoloCoeffs.J.asn(arg,lkappa);
c      print *,"asn:",asn
      c1 = M*asn;
      c2 = 2*zolo%bigKdash/Ns;  
c      print *,"c1:",c1
c      print *,"c2:",c2
      do s=0,Ns-1
        so2 = (s+1)/2
c        print *,"so2:",so2
        vs(s+1) = (-one)**(s)*c1 + so2*c2
c        print *,"plus/minus:", (-one)**(s)*c1
c        print *,"base:",so2*c2
c        print *,s,"vs:",vs(s+1)
      enddo

      do s=0,Ns-1
        call jelp(vs(s+1),zolo%kdash, snda, t1, t2, t3 )
        sn2=snda*snda
c       double sn2 = zoloCoeffs.J.snd2(vs[s]);
c       double kdash = zoloCoeffs.kdash;
        arg = 1 - zolo%kdash*zolo%kdash*sn2;
        omega = sqrt(arg)/zolo%xmin
        zolo%roots(s+1) = 1/omega;
c        print *,s,"roots(s)",1/omega
      end do
      deallocate(vs)
      return
      end subroutine getRoots
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function zsgnrts(x,zolo) 
!     calculate the sign function from the roots
!      use numbers
      implicit none
      real(prc) x
      type(zolotarev) zolo
      real(prc) T,num,denom
      integer Nrts,j

      Nrts=size(zolo%roots)
      T=one
      do j=1,Nrts
        num=zolo%roots(j)-x
        denom=zolo%roots(j)+x
        T=T*num/denom
      end do
      zsgnrts=(1-T)/(1+T)
      return
      end function zsgnrts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testZolo(Nr)
!      use numbers
      implicit none
      integer Nr
      type(zolotarev) zolo
      integer Nz,j
      real(prc) dx,x,xdash,Sest,SestRT,lmin,lmax
     
      lmin=one/1000 ; lmax=20*one
      call setZolo(lmin,lmax,Nr,zolo)
      call getRoots(zolo)
      open(unit=10,file='zolo.dat',form='formatted',status='unknown')
      Nz=200
      dx=(lmax-lmin)/(Nz-1)
      dx=(log(lmax)-log(lmin))/(Nz-1)
      x=lmin
      xdash=0
      do j=1,Nz
        Sest=zsgn(x,zolo)
        SestRT=zsgnrts(x,zolo)
        print *,x,Sest,Sest-one,SestRT,SestRT-one
        write(10,'(5E13.5E2)') x,Sest,Sest-one,SestRT,SestRT-one
        x=x+dx
        xdash=xdash+dx
        x=lmin*exp(xdash)
      end do
      close(10)
      open(unit=10,file='zrts.dat',form='formatted',status='unknown')
      do j=1,Nr
        write(10,*) zolo%roots(j),1,0
      end do
      close(10)
      return
      end subroutine testZolo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module zolomodule

