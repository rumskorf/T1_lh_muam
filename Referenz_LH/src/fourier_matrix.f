      subroutine fourier_matrix
      include 'com_main.fc'
      include 'com_namelist.fc' !FL
c      parameter(igit=64,nb=36,nfmax=igit/2+1,nfmaxl=10,nfleg=10)
      parameter(nfmax=igit/2+1,nfmaxl=10,nfleg=10)
      double precision fun(nfmax,igit),snorm
      double precision theta(nb),pleg(nfmaxl,nb)
      double precision pf1(igit,igit),pf3(igit,igit),pf5 (igit,igit),
     $                 pf7(igit,igit),pf9(igit,igit),pf11(igit,igit),
     $                 pf2(igit,igit),pf12(igit,igit), pfleg(nb,nb)
c      common/legandr/pfleg
c      common/fourier/pf1,pf3,pf5,pf7,pf9,pf11,pf2,pf12
      double precision a(nfmax,nfmax),am1(nfmax,nfmax),u(nfmax,nfmax),
     $     v(nfmax,nfmax),sigma(nfmax),wk1(nfmax),wk2(nfmax)
      double precision al(nfmaxl,nfmaxl),am1l(nfmaxl,nfmaxl),
     $	   ul(nfmaxl,nfmaxl),
     $     vl(nfmaxl,nfmaxl),sigmal(nfmaxl),wk1l(nfmaxl),wk2l(nfmaxl)
      
      !open(1,file='U_test_1.dat')
c      open(2,file='U_test_2.dat')
      
      pi=2.*asin(1.)
      dlon=2.*pi/float(igit) !abstand gitterpunkte breite in rad
      do i=1,igit
        xlon=float(i-1)*dlon !breitengrad in rad
        fun(1,i)=1.0D0 
        do if=1,igit/4 !1,16 (0.90deg) -> zuweisung von 33 werten zu jeder breite
        iif=2*if
        fun(iif  ,i)=dsin(dble(if)*xlon)
        fun(iif+1,i)=dcos(dble(if)*xlon)
        end do
      end do
      DO IF=1,IGIT/2+1 !1,33
       SNORM=0.0D0
       DO I=1,IGIT
          SNORM=SNORM+FUN(IF,I)*FUN(IF,I) !snorm=aufaddieren von fun**2 ueber alle breiten
       END DO
        DO I=1,IGIT
          FUN(IF,I)=FUN(IF,I)/DSQRT(SNORM) !Normierung
        END DO
      END DO
c	harmonic filtration: 1 harmonic (on poles) to 12 harmonics ( on equator)      
      relerr=1.e-4
      print*,'1 harmonic'
      call filter(igit,igit,igit,nfmax, 3,fun,fun,a,am1,
     $            u,v,sigma,wk1,wk2, pf1,relerr)
      print*,'2 harmonics'
      call filter(igit,igit,igit,nfmax, 5,fun,fun,a,am1,
     $            u,v,sigma,wk1,wk2, pf2,relerr)
      print*,'3 harmonics'
      call filter(igit,igit,igit,nfmax, 7,fun,fun,a,am1,
     $            u,v,sigma,wk1,wk2, pf3,relerr)
      print*,'5 harmonics'
      call filter(igit,igit,igit,nfmax,11,fun,fun,a,am1,
     $            u,v,sigma,wk1,wk2, pf5,relerr)
      print*,'7 harmonics'
      call filter(igit,igit,igit,nfmax,15,fun,fun,a,am1,
     $            u,v,sigma,wk1,wk2, pf7,relerr)
      print*,'9 harmonics'
      call filter(igit,igit,igit,nfmax,19,fun,fun,a,am1,
     $            u,v,sigma,wk1,wk2, pf9,relerr)
      print*,'11 harmonics'
      call filter(igit,igit,igit,nfmax,23,fun,fun,a,am1,
     $            u,v,sigma,wk1,wk2,pf11,relerr)
      print*,'12 harmonics'
      call filter(igit,igit,igit,nfmax,25,fun,fun,a,am1,
     $            u,v,sigma,wk1,wk2,pf12,relerr)
      do j=1,nb 
c------------------------------------------- co-latitude COMMA
        theta(j)=pi/180.*(2.5+ 5.*FLOAT(j-1)) !lat in rad 0-pi
      end do
      CALL POLLEG(nb,nfmaxl,nfleg,1,theta,pleg) !P0L Legendre: berechnet pleg(1-10,j)
c      CALL PLLEG0(NJ,NFMAX,NF,ARG,PL0) 
      print*,'filter 10 legendre'
      !filter 10 -> 4.5 harmonics?
      CALL FILTER(nb,nb,nb,nfmaxl,nfleg,pleg,pleg,
     $            Al,AM1l,Ul,Vl,SIGMAl,WK1l,WK2l,pfleg,RELERR)
      print*,pfleg(10,30),pf1(1,1),pf2(1,1),pf3(1,1),pf5(1,1),pf7(1,1),
     +	     pf9(1,1),pf11(1,1),pf12(1,1)	 
      return
      end 
c	------------------------------------------------
c      call filter(njmax=igit,njs=igit,njd=igit,nfmax=33,nf=25,pls=fun,pld=fun,a,am1,u,v,sigma,wk1,wk2,pf12,relerr=1.e-4)
      SUBROUTINE FILTER(NJMAX,NJS,NJD,NFMAX,NF,
     $ PLS,PLD,A,AM1,U,V,SIGMA,WK1,WK2,PF,RELERR)
      double precision PLD(NFMAX,NJD),PLS(NFMAX,NJS),
     $     AM1(NFMAX,NF),
     $     V(NFMAX,NF),SIGMA(NF),WK1(NF),WK2(NF),
     $     PF(NJMAX,NJD)
      double precision A(NFMAX,NF),U(NFMAX,NF)
      integer NFMAX
      DO 10 IF=1,NF
	DO 20 KF=1,NF
	  S=0.0D0
	  DO 30 JD=1,NJD
   30 	    S=S+PLD(IF,JD)*PLD(KF,JD)
	  A(KF,IF)=S
   20   CONTINUE
   10 CONTINUE
      CALL INVSVD(NFMAX,NF,NF,A,AM1,U,V,SIGMA,WK1,RELERR) !call INVSVD(nfmax=33,nf=25,nf=25,A,AM1,U,V,SIGMA,WK1,RELERR)
      DO 40 JD=1,NJD
	DO 50 IF=1,NF
	  S=0.0D0
	    DO 55 KF=1,NF
   55 	      S=S+AM1(KF,IF)*PLD(KF,JD)
	    WK1(IF)=S
   50 	CONTINUE
	DO 60 JS=1,NJS
	  DO 70 KF=1,NF
c---- 		filter out the wave with zonal wave number 6, 1, 2
c      		if(kf.eq.12.or.kf.eq.13) then
c      		if(kf.eq.2.or.kf.eq.3) then
c      		if(kf.eq.4.or.kf.eq.5) then
c        	wk2(kf)=0.
c                             else
	    WK2(KF)=PLS(KF,JS)
c      		end if
   70 	  continue
	  CALL PRVV(WK1,NF,WK2,RESULT)
	  PF(JS,JD)=RESULT
   60 	CONTINUE
   40 CONTINUE
      RETURN
      END
c------------------------------------------------------------
c	call INVSVD(ndim=33,m=25,n=25,A,AM1,U,V,SIGMA,WK1,RELERR) 
      SUBROUTINE INVSVD(NDIM,M,N,A,AM1,U,V,SIGMA,WORK,RELERR)
      DOUBLE PRECISION S
      double precision V(NDIM,N),
     *     SIGMA(N),WORK(N),AM1(NDIM,M)
      double precision A(NDIM,N),U(NDIM,N)
      integer NDIM
      CALL SVD(NDIM,M,N,A,SIGMA,.TRUE.,U,.TRUE.,V,IERR,WORK) !single value decomposition: call svd(ndim=33,m=25,n=25,A,W=Sigma,.TRUE.,U,.TRUE.,V,IERR,WORK)
      IF(IERR.NE.0) WRITE(99,1) IERR !trouble shooting
   1  FORMAT(/'   /INVSVD/    TROUBLE.IERR = ',I4/)
      SIGMA1=0.0D0
      DO 2 I=1,N
	IF(SIGMA(I).GT.SIGMA1) SIGMA1=SIGMA(I)
   2  CONTINUE
      TAU=RELERR*SIGMA1
      DO 4 I=1,N
	DO 4 K=1,M
	  S=0.0D0
	  DO 3 J=1,N
	    IF(SIGMA(J).LE.TAU) GO TO 3
	    S=S+DBLE(V(I,J)*U(K,J)/SIGMA(J))
   3  	  CONTINUE
	  AM1(I,K)=S
   4  CONTINUE
      RETURN
      END
c------------------------------------------------------
c	call svd(nm=33,m=25,n=25,A,W=Sigma,MATU=.TRUE.,U,MATV=.TRUE.,V,IERR,RV1)
      SUBROUTINE SVD(NM,M,N,A,W,MATU,U,MATV,V,IERR,RV1)
c
c     this subroutine is a translation of the algol procedure svd,
c     numer. math. 14, 403-420(1970) by g.h. golub and c. reinsch.
c     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).

c     this subroutine determines the singular value decomposition: 
c     computation of the singular values and complete orthogonal 
c     decomposition of a real rectangular matrix A:
c     A = U diag(RV1) V**T  of a real m*n rectangular matrix.  
c     Householder' Bidiagonalization and a variant of the qr algorithm are used.

c the actual parameters corresponding to a,u,v may all be identical unless matu=matv=true
c in this case, the actual parameters corresponding to u and v must differ. m>=n assumed
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.  note that nm must be at least
c          as large as the maximum of m and n.
c
c        m is the number of rows of a (and u).
c
c        n is the number of columns of a (and u) and the order of v.
c
c        a contains the rectangular input matrix to be decomposed.
c
c        matu should be set to .true. if the u matrix in the
c          decomposition is desired, and to .false. otherwise.
c
c        matv should be set to .true. if the v matrix in the
c          decomposition is desired, and to .false. otherwise.
c
c     on output
c
c        a is unaltered (unless overwritten by u or v).
c
c        w contains the n (non-negative) singular values of a (the
c          diagonal elements of s).  they are unordered.  if an
c          error exit is made, the singular values should be correct
c          for indices ierr+1,ierr+2,...,n.
c
c        u contains the matrix u (orthogonal column vectors) of the
c          decomposition if matu has been set to .true.  otherwise
c          u is used as a temporary array.  u may coincide with a.
c          if an error exit is made, the columns of u corresponding
c          to indices of correct singular values should be correct.
c
c        v contains the matrix v (orthogonal) of the decomposition if
c          matv has been set to .true.  otherwise v is not referenced.
c          v may also coincide with a if u is not needed.  if an error
c          exit is made, the columns of v corresponding to indices of
c          correct singular values should be correct.
c
c        ierr is set to
c          zero       for normal return,
c          k          if the k-th singular value has not been
c                     determined after 30 iterations.
c
c        rv1 is a temporary storage array.

C	   A(33,25),W(25),U(33,25),V(33,25),RV1(25)	

      integer i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its,ierr
      double precision A(NM,N),W(N),U(NM,N),V(NM,N),RV1(N)
      double precision C,F,G,H,S,Y,Z,tst1,tst2,SCALE,ANORM      
     
      LOGICAL MATU,MATV
C
      IERR=0
      DO 100 I=1,M !1,25
	DO 100 J=1,N !1,25
  100 	  U(I,J)=A(I,J) !in I nur bis 25 statt bis 33
c     householder reduction to bidiagonal form  
      G=0.0D0
      SCALE=0.0D0
      ANORM=0.0D0
      DO 300 I=1,N !1,25
	L=I+1
	RV1(I)=SCALE*G
	G=0.0D0
	S=0.0D0
	SCALE=0.0D0
	IF(I.GT.M) then !wenn I groesser 25 wird
c	  print*,'l152 I.GT.M',i,m
	  GO TO 210
	endif  
	DO 120 K=I,M !von laufindex I bis 25
	  SCALE=SCALE+DABS(U(K,I)) !aufaddition der U von I bis 25 -> je groesser I desto kleiner wird scale
  120 	CONTINUE
	
	IF(SCALE.EQ.0.0D0) then
c	  print*,'l160 scale.eq.0',scale
	  GO TO 210
	endif
	
	DO 130 K=I,M !von laufindex I bis 25
	  U(K,I)=U(K,I)/SCALE !normierung?
c	  if(n.eq.25) print*,k,i,u(k,i)
c	  if(i.eq.24.and.k.eq.25) print*,'4',k,i,u(i,k),scale
	  S=S+U(K,I)**2
  130 	CONTINUE
	F=U(I,I)	
	G=-DSIGN(DSQRT(S),F)
	H=F*G-S
	U(I,I)=F-G	
	IF(I.EQ.N) then
c	  print*,'l173 i.eq.n',i,n
	  GO TO 190
	endif  
	DO 150 J=L,N
	  S=0.0D0
	  DO 140 K=I,M
	    S=S+U(K,I)*U(K,J)
  140 	  CONTINUE
	  F=S/H
	  DO 150 K=I,M !I-25
c	    if(i.eq.24.and.k.eq.25) print*,'5',j,k,i,u(i,k),F,scale
c	    if(k.eq.24.and.j.eq.25) print*,'6',j,k,i,u(k,j),u(k,i),F
c	    if(n.eq.10) print*,n,i,j,k,u(k,j),F,u(k,i),F*U(K,I)
	    if(u(k,j).eq.0.) then
	      pr=1
	      print*,'v',n,i,j,k,u(k,j),F,u(k,i),F*U(K,I) !U(24,25) falsch!
	    else 
	      pr=0
	    endif  
	    U(K,J)=U(K,J)+F*U(K,I)
	    if(pr.eq.1) print*,'n',n,i,j,k,u(k,j)
c	    if(i.eq.24.and.k.eq.25) print*,'0',k,i,u(i,k),F,scale
cFL 	    U(K,J) wird hier 0 -> vorlaeufig ersetzt		    	    	    
c	    if (abs(u(k,j)).lt.1.e-10) then
c	      write(1,*) k,j,u(k,j)
c	    endif	     
c	    if(u(k,j).eq.0) u(k,j)=0.1
  150 	CONTINUE
  190	continue

   	DO 200 K=I,M !I,25
   	  U(K,I)=SCALE*U(K,I)
c   	  if(k.eq.24.and.i.eq.25) print*,'1',k,i,u(k,i),scale
  200	continue
  210 	W(I)=SCALE*G
	G=0.0D0
	S=0.0D0
	scale=0.0D0 ! forgotten?
	IF(I.GT.M.OR.I.EQ.N) GO TO 290
	DO 220 K=L,N !(I+1),25
  220 	  SCALE=SCALE+DABS(U(I,K))
	IF(SCALE.EQ.0.0D0) GO TO 290
	DO 230 K=L,N!(I+1),25
c	  if(k.eq.24.and.i.eq.25) print*,'2',i,k,u(k,i),scale
c	  if(n.eq.25) print*,i,k,u(i,k)
	  U(I,K)=U(I,K)/SCALE
	  S=S+U(I,K)**2
  230 	CONTINUE
	F=U(I,L)
	G=-DSIGN(DSQRT(S),F)
	H=F*G-S
	U(I,L)=F-G
	DO 240 K=L,N
c	  print*,i,k,H 
c	  if (H.eq.0) then
c	    write(1,*) i,k,u(i,k),H,F,G,S
c	    H=-2.3391174E-13
c	  endif
c	  if(n.eq.25) print*,i,k,u(i,k),H
c	  hilf1=real(H,8)
c	  if(n.eq.25.and.k.eq.25) print*,i,k,u(i,k),hilf1
c	  hilf=U(I,K)/hilf1	  
c	  RV1=sngl(hilf)
  240 	  RV1(K)=U(I,K)/H
c  240	  continue	
c240	  if(n.eq.25) print*,i,k,RV1(K)
	IF(I.EQ.M) GO TO 270
	DO 260 J=L,M
	  S=0.0D0
	  DO 250 K=L,N
  250 	    S=S+U(J,K)*U(I,K)
	  DO 260 K=L,N
  260 	    U(J,K)=U(J,K)+S*RV1(K)
  
  270 	DO 280 K=L,N
  280 	  U(I,K)=SCALE*U(I,K)
  
  290 	ANORM=DMAX1(ANORM,DABS(W(I))+DABS(RV1(I)))
  
  300 CONTINUE
  
c	accumulation of right-hand transformations 
      IF(.NOT.MATV) GO TO 410
c	for i=n step -1 until 1 do --      
      DO 400 II=1,N
	I=N+1-II
	IF(I.EQ.N ) GO TO 390
	IF(G.EQ.0.0D0) GO TO 360
	DO 320 J=L,N
  320 	  V(J,I)=(U(I,J)/U(I,L))/G
	DO 350 J=L,N
	  S=0.0D0
	  DO 340 K=L,N
  340 	    S=S+U(I,K)*V(K,J)
	  DO 350 K=L,N
  350 	    V(K,J)=V(K,J)+S*V(K,I)
  
  360 	DO 380 J=L,N
	  V(I,J)=0.0D0
	  V(J,I)=0.0D0
  380 	CONTINUE
  
  390 	V(I,I)=1.0D0
	G=RV1(I)
	L=I
      
  400 CONTINUE
  
c  accumulation of left-hand transformations  
  410 IF(.NOT.MATU) GO TO 510
c	for i=min(m,n) step -1 until 1 do  
      MN=N
      IF(M.LT.N)MN=M
      DO 500 II=1,MN
	I=MN+1-II
	L=I+1
	G=W(I)
	IF(I.EQ.N)  GO TO 430
	DO 420 J=L,N
  420	  U(I,J)=0.0D0
  
  430 	IF(G.EQ.0.0D0) GO TO 475
	IF(I.EQ.MN) GO TO 460
	DO 450 J=L,N
	  S=0.0D0
	  DO 440 K=L,M
  440 	    S=S+U(K,I)*U(K,J)
	  F=(S/U(I,I))/G
	  DO 450 K=I,M
	    U(K,J)=U(K,J)+F*U(K,I)
  450 	CONTINUE
  
  460 	DO 470 J=I,M
  470 	  U(J,I)=U(J,I)/G
	GO TO 490
  475 	DO 480 J=I,M
  480 	  U(J,I)=0.0D0
  
  490 	U(I,I)=U(I,I)+1.
  
  500 CONTINUE
  
c	diagonalization of the bidiagonal form  
      !510 tst1=anorm
  510 DO 700 KK=1,N
	K1=N-KK
	K=K1+1
	ITS=0 !no. of iterations
c	test for splitting
c	for l=k step -1 until do --
  520 	DO 530 LL=1,K
	  L1=K-LL
	  L=L1+1
	  !tst2=tst1+dabs(rv1(l))
	  !if (tst2.eq.tst1) go to 565
c	  rv1(1) is always zero, so there is no exit through the bottom of the loop	  
	  IF(DABS(RV1(L))+ANORM.EQ.ANORM)  GO TO 565
	  !tst2=tst1+dabs(w(l1))
	  !if (tst2.eq.tst1) go to 540
	  IF(DABS(W(L1)) +ANORM.EQ.ANORM)  GO TO 540
  530 	CONTINUE
  
c	cancellation of rv1(l) if L greater than 1  
  540 	C=0.0D0
	S=1.0D0
	DO 560 I=L,K
	  F=S*RV1(I)
	  RV1(I)=C*RV1(I)
	  !tst2=tst1+dabs(f)
	  !if (tst2 .eq. tst1) go to 565
	  IF(DABS(F)+ANORM.EQ.ANORM) GO TO 565
	  G=W(I)
	  H=DSQRT(F*F+G*G)
	  W(I)=H
	  C= G/H
	  S=-F/H
	  IF(.NOT.MATU) GO TO 560
	  DO 550 J=1,M
	    Y=U(J,L1)
	    Z=U(J,I)
	    U(J,L1)=Y*C+Z*S
	    U(J,I)=-Y*S+Z*C
  550 	  CONTINUE
  
  560   CONTINUE

c	test for converegence  
  565 	Z=W(K)
	IF( L .EQ. K) GO TO  650
c	shift from bottom 2 by 2 minor	
	!IF(ITS.EQ.30) GO TO 1000
	IF(ITS.EQ.300) GO TO 1000 !if no. of iterations exceeds this value
	ITS=ITS+1 !next iteration step
	X=W(L)
	Y=W(K1)
	G=RV1(K1)
	H=RV1(K)
	F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.*H*Y)
	G=DSQRT(F*F+1.0D0)
	F=((X-Z)*(X+Z)+H*(Y/(F+DSIGN(G,F))-H))/X
c	next qr transformation	
C
	C=1.0D0
	S=1.0D0
	DO 600 I1=L,K1
	  I=I1+1
	  G=RV1(I)
	  Y=W(I)
	  H=S*G
	  G=C*G
	  Z=DSQRT(F*F+H*H)
	  RV1(I1)=Z
	  C=F/Z
	  S=H/Z
	  F= X*C+G*S
	  G=-X*S+G*C
	  H=Y*S
	  Y=Y*C
	  IF(.NOT.MATV) GO TO 575
	  DO 570 J=1,N
	    X=V(J,I1)
	    Z=V(J,I)
	    V(J,I1)=X*C+Z*S
	    V(J,I)=-X*S+Z*C
  570 	  CONTINUE
  
  575 	  Z=DSQRT(F*F+H*H)
	  W(I1)=Z
c	  rotation can be arbitrary ig z is zero	  
	  IF(Z.EQ.0.0D0)   GO TO 580
	  C=F/Z
	  S=H/Z
	  
  580 	  F=C*G+S*Y
	  X=-S*G+C*Y
	  IF(.NOT.MATU) GO TO 600
	  DO 590 J=1,M
	    Y=U(J,I1)
	    Z=U(J,I)
	    U(J,I1)=Y*C+Z*S
	    U(J,I)=-Y*S+Z*C
  590 	  CONTINUE
  
  600 	CONTINUE
  
	RV1(L)=0.0D0
	RV1(K)=F
	W(K)  =X
	GO TO 520
c	convergence	
  650 	IF(Z.GE.0.0D0)   GO TO 700
c	w(k) is made non-negative  
	W(K)=-Z
	IF(.NOT.MATV) GO TO 700
	DO 690 J=1,N
  690	  V(J,K)=-V(J,K)
  
  700 CONTINUE
      GO TO 1001
c	set error -- no converegence to a sungular value after 30 iterations      
 1000 IERR=K
C
 1001 RETURN
      END
C	-------------------<POLLEG>---------------------
c	POLLEG(nb,nfmaxl=10,nfleg=10,1,theta,pleg)      
      SUBROUTINE POLLEG(NJ,NFMAX,NF,M,ARG,PL)
      double precision ARG(NJ),PL(NFMAX,NJ) !ARG=breite in rad, PL(10,nb) ist ergebnis
      double precision X,XX,XSQ
      NFM1=NF-1
      DO 10 J=1,NJ !ueber alle breiten
	X=ARG(J) 
	X=DCOS(X)
	XX=X*X
	IF(M.GT.1) GO TO 20 !nicht erfuellt, M=1
	  
	  XSQ=DSQRT(1.-XX)
	  PL(1,J)=XSQ !Anfangsbedingung 1
	  PL(2,J)=3.0D0*X*XSQ !AB 2
	  GO TO 30
   
   20 	  XSQ=1.-XX 		!uebersprungen
	  PL(1,J)=3.0D0*XSQ	!ueberspringen
	  PL(2,J)=15.0D0*X*XSQ	!ueberspringen	
	  
   30 	  CONTINUE
	  DO 40 IF=2,NFM1 !2,9 -> PL(3-10,j)
	    KF=IF-M+1 !KF=IF (M=1)
   40 	    PL(IF+1,J)=((2.0D0*IF+1.)*X*PL(IF,J)-(IF+M)*PL(IF-1,J))/KF
   10 CONTINUE
      RETURN 
      END
c---------------------------------------                  
      SUBROUTINE  PRMV(NDIM,NI,NJ,A,B,C)
      double precision  B(NJ), C(NI)
      double precision A(NDIM,NJ) 
      DO  2  I = 1, NI
          SUM  = 0.0D0
      DO  1  J = 1, NJ
          SUM  = SUM + A(I,J) * B(J)
   1  CONTINUE
          C(I) = SUM
   2  CONTINUE
      RETURN
      END
c--------------------------------------
      SUBROUTINE  PRVV(A,NI,B,RESULT)
      double precision  A(NI), B(NI)
        RESULT = 0.0D0
      DO  1  I = 1, NI
        RESULT = RESULT + A(I) * B(I)
   1  CONTINUE
      RETURN
      END

