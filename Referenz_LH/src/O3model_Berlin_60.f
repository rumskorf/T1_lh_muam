      subroutine O3model_Berlin_60(nmon)
c
c     rearranged for Berlin data 24.01.2006, Alex Pogoreltsev
c
c----------------- nb=36 is fixed for given O3 model
c -----------------and kgit has to be set (arbitrary)
c---------------------------------------------------------------------
      include 'com_main.fc'
csl
      include 'com_namelist.fc' !FL
csl
      parameter (klat=37,klev=50,kmon=12)
c      parameter (nb=36,kgit=48)
      real       o3vmr (klat,klev,kmon),blat(klat),blev(klev)
      real       press(klev),temp(klev)
      real       o3lat(klat),B(klat),C(klat),D(klat)
      real       comlat,o3out1,o3out2(nb,klev)
      
      real	 faktor,tave,tavekm1,tempjk
      character*10 text(klat,klev,kmon)
      
      !by FL 2014
      integer jahr_1000,jahr_100, jahr_10, jahr_1
      character*4 jahr
      
      common/o3_Berl/blev,o3out2,h_km 
c --------------------------------------------------------------------
c------ CO3(45,nb) will be used to calculate O3 cooling
c
c	O3ppmv - O3 part per million volume
c	CO3    - for O3 cooling
c	SO3    - Integral of O3 number density
c	VMOINT - integralnaya koncentraciya O3
c	VMO    - koncentraciya O3
c---------------------------------------------------------------------
c      DIMENSION O3PPMV(nb,kgit), CO3(45,nb), SO3(nb,kgit)
      EXTERNAL VMOINT
      REAL Z_km(kgit+1)
      data press/899.,795.,701.,616.,540.,472.,411.,357.,308.,265.,
     $ 		 226.,193.,165.,141.,120.,103., 88., 75., 64., 55.,
     $		  47., 40., 35., 30., 25., 22., 19., 16., 14., 12.,
     $            10., 8.9, 7.7, 6.6, 5.7, 5.0, 4.3, 3.8, 3.3, 2.9,
     $		  2.5, 2.2, 1.9, 1.7, 1.5, 1.3, 1.2, 1.0,0.90,0.78/ !taken from US standard atmosphere 1-50km
      data temp/
     $ 	    281.7,275.2,268.7,262.2,255.7,249.2,242.7,236.2,229.7,223.3,
     $	    216.8,9*216.7,
     $	    217.6,218.6,219.6,220.6,221.6,222.5,223.5,224.5,225.5,226.5,
     $	    227.5,228.5,231.0,233.7,236.5,239.3,242.1,244.8,247.6,250.4,
     $	    253.1,255.9,258.6,261.4,264.2,266.4,269.7,3*270.7/
      h_km=7.
      
      !by FL 2014
      jahr_1000=YEAR/1000
      jahr_100 =(YEAR-jahr_1000*1000)/100
      jahr_10  =(YEAR-jahr_1000*1000-jahr_100*100)/10
      jahr_1   =(YEAR-jahr_1000*1000-jahr_100*100-jahr_10*10)
      jahr=char(jahr_1000+48)//char(jahr_100+48)//char(jahr_10+48)//
     &	     char(jahr_1+48)
      !reading data
      open  (10,file=FO3)
      read  (10,*)
      read  (10,*)
      read  (10,*)
      do jt=1,kmon
         read(10,*)
         read(10,*)
         read(10,*)
         write (6,*) 'SPARC / month : ',jt
         do jk=1,klev
	    read  (10,'(I3,37(F10.2))') n,(o3vmr(jl,jk,jt),jl=1,klat) !in DU/km
	    do jl=1,klat
	      if (o3vmr(jl,jk,jt).GT.999.) o3vmr(jl,jk,jt)=0.	   
	      o3vmr(jl,jk,jt)=o3vmr(jl,jk,jt)*real( 2.7d20*1.38d-23
     &		     	      *temp(jk)/(press(jk)*1.d2)  !!V(DU)=NkT/p
     &		     	      /1.d3 *1.d6 ,4)! /V(1m*1m*1km)*10^6 <- *10^6 for units in ppm
	    enddo  
c	    print*, o3vmr(1,jk,jt)
	 enddo             
      end do
      close (10)         
      
      do k=1,klev
        blev(k)=-h_km*log(press(k)/press(1)) !log-p height for all 50 levels
        do j=1,klat
csl---------------------------------------------------------------------------------------        
csl       o3lat(j)=o3vmr(j,k,nmon)                                      ! original
csl	  o3lat(j)=o3vmr(j,k,nmon)*1.02                                 ! + 2%
Csl	  o3lat(j)=o3vmr(j,k,nmon)-0.5*o3vmr(j,k,nmon)                  ! Halbierung der Ozonkonzentration
Csl	  o3lat(j)=o3vmr(j,k,nmon)*OZON                                 ! mit Faktor aus Namelist
Csl       factor = sin(float(j-1)*PI/float(klat-1))**2                  ! wichtungsfaktor sin quadrat
c          
	    pi=2.*asin(1.) !FL

c	Wichtungsfunktion: 
c	OZON=1 keine Wichtung (OZON-1=0) (Ozon bleibt Original)
c	OZON=0 maximale Wichtung (OZON-1=-1)
c	sin(float(j-1)*PI/float(klat-1))**2: =1 am ??quator, =0 an den Polen
c	d.h. bei OZON=0 verschwindet Ozon am ??quator und bleibt Original an den Polen
	    o3lat(j)= o3vmr(j,k,nmon)!+   
c	    print*,o3lat(j), o3vmr(j,k,nmon)
!      &             ((OZON-1.0) * sin(float(j-1)*PI/float(klat-1))**2)   
!      &            * o3vmr(j,k,nmon)     
	    
csl---------------------------------------------------------------------------------------
	    blat(j)=-90.+float(j-1)*180./float(klat-1) !latitudes in degree -90...90
          end do
c          print*,o3lat(10)
          CALL SPL3(klat,blat,o3lat,B,C,D) !anstiege 1.,2.,3. Grades der splines?
          do jc=1,nb
	    comlat = -87.5+float(jc-1)*180./float(nb) !comma latitudes
	    CALL SEVAL0(klat,comlat,blat,o3lat,B,C,D,o3out1) !umrechnen von ozonwerten von ursprungsgitter in commagitter
	    o3out2(jc,k)=o3out1
          end do
      end do
      
c-----------------------------------------------------------------
      DZKM=135./(48.-0.5)
      ZUNTKM=0.5*DZKM
      Z_km(1)=0.
      DO 20 I=1,NB
         DO 10 K=1,KGIT
                Z_km(K+1)=ZUNTKM+FLOAT(K-1)*DZKM
                OZPPMV(I,K)=VMO(Z_km(K+1),I)
c                write(*,'(A10,2I3,F10.4)') 'ergebnis:',i,k,ozppmv(j,k)                
   10    CONTINUE      
c------------------------------ O3 ppmv for cooling
         DO 11 K=1,45
                ZINT=FLOAT(K-1)*0.25*h_km
                CO3(K,I)=VMO(ZINT,I)
 11      CONTINUE
 20   CONTINUE 
      
C
      ZOBEN1=Z_km(11)
      ZOBEN2=200.
      EPSMAX=1.E-6
      DO 42 I=1,NB
      DO 41 K=1,KGIT
       IF(K.LT.10) THEN
        CALL ADAPT(I,VMOINT,Z_km(K+1),ZOBEN1,EPS,EPSMAX,RES1)
        CALL ADAPT(I,VMOINT,ZOBEN1   ,ZOBEN2,EPS,EPSMAX,RES2)
         RES=RES1+RES2
       ELSE
         CALL ADAPT(I,VMOINT,Z_km(K+1),ZOBEN2,EPS,EPSMAX,RES)
       END IF 
c         CALL ADAPT(PHI(I),I,VMOINT,Z_km(K+1),ZOBEN,EPS,EPSMAX,SW)
c--------------   N_O3*1e5(km->cm)*P_0/K_B/T_0
         SO3(I,K)=RES*1.e5*1.013/1.38066e-16/273.
 41   CONTINUE
 42   CONTINUE       
      
      RETURN
      END 
   
C_______________________________________________________________________   
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C--------------------------SUBROUTINES----------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C_______________________________________________________________________

C
      SUBROUTINE ADAPT (II,FCT,A1,B1,EPS,EPSMAX,XINT)
C      EXTERNAL FCT
      DIMENSION QI(25),BI(25)
      NMAX = 25
      XFF  = 8.192E3
      XFAK = 1./(XFF - 1.)
C
      BI(1) = B1
      XA    = A1
      CALL GAUSS (II,FCT,QI(1),XA,B1)
C
      XINT  = 0.
      EPS   = 0.
      J     = 1
C
 100  IF (J .LT. 1)                RETURN
      QI1  = QI(J)
      BIE  = BI(J)
      XAP  = 0.5*(XA + BIE)
      CALL GAUSS (II,FCT,QI2,XA,XAP)
C
      CALL GAUSS (II,FCT,QI3,XAP,BIE)
C
      QI0  = QI2 + QI3
      DELT = XFAK*ABS(QI0 - QI1)
      DREL = DELT
      IF (QI0 .NE. 0.)             DREL = DELT/ABS(QI0)
      IF (DREL .LE. EPSMAX)        GO  TO  10
C
      JP1 = J + 1
      IF (JP1 .GT. NMAX)           GO  TO  10
      QI(J)  = QI3
C
      J = JP1
      BI(J) = XAP
      QI(J) = QI2
      GO  TO  100
C
 10   XINT = XFAK*(XFF*QI0 - QI1) + XINT
      EPS  = EPS + DELT
      XA   = BIE
      J    = J - 1
      GO  TO  100
C
      END
C
C-----------------------------------------------------------
C
      SUBROUTINE GAUSS (II,FCT,VALUE,AF,BF)
C
      DIMENSION STUETZ(6) , GEW(6)
      DATA STUETZ /-0.9324695,-0.6612094,-0.2386192,0.2386192
     *,0.6612094,0.9324695/
      DATA GEW /0.1713245,0.3607616,0.4679139,0.4679139,
     * 0.3607616,0.1713245/
C
      HH1   = 0.5*(BF - AF)
      VALUE = 0.
      XXM   = 0.5*(AF + BF)
      DO  100  I = 1 , 6
          VAR = XXM + HH1*STUETZ(I)
          SUM = FCT (VAR,II)*GEW(I)
 100      VALUE = VALUE + SUM
C
      VALUE = HH1*VALUE
      RETURN
      END
C
C-----------------------------------------------------------
      FUNCTION VMOINT(ZZ,I)
      Parameter(klev=50,nb=36)
      Real blev(klev),o3lev(klev),o3out2(nb,klev),
     $     B(klev),C(klev),D(klev) 
      common/o3_Berl/blev,o3out2,h_km 
        do k=1,klev
          o3lev(k)=o3out2(I,k)
        end do
      CALL SPL3(klev,blev,o3lev,B,C,D)
      if(zz.le.blev(klev)) then
        CALL SEVAL0(klev,zz,blev,o3lev,B,C,D,O3ppmv)
      else
        O3ppmv=o3lev(klev)*exp((blev(klev)-zz)/h_km)
      end if
      VMOINT=O3ppmv*EXP(-ZZ/h_km)
      RETURN
      END
C--------------------------------------------------------------
      FUNCTION VMO(ZZ,I)
      Parameter(klev=50,nb=36)
      real blev(klev),o3lev(klev),o3out2(nb,klev),
     $     B(klev),C(klev),D(klev) 
      common/o3_Berl/blev,o3out2,h_km 
        do k=1,klev
          o3lev(k)=o3out2(I,k)
        end do
      CALL SPL3(klev,blev,o3lev,B,C,D)
      if(zz.le.blev(klev)) then
        CALL SEVAL0(klev,zz,blev,o3lev,B,C,D,O3ppmv)
      else
        O3ppmv=o3lev(klev)*exp((blev(klev)-zz)/h_km)
      end if
      VMO=O3ppmv
      RETURN
      END
c--------------------------------------------------------------
      SUBROUTINE  SPL3(N,X,Y,B,C,D)
      REAL  T, X(N), Y(N), B(N), C(N), D(N)
         NM1 = N - 1
      IF( N .LT. 2 )  RETURN
      IF( N .LT. 3 )  GO TO 50
        D(1)   =  X(2)-X(1)
        C(2)   = (Y(2)-Y(1))/D(1)
      DO 10  I = 2, NM1
        D(I)   =  X(I+1)-X(I)
        B(I)   = 2.*(D(I-1)+D(I))
        C(I+1) = (Y(I+1)-Y(I))/D(I)
        C(I)   =  C(I+1)-C(I)
  10  CONTINUE
       B(1)  = -D(1)
       B(N)  = -D(N-1)
       C(1)  =  0.
       C(N)  =  0.
      IF( N .EQ. 3 )  GO TO 15
       C(1)  = C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
       C(N)  = C(N-1)/(X(N)-X(N-2))
     *                         -C(N-2)/(X(N-1)-X(N-3))
       C(1)  = C(1)*D(1)**2/(X(4)-X(1))
       C(N)  =-C(N)*D(N-1)**2/(X(N)-X(N-3))
  15  CONTINUE
      DO 20  I = 2, N
          T    = D(I-1)/B(I-1)
          B(I) = B(I)-T*D(I-1)
          C(I) = C(I)-T*C(I-1)
  20  CONTINUE
           C(N) = C(N)/B(N)
      DO 30  IB = 1, NM1
           I    = N - IB
           C(I) = (C(I)-D(I)*C(I+1))/B(I)
  30  CONTINUE
       B(N)    = (Y(N)-Y(NM1))/D(NM1)
     *          + D(NM1)*(C(NM1)+2.*C(N))
      DO 40  I = 1, NM1
       B(I)    = (Y(I+1)-Y(I))/D(I)
     *                        -D(I)*(C(I+1)+2.*C(I))
       D(I)    = (C(I+1)-C(I))/D(I)
       C(I)    = 3.*C(I)
  40  CONTINUE
       C(N)    = 3.*C(N)
       D(N)    = D(N-1)
      RETURN
  50  CONTINUE
       B(1)  =  (Y(2)-Y(1))/(X(2)-X(1))
       C(1)  =  0.
       D(1)  =  0.
       B(2)  =  B(1)
       C(2)  =  0.
       D(2)  =  0.
      RETURN
      END
c----------------------------------------------------------
c	SEVAL0(klat,comlat,blat,o3lat,B,C,D,o3out1)
c	N - No. grid points in old grid
c	U - current latitude (deg) of new grid (comma)
c	X - latitude field of old grid (deg) with N grid points
c	Y - data to be interpolated (ozone) with N grid points
c	B - aus SPL3
c	C - aus SPL3
c	D - aus SPL3
c	S - field of data (ozone) with new grid (comma latitudes)
      SUBROUTINE  SEVAL0(N,U,X,Y,B,C,D,S)
      REAL  U, X(N), Y(N), B(N), C(N), D(N), S
      DATA  I/1/
      IF( I .GE. N      )    I = 1 !if N is only 1 data point
      IF( U .LT. X(I)   )  GO TO 10 !extrapolieren
      IF( U .LE. X(I+1) )  GO TO 30 !interpolieren
  10  CONTINUE
      I = 1
      J = N + 1
  20  CONTINUE
      K = (I+J) / 2
      IF( U .LT. X(K) )  J = K
      IF( U .GE. X(K) )  I = K
      IF( J .GT. I+1  )  GO TO 20
  30  CONTINUE
      DX  =  U - X(I)
      S   =  Y(I) + DX*(B(I)+DX*(C(I)+DX*D(I)))
      RETURN
      END

