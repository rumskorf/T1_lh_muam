C****************************************************************************
C*                                                                          *
C*                       SUBROUTINE CCO2GR
C*                                                                          *
C****************************************************************************
C
C  to calculate parameterization coefficients for 15 um CO2 band cooling
C  rate calculation for arbitrary CO2 volume mixing ratio (without height
C  variation up to x=12.5) in a range of 150-720ppm.
C  Coefficients are calculated for atmosphere layer from the pressure scale
C  height x=2 up to x=16.5. Irregular vertical grid is used to account for
C  internal heat exchange within "LTE layer" (x=0-12.5)
C  ***** Important!!!!!*****
C  From x=2 up to x=13.75 coefficients can be calculated only for a vertical
C  grid with a step of 0.25. Eventually, any vertical step could be used 
C  within the atmospheric layer of x=14-16.5. However, to do it,  some
C  modification is needed. No coefficients are needed to calculate cooling
C  rates above x=16.5 (see driving programm).
C
C                                    May, 1996.  V.I. Fomichev.

C Called by a driving programm
C Calls A18LIN (linear interpolation)
C Calls A18INT (third order spline)

      SUBROUTINE CCO2GR(CO2,AM,UTOP,RCO2)

C input:
C  CO2(17) - CO2 vmr at a grid of x=12.5-16.5 with a step = 0.25
C  AM(17)  - air molecular weight (g/mole) at the same grid
C  UTOP    -  CO2 column amount (cm-2) above the level of x=16.5
C  RCO2 - CO2 volume mixing in the region below x=12.5
C  initial data from BLOCK DATA PCO2O3 come through common blocks

      REAL CO2(17), AM(17)
      real uref(4), co2int(4), uco2(17)
C output: parameterization coefficients for both, the matrix parameterization 
C         (is used between x=2 and 12.5) and reccurence formula. 
C         Passing on to other subroutines through common block CO2CFG.
C  AMAT,BMAT(43,9) - coefficients for the matrix parameterization
C  AL(17) - coeficients for the reccurence formula. Note that starting up
C           from x=14.00, these coefficients are identical to escape functions
C           which can be calculated at any arbitrary vertical grid.
C           So that starting up from this level any arbitrary vertical grid
C           can be used.

      COMMON /CO2CFG/ AMAT(43,9),BMAT(43,9),AL(17)

C data from BLOCK DATA PCO2O3
      include 'pco2o3.f'

c      common /PIRGRD/ xr(67), xrh(59), xrf(17)
c      common /PIRCO2/ co2vp(67), co2cl(67)
c      common /PIRLTE/ a150(43,9),b150(43,9), a360(43,9),b360(43,9),
c     1                a540(43,9),b540(43,9), a720(43,9),b720(43,9),
c     2                co2o(4)
c      common /PIRMGR/ IG(9)
c      common /PIRNE/ uco2ro(51),alo(51),
c     1               cor150(6),cor360(6),cor540(6),cor720(6),uco2co(6)


c auxiliary arrays

c      real uref(4), co2int(4), uco2(17)

c calculate coefficients for the matrix paramerization:

      do 1 i = 1,43
      do 1 j = 1,9
      if((i.le.5).and.(j.eq.2)) goto 1
      isgn = int(sign(1.,a150(i,j))+sign(1.,a360(i,j))+
     +           sign(1.,a540(i,j))+sign(1.,a720(i,j)))
      co2int(1)=a150(i,j)/co2o(1)
      co2int(2)=a360(i,j)/co2o(2)
      co2int(3)=a540(i,j)/co2o(3)
      co2int(4)=a720(i,j)/co2o(4)
      if(isgn.eq.-4) then
      co2int(1) = alog(-co2int(1))
      co2int(2) = alog(-co2int(2))
      co2int(3) = alog(-co2int(3))
      co2int(4) = alog(-co2int(4))
       a = -exp(a18lin(RCO2,co2o,co2int,1,4))
      else if (isgn.eq.4) then
      co2int(1) = alog(co2int(1))
      co2int(2) = alog(co2int(2))
      co2int(3) = alog(co2int(3))
      co2int(4) = alog(co2int(4))
       a = exp(a18lin(RCO2,co2o,co2int,1,4))
      else
      call a18int(co2o,co2int,RCO2,a,4,1)
      end if
       AMAT(i,j)=a*RCO2

      isgn = int(sign(1.,b150(i,j))+sign(1.,b360(i,j))+
     +           sign(1.,b540(i,j))+sign(1.,b720(i,j)))
      co2int(1)=b150(i,j)/co2o(1)
      co2int(2)=b360(i,j)/co2o(2)
      co2int(3)=b540(i,j)/co2o(3)
      co2int(4)=b720(i,j)/co2o(4)
      if(isgn.eq.-4) then
      co2int(1) = alog(-co2int(1))
      co2int(2) = alog(-co2int(2))
      co2int(3) = alog(-co2int(3))
      co2int(4) = alog(-co2int(4))
       a = -exp(a18lin(RCO2,co2o,co2int,1,4))
      else if (isgn.eq.4) then
      co2int(1) = alog(co2int(1))
      co2int(2) = alog(co2int(2))
      co2int(3) = alog(co2int(3))
      co2int(4) = alog(co2int(4))
       a = exp(a18lin(RCO2,co2o,co2int,1,4))
      else
      call a18int(co2o,co2int,RCO2,a,4,1)
      end if
       BMAT(i,j)=a*RCO2
    1 continue

c calculate uco2(17) (CO2 column amount) for CO2 and AM profiles at XRF grid:

      uco2(17)=utop

      const = 1.e6*6.024e23/970.*0.25/2.
      u = co2(17)*exp(-xrf(17))/am(17)
      do 2 ii=1,16
      i=17-ii
      u1=co2(i)*exp(-xrf(i))/am(i)
      uco2(i) = uco2(i+1) + const*(u1+u)
      u=u1
    2 continue

c calculate coeficients for the reccurence formula:
c between x=12.5 and 13.75 these coefficients (al) are calculated using
c correction to escape function. Starting up from x=14.00 parameterization
c coeficients equal escape function.

      do 3 i=1,6
      a = a18lin(uco2(i),uco2ro,alo,1,51)
      co2int(1)=cor150(i)
      co2int(2)=cor360(i)
      co2int(3)=cor540(i)
      co2int(4)=cor720(i)
      uref(1) =uco2co(i)*150./360.
      uref(2) =uco2co(i)
      uref(3) =uco2co(i)*540./360.
      uref(4) =uco2co(i)*720./360.
      cor = a18lin(uco2(i),uref,co2int,1,4)
      AL(i)=exp(cor+a)
    3 continue

      do 4 i=7,17
      a = a18lin(uco2(i),uco2ro,alo,1,51)
      AL(i)=exp(a)
    4 continue

      return
      end

C****************************************************************************
C*                                                                          *
C*                            FUNCTION A18LIN                               *
C*                                                                          *
C****************************************************************************
C
C linear interpolation
C
C called by cco2gr
C calls nothing
C
C input:
C  X - argument for which a value of function should be found
C  XN(N),YN(N) - values of function YN(N) at XN(N) grid. X(N) should be
C                ordered so that X(I-1) < X(I).
C output:
C  A18LIN - value of function for X

      FUNCTION A18LIN(X,XN,YN,M,N)

      dimension XN(N),YN(N)

      k=m-1
      do 1 i=m,n
      k=k+1
      if(x-xn(i)) 2,2,1
    1 continue
    2 if(k.eq.1) k=2

c k has been found so that xn(k).le.x.lt.xn(k+1)

      A18LIN=(yn(k)-yn(k-1))/(xn(k)-xn(k-1))*(x-xn(k))+yn(k)
      return
      end

C****************************************************************************
C*                                                                          *
C*                    SUBROUTINE A18INT                                     *
C*                                                                          *
C****************************************************************************
C
C         third order spline interpolation
C input argument and function:  X1(1:N1),Y1(1:N1)
C output argument and function: X2(1:N2),Y2(1:N2)
C the necessary conditions are: X1(I) < X1(I+1), and the same for X2 array.
C
C called by cco2gr
C calls nothing

      SUBROUTINE A18INT(X1,Y1,X2,Y2,N1,N2)
      DIMENSION X1(N1),X2(N2),Y1(N1),Y2(N2)
     *,A(150),E(150),F(150),H(150)
      H2=X1(1)
      NVS=N1-1
      DO 1 K=1,NVS
      H1=H2
      H2=X1(K+1)
      H(K)=H2-H1
    1 CONTINUE
      A(1)=0.
      A(N1)=0.
      E(N1)=0.
      F(N1)=0.
      H1=H(N1-1)
      F1=Y1(N1-1)
      F2=Y1(N1)
      DO 2 KR=2,NVS
      K=NVS+2-KR
      H2=H1
      H1=H(K-1)
      F3=F2
      F2=F1
      F1=Y1(K-1)
      G=1/(H2*E(K+1)+2.*(H1+H2))
      E(K)=-H1*G
      F(K)=(3.*((F3-F2)/H2-(F2-F1)/H1)-H2*F(K+1))*G
    2 CONTINUE
      G=0.
      DO 3 K=2,NVS
      G=E(K)*G+F(K)
      A(K)=G
    3 CONTINUE
      L=1
      DO 4 I=1,N2
      G=X2(I)
      DO 6 K=L,NVS
      IF(G.GT.X1(K+1))GOTO6
      L=K
      GOTO 5
    6 CONTINUE
      L=NVS
    5 G=G-X1(L)
      H2=H(L)
      F2=Y1(L)
      F1=H2**2
      F3=G**2
      Y2(I)=F2+G/H2*(Y1(L+1)-F2-(A(L+1)*(F1-F3)+
     *               A(L)*(2.*F1-3.*G*H2+F3))/3.)
    4 CONTINUE
      RETURN
      END















