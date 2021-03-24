C****************************************************************************
C*                                                                          *
C*                       SUBROUTINE  PCOOLG                                 *
C*                                                                          *
C****************************************************************************
C
C  to calculate the heating rate in both, 15 um CO2 and 9.6 um O3 bands.
C  NLTE conditions are taken into account for CO2 band. Irregular vertical
C  grid to account for heat exchange in the "LTE-layer" (x=2-12.5), is used.
C  The subroutine can be used for an arbitrary CO2 volume mixing ratio
C  with the parameterization coefficients calculated by the CCO2GR
C
C                                    December, 1998.  V.I. Fomichev.

C Called by a driving programm
C Calls nothing

      SUBROUTINE PCOOLG(H,X,T,O3,CO2,O2,SN2,O,AM,ZCO2O,FLUX)

C input:
C T(67)  - temperature (K) at x(=0-16.5, 0.25) levels
C X(67)  - pressure scale height array; X from 0 up to 16.5, step 0.25
C O3(45) - ozone volume mixing ratio (vmr) at x=0-11., 0.25
C CO2,O2,N2,O(17) - volume mixing ratios (vmr's) for corresponding gases
C                   at X = 12.5 - 16.5, with a step of 0.25
C AM(17) - air molecular weight at X = 12.5 - 16.5, 0.25
C ZCO2O  - CO2-O collisional deactivation rate constant (in cm3/sec)
C parameterization coeficients for CO2 (calculated by CCO2GR) come from
C the CO2CFG and PIRO3 common blocks. Irregular vertical grid (ig(9)) 
c comes from the BLOCK DATA PCO2O3 (PIRMGR common block)

C output: 
C  H(59) - heating rates in erg/g/sec at X = 2-16.5, 0.25
C  FLUX - upward flux at X=16.5 to calculate the heating rate above this level.
C
C Note: 1) Actually, any arbitrary vertical grid could be utilized starting
C          up from x=13.75. To do so, the coefficients for the reccurence 
C          formula must be determined at a proper grid (using CCO2GR) and all
C          input information should be given at this grid too.
C       2) As of now, ZCO2O is recommended to vary in a range of (1.5-6)e-12
C          with a mean value of 3e-12 cm3/sec, without T-dependence
C       3) To get the heating rate values in K/sec, the heating rates 
C          calculated here should be divided by Cp - specific heat at constant
C          pressure. Where Cp = R/AM * {7/2*SUM[Rv,i - for 2-atom gases] +
C          5/2*SUM[Rv,i - for 1-atom gases]}, R - gas constant, Rv - vmr's.
C          

      REAL H(59),X(67),T(67),O3(45),CO2(17),O2(17),SN2(17),O(17),AM(17)

      common /CO2CFG/ amat(43,9),bmat(43,9),al(17)
      common /PIRMGR/ ig(9)
      common /PIRO3/  ao3(35,9)

C INTERNAL arrays:

C su(67) - exponential part of the Planck function for CO2 band at X(67) grid
C so3(49) - the same for O3 band at x=(0-11,0.25) and 0 for x=11.25-11.5
C lambda(17) - quantum survival probability at the X grid for levels 
C              starting up from x=12.5 (indexes 51-67)
 
      real su(67), so3(49), lambda(17)

      data so3/49*0./

c const - constant used to determine the heating rate:
c         const = Na*h*c*V*Gv'/Gv*Avv, where
c         Na - Avogadro number, h- Planck's constant, c- light speed, 
c         V - frequency of the fundamental vibrational transition for the
c         main isotope, Gv'/Gv = 2 - ratio of the statistical weights for the
c         fundamental transition, Av'v - Einstine coefficient for the
c         fundamental band of the main isotope.
c constb = const/28.96 - to determine a boundary condition for the reccurence
c         formula at x=12.5 (air molecule weight is supposed to be 28.96
c         below x=12.5... could be changed, if neccessary)
c boltz - Boltzmann's constant
c a10   - Einstine's coefficient for the fundamental band for the main isotope

      data const/2.55521e11/, constb/8.82325e9/,

     :     boltz/1.38066e-16/, a10/1.5988/

c aku = h*c/k*Vco2 - to determine the Planck function, k - Boltzmann's constant
c ako3 - the same as aku, but for O3

      data aku/-960.217/, ako3/-1500./

c to determine su(67), so3(45):

      do 1 i = 1,67
      su(i)=exp(aku/T(i))
      if(i.ge.46) goto 1
      so3(i)=exp(ako3/T(i))
    1 continue

c to determine lambda(17)

      do 2 i=1,17

c --- CO2-O2 and CO2-N2 V-T constants:
      tt=T(i+50)
      y=tt**(-1./3.)
      zn2=5.5e-17*sqrt(tt)+6.7e-10*exp(-83.8*y)
      zo2=1.e-15*exp(23.37-230.9*y+564.*y*y)

c --- c air number density:
      d=1.e6*exp(-X(i+50))/(tt*boltz)

c --- collisional deactivation rate:
      z=(SN2(i)*zn2+O2(i)*zo2+O(i)*ZCO2O)*d

      lambda(i) = a10/(a10+z)

    2 continue

c cooling rate in both O3 and CO2 bands (x=2-10.5, the matrix approach is used)

      do 3 i=1,5
      is = i+8
      h2=(amat(i,1)+bmat(i,1)*su(is))*su(1)
      h3=ao3(i,1)*so3(1)
      do 4 j=3,9
      js=is+ig(j)
      h2=h2+(amat(i,j)+bmat(i,j)*su(is))*su(js)
    4 h3=h3+ao3(i,j)*so3(js)
    3 H(i) = h2+h3*O3(is)

      do 5 i=6,18
      is = i+8
      h2=(amat(i,1)+bmat(i,1)*su(is))*su(1)
      h3=ao3(i,1)*so3(1)
      do  6 j=2,9
      js=is+ig(j)
      h2=h2+(amat(i,j)+bmat(i,j)*su(is))*su(js)
    6 h3=h3+ao3(i,j)*so3(js)
    5 H(i) = h2+h3*O3(is) 

      do 7 I=19,35
      is=i+8
      h2=0.
      h3=0.
      do 8 j=1,9
      js=is+ig(j)
      h2=h2+(amat(i,j)+bmat(i,j)*su(is))*su(js)
    8 h3=h3+ao3(i,j)*so3(js)
    7 H(i) = h2+h3*O3(is)

c cooling rate in CO2 bands (x=10.75-12.5, the matrix approach is used)

      do 9 I=36,43
      is=i+8
      h2=0.
      do 10 j=1,9
      js=is+ig(j)
   10 h2=h2+(amat(i,j)+bmat(i,j)*su(is))*su(js)
    9 H(i) = h2

c calculate the heating rate for x=12.75-16.5 (the reccurence formula is used)

c --- to form the boundary condition at x=12.5
      h1=H(43)/(CO2(1)*(1.-lambda(1))*constb)

c --- the reccurence formula
      do 11 I=2,17
      im=i-1
      aa1=1.-lambda(im)*(1.-.25*AL(i)-.75*AL(im))
      aa2=1.-lambda(i)*(1.-.75*AL(i)-.25*AL(im))
      d1=-.25*(AL(i)+3.*AL(im))
      d2=.25*(3.*AL(i)+AL(im))
      h2=(aa1*h1-d1*su(im+50)-d2*su(i+50))/aa2
      H(i+42)=h2*CO2(i)*(1.-lambda(i))/AM(i)*const
      h1=h2
   11 continue

c to determine FLUX
c cooling rate above x=16.5 is suggested to be calculated by the formula
c           H(i) = const/AM(i)*CO2(i)*(1.-lambda(i))*(FLUX-su(i))
c no any parameterization coefficients are needed and an arbitrary hight grid
c can be used above x=16.5 level

      FLUX = h2 + su(67)

      return
      end
