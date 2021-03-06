C!!!!!!!!!!! - in this version the cooling does not depend on longitude????
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     ttau is introduced
c     Alex Pogoretsev, January 2006
c--------------------------------------------------------------------------
      subroutine ircool_60
 
c     XX-level version feb'95 (quelle: 118-level version'91)
c     (beachte aenderungen gegenueber alter version !!!!)
c
      parameter (kh2o=25, ngit=118)
      include 'com_main.fc'
c      include 'com_namelist.fc' !FL
      include 'data_irc_LIM.f'
      include 'udssr.f'
c--------------------------------------------------------------------
c   included by lange(october 1999)
c      real rmrfdmol(nb,73)
c--------------------------------------------------------------------
c   included by VPO(September 1999)
C
      real Tkinet(81),XGRDFUL(81),XXHELP(ngit), COOLCO2(ngit)
      real su(81),lambda(17)
c--------------------------- new parameterization for H2O cooling 
c--------------- Chou et al., JAS, 1993, Vol. 50, No. 14, p. 2294
      REAL h2oA(10),h2oB(10),h2oK1(10),h2oN(10),CSLAVE(13,6),
     $h2oC1(10),h2oC2(10),h2oC3(10),h2oC4(10),h2oC5(10),h2oC6(10),
     $     h2onue(10),dh2onu(10),e_co2(kh2o)
      REAL co2A(3),co2B(3),co2K1(3),co2N(3),co2M(3),co2Pr(3),
     $ co2C1(3),co2C2(3),co2C3(3),co2C4(3),co2C5(3),co2C6(3)
      DATA h2onue/170.,440.,580.,670.,760.,890.,1040.,1240.,
     $            1640.,2450./
      DATA dh2onu/340.,200.,80.,100.,80.,180.,120.,280.,520.,1100./
      DATA h2oTr/250./,co2Pr/300.,30.,300./,
     $                               co2M/0.5,0.85,0.5/,co2N/3*8./
      DATA h2oA /.0021,.0140,.0149,.0199,.0231,.0302,.0307,.0154,
     $           .0008,.0096/
      DATA h2oB /-1.01e-5,5.57e-5,6.20e-5, 5.57e-5,1.70e-4,
     $            2.96e-4,2.86e-4,7.53e-5,-3.52e-6,1.64e-5/ 
      DATA h2oK1/1.78e+1,2.51e-1,6.40e-2,8.00e-3,1.00e-3,
     $           3.16e-4,3.16e-4,1.41e-3,7.95e-1,3.16e-4/
      DATA h2oN /2*6.,3*8.,2*6.,8.,6.,16./
      DATA h2oC1/.2747,.1521,.3153,.2375,.2180,.4654,.5543,.1845,
     $           .0740,.1436/
      DATA h2oC2/.2717,.3974,.4604,.4309,.4468,.2991,.2723,.2732,
     $           .1636,.2197/
      DATA h2oC3/.2752,.1778,.1326,.2376,.2214,.1343,.1131,.2353,
     $           .4174,.3185/
      DATA h2oC4/.1177,.1826,.0798,.0482,.0802,.0646,.0443,.1613,
     $           .1783,.2351/
      DATA h2oC5/.0352,.0374,.0119,.0458,.0253,.0226,.0160,.1146,
     $           .1101,.0647/
      DATA h2oC6/.0255,.0527,.0000,.0000,.0083,.0140,.0000,.0311,
     $           .0566,.0184/

      DATA co2A/.0179,.0042,.0184/,co2B/1.02e-4,2.e-5,1.12e-4/
      DATA co2K1/1.6e-5,1.6e-3,1.6e-5/
      DATA co2C1/.2673,.1970,.1784/,co2C2/.2201,.3528,.2432/,
     $     co2C3/.2106,.3056,.3086/,co2C4/.2409,.0861,.1983/,
     $     co2C5/.0196,.0434,.0424/,co2C6/.0415,.0151,.0291/ 
c----------- END of DATA for new paprameterization of H2O cooling
      common/CO2_PARA/ pco2(17),pam(17),pxnlte(17),
     &  HLOGPRES(93),HCO2VMR(73),HN2VMR(93),HO2VMR(93),HOVMR(93),
     &  HAVWEI(73),RFDMOL(73)
      common /PIRGRD/ XR(67),XRH(59),XRF(17)
      common /PIRFEC/ IREP
      common /CO2CFG/ amat(43,9),bmat(43,9),al(17)
      common /PIRMGR/ ig(9)
C---------------------------------------------------------------------
c....FUNKTIONEN
c     formel-funktion 1: h2o-saeule [g/cm*cm]
c
c------------------ 2. = (1000mb/500mb)
c      sh2o (afkt1,bfkt1,ajfkt1,e1fkt1,hfkt1) = 2.*(
c     *     ajfkt1/e1fkt1*(exp(afkt1*e1fkt1) - exp(bfkt1*e1fkt1)) +
c     * 3.819e-4*(-hfkt1)*(exp(-afkt1/hfkt1) - exp(-bfkt1/hfkt1)))
c     * 2.542e-4*(-hfkt1)*(exp(-afkt1/hfkt1) - exp(-bfkt1/hfkt1)))
c
c     formel-funktion 2
c
c       function planck(x,t)
c
c      einheit  w/m**2/sr*cm
c

c       c1=1.1911e-8
c       c2=1.439
c
c       planck=c1*x*x*x/(exp(c2*x/t)-1.)
c
c       return
c       end
c
c      planck (xfkt2,tfkt2) = 1.1911e-8*xfkt2*xfkt2*xfkt2/
c     *                       (exp(1.439*xfkt2/tfkt2)-1.)
c
c
c     formel-funktion 3
c
c      soloa (afkt3,cfkt3,xfkt3,yfkt3,abfkt3) =
c     *      2.*afkt3*cfkt3*abfkt3/(1.013e5*(sqrt(xfkt3*xfkt3-
c     *      yfkt3*yfkt3)+cfkt3*(xfkt3*xfkt3-yfkt3*yfkt3)))
c 
c     formel-funktion 4
c
c      Tipol(xz1,xz2,xzi,Tz1,Tz2)=Tz1+(xzi-xz1)*(Tz2-Tz1)/(xz2-xz1)
c
c-----------------     time for initialization
      ttau = float(nsec)/5./86400.
      tsec=(1.-exp(-ttau*ttau))/86400.
c--      print*,'ircool: tsec umschreiben'
c     intervalle: o3 
      DO iband=1,10
        CSLAVE(iband,1)=h2oC1(iband)
        CSLAVE(iband,2)=h2oC2(iband)
        CSLAVE(iband,3)=h2oC3(iband)
        CSLAVE(iband,4)=h2oC4(iband)
        CSLAVE(iband,5)=h2oC5(iband)
        CSLAVE(iband,6)=h2oC6(iband)
      END DO
      DO iband=11,13
        ibandm=iband-10
        CSLAVE(iband,1)=co2C1(ibandm)
        CSLAVE(iband,2)=co2C2(ibandm)
        CSLAVE(iband,3)=co2C3(ibandm)
        CSLAVE(iband,4)=co2C4(ibandm)
        CSLAVE(iband,5)=co2C5(ibandm)
        CSLAVE(iband,6)=co2C6(ibandm)
      END DO
      k1_o3 =int((22000.-zunten)/dz) +1   !1.level ! auf Modellgitter
      k2_o3 =int((71000.-zunten)/dz) +1   !n.level
c-------------------------------------------------------------------------------
c included by VPO (September 1999)

c const - constant used to determine the heating rate:
c         const = Na*h*c*V*Gv'/Gv*Avv, where
c         Na - Avogadro number, h- Planck's constant, c- light speed, 
c         V - frequency of the fundamental vibrational transition for the
c         main isotope, Gv'/Gv = 2 - ratio of the statistical weights for the
c         fundamental transition, Av'v - Einstein coefficient for the
c         fundamental band of the main isotope.
c constb = const/28.96 - to determine a boundary condition for the reccurence
c         formula at x=12.5 (air molecule weight is supposed to be 28.96
c         below x=12.5... could be changed, if neccessary)
c boltz - Boltzmann's constant
c a10   - Einstein coefficient for the fundamental band for the main isotope

      data const/2.55521e11/, constb/8.82325e9/,
     *     boltz/1.38066e-16/, a10/1.5988/
c aku = h*c/k*Vco2 - to determine the Planck function, k - Boltzmann's constant
      data aku/-960.217/

c----------------------------------------------------------------------
      rgit  = float(igit)
c      OPEN(1111,FILE='tst_cool.dat')
      do 720 il=1,igit    !***********     laengenschleife    ***********
      do 720 ib=1,nb      !***********     breitenschleife    ***********
         tbott    = T_1000(ib,il)
         a1h2o_ib = a1h2o (il,ib)
c--     ZONAL MEANS / LONGITUDE 
c      tbott    =0.
c      a1h2o_ib =0.
c      do il=1,igit
c         if(msurf.eq.2) tunten=tbot1(ib,il) !prognose
c         tbott    = tunten
c         tbott    =tbott    +T_1000(ib,il)/rgit !surface (z=0)
c         a1h2o_ib =a1h2o_ib +a1h2o (il,ib)/rgit
c      enddo
      do k=2,kgit+1
      tequib(k)= tn1(ib,k-1,il)
c         tequib(k)=0.
c         do il=1,igit
c            tequib(k)=tequib(k)+tn1(ib,k-1,il)/rgit !k=1,kgit
c         enddo
      enddo

      tequib(1)=tbott*2.-tequib(2)

      do 500 k =1,kgit+1
         x(k)=zunten+float(k-2)*dz !log-press.hoehe/m
         y(k)=tequib(k)
 500  continue

c.....INTERPOLATION -->strahlungs-level

c.....XX ->118 (CO2,COXO,NO)

      do 509 k=1,ngit
	     k1=ip118var(k)
c            zi aus irc118
	     T118(k)=Tipol(x(k1),x(k1+1),zi(k)*1000.,y(k1),y(k1+1))
 509  continue

c....XX ->72-H2O

      zh2o(1)=0. !surface
      do 1000 k=1,72
	      k1=ipolH2O(k)
              zh2o(k)=zh2o(1)+float(k-1)*zi(1)
	      TH2O(k)=Tipol(x(k1),x(k1+1),zh2o(k)*1000.,y(k1),y(k1+1))
 1000         continue

c....XX ->72-COX (O3)

      zint(1)=0.
      do 1002 k=1,72
	      k1=ipolCOX(k)
              zint(k)=zint(1)+float(k-1)*0.25*skh
	      T(k)   =Tipol(x(k1),x(k1+1),zint(k)*1000.,y(k1),y(k1+1))
 1002    continue
c
c....PARAMETRISIERUNGEN
c....O3
         do 11 i=1,33
             mi2=mina2(i)
            do 2 j=mi2,7
            sum=0.
              do 1 j1=j,7
             sum=sum+(co3(index2(i,j1+1),ib)*exp(-xlnp(index2(i,j1+1)))
     *              + co3(index2(i,j1),ib)*exp(-xlnp(index2(i,j1))))
     *              *(xlnp(index2(i,j1+1))-xlnp(index2(i,j1)))
 1            continue
            sum=sum*0.5
            bl(j)=ao3(i,j)*exp(-bo3(i,j)*sum)
 2          continue
         bl(8)=1.
            ma2=maxa2(i)
            do 4 j=9,ma2
            sum=0.
              do 3 j1=9,j
              sum=sum+(co3(index2(i,j1),ib)*exp(-xlnp(index2(i,j1)))
     *             + co3(index2(i,j1-1),ib)*exp(-xlnp(index2(i,j1-1))))
     *             * (xlnp(index2(i,j1))-xlnp(index2(i,j1-1)))
 3            continue
            sum=sum*0.5
            bl(j)=ao3(i,j)*exp(-bo3(i,j)*sum)
 4          continue
c
c
c
             cpara(mina2(i))=bl(mina2(i)+1)+bl(mina2(i))
             mi2=mina2(i)+1
             do 6 j=mi2,7
                cpara(j)=bl(j+1)-bl(j-1)
 6           continue
             if(maxa2(i).eq.8)  cpara(8)=-3.-bl(7)
             if(maxa2(i).eq.8)                           goto 9
             cpara(8)=-2.-bl(9)-bl(7)
             if(maxa2(i).eq.9)                           goto 8
             ma2=maxa2(i)-1
             do 7 j=9,ma2
                cpara(j)=bl(j-1)-bl(j+1)
 7           continue
 8           cpara(maxa2(i))=bl(maxa2(i)-1)-bl(maxa2(i))
 9         continue
c
c
              sum=0.
              mi2=mina2(i)
              ma2=maxa2(i)
              do 10 j=mi2,ma2
                 aexp=exp(-1500./t(index2(i,j)))
                 sum=sum+cpara(j)*aexp
 10           continue
              epso3(i)=sum*2.37e4*co3(index2(i,8),ib)*0.5*86400./
     *                    1.005e7  ! [cp] = cm**2/(s**2 K)
 11        continue

c
c
c......CO2
c (New parameterization for 15-micron CO2 cooling for both LTE and non-LTE)

c----------------------------------------------------------------
c included by VPO(September 1999)

C  To calculate the heating rate in 15 um CO2 band.
C  NLTE conditions are taken into account for CO2 band. Irregular vertical
C  grid to account for heat exchange in the "LTE-layer" (x=2-12.5), is used.
C  The subroutine can be used for an arbitrary CO2 volume mixing ratio
C  with the parameterization coefficients calculated by the CCO2GR
C
C!!!!!!!!!!!!!!!!! WAS 3.0 !!!!!!!!!!!!!!!!
      ZCO2O= 1.5e-12 ! quenching rate constant [cm^3/s]

C--  Determination the full grid XGRDFUL
      do i=1,8
        XGRDFUL(i)=XR(i)
      enddo
      do i=9,81
        ii=i-8
        XGRDFUL(i)=HLOGPRES(ii)
      enddo

      do k=1,ngit
        xxhelp(k)=zi(k)/(H/1000.0)
c        print*,'k=',k,'  zi=',zi(k),'  xxhelp=',xxhelp(k)
      enddo

c...................  To interpolate T118 at XGRDFUL to get Tkinet(81)
      do i=1,81
        xxx=XGRDFUL(i)
        if(xxx.ge.xxhelp(1).and.xxx.le.xxhelp(ngit)) then
          Tkinet(i)=A18LIN(xxx,xxhelp,T118,1,NGIT)
          KBD2=i
        else if(xxx.lt.xxhelp(1)) then
          Tkinet(i)=T118(1)
        else
          Tkinet(i)=T118(ngit)
        endif
c        print*,'i=',i,'  XGRDFUL=',XGRDFUL(i),'  Tkinet=',Tkinet(i)
      enddo
c      print*, 'KBD2=',KBD2

      KBD2P1=KBD2+1
      KBD2M8=KBD2-8
      
      do i=1,ngit
        KBD1=i
        if(xxhelp(i).gt.16.5) go to 202
      enddo
  202  continue

c       print*, 'KBD1=',KBD1
       KBD1P1=KBD1+1
          
c...... to determine su:

      do  i = 1,KBD2P1
        su(i)=exp(aku/Tkinet(i))
      enddo

      do i=1,73
        RFDMOL(i)=0.0
      enddo

c...... to determine lambda(17)

      do i=1,17

c --- CO2-O2 and CO2-N2 V-T constants:
      tt=Tkinet(i+50)
      yhlp=tt**(-1./3.)
      zn2=5.5e-17*sqrt(tt)+6.7e-10*exp(-83.8*yhlp)
      zo2=1.e-15*exp(23.37-230.9*yhlp+564.*yhlp*yhlp)

c --- c air number density:
      dpres=1.013e6*exp(-HLOGPRES(i+42))/(tt*boltz)

c --- collisional deactivation rate:
      zhlp=(HN2VMR(i+42)*zn2+HO2VMR(i+42)*zo2+HOVMR(i+42)*ZCO2O)*dpres
      lambda(i) = a10/(a10+zhlp)
      enddo
C
c................ cooling rate  (x=2-10.5, the matrix approach is used)

      do  i=1,5
          is = i+8
          RFDCO2=(amat(i,1)+bmat(i,1)*su(is))*su(1)
          do  j=3,9
              js=is+ig(j)
               RFDCO2=RFDCO2+(amat(i,j)+bmat(i,j)*su(is))*su(js)
          enddo
          RFDMOL(i) = RFDCO2
      enddo

      do  i=6,18
          is = i+8
          RFDCO2=(amat(i,1)+bmat(i,1)*su(is))*su(1)
          do  j=2,9
              js=is+ig(j)
              RFDCO2=RFDCO2+(amat(i,j)+bmat(i,j)*su(is))*su(js)
          enddo
          RFDMOL(i) = RFDCO2
      enddo

c.......cooling rate in CO2 bands (x=10.75-12.5, the matrix approach is used)
      do  I=19,43
          is=i+8
          RFDCO2=0.
          h3=0.
          do  j=1,9
              js=is+ig(j)
              RFDCO2=RFDCO2+(amat(i,j)+bmat(i,j)*su(is))*su(js)
          enddo
          RFDMOL(i) = RFDCO2
      enddo

c      do i=1,43
c        print*, 'i=',i,'  HLOGPRES=',HLOGPRES(i),'  RFDMOL=', RFDMOL(i)
c      enddo

c calculate the heating rate for x=12.75-16.5 (the reccurence formula is used)

c --- to form the boundary condition at x=12.5
      RFDBORD=RFDMOL(43)/(PCO2(1)*(1.-lambda(1))*constb)

c --- the reccurence formula
      do  I=2,17
        im=i-1
        aa1=1.-lambda(im)*(1.-.25*AL(i)-.75*AL(im))
        aa2=1.-lambda(i)*(1.-.75*AL(i)-.25*AL(im))
        ddd1=-.25*(AL(i)+3.*AL(im))
        ddd2=.25*(3.*AL(i)+AL(im))
        RFDCO2=(aa1*RFDBORD-ddd1*su(im+50)-ddd2*su(i+50))/aa2
        RFDMOL(i+42)=RFDCO2*PCO2(i)*(1.-lambda(i))/PAM(i)*const
        RFDBORD=RFDCO2
c        print*, '(i-42)=',i,'  HLOGPRES=',HLOGPRES(i),
c     &    '  RFDMOL=', RFDMOL(i+42)
      enddo

c to determine FLUX
c cooling rate above x=16.5 is suggested to be calculated by the formula

c           RFDMOL(i) = const/PAM(i)*PCO2(i)*(1.-lambda(i))*(FLUX-su(i))
c no any parameterization coefficients are needed and an arbitrary hight grid
c can be used above x=16.5 level

      FLUX = RFDBORD+su(67)

C
      do  i=1,73
c --- CPCONST - specific heat at constant pressure
      if(i.le.42) then
        CPCONST = 7./2.*8.31441e7/HAVWEI(i)
      else
        CPCONST = 8.31441e7/HAVWEI(i)*
     &    (7./2.*(1.-HOVMR(i))+5./2.*HOVMR(i))
      endif
c --- to convert heating rate (at x=2-16.5) into "K/day"
      if(i.le.59) then
           RFDMOL(i) = RFDMOL(i)*86400./CPCONST
      else
c ***** to calculate heating rate above x=16.5
C
c determine quantum survival probability -- alambda:
C
c --- V-T constants for O2 and N2
         tt=Tkinet(i+8)
         yhlp=tt**(-1./3.)
         zn2=5.5e-17*sqrt(tt)+6.7e-10*exp(-83.8*yhlp)
         zo2=1.e-15*exp(23.37-230.9*yhlp+564.*yhlp*yhlp)
C
c --- air number density:
         ddhlp=1.013e6*exp(-HLOGPRES(i))/(tt*boltz)
c --- collisional deactivation rate:
         zcollis=(HN2VMR(i)*zn2+HO2VMR(i)*zo2+HOVMR(i)*ZCO2O)*ddhlp
C
         alambda=a10/(a10+zcollis)

c --- source function
         sour=exp(aku/Tkinet(i+8))

c --- cooling rate above x=16.75 in K/day
         RFDCO2=const/HAVWEI(i)*HCO2VMR(i)*(1.-alambda)*(FLUX-sour)
         RFDMOL(i)=RFDCO2*86400./CPCONST
       end if
c ---  print*, 'i=',i,'  HLOGPRES=',HLOGPRES(i),'  RFDMOL=', RFDMOL(i)
      enddo

c --- now to interpolate RFDMOL(73) to the grid of given problem
c --- interpolation on the internal grid
c      IF(ib.eq.1.and.il.eq.1) then
c      do kint=5,34
c      kk=kint-4
c      write(1111,1112) HLOGPRES(kint),tkinet(kint+8),
c     $                 epso3(kk),RFDMOL(kint)
c 1112 format(1x,4e12.3)
c      end do
c      close(1111)
c      stop
c      end if
      do i=1,ngit
        COOLCO2(i)=0.0
      enddo

      do i=1,ngit   !! interpolation in the range x=2-20,0
        xxx=xxhelp(i)
        if(xxx.ge.HLOGPRES(1).and.xxx.le.HLOGPRES(73)) then
          COOLCO2(i)=A18LIN(xxhelp(i),HLOGPRES,RFDMOL,1,73)
        endif
      enddo

c      do i=1,73
c         rmrfdmol(ib,i)=rfdmol(i)
c      enddo
c      do i=1,ngit
c        print*, 'i=',i,'  XXHELP=',xxhelp(i),'  COOLCO2=', COOLCO2(i)
c      enddo
      
C if necessary to convert RFDMOL(i) from K/day to another units

cc ******* to print out results *******
c  605 format((1x,i2,1x,f5.2, 1x, 12(1x,f7.2)))
c      do 630 i=1,73
c        print 605, i,HLOGPRES(i), RFDMOL(i)
c  630 continue

c------------------------------------------------------------------------------
c
c.....H2O
             diffus = 1.66
             satmcm =7.*0.1*co2pmv
c--- PLANCK
	    do 400 inu=1,10  ! nu
	    do 400 ih2o=1,57  ! T
 400        bnue(inu,ih2o) = pi*planck(h2onue(inu),th2o(ih2o))  ! pi*B(nu,T)
      nh2o=kh2o+kh2o
      do 440 i=1,kh2o    ! H2O-level
              ii=min(4,1+i/5)
              ih2o=i+i
           h2o_LW = 0.
	     h2o_UP = 0.
	     co2_LW = 0.
	     co2_UP = 0.
c---  I)  COOL TO EARTH = -[B(0)-B(dz/2)]*dT(z,bot)/dz
c---------------------------------------------------- H2O
      u = diffus*sh2o(zh2o(ih2o),0.,a1h2o_ib,e1nir(ii,ib),skhnir(ii))
       deltuT=th2o(1   )-h2oTr
       deltaT=th2o(ih2o)-h2oTr
      D_1000=1.013e6*rm(1)/RgSGS/tbott
      dudz = 2.*diffus*1.273*roi(i)*(a1h2o_ib/D_1000/1.e5*
     $                 exp(-zi(i)/H_wv(ib))+3.e-6)*roi(i)**xnc(ii)
      do 403 inu=1,10
      if(inu.lt.3.or.inu.gt.5) go to 404
c----------------------------------------------------- CO2
      inuco2=inu-2
      hkmco2 = 7./(1.+co2M(inuco2))
      uco2 =diffus*satmcm/(1.+co2M(inuco2))*(1.-EXP(-zh2o(ih2o)/hkmco2))
      dudz_c=diffus*satmcm/7.e3*EXP(-zh2o(ih2o)/hkmco2)
      F_co2U=(1000./co2Pr(inuco2))**co2M(inuco2)*
     $       (1.+co2A(inuco2)*deltuT+co2B(inuco2)*deltuT*deltuT)
      F_co2Z=(1000./co2Pr(inuco2))**co2M(inuco2)*
     $       (1.+co2A(inuco2)*deltaT+co2B(inuco2)*deltaT*deltaT)
      W_co2=uco2*F_co2U
       sum=0.
       co2KN=co2K1(inuco2) 
      DO in=1,6
       sum=sum-
     $     CSLAVE(inuco2+10,in)*EXP(-co2KN*W_co2)*co2KN*F_co2Z*dudz_c
       co2KN=co2N(inuco2)*co2KN
      END DO
       co2_LW=co2_LW+sum*0.5*(bnue(inu,2)-bnue(inu,1))*dh2onu(inu)
      GO TO 403
 404  continue
       F_h2oU=1.+h2oA(inu)*deltuT+h2oB(inu)*deltuT*deltuT
       F_h2oZ=1.+h2oA(inu)*deltaT+h2oB(inu)*deltaT*deltaT
       W_h2o=u*F_h2oU
       sum=0.
       h2oKN=h2oK1(inu) 
      DO in=1,6
       sum=sum-CSLAVE(inu,in)*EXP(-h2oKN*W_h2o)*h2oKN*F_h2oZ*dudz
       h2oKN=h2oN(inu)*h2oKN
      END DO
       h2o_LW=h2o_LW+sum/10.*0.5*(bnue(inu,2)-bnue(inu,1))*dh2onu(inu)
 403  continue
c------------------------------------------------------ II)  LOWER SIDE
       if(ih2o.eq.2) go to 412
       do 410 ifo = 2,ih2o-2,2   ! 0 < z' < z
      u=diffus*tbott/th2o(ifo)*
     $  sh2o(zh2o(ih2o),zh2o(ifo),a1h2o_ib,e1nir(ii,ib),skhnir(ii))
       deltuT=th2o(ifo)-h2oTr
       sum_h2o=0.
       sum_co2=0.
      do 411 inu=1,10
       bnuelw=bnue(inu,ifo-1)
       if(ifo.eq.2) bnuelw=(bnue(inu,1)+bnue(inu,2))/2.
      if(inu.lt.3.or.inu.gt.5) go to 413
c----------------------------------------------------- CO2
      inuco2=inu-2
      hkmco2 = 7./(1.+co2M(inuco2))
      uco2 = diffus*satmcm/(1.+co2M(inuco2))*273./th2o(ifo)*
     $              (EXP(-zh2o(ifo)/hkmco2)-EXP(-zh2o(ih2o)/hkmco2))
      dudz_c=diffus*satmcm/7.e3*EXP(-zh2o(ih2o)/hkmco2)
      F_co2U=(1000./co2Pr(inuco2))**co2M(inuco2)*
     $       (1.+co2A(inuco2)*deltuT+co2B(inuco2)*deltuT*deltuT)
      F_co2Z=(1000./co2Pr(inuco2))**co2M(inuco2)*
     $       (1.+co2A(inuco2)*deltaT+co2B(inuco2)*deltaT*deltaT)
      W_co2=uco2*F_co2U
       sum=0.
       co2KN=co2K1(inuco2) 
      DO in=1,6
       sum=sum-
     $     CSLAVE(inuco2+10,in)*EXP(-co2KN*W_co2)*co2KN*F_co2Z*dudz_c
       co2KN=co2N(inuco2)*co2KN
      END DO
       sum_co2=sum_co2+sum*(bnue(inu,ifo+1)-bnuelw)*dh2onu(inu)
      GO TO 411
 413  continue
       F_h2oU=1.+h2oA(inu)*deltuT+h2oB(inu)*deltuT*deltuT
       F_h2oZ=1.+h2oA(inu)*deltaT+h2oB(inu)*deltaT*deltaT
       W_h2o=u*F_h2oU
       sum=0.
       h2oKN=h2oK1(inu) 
      DO in=1,6
       sum=sum-CSLAVE(inu,in)*EXP(-h2oKN*W_h2o)*h2oKN*F_h2oZ*dudz
       h2oKN=h2oN(inu)*h2oKN
      END DO
       sum_h2o=sum_h2o+sum/10.*(bnue(inu,ifo+1)-bnuelw)*dh2onu(inu)
 411   continue
      h2o_LW = h2o_LW +  sum_h2o
      co2_LW = co2_LW +  sum_co2
 410   continue
 412   continue
c--- III)  UPPER SIDE (including cooling to space)
      if(ih2o+2.gt.nh2o) go to 422
      do 420 ifo = ih2o+2,nh2o,2  ! z < z' < top
      u=diffus*tbott/th2o(ih2o)*
     $  sh2o(zh2o(ifo),zh2o(ih2o),a1h2o_ib,e1nir(ii,ib),skhnir(ii))
       sum_h2o=0.
       sum_co2=0.
      do 421 inu=1,10  ! nu
       bnueup=bnue(inu,ifo+1)
       if(ifo.eq.nh2o) bnueup=0.
      if(inu.lt.3.or.inu.gt.5) go to 423
c----------------------------------------------------- CO2
      inuco2=inu-2
      hkmco2 = 7./(1.+co2M(inuco2))
      uco2 = diffus*satmcm/(1.+co2M(inuco2))*273./th2o(ih2o)*
     $              (EXP(-zh2o(ih2o)/hkmco2)-EXP(-zh2o(ifo)/hkmco2))
      dudz_c=diffus*satmcm/7.e3*EXP(-zh2o(ih2o)/hkmco2)
      F_co2U=(1000./co2Pr(inuco2))**co2M(inuco2)*
     $       (1.+co2A(inuco2)*deltaT+co2B(inuco2)*deltaT*deltaT)
      W_co2=uco2*F_co2U
       sum=0.
       co2KN=co2K1(inuco2) 
      DO in=1,6
       sum=sum+
     $     CSLAVE(inuco2+10,in)*EXP(-co2KN*W_co2)*co2KN*F_co2U*dudz_c
       co2KN=co2N(inuco2)*co2KN
      END DO
       sum_co2=sum_co2+sum*(bnueup-bnue(inu,ifo-1))*dh2onu(inu)
      GO TO 421
 423  continue
       F_h2oU=1.+h2oA(inu)*deltaT+h2oB(inu)*deltaT*deltaT
       W_h2o=u*F_h2oU
       sum=0.
       h2oKN=h2oK1(inu) 
      DO in=1,6
       sum=sum+CSLAVE(inu,in)*EXP(-h2oKN*W_h2o)*h2oKN*F_h2oU*dudz
       h2oKN=h2oN(inu)*h2oKN
      END DO
       sum_h2o=sum_h2o+sum/10.*(bnueup-bnue(inu,ifo-1))*dh2onu(inu)
 421   continue
      h2o_UP = h2o_UP +  sum_h2o
      co2_UP = co2_UP +  sum_co2
 420   continue
 422   continue 
c--- IV)  contribution from the nearest layers
c------------------------------------------------------------------
      zh2om1=zh2o(ih2o-1)
      th2om1=th2o(ih2o-1)
      ih2om2=ih2o-2 
      if(ih2o.eq.2) then
        zh2om1=zi(1)/2. 
        th2om1=(th2o(1)+th2o(2))/2. 
        ih2om2=ih2o-1
      end if
      u=diffus*tbott/th2om1*
     $  sh2o(zh2o(ih2o),zh2om1,a1h2o_ib,e1nir(ii,ib),skhnir(ii))
      deltuT=th2om1-h2oTr
      sum_h2o=0.
      sum_co2=0.
      do 431 inu=1,10  ! nu
      if(inu.lt.3.or.inu.gt.5) go to 433
c----------------------------------------------------- CO2
      inuco2=inu-2
      hkmco2 = 7./(1.+co2M(inuco2))
      uco2 = diffus*satmcm/(1.+co2M(inuco2))*273./th2om1*
     $              (EXP(-zh2om1/hkmco2)-EXP(-zh2o(ih2o)/hkmco2))
      dudz_c=diffus*satmcm/7.e3*EXP(-zh2o(ih2o)/hkmco2)
      F_co2U=(1000./co2Pr(inuco2))**co2M(inuco2)*
     $       (1.+co2A(inuco2)*deltuT+co2B(inuco2)*deltuT*deltuT)
      F_co2Z=(1000./co2Pr(inuco2))**co2M(inuco2)*
     $       (1.+co2A(inuco2)*deltaT+co2B(inuco2)*deltaT*deltaT)
      W_co2=uco2*F_co2U
       sum=0.
       co2KN=co2K1(inuco2) 
      DO in=1,6
       sum=sum-
     $ CSLAVE(inuco2+10,in)*EXP(-co2KN*W_co2)*co2KN*F_co2Z*dudz_c
       co2KN=co2N(inuco2)*co2KN
      END DO
      sum_co2=sum_co2+sum*(bnue(inu,ih2o)-bnue(inu,ih2om2))*dh2onu(inu)
      GO TO 431
 433  continue
       F_h2oU=2.*(1.+h2oA(inu)*deltuT+h2oB(inu)*deltuT*deltuT)
       F_h2oZ=2.*(1.+h2oA(inu)*deltaT+h2oB(inu)*deltaT*deltaT)
       W_h2o=u*F_h2oU
       sum=0.
       h2oKN=h2oK1(inu) 
      DO in=1,6
       sum=sum-CSLAVE(inu,in)*EXP(-h2oKN*W_h2o)*h2oKN*F_h2oZ*dudz
       h2oKN=h2oN(inu)*h2oKN
      END DO
      sum_h2o=sum_h2o+
     $       sum/10.*(bnue(inu,ih2o)-bnue(inu,ih2om2))*dh2onu(inu)
  431 continue
      h2o_LW=h2o_LW + sum_h2o
      co2_LW=co2_LW + sum_co2
c---------------------------------------------------------------------
      u=diffus*tbott/th2o(ih2om2)*
     $  sh2o(zh2o(ih2o),zh2o(ih2om2),a1h2o_ib,e1nir(ii,ib),skhnir(ii))
      deltuT=th2o(ih2om2)-h2oTr
      sum_h2o=0.
      sum_co2=0.
      do 441 inu=1,10  ! nu
       bnuelw=bnue(inu,ih2o-1)
       if(ih2o.eq.2) bnuelw=(bnue(inu,1)+bnue(inu,2))/2.
      if(inu.lt.3.or.inu.gt.5) go to 443
c----------------------------------------------------- CO2
      inuco2=inu-2
      hkmco2 = 7./(1.+co2M(inuco2))
      uco2 = diffus*satmcm/(1.+co2M(inuco2))*273./th2o(ih2om2)*
     $              (EXP(-zh2o(ih2om2)/hkmco2)-EXP(-zh2o(ih2o)/hkmco2))
      dudz_c=diffus*satmcm/7.e3*EXP(-zh2o(ih2o)/hkmco2)
      F_co2U=(1000./co2Pr(inuco2))**co2M(inuco2)*
     $       (1.+co2A(inuco2)*deltuT+co2B(inuco2)*deltuT*deltuT)
      F_co2Z=(1000./co2Pr(inuco2))**co2M(inuco2)*
     $       (1.+co2A(inuco2)*deltaT+co2B(inuco2)*deltaT*deltaT)
      W_co2=uco2*F_co2U
       sum=0.
       co2KN=co2K1(inuco2) 
      DO in=1,6
       sum=sum-
     $ CSLAVE(inuco2+10,in)*EXP(-co2KN*W_co2)*co2KN*F_co2Z*dudz_c
       co2KN=co2N(inuco2)*co2KN
      END DO
      sum_co2=sum_co2+sum*(bnue(inu,ih2om2)-bnuelw)*dh2onu(inu)
      GO TO 441
 443  continue
       F_h2oU=2.*(1.+h2oA(inu)*deltuT+h2oB(inu)*deltuT*deltuT)
       F_h2oZ=2.*(1.+h2oA(inu)*deltaT+h2oB(inu)*deltaT*deltaT)
       W_h2o=u*F_h2oU
       sum=0.
       h2oKN=h2oK1(inu) 
      DO in=1,6
       sum=sum-CSLAVE(inu,in)*EXP(-h2oKN*W_h2o)*h2oKN*F_h2oZ*dudz
       h2oKN=h2oN(inu)*h2oKN
      END DO
      sum_h2o=sum_h2o+sum/10.*(bnue(inu,ih2om2)-bnuelw)*dh2onu(inu)
  441 continue
      h2o_LW=h2o_LW + sum_h2o
      co2_LW=co2_LW + sum_co2
c---------------------------------------------------------------------
      if(ih2o.ge.nh2o-2) go to 555 
      u=diffus*tbott/th2o(ih2o)*
     $  sh2o(zh2o(ih2o+1),zh2o(ih2o),a1h2o_ib,e1nir(ii,ib),skhnir(ii))
      sum_h2o=0.
      sum_co2=0.
      do 451 inu=1,10  ! nu
      if(inu.lt.3.or.inu.gt.5) go to 453
c----------------------------------------------------- CO2
      inuco2=inu-2
      hkmco2 = 7./(1.+co2M(inuco2))
      uco2 = diffus*satmcm/(1.+co2M(inuco2))*273./th2o(ih2o)*
     $              (EXP(-zh2o(ih2o)/hkmco2)-EXP(-zh2o(ih2o+1)/hkmco2))
      dudz_c=diffus*satmcm/7.e3*EXP(-zh2o(ih2o)/hkmco2)
      F_co2U=(1000./co2Pr(inuco2))**co2M(inuco2)*
     $       (1.+co2A(inuco2)*deltaT+co2B(inuco2)*deltaT*deltaT)
      W_co2=uco2*F_co2U
       sum=0.
       co2KN=co2K1(inuco2) 
      DO in=1,6
       sum=sum+
     $ CSLAVE(inuco2+10,in)*EXP(-co2KN*W_co2)*co2KN*F_co2U*dudz_c
       co2KN=co2N(inuco2)*co2KN
      END DO
      sum_co2=sum_co2+sum*(bnue(inu,ih2o+2)-bnue(inu,ih2o))*dh2onu(inu)
      GO TO 451
 453  continue
       F_h2oU=2.*(1.+h2oA(inu)*deltaT+h2oB(inu)*deltaT*deltaT)
       W_h2o=u*F_h2oU
       sum=0.
       h2oKN=h2oK1(inu) 
      DO in=1,6
       sum=sum+CSLAVE(inu,in)*EXP(-h2oKN*W_h2o)*h2oKN*F_h2oU*dudz
       h2oKN=h2oN(inu)*h2oKN
      END DO
      sum_h2o=sum_h2o+
     $       sum/10.*(bnue(inu,ih2o+2)-bnue(inu,ih2o))*dh2onu(inu)
  451 continue
      h2o_UP=h2o_UP + sum_h2o
      co2_UP=co2_UP + sum_co2
c-----------------------------------------------------------------------
      u=diffus*tbott/th2o(ih2o)*
     $   sh2o(zh2o(ih2o+2),zh2o(ih2o),a1h2o_ib,e1nir(ii,ib),skhnir(ii))
      sum_h2o=0.
      sum_co2=0.
      do 461 inu=1,10  ! nu
      if(inu.lt.3.or.inu.gt.5) go to 463
c----------------------------------------------------- CO2
      inuco2=inu-2
      hkmco2 = 7./(1.+co2M(inuco2))
      uco2 = diffus*satmcm/(1.+co2M(inuco2))* 273./th2o(ih2o)*
     $              (EXP(-zh2o(ih2o)/hkmco2)-EXP(-zh2o(ih2o+2)/hkmco2))
      dudz_c=diffus*satmcm/7.e3*EXP(-zh2o(ih2o)/hkmco2)
      F_co2U=(1000./co2Pr(inuco2))**co2M(inuco2)*
     $       (1.+co2A(inuco2)*deltaT+co2B(inuco2)*deltaT*deltaT)
      W_co2=uco2*F_co2U
       sum=0.
       co2KN=co2K1(inuco2) 
      DO in=1,6
       sum=sum+
     $ CSLAVE(inuco2+10,in)*EXP(-co2KN*W_co2)*co2KN*F_co2U*dudz_c
       co2KN=co2N(inuco2)*co2KN
      END DO
      sum_co2=sum_co2+
     $        sum*(bnue(inu,ih2o+1)-bnue(inu,ih2o+2))*dh2onu(inu)
      GO TO 461
 463  continue
       F_h2oU=2.*(1.+h2oA(inu)*deltaT+h2oB(inu)*deltaT*deltaT)
       W_h2o=u*F_h2oU
       sum=0.
       h2oKN=h2oK1(inu) 
      DO in=1,6
       sum=sum+CSLAVE(inu,in)*EXP(-h2oKN*W_h2o)*h2oKN*F_h2oU*dudz
       h2oKN=h2oN(inu)*h2oKN
      END DO
      sum_h2o=sum_h2o+
     $       sum/10.*(bnue(inu,ih2o+1)-bnue(inu,ih2o+2))*dh2onu(inu)
  461 continue
      h2o_UP=h2o_UP + sum_h2o
      co2_UP=co2_UP + sum_co2
  555 continue
c------------------------------------------------------------------
      dicpt2=86400./(1005.*dichti(i)*1000.)   !=tag/cp*rho(0)

c---         ERWAERMUNGSRATEN I+II+III+IV
          epsh1= dicpt2 * h2o_LW
          epsh2= dicpt2 * h2o_UP 
          epsh3= dicpt2 * co2_LW 
          epsh4= dicpt2 * co2_UP  
          epsh2o(i)= epsh1+epsh2
          e_co2 (i)= epsh3+epsh4
 440   continue
c.....INTERPOLATION --->XX level
      do k=1,kgit
	 cooli(k)=0.
      enddo
      do kint=1,ngit
	 cool(kint)=0.
      enddo

      do kint=13,42                  ! O3
	 kk=kint-12
	 cool(kint)=epso3(kk)
      enddo
 
c.....72-COX (O3) -> XX


c      do 503 k=19,62       ! O3
c--- check the altitude range!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do 503 k=k1_o3,k2_o3
             k1=ipCOXvar(k)
       cooli(k)=Tipol(zint(k1),zint(k1+1),z(k),cool(k1),cool(k1+1))
 503  continue

c.....118-->XX

c--------------------------------------------
c  modified by Alex Pogoreltsev (January, 2006)
      do k=1,20
        if(k.le.15) then
          cool(k)=e_co2(k)
        else
          cool(k)=(float(20-k)*e_co2(k)+float(k-15)*COOLCO2(k))/5.
        end if
      end do
      do k=21,ngit  
          cool(k)=COOLCO2(k) 
      enddo
c--------------------------------------------
       do k=1,kh2o
        cool(k)=cool(k)+epsh2o(k)      ! H2O
       enddo

c.....118 ->XX

c      do 508 k=1,kgit
      do 508 k=1,48
	     k1=ipVar118(k)
	     cooli(k)=cooli(k)+
     +              Tipol(zi(k1),zi(k1+1),z(k),cool(k1),cool(k1+1))
 508  continue
      do 519 k=48,kgit
             cooli(k)=cooli(48)*EXP((z(48)-z(k))/skh)
 519  continue
c      do 520 il=1,igit
      do 520 k=1,kgit
      fc(k,il,ib)=cooli(k)*tsec
 520  continue
 720  continue   !********** ende breitenschleife  **************
      return
      end
c------------------------------------------------------ Functions
c     formel-funktion 1: h2o-saeule [g/cm*cm]
c-------------------------------     see strobel_60.f for details 
c------------------ 2. = (1000mb/500mb)
c      module FUN_cool
      FUNCTION sh2o (afkt1,bfkt1,ajfkt1,e1fkt1,hfkt1)
      sh2o = 2.*(
     *     ajfkt1/e1fkt1*(exp(afkt1*e1fkt1) - exp(bfkt1*e1fkt1)) +
     * 3.819e-4*(-hfkt1)*(exp(-afkt1/hfkt1) - exp(-bfkt1/hfkt1)))
c     * 2.542e-4*(-hfkt1)*(exp(-afkt1/hfkt1) - exp(-bfkt1/hfkt1)))
      RETURN
      END
c----------------------------------------------------------------
c     formel-funktion 2
c
c       function planck(x,t)
c
c      einheit  w/m**2/sr*cm
c

c       c1=1.1911e-8
c       c2=1.439
c
c       planck=c1*x*x*x/(exp(c2*x/t)-1.)
c
c       return
c       end
c----------------------------------------------------------------
      FUNCTION planck (xfkt2,tfkt2)
      planck = 1.1911e-8*xfkt2*xfkt2*xfkt2/
     *                       (exp(1.439*xfkt2/tfkt2)-1.)
      RETURN 
      END
c
c--------------------------------- I didn't find where it is used
c--------------------------------- Alex Pogoreltsev 29.11.2007
c     formel-funktion 3
c
c      soloa (afkt3,cfkt3,xfkt3,yfkt3,abfkt3) =
c     *      2.*afkt3*cfkt3*abfkt3/(1.013e5*(sqrt(xfkt3*xfkt3-
c     *      yfkt3*yfkt3)+cfkt3*(xfkt3*xfkt3-yfkt3*yfkt3)))
c---------------------------------------------------------------- 
c     formel-funktion 4
      FUNCTION Tipol(xz1,xz2,xzi,Tz1,Tz2)
      Tipol=Tz1+(xzi-xz1)*(Tz2-Tz1)/(xz2-xz1) 
      RETURN
      END
c      end module FUN_cool


 


