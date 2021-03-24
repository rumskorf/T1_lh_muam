c--------------------------------------------------------------------
c     ttau is introduced to calculate eff_t, i.e., establishment
c     KGIT --> k48fix
c     Alex Pogoretsev, February 2006
c-------------------------------------------------------------------- 
      subroutine strobel_60
      include 'com_main.fc'
      dimension secchi(kgit),fsolar(kgit),hs_tide(kgit,igit,nb),
     * secnox(kgit),secno2(kgit),secno3(kgit),secnn2(kgit)
c---------  dimension phgdk(9,48),x0li(9),dnueli(9),fracd(9)
      dimension cli(9),dli(9),dkli(9),fracli(9),uliou(9),dudz(9),
     $          aliou0(9),uliou0(9),ulioub(9) 
      real sech2o(kgit,nb)
      double precision f_chem

csl........ ef_lange faktor mit dem die uv routine stimuliert wird
      real ef_lange
      include 'com_namelist.fc'
csl
      common /EUVFCS/ EUVflx(37),sEUVo2(37),sEUVn2(37),sEUVo(37)
      common /EUVI/ IEUV
      include 'data_strobel.f'    
CPH      ef_lange=1.05
         ef_lange=SOLE      
c---------------------------------     lokal:
c---- data for Ly alpha heating 
c     (flc_Ly is 5 for low and 10 for)
c--------------------------------------------------------------------
c
csl   effizienz geändert für solarmax 1. -> ef_s *ef_lange 2. -> ef_m *ef_lange (mit faktor)
csl      data ef_s/1./,ef_m/1.9/,flc_Ly/5./					(verdoppelung)
csl      data ef_s/0.5/,ef_m/0.95/,flc_Ly/5./					(original quellcode)

      ef_s   = 0.5  *ef_lange
      ef_m   = 0.95 *ef_lange
      flc_Ly = 5.0  *ef_lange
c
csl---------------------------------------------------------------------------------
c
c---- data for Hz band heating
      data sO21/5.5e-24/,sO22/1.34e-26/,sO31/6.22e-18/,sO32/1.66e-18/
      data fHz1/913./,fHz2/631./
c---- data for Ha band heating
      data sHa1/9.65e-18/,sHa2/4.79e-18/,fHa1/3680./,fHa2/1000./
c---- data for Hu band heating
      data sHu1/8.94e-19/,sHu2/1.09e-19/,sHu3/3.07e-18/
      data fHu1/6110.   /,fHu2/13250.  /,fHu3/2280.   /
c---- data for Ch band heating
      data sCh1/3.16e-21/,fCh1/364900./
      data cli   /16.,9.,-135.,-292.,202.,127.,337.,-144.,-137./
      data dnueli/1375.,1525.,1400.,1000.,1500.,1100.,1000.,540.,320./
      data x0li  /0.3047,0.3548,3.86,7.02,0.36,0.28,0.04,3.25,61.66/
      data dli   /31.,20.,230.,345.,460.,232.,246.,295.,77./
C------------------------------------------------ K_i/D_i, Liou (1992)
C---------------- pressure factors
      data dkli  /0.226,0.25,0.54,0.52,0.43,0.62,0.62,0.51,0.87/
      data fracli/0.0694,0.0801,0.1346,0.0892,0.1021,0.0622,0.03,
     $                                                0.0218,0.03/
c--- The heating efficiency is fixed for this version at efEUV = 0.366:
      data efEUV/0.366/
C.....DEFINITIONEN DER FORMEL - FUNKTIONEN:
c       function sh2o(a,b,ib,xn)
c      aa,bb   zoben,zunten [km]
c      sh2o  saueleninhalt wasserdampf [g/cm**2]
c      xn    druckkorrektur
c       skh  = 7.
c       skhn = skh/(1.+xn)
c       gphi = abs(- 87.5 + float(ib-1)*5.)
c       a1   = exp(-gphi*gphi/1700.)
c       dich0=1.273e-3
c       e1   = -1./skhn - 1./2.5
c       a1   = a1*0.02*dich0*1.e5/e1
c       a2   = 3.e-6*dich0*(-skhn)*1.e5
c       sh2o = a1* (exp(aa*e1) - exp(bb*e1)) +
c     *       a2* (exp(-aa/skhn) - exp(-bb/skhn))
c       return
c       end
c
c      sh2o1 (afkt1,bfkt1,ajfkt1,e1fkt1) =
c     *     ajfkt1/e1fkt1*(exp(afkt1*e1fkt1) - exp(bfkt1*e1fkt1))
c      sh2o2 (afkt1,bfkt1,hfkt1) =
c     *  3.819e-4*(-hfkt1)*(exp(-afkt1/hfkt1) - exp(-bfkt1/hfkt1))
c---------------- what is this????????????????????????????????????????
c     *  2.542e-4*(-hfkt1)*(exp(-afkt1/hfkt1) - exp(-bfkt1/hfkt1))
c---------------------------------------------------------------------
      tau_h=float(nsec)/10./86400.
      tsec=1.-exp(-tau_h*tau_h)   ! establishment of the solar heating
c------------------------------------------------ just for information,
c                                 these values have to be done in iNIT
c         co2pmv=360.*CO2
c         dich0=1.273e-3
c         h=7000.
c         skh=7.
c         RgSGS=8.31441e7
c---------------------------------------------------------------------          
      co2pmm=co2pmv*44./28.9644
      t0=1.013e6*rm(1)/RgSGS/dich0
      rgit=float(igit)
      rb=float(nb)
      do k=1,kgit
         tz(k)=0.
         do j=1,nb
            do i=1,igit     
              tz(k)=tz(k)+tn1(j,k,i)/rgit/rb
            end do
         end do
      Hsclkm=RgSGS/rm(k)*tz(k)*1.e-7/9.81
         do j=1,nb
c--- numerical densities of O3, O2, O - are INdependent on latitude
            xnumo3(j,k)=1.013  *rou(k)*ozppmv(j,k)/1.38066e-16/tz(k)
            xnumo2(j,k)=1.013e6*rou(k)*ro2vmr(  k)/1.38066e-16/tz(k)
            xnumox(j,k)=1.013e6*rou(k)*roxvmr(  k)/1.38066e-16/tz(k)
            xnumn2(j,k)=1.013e6*rou(k)*rn2vmr(  k)/1.38066e-16/tz(k)
c-------------------------------------------------------------------
            sox(j,k)=Hsclkm*rm(k)/16.*xnumox(j,k)*1.e5
            so2(j,k)=Hsclkm*rm(k)/32.*xnumo2(j,k)*1.e5
            sn2(j,k)=Hsclkm*rm(k)/28.*xnumn2(j,k)*1.e5
          end do
c---------------------------------- density and related profiles
       dichte(k)=dich0*rou(k)*rm(k)/rm(1)*t0/tz(k)
c
c----------------- Cp(z) - approx. 1.e7 (SGS) in lower atmosphere
c
       cp_z=RgSGS/rm(k)*(7./2.*(1.-roxvmr(k))+5./2.*roxvmr(k))
       dichcp(k)=dichte(k)*cp_z
      end do
c----------------------------------------------------------------
cfl      if(ncom.eq.0) print*,'heating rates: strobel_60.f'
         kf_h2o=int((22000.-zunten)/dz) +1
         kf_co2=int((55000.-zunten)/dz) +1
         do 15 ib=1,8
            fracd (ib) = fracli(ib) * dli(ib) * 0.4343
            do k=1,kf_h2o
              phgdk(ib,k) = (760.*rou(k)) ** dkli(ib)
            end do
 15      continue
            fracd (9) = fracli(9) * dli(9) * 0.4343
            do k=1,kf_co2
              phgdk(9,k) = (760.*rou(k)) ** dkli(9)
            end do
c     -----------------------------------------------------------------
      nday=ncom-nsun           !zeitschritt mit bewegter sonne
      dayl=float(ntime*24)     !zeitschritte/tag
      if(nday.lt.0) nday=0 
c ..............................DEKLINATION
c
c-------------------------------deklination in jul.tagen nach 
c                               31.dezember (12:00 Greenwich)
c                               ndek=1 -->1.jan 12:00
      day=float(ndek) !+float(nday)/dayl
c     day=80 --> fruehlings-equinoktium --> dek=0 (12:00 Greenwich)
      dek=sin(2.*pi*(day-82.)/365.)*23.5*pi/180.
      day64=dayl/float(igit)
c      day2=dayl/2.
      rncom=float(ncom)
      zeit=amod(rncom,dayl)
c
      do 160 j=1,nb
      do 150 i=1,igit
      iflag=1
      kflag=1
      ortzeit=zeit+float(i-1)*day64
      ortzeit=amod(ortzeit,dayl)
      stuwi=2.*pi*ortzeit/dayl
c------------------------------------------------------------------------
C     coschi= - cos(phi(j))*cos(dek)*cos(stuwi)+sin(phi(j))*sin(dek)
C                              ! coschi=-1 -->chi=180 ncom=0  -->24:00 UT
c------------------------------------------------------------------------
      coschi= + cos(phi(j))*cos(dek)*cos(stuwi)+sin(phi(j))*sin(dek)
c------------------------------------------------------------------------
C                              ! coschi=+1 -->chi=0   ncom=0  -->12:00 UT
c------------------------------------------------------------------------
      if(abs(coschi).gt.1.) coschi=sign(1.,coschi)
      chi=acos(coschi)*180./pi
      if(chi.le.75.) then
         do k=1,kgit
            secchi(k)  =1./coschi
            sech2o(k,j)=1./coschi
         end do
      else if (chi.gt.75.and.chi.lt.90.)  then
        do 131 k=1,kgit
            x=(6373.+z(k))/skh
            y=sqrt(.5*x)*abs(cos(chi*pi/180.))
            erfy=(aa+bb*y)/(ccc+dd*y+y*y)
            if (y.gt.8.)     erfy=f/(gs+y)
            secchi(k)=sqrt(pi/2.*x)*erfy
c            secchi(k)=atm_chapman(x,chi)
c------------------------------------------- H2O
            x=(6373.+z(k))*(skh+H_wv(j))/skh/H_wv(j)
            y=sqrt(.5*x)*abs(cos(chi*pi/180.))
            erfy=(aa+bb*y)/(ccc+dd*y+y*y)
            if (y.gt.8.)     erfy=f/(gs+y)
            sech2o(k,j)=sqrt(pi/2.*x)*erfy
c            sech2o(k,j)=atm_chapman(x,chi)
 131    continue
      else if (chi.gt.sunris(kgit))         then
        iflag=0
      else
        do 133 k=1,kgit
            if (chi.gt.sunris(k))         then
               secchi(k)  =0.
               sech2o(k,j)=0.
               kflag=k
            else
               x=(6373.+z(k))/skh
               y=sqrt(.5*x)*abs(cos(chi*pi/180.))
               erfy=(aa+bb*y)/(ccc+dd*y+y*y)
               if (y.gt.8.)     erfy=f/(gs+y)
               arg=x*(1.-sin (chi*pi/180.))
               secchi(k)=sqrt(2.*pi*x)*(sqrt(sin(chi*pi/180.))
     *              *exp(arg)-0.5*erfy)
c               secchi(k)=atm_chapman(x,chi)
c------------------------------------------- H2O
               x=(6373.+z(k))*(skh+H_wv(j))/skh/H_wv(j)
c------------- to avoid owerflow for large X
               if(x.gt.3000.) x=3000.
               y=sqrt(.5*x)*abs(cos(chi*pi/180.))
               erfy=(aa+bb*y)/(ccc+dd*y+y*y)
               if (y.gt.8.)     erfy=f/(gs+y)
               arg=x*(1.-sin (chi*pi/180.))
               sech2o(k,j)=sqrt(2.*pi*x)*(sqrt(sin(chi*pi/180.))
     *              *exp(arg)-0.5*erfy)
c            sech2o(k,j)=atm_chapman(x,chi)
            endif
 133    continue
         kflag=kflag+1
      endif
       if (iflag.eq.0)                   then
         do k=1,kgit
            fsolar(k)=0.
         end do
      else   
c         mark6 =int(( 16000.-zunten)/dz) + 1
         mark6 =int(( 6000.-zunten)/dz) + 1
         mark12=int(( 32000.-zunten)/dz) + 1
         mark26=int(( 73000.-zunten)/dz) + 1
         mark32=int(( 90000.-zunten)/dz) + 1
         mark39=int((110000.-zunten)/dz) + 1
c
         k1flg  = max0 (1,kflag)
         k6flg  = max0 (mark6 ,kflag)
         k11flg = max0 (mark12,kflag)
         k26flg = max0 (mark26,kflag)
         k32flg = max0 (mark32,kflag)
c
         k26fix = mark26
         k32fix = mark32
         k39fix = mark39
c     48-level version:
c         k1flg  = max0 (1,kflag)
c         k6flg  = max0 (6,kflag)
c         k11flg = max0 (12,kflag)
c         k32flg = max0 (32,kflag)
c         k26fix = 26
c         k32fix = 32
c         k39fix = 39
         do k=1,kgit
            fsolar(k)=0.
            secnox(k)=sox(j,k)*secchi(k)
            secno2(k)=so2(j,k)*secchi(k)
            secnn2(k)=sn2(j,k)*secchi(k)
            secno3(k)=so3(j,k)*secchi(k)*273./tn1(j,k,i)
         end do
C**********************************************************************
C     COMMA_LIM version of H2O heating,  October 2002, Alex Pogoreltsev
C     Rayleygh scattering and surface reflection are taken into account
C     Liou K.N. "Radiation and Cloud Processes in the Atmosphere",
C     Oxford Univ. Press, Ney York, Oxford, 1992, 487 pp.
C     Lambertian albedo =0.25, Strugnell et al., GRL, 2001, 28, 191-194 
C**********************************************************************
C--- the upper boundary for H2O heating calculation
      zh2oa = 200.
C--- The diffusivity factor [1/sec(mu_av), Liou, 1992]
      sec_ef= 1.66
C--- reflection by the Rayleigh layer (=0 for H2O, Liou, 1992)
      r_a   = 0.
C--- The global reflection due to Rayleigh scattering
      ra_ef = 0.144
C--- Lambertian albedo (globally averaged)
      r_s   = 0.25
C--- The combine reflection due to Rayleigh leyer and the surface
      rmu_0 = r_a + (1.-r_a)*(1.-ra_ef)*r_s/(1.-ra_ef*r_s)
C-----------------------------------------------------------------------
      DO ib=1,8 
C--- H2O absorbing path length for atmosphere (z=0 -- > z=200 km) 
      uliou0(ib)=760.**dkli(ib)* 
     $       (sech2o(1,j)*sh2o1(zh2oa,0.,a1h2o(i,j),expnuv(ib,j))+
     $        secchi(1  )*sh2o2(zh2oa,0.,skhnuv(ib)))
C--- H2O absorptance at the ground
      aliou0(ib)=(cli(ib)+
     $            dli(ib)*ALOG10(uliou0(ib)+x0li(ib)))/dnueli(ib)
      END DO
      do 142 k=k1flg,kf_h2o
      abliou = 0.
      DO ib=1,8
C*********************************************************************
C integral(Z -- > Z_inf) [rho_H2O*p**(D/K)dZ']/COS(mu_0), dkli(ib)=D/K 
C "uliou(ib)"  has to be in (g/cm^2)
C*********************************************************************
C--- H2O absorbing path length (z   -- > 200 km)
      uliou (ib)=760.**dkli(ib)*(T_1000(j,i)/tn1(j,k,i))*
     $                          (273.     /tn1(j,k,i))**0.45*
     $       (sech2o(k,j)*sh2o1(zh2oa,z(k),a1h2o(i,j),expnuv(ib,j))+
     $        secchi(k  )*sh2o2(zh2oa,z(k),skhnuv(ib)))
C--- H2O absorbing path length (z=0 -- > z)
      ulioub(ib)=sec_ef   *760.**dkli(ib)*
     $           (sh2o1(z(k) ,0.  ,a1h2o(i,j),expnuv(ib,j))+
     $            sh2o2(z(k) ,0.  ,skhnuv(ib)))
C*********************************************************************
C H2O density (g/cm^3),(dU/dz= - rho_H2O), D_1000 = density at 1000 mb 
C********************************************************************* 
      D_1000=1.013e6*rm(1)/RgSGS/T_1000(j,i)
      dudz(ib)=-dichte(k)*
     $          (a1h2o(i,j)/D_1000/1.e5*exp(-z(k)/H_wv(j))+3.e-6)*
     $                        phgdk(ib,k)*(273./tn1(j,k,i))**0.45
C--- nondimansional
      abliou = abliou + fracd(ib)*dudz(ib)/dnueli(ib)*
     $                                  (1./(uliou (ib)+x0li(ib))+
     $                     sec_ef*rmu_0/secchi(k)*(1.-aliou0(ib))*
     $                                   1./(ulioub(ib)+x0li(ib)))
      END DO 
C**********************************************************************
C           f0sun(J/m^2/s)*1.e7/1.e4 -- > (erg/cm/s)
C**********************************************************************
      fsolar(k)=-f0sun*1.e3*abliou/dichcp(k)
      go to 5555
C----------------------------------------------------------------------
C    new parameterization of H2O heating by Freidenreich and Ramaswamy, 
C    JGR, 1999, Vol. 104, No. D24, P. 31,389-31,409
C    Alex Pogoreltsev (October, 2003)
C----------------------------------------------------------------------
      ht_h2o=0.     
      DENS_0=1.e6*rm(k)/RgSGS/tn1(j,k,i)
      do 1111 kk=1,53
      e_h2o_1=EXP(-z(k)/H_h2o_1(kk,j)) 
      e_h2o_2=EXP(-z(k)/H_h2o_2(kk)) 
      w1_h2o = s1_h2o(kk)*DENS_0*1.e5*
     $                            c1_h2o(i,j)*H_h2o_1(kk,j)*e_h2o_1
      w2_h2o = s1_h2o(kk)*DENS_0*1.e5*
     $                            c2_h2o     *H_h2o_2(kk)  *e_h2o_2
c---------------------------------------------------------------------
      s2_h2o(kk)=F_h2o(kk)/1357.96*a_h2o(kk)*b_h2o(kk)*
     $                                 s1_h2o(kk)*1.273e-3*t0/tz(k)
      ht_h2o = ht_h2o+
     $                  s2_h2o(kk)*(c1_h2o(i,j)*e_h2o_1+
     $                              c2_h2o     *e_h2o_2)*
     $        EXP(-b_h2o(kk)*(sech2o(k,j)*w1_h2o+secchi(k)*w2_h2o))
 1111 continue
c      fsolar(k)= f0sun*1.e3*ht_h2o/dichcp(k)
5555  continue
  142 continue 
C--------------------------------------------------- END of H2O heating
C**********************************************************************
C     COMMA_LIM version of CO2 heating,  October 2002, Alex Pogoreltsev
C     Liou K.N. "Radiation and Cloud Processes in the Atmosphere",
C     Oxford Univ. Press, Ney York, Oxford, 1992, 487 pp.
C**********************************************************************
      do 1142 k=k1flg,kf_co2
C----------------------------------------------------------------------
C The absorbing path (cm-atm for CO2)
C     satmcm =skh*1.e5*6.0221e23*co2pmm*1.e-6*dichte(k)/2.68681.e19/44.
C     satmcm = H(cm)*N_A*co2pmm*1.e-6*rho(z)/n_0/44,
C     where N_A is Avogadro number, n_0 is Loschmidt's number
C---------------------------------------------------------------------- 
c      satmcm =skh*0.051e3*co2pmm*dichte(k)
c----------------------------------------------------------------------
c     skh(cm)*co2pmv*1.e-6*rou(k) -- > skh(km)*1.e5*co2pmv*1.e-6*rou(k)
c----------------------------------------------------------------------
      satmcm =skh*0.1*co2pmv*rou(k)
C--- cm-atm * 760.^dkli(9)
c      uco2 = secchi(k)*satmcm*phgdk(9,k)/(1.+dkli(9))
C--- cm-atm/cm * 760.^dkli(9)
c      dudzco = - satmcm*phgdk(9,k)/skh/1.e5
c---------------------------------------------- T dependence 
      uco2 = secchi(k)*satmcm*phgdk(9,k)*(273./tn1(j,k,i))**10./
     $                                                     (1.+dkli(9))
      dudzco = - satmcm*phgdk(9,k)*(273./tn1(j,k,i))**10./skh/1.e5
C--- nondimantional
      abliou=fracd(9)*dudzco/dnueli(9)/(uco2+x0li(9)) 
C**********************************************************************
C           f0sun(J/m^2/s)*1.e7/1.e4 -- > (erg/cm/s)
C**********************************************************************
      fsolar(k)  = fsolar(k) - f0sun*1.e3*abliou/dichcp(k)
C--------------------------------------------------- END of CO2 heating
1142  continue
c     ****************************************************************
c     **                    chappius                                **
c     ****************************************************************
      do 143 k=k6flg,k39fix  !!!=109.4km!!!
c---------------------------------------------------------------------
c     account of the diffuse scattered solar radiation, Strobel (1978)
c---------------------------------------------------------------------
c      ef_ch=0.5 
csl    effizienz geändert für solarmax 1. -> 1. *ef_lange
csl      ef_ch=1. 
Csl      ef_ch=1.*ef_lange

      ef_ch=0.5 *ef_lange
        expqc  =-sCh1*secno3(k)
         expqc0 =-sCh1*secno3(1)
         E2     =exp(-1.9*sCh1*(so3(j,1)-so3(j,k)))         
         fsolar(k)=ef_ch*xnumo3(j,k)*fCh1*sCh1*
     $          (exp(expqc)+2.*0.25*exp(expqc0)*E2/secchi(k))/dichcp(k)
     $         + fsolar(k)
c          fsolar(k)=ef_ch*xnumo3(j,k)*fCh1*sCh1*exp(expqc)/dichcp(k)
c     $         + fsolar(k)
 143  continue
c     ****************************************************************
c     **                    Ly alpha O2                             **
c     ****************************************************************
      do 243 k=k26flg,kgit
      if(secno2(k).le.1.e19) then
       fac_Ly= 1.e-20
       expqc =-1.e-20*secno2(k)
      else
       fac_Ly= 4.17e-19/(secno2(k)**0.083)
       expqc =-4.17e-19*(secno2(k)**0.917)
      end if
      fsolar(k)=ef_s*ef_m*xnumo2(j,k)*flc_Ly*fac_Ly*exp(expqc)/dichcp(k)
     $         + fsolar(k)
 243  continue
c     ****************************************************************
c     **                    huggins                                 **
c     ****************************************************************
c--- efficiencies according Strobel (1981)      
c      ef_hu = 0.75
csl   effizienz geändert für solarmax 1. -> 1.*ef_lange
csl      ef_hu = 1.
      ef_hu = 0.75*ef_lange
      do 144 k=k6flg,k39fix  !!!=109.4km, 26=72.5km, 31=86.7!!!
         exphu1= -sHu1*secno3(k) 
         exphu2= -sHu2*secno3(k) 
         exphu3= -sHu3*secno3(k) 
         fsolar(k)=xnumo3(j,k)*ef_hu*
     $    (sHu1*fHu1*exp(exphu1)+sHu2*fHu2*exp(exphu2)+
     $     sHu3*fHu3*exp(exphu3))/dichcp(k)+fsolar(k)
 144  continue
c     ****************************************************************
c     **                    herzberg                                **
c     ****************************************************************
c      ef_O2 = 0.055
c      ef_O3 = 0.806 
c
csl   effizienz geändert für solarmax 1. -> 1.*ef_lange
csl      ef_O2 = 1.
csl      ef_O3 = 1.

   
      ef_O2 = 0.055*ef_lange
      ef_O3 = 0.806*ef_lange     
      do 145 k=k6flg,k39fix
         exphz1 = -sO21*secno2(k)-sO31*secno3(k)
         exphz2 = -sO22*secno2(k)-sO32*secno3(k)
         fsolar(k)=(fHz1*(ef_O2*sO21*xnumo2(j,k)+
     $                    ef_O3*sO31*xnumo3(j,k))*exp(exphz1)+
     $              fHz2*(ef_O2*sO22*xnumo2(j,k)+
     $                    ef_O3*sO32*xnumo3(j,k))*exp(exphz2))
     $              /dichcp(k) + fsolar(k)
 145  continue
c     ****************************************************************
c     **                    hartley                                 **
c     ****************************************************************
c--- efficiencies according Strobel (1981), Mlynczak and Solomon (1993)

csl   effizienz geändert für solarmax 1. -> 1.*ef_lange
c        ef_ha = 0.786         
Csl      ef_ha = 1.*ef_lange

      ef_ha = 0.786*ef_lange
      e_ms=1.
      do 146 k=k6flg,k39fix
c----------------------------------------------------------------------
c      The efficiency is calculated using the expression proposed by
c      Mlynczak and Solomon, JGR, 1993, v. 98, No. D6, p. 10517-10541. 
c----------------------------------------------------------------------
       pressz = 1013.*EXP(-z(k)/skh)
       IF(pressz.gt.1.) GO TO 10
        if(pressz.le.0.01) then
         xpr = ALOG10(pressz)+3.
         xpr2= xpr*xpr
         e_ms = 0.669650-0.009682*xpr+0.033093*xpr2+0.017938*xpr*xpr2
        else
         xpr = ALOG10(pressz)+1.
         xpr2= xpr*xpr
         e_ms = 0.926210+0.133960*xpr-0.076863*xpr2+0.006897*xpr*xpr2
        end if
   10 CONTINUE
         expha1 = -sHa1*secno3(k)
         expha2 = -sHa2*secno3(k)
c---------------------------- looses of the energy only due to airglow
         fsolar(k)=(1.-ef_ha*(1.-e_ms))*xnumo3(j,k)*
     $             (fHa1*sHa1*exp(expha1)+fHa2*sHa2*exp(expha2))
     $             /dichcp(k) + fsolar(k)
  146 continue
c     ****************************************************************
c     **              schumann runge bande                          **
c     ****************************************************************
      do 147 k=k26flg,kgit
         if(secno2(k).gt.1.e18)      then
            fsolar(k)=xnumo2(j,k)/(0.67*secno2(k)
     *           +3.44e9*sqrt(secno2(k)))/dichcp(k) + fsolar(k)
         else
            fsolar(k)=xnumo2(j,k)*2.43e-19/dichcp(k) + fsolar(k)
         endif
 147  continue
c     ****************************************************************
c     **      schumann runge kont. (1250 bis 1520 ang)              **
c     ****************************************************************
      do 148 k=k32flg,kgit
         exprs=-abscr*secno2(k)
         exa=exp(exprs)
         qscr1=escr*fscr*abscr*exa*xnumo2(j,k)/dichcp(k)
c     ****************************************************************
c     **      schumann runge kont.  (1520 bis 1750 ang)             **
c     ****************************************************************
         exprs1=-abl*secno2(k)
         exprs2=-abm*secno2(k)
         exprs3=-abk*secno2(k)
         ex1=exp(exprs1)
         ex2=exp(exprs2)
         ex3=exp(exprs3)
c----------------------------------------------------------------------
c      The efficiency is calculated using the expression proposed by
c      Mlynczak and Solomon, JGR, 1993, v. 98, No. D6, p. 10517-10541. 
c----------------------------------------------------------------------
        pressz = 1013.*EXP(-z(k)/skh)
        if(pressz.ge.0.0001) then
         xpr = ALOG10(pressz)+3.
        else
         xpr = ALOG10(0.0001)+3.
        end if
         xpr2= xpr*xpr
         xeps= 0.75349+0.0036*xpr+0.059468*xpr2-0.022795*xpr*xpr2
        if(xeps.gt.1.) xeps=1.
         fsolar(k)=xeps*xnumo2(j,k)*(eilm*ex1+(eism-eilm)*ex2-eism*ex3)
     &                  /dichcp(k)/secno2(k)+qscr1 + fsolar(k)
 148  continue
C***********************************************************************
C                               EUV heating
c to calculate the solar heating due to absorption in Extreme UV region
c (5-105nm). Data on fluxes and absorption coefficients are taken from
c Richards et al., 1994, JGR, vol 99, No A7, 13,283
c The heating efficiency is fixed for this version at efEUV = 0.366:
c Fluxes (erg/cm2/s) and absorption coefficients (cm2) in 37 intervals
c (received from  Richards et al)
c
c ---- to determine the flux and absorption cross-sections, if first entry:
c
csl ---- F107 und F107A max zwischen 200 und 250 (im jahr 2002) mittel 120
csl ---- F107 und F107A max zwischen 200 und 250 (im jahr 2002) mittel 120, min zwischen 80 und 60 (im Jahr 2009)**http://www.swpc.noaa.gov/SolarCycle/
CPH      F107  = 240.
CPH      F107A = 240.
c
      if(IEUV.eq.0) then
      IEUV = 1
      F107A = F107
c
c
      call EUVflcs(F107,F107A)
      end if
C-------------------------------------------------- EUV heating rate (K/s):
      do 1148 k=k32flg,kgit
c---------------------------------------- loop  over the spectral intervals
      do it=1,37
      eEUV = sEUVo2(it)*secno2(k)+sEUVn2(it)*secnn2(k)+
     $                            sEUVo (it)*secnox(k)
      sEUV= EUVflx(it)*(sEUVo2(it)*xnumo2(j,k)+
     +                  sEUVn2(it)*xnumn2(j,k)+sEUVo(it)*xnumox(j,k))
      fsolar(k)=efEUV*sEUV*EXP(-eEUV)/dichcp(k) + fsolar(k)
      end do
 1148 continue
      end if
c--------------------------------------------------------- chemical heating
      do 149 k=26,kgit
      cp_z=RgSGS/rm(k)*(7./2.*(1.-roxvmr(k))+5./2.*roxvmr(k))
      
      !included by F. Lilienthal (overflow without double: 1.6022d-12*9.59d-34=0.0e0 )
      f_chem=(5.12*1.6022d-12*9.59d-34*DEXP(480.d0/dble(tn1(j,k,i)))
     $	     *dble(xnumox(j,k))**2
     $	     +1.d0*1.05d0*1.6022d-12*6.0d-34*
     $	     (300.d0/dble(tn1(j,k,i)))**2.4d0*dble(xnumox(j,k))**2)
     $	     *6.0221d23/dble(rm(k))/dble(cp_z)  
c        print*,f_chem !check
        fsolar(k)=fsolar(k)+real(f_chem,4)
!----------------earlier solution with real:----------------------------		
        !       fsolar(k)=fsolar(k)+
!      $      (5.12*1.6022e-12*
! c------------------ Riese et al., JGR, 1994, 99, D7, 14,585-14,593. 
! c     $     4.7e-33*(300./tn1(j,k,i))**2  *xnumox(j,k)*xnumox(j,k)+
! c------------------ Roble in The Upper Mesosphere an Lower Thermosphere
! c                   Geophys. Monograph 87, pp. 1-21.
!      $     9.59e-34*EXP(480./tn1(j,k,i)) *xnumox(j,k)*xnumox(j,k)+
! c------------------------- efficiency = 0.0 OR 1.0
!      $     1.0*1.05*1.6022e-12*
!      $     6.0e-34*(300./tn1(j,k,i))**2.4*xnumox(j,k)*xnumo2(j,k))*
!      $                      6.0221e23/rm(k)/cp_z  
!-----------------------------------------------------------------------       
 149  continue
      do k=1,kgit
        hs(k,i,j)=fsolar(k)*tsec
      end do
 150  continue  !! -i-schleife, longitude
 160  continue  !! -j-schleife, latitude
c---------------------------------------------------- tidal forcing     
      rgit  = float(igit)
c      nfor=max0(ncom-nphi,0)
c----------------------------------------------------------------------
c------------------------------- without diurnal variability of heating
c                                until ncom > 207360 
c                               (240 days with 2-hourly output, 100sec)  
c----------------------------------------------------------------------
      nfor=max0(ncom-nphi,0) 
      xsec   = 3600./float(ntime)          ! seconds per step
      tdfor  = float(nfor)*xsec/10./86400. 
c--------------------- efficiency of the generation of tides after nphi
c                      NOW after 84489 steps!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      t_tides= 1.- exp(-tdfor*tdfor) 
      do j=1,nb
        do k=1,kgit
        h_ave=0.     
          do i=1,igit
            h_ave=h_ave+hs(k,i,j)/rgit
          end do
            do i=1,igit
            hs_tide(k,i,j)=hs(k,i,j)-h_ave
            hs(k,i,j)=h_ave+t_tides*hs_tide(k,i,j)
            end do
        end do
      end do
      print *,'hs:',hs(15,1,26),hs(32,1,26) !FL hs(k,i,j) k=15 44km,k=32 92km, i=0 0deg, j=26 42.5degN
      return
      end
c------------------------------------------------------------------------
      FUNCTION sh2o1 (afkt1,bfkt1,ajfkt1,e1fkt1)
      sh2o1 =
     *     ajfkt1/e1fkt1*(exp(afkt1*e1fkt1) - exp(bfkt1*e1fkt1))
      RETURN
      END
c-------------------------------------------------------------
      FUNCTION sh2o2 (afkt1,bfkt1,hfkt1)
      sh2o2 =
     *  3.819e-4*(-hfkt1)*(exp(-afkt1/hfkt1) - exp(-bfkt1/hfkt1))
      RETURN
      END
c -----------------------------------------------------------------------
      BLOCK DATA EUVREP
      common /EUVI/ IEUV
      data IEUV/0/
      end

C***********************************************************************
C*                                                                     *
C*              SUBROUTINE EUVheat                                     *
C*                                                                     *
C***********************************************************************
c
c to calculate the solar flux (erg/cm2/s) and absorption cross-sections (cm2)
c in 37 intervals, using Richards et al. (1994, JGR, vol 99, No A7, 13283)
c data
c                                           December 1997, V. Fomichev

      SUBROUTINE EUVflcs(F107,F107A)

      common /EUVFCS/ EUVflx(37),sEUVo2(37),sEUVn2(37),sEUVo(37)

      real al(37), so2o(37),sn2o(37),soo(37),
     :     F74113(37), Ai(37), photon(37)

c centers of the bins (nm):

      data (al(i), i=1,37)/
     & 7.5,   12.5,  17.5,  22.5,   25.63, 28.415, 27.5, 30.331,
     & 30.378,32.5,  36.807,37.5,   42.5,  46.522, 47.5, 52.5,
     & 55.437,58.433,57.5,  60.976, 62.973,62.5,   67.5, 70.336,
     & 72.5,  76.515,77.041,78.936, 77.5,  82.5,   87.5, 92.5,
     & 97.702,97.5, 102.572,103.191,102.5/

c cross sections (in 1.e18 cm2)

      data (so2o(i), i=1,37)/
     &  1.316, 3.806, 7.509,10.900,13.370,15.790,14.387,16.800,
     & 16.810,17.438,18.320,18.118,20.310,21.910,23.101,24.606,
     & 26.040,22.720,26.610,28.070,32.060,26.017,21.919,27.440,
     & 28.535,20.800,18.910,26.668,22.145,16.631, 8.562,12.817,
     & 18.730,21.108, 1.630, 1.050, 1.346/ 
      data (sn2o(i), i=1,37)/
     &  0.720, 2.261, 4.958, 8.392,10.210,10.900,10.493,11.670,
     & 11.700,13.857,16.910,16.395,21.675,23.160,23.471,24.501,
     & 24.130,22.400,22.787,22.790,23.370,23.339,31.755,26.540,
     & 24.662,120.49,14.180,16.487,33.578,16.992,20.249, 9.680,
     &  2.240,50.988, 0.0,   0.0,  0.0/
      data (soo(i), i=1,37)/
     &  0.730, 1.839, 3.732, 5.202, 6.050, 7.080, 6.461, 7.680,
     &  7.700, 8.693, 9.840, 9.687,11.496,11.930,12.127,12.059,
     & 12.590,13.090,13.024,13.400,13.400,13.365,17.245,11.460,
     & 10.736, 4.000, 3.890, 3.749, 5.091, 3.498, 4.554, 1.315,
     &  0.0,   0.0,   0.0,   0.0,   0.0/

c parameters for the solar EUV flux model (flux in 1.e-9 photon/cm2/s)

      data (F74113(i), i=1,37)/
     & 1.200, 0.450, 4.800, 3.100, 0.460, 0.210, 1.679, 0.800,
     & 6.900, 0.965, 0.650, 0.314, 0.383, 0.290, 0.285, 0.452,
     & 0.720, 1.270, 0.357, 0.530, 1.590, 0.342, 0.230, 0.360,
     & 0.141, 0.170, 0.260, 0.702, 0.758, 1.625, 3.537, 3.000,
     & 4.400, 1.475, 3.500, 2.100, 2.467/
      data (Ai(i), i=1,37)/
     & 1.0017e-02, 7.1250e-03, 1.3375e-02, 1.9450e-02, 2.7750e-03,
     & 1.3768e-01, 2.6467e-02, 2.5000e-02, 3.3333e-03, 2.2450e-02,
     & 6.5917e-03, 3.6542e-02, 7.4083e-03, 7.4917e-03, 2.0225e-02,
     & 8.7583e-03, 3.2667e-03, 5.1583e-03, 3.6583e-03, 1.6175e-02,
     & 3.3250e-03, 1.1800e-02, 4.2667e-03, 3.0417e-03, 4.7500e-03,
     & 3.8500e-03, 1.2808e-02, 3.2750e-03, 4.7667e-03, 4.8167e-03,
     & 5.6750e-03, 4.9833e-03, 3.9417e-03, 4.4167e-03, 5.1833e-03,
     & 5.2833e-03, 4.3750e-03/

c --- to determine cross sections:

      do 1 i = 1,37
      sEUVo2(i) = so2o(i)*1.e-18
      sEUVn2(i) = sn2o(i)*1.e-18
      sEUVo (i) = soo(i)*1.e-18
    1 continue

c --- to calculate the flux for given F107 and F107A:

      P = (F107 + F107A)/2.

c --- fluxes (in photon/cm2/s) for given F107 and F107A:

      do 2 i = 1,37
      photon(i) = F74113(i) * (1 + Ai(i)*(P-80)) *1.e9
    2 continue

c --- fluxes in erg/cm2/s, flux = h*c/al * photon (1nm = 1.e-7 cm):

      do 3 i = 1,37
      EUVflx(i) = 1.98648e-09/al(i) * photon(i)
    3 continue

      return
      end











