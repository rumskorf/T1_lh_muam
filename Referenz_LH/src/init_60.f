*********************************************************************
c     Pogoreltsev 2006
C     Rf_c=0.4
C     SPW1 and SPW2 inserted from the UKMO or NCEP/NCAR data
C
C --- INPUT: nb, igit, zoben, rm(kgit), rhos (kg/m^3), skh
C--------------------------------------------------------------------
C --- OUTPUT: 
C         dz,defi,zoben,zunten,nstep,dt,dt2,cor(nb+1),tgfia(nb+1),...
C         dx(nb+1),phi(nb+1)  
C      rou(kgit),row(kgit),alfa(kgit),drag(nb,kgit,1),drag(nb,kgit,2) 
C         dichte(kgit),
C         a1h2o(igit,nb),
C         skhnuv(9),skhnir(4),e1nuv(9),e1nir(4) --> strobel.f ircool.f
C-------------------------------------------------------------------
        subroutine init_60
! c
! c--------------------------------------    INITIALIZATION FOR COMMA_64 
! c
        include 'com_main.fc'
        include 'com_namelist.fc'
        include 'data_TMP_SPW.f'
        include 'data_GH_SPW.f'
        include 'data_H2O.f'
! c
! c     --------------------------------------------------------------
! c     ntime        anzahl der zeitschritte pro stunde (dt=3600/ntime
! c                  ntime ist so zu  waehlen, dass dt ganz zahl)
! c     ncom         aktueller zeitschritt
! c     nsec         aktueller zeitpunkt
! c     nend         letzter prognostizierter zeitschritt
! c     nout         anzahl der zeitschritte bei denen jeweils auslagerung
! c                  von zwischenergebnissen auf band oder platte
! c     nprint       anzahl der zeitschritte  bei denen jeweils  ausdruck
! c                  von zwischenergebnissen erfolgt
! c     igit         anzahl der gitterpunkte laengs eines breitenkreises
! c     kgit         anzahl der modellflaechen
! c     nwelle(1-3)  gibt an, welche wellen im forcing enthalten sind
! c     nt0          t0 in forcing function
! c     nfase(1-3)   phase der wellen 1-3
! c     jforc1       breitenkreis bei dem forcing beginnt
! c     jforc2       breitenkreis bei dem forcing endet
! c     nfinul(1-3)  max. amplitude des forcing
! c     rnn        : vert. profile von n2, o!!!!!!!!!!!!!!!!!!!!!!!!!, o2, no
! c     rm         : vert. profil des molekulargewichts
! c     rl (kgit+1)  molekulare waermeleitfaehigkeit
! c     rx (kgit+1)  dynamische viskositaet
! c     xmt (kgit+1)   xmt=m'/tau'
! c     rou (kgit)     z-koordinate im niveau k
! c     row (kgit)     z-koordinate im niveau k-0.5
! c----------------------------------------------------------------------
      REAL T0_NCEP(36),G0_NCEP(36) 
      REAL gh1a(36),gh1f(36),gh2a(36),gh2f(36),gh3a(36),gh3f(36)
      REAL tp1a(36),tp1f(36),tp2a(36),tp2f(36),tp3a(36),tp3f(36)
      REAL dkli(9),xncir(4),rh2o(36)
      REAL alfaNC(60),sigma1(60),sigma2(60)
      REAL*4 ht_m1 (3,nb,10,4),ht_m2 (3,nb,10,4),
     &       ht_m3 (3,nb,10,4),ht_m4 (3,nb,10,4),
     &       ht_m0 (2,nb,10,2)
      REAL*4 ht_a(nb,10),ht_a0(nb,10)



!       INTEGER NPHI,NEND,NPRINT,NOUT,NTIME,KRET,NSUN,NDEK
!       CHARACTER*6 FPHI
!       CHARACTER*7 FVER,FGFU,FGFV,FGAU,FGAV
!       CHARACTER*12 FOLD
!       CHARACTER*14 FNEW,FO3
!       CHARACTER*17 FUVT
!       CHARACTER*11 TLB
!       CHARACTER*13 T00,T11,T22,T33,G00,G11,G22,G33
!       INTEGER MFORC,MTURB,MCOOL,MUVIR,MIOND,MMCON,MSURF
!       INTEGER I_2d,I_SK,I_KW,I_FK,I_NK,I_10,I_11,I_12,I_13,I_21,I_23
!       REAL Xm2d,XmSK,XmFK,XmKW,XmNK,Xm10,Xm11,Xm12,Xm13,Xm21,Xm23
!       REAL Per_2d,Per_SK,Per_FK,Per_KW,Per_NK,Per_10,Per_11,Per_12,
!      #     Per_13,Per_21,Per_23
!       REAL fi0_2d,fi0_SK,fi0_FK,fi0_KW,fi0_NK,fi0_10,fi0_11,fi0_12,
!      #     fi0_13,fi0_21,fi0_23
!       REAL XL,WAMPL,CPHASE(6),TETA(8)
!       REAL t0
! 
!       NAMELIST /RUN/ NPHI,NEND,NPRINT,NOUT,NTIME,KRET,NSUN,NDEK
!       NAMELIST /INI/ FOLD,FNEW,FO3
!       NAMELIST /FDX/ FUVT,FPHI,FVER,FGFU,FGFV,FGAU,FGAV
!       NAMELIST /LOB/ TLB,T00,T11,T22,T33,G00,G11,G22,G33
!       NAMELIST /FOR/ MFORC,MTURB,MCOOL,MUVIR,MIOND,MMCON,MSURF
!       NAMELIST /PWS/ I_2d,I_SK,I_KW,I_FK,I_NK,I_10,I_11,I_12,
!      #               I_13,I_21,I_23
!       NAMELIST /PWX/ Xm2d,XmSK,XmFK,XmKW,XmNK,Xm10,Xm11,Xm12,
!      #               Xm13,Xm21,Xm23
!       NAMELIST /PWP/ Per_2d,Per_SK,Per_FK,Per_KW,Per_NK,Per_10,
!      #		     Per_11,Per_12,Per_13,Per_21,Per_23
!       NAMELIST /PWF/ fi0_2d,fi0_SK,fi0_FK,fi0_KW,fi0_NK,fi0_10,
!      #		     fi0_11,fi0_12,fi0_13,fi0_21,fi0_23
!       NAMELIST /GRW/ XL,WAMPL,CPHASE,TETA
! !      NAMELIST /SOL/

cfl:
!       OPEN(UNIT=10,FILE='muam_mod.txt')
! 
!       READ(10,RUN)
!       READ(10,INI)
!       READ(10,FDX)
!       READ(10,LOB)
!       READ(10,FOR)
!       READ(10,PWS)
!       READ(10,PWX)
!       READ(10,PWP)
!       READ(10,PWF)
!       READ(10,GRW)
!	READ(10,SOL)
! 
!       CLOSE(10)

      print*,'check parameters'
      print*,'RUN ',NPHI,NEND,NPRINT,NOUT,NTIME,KRET,NSUN,NDEK,YEAR
      print*,'INI ',FOLD,' ',FNEW,' ',FO3
      print*,'FDX ',FUVT,' ',FPHI,' ',FVER,' ',FGFU,' ',FGFV,' ',FGAU,
     +' ',FGAV,' ',FGAT
      print*,'LOB ',TLB,' ',T00,' ',T11,' ',T22,' ',T33,' ',G00,' ',
     +G11,' ',G22,' ',G33
      print*,'FOR ',MFORC,MUVIR,MIOND,LODEF
      print*,'PWS ',I_2d,I_SK,I_KW,I_FK,I_NK,I_10,I_11,I_12,I_13,I_21,
     +	I_23
      print*,'PWX ',Xm2d,XmSK,XmFK,XmKW,XmNK,Xm10,Xm11,Xm12,Xm13,Xm21,
     +  Xm23
      print*,'PWP ',Per_2d,Per_SK,Per_FK,Per_KW,Per_NK,Per_10,Per_11,
     + Per_12,Per_13,Per_21,Per_23
      print*,'PWF ',fi0_2d,fi0_SK,fi0_FK,fi0_KW,fi0_NK,fi0_10,fi0_11,
     + fi0_12,fi0_13,fi0_21,fi0_23
      print*,'GRW ',XL,WAMPX,CPHASE,TETA
      print*,'SOL ',F107,SOLE,SOLC,OZON,CO2

      data a /6366197./
C-------------------------------------------------- K_i/D_i, Liou (1992)
C---------------- pressure factors
      data dkli  /0.226,0.25,0.54,0.52,0.43,0.62,0.62,0.51,0.87/
c------------------------------ coeff. For h2o column in ircool_newCO2.f
c      data xncir /0.5,0.5,0.5,0.5/ ! xnc
      data xncir /1.,1.,1.,1./ ! xnc
c--------------------------------------------------- PHYSICAL PARAMETERS:
c
c---------- daily averaged - /1.5
      data sigma1/  35*0.0, !60 levels coefficient for ion drag
     $0.02,0.024,0.08,0.29,0.48,1.22,4.58,11.2,23.0,38.2,47.9,54.8,
     $13*60.5/
c---------- daily averaged - /1.5
      data sigma2/  35*0.0, !60 level coefficient for Lorentz deflection
     $0.2,0.8,1.8,3.2,4.9,7.7,13.3,20.0,23.7,21.1,14.5,9.4,
     $13*6.0/
c-------------------- Newtonian cooling profile 
c  Zhu, (1993) JAS, v. 50, No 17, p. 3008-3021.    
      data alfaNC/3*0.070,0.08,0.110,0.124,0.176,0.277,0.390,0.533,
     $            0.674,0.859,1.140,1.430,1.690,1.950,2.200,2.360,
     $            2.360,2.350,2.390,2.520,2.650,2.750,2.760,2.760,
     $            2.750,2.660,2.390,2.100,2.230,2.810,3.660,4.650,
     $            5.180,6.700,7.310,7.210,22*6.890/
      pi=2.*asin(1.)

       OPEN ( 11,FILE=TLB)
       OPEN ( 12,FILE=G11)
       OPEN ( 13,FILE=G22)
       OPEN (313,FILE=G33)
       OPEN ( 14,FILE=T11)
       OPEN ( 15,FILE=T22)
       OPEN (315,FILE=T33)
       OPEN ( 16,FILE=G00)
       OPEN ( 17,FILE=T00)

      do k=1,11
        do j=1,36
         jj=36-j+1
         read(11,211)rNCEP,Z_COMMA,T_NCEP(jj,k) !TLB: breite, hoehe, zonalmittel -> used only in ncool_PWs_60.f
        end do
      end do
  211 format(1x,3e13.4)
      CLOSE(11)
      do j=1,36
        jj=36-j+1
        read(16,211)rNCEP,Z_COMMA,G0_NCEP(jj) !G0: Breite, Hoehe (=0), zonalmittel
        read(17,211)rNCEP,Z_COMMA,T0_NCEP(jj) !T0: -"-
	write(*,*) rNCEP,Z_COMMA,T0_NCEP(jj)
      end do
      do j=1,36
	jj=36-j+1
	read(12,211)rNCEP,gh1a(jj),gh1f(jj) !G1 	BREITE, AMPLITUDE, PHASE
	read(13,211)rNCEP,gh2a(jj),gh2f(jj) !G2           "
	read(313,211)rNCEP,gh3a(jj),gh3f(jj) !G3		"	
	read(14,211)rNCEP,tp1a(jj),tp1f(jj) !T1		"	
	read(15,211)rNCEP,tp2a(jj),tp2f(jj) !T2		"	
	read(315,211)rNCEP,tp3a(jj),tp3f(jj) !T3		"
      end do
      
c-------------------- read latent-heating input data

        open (40, FILE='lh_m0.dx' ,
     $      form='unformatted',
     *      access='direct', status='old',recl=2*nb*10*2)
        open (41, FILE='lh_m1.dx' ,
     $      form='unformatted',
     *      access='direct', status='old',recl=3*nb*10*4)	
        open (42, FILE='lh_m2.dx' ,
     $      form='unformatted',
     *      access='direct', status='old',recl=3*nb*10*4)
        open (43, FILE='lh_m3.dx' ,
     $      form='unformatted',
     *      access='direct', status='old',recl=3*nb*10*4)	
        open (44, FILE='lh_m4.dx' ,
     $      form='unformatted',
     *      access='direct', status='old',recl=3*nb*10*4)
     
        READ(40,REC=1) ht_m0
	READ(41,REC=1) ht_m1
	READ(42,REC=1) ht_m2
	READ(43,REC=1) ht_m3
	READ(44,REC=1) ht_m4
	
	do k=1,10
          do j=1,nb
c-------------------------------------
          Alhm0_24(j,k)=ht_m0(1,j,k,1)
          Flhm0_24(j,k)=ht_m0(1,j,k,2)
          Alhm0_12(j,k)=ht_m0(2,j,k,1)
          Flhm0_12(j,k)=ht_m0(2,j,k,2)
c-------------------------------------
          Alhm1w24(j,k)=ht_m1(1,j,k,1)
          Flhm1w24(j,k)=ht_m1(1,j,k,3)
          Alhm1w12(j,k)=ht_m1(2,j,k,1)
          Flhm1w12(j,k)=ht_m1(2,j,k,3)
          Alhm1e24(j,k)=ht_m1(1,j,k,2)
          Flhm1e24(j,k)=ht_m1(1,j,k,4)
          Alhm1e12(j,k)=ht_m1(2,j,k,2)
          Flhm1e12(j,k)=ht_m1(2,j,k,4)
c          Alhm1spw(j,k)=ht_m1(5,j,k,1)
c          Flhm1spw(j,k)=ht_m1(5,j,k,3)
c--------------------------------------
          Alhm2w24(j,k)=ht_m2(1,j,k,1)
          Flhm2w24(j,k)=ht_m2(1,j,k,3)
          Alhm2w12(j,k)=ht_m2(2,j,k,1)
          Flhm2w12(j,k)=ht_m2(2,j,k,3)
c          Alhm2spw(j,k)=ht_m2(5,j,k,1)
c          Flhm2spw(j,k)=ht_m2(5,j,k,3)
c--------------------------------------
          Alhm1e24(j,k)=ht_m2(1,j,k,2)
          Flhm1e24(j,k)=ht_m2(1,j,k,4)
          Alhm1e12(j,k)=ht_m2(2,j,k,2)
          Flhm1e12(j,k)=ht_m2(2,j,k,4)
c-------------------------------------
          Alhm3w24(j,k)=ht_m3(1,j,k,1)
          Flhm3w24(j,k)=ht_m3(1,j,k,3)
          Alhm3w12(j,k)=ht_m3(2,j,k,1)
          Flhm3w12(j,k)=ht_m3(2,j,k,3)
          Alhm3e24(j,k)=ht_m3(1,j,k,2)
          Flhm3e24(j,k)=ht_m3(1,j,k,4)
          Alhm3e12(j,k)=ht_m3(2,j,k,2)
          Flhm3e12(j,k)=ht_m3(2,j,k,4)
c-------------------------------------
          Alhm4w24(j,k)=ht_m4(1,j,k,1)
          Flhm4w24(j,k)=ht_m4(1,j,k,3)
          Alhm4w12(j,k)=ht_m4(2,j,k,1)
          Flhm4w12(j,k)=ht_m4(2,j,k,3)
          Alhm4e24(j,k)=ht_m4(1,j,k,2)
          Flhm4e24(j,k)=ht_m4(1,j,k,4)
          Alhm4e12(j,k)=ht_m4(2,j,k,2)
          Flhm4e12(j,k)=ht_m4(2,j,k,4)
c######################################
	   end do
        end do
    
        open (40, FILE='lh_a0.dx' ,
     $      form='unformatted',
     *      access='direct', status='old',recl=nb*10)
    
        do k=1,10
                do j=1,nb
                        ht_a0(j,k) = 0.
                end do
        end do
c------------------- zonally averaged LH (averaged over January 31 days)
        do l=1,248
                READ(40,REC=l) ht_a
                do k=1,10
                        do j=1,nb
                                ht_a0(j,k)=ht_a0(j,k)+ht_a(j,k)/248.
                        enddo
                enddo
	enddo
c------------------- SPW1 and SPW2 in the LH
        do k=1,10
                do j=1,nb
                        do i=1,igit
                                rlon=float(i-1)*5.625
                                heat_LH0(i,j,k)=ht_a0(j,k)
     *  +ht_m1(3,j,k,1)*cos(1.*(rlon-ht_m1(3,j,k,3))*VRAD)
     *  +ht_m2(3,j,k,1)*cos(2.*(rlon-ht_m2(3,j,k,3))*VRAD)
     *  +ht_m3(3,j,k,1)*cos(3.*(rlon-ht_m3(3,j,k,3))*VRAD)
     *  +ht_m4(3,j,k,1)*cos(4.*(rlon-ht_m4(3,j,k,3))*VRAD)
                        end do
                end do
        end do   
     
     
c-------------------- read input data on harmonics
      
      day=16.
      td1=2.*pi*day/365.
      td2=2.*td1
      td3=3.*td1
      td4=4.*td1
      
      do 200 j=1,nb
        tp1000=tpav_00(j)+tpav_s1(j)*SIN(td1)+tpav_c1(j)*COS(td1)+
     $                    tpav_s2(j)*SIN(td2)+tpav_c2(j)*COS(td2)+
     $                    tpav_s3(j)*SIN(td3)+tpav_c3(j)*COS(td3)+
     $                    tpav_s4(j)*SIN(td4)+tpav_c4(j)*COS(td4)
        tpr_m1=tpm1_r00(j)+tpm1_rs1(j)*SIN(td1)+tpm1_rc1(j)*COS(td1)+
     $                     tpm1_rs2(j)*SIN(td2)+tpm1_rc2(j)*COS(td2)+
     $                     tpm1_rs3(j)*SIN(td3)+tpm1_rc3(j)*COS(td3)+
     $                     tpm1_rs4(j)*SIN(td4)+tpm1_rc4(j)*COS(td4)
        tpi_m1=tpm1_i00(j)+tpm1_is1(j)*SIN(td1)+tpm1_ic1(j)*COS(td1)+
     $                     tpm1_is2(j)*SIN(td2)+tpm1_ic2(j)*COS(td2)+
     $                     tpm1_is3(j)*SIN(td3)+tpm1_ic3(j)*COS(td3)+
     $                     tpm1_is4(j)*SIN(td4)+tpm1_ic4(j)*COS(td4)
c--------------------------------------------------------------------
        tpr_m2=tpm2_r00(j)+tpm2_rs1(j)*SIN(td1)+tpm2_rc1(j)*COS(td1)+
     $                     tpm2_rs2(j)*SIN(td2)+tpm2_rc2(j)*COS(td2)+
     $                     tpm2_rs3(j)*SIN(td3)+tpm2_rc3(j)*COS(td3)+
     $                     tpm2_rs4(j)*SIN(td4)+tpm2_rc4(j)*COS(td4)
        tpi_m2=tpm2_i00(j)+tpm2_is1(j)*SIN(td1)+tpm2_ic1(j)*COS(td1)+
     $                     tpm2_is2(j)*SIN(td2)+tpm2_ic2(j)*COS(td2)+
     $                     tpm2_is3(j)*SIN(td3)+tpm2_ic3(j)*COS(td3)+
     $                     tpm2_is4(j)*SIN(td4)+tpm2_ic4(j)*COS(td4)
c--------------------------------------------------------------------
        gh1000=ghav_00(j)+ghav_s1(j)*SIN(td1)+ghav_c1(j)*COS(td1)+
     $                    ghav_s2(j)*SIN(td2)+ghav_c2(j)*COS(td2)+
     $                    ghav_s3(j)*SIN(td3)+ghav_c3(j)*COS(td3)+
     $                    ghav_s4(j)*SIN(td4)+ghav_c4(j)*COS(td4)
        ghr_m1=ghm1_r00(j)+ghm1_rs1(j)*SIN(td1)+ghm1_rc1(j)*COS(td1)+
     $                     ghm1_rs2(j)*SIN(td2)+ghm1_rc2(j)*COS(td2)+
     $                     ghm1_rs3(j)*SIN(td3)+ghm1_rc3(j)*COS(td3)+
     $                     ghm1_rs4(j)*SIN(td4)+ghm1_rc4(j)*COS(td4)
        ghi_m1=ghm1_i00(j)+ghm1_is1(j)*SIN(td1)+ghm1_ic1(j)*COS(td1)+
     $                     ghm1_is2(j)*SIN(td2)+ghm1_ic2(j)*COS(td2)+
     $                     ghm1_is3(j)*SIN(td3)+ghm1_ic3(j)*COS(td3)+
     $                     ghm1_is4(j)*SIN(td4)+ghm1_ic4(j)*COS(td4)
c--------------------------------------------------------------------
        ghr_m2=ghm2_r00(j)+ghm2_rs1(j)*SIN(td1)+ghm2_rc1(j)*COS(td1)+
     $                     ghm2_rs2(j)*SIN(td2)+ghm2_rc2(j)*COS(td2)+
     $                     ghm2_rs3(j)*SIN(td3)+ghm2_rc3(j)*COS(td3)+
     $                     ghm2_rs4(j)*SIN(td4)+ghm2_rc4(j)*COS(td4)
        ghi_m2=ghm2_i00(j)+ghm2_is1(j)*SIN(td1)+ghm2_ic1(j)*COS(td1)+
     $                     ghm2_is2(j)*SIN(td2)+ghm2_ic2(j)*COS(td2)+
     $                     ghm2_is3(j)*SIN(td3)+ghm2_ic3(j)*COS(td3)+
     $                     ghm2_is4(j)*SIN(td4)+ghm2_ic4(j)*COS(td4)
c--------------------------------------------------------------------
      rh2o(j)=h2oav_00(j) + h2oav_s1(j)*SIN(td1)+ h2oav_c1(j)*COS(td1)+
     $                      h2oav_s2(j)*SIN(td2)+ h2oav_c2(j)*COS(td2)+
     $                      h2oav_s3(j)*SIN(td3)+ h2oav_c3(j)*COS(td3)+
     $                      h2oav_s4(j)*SIN(td4)+ h2oav_c4(j)*COS(td4)
c--------------------------------------------------------------------
      H_wv(j)=hh2oav_00(j)+hh2oav_s1(j)*SIN(td1)+hh2oav_c1(j)*COS(td1)+
     $                     hh2oav_s2(j)*SIN(td2)+hh2oav_c2(j)*COS(td2)+
     $                     hh2oav_s3(j)*SIN(td3)+hh2oav_c3(j)*COS(td3)+
     $                     hh2oav_s4(j)*SIN(td4)+hh2oav_c4(j)*COS(td4)
c--------------------------------------------------------------------
        do i=1,igit
          dlong1    = 2.*pi*float(i-1)/float(igit) !längen in rad
          dlong2    = 2.*dlong1
c----------------------- UKMO data
c          T_1000(j,i)= tp1000+
c     $                 tpr_m2*COS(dlong2)+tpi_m2*SIN(dlong2)+
c     $                 tpr_m1*COS(dlong1)+tpi_m1*SIN(dlong1)
c----------------------- NCEP/NCAR data
         T_1000(j,i)= T0_NCEP(j)!+
!      $             1.0*tp1a(j)*COS(1.*(dlong1-tp1f(j)*pi/180.))+
!      $             1.0*tp3a(j)*COS(3.*(dlong1-tp3f(j)*pi/180.))+
!      $             1.0*tp2a(j)*COS(2.*(dlong1-tp2f(j)*pi/180.))
c--------------------------------------------------------------------
c          G_1000(j,i)= 9.81*(0.+
c          G_1000(j,i)= 9.81*gh1000
c----------------------- UKMO data
c          G_1000(j,i)= 9.81*(gh1000+
c     $                 ghr_m2*COS(dlong2)+ghi_m2*SIN(dlong2)+
c     $                (ghr_m1*COS(dlong1)+ghi_m1*SIN(dlong1)))
c----------------------- NCEP/NCAR data
          G_1000(j,i)= 9.81*G0_NCEP(j)!+
!      $    	   1.0*gh1a(j)*COS(1.*(dlong1-gh1f(j)*pi/180.))+
!      $             1.0*gh3a(j)*COS(3.*(dlong1-gh3f(j)*pi/180.))+
!      $             1.0*gh2a(j)*COS(2.*(dlong1-gh2f(j)*pi/180.)))
        end do 
 200  continue
      rb  =float(nb) 
      rgit=float(igit)
c
c.....GEOMETRIC PARAMETERS:
c     this gives altitude of highest level
      z_48 =135000.  !!!recommended!!!
      dz    = z_48/(48.-0.5)
      defi  = 180./rb   ! latitudinal increment
      dy    = a*pi/rb
      dela  = 360./rgit ! longitudinal increment
c--------------------     zunten  is the lowest level of dynamic fields 
c-----     (!!!vertical wind is calculated at ground, not at zunten!!!)
      zunten=dz*0.5
      dzn   = z_48/(float(ngit)-0.5)   !!! ngit - dzi in COMMON-BLOCK
      zunteni=dzn*0.5
C------------------------------- (g/sm^3) 
      dich0  =1.273e-3
      rhos   =1.273  
      RgSGS=8.31441e+7
      RgSI =8314.41
      C_p  =1005.
      G    =9.81
      Rf_c =0.4
cfl      co2pmv=360.
cfl      co2pmv=360.*CO2
      co2pmv=co2yr
      print*,'co2pmv=',co2pmv
      H  =7000.
      skh=7.
c--------------------------------------------------     time step
      nstep  =24*ntime          ! steps/day
      dt     =3600./float(ntime)
      dt2    =dt/2. 
      print 1045,dt
 1045 format(10x,'one time step = ',f6.1,'sec',//)
      do 1 j=1,nb+1
         gphi  = float(j-1)*defi - 87.5        ! -87.5,...,92.5 North
         fi    =(float(j)-19.5)*defi/180.*pi   ! -92.5,...,87.5 North
         fi2   =(float(j)-19.0)*defi/180.*pi   ! -90.0,...,90.0 North
c--------------------------------------------------------------------
         sinfi(j)=sin(fi)
         cosfi(j)=cos(fi)
         cosf2(j)=cos(fi2)
         sinf2(j)=sin(fi2) 
         tgfia(j)=tan(fi)/a
         cor(j)=2.*7.292116e-5*sinfi(j) !corioliskraft 2*2pi/86164*sin(phi) 
         dx(j)=2.*pi/rgit*a*cosfi(j)
         if(j.eq.nb+1) goto 1
         phi(j)=gphi*pi/180.
  1   continue
      cosf2(nb+1)=0.
      cosf2(1)=0.
         fi    =(float(nb+2)-19.5)*defi/180.*pi   ! 92.5
      cosfi(nb+2)=cos(fi)
c.......................................NEWTONIAN-COOLING-COEFFICIENTS
      dmax=100.
      d0=dmax/3.
      do k=1,kgit
         z(k)=zunten/1000.+float(k-1)*dz/1000.
         zw     = z(k)*1000.-dz/2.
         rou(k) = exp(-z(k)*1000./h)
         row(k) = exp(-zw /h)   
c...........................................     PROFILE OF NEWTONIAN COOLING:
          alfa(k)=alfaNC(k)*1.e-6
c         alfa(k)=(2.3+1.7*tanh((zu-35000.)/8000.))  *1.e-6 ! schoeberl'78
c         alfa(k)=(1.5+1.0*tanh((zu-35000.)/7000.))  *1.e-6 ! holton'76
c.........................analytical profile of the eddy diffusion coefficient
         IF(z(k).GE.90.)geddy0(k)=
     $                   dmax    *EXP(-(z(k)-90.)/20.*(z(k)-90.)/20.)
         IF(z(k).LT.90.)geddy0(k)=
     $                  (dmax-d0)*EXP(-(z(k)-90.)/20.*(z(k)-90.)/20.)+
     $                        d0 *EXP( (z(k)-90.)/20.)
      end do
      
      
! c...................................................................ION-DRAG:
! c	after Dissertation Fröhlich Eqns 3.28, 3.29, 3.30
!       if (miond.gt.0)  then
!        do k=1,kgit
! c         sigma1(k)=sigma1(k)+300.*exp((z(k)-z(kgit))/10.)
!       print*,'RAYLEIGH-ION-DRAG',' z(',k,')=', z(k),' sigma1=',sigma1(k)
!        end do
!          do j=1,nb
! 	    sinphi=sin(phi(j))
!             do k=1,kgit
!        drag(j,k,1)=sigma1(k)*1.e-6						!second part of beta_(r,lambda)
!        drag(j,k,2)=sigma1(k)*4.e-6*sinphi*sinphi/(1.+3.*sinphi*sinphi)		!second part of beta_(r,phi)
!        drag(j,k,3)=sigma2(k)*2.e-6*sinphi   /sqrt(1.+3.*sinphi*sinphi)		!Coriolos+Lorentz	
!             enddo
!          enddo
!       endif
! c.....RAYLEIGH-FRICTION:
!       if (miond.gt.0)     then
!          do j=1,nb
!          do k=1,kgit
! cfl	    Brayl0=1.E-8   ! why fix???
! c	    z(k)-z0 -> z0=40 or 60km referring to strong or weak coefficient of the horizontal diffusion in the stratosphere
! cfl	    earlier: BRayl0 = (1.+0.5*tanh((z(k)-50.)/20.))*1.e-6	
!             BRayl0 = (1.25+0.75*tanh((z(k)-60.)/20.))*1.e-6	!background rayleigh frisction coefficient beta_r    	
!             drag(j,k,1)=drag(j,k,1)+BRayl0 !beta_(r,lambda)
!             drag(j,k,2)=drag(j,k,2)+BRayl0 !beta_(r,phi)      
!          enddo
!          enddo
!       endif

c...................................................................ION-DRAG:
c	FL after Dissertation Fröhlich Eqns 3.28, 3.29, 3.30
      !Initialisierung:
      do j=1,nb
	do k=1,kgit
	  drag(j,k,1)=0.						
	  drag(j,k,2)=0.		
	  drag(j,k,3)=0.
	enddo
      enddo	
c.....BACKGROUND RAYLEIGH-FRICTION:
      if (miond.ge.1)  then	
        do j=1,nb
	    sinphi=sin(phi(j))
            do k=1,kgit
cfl	      Brayl0=1.E-8   ! why fix???
c	      z(k)-z0 -> z0=40 or 60km referring to strong or weak coefficient of the horizontal diffusion in the stratosphere
cfl	      earlier: BRayl0 = (1.+0.5*tanh((z(k)-50.)/20.))*1.e-6	              
	      BRayl0 = (1.25+0.75*tanh((z(k)-60.)/20.))*1.e-6	!background rayleigh friction coefficient beta_r  
	      !BRayl0 = 1.e-8
	      drag(j,k,1)=drag(j,k,1)+BRayl0 !beta_(r,lambda)
	      drag(j,k,2)=drag(j,k,2)+BRayl0 !beta_(r,phi) 
            enddo
         enddo
      endif
C.....ION DRAG           
      if (miond.ge.2)     then
        do k=1,kgit
c         sigma1(k)=sigma1(k)+300.*exp((z(k)-z(kgit))/10.)
	  print*,'ION-DRAG:',' z(',k,')=', z(k),' sigma1=',
     + 	  sigma1(k)
        end do 
        do j=1,nb
          do k=1,kgit
            drag(j,k,1)=drag(j,k,1)+sigma1(k)*1.e-6			!second part of beta_(r,lambda)
	    drag(j,k,2)=drag(j,k,2)+sigma1(k)*4.e-6*sinphi/
     +	      (1.+3.*sinphi*sinphi)		!second part of beta_(r,phi)	
          enddo
        enddo
      endif
c......LORENTZ DEFLECTION (to be subtracted from Coriolis force)     
      if (lodef.ge.1) then
	do k=1,kgit
          print*,'LORENTZ DEFLECTION:',' z(',k,')=', z(k),' sigma2=',
     +    sigma2(k)
        end do 
	do j=1,nb
          do k=1,kgit
	    drag(j,k,3)=drag(j,k,3)+sigma2(k)*2.e-6*sinphi/
     +	      sqrt(1.+3.*sinphi*sinphi)	*sinphi	!Lorentz deflection term multiplied by phi due to cor(j)=2*omega*sin(phi) 
          enddo
        enddo  
      endif
      
      
c................................GLOBAL MEAN TEMPERATURE tz(kgit)
c----------------- initial setting, will be improved in strobel.f
      do k=1,kgit
         tz(k)=0.
         do j=1,nb
            do i=1,igit     
               tz(k)=tz(k)+tn1(j,k,i)/rgit/rb
            enddo
         enddo
      end do
C................................... DENSITY and related PROFILES
c----------------- initial setting, will be improved in strobel.f
c---------- T(0), pressure at the lower boundary is 1.013e6 - SGS
         t0=1.013e6*rm(1)/RgSGS/dich0
      do k=1,kgit
      	 print*,k,RgSGS,rm(k),tz(k)
      	 Hsclkm=RgSGS/rm(k)*tz(k)*1.e-7/9.81
         dichte(k)=dich0*rou(k)*rm(k)/rm(1)*t0/tz(k)
c------------------------------------- xmt will be used in grav.f
         xmt(k)=rm(k)/rm(1)*t0/tz(k)
c------------------------------------------     for UV radiation:
c
c----------------- Cp(z) - approx. 1.e7 (SGS) in lower atmosphere
c
         cp_z=RgSGS/rm(k)*(7./2.*(1.-roxvmr(k))+5./2.*roxvmr(k))
         dichcp(k)=dichte(k)*cp_z
         sunris(k)=acos(6373./(6373.+z(k)))*180./pi + 90.
         do j=1,nb
c
c------------- initial setting, will be recalculated in strobel.f 
c
            xnumo3(j,k)=1.013  *rou(k)*ozppmv(j,k)/1.38066e-16/tz(k)
            xnumn2(j,k)=1.013e6*rou(k)*rn2vmr(  k)/1.38066e-16/tz(k)
            xnumo2(j,k)=1.013e6*rou(k)*ro2vmr(  k)/1.38066e-16/tz(k)
            xnumox(j,k)=1.013e6*rou(k)*roxvmr(  k)/1.38066e-16/tz(k)
            sox(j,k)=Hsclkm*rm(k)/16.*xnumox(j,k)*1.e5
            so2(j,k)=Hsclkm*rm(k)/32.*xnumo2(j,k)*1.e5
            sn2(j,k)=Hsclkm*rm(k)/28.*xnumn2(j,k)*1.e5
         end do
      end do
c----     M/M0*T0/T
      xmt(kgit+1)=(8.*xmt(kgit)-3.*xmt(kgit-1))/5.
c
c....................................................H2O:
c      ratio=28.9644/18.016 !mol_luft/mol_h2o
c------------------------     this is the default value at ground
c--------- new value Forbes and Garrett (1979)
c      ppm=0.0175 
c      ppm=0.02 !max.massemmisch.verh.(=0.032ppv)
c      ppm=0.03 !max.massemmisch.verh.(=0.048ppv)
c
c      c2_h2o=2.e-6
      do 12 i=1,igit
      do 12 j=1,nb
      D_1000=1.e6*rm(1)/RgSGS/T_1000(j,i)
c      D_1000=1.e6*rm(1)/RgSGS/273.
c         gphi      = -87.5+float(j-1)*5.       ! -87.5,...,87.5 North
c---------------------ANALITICAL H2O(phi)
c         c1_h2o(i,j)= ppm*exp (-gphi*gphi/1600.) 
c        a1h2o(i,j) = exp (-gphi*gphi/1600.)*ppm*dich0 *1.e5 
c        a1h2o(i,j) = exp (-gphi*gphi/1600.)*ppm*D_1000*1.e5 
c------------------- rh2o(j) is the empirical (NCEP/NCAR) distribution
c              D_1000 (g/cm^3), 1.e5 (km - cm), a1h2o (g/cm^2/km)
c------------------NCEP/NCAR H2O DATA---------------------------------
       a1h2o(i,j)= rh2o(j)*D_1000*1.e5
 12   continue
c---------------------------------------------  factors for h2o column
c------------------- will be used in ircool_newCO2.f and strobel.f
      do 13 kk=1,4
            skhnir(kk) = skh/(1.+xncir(kk))
            do j=1,nb
              e1nir (kk,j) = -1./skhnir(kk) - 1./H_wv(j)
            end do
 13   continue
      do 14 kk=1,9
            skhnuv(kk) = skh/(1.+dkli(kk))
            do j=1,nb
              expnuv(kk,j)  = -1./skhnuv(kk) - 1./H_wv(j)
            end do 
 14   continue
c---------------------------------- new H2O heating???????
c      do 15 kk=1,53
c        H_h2o_2(kk) = skh/(1.+cm_h2o(kk))
c          do j=1,nb
c          H_h2o_1(kk,j)  = 1./(1./H_h2o_2(kk)+1./H_wv(j))
c          end do 
c----------------- scaling factors for each interval
c      s1_h2o(kk)=(1013./p0_h2o(kk))**cm_h2o(kk) 
c 15   continue
       return
      end









