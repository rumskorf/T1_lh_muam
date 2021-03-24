      subroutine plwaves_LIM
      include 'com_main.fc'
      include 'com_namelist.fc'
      real fphiKW(nb),fphi2d(nb),fphi11(nb),fphi12(nb),fphi13(nb),
     $     fphi21(nb),fphi23(nb),fphiFK(nb),fphiSK(nb),fphiNK(nb),
     $     fphi10(nb)
c      real fi0_2d,fi0_KW,fi0_FK,fi0_11,fi0_12,fi0_13,fi0_21,fi0_23,
      real cos_2d,cos_KW,cos_FK,cos_11,cos_12,cos_13,cos_21,cos_23,
     $     cos_SK
c
c      integer I_2d,I_SK,I_KW,I_FK,I_NK,I_10,I_11,I_12,I_13,I_21
c
c     data I_2d/0/,I_SK/0/,I_KW/0/,I_FK/0/,I_NK/0/,
c     $     I_10/0/,I_11/1/,I_12/1/,I_13/1/,
c     $     I_21/1/,I_23/0/
c----------------------- Hough functions - latitudinal structure
      data fphi2d/3.49e-05,9.34e-04,4.26e-03,1.14e-02,2.34e-02,4.07e-02,
     $            6.33e-02,9.00e-02,1.19e-01,1.48e-01,1.73e-01,1.90e-01,
     $            1.97e-01,1.90e-01,1.69e-01,1.33e-01,8.51e-02,2.93e-02,
     $      -2.93e-02,-8.51e-02,-1.33e-01,-1.69e-01,-1.90e-01,-1.97e-01,
     $      -1.90e-01,-1.73e-01,-1.48e-01,-1.19e-01,-9.00e-02,-6.33e-02,
     $      -4.07e-02,-2.34e-02,-1.14e-02,-4.26e-03,-9.34e-04,-3.49e-05/
      data fphiSK/-.29124E-04,.24076E-04,.78733E-05,.21911E-05,
     $             .39120E-05,.21749E-05,.56963E-05,.28159E-04,
     $             .17179E-03,.92860E-03,.42880E-02,.16833E-01,
     $             .55505E-01,.15226E+00,.34469E+00,.63982E+00,
     $             .96901E+00,.11935E+01,.11935E+01,.96901E+00,

     $             .63982E+00,.34469E+00,.15226E+00,.55505E-01,
     $             .16833E-01,.42880E-02,.92859E-03,.17182E-03,
     $             .28112E-04,.57196E-05,.21461E-05,.38822E-05,
     $             .22303E-05,.78517E-05,.24066E-04,-.29127E-04/
c--------------- UFKW, T=3.5 days
c      data fphiKW/0.00053643,0.0017932,0.0036518,0.0067036,0.011848,
c     $            0.020452,0.034507,0.056729,0.0905,0.13953,0.20708,
c     $            0.29483,0.40145,0.52135,0.64432,0.75644,0.84247,
c     $            0.88928,0.88928,0.84247,0.75644,0.64432,0.52135,
c     $            0.40145,0.29483,0.20708,0.13953,0.0905,0.056729,
c     $            0.034507,0.020452,0.011848,0.0067036,0.0036518,
c     $            0.0017932,0.00053643/
c--------------- KW, T=5 days
      data fphiKW/.50352E-04,.17973E-03,.43180E-03,.97296E-03,
     $            .21420E-02,.46122E-02,.96556E-02,.19510E-01,
     $            .37795E-01,.69753E-01,.12197E+00,.20107E+00,
     $            .31113E+00,.45023E+00,.60739E+00,.76200E+00,
     $            .88729E+00,.95777E+00,.95777E+00,.88729E+00,
     $            .76200E+00,.60739E+00,.45023E+00,.31113E+00,
     $            .20107E+00,.12197E+00,.69753E-01,.37795E-01,
     $            .19510E-01,.96555E-02,.46122E-02,.21420E-02,
     $            .97297E-03,.43180E-03,.17973E-03,.50352E-04/
      data fphiFK/-.46899E-05,-.54356E-05,-.17457E-04,-.57160E-04,
     $            -.17727E-03,-.52825E-03,-.15137E-02,-.41173E-02,
     $            -.10539E-01,-.25168E-01,-.55638E-01,-.11308E+00,
     $            -.21004E+00,-.35470E+00,-.54224E+00,-.74776E+00,
     $            -.92779E+00,-.10339E+01,-.10339E+01,-.92779E+00,
     $            -.74776E+00,-.54224E+00,-.35470E+00,-.21004E+00,
     $            -.11308E+00,-.55638E-01,-.25168E-01,-.10539E-01,
     $            -.41174E-02,-.15136E-02,-.52826E-03,-.17727E-03,
     $            -.57147E-04,-.17473E-04,-.54179E-05,-.46896E-05/
      data fphiNK/.10557E-01,.32257E-01,.55699E-01,.82014E-01,
     $            .11227E+00,.14743E+00,.18823E+00,.23512E+00,
     $            .28810E+00,.34664E+00,.40951E+00,.47480E+00,
     $            .53993E+00,.60177E+00,.65690E+00,.70192E+00,
     $            .73382E+00,.75037E+00,.75037E+00,.73382E+00,
     $            .70192E+00,.65690E+00,.60177E+00,.53993E+00,
     $            .47480E+00,.40951E+00,.34664E+00,.28810E+00,
     $            .23512E+00,.18823E+00,.14743E+00,.11227E+00,
     $            .82014E-01,.55699E-01,.32257E-01,.10557E-01/
      data fphi10/-.17989E-01,-.54185E-01,-.90979E-01,-.12856E+00,
     $            -.16682E+00,-.20521E+00,-.24270E+00,-.27769E+00,
     $            -.30806E+00,-.33128E+00,-.34458E+00,-.34532E+00,
     $            -.33134E+00,-.30132E+00,-.25522E+00,-.19442E+00,
     $            -.12180E+00,-.41487E-01, .41487E-01, .12180E+00,
     $             .19442E+00, .25522E+00, .30132E+00, .33134E+00,
     $             .34532E+00, .34458E+00, .33128E+00, .30806E+00,
     $             .27769E+00, .24270E+00, .20521E+00, .16682E+00,
     $             .12856E+00, .90979E-01, .54185E-01, .17989E-01/
      data fphi11/0.043283,0.12906,0.21237,0.29134,0.36378,0.42721,

     $            0.47895,0.51639,0.53734,0.54047,0.52574,0.49475,
     $            0.45094,0.39939,0.34643,0.29885,0.26298,0.24371,
     $            0.24371,0.26298,0.29885,0.34643,0.39939,0.45094,

     $            0.49475,0.52574,0.54047,0.53734,0.51639,0.47895,
     $            0.42721,0.36379,0.29134,0.21237,0.12906,0.043283/
      data fphi12/0.074316,0.21891,0.35153,0.46471,0.55183,0.60769,
     $            0.62906,0.61531,0.56886,0.49534,0.40328,0.30323,
     $            0.20637,0.12278,0.059682,0.020024,0.0018908,
     $           -0.0011728,0.0011728,-0.0018908,-0.020024,-0.059682,
     $           -0.12278,-0.20637,-0.30323,-0.40328,-0.49534,
     $           -0.56886,-0.61531,-0.62906,-0.60769,-0.55184,
     $           -0.46471,-0.35153,-0.21891,-0.074316/
      data fphi13/-0.1012,-0.29359,-0.45698,-0.57553,-0.63838,-0.64115,
     $            -0.58668,-0.48502,-0.35217,-0.2078,-0.071996,
     $             0.038239,0.11222,0.14749,0.14994,0.13184,0.10822,
     $             0.092488,0.092488,0.10822,0.13184,0.14994,0.14749,
     $             0.11222,0.038239,-0.071996,-0.2078,-0.35217,
     $            -0.48502,-0.58668,-0.64115,-0.63838,-0.57553,
     $            -0.45698,-0.29359,-0.1012/
      data fphi21/-0.0022016,-0.019591,-0.053164,-0.10051,-0.15802,
     $            -0.22101,-0.2839,-0.34064,-0.38534,-0.41296,-0.42022,
     $            -0.40632,-0.37339,-0.32654,-0.27329,-0.22245,
     $            -0.18268,-0.16089,-0.16089,-0.18268,-0.22245,
     $            -0.27329,-0.32654,-0.37339,-0.40632,-0.42022,

     $            -0.41296,-0.38534,-0.34064,-0.2839,-0.22101,
     $            -0.15802,-0.10051,-0.053165,-0.019591,-0.0022017/
      data fphi23/.64766E-02, .56259E-01, .14543E+00, .25519E+00,
     $            .36205E+00, .44259E+00, .47845E+00, .46043E+00,
     $            .39079E+00, .28285E+00, .15786E+00, .39635E-01,
     $           -.51618E-01,-.10459E+00,-.11949E+00,-.10707E+00,
     $           -.84259E-01,-.67737E-01,-.67737E-01,-.84259E-01,
     $           -.10707E+00,-.11949E+00,-.10459E+00,-.51618E-01,
     $            .39635E-01, .15786E+00, .28285E+00, .39079E+00,
     $            .46043E+00, .47845E+00, .44259E+00, .36205E+00,
     $            .25519E+00, .14543E+00, .56260E-01, .64766E-02/


c      OPEN(UNIT=10,FILE='muam_mod.txt')

c      READ(10,PWS)
c      READ(10,PWX)
c      READ(10,PWP)
c      READ(10,PWF)

c      CLOSE(10)

c      print*,I_2d,I_SK,I_KW,I_FK,I_NK,I_10,I_11,I_12,I_13,I_21,I_23
c      print*,Xm2d,XmSK,XmFK,XmKW,XmNK,Xm10,Xm11,Xm12,Xm13,Xm21,Xm23
c      print*,Per_2d,Per_SK,Per_FK,Per_KW,Per_NK,Per_10,Per_11,Per_12,
c     #       Per_13,Per_21,Per_23
c      print*,fi0_2d,fi0_SK,fi0_FK,fi0_KW,fi0_NK,fi0_10,fi0_11,fi0_12,
c     #       fi0_13,fi0_21,fi0_23

cFL for DAU: all number greater than 1 are set to 1
      if(I_2d.gt.1) I_2d=1
      if(I_SK.gt.1) I_SK=1
      if(I_FK.gt.1) I_FK=1 
      if(I_KW.gt.1) I_KW=1
      if(I_NK.gt.1) I_NK=1
      if(I_10.gt.1) I_10=1
      if(I_11.gt.1) I_11=1
      if(I_12.gt.1) I_12=1
      if(I_13.gt.1) I_13=1
      if(I_21.gt.1) I_21=1
      if(I_23.gt.1) I_23=1
c-------------------------------------------------------------------------
c      pi=2.*asin(1.)
      do j=1,nb
        do i=1,igit
         do k=1,kgit
           h_PWs(k,i,j)=0.
         end do 
        end do
      end do
c-----------------------------------------------------------------
      nfor=max0(ncom-nphi,0) 
      IF(nfor.gt.0) then
      xsec  = 3600./float(ntime)      ! seconds per step (ntime=8 -> xsec=450s)
c-----------------------------------------------------------------
c        *****  the time dependence of amplitudes  *****
c-----------------------------------------------------------------
      tdfor = float(nfor)*xsec/86400. ! time (days) of the forcing
      tPWs  = 1.- exp(-tdfor/10.)     ! travelling PWs 
C--------------------- MAX of Travelling PWs on July 10!!!!!!!!!!!
c      nmax=21120-nphi
c      tdmax = float(nmax)*xsec/86400.
c------------     a smooth start of forcing ....
c      tPWs_1 = 1.+EXP(-(tdfor-tdmax)*(tdfor-tdmax)/10./10.)
c      fi0_2d = fi0_2d*(1.-0.5*COS(2.*pi*tdfor/12.))
c      fi0_KW = fi0_KW*(1.-0.5*COS(2.*pi*tdfor/15.))
c      fi0_11 = fi0_11*(1.-0.5*COS(2.*pi*tdfor/30.))
c
c--------------- the 2-day wave, (3,0) mode
c--------------- amplitude of mer. wind is about 8 m/s in stratosphere 
!       Xm2d   = 3.
!       Per_2d =-52.5*3600. 
! 
! c      Per_2d =-56.*3600. 
! 
! c      fi0_2d = 1500.*float(I_2d)     
! c-----------  ampilude too strong, to investigate the resonant response
!       fi0_2d = 2.9e-4*float(I_2d)     
! c---- weak 2-day wave - has been used to simulate the resonant response
! c      fi0_2d = 1.e-4*float(I_2d)     
! c--------------- slow Kelvin wave
! c--------------- amplitude of geop. Height is about 50 m at 50 km
!       XmSK   = 1.
!       Per_SK = 14.0*24.*3600.
!       fi0_SK =  200.*float(I_SK)
! c--------------- fast Kelvin wave
! c--------------- amplitude of geop. Height is about 50 m at 50 km
!       XmFK   = 1.
!       Per_FK = 7.0*24.*3600.
!       fi0_FK = 100.*float(I_FK)
! c--------------- ultra-fast Kelvin wave
! c--------------- amplitude of geop. Height is about 50 m at 50 km
!       XmKW   = 1.
!       Per_KW = 3.75*24.*3600.
! c      fi0_KW =   30.*float(I_KW)
! c      Per_KW = 5.*24.*3600.
! c      fi0_KW =   10.*float(I_KW)
!       fi0_KW =   0.024e-4*float(I_KW)
! c--------------- Normal mode Kelvin wave, T=1.367*24=32.808 (h_n=9.8 km)
! c--------------- amplitude of geop. height is about ?? m at 50 km
!       XmNK   = 1.
!       Per_NK = 1.367*24.*3600.*0.8
!       fi0_NK =   0.02e-4*float(I_NK)
! c-----------------------------------------------------------------------
! c--------------- the 30h wave, (1,0) mode, T=-1.19*24=28.56 (h_n=9.8 km)
! c--------------- amplitude of geop. height is about ?? m at 50 km
!       Xm10   = 1.
!       Per_10 =-1.19*24.*3600.*0.8
!       fi0_10 =  0.04e-4*float(I_10)
! c----------------------------------------------------------------------- 
! c--------------- the 5-day wave, (1,1) mode
! c--------------- amplitude of geop. Height is about 150 m at 50 km
!       Xm11   = 1.
!       Per_11 =-120.*3600.*0.8
! c      fi0_11 =  0.5*400.*float(I_11) 
!       fi0_11 =  0.20e-4*float(-1*I_11)
! c--------------- the 10-day wave, (1,2) mode
!       Xm12   = 1.
!       Per_12 =-220.*3600.
! c      fi0_12 = 2.*300.*float(I_12) 
! 
!       fi0_12 = 0.37e-4*float(-1*I_12)
! Cmpw      fi0_12 = 0.37e-4*float(I_12) 
! Clpw      fi0_12 = 0.37e-4*float(I_12)/2.
! c--------------- the 16-day wave, (1,3) mode
!       Xm13 = 1.
!       Per_13 =-360.*3600.
! c      fi0_13 = 2.*600.*float(I_13) 
! 
!       fi0_13 = 0.5e-4*float(-1*I_13)
! cmpw      fi0_13 = 0.5e-4*float(I_13) 
! Clpw      fi0_13 = 0.5e-4*float(I_13)/4.
! c--------------- the 4-day wave, (2,1) mode
!       Xm21   = 2.
!       Per_21 =-96.*3600.
! 
! c      fi0_21 =  400.*float(I_21) 
!       fi0_21 =  0.55e-4*float(-1*I_21) 
! c--------------- the 15-day wave, (2,3) mode
!       Xm23   = 2.
!       Per_23 = 360.*3600.
!       fi0_23 =  1000.*float(I_23) 
c-----------------------------------------------------------------
        if(ncom.eq.null+1)then
           if(I_2d.eq.1)print*,'QDTW, Heating_Amp=', fi0_2d
           if(I_SK.eq.1)print*,'SLow KW, Heating_Amp=', fi0_SK
           if(I_FK.eq.1)print*,'Fast KW, Heatimg_Amp=', fi0_FK 
           if(I_KW.eq.1)print*,'UltraFast KW,  Heat.Amp=', fi0_KW 
           if(I_NK.eq.1)print*,'Norm. mode KW, Heat.Amp=', fi0_NK 
           if(I_10.eq.1)print*,'30h m=1 NM, Heating_Amp=', fi0_10
           if(I_11.eq.1)print*,'5 DW, Heating_Amp=', fi0_11
           if(I_12.eq.1)print*,'10 DW, Heating_Amp=', fi0_12
           if(I_13.eq.1)print*,'16 DW, Heating_Amp= ', fi0_13
           if(I_21.eq.1)print*,'4 DW, Heating_Amp= ', fi0_21
           if(I_23.eq.1)print*,'16 DW, Heating_Amp= ', fi0_13
        endif
c
c     *****  the phases = omega*t  ********************
c
        fase2d = float(nfor)*xsec*2.*pi/per_2d
        faseSK = float(nfor)*xsec*2.*pi/per_SK
        faseFK = float(nfor)*xsec*2.*pi/per_FK
        faseKW = float(nfor)*xsec*2.*pi/per_KW
        faseNK = float(nfor)*xsec*2.*pi/per_NK
        fase10 = float(nfor)*xsec*2.*pi/per_10
        fase11 = float(nfor)*xsec*2.*pi/per_11
        fase12 = float(nfor)*xsec*2.*pi/per_12
        fase13 = float(nfor)*xsec*2.*pi/per_13
        fase21 = float(nfor)*xsec*2.*pi/per_21
        fase23 = float(nfor)*xsec*2.*pi/per_23
c
c     *****  geopotential height perturbations at the lower boundary 
c
      do 100 j= 1,nb
         do i=1,igit
             cos_2d= cos(xm2d*2.*pi*float(i-1)/float(igit)+fase2d)
             cos_SK= cos(xmSK*2.*pi*float(i-1)/float(igit)+faseSK)
             cos_FK= cos(xmFK*2.*pi*float(i-1)/float(igit)+faseFK)
             cos_KW= cos(xmKW*2.*pi*float(i-1)/float(igit)+faseKW)
             cos_NK= cos(xmNK*2.*pi*float(i-1)/float(igit)+faseNK)
             cos_10= cos(xm10*2.*pi*float(i-1)/float(igit)+fase10)
             cos_11= cos(xm11*2.*pi*float(i-1)/float(igit)+fase11)
             cos_12= cos(xm12*2.*pi*float(i-1)/float(igit)+fase12)
             cos_13= cos(xm13*2.*pi*float(i-1)/float(igit)+fase13)
             cos_21= cos(xm21*2.*pi*float(i-1)/float(igit)+fase21)
             cos_23= cos(xm23*2.*pi*float(i-1)/float(igit)+fase23)
c-----------------------------------------------------------------
c-------- tPWs   is factor for the establishment of travelling Pws
c-------- tPWs_1 is factor for increase in ampl. of travelling Pws
c-----------------------------------------------------------------

c------------------------------------------------------------------------
c FL: previous version all waves on, no SK,Fk
c------------------------------------------------------------------------

c	    do k=1,kgit
c               vertst=EXP(-(z(k)-10.)*(z(k)-10.)/25.)
c               h_PWs(k,i,j)=tPWs*vertst*(fi0_2d*fphi2d(j)*cos_2d+
c     $                                   fi0_10*fphi10(j)*cos_10+
c     $                                   fi0_11*fphi11(j)*cos_11+
c     $                                   fi0_12*fphi12(j)*cos_12+
c     $                                   fi0_13*fphi13(j)*cos_13+
c     $                                   fi0_NK*fphiNK(j)*cos_NK+
c     $                                   fi0_KW*fphiKW(j)*cos_KW+
c     $                                   fi0_21*fphi21(j)*cos_21+
c     $                                   fi0_23*fphi23(j)*cos_23)
c             end do 

c------------------------------------------------------------------------
c FL new version with on/off buttons
c------------------------------------------------------------------------
	    
             do k=1,kgit
               vertst=EXP(-(z(k)-10.)*(z(k)-10.)/25.) !vertst increases until 10km (to 1) and decreases above (at 20km close to 0)
               h_PWs(k,i,j)=tPWs*vertst*(fi0_2d*fphi2d(j)*cos_2d*I_2d+
     $                                   fi0_10*fphi10(j)*cos_10*I_10+
     $                                   fi0_11*fphi11(j)*cos_11*I_11+
     $                                   fi0_12*fphi12(j)*cos_12*I_12+
     $                                   fi0_13*fphi13(j)*cos_13*I_13+
     $					 fi0_SK*fphiSK(j)*cos_SK*I_SK+ 
     $                                   fi0_NK*fphiNK(j)*cos_NK*I_NK+
     $					 fi0_FK*fphiFK(j)*cos_FK*I_FK+
     $                                   fi0_KW*fphiKW(j)*cos_KW*I_KW+
     $                                   fi0_21*fphi21(j)*cos_21*I_21+
     $                                   fi0_23*fphi23(j)*cos_23*I_23)
             end do 
         enddo 
 100  continue
      END IF
      return
      end



