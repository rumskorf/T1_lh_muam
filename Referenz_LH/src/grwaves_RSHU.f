c ---------------------------------------------------------------------
c      version by Alex Pogoreltsev, 2005 - see comments in the text
c----------------------------------------------------------------------
      subroutine grwaves_RSHU
C     =================================================================
C     GRAVITY WAVE DRAG parameterization based on an analytical solution
C     (WKB approximation) of the vertical structure equation for the GWs
C     in the atmosphere with arbitrary background wind U0(z), V0(z) and 
C     realistic radiative damping (Newtonian cooling).
C----------------------------------------------------------------------  
C     Eddy diffusion coefficient is estimated using the idea of GW breaking
C     due to instability proposed by Lindzen (1981). 
C     Heating/cooling rates are calculated following 
C     Medvedev and Klaassen (2002), however, the efficiency of mechanical
C     energy conversion into heat is taken into account, and turbulent
C     Prandl number = 3 (Huang und Smith, 1991).
C     Heat flux and its divergence are calculated explicitly,
C     using the polarization relations for GWs and analytical solution.  
C     ================================================================
      include 'com_main.fc'
      include 'com_namelist.fc'
      real    u0(kgit),v0(kgit),t0(kgit),ac0(kgit),ac1(kgit),
     $        flux0(kgit), flux1(kgit),      
     $        Fheat1(kgit), heatc1(kgit), Fheatc(kgit), heatco(kgit),
     $        u0z(kgit),v0z(kgit),t0z(kgit), u0zz(kgit), v0zz(kgit),
     $        cphint(kgit), omb2(kgit),derkz2(kgit), 
     $        debg0(kgit),debgm(kgit),geddym(kgit,igit,nb),
     $        degw(kgit),alfaNC(kgit),w0(48),
     $        weight(6),ww(6,36) 
      real*8  forceU(kgit), forceV(kgit), forceT(kgit), 
     $        work(6,igit), fluxU(kgit), fluxV(kgit)                   ! by Rliu 
      
      real ep, lat
      character*2 mon
      character*17 datei
	  integer nmon
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c             cphase = phase speed [m/s]
c             xl     = horizontal wave length [m]
c             teta   = azimuth of GW [deg]
c             ww     = vertical velocity perturbations at lower bondary [m/s]
c             pr     = 3! Prandl number (Huang and Smith, 1991)
c             eff    = efficiency of impuls fluxes 
c             Rg     = gas constant
c             gam    = Cp/Cv
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      data nwaves /48/, eff0/1./
CPH      data cphase/5.,10.,15.,20.,25.,30./,
CPH     $     teta/0.,45.,90.,135.,180.,225.,270.,315./
     
c      data cphase/2.5,5.,7.5,10.0,12.5,15.0/,
c     $     teta/0.,45.,90.,135.,180.,225.,270.,315./
           
    
c----------------------------------------------------- the amplitudes
CPH      data Rg/287.04/,gam/1.4/,xl/300000./,pr/3./
      data Rg/287.04/,gam/1.4/,pr/3./

      data omega0/1.e-4/,alpha/0.8333/,gamma/1.5/,c1/10./,c2/4./

c      NAMELIST /GRW/ XL,WAMPL,CPHASE,TETA
c      OPEN(UNIT=10,FILE='muam_mod.txt')

c      READ(10,GRW)

c      CLOSE(10)

c      print*,XL,WAMPL,CPHASE,TETA


c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      rkx=2.*pi/xl
      rkzave=0. !FL: only initializing if not set for first steps      
c------------------------------ Gavrilov
      do i=1,6
        omegai=cphase(i)*rkx
        omrel = omegai/omega0
        c1rel = cphase(i)/c1
        c2rel = c2/cphase(i)
        weight(i)=cphase(i)/(1.+omrel**alpha)/(1.+c1rel+c2rel**gamma)
      end do
      
!Open the file of the potential energy

      if(NDEK.le.31)                  mon='01'
      if(NDEK.ge. 32.and.NDEK.le. 59) mon='02'
      if(NDEK.ge. 60.and.NDEK.le. 90) mon='03'
      if(NDEK.ge. 91.and.NDEK.le.120) mon='04'
      if(NDEK.ge.121.and.NDEK.le.151) mon='05'
      if(NDEK.ge.152.and.NDEK.le.181) mon='06'
      if(NDEK.ge.182.and.NDEK.le.212) mon='07'
      if(NDEK.ge.213.and.NDEK.le.243) mon='08'
      if(NDEK.ge.244.and.NDEK.le.273) mon='09'
      if(NDEK.ge.274.and.NDEK.le.304) mon='10'
      if(NDEK.ge.305.and.NDEK.le.334) mon='11'
      if(NDEK.ge.335.and.NDEK.le.365) mon='12'

	  datei='zonmean'//mon//'.'//GWAM//'.dat'
      !datei='zonmean01.'//GWAM//'.dat'
      print*,datei
      open(1,file='../ep_zon/'//datei)
      

!Calculation of the GW momentum flux via midfrequency analysis

      do j=1,nb
      
       read(1,197) lat,ep    
       !print*,lat,ep
  197  format(E10.7,1x,E10.7)
         
       do i=1,6
	 
csl-----------------------------------------------------------------------------
csl	schwerewellen amplitude einstellen -> verschiebung der windumkehr faktor 0.02 kleiner -> verschiebung nach oben

c--------------------------- for January
c        wampl0 = 0.01*(1.+1./3.*tanh((float(j)-18.5)/6.))  !lgw            ! by Rliu
c	     wampl0 = 0.015*(1.+1./3.*tanh((float(j)-18.5)/6.)) !mgw
c        wampl0 = 0.02*(1.+1./3.*tanh((float(j)-18.5)/6.)) !Jan
c        wampl0 = 0.022*(1.+1./3.*tanh((float(j)-18.5)/6.)) !hgw
c        wampl0 = 0.03*(1.+1./3.*tanh((float(j)-18.5)/6.)) !egw
       
c--------------------------- for July
c        wampl0 = 0.025-0.005*tanh((float(j)-18.5)/9.) 
c        wampl0 = 0.015-0.005*tanh((float(j)-18.5)/9.)
c        wampl0 = 0.02-0.005*tanh((float(j)-18.5)/9.) 

c--------------------------- for April
c        wampl0 = 0.015
	     wampl0 = ep*WAMPX

! 	  if(GWAM.eq.'Jan') then 
! 	    wampl0 = WAMPX*(1.+1./3.*tanh((float(j)-18.5)/6.))
! 	  elseif(GWAM.eq.'Feb') then
! 	    wampl0 = WAMPX*(1.+1./5.*tanh((float(j)-18.5)/6.))
! 	  elseif(GWAM.eq.'Mar') then
! 	    wampl0 = WAMPX*(1.+1./12.*tanh((float(j)-18.5)/6.))
! 	  elseif(GWAM.eq.'Apr') then
! 	    wampl0 = WAMPX
! 	  elseif(GWAM.eq.'May') then 
! 	    wampl0 = WAMPX-0.0015*tanh((float(j)-18.5)/9.)
! 	  elseif(GWAM.eq.'Jun') then
! 	    wampl0 = WAMPX-0.003*tanh((float(j)-18.5)/9.)
! 	  elseif(GWAM.eq.'Jul') then
! 	    wampl0 = WAMPX-0.005*tanh((float(j)-18.5)/9.)
! 	  elseif(GWAM.eq.'Aug') then
! 	    wampl0 = WAMPX-0.003*tanh((float(j)-18.5)/9.)
! 	  elseif(GWAM.eq.'Sep') then
! 	    wampl0 = WAMPX-0.0015*tanh((float(j)-18.5)/9.)
! 	  elseif(GWAM.eq.'Oct') then
! 	    wampl0 = WAMPX
! 	  elseif(GWAM.eq.'Nov') then
! 	    wampl0 = WAMPX*(1.+1./12.*tanh((float(j)-18.5)/6.))
! 	  elseif(GWAM.eq.'Dec') then
! 	    wampl0 = WAMPX*(1.+1./5.*tanh((float(j)-18.5)/6.))
! 	  endif

csl--------------------------------------------------------------------------------------------------------------

	    ww(i,j)=wampl0*weight(i)
        end do
      end do
      close(1)
c--------------------------------------------------------------------
c      h = 7000.
c      Cp = Rg*gam/(gam-1.)
c------------------- Prandtl numbers for molecular duffusion
      Pr_mol=4.*gam/(9.*gam-5.)
      eff  = eff0/float(nwaves)*(1.-exp(-float(nsec)/5./86400.))
      efheat=0.5
      fdec = EXP(-dz/h)
      kgitm1=kgit-1
      kgitm2=kgit-2
c------------------------------------------ COMMA to GW levels
      do k=1,kgitm1
        kp1=k+1 
        alfaNC(k)=(alfa  (k)+alfa  (kp1))/2.
        debg0 (k)=(geddy0(k)+geddy0(kp1))/2.
      end do
c-------------------------------------------------- LATITUDES
      do 1000 j=1,nb
        do nt=1,8
         w0(nt)    = ww(1,j) 
         w0(nt+ 8) = ww(2,j) 
         w0(nt+16) = ww(3,j) 
         w0(nt+24) = ww(4,j) 
         w0(nt+32) = ww(5,j) 
         w0(nt+40) = ww(6,j) 
        end do
c-------------------------------------------------- LONGITUDES
	do 1000 i= 1,igit
	   do k= 1,kgit
	      fluxU(k) = 0.                                               !by Rliu
	      fluxV(k) = 0.                                               !by Rliu
	      forceU(k)= 0.
	      forceV(k)= 0.
	      forceT(k)= 0.
	      geddy(k,i,j)=geddy0(k)
	      rheat=0.015/rm(k)*tn1(j,k,i)**0.667
	      rvisk=Pr_mol/C_p*rheat
	      geddym(k,i,j)=rvisk/rhos/rou(k)
	   end do
c-------------------------------         background profiles
	   do k=1,kgitm1
	      kp1=k+1
	      wzonk =un1(j,k,i)
	      wmerk =vn1(j,k,i)
	      tempk =tn1(j,k,i)
	      wzonk1=un1(j,kp1,i)
	      wmerk1=vn1(j,kp1,i)
	      tempk1=tn1(j,kp1,i)
	      u0 (k)=(wzonk1+wzonk)/2.
	      u0z(k)=(wzonk1-wzonk)/dz
	      v0 (k)=(wmerk1+wmerk)/2.
	      v0z(k)=(wmerk1-wmerk)/dz
	      t0 (k)=(tempk1+tempk)/2.
	      t0z(k)=(tempk1-tempk)/dz
	      omb2(k) = Rg/h*((gam-1.)/gam*t0(k)/h+t0z(k))
	      debgm(k)= (geddym(k,i,j)+geddym(kp1,i,j))/2.
	    enddo
	   do k=2,kgit-2
	      km1=k-1
	      kp1=k+1
	      u0zz(k)=(u0(km1)-2.*u0(k)+u0(kp1))/dz/dz
	      v0zz(k)=(v0(km1)-2.*v0(k)+v0(kp1))/dz/dz
	    end do
c-------------------------------        over all waves considered
	   do 2000 nw=1,nwaves
cfl	      print*, 'nw=',nw
	      do k= 1,kgit
		  degw(k)=0.
		  Fheat1(k)=0.
		  heatc1(k)=0.
		  Fheatc(k)=0.
		  heatco(k)=0.
	      end do
	      mc= int((nw-1)/7.9) + 1
	      mpha= nw - (mc-1)*8
	      cost= cos(pi*teta(mpha)/180.)
	      sint= sin(pi*teta(mpha)/180.)
c-----------------------------            critical level (kcrit) 
	      do k=4,kgitm1
		 if (omb2(k).le.0.) go to 10
c-----------------------------            intrinsic phase speed 
		 if(k.eq.4) cphkm1=cphase(mc)-u0(3)*cost-v0(3)*sint
		 cphint(k)=cphase(mc)-u0(k)*cost-v0(k)*sint
		 if(abs(cphint(k)).le.1.) go to 10 
		 sq= cphint(k)*cphkm1
C--------------------------------------------------------------
C    the perturbation approach is applicable if only the radiative 
C    damping is weak, so "alfaNC*SQRT(omb2)*h/rkx/cphint/cphint"
C    has to be small (<< 1) in other cases we assume that we are
C    close to the critical level  
C--------------------------------------------------------------
		 geddef=debg0(k)*(1.+1./pr)+debgm(k)*(1.+1./Pr_mol)
		 alfaEF=alfaNC(k)+geddef*omb2(k)/cphint(k)/cphint(k)
		 dcrit = alfaEF*SQRT(omb2(k))*h/rkx 
		 if (sq.le.dcrit) go to 10
		 cphkm1=cphint(k)
	      end do
   10         kcrit  = k-1
cfl	      if (kcrit.le.4) print*,'kcrit=',kcrit,'.le.4 - goto 2000'
cfl	      if (kcrit.gt.4) print*,'kcrit=',kcrit,'.gt.4 - ok!'
	      if (kcrit.le.4) go to 2000
	      kcrtp1=kcrit+1
	      kcrtp2=kcrit+2
	      wamkm1 = w0(nw)
C----------------  the vertical wave number Kz at the lower boundary
C                  (altitude about of 10 km)
	      rkz0   = SQRT(omb2(4)/cphint(4)/cphint(4))
	      rkzkm1 = rkz0
	      ompkm1 = rkx*cphint(4)
	      realL1 = -1./h
	      do 3000 k=5,kcrit
		  omb   = SQRT(omb2(k))
		  omega = rkx*cphase(mc)
		  ompls = rkx*cphint(k )
C--------------------------------------------------------------------
C                                   Eqn. [d^2/dz^2+L(z)*d/dz+M(z)]W=0
C--------------------------------------------------------------------
		  realL = -1./h
		  omplsz= -rkx*(u0z (k)*cost+v0z (k)*sint)
		  omplzz= -rkx*(u0zz(k)*cost+v0zz(k)*sint)
		  realM = rkx*rkx*omb2(k)/ompls/ompls-omplsz/ompls/h-
     &			  omplzz/ompls
C---------------------- W.K.B. solution -----------------------------
C      rkz2  = realM-realL*realL/4.-d(realL)/dz/2.  
C--------------------------------------------------------------------
		  rkz2  = realM-realL*realL/4.
cfl		  print*,'rkz2.gt.0? rkz2=',rkz2
		  if(rkz2.gt.0.) go to 50 
		  wampl  = wamkm1
		  rkz=rkzkm1
		  go to 40
   50 		  continue 
		  rkz   = sqrt(rkz2)
C-------------------  correction to the L(z) due to Newtonian cooling 
C                     and background diffusion      (has to be < 1/h) 
		  geddef=debg0(k)*(1.+1./pr)+debgm(k)*(1.+1./Pr_mol)
		  realL = realL+
     &	          	  rkz*(alfaNC(k)+geddef*rkz*rkz)/abs(ompls)
C------------------------- W.K.B. solution for amplitude !W'(z)!
		  wampl  = wamkm1*SQRT(rkzkm1/rkz)*EXP(-(realL+realL1)
     &		           *dz/4.)     
c-----------------GW breaking-------------------------------------------     
		  wbreak = abs(ompls)/rkz !amplitude wo GW bricht		  
cfl		  print*,'wampl.le.wbreak?',wampl,wbreak		  
		  if(wampl.le.wbreak) go to 40 
		    wampl = wbreak !wampl cannot exceed wbreak due to breaking
c--------------- Eddy diffusion coefficient due to breaking of GW 
cfl		    print*,'berechne jetzt rkzave. vorher:',rkzave
		    rkzave = (rkz  +rkzkm1)/2. !average vertical wave length between level k and k-1
cfl		    print*,'hinterher:',nw,i,j,k,rkz,rkzkm1,rkzave
		    omplav = (ompls+ompkm1)/2. !average intrinsic frequency between level k and k-1
C-------------------------------------------------- COMMA levels
		    geddef=geddy0(k)*(1.+1./pr)+geddym(k,i,j)*
     &		  	   (1.+1./Pr_mol)
		    degw(k) = abs(omplav)/(1.+1./pr)/rkzave/rkzave/
     &			      rkzave*(1./h-rkzave*(alfa(k)+geddef*
     &			      rkzave*rkzave)/abs(omplav)-2.*
     &			      ALOG(wampl/wamkm1*SQRT(rkz/rkzkm1))/dz)
		    if(degw(k).gt.1000.) degw(k)=1000.
   40 		  continue
		  !print*,'nach 40:',nw,i,j,k,rkz,rkzkm1,rkzave
		  rkzave = (rkz  +rkzkm1)/2. !FL hier gesetzt (Vermutung!)
cfl		  print*,'nach 40:',nw,i,j,k,rkz,rkzkm1,rkzave
C---------------------------------------------------------------
C     Force per unit mass due to divergency of  GW momentum flux
C     The simplified Eqn. U'=i/Kx*(d/dz-1/h)W' has been used,
C     0.5 is due to averaging over period ave.(W'W')= 0.5*!W'!^2
C--------------------------------------------------------------- 
C-------------------------------------------------- COMMA levels
		  flux0(k)=0.5*(rkz*wampl*wampl+rkzkm1*wamkm1*wamkm1)/
     &			   2./rkx        !by Rliu
		  ac0(k) = -0.5/rkx*
     &		          ((rkz*wampl*wampl-rkzkm1*wamkm1*wamkm1)/dz-
     &		          (rkz*wampl*wampl+rkzkm1*wamkm1*wamkm1)/2./h)
		  if(cphint(3).lt.0.) ac0(k)=-ac0(k)
C---------------------------------------------------------------
C               Heating/Cooling terms and flux of heat due to GW
C---------------------------------------------------------------
		  wampla = (wampl+wamkm1)/2.   
		  geddef = degw(k)+geddy0(k)+geddym(k,i,j) !degw(k)=0 initialisiert falls nicht unter 50 (oben) geändert
		  if(geddef.gt.1000.) geddef=1000.
		  !print*,nw,i,j,k,efheat,wampla,rkzave,C_p,rkx,geddef !abbruch: rkzave(nw=9,i=1,j=1,k=5)=e+35, kein abbruch: rkzave=-0.2, e-12, e-18 ???
		  heatc1(k) = efheat*wampla*wampla*rkzave*rkzave/2./
     &		              C_p/rkx/rkx*geddef*rkzave*rkzave
		  derkz2(k) = geddef*rkzave*rkzave
C---------------------------------------------------- GWs levels
		 Fheat1(k) = - wampl*wampl*h*omb2(k)/2./Rg/ompls/ompls
		  wamkm1 = wampl
		  rkzkm1 = rkz
		  ompkm1 = ompls
		  realL1 = realL
 3000 	      continue
	      do k=1,kgit-1
C--------------------------- COMMA to GWs levels
		  Fheat1(k) = Fheat1(k)*(alfaNC(k)+
     &			      (derkz2(k+1)+derkz2(k))/2./pr)
	      end do
	      do k=2,kgit
C--------------------------- GWs to COMMA levels
		heatc1(k) = heatc1(k)+(Fheat1(k)+Fheat1(k-1))/2./gam/h
	      end do
C------------------- Smoothing of heating/cooling and flux terms
	      if(kcrit.le.kgitm1) then
		flux0 (kcrtp1)=flux0 (kcrit)*fdec                               !by Rliu  
		ac0   (kcrtp1)=ac0   (kcrit)*fdec
		Fheat1(kcrtp1)=Fheat1(kcrit)*fdec
		heatc1(kcrtp1)=heatc1(kcrit)*fdec
		degw  (kcrtp1)=degw  (kcrit)*fdec
	      end if
	      if(kcrit.le.kgitm2) then
		flux0 (kcrtp2)=flux0 (kcrtp1)*fdec                              !by Rliu  
		ac0   (kcrtp2)=ac0   (kcrtp1)*fdec
		Fheat1(kcrtp2)=Fheat1(kcrtp1)*fdec
		heatc1(kcrtp2)=heatc1(kcrtp1)*fdec
		degw  (kcrtp2)=degw  (kcrtp1)*fdec
	      end if
	      do k=4,kgitm2
		km1=k-1
		km2=k-2
		kp1=k+1
		kp2=k+2
c      ac1(k)=ac0(k)-(ac0(km2)-4.*ac0(km1)+6.*ac0(k)+
c     $               ac0(kp2)-4.*ac0(kp1))/16. 
c      Fheatc(k)=Fheat1(k)-(Fheat1(km2)-4.*Fheat1(km1)+6.*Fheat1(k)+
c     $                     Fheat1(kp2)-4.*Fheat1(kp1))/16. 
c      heatco(k)=heatc1(k)-(heatc1(km2)-4.*heatc1(km1)+6.*heatc1(k)+
c     $                     heatc1(kp2)-4.*heatc1(kp1))/16. 
c      geddy(k,i,j)=geddy(k,i,j)+
c     $        eff*(degw(k)-(degw(km2)-4.*degw(km1)+6.*degw(k)+
c     $                      degw(kp2)-4.*degw(kp1))/16.) 
		flux1(k)=(flux0(km2)+flux0(km1)+flux0(k)+flux0(kp2)+      !by Rliu
     &			 flux0(kp1))/5.
		ac1(k)=(ac0(km2)+ac0(km1)+ac0(k)+ac0(kp2)+ac0(kp1))/5. 
		Fheatc(k)=(Fheat1(km2)+Fheat1(km1)+Fheat1(k)+
     &	                  Fheat1(kp2)+Fheat1(kp1))/5. 
		heatco(k)=(heatc1(km2)+heatc1(km1)+heatc1(k)+
     &	                  heatc1(kp2)+heatc1(kp1))/5. 
		geddy(k,i,j)=geddy(k,i,j)+
     &			     eff*(degw(km2)+degw(km1)+degw(k)+
     &	                     degw(kp2)+degw(kp1))/5. 
		if(geddy(k,i,j).gt.1000.) geddy(k,i,j)=1000.
	      end do
c      ac1(kgitm1)=ac0(kgitm1)+
c     $           (ac0(kgitm2)-2.*ac0(kgitm1)+ac0(kgit))/4. 
c      Fheatc(kgitm1)=Fheat1(kgitm1)+
c     $           (Fheat1(kgitm2)-2.*Fheat1(kgitm1)+Fheat1(kgit))/4. 
c      heatco(kgitm1)=heatc1(kgitm1)+
c     $           (heatc1(kgitm2)-2.*heatc1(kgitm1)+heatc1(kgit))/4. 
	    flux1(kgitm1)=(flux0(kgitm2)+flux0(kgitm1)+flux0(kgit))/3.        !by Rliu
	      ac1(kgitm1)=(ac0(kgitm2)+ac0(kgitm1)+ac0(kgit))/3. 
	      Fheatc(kgitm1)=(Fheat1(kgitm2)+Fheat1(kgitm1)+
     &			     Fheat1(kgit))/3. 
	      heatco(kgitm1)=(heatc1(kgitm2)+heatc1(kgitm1)+
     &			     heatc1(kgit))/3. 
	  do k=5,kcrit
	    km1=k-1
		
		heatGW = heatco(k)-(Fheatc(k)-Fheatc(km1))/dz
		fluxU(k) = fluxU(k) + eff*flux1(k)*cost                          !by Rliu
		fluxV(k) = fluxV(k) + eff*flux1(k)*sint                          !by Rliu
		forceU(k) = forceU(k) + eff*ac1(k)*cost
		forceV(k) = forceV(k) + eff*ac1(k)*sint
		forceT(k) = forceT(k) + eff*heatGW
	      end do
	      fluxU(kcrtp1)= fluxU(kcrit)*fdec                                 !by Rliu
	      fluxV(kcrtp1)= fluxV(kcrit)*fdec                                 !by Rliu
	      forceU(kcrtp1)= forceU(kcrit)*fdec
	      forceV(kcrtp1)= forceV(kcrit)*fdec
	      forceT(kcrtp1)= forceT(kcrit)*fdec
	      
 2000 	   continue
	   do k=1,kgit
	      fluxgru(k,i,j)= fluxU(k)                                         !by Rliu
	      fluxgrv(k,i,j)= fluxV(k)                                         !by Rliu
	      fgru(k,i,j)= forceU(k)
	      fgrv(k,i,j)= forceV(k)
	      fgrt(k,i,j)= forceT(k)
	      
	   end do
 1000 continue
 
      do j=1,nb
	do k=4,kgit
	  do i=1,igit
	      work(1,i)=geddy(k,i,j)
	      work(2,i)=fgru (k,i,j)
	      work(3,i)=fgrv (k,i,j)
	      work(4,i)=fgrt (k,i,j)
	      work(5,i)=fluxgru(k,i,j)		                                  !by Rliu
	      work(6,i)=fluxgrv(k,i,j)			                              !by Rliu
	      
	  end do
	  do i=1,igit
	      ip1=i+1
	      ip2=i+2
	      im1=i-1
	      im2=i-2
c---------------------------
	      if(i.eq.1) then
		im1=igit
		im2=igit-1
	      end if
	      if(i.eq.2) im2=igit
c---------------------------
	      if(i.eq.igit) then
		ip1=1
		ip2=2
	      end if
	      if(i.eq.igit-1) ip2=1
	      geddy(k,i,j) = (work(1,im2)+work(1,im1)+work(1,i)+
     &	                     work(1,ip2)+work(1,ip1))/5.  
	      fgru (k,i,j) = (work(2,im2)+work(2,im1)+work(2,i)+
     &	                     work(2,ip2)+work(2,ip1))/5.  
	      fgrv (k,i,j) = (work(3,im2)+work(3,im1)+work(3,i)+
     &	                     work(3,ip2)+work(3,ip1))/5.  
	      fgrt (k,i,j) = (work(4,im2)+work(4,im1)+work(4,i)+
     &	                     work(4,ip2)+work(4,ip1))/5.  
	      fluxgru(k,i,j) = (work(5,im2)+work(5,im1)+work(5,i)+
     &	                       work(5,ip2)+work(5,ip1))/5.                      !by Rliu
	      fluxgrv(k,i,j) = (work(6,im2)+work(6,im1)+work(6,i)+
     &	                       work(6,ip2)+work(6,ip1))/5.                      !by Rliu
     
          end do
        end do
      end do
      
      return
      
      end

