      SUBROUTINE EYGwave_RSHU(ii,jj,geddy,rho,vx1d,vy1d,temp1d,cp1d,h,
     &   ut_gwd, vt_gwd, pres, gwh_ir, gwh_dif, var_tot, flux_tot)
       IMPLICIT NONE
      integer ii,jj 
      integer, parameter  :: kgit=56
!AM Change this later!!!      
      REAL*8,  PARAMETER	:: RGAS=287.
      REAL*8,  PARAMETER  :: GRAV=9.81
      REAL*8,  PARAMETER  :: PPI=3.1415926
!AM end for future change

      INTEGER, PARAMETER  :: ht_dim= kgit    !Number of vertical levels
      INTEGER, PARAMETER  :: SLEV = 4        !Source level
      INTEGER, PARAMETER  :: nh = 30         !Number of harmonics
! Spectral parameters
!      REAL*8, PARAMETER   :: flux0    = 0.00026
      REAL*8, PARAMETER   ::omega0=1.e-4

      REAL*8, PARAMETER   ::alpha1=5./6.

      REAL*8, PARAMETER   ::betta=1.

      REAL*8, PARAMETER   ::gamma=1.5 

      REAL*8, PARAMETER   ::c1=10.

      REAL*8, PARAMETER   ::c2=4.
      REAL*8, PARAMETER   :: flux00    = 0.000005 ! 0.000015
      REAL*8, PARAMETER   :: cw       = 75.
      REAL*8, PARAMETER   :: max_cp_y = 120.0 ! 61.47 (nh=28) ! 104.12  !80.0 (nh=30)
      REAL*8, PARAMETER   :: kx       = 2.*PPI/300.e3
      REAL*8, PARAMETER   :: max_ht = 500.e3 !No calculations above this height
c-----------------------------------------------------------------------------
      REAL*8, PARAMETER :: dmax=300.
      REAL*8, PARAMETER :: d0  =100.
      REAL*8 zz,flux0,omega1,omrel,c1rel,c2rel,uw_Gavrilov,persek

! Input variables
      REAL*8, INTENT(IN)  :: vy1d(ht_dim)   !v-wind
      REAL*8, INTENT(IN)  :: vx1d(ht_dim)   !u-wind
      REAL*8, INTENT(IN)  :: temp1d(ht_dim) ! 
      REAL*8, INTENT(IN)  :: cp1d(ht_dim)   ! Cp
      REAL*8, INTENT(IN)  :: h(ht_dim)      !geopotential height (in m)
      REAL*8, INTENT(IN)  :: pres(ht_dim)   !pressure (in Pa)
      REAL*8, INTENT(IN)  :: rho(ht_dim)    !density
      REAL*4, INTENT(IN)  :: geddy(kgit,64,36)    !eddy diffusion
      REAL*8              :: v_eddy(ht_dim)=0.0  !eddy viscosity (=null)
      REAL*8              :: eden(ht_dim)=0.0 ! electron density
      
! Output variables
      REAL*8, INTENT(OUT) :: ut_gwd(ht_dim)
      REAL*8, INTENT(OUT) :: vt_gwd(ht_dim)

! GW heating terms
      REAL*8, INTENT(OUT) :: gwh_ir(ht_dim) !irreversible heating (K/s)
      REAL*8, INTENT(OUT) :: gwh_dif(ht_dim)!Differential heating/cooling (K/s)
      REAL*8 :: scht1d(ht_dim)              ! 1D scale height

      REAL*8 :: m_vis(ht_dim),              ! molecular kinematic viscosity
     &          var_tot(ht_dim),
     &          flux_tot(ht_dim)

! mid-point variables
      REAL*8 :: temp1d_m(ht_dim)
      REAL*8 :: u_source_m(ht_dim)
      REAL*8 :: pres_m(ht_dim)
      REAL*8 :: rho_m(ht_dim)
      REAL*8 :: h_m(ht_dim)
      REAL*8 :: eden_m(ht_dim)
      REAL*8 :: scht1d_m(ht_dim)

      REAL*8 :: drag(ht_dim,nh)

      REAL   :: phasespeed(nh)              ! Phase speed
      REAL*8 :: uw_mom(nh) = 0.
      REAL*8 :: brunt(ht_dim) = 0.02        ! Buoyancy frequency
      REAL*8 :: theta(ht_dim), thetap(ht_dim)
      INTEGER :: b_lev(nh), c_lev(nh)
      REAL*8 :: dz(ht_dim), rho_pr(ht_dim)

! Dissipation variables
      REAL*8 :: tau(ht_dim, nh)     
      REAL*8 :: beta(ht_dim, nh)           
      REAL*8 :: beta_dif(ht_dim,nh)  = 0.   
      REAL*8 :: beta_mol(ht_dim,nh)  = 0.   
      REAL*8 :: beta_nc(ht_dim,nh)   = 0.   
      REAL*8 :: beta_non(ht_dim,nh) ! = 0.   
      REAL*8 :: beta_cond(ht_dim,nh) = 0.  
      REAL*8 :: beta_eddy(ht_dim,nh) = 0.  
      REAL*8 :: beta_ion(ht_dim,nh)  = 0.  

      REAL*8 :: fac1, fac2  ! Factors to simplify the dissipation variable

      REAL*8 :: alpha(ht_dim)        = 0. ! Newtonian cooling coefficient
      REAL*8 :: vin(ht_dim)          = 0. ! ion-neutral collision frequency


      INTEGER :: i, j, k, n, nd, m

      REAL*8 :: c_int(ht_dim,nh) 
      REAL*8 :: sign_c(ht_dim,nh)
      REAL*8 :: flux(ht_dim,nh)
      REAL*8 :: up(ht_dim,nh)
      REAL*8 :: upSq(ht_dim,nh)

! Total horizontal wind variance at the source level
      REAL*8 :: sigmaSq_s=0., rms

      INTEGER :: sgn(nh)

! Nonlinearity parameters
      REAL*8 :: sigma(ht_dim,nh)=0.
      REAL*8 :: sigmaSq(ht_dim,nh)=0.
      REAL*8 :: alpha_ins(ht_dim,nh)=0.

      REAL*8, PARAMETER   :: S2  = 1.4142136d0 ! (2.**0.5
      REAL*8, PARAMETER   :: S2P = 2.5066283d0 ! (2.*PPI)**0.5

! Anisotropy variables
      REAL*8  :: u_source(ht_dim), gwd(ht_dim), 
     &           gwh(ht_dim), gwhd(ht_dim)
      REAL*8  :: xv, yv, H_tld, EPS 

      tau(:,:)     = 1.;    gwd(:)       = 0.;   gwh(:)       = 0.
      gwhd(:)      = 0.;    drag(:,:)    = 0.
      sigmaSq(:,:) = 0.;    beta(:,:)    = 0.;   flux_tot(:)  = 0.
      var_tot(:)   = 0.;    upSq(:,:)    = 0.;   up(:,:)      = 0.
      flux(:,:)    = 0.;    beta_non(:,:)= 0.

      flux0    = flux00*(1.+1./3.*tanh((float(jj)-18.5)/6.))    
! Project the upper level winds onto the source level winds
      PROJECTION : DO n = SLEV, ht_dim
       IF(n.eq.SLEV) THEN
          u_source(n) = SQRT(vy1d(n)*vy1d(n) + vx1d(n)*vx1d(n))
          IF(u_source(n).eq.0.) u_source(n) = 1.
          yv = vy1d(n)/u_source(n)
          xv = vx1d(n)/u_source(n)
       ENDIF
          u_source(n) = vy1d(n)*yv + vx1d(n)*xv
      ENDDO PROJECTION

      theta(:) = temp1d(:)*(pres(1)/pres(:))**(RGAS/cp1d(:))
      scht1d(:)= temp1d(:)*RGAS/GRAV
c-----------------------------------------------------------------------
      do n=1,ht_dim
         zz=h(n)/1000.
c.........................analytical profile of the eddy diffusion coefficient
c         IF(zz.GE.90.)v_eddy(n)=
c     $                   dmax    *EXP(-(zz-90.)/20.*(zz-90.)/20.)
c         IF(zz.LT.90.)v_eddy(n)=
c     $                  (dmax-d0)*EXP(-(zz-90.)/20.*(zz-90.)/20.)+
c     $                        d0 *EXP( (zz-90.)/20.)
       v_eddy(n)=geddy(n,ii,jj)
c       v_eddy(n)=0.
      end do
c------------------------------------------------------------------------
! Mid-point values of model variables
      MID_POINT : DO n = 1, ht_dim-1
         temp1d_m(n)   = 0.5*(temp1d(n)   + temp1d(n+1))
         u_source_m(n) = 0.5*(u_source(n) + u_source(n+1))
         pres_m(n)     = 0.5*(pres(n)     + pres(n+1))
         rho_m(n)      = 0.5*(rho(n)      + rho(n+1))
         h_m(n)        = 0.5*(h(n)        + h(n+1))
         eden_m(n)     = 0.5*(eden(n)     + eden(n+1))
         scht1d_m(n)   = 0.5*(scht1d(n)   + scht1d(n+1))
      ENDDO MID_POINT

! Brunt frequency and background parameters
      BACKGROUND1: DO n = 1, ht_dim-1
         dz(n)    = h(n+1) - h(n)                ! Midpoint value
         rho_pr(n)= (rho(n+1) - rho(n))/dz(n)    ! Midpoint value
!         m_vis(n) = 3.87e-5*temp1d_m(n)**0.69/rho_m(n)    
!           m_vis(n) = 3.563E-7*temp1d_m(n)**0.69/rho_m(n)!Sensitivity paper visc.
         m_vis(n) = 4.8E-7*temp1d_m(n)**0.667/rho_m(n)
!       m_vis(n) = 3.77E-7*temp1d_m(n)**0.667/pres_m(n)*temp1d_m(n)*RGAS
!    &                                       /H_tld/H_tld
!Mars    m_vis(n) = 3.128e-7*temp1d_m(n)**0.69/rho_m(n) +1.e-1!For CO2 atmos
!AM      m_vis(n) = 4.5E-05*(temp1d_m(n)/1000.0)**0.71/rho_m(n)!CMAT2 viscousity
!         vin(n)  = 7.22e-17*temp1d_m(n)**0.37*eden_m(n)
!Newtonian cooling
!        alpha(n) = 3.e-6*(1.2+TANH(-7.*LOG(pres_m(n)/100000.)-50.)*0.833)
!        alpha(n) = 3.e-6*(1.2+TANH(-7.*LOG(pres_m(n)/100000.)-50.))
!        alpha(n)=0. 
!       print*,real(h(n)), real(m_vis(n)),real(m_vis(n))
      ENDDO BACKGROUND1
       
!	       STOP

      BACKGROUND2: DO n=1, ht_dim-1
         thetap(n)= (theta(n+1) - theta(n))/dz(n) ! gradient of theta
         brunt(n) = (ABS((2.*GRAV/(theta(n)+theta(n+1)))*thetap(n)))
         brunt(n) = sqrt(brunt(n))
      ENDDO BACKGROUND2

!     WRITE(6,*) "EY Gravity Wave Spectrum"
!     WRITE(6,*) " i      Phasespeed      Flux         U_prime       Sigma_sq_tot    Hwhm"
!     WRITE(6,*) "          [ms-1]       [m2s-2]        [m/s]           [m2s-2]"
!     WRITE(6,*)
!      open(66,file='spectr.dat')

      SPECTRUM : DO i = 1, nh
       if(i.le.nh/2) then
!         persek=(40.+float(i-1)*10.)*60.
!         persek=(60.+float(i-1)*5.)*60.
!         phasespeed(i)=300.e3/persek
         phasespeed(i)=35.+float(i-1)*5.
         persek=300.e3/phasespeed(i)  
       else
         phasespeed(i)=-phasespeed(i-nh/2)
       end if

!         IF(i.eq.1)  phasespeed(i) = -max_cp_y  !80.
!         IF(i.GE.2)  phasespeed(i) = phasespeed(i-1)
!     &                  *((max_cp_y/2)**(-1./((nh/2.)-1.)))
!          
!         IF(i.eq.(nh/2+1))   phasespeed(i) = 2. 
!         IF(i.GE.(nh/2+2)) phasespeed(i) = phasespeed(i-1)
!     &                       *((max_cp_y/2)**(1./((nh/2.)-1.)))
!--------------------------------------------------------------
!          IF(i.eq.1)  phasespeed(i) = -max_cp_y  !80.
!          IF(i.GE.2)  phasespeed(i) = phasespeed(i-1)
!     &               + (max_cp_y-2.)/(nh/2-1)
!          IF(i.eq.(nh/2+1))   phasespeed(i) = 2.
!          IF(i.GE.(nh/2+1)) phasespeed(i) = phasespeed(i-1)
!     &               + (max_cp_y-2.)/(nh/2-1)
!      write(66,*) i,persek/60.,phasespeed(i)
      ENDDO SPECTRUM

!     phasespeed(:) = phasespeed(:) + u_source_m(SLEV)     !AM
c--------------------------------------------------------------


       SPECTRUM2: DO i = 1, nh
        omega1=abs(phasespeed(i))*kx

        omrel = omega1/omega0

        c1rel = abs(phasespeed(i))/c1

        c2rel = c2/abs(phasespeed(i))
        sgn(i)    = phasespeed(i)/ABS(phasespeed(i))
        uw_Gavrilov = sgn(i)*phasespeed(i)/
     $             (1.+omrel**alpha1)/(1.+c1rel**betta+c2rel**gamma)
!        uw_mom(i)=flux0*uw_Gavrilov*uw_Gavrilov
!      write(66,*) i,uw_mom(i)
! --- Aymmetric/Shifted spectrum ---------------------------------------------
!         sgn(i) = (phasespeed(i)-u_source_m(SLEV))
!     &            /(ABS(phasespeed(i)-u_source_m(SLEV)))  
!         uw_mom(i) = sgn(i)
!     &      *(flux0*EXP(-((phasespeed(i)-u_source_m(SLEV))/cw)**2))
! ----------------------------------------------------------------------------
! ---- Symmetric spectrum ---------------------------
        sgn(i)    = phasespeed(i)/ABS(phasespeed(i))
        uw_mom(i) = sgn(i)*(flux0*EXP(-(phasespeed(i)/cw)**2))
!      write(66,*) i,uw_mom(i)
!----------------------------------------------------
!Half spectrum:
!        IF(i.LE.nh/2) uw_mom(i)=0.
!-----------------------------------


         tau(SLEV,i)  = 1.
         flux(SLEV,i) = uw_mom(i)
         flux_tot(SLEV)=flux_tot(SLEV) + flux(SLEV,i)
         c_int(SLEV,i)  = phasespeed(i) - u_source_m(SLEV)
         IF(abs(c_int(SLEV,i)) .ge. brunt(SLEV)*100.) upSq(SLEV,i)     !AM
     &            = ABS(uw_mom(i))*brunt(SLEV)/(kx*ABS(c_int(SLEV,i)))
         IF(ABS(upSq(SLEV,i)).GT.0.) up(SLEV,i) = SQRT(upSq(SLEV,i))
         IF (c_int(SLEV,i).NE.0) sign_c(SLEV,i) = 
     &                          c_int(SLEV,i)/ABS(c_int(SLEV,i))
         IF (c_int(SLEV,i).EQ.0) sign_c(SLEV,i) = 1.

         var_tot(SLEV) = var_tot(SLEV) + upSq(SLEV,i)
      ENDDO SPECTRUM2

      c_lev(:) = -1
      b_lev(:) = -1

! Altitude loop starts
      HT_LOOP : DO n = SLEV+1, ht_dim-1
       ! phasespeed loop
        PS_LOOP : DO i = 1, nh
          
          IF(h(n) > max_ht) CYCLE
          
          c_int(n,i)  = phasespeed(i) - u_source_m(n)
          IF(c_int(n,i).NE.0) sign_c(n,i) = c_int(n,i)/ABS(c_int(n,i))
          IF(c_int(n,i).EQ.0) sign_c(n,i) = 1.
      
	!AM Added for trying only. Remove it!
		if( abs(c_int(n,i)) .le. 1.) 
     &         c_int(n,i) = 1.*(c_int(n,i)/abs(c_int(n,i)))
	!AM End of remove
          
          CRITICAL : IF(c_lev(i).LT.0 
     &                 .AND. (sign_c(n,i).NE.sign_c(n-1,i)
     &                 .AND. c_int(n,i).NE.0.)) THEN
 !AM         .AND. abs(c_int(n,i)).GE.brunt(n)*90.))  THEN    !AM
                       c_lev(i) = n
                       flux(n,i)  = 0.

           ELSE IF (c_lev(i).LT.0) THEN
	!AM Added for trying only. Remove it!
		if(abs(c_int(n,i)) .le. 1.) 
     &         c_int(n,i) = 1.*(c_int(n,i)/abs(c_int(n,i)))
	!AM End of remove
             fac1      = 2.*brunt(n)*brunt(n)*brunt(n)
     &                  /(kx*c_int(n,i)**4)
             fac2      = 2.*brunt(n)/(kx*c_int(n,i)*c_int(n,i))
!################################################################
       EPS=2.*(m_vis(n)+v_eddy(n))*
     &         brunt(n)*brunt(n)/kx/abs(c_int(n,i))**3
       IF(EPS.ge.1.) fac1 = fac1/EPS/EPS
!       if(i.eq.nh.and.jj.eq.6.and.ii.eq.1) then
!        print *, EPS,brunt(n),c_int(n,i),rho(n)
!        If(n.eq.ht_dim-1) stop
!       end if
!################################################################
             beta(n,i) = fac1*(m_vis(n)+v_eddy(n)) 
     &                  + fac2*(vin(n)+alpha(n))
!             beta(n,i) = 0.

             NON_LINEARITY : DO j = 1, nh
                c_int(n,j) = phasespeed(j) - u_source_m(n)

	           !AM Added for trying only. Remove it!
			    if(abs(c_int(n,j)) .le. 1.) 
     &              c_int(n,j) = 1.*(c_int(n,j)/abs(c_int(n,j)))
	           !AM End of remove
                
			  IF (ABS(c_int(n,i)) .GE. ABS(c_int(n,j))) THEN
                   sigmaSq(n,i) = sigmaSq(n,i) + upSq(n-1,j)
                ENDIF
             ENDDO NON_LINEARITY
                
             IF(sigmaSq(n,i).ge.1e-36) THEN
                sigma(n,i)     = SQRT(sigmaSq(n,i))
                alpha_ins(n,i) = ABS(c_int(n,i))/S2/sigma(n,i)
                IF(alpha_ins(n,i).GE.1.e10) THEN
                   beta_non(n,i) = 0.
                ELSE
                   beta_non(n,i)  = S2P*brunt(n)/sigma(n,i) 
     &                            *EXP(-alpha_ins(n,i)*alpha_ins(n,i))
!                   beta_non(n,i)  = 0.
                ENDIF
             ENDIF
!##################################################################
!       if(i.eq.nh/2+8.and.jj.eq.6.and.ii.eq.1) then
!        print *, v_eddy(n), h(n), beta(n,i), beta_non(n,i)
!        If(n.eq.ht_dim-1) stop
!       end if
!##################################################################              
             beta(n,i) = beta(n,i) + beta_non(n,i)
             tau(n,i)  = tau(n-1,i)
     &                    *EXP(-dz(n)*(beta(n,i)+beta(n-1,i))*0.5)
             flux(n,i) = uw_mom(i)*rho(SLEV)/rho(n)*tau(n,i)
             upSq(n,i) = ABS(flux(n,i))*brunt(n)/kx/ABS(c_int(n,i))
             IF(ABS(upSq(n,i)).GT.0.) up(n,i) = SQRT(upSq(n,i))
                
             drag(n,i) = beta(n,i)*flux(n,i)

             IF(alpha_ins(n,i).LT.0.75) drag(n,i) = 0.

             flux_tot(n) = flux_tot(n) + flux(n,i) 
             var_tot(n)  = var_tot(n)  + upSq(n,i)
             gwd(n)    = gwd(n) + drag(n,i)
             gwh(n)    = gwh(n) + drag(n,i)*c_int(n,i)/cp1d(n)
             gwhd(n)   = gwhd(n) + ((scht1d_m(n)/2./RGAS/rho_m(n))*!EY 10.10.08 
     &            ( rho_pr(n)*c_int(n,i)*drag(n,i)  + 
     &            rho_m(n)*((c_int(n,i)-c_int(n-1,i))/dz(n))*drag(n,i)+ 
     &            rho_m(n)*c_int(n,i)*((drag(n,i)-drag(n-1,i))/dz(n))))   
          ENDIF CRITICAL
        ENDDO PS_LOOP
      ENDDO HT_LOOP
    
! Back-interpolate the drag onto the model full-levels.
      DO n=2, ht_dim-1
         gwd(n) = 0.5*(gwd(n) + gwd(n-1))
         gwh(n) = 0.5*(gwh(n) + gwh(n-1))
      ENDDO

      SMOOTH : DO n = 2, ht_dim-1
        ! Do 2 delta smoothing of drag at full-levels      
        gwd(n)  = (gwd(n-1)  + 2.*gwd(n)  + gwd(n+1))*0.25
        gwh(n)  = (gwh(n-1)  + 2.*gwh(n)  + gwh(n+1))*0.25
        gwhd(n) = (gwhd(n-1) + 2.*gwhd(n) + gwhd(n+1))*0.25
      ENDDO SMOOTH
      gwd(ht_dim) = gwd(ht_dim-1)
    
      BACK_PROJECT : DO n=SLEV, ht_dim
        ut_gwd(n) = xv * gwd(n) 
        vt_gwd(n) = yv * gwd(n)
        gwh_ir(n) = gwh(n)
        gwh_dif(n)= gwhd(n)
      ENDDO BACK_PROJECT
   
      RETURN
      END SUBROUTINE EYGwave_RSHU




