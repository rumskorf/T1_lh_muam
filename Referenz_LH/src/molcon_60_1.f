c----------------------------------- molecular and eddy duffusion terms  
      subroutine molcon_60_1(i,j)
c----------------------------------------------------------------------
c     berechnung von : frunp1 : innere reibung beim zonalwind
c                      frvnp1 : innere reibung beim meridionalwind
c                      fmcnp1 : molekulare waermeleitung
c----------------------------------------------------------------------
      include 'com_main.fc'
      include 'com_namelist.fc' !FL
c
      real rho (kgit),cp_z(kgit) 
         rgit=real(igit)
         rb=real(nb)
c------- Pr_eddy is the turbulent Prandtl number (Huang+Smith, 1991)
         Pr_eddy=3.
c------------------------------------------------- just for information  
c        Critical Richardson number = 0.2-0.6 (Izakov, 1978)
c        0.25 (Ebel, et al., 1983)
c         are given in "init_60"
c         Rf_c=0.4
c         RgSI=8314.41
c         g=9.81 
c         h=7000.
c----------------------------------------------------------------------
         t0=273. 
         dz5  =.5/dz
         eff_u=1.
         do k=1,kgit
         fmcnp1(k)=0.
         frunp1(k)=0.
         frvnp1(k)=0.
         cp_z(k)=RgSI/rm(k)*(7./2.*(1.-roxvmr(k))+5./2.*roxvmr(k))
         end do
c-------------------------------------------------- initialization TIME
         t_mol= float(nsec)/86400.
         eff_mol= 1.-0.9*exp(-t_mol*t_mol)
c---------------------------------------------------------------------------
c   Molecular+Eddy (if mturb=1) VISCOUSITY rx(k) & Thermal CONDUCTION rl(k)
c   The molecular thermal conduction empirical coefficient, expressions, and 
c   molecular Prandtl number are taken from
c   Forbes and Garrett, Rev. Geophys. Space Phys., 1979, p. 1951-1981.  
c---------------------------------------------------------------------------
      do 20 k=1,kgit
         rho(k)=rhos*rou(k)*rm(k)/rm(1)*t0/tn1(j,k,i)
c------------------- Prandtl numbers for molecular duffusion
         gam_z=cp_z(k)/(cp_z(k)-RgSI/rm(k)) 
         Pr_mol=4.*gam_z/(9.*gam_z-5.)
c         if (ncom.ge.34.and.k.ge.48.and.j.eq.1) then
c         	print*,i,j,k,rm(k),tn1(j,k,i)
c         endif	
         rl(k)=0.015/rm(k)*tn1(j,k,i)**0.667
         rx(k)=Pr_mol/cp_z(k)*rl(k)
         hscale=RgSI*tn1(j,k,i)/g/rm(k)
         veddy(k)=geddy(k,i,j)
         teddy(k)=cp_z(k)*veddy(k)/Pr_eddy
c         rx(k)=(rx(k)+veddy(k)*rho(k))*h/hscale
c         rl(k)=(rl(k)+teddy(k)*rho(k))*h/hscale
         rx(k)=(eff_mol*rx(k)+veddy(k)*rho(k))*h/hscale
         rl(k)=(eff_mol*rl(k)+teddy(k)*rho(k))*h/hscale
   20 continue
         kgp1=kgit+1 
         kgm1=kgit-1
         rl (kgp1)=(8.*rl(kgit)-3.*rl(kgm1))/5.
         rx (kgp1)=(8.*rx(kgit)-3.*rx(kgm1))/5.
         do 10 k=4,kgit-2
           km1=k-1
           kp1=k+1
           km2=k-2
           kp2=k+2
           hscale=RgSI*tn1(j,k,i)/g/rm(k)
           det=rho(k)*hscale/h
c
            duz0=(un1(j,k,i)-un1(j,km2,i))*dz5
            dvz0=(vn1(j,k,i)-vn1(j,km2,i))*dz5
            dtz0=(tn1(j,k,i)-tn1(j,km2,i))*dz5
c
            duz2=(un1(j,kp2,i)-un1(j,k,i))*dz5
            dvz2=(vn1(j,kp2,i)-vn1(j,k,i))*dz5
            dtz2=(tn1(j,kp2,i)-tn1(j,k,i))*dz5
c
            duz =(un1(j,kp1,i)-un1(j,km1,i))*dz5
            dvz =(vn1(j,kp1,i)-vn1(j,km1,i))*dz5
            dtz =(tn1(j,kp1,i)-tn1(j,km1,i))*dz5
c
c--------------- additional terms (eddy diff.) only up to kgit-2
c
      fmcnp1(k)=(rl(kp1)*dtz2-rl(km1)*dtz0+
     $          g*(rho(kp1)*teddy(kp1)/cp_z(kp1)-
     $             rho(km1)*teddy(km1)/cp_z(km1)))*
     $                                            dz5/det/cp_z(k)
c----------------------------------------------------------------
c                      additional terms due to eddy diffusion and
c                    conversion of the mechanical energy to heat,
c                      in general, parameter eff_u has to be < 1. 
c                                    Alex Pogoreltsev, April 2002
c----------------------------------------------------------------  
       cp_z2 = cp_z(k)*cp_z(k)   
       fmcnp1(k)=fmcnp1(k)+
     $  g*teddy(k)*(dtz*h/hscale+g/cp_z(k))/cp_z2/tn1(j,k,i)/Rf_c+
     $  eff_u*(rx(k)*h/hscale*(duz*duz+dvz*dvz)/rho(k)+
     $               drag(j,k,1)*un1(j,k,i)*un1(j,k,i)+
     $               drag(j,k,2)*vn1(j,k,i)*vn1(j,k,i))/cp_z(k)
c----------------------------------------------------------------
           frunp1(k)=(rx(kp1)*duz2-rx(km1)*duz0)*dz5/det
           frvnp1(k)=(rx(kp1)*dvz2-rx(km1)*dvz0)*dz5/det
   10 continue
c------------------------------------------------------- k=kgit-1
      kgm2=kgit-2
      kgm3=kgit-3
        hscale=RgSI*tn1(j,kgm1,i)/g/rm(kgm1)
        det=rho(kgm1)*hscale/h
      duz0=(un1(j,kgm1,i)-un1(j,kgm3,i))*dz5
      dvz0=(vn1(j,kgm1,i)-vn1(j,kgm3,i))*dz5
      dtz0=(tn1(j,kgm1,i)-tn1(j,kgm3,i))*dz5
c---------------- assumption !!!! dT/dz(kgit)=1./3.*dT/dz(kgit-1)
c      duz2=1./3.*(un1(j,kgit,i)-un1(j,kgm2,i))*dz5
c      dvz2=1./3.*(vn1(j,kgit,i)-vn1(j,kgm2,i))*dz5
c      dtz2=1./3.*(tn1(j,kgit,i)-tn1(j,kgm2,i))*dz5
      duz2=1./4.*(7.*un1(j,kgit,i)-8.*un1(j,kgm1,i)+un1(j,kgm2,i))*dz5
      dvz2=1./4.*(7.*vn1(j,kgit,i)-8.*vn1(j,kgm1,i)+vn1(j,kgm2,i))*dz5
      dtz2=1./4.*(7.*tn1(j,kgit,i)-8.*tn1(j,kgm1,i)+tn1(j,kgm2,i))*dz5
        duz =(un1(j,kgit,i)-un1(j,kgm2,i))*dz5
        dvz =(vn1(j,kgit,i)-vn1(j,kgm2,i))*dz5
c
      fmcnp1(kgm1)=(rl(kgit)*dtz2-rl(kgm2)*dtz0)*dz5/det/cp_z(kgm1)
c
      fmcnp1(kgm1)=fmcnp1(kgm1)+
     $   eff_u*(rx(kgm1)*h/hscale*(duz*duz+dvz*dvz)/rho(kgm1)+
     $   drag(j,kgm1,1)*un1(j,kgm1,i)*un1(j,kgm1,i)+
     $   drag(j,kgm1,2)*vn1(j,kgm1,i)*vn1(j,kgm1,i))/cp_z(kgm1)
c
      frunp1(kgm1)=(rx(kgit)*duz2-rx(kgm2)*duz0)*dz5/det
      frvnp1(kgm1)=(rx(kgit)*dvz2-rx(kgm2)*dvz0)*dz5/det
c---------------------------------------------------------- k=kgit
        hscale=RgSI*tn1(j,kgit,i)/g/rm(kgit)
        det=rho(kgit)*hscale/h
      duz0=(un1(j,kgit,i)-un1(j,kgm2,i))*dz5
      dvz0=(vn1(j,kgit,i)-vn1(j,kgm2,i))*dz5
      dtz0=(tn1(j,kgit,i)-tn1(j,kgm2,i))*dz5
c---------------------------------------------------- dT/dz(kgit+1)=0
      duz2=0.
      dvz2=0.
      dtz2=0.
c        duz=1./3.*(un1(j,kgit,i)-un1(j,kgm2,i))*dz5
c        dvz=1./3.*(vn1(j,kgit,i)-vn1(j,kgm2,i))*dz5
      duz=1./4.*(7.*un1(j,kgit,i)-8.*un1(j,kgm1,i)+un1(j,kgm2,i))*dz5
      dvz=1./4.*(7.*vn1(j,kgit,i)-8.*vn1(j,kgm1,i)+vn1(j,kgm2,i))*dz5
c--------------------------------------------------------------------
      fmcnp1(kgit)=(rl(kgp1)*dtz2-rl(kgm1)*dtz0)*dz5/det/cp_z(kgit)
c
      fmcnp1(kgit)=fmcnp1(kgit)+
     $      eff_u*(rx(kgit)*h/hscale*(duz*duz+dvz*dvz)/rho(kgit)+
     $      drag(j,kgit,1)*un1(j,kgit,i)*un1(j,kgit,i)+
     $      drag(j,kgit,2)*vn1(j,kgit,i)*vn1(j,kgit,i))/cp_z(kgit)
c
      frunp1(kgit)=(rx(kgp1)*duz2-rx(kgm1)*duz0)*dz5/det
      frvnp1(kgit)=(rx(kgp1)*dvz2-rx(kgm1)*dvz0)*dz5/det
c---------------------------------------------------- spin up interval
        ttau = float(nsec)/10./86400.
c      if(ncom.le.nphi) then
c------------------------------------- strong damping at the beginning
        t_diss = exp(-ttau*ttau)
c      else
c----------------------------------------------------------------------
c        nfor=max0(ncom-nphi,0) 
c        xsec = 3600./float(ntime)          ! seconds per step
c        td1  = float(nfor)*xsec/ 1./86400. 
c        td2  = float(nfor)*xsec/30./86400. 
c------------------------------- damping during establishment of tides
c        t_diss = exp(-td2*td2) - exp(-td1*td1) 
c        t_diss = exp(-td2*td2) 
c      end if
      do 30 k=1,kgit
       alfaNC = 1./86400.*(0.5+0.5*tanh((z(k)-100.)/20.))*t_diss 
       betaRF = 1./86400.*(0.5+0.5*tanh((z(k)-100.)/20.))*t_diss 
       fmcnp1(k)=fmcnp1(k)-alfaNC*(tn1(j,k,i)-tz(k))
       frunp1(k)=frunp1(k)-betaRF*un1(j,k,i)
       frvnp1(k)=frvnp1(k)-betaRF*vn1(j,k,i)
   30 continue
      return
      end




