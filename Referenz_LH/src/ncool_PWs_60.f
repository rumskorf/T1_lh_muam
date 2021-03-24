c     Alex Pogoreltsev, February 2006
c
c     RETURNED TO THE NCEP/NCAR DATA, i.e. adjustment starts
c     at second level !!!!!!!!!! (see bellow)
c
c     convective adjustment at the level k=1 is included,
c     adjustment to the T(NCEP/NCAR) is removed
c     tau = 5 days and ttau is introduced
c-----------------------------------------------------------
      subroutine ncool_PWs_60
c 
      include 'com_main.fc'
      include 'com_namelist.fc' !FL
c
c      g=9.81
c      RgSI=8314.41
c      h=7000.
c      c_p=1005.
      rgit=float(igit)
      gradT0=-g/c_p
      tau=5.*86400.
      ttau = float(nsec)/tau
      eff_t  = 1.-exp(-ttau*ttau)
      
      do 10 k=1,kgit
      do 20 j=1,nb
         tav=0.
         do i=1,igit     
            tav=tav+tn1(j,k,i)/rgit !zonal mean temperature
         enddo
         do i=1,igit
c     -----
c     NEWTONIAN COOLING
c     -----
c-------------------------- adaptation to the NCEP/NCAR temperature
c
            fnt(k,i,j)=0.
	    if(k.le.4) fnt(k,i,j) = -eff_t/tau*(tav-T_NCEP(j,k))
c--------------------------------------------- Latent heating
	    if(k.le.10) fnt(k,i,j) = fnt(k,i,j) + eff_t*heat_LH(i,j,k)/86400.
         enddo 
   20 continue
   10 continue
   
      do i=1,igit     
         do k=2,kgit
             do j=1,nb
c     -----
c     CONVECTIVE ADJUSTMENT
c     -----
		if(k.eq.1) then
		  tav  =   (tn1(j,k,i)+T_1000(j,i))/2.
		  gradT=2.*(tn1(j,k,i)-T_1000(j,i))/dz
		else
		  km1=k-1 
		  tav  = (tn1(j,k,i)+tn1(j,km1,i))/2.
		  gradT= (tn1(j,k,i)-tn1(j,km1,i))/dz
		end if
		hscale=RgSI*tav/g/rm(k)
		gradTh=gradT0*hscale/h
		if(gradT.le.gradT0) then
		  conv_c=(gradTh-gradT)*dz/tau
		  fnt(k,i,j)=fnt(k,i,j)+eff_t*conv_c
		  if(k.ge.2) fnt(km1,i,j)=fnt(km1,i,j)-eff_t*conv_c
		end if
             enddo 
         end do
      enddo
      
      return
      end




