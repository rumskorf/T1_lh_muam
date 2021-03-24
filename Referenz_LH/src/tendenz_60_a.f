      subroutine tendenz_60_a
      include 'com_main.fc'
      include 'com_namelist.fc' !FL
      dimension aux(kgit),avx(kgit),atx(kgit),                  !
     &          auy(kgit),avy(kgit),aty(kgit),
     &          auz(kgit),avz(kgit),atz(kgit),ft1(kgit),
     &          fu1(kgit),fv1(kgit),fu2(kgit),fv2(kgit),
     &          au(kgit),av(kgit),at(kgit)
      dimension xm(kgit)
      real un1_00,vn1_00,tn1_00,fu(kgit),fv(kgit),ft(kgit),
c     $     frun(kgit),frvn(kgit),fmcn(kgit),
     $     un1pl,un1mn,vn1pl,vn1mn,tn1pl,tn1mn,wn1pl,wn1mn
      real fiave(kgit,nb),fiprt(kgit,igit,nb),fv2ave(kgit),fiplat
      real un1_w(kgit),vn1_w(kgit),tn1_w(kgit)
      real un1_0(nb,igit),vn1_0(nb,igit),tn1_0(nb,igit),fin1_0(nb,igit)    
      
c----------------------------------------------------------------------
      dt2=dt/2.
      rgit  =float(igit)
      kgm1=kgit-1
      Rh    = 287.04/h      ! R/M/H
      do 333 istep=1,2
c-----------------------     GEOPOTENTIAL (from hydrostatic equation)
      call geopot_60 
c-----------------------     VERTICAL WIND (from continuity equation)
      call vbew
c---------------     NC and heating due to PWs and tidal oscillations
      call ncool_PWs_60
c
c-------------------     smoothing, bi_harmonic diffusion of the momentum
c                                   and temperature
      call bi_diff_60
c--------------------------------------------------   advectiv terms
      do 11 j=1,nb
        if(j.gt.2.and.j.lt.nb-1) go to 11
c---- zonally averaged geopotential and perturbations near the poles
        do k=1,kgit 
           fiave(k,j)=0.
           do i=1,igit     
            fiave(k,j)   = fiave(k,j)+fin1(j,k,i)/rgit
           enddo
           do i=1,igit
            fiprt(k,i,j) = fin1(j,k,i)-fiave(k,j)
           enddo
        end do
   11 continue 
c-----------------------------------------------LATITUDE
      do 1 j=1,nb
         jplus1=j+1           
         jm1 = max(1 ,j - 1)
         jp1 = min(nb,jplus1)
         dey = 1./dy/cosfi(jplus1)
        do 12 k=1,kgit
        if(j.gt.1.and.j.lt.nb) go to 12
c-------- zonally averaged latitudinal gradient of geopotential
c         tends to zero on the poles 
c         PHIave_varphi(1 )=PHIave_varphi(1 +1/2)/2.
c         PHIave_varphi(nb)=PHIave_varphi(nb-1/2)/2.
c--------------------------------------------------------------
          fv2ave(k) = -(fiave(k,jp1) - fiave(k,jm1))/dy/2.
c-----------------------------------------------------------------------
c      see bellow - leads to the temperature oscillations near the poles
c
c      if(j.eq. 1) fv2ave(k) = -1./3.*(fiave(k, 3)-fiave(k,   1))/dy/2.
c      if(j.eq.nb) fv2ave(k) = -1./3.*(fiave(k,nb)-fiave(k,nb-2))/dy/2.
c-----------------------------------------------------------------------
   12 continue
c-----------------------------------------------LONGITUDE
      do 2 i=1,igit
	    i2=i+1
	    i0=i-1
	    if(i.eq.1   ) i0=igit
	    if(i.eq.igit) i2=1
        do k=1,kgit
        fu2(k) = -(fin1(j,k,i2) - fin1(j,k,i0))*0.5/dx(jplus1)
      if(j.gt.1.and.j.lt.nb) then
        fv2(k) = -(fin1(jp1,k,i) - fin1(jm1,k,i))/dy/2.
      else
c-------------------- geopotential gradient terms for perturbations
c         (perturbations of the geopotential are zero on the poles),
c          i.e.,
c          PHI_varphi(j=1 ) = (PHI(1   )+PHI(2 ))/2/dy
c          PHI_varphi(j=nb) =-(PHI(nb-1)+PHI(nb))/2/dy
c-----------------------------------------------------------------
         if(j.eq.1 ) fiplat= (fiprt(k,i,2   )+fiprt(k,i,1 ))/2./dy  
         if(j.eq.nb) fiplat=-(fiprt(k,i,nb-1)+fiprt(k,i,nb))/2./dy
        fv2(k)=fv2ave(k)-fiplat
      end if
        end do
      do 10 k=1,kgit
      xm(k)=rm(k)/rm(1) 
      un1_00 = un1(j,k,i)
      vn1_00 = vn1(j,k,i)
      tn1_00 = tn1(j,k,i)
c----------------------------------------     x-anteil/ro
      un1pl = (un1(j,k,i2)+un1_00)/2.
      un1mn = (un1(j,k,i0)+un1_00)/2.
      aux(k)= (un1pl*un1pl-un1mn*un1mn)/dx(jplus1)
c
      vn1pl = (vn1(j,k,i2)+vn1_00)/2.
      vn1mn = (vn1(j,k,i0)+vn1_00)/2.
      avx(k)= (un1pl*vn1pl-un1mn*vn1mn)/dx(jplus1)
c
      tn1pl = (tn1(j,k,i2)+tn1_00)/2.
      tn1mn = (tn1(j,k,i0)+tn1_00)/2.
      atx(k)= (un1pl*tn1pl-un1mn*tn1mn)/dx(jplus1)
c----------------------------------------     y-anteil/ro
      un1pl = (un1(jp1,k,i)+un1_00)/2.
      un1mn = (un1(jm1,k,i)+un1_00)/2.
      vn1pl = (vn1(jp1,k,i)+vn1_00)/2.
      vn1mn = (vn1(jm1,k,i)+vn1_00)/2.
      auy(k)= (un1pl*vn1pl*cosf2(jplus1)-un1mn*vn1mn*cosf2(j))*dey 
      avy(k)= (vn1pl*vn1pl*cosf2(jplus1)-vn1mn*vn1mn*cosf2(j))*dey 
c
      tn1pl = (tn1(jp1,k,i)+tn1_00)/2.
      tn1mn = (tn1(jm1,k,i)+tn1_00)/2.
      aty(k)= (vn1pl*tn1pl*cosf2(jplus1)-vn1mn*tn1mn*cosf2(j))*dey 
10    continue
c------------------------     z-anteil fuer die flaechen 'im' modell
      do 20 k=1,kgm1
      un1_00 = un1(j,k,i)
      vn1_00 = vn1(j,k,i)
      tn1_00 = tn1(j,k,i)
      kp1=k+1
      km1=k-1
      un1pl = (un1(j,kp1,i)+un1_00)/2.
      vn1pl = (vn1(j,kp1,i)+vn1_00)/2.
      tn1pl = (tn1(j,kp1,i)+tn1_00)/2.
      wn1pl =  wn1(j,kp1,i)
      wn1mn =  wn1(j,k  ,i)
      if(k.eq.1) then
c       tn1mn = (2.*tn1_00+0.64e-2*dz)/2.
c---------------------------------------- interpolation to the ground
c                                         A.Pogoreltsev, February 2006 
       un1mn = (3.*un1(j,1,i)-un1(j,2,i))/2.
       vn1mn = (3.*vn1(j,1,i)-vn1(j,2,i))/2.
       tn1mn = (3.*tn1(j,1,i)-tn1(j,2,i))/2.
      else
       un1mn = (un1(j,km1,i)+un1_00)/2.
       vn1mn = (vn1(j,km1,i)+vn1_00)/2.
       tn1mn = (tn1(j,km1,i)+tn1_00)/2.
      end if
      auz(k)= (un1pl*wn1pl-un1mn*wn1mn)/dz
      avz(k)= (vn1pl*wn1pl-vn1mn*wn1mn)/dz
      atz(k)= (tn1pl*wn1pl-tn1mn*wn1mn)/dz
c----------------------------------------------------- original COMMA
c      ft1(k) =-.28577*(wn1pl*tn1pl+wn1mn*tn1mn)/2./rou(k)/xm(k)/h
c
c------------------- has been changed by A.Pogoreltsev, February 2006
c see also Eqn. (3.6) in Dissertation Fröhlich
c		R/c_p		        w		T	m'  H 
c		  !			!		!	!   !
      ft1(k) =-.28577*(wn1pl/row(kp1)+wn1mn/row(k))/2.*tn1_00/xm(k)/h
 20   continue

c----------------------------------------------     z-terms at k=kgit
      un1pl = un1(j,kgit,i)
      vn1pl = vn1(j,kgit,i)
      tn1pl = tn1(j,kgit,i)
      wn1pl = wn1(j,kgit+1,i)
        un1mn = (un1(j,kgit,i)+un1(j,kgm1,i))/2.
        vn1mn = (vn1(j,kgit,i)+vn1(j,kgm1,i))/2.
        tn1mn = (tn1(j,kgit,i)+tn1(j,kgm1,i))/2.
        wn1mn =  wn1(j,kgit,i)

c        auz(kgit) = -un1mn*wn1mn/dz
c        avz(kgit) = -vn1mn*wn1mn/dz
c        atz(kgit) = -tn1mn*wn1mn/dz
      auz(k)= (un1pl*wn1pl-un1mn*wn1mn)/dz
      avz(k)= (vn1pl*wn1pl-vn1mn*wn1mn)/dz
      atz(k)= (tn1pl*wn1pl-tn1mn*wn1mn)/dz
c
c------------------- has been changed by A.Pogoreltsev, February 2006
c
      ft1(kgit) = -.28577*tn1(j,kgit,i)*wn1mn/2./row(kgit)/xm(kgit)/h
c---------------------------------------------------------------------
      do 30 k=1,kgit
      au(k)=aux(k)+auy(k)+auz(k)/rou(k)
      av(k)=avx(k)+avy(k)+avz(k)/rou(k)
      at(k)=atx(k)+aty(k)+atz(k)/rou(k)
c
c------ the correction to the Coriolis terms due to Hall conductivity
c------ Alex Pogoreltsev, March 2002
c
      un1_00 = un1(j,k,i)
      vn1_00 = vn1(j,k,i)
      
c      	         Coriolis force   Lorentz deflection   tan(phi)/a
c			!		!               !
      fu1(k)= vn1_00*(cor(jplus1)-drag(j,k,3)+un1_00*tgfia(jplus1))
      fv1(k)=-un1_00*(cor(jplus1)-drag(j,k,3)+un1_00*tgfia(jplus1))
 30   continue
c--------------------------------------------------------------------
      do 50 k=1,kgit

c-----------------------------------------------------------------------
c				Tendenzterme      
c-----------------------------------------------------------------------
      
c     NETTO ADVEC CORIOL PRESS        bi_DIFF        GWs
c        !    !     !      !           !             !     
      fu(k)=-au(k)+fu1(k)+fu2(k)    +fut(k,i,j)+fgru(k,i,j)   
      fv(k)=-av(k)+fv1(k)+fv2(k)    +fvt(k,i,j)+fgrv(k,i,j)
c
c             ADVEC ADIAB  New.cool     bi_diff     abs.sol.rad. (strobel)
c              !     !      !           !           !
      ft(k)=-at(k)+ft1(k)+fnt(k,i,j)+ftt(k,i,j)+hs  (k,i,j)+
     $                               fc (k,i,j)+fgrt(k,i,j)
c				        !           !
c				      ircool       GWs      
      
c--------------------------------------------------------- first step
c
      if(istep.eq.1) then
       un0(j,k,i)=un1(j,k,i)  
       vn0(j,k,i)=vn1(j,k,i)  
       tn0(j,k,i)=tn1(j,k,i)
      end if 
c					     ion drag
c						!
       un2(j,k,i)=un0(j,k,i) + dt*(fu(k)-drag(j,k,1)*un1(j,k,i))
       vn2(j,k,i)=vn0(j,k,i) + dt*(fv(k)-drag(j,k,2)*vn1(j,k,i))
       
c       if(i.eq.1.and.j.eq.1.and.k.ge.48) then
c       	  print*,'vor tn2 neu',k,i,j,tn2(j,k,i),tn0(j,k,i),dt,ft(k),h_PWs(k,i,j)       
c       endif
                                                                  
c					Plan. Waves
c					   !	
       tn2(j,k,i)=tn0(j,k,i) + dt*(ft(k)+h_PWs(k,i,j))       
       
c       if(i.eq.1.and.j.eq.1.and.k.ge.48) print*,k,i,j,'nach tn2 neu',tn2(j,k,i)
       
!        open(557,file='testreihe2.txt')
!        if (istep.eq.1.and.j.eq.7.and.k.eq.43.and.i.eq.26) then
! 	  write(557,'(5I5,1x,100(E10.3,1x))') ncom,ncom/nstep,j,k,i,     
!      +	  un1(j,k,i),vn1(j,k,i),tn1(j,k,i), 
!      +    fu(k),fv(k),ft(k), !netto
!      +	  drag(j,k,1),drag(j,k,2),drag(j,k,3), !ion drag and Lorentz deflection
!      +	  h_PWs(k,i,j), !daempfung/ planetary waves
!      +    au(k),av(k),at(k), !advection
!      +    fu1(k),fv1(k),ft1(k), !coriolis
!      +    fu2(k),fv2(k),fnt(k,i,j), !pressure/ newton cooling 
!      +    fut(k,i,j),fvt(k,i,j),ftt(k,i,j), !bi_DIFF
!      +	  fgru(k,i,j),fgrv(k,i,j),fgrt(k,i,j),hs(k,i,j),fc(k,i,j)    
!       endif
   50 continue
  2   continue !i-schleife
  1   continue !j-schleife
c-----------------------------------------------------------------
      do j=1,nb
        do i=1,igit
        if(ncom.le.nphi) philb1(j,i)=0.
c-----------------------------------------------------------------
c        to simulate with old lower boundary condition
c        if(ncom.le.nend) philb1(j,i)=0.
c----------------------------------------- U,V,T,PHI on the ground
        un1_0 (j,i)=(3.*un1(j,1,i)-un1(j,2,i))/2.
        vn1_0 (j,i)=(3.*vn1(j,1,i)-vn1(j,2,i))/2.
        tn1_0 (j,i)=(3.*tn1(j,1,i)-tn1(j,2,i))/2.
        fin1_0(j,i)=G_1000(j,i)+philb1(j,i)
        end do
      end do
c----------------------------------------------------------------
       do j=1,nb
         jplus1=j+1
         jm1 = max(1, j - 1)
         jp1 = min(nb,jplus1)
         do i=1,igit
            ip1=i+1
            im1=i-1
            if(i.eq.1   )   im1=igit
            if(i.eq.igit)   ip1=1
               phix = fin1_0(j,ip1)-fin1_0(j,im1)
               phiy1=(vn1_0(jp1,i)*fin1_0(jp1,i)+
     $                vn1_0(j  ,i)*fin1_0(j  ,i))*cosf2(jplus1)-
     $               (vn1_0(j  ,i)*fin1_0(j  ,i)+
     $                vn1_0(jm1,i)*fin1_0(jm1,i))*cosf2(j)
               phiy2=(vn1_0(jp1,i)+vn1_0(j  ,i))*cosf2(jplus1)-
     $               (vn1_0(j  ,i)+vn1_0(jm1,i))*cosf2(j)
         if(ncom.le.nphi) then
           philb0(j,i)=0.
            fphi_0(j,i)=-(un1_0(j,i)*phix/dx(jplus1)+
     $                   (phiy1-phiy2*fin1_0(j,i))/dy/cosfi(jplus1))*.5-
     $                    Rh*tn1_0(j,i)*wn1(j,1,i)
            f_phi = fphi_0(j,i)
         else  
            f_phi       =-(un1_0(j,i)*phix/dx(jplus1)+
     $                   (phiy1-phiy2*fin1_0(j,i))/dy/cosfi(jplus1))*.5-
     $                    Rh*tn1_0(j,i)*wn1(j,1,i)
         end if
      if(istep.eq.1) philb0(j,i)=philb1(j,i)
         philb1(j,i) = philb0(j,i)+dt*(f_phi-fphi_0(j,i))
        enddo  ! i - longitude 
       enddo   ! j - latitude
      
c      do k=48,56 
c      	print*,'vor 2201',k,tn1(1,k,1),tn2(1,k,1)
c      enddo
       
      do 2201 i=1,igit
      do 2201 k=1,kgit
      do 2201 j=1,nb
       un1(j,k,i)=un2(j,k,i)
       vn1(j,k,i)=vn2(j,k,i)
       tn1(j,k,i)=tn2(j,k,i)
 2201 continue
  333 continue
      if(ncom.ge.nphi) call legandr_philb
c--------- write just to check the behaviour of the lower bound. cond. 
      write(70,170)ncom,philb1(17,1),philb1(28,1),
     $                  fin1_0(28,1),fphi_0(28,1),wn1(17,1,1) 
  170 format(1x,i5,5e13.4)
c----------------------------------------------------------------------
      do 3 j=1,nb
      do 3 i=1,igit           
c---------------------------------------------------------------------
      do 100 istep=1,2
c--------------------------------------------------- molec. diffusion
c      if(i.eq.1.and.j.eq.1) then       
c             do k=48,56       
c                            print*,'jetzt molcon',i,j,k,tn1(j,k,i)            
c             enddo       
c      endif  
      
      call molcon_60_1(i,j)
      
      do 60 k=1,kgit
      if(istep.eq.1) then
c        fu(k) = fgru(k,i,j)
c        fv(k) = fgrv(k,i,j)
c        ft(k) = hs(k,i,j)+fc(k,i,j)+fgrt(k,i,j)
c
        un1_w(k)=un1(j,k,i)
        vn1_w(k)=vn1(j,k,i)
        tn1_w(k)=tn1(j,k,i)
      end if
c------------------------------------------------------- step in time 
c      un1(j,k,i)=un1_w(k)+dt2*(frunp1(k)+fu(k))
c      vn1(j,k,i)=vn1_w(k)+dt2*(frvnp1(k)+fv(k))
c      tn1(j,k,i)=tn1_w(k)+dt2*(fmcnp1(k)+ft(k))
      un1(j,k,i)=un1_w(k)+dt2*frunp1(k)
      vn1(j,k,i)=vn1_w(k)+dt2*frvnp1(k)
      tn1(j,k,i)=tn1_w(k)+dt2*fmcnp1(k)
   60 continue
  100 continue
    3 continue
cFL      print *,'u,v,t',un0(26,15,0),vn0(26,15,0),tn0(26,15,0)
cFL      print *,'u,v,t',un1(26,15,0),vn1(26,15,0),tn1(26,15,0)
cFL      print *,'u,v,t',un2(26,15,0),vn2(26,15,0),tn2(26,15,0)
      !hs(k,i,j) k=15 44km, i=0 0deg, j=26 42.5degN
      return
      end

 
