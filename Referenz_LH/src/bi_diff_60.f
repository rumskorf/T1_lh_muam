c-----------------------------------------------------------------------
c     Alex Pogoreltsev, January 2006
c
c     1. vertical bi-harmonic diffusion at the level k=1 is removed
c
c     2. the influence of diff_z on solution was tested using diff_z/5,
c        there were not significant changes in zonal mean wind and SPW
c        additionally, diff_z to smooth the temperature is the same,
c        i.e., without of /3.
c        FINALLY, REMAIN diff_z WITHOUT /5 (TO HAVE A STABLE SOLUTION NEAR
c        THE EQUATOR WITHOUT FRICTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c 
c     3. diff_H is increased (*2e6) (probably this is not necessary?)
c
c-----------------------------------------------------------------------
      subroutine bi_diff_60
c------------------------- THIS HAS TO BE REMOVED IF WE USE NCEP/NCAR
c                          TEMPERATURE IN THE TROPOSPHERE!!!!!!!!!!!!
c------------------------- includes diffusion for the first level, i.e.,
c                          temperature is connected with T_1000
c-----------------------------------------------------------------------     
      include 'com_main.fc'
      include 'com_namelist.fc' !FL
      real uxy,vxy,txy
      data a /6366197./
      cappa=RgSI/rm(1)/c_p
      dif2_z=30.
      dzdz = dz*dz  
      dz4  = dzdz*dzdz
      kgm1=kgit-1
      kgm2=kgit-2
c------------------ HORIZONTAL DIFFUSION,  based on the equations suggested by
c Marchuk, G.I., V. Dymnikov, V. Zalesny, V. Lykossov, and V. Ya. Galin, 1984:
c Mathematical Modeling of the General Circulation of the Atmosphere and Ocean. 
c Gidrometeoizdat, 320 pp. (in Russian).
c-----------------------------------------------------------------------------
c      d_bh   = 1.e18
c      ttau = float(nsec)/10./86400.
      do 2 k=1,kgit
c      diff_H =  1.e6 *(1.5+0.5*tanh((z(k)-100.)/20.))
c     $         +1.e6 *EXP((z(k)-z(kgit))/20.)
c      d_bh   =  1.e18*(1.5+0.5*tanh((z(k)-100.)/20.))
c     $         +1.e18*EXP((z(k)-z(kgit))/20.)
c------------------------ two times weaker in the lower atmosphere
c                         and very strong smoothing at the beginning
      diff_H =  0.5e6 *(1.5+0.5*tanh((z(k)-120.)/20.))
c     $         +1.0e6 *exp(-ttau*ttau)*(0.5+0.5*tanh((z(k)-100.)/20.))
c     $  + 1.0e6 *exp((z(k)-120.)/20.)
c      d_bh   =  1.0e17*(1.5+0.5*tanh((z(k)-120.)/20.))
c     $         +1.0e18*exp(-ttau*ttau)*(0.5+0.5*tanh((z(k)-100.)/20.))
c     $  + 1.0e18*exp((z(k)-120.)/20.)
c-------------------------------------- latitude
      do 2 j=1,nb
         jplus1=j+1
c      diff_H =  1.0e6 *(1.5+0.5*tanh((z(k)-120.)/20.))*cosfi(jplus1)
      d_bh   =  0.5e18*(1.5+0.5*tanh((z(k)-120.)/20.))*cosfi(jplus1)
         jm1 = max(1 ,j - 1)
         jp1 = min(nb,jplus1)
         dey = 1./dy/cosfi(jplus1)
         dx0 = a*cosfi(jplus1)
         dx2 = dx(jplus1)*dx(jplus1)/cosfi(jplus1)
         dx4 = dx2*dx2 
c-------------------------------------- longitude
      do 2 i=1,igit
	    ip1 = i+1
	    im1 = i-1
	    if(i.eq.1   ) im1=igit
	    if(i.eq.igit) ip1=1
	    ip2 = i+2
	    im2 = i-2
	    if(i.eq.1   ) im2=igit-1
	    if(i.eq.igit) ip2=2
	    if(i.eq.2     ) im2=igit
	    if(i.eq.igit-1) ip2=1
c--------------------------------------
        un_00  = un1(j,k,i)
        un1ip  = un1(j,k,ip1)
        un1im  = un1(j,k,im1)
        un1ip2 = un1(j,k,ip2)
        un1im2 = un1(j,k,im2)
        un1jp  = un1(jp1,k,i)
        un1jm  = un1(jm1,k,i)
c        un1p05 = (un1jp+un_00)/2.
c        un1m05 = (un1jm+un_00)/2.
c--------------------------------------
        vn_00  = vn1(j,k,i)
        vn1ip  = vn1(j,k,ip1)
        vn1im  = vn1(j,k,im1)
        vn1ip2 = vn1(j,k,ip2)
        vn1im2 = vn1(j,k,im2)
        vn1jp  = vn1(jp1,k,i)
        vn1jm  = vn1(jm1,k,i)
c        vn1p05 = (vn1jp+vn_00)/2.
c        vn1m05 = (vn1jm+vn_00)/2.
c----------------------------------------------------------------------
c     dx(j+1) = a*COS(varphi)*d_lambda
c     dx0     = a*COS(varphi)
c     dx2     = a^2*COS(varphi)*(d_lambda)^2
c     dx4     = a^4*COS^2(varphi)*(d_lambda)^4
c     dy      = a*d_varphi
c     dey     = 1./(a*COS(varphi)*d_varphi)
c
c         fi    =(float(j)-19.5)*defi/180.*pi   ! -92.5,...,87.5 North
c         fi2   =(float(j)-19.0)*defi/180.*pi   ! -90.0,...,90.0 North
c         i.e., sinfi & cosfi has to be used for (j+1) point and
c               sinf2 & cosf2 related to the point between j and j+1  
c--------------------------------------------------------------------- 
      uxy = (un1ip-2.*un_00+un1im)/dx(jplus1)/dx(jplus1)-
     $      2.*sinfi(jplus1)*(vn1ip-vn1im)/2./dx(jplus1)/dx0
c      uxy = uxy + dey*dey*
c     $      ((un1jp - un_00)*cosf2(jplus1)*cosf2(jplus1)-
c     $       (un_00 - un1jm)*cosf2(j     )*cosf2(j     ))+
c     $            dey/dx0*
c     $  (un1p05*sinf2(jplus1)*cosf2(jplus1)-un1m05*sinf2(j)*cosf2(j))
      uxy = uxy + dey*dey*
     $      ((un1jp/cosfi(j+2   )-un_00/cosfi(jplus1))*
     $       cosf2(jplus1)*cosf2(jplus1)*cosf2(jplus1)-
     $       (un_00/cosfi(jplus1)-un1jm/cosfi(j     ))*
     $       cosf2(j     )*cosf2(j     )*cosf2(j     ))
      fut(k,i,j) = diff_H*uxy-
     $      d_bh*(un1im2-4.*un1im+6.*un_00-4.*un1ip+un1ip2)/dx4
c----------------------------------------------------------------------
      vxy = (vn1ip-2.*vn_00+vn1im)/dx(jplus1)/dx(jplus1)+
     $      2.*sinfi(jplus1)*(un1ip-un1im)/2./dx(jplus1)/dx0
c      vxy = vxy + dey*dey*
c     $      ((vn1jp -vn_00)*cosf2(jplus1)*cosf2(jplus1)-
c     $       (vn_00 -vn1jm)*cosf2(j     )*cosf2(j     ))+
c     $            dey/dx0*
c     $  (vn1p05*sinf2(jplus1)*cosf2(jplus1)-vn1m05*sinf2(j)*cosf2(j))
      vxy = vxy + dey*dey*
     $      ((vn1jp/cosfi(j+2   )-vn_00/cosfi(jplus1))*
     $       cosf2(jplus1)*cosf2(jplus1)*cosf2(jplus1)-
     $       (vn_00/cosfi(jplus1)-vn1jm/cosfi(j     ))*
     $       cosf2(j     )*cosf2(j     )*cosf2(j     ))
      fvt(k,i,j) = diff_H*vxy-
     $      d_bh*(vn1im2-4.*vn1im+6.*vn_00-4.*vn1ip+vn1ip2)/dx4
c----------------------------------------------------------------------
        tn_00 = tn1(j,k,i)
        tn1ip = tn1(j,k,ip1)
        tn1im = tn1(j,k,im1)
        tn1ip2= tn1(j,k,ip2)
        tn1im2= tn1(j,k,im2)
        tn1jp = tn1(jp1,k,i)
        tn1jm = tn1(jm1,k,i)
c        tn1p05= (tn1jp+tn_00)/2.
c        tn1m05= (tn1jm+tn_00)/2.
      txy = (tn1ip-2.*tn_00+tn1im)/dx(jplus1)/dx(jplus1)
      txy = txy + dey/dy*
     $      ((tn1jp-tn_00)*cosf2(jplus1)-(tn_00-tn1jm)*cosf2(j))
      ftt(k,i,j) = diff_H*txy-
     $      d_bh*(tn1im2-4.*tn1im+6.*tn_00-4.*tn1ip+tn1ip2)/dx4
    2 continue
c
c-------------------------   vertical bi-harmonic diffusion smoothing
c      g=9.81
c      RgSI=8314.41
c      h=7000.
c      c_p=1005.
c         Rf_c=0.4
c-(0.5e8) - lamb_z = 20 km, effective friction coeff. is about of 0.5e-6
c      diff_z=1.e8
c- diffusion near the lower boundary, lamb_z = 20 km, friction is 0.3e-5
c------------------------------------------------- dif2_z = Pr * 10 m2/s
      do 3 j=1,nb
      do 3 i=1,igit     
      do k=3,kgm2
      diff_z = (0.75+0.25*tanh((z(k)-50.)/10.))*1.e8
      fut(k,i,j)=fut(k,i,j)-diff_z*(un1(j,k-2,i)-4.*un1(j,k-1,i)+
     $          6.*un1(j,k,i)+un1(j,k+2,i)-4.*un1(j,k+1,i))/dz4
      fvt(k,i,j)=fvt(k,i,j)-diff_z*(vn1(j,k-2,i)-4.*vn1(j,k-1,i)+
     $          6.*vn1(j,k,i)+vn1(j,k+2,i)-4.*vn1(j,k+1,i))/dz4
      ftt(k,i,j)=ftt(k,i,j)-diff_z*(tn1(j,k-2,i)-4.*tn1(j,k-1,i)+
     $          6.*tn1(j,k,i)+tn1(j,k+2,i)-4.*tn1(j,k+1,i))/dz4
      end do
c-----------------------------------------------------------------------
      fut(kgm1,i,j) = fut(kgm1,i,j)+dif2_z*(un1(j,kgm2,i)-
     $                          2.*un1(j,kgm1,i)+un1(j,kgit,i))/dzdz
      fvt(kgm1,i,j) = fvt(kgm1,i,j)+dif2_z*(vn1(j,kgm2,i)-
     $                          2.*vn1(j,kgm1,i)+vn1(j,kgit,i))/dzdz
      ftt(kgm1,i,j) = ftt(kgm1,i,j)+dif2_z*(tn1(j,kgm2,i)-
     $                          2.*tn1(j,kgm1,i)+tn1(j,kgit,i))/dzdz
c------  an assumption that wind and T near the upper boundary are CONST
      fut(kgit,i,j)=fut(kgit,i,j)+
     $                     dif2_z*(un1(j,kgm1,i)-un1(j,kgit,i))/dzdz
      fvt(kgit,i,j)=fvt(kgit,i,j)+
     $                     dif2_z*(vn1(j,kgm1,i)-vn1(j,kgit,i))/dzdz

      ftt(kgit,i,j)=ftt(kgit,i,j)+
     $                     dif2_z*(tn1(j,kgm1,i)-tn1(j,kgit,i))/dzdz
c-----------------------------------------------------------------------

      do k=1,2 
      U_k   =  un1(j,k,i)
      U_pl  =  un1(j,k+1,i)
      V_k   =  vn1(j,k,i)
      V_pl  =  vn1(j,k+1,i)
      T_k   =  tn1(j,k,i)
      T_pl  =  tn1(j,k+1,i)
      if(k.gt.1) then
        U_mn  =  un1(j,k-1,i)
        V_mn  =  vn1(j,k-1,i)
        T_mn  =  tn1(j,k-1,i)
      else
        U_mn  = (3.*un1(j,1,i)-un1(j,2,i))/2.
        V_mn  = (3.*vn1(j,1,i)-vn1(j,2,i))/2.
c-------------------  T increases to the ground with the gradient T_grad
        T_grad = 2.*(T_k-T_1000(j,i))/dz 
        T_mn   = T_k-T_grad*dz
      end if
      U_zz  = (U_mn-2.*U_k+U_pl)/dzdz 
      U_grad= (U_pl-U_mn)/2./dz
      V_zz  = (V_mn-2.*V_k+V_pl)/dzdz 
      V_grad= (V_pl-V_mn)/2./dz
      T_zz  = (T_mn-2.*T_k+T_pl)/dzdz 
      T_grad= (T_pl-T_mn)/2./dz
      hscale=RgSI*T_k/g/rm(k)
      fut(k,i,j)=fut(k,i,j)+dif2_z*(U_zz-(1./h+2.*T_grad/T_k)*U_grad)
      fvt(k,i,j)=fvt(k,i,j)+dif2_z*(V_zz-(1./h+2.*T_grad/T_k)*V_grad)
      ftt(k,i,j)=ftt(k,i,j)+dif2_z/3.* 
     $           (T_zz-2.*T_grad*T_grad/T_k-g/c_p*hscale/h*T_grad/T_k-
     $           (1.-cappa/Rf_c)/h*(T_grad+hscale/h*g/c_p))
      end do
    3 continue
      RETURN
      END






