      subroutine acc_ASM
      include 'com_main.fc'
      include 'com_namelist.fc'
! Input parameters for egwd
      REAL*8 :: vy1d(kgit)     !v-wind
      REAL*8 :: vx1d(kgit)     !u-wind
      REAL*8 :: temp1d(kgit)   ! 
      REAL*8 :: cp1d(kgit)     ! Cp
      REAL*8 :: hm(kgit)        !geopotential height (in m)
      REAL*8 :: pres(kgit)     !pressure (in Pa)
      REAL*8 :: rho(kgit)      !density
! Output variables for egwd
      REAL*8 :: ut_gwd(kgit)
      REAL*8 :: vt_gwd(kgit)
     &         ,gwh_ir(kgit), gwh_dif(kgit)
     &         ,var_tot(kgit), flux_tot(kgit)
      real acc(igit,nb,kgit,3)
      real fgru_AM(kgit,igit,nb),fgrv_AM(kgit,igit,nb),
     &     fgrt_AM(kgit,igit,nb)
      real work1(nb),work2(nb),work3(nb)
      
      character*6 ncomchar
     
      
      eff  = 1.-exp(-float(nsec)/30./86400.)
       do i=1,igit
         do j=1,nb
           do k=1,kgit
             vx1d(k)=an1(j,k,i,1)
             vy1d(k)=an1(j,k,i,2)
             temp1d(k)=an1(j,k,i,3)

      !AM  Below must be the GEOPOTENTIAL HEIGHT!! We put the log-pressure
	!    height for trial only. Changde this!!! 
c-------------------------------------------------------------------------
c        can be calculated using the geopotential height from the solution
c-------------------------------------------------------------------------
	     hm(k)=fin1(j,k,i)/9.81
             cp1d(k)=1005.0
             pres(k)=1.e5*exp(-z(k)/7.)
             rho(k)=pres(k)*rm(k)/an1(j,k,i,3)/RgSI
           enddo
	
           call EYGwave_RSHU(i,j,geddy,rho,vx1d,vy1d,temp1d,cp1d,hm, 
     &       ut_gwd, vt_gwd, pres, gwh_ir, gwh_dif,
     &       var_tot, flux_tot)
	     
           do k=1,kgit
	          acc(i,j,k,1)=real(ut_gwd(k))
		   acc(i,j,k,2)=real(vt_gwd(k))
		   acc(i,j,k,3)=real(gwh_ir(k)+gwh_dif(k))
	   enddo
         enddo
       enddo
c------------------------------------ smoothing in longitude
       do j=1,nb
        do k=1,kgit
         do i=1,igit
           ip1=i+1
           ip2=i+2
           im1=i-1
           im2=i-2
           if(i.eq.1) then
            im1=igit
            im2=igit-1
           end if
           if(i.eq.2)im2=igit
           if(i.eq.igit) then
            ip1=1
            ip2=2
           end if
           if(i.eq.igit-1)ip2=1
           fgru_AM(k,i,j) = eff*(acc(im2,j,k,1)+acc(im1,j,k,1)+
     &      acc(i,j,k,1)+acc(ip2,j,k,1)+acc(ip1,j,k,1))/5.  
           fgrv_AM(k,i,j) = eff*(acc(im2,j,k,2)+acc(im1,j,k,2)+
     &      acc(i,j,k,2)+acc(ip2,j,k,2)+acc(ip1,j,k,2))/5. 
           fgrt_AM(k,i,j) = eff*(acc(im2,j,k,3)+acc(im1,j,k,3)+
     &      acc(i,j,k,3)+acc(ip2,j,k,3)+acc(ip1,j,k,3))/5.
c           geddy(k,i,j)=geddy0(k)
         enddo
        enddo
       enddo
c------------------------------------ smoothing in latitude
c                               Alex Pogoreltsev 20.03.2013
c----------------------------------------------------------
      do i=1,igit
       do k=8,kgit
        do j=1,nb
         work1(j)=fgru_AM (k,i,j)
         work2(j)=fgrv_AM (k,i,j)
         work3(j)=fgrt_AM (k,i,j)
        end do
        do j=2,nb-1
          jp1=j+1
          jm1=j-1
      fgru_AM (k,i,j)= (work1(jm1)+work1(j)+work1(jp1))/3.  
      fgrv_AM (k,i,j)= (work2(jm1)+work2(j)+work2(jp1))/3.  
      fgrt_AM (k,i,j)= (work3(jm1)+work3(j)+work3(jp1))/3.  
        end do
       end do
      end do
      do j=1,nb
        do k=3,kgit-2
         do i=1,igit
           acc(i,j,k,1)=fgru_AM(k,i,j)
           acc(i,j,k,2)=fgrv_AM(k,i,j)
           acc(i,j,k,3)=fgrt_AM(k,i,j)
         end do
        end do
      end do
      do j=1,nb
        do k=3,kgit-2
         do i=1,igit
           kp1=k+1
           kp2=k+2
           km1=k-1
           km2=k-2
           fgru(k,i,j)= fgru(k,i,j)+(acc(i,j,km2,1)+acc(i,j,km1,1)+
     &      acc(i,j,k,1)+acc(i,j,kp1,1)+acc(i,j,kp2,1))/5.  
           fgrv(k,i,j)= fgrv(k,i,j)+(acc(i,j,km2,2)+acc(i,j,km1,2)+
     &      acc(i,j,k,2)+acc(i,j,km1,2)+acc(i,j,kp2,2))/5. 
           fgrt(k,i,j)= fgrt(k,i,j)+(acc(i,j,km2,3)+acc(i,j,km1,3)+
     &      acc(i,j,k,3)+acc(i,j,kp1,3)+acc(i,j,kp2,3))/5. 
         enddo
        enddo
       enddo
       
       print*,ncom,fgru(19,18,28),'acc'
       
       
      return
      end


