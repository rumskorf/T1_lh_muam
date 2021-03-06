c     Alex Pogoreltsev, 2006
c     tlbd has been introduced (SPWs started after 30 days) and
c     tsec is calculated with the characteristic time 30 days
c-----------------------------------------------------------------
      subroutine geopot_60
      include 'com_main.fc'
      include 'com_namelist.fc' !FL
c-----------------------------------------------------------------
      Rh    = 287.04/h      ! R/M/H
      dzrh  = Rh*dz         ! dz*R/M/H
      dzrm  = dzrh*rm(1)    ! dz*R/H 
      tlbd=float(nsec)-30.*86400.
      if(tlbd.gt.0.) then
      ttau = tlbd/30./86400.
      tsec=1.-exp(-ttau*ttau)
      do j=1,nb
        do i=1,igit
        if(ncom.le.nphi) philb1(j,i)=0.
c-----------------------------------------------------------------
c        to simulate with old lower boundary condition
c        if(ncom.le.nend) philb1(j,i)=0.
c----------------------------------------- T and PHI on the ground
        tn1_0  = (3.*tn1(j,1,i)-tn1(j,2,i))/2.
        fin1_0 = G_1000(j,i)+philb1(j,i)
c-----------------------------------------------------------------
c             G_1000(j,i) includes zonally averaged filed and SPWs

c------------------------- SPWs and zonally averaged geopotential,
c                i.e., G_1000(j,i) or fin1_0(j,i) are given at z=0
c-----------------------------------------------------------------
      fin1(j,1,i)=tsec*(fin1_0+dzrh/4.*(tn1(j,1,i)+tn1_0))
        end do
      end do
      else
      do j=1,nb
        do i=1,igit
          fin1(j,1,i)=0.
        end do
      end do
      end if
c------------------------------------------------------------------
       do j=1,nb
        do i=1,igit 
         do k =1,kgit-1
          k1=k+1
          fin1(j,k1,i)= fin1(j,k,i)+
     $               (tn1(j,k,i)/rm(k)+tn1(j,k1,i)/rm(k1))/2.*dzrm 
         enddo ! k - height
        enddo  ! i - longitude 
       enddo   ! j - latitude
      rgit=float(igit)
      do k=1,1
        sum1=0.
        sum2=0.
        sum3=0.
        sum4=0.
        do i=1,igit
         sum1=sum1+fin1(1   ,k,i)/rgit
         sum2=sum2+fin1(2   ,k,i)/rgit
         sum3=sum3+fin1(nb-1,k,i)/rgit
         sum4=sum4+fin1(nb  ,k,i)/rgit
        end do
        do i=1,igit
         fin1(1 ,k,i)=sum1+(fin1(2   ,k,i)-sum2)/3. 
         fin1(nb,k,i)=sum4+(fin1(nb-1,k,i)-sum3)/3.
        end do  
      end do
      return
      end




