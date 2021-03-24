      subroutine vbew
c     berechnung der vertikalgeschwindigkeit
c
      include 'com_main.fc'
      include 'com_namelist.fc' !FL
      roukp1=exp(-(z(kgit)*1000.+dz)/h)
      dz5=dz*0.5
      do 500 j=1,nb
         jm1 = max(1, j - 1)
         jp1 = min(nb,j + 1)
         do 500 i=1,igit
            ip1=i+1
            im1=i-1
            if(i.eq.1   )   im1=igit
            if(i.eq.igit)   ip1=1
               dvx= un1(j,  kgit,ip1)-un1(j,  kgit,im1)
               dvy=(vn1(jp1,kgit,i  )+vn1(j  ,kgit,i  ))*cosf2(j+1)
     +            -(vn1(j,  kgit,i  )+vn1(jm1,kgit,i  ))*cosf2(j)
               div=(dvx/dx(j+1)+dvy/dy/cosfi(j+1))*dz5*roukp1
               wn1(j,kgit+1,i)= div 
c            wn1(j,kgit,i)=0.
            do k=kgit,1,-1
c               kp1   =min(k+1,kgit)
               kp1=k+1
               dvx= un1(j,  k,ip1)-un1(j,  k,im1)
               dvy=(vn1(jp1,k,i  )+vn1(j  ,k,i  ))*cosf2(j+1)
     +            -(vn1(j,  k,i  )+vn1(jm1,k,i  ))*cosf2(j)
               div=(dvx/dx(j+1)+dvy/dy/cosfi(j+1))*dz5*rou(k)
               wn1(j,k,i)= wn1(j,kp1,i) + div 
            enddo
  500 continue
      rgit=float(igit)
      do k=1,1
        sum1=0.
        sum2=0.
        sum3=0.
        sum4=0.
        do i=1,igit
         sum1=sum1+wn1(1   ,k,i)/rgit
         sum2=sum2+wn1(2   ,k,i)/rgit
         sum3=sum3+wn1(nb-1,k,i)/rgit
         sum4=sum4+wn1(nb  ,k,i)/rgit
        end do
        do i=1,igit
         wn1(1 ,k,i)=sum1+(wn1(2   ,k,i)-sum2)/3. 
         wn1(nb,k,i)=sum4+(wn1(nb-1,k,i)-sum3)/3.
        end do  
      end do
      return
      end


