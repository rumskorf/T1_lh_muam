      subroutine levy

      include 'com_main.fc'
      include 'com_namelist.fc' !FL

      dimension xmin(nb),kmin(nb)

c      nstep  =24*ntime          ! zeitschritte/tag
      tleps=1.e-9

      do j=1,nb
         amax=1.e6
         bmax=1.e6
         cmax=1.e6
         dtx =dt*1.
         do i=1,igit
            do k=1,kgit
               tlevyx=dx(j+1)/(abs(un1(j,k,i)*dtx)+tleps)
               tlevyy=dy     /(abs(vn1(j,k,i)*dtx)+tleps)
               tlevyz=dz     /(abs(wn1(j,k,i)*dtx)+tleps)
               if (tlevyx.lt.amax) then
                  amax=tlevyx
                  kk=k
               endif
               if (tlevyy.lt.bmax)  then
                  bmax=tlevyy
c                  kk=k
               endif
               if (tlevyz.lt.cmax)  then
                  cmax=tlevyz
c                  kk=k
               endif
            enddo
         enddo
         xmin(j)=min(amax,bmax,cmax)
         kmin(j)=kk
      enddo

c     kriterium erfuellt? (levy =2 --> c =500 m/s)

      amin=1.e6
      do j=1,nb
         if (xmin(j).lt.amin) then
            amin =xmin(j)
            kkmin=kmin(j)
            jjmin=j
         endif
      enddo

      print 348,ncom/nstep,ncom,amax,bmax,cmax,jjmin,kkmin
348   format 
     +     (' day ',i3,' step ',i10,' Levy-C. ',f9.3,'X',f7.0,'X',f7.0,
     +     ' at ',2i3)    

      return
      end

