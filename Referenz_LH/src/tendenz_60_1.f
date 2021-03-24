      subroutine tendenz_60_1
      include 'com_main.fc'
      include 'com_namelist.fc' !FL
      real un1_w(kgit),vn1_w(kgit),tn1_w(kgit)
c--------------------------------------------------------------------
      dt2=dt/2.
c-------------------     smoothing, bi_hamonic duffusion of the momentum
c                                   and temperature
      do 1 j=1,nb
      do 1 i=1,igit
      do 2 istep=1,2
c--------------------------------------------------- molec. diffusion
      
c      if(i.eq.1.and.j.eq.1.) then
c      	do k=48,56
c      		print*,'vorher',tn1(j,k,i)
c      	enddo	
c      endif
      
      call molcon_60_1(i,j)
      
c      if(i.eq.1.and.j.eq.1) then
c      	do k=48,56
c      		print*,'nachher',tn1(j,k,i)	
c      	enddo
c      endif
c
      do k=1,kgit
      if(istep.eq.1) then
        un1_w(k)=un1(j,k,i)
        vn1_w(k)=vn1(j,k,i)
        tn1_w(k)=tn1(j,k,i)
      end if
c------------------------------------------------------- step in time 
      un1(j,k,i)=un1_w(k)+dt2*frunp1(k)
      vn1(j,k,i)=vn1_w(k)+dt2*frvnp1(k)
      tn1(j,k,i)=tn1_w(k)+dt2*fmcnp1(k)
      end do
    2 continue
    1 continue
      return
      end

 
