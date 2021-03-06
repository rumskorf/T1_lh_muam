      program MUAM_48  
c.....COMMON-BLOCKS
      include 'com_main.fc'
      include 'com_namelist.fc' !FL
c-----------------------------------------------------------------------
      real rmyfin1(igit,nb,kgit),f107men
c      COMMON/adams/funm1(nb,kgit,igit),fvnm1(nb,kgit,igit),
c     $            ftnm1(nb,kgit,igit)      
      integer pmac,numb
csl
csl      COMMON JAHR,NDEK,f107mean
csl      
c-----------------------------------------------------------------------
      print*
      print*,'*******************************************************'
      PRINT*,'*                                                     *'
      print*,'*                        M U A M                      *'
      PRINT*,'*                                                     *'
      PRINT*,'*           Middle  and  Upper  Atmosphere  Model     *'
      print*,'*                                                     *'
      print*,'*                                                     *'
      print*,'*******************************************************'
c------------------------------------------------ first record               
       pmac=0.
       nc10=0.
       numb=0.
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
C.....READ EXTERNAL DATA:
      call open_60      
      print*
      print*,'ALL DATA READ'
      print*
c
c.....INITIALIZE MODEL:
      call init_60    
      call fourier_matrix
      call fourier_filter_lb
      print*
      print*,'MODEL INITIALIZED'
      print*
c
c.....ADJUSTMENTS FOR IR RADIATION (Berger'92)
      call irc_60
c
      print*,'step5'
C
************************************************************************ 
*                                                                      *
*     START OF TIME LOOP . . .                                         *
*                                                                      *
************************************************************************
c
      print*
      print*,'++++++ STARTING TIME INTEGRATION +++++++'
      print*
      null=ncom                 ! initial step from pre-run
CPH      null=ncom-5760
CPH      null=ncom-192
CPH      null=ncom-192*5
c      nsec=ncom*int(dt)
      print*
      print*,'first step is ',null,' last is ',nend
      print*
c
      print*,'dz=',dz
c
c------------------------------------- FIRST CALL
      do 100 ncom=null,nend
c
      do j=1,nb
       do i=1,igit
        firn1(j,i)=0.
        do k=1,kgit
        fnt(k,i,j)=0.
        frl(k,i,j)=0.
        heat_LH (i,j,k) = 0.
        if(k.gt.10) heat_LH0(i,j,k) = 0.
        end do 
       end do
      end do
      IF(ncom.eq.null) then
c-----------------------     GEOPOTENTIAL (from hydrostatic equation)
      call geopot_60 
c-----------------------     VERTICAL WIND (from continuity equation)
      call vbew
c	
c.....OUTPUT -II- OF BINARIES in DIRECT-ACCESS-MODE:
c     the analyse package can read this
c     !record-Nr.
	pmac=pmac+1
        do k=1,kgit
           do j=1,nb
              do i=1,igit
                 rmyfin1(i,j,k)=fin1(j,k,i)/9.81
                 rmywn1 (i,j,k)=wn1(j,k,i)/row(k)
                 do m=1,3
                    rmyan1(i,j,k,m)=an1(j,k,i,m)
                 enddo
              enddo
           enddo
        enddo
        write(72,rec=pmac) rmyan1 
        write(24,rec=pmac) rmyfin1 
        write(25,rec=pmac) rmywn1
c--------------------------------------------------------------------
C.....PARAMETRIZATIONS:
c     some calculations are not repeated every time step to
c     save cpu-time (e.g. radiation and gravity-wave modules)
c 
c     But THIS MAY BE NECCESSARY if you need more accurate results.
c
c-------------------------------     GRAVITY WAVES
      call grwaves_RSHU

        write(26,rec=pmac) fluxgru            !(kgit,igit,nb) by Rliu 
        write(27,rec=pmac) fluxgrv 			!(kgit,igit,nb) by Rliu
        write(28,rec=pmac) fgru 				!(kgit,igit,nb) by Rliu
        write(29,rec=pmac) fgrv 				!(kgit,igit,nb) by Rliu
        write(30,rec=pmac) fgrt 
        
csl
csl       CALL  sub_F107mean(JAHR,NDEK,f107mean)
csl      
c          call fourier_filter_GW       
c-------------------------------     UV, IR
        if(muvir.ge.1) call strobel_60
        if(muvir.ge.2) call ircool_60
c-------------------------------     Latent heating
        call heating_LH
c-----------------------     FORCING of transient waves
      if(mforc.eq.2.and.ncom.ge.nphi) call plwaves_LIM
      go to 200
      END IF
c--------------------------------------------------------------------
C.....PARAMETRIZATIONS:
c     some calculations are not repeated every time step to
c     save cpu-time (e.g. radiation and gravity-wave modules)
c 
c     But THIS MAY BE NECCESSARY if you need more accurate results.
c
c-------------------------------     GRAVITY WAVES
      IF(mod(ncom,2*ntime).eq.0) then
          call grwaves_RSHU
          call acc_ASM
c          call fourier_filter_GW
      END IF
c
c-------------------------------     UV, IR
      IF(ncom.lt.nphi) then
        if(muvir.ge.1.and.mod(ncom,ntime).eq.0) call strobel_60
        if(muvir.ge.2.and.mod(ncom,ntime).eq.0) call ircool_60
      ELSE
        if(muvir.ge.1.and.mod(ncom,ntime).eq.0) call strobel_60
        if(muvir.ge.2.and.mod(ncom,ntime).eq.0) call ircool_60
      END IF
c-----------------------     FORCING of transient waves
      if(mforc.eq.2.and.ncom.ge.nphi) call plwaves_LIM
c
c...........................PROGNOSTIC EQUATIONS / INTEGRATION STEP 1
      call tendenz_60_1
c......................................................FOURIER-FILTER
c      call fourier_filter
c
      if(ncom.eq.null+240*ntime) then
      do k=1,kgit
      print*,'z=',z(k),' h=',hs(k,1,1),' c=',fc(k,1,1),' g=',fgrt(k,1,1)
      enddo
      end if
c
c...........................PROGNOSTIC EQUATIONS / INTEGRATION STEP 2
      call tendenz_60_a
c
c.....INCREASE TIME-INCREMENT:
      nsec=nsec+int(dt)
c--------------- correction near the poles, A.Pogoreltsev, March 2006
      rgit=float(igit)
      do k=1,kgit
        sum1=0.
        sum2=0.
        sum5=0.
        sum6=0.
        do i=1,igit
         sum1=sum1+tn1(1,k,i)/rgit
         sum2=sum2+tn1(2,k,i)/rgit
         sum5=sum5+tn1(nb-1,k,i)/rgit
         sum6=sum6+tn1(nb  ,k,i)/rgit
        end do
        do i=1,igit
         tn1(1,k,i)=sum1+(tn1(2,k,i)-sum2)/3. 
         tn1(nb  ,k,i)=sum6+(tn1(nb-1,k,i)-sum5)/3.
        end do  
      end do
c......................................................FOURIER-FILTER
      call fourier_filter
c.....OUTPUT ZONALLY AVERAGED VALUES
c     two ascii formats possible
c      if(mod(ncom,nprint).eq.0) call ausgaz
c      if(ncom.eq.nend) call ausgaz
c
c.....SAVE FORT.2=an0,an1,an2,...
c
c     a rewind of fort.2 saves disk space
c
 80   if(kret)    89,86,87   !kret<0 NO SAVE
			     !    =0 SAVE ONLY VERY LAST STEPS
			     !    >0 ONE SAVE PER DAY
c
 87   if(mod(ncom,nstep)) 89,88,89   
 86   if(    ncom-nend)   89,88,89   
c
 88   continue
      write (2) an0
      write (2) an1
      write (2) an2
      write (2) philb0
      write (2) philb1
      write (2) fphi_0
      write (2) ken1 
      close (2)
      print*, ' f2_new saved!'
c      write (113,1013) an0
c      write (113,1013) an1
c      write (113,1013) an2
c      write (113,1013) philb0
c      write (113,1013) philb1
c      write (113,1013) fphi_0
c      write (113,2013) ken1 
c      close (113)
 1013 format(10e13.6)
 2013 format(10i6)
      print*, ' f2_nfm saved!'
c
 89   if(mod(ncom,nout)) 90,92,90     
 92   continue
c	
c.....OUTPUT -II- OF BINARIES in DIRECT-ACCESS-MODE:
c     the analyse package can read this
c     !record-Nr.
	pmac=pmac+1
c-----------------------     GEOPOTENTIAL (from hydrostatic equation)
      call geopot_60 
c-----------------------     VERTICAL WIND (from continuity equation)
      call vbew
        do k=1,kgit
           do j=1,nb
              do i=1,igit
                    rmyfin1(i,j,k)=fin1(j,k,i)/9.81
                    rmywn1 (i,j,k)= wn1(j,k,i)/row(k)
                 do m=1,3
                    rmyan1(i,j,k,m)=an1(j,k,i,m)
                 enddo
              enddo
           enddo
        enddo
        write(72,rec=pmac) rmyan1 
        write(24,rec=pmac) rmyfin1 
        write(25,rec=pmac) rmywn1

        write(26,rec=pmac) fluxgru            !(kgit,igit,nb) by Rliu 
        write(27,rec=pmac) fluxgrv 			!(kgit,igit,nb) by Rliu
        write(28,rec=pmac) fgru 				!(kgit,igit,nb) by Rliu
        write(29,rec=pmac) fgrv 				!(kgit,igit,nb) by Rliu
        write(30,rec=pmac) fgrt 

   90 continue
c
      nc10 = nc10+1
  200 continue
c
c     LEVY-COURANT-CHECK (dX/c/dt > 1 !)
      ncase=mod(ncom,nstep/4)*(mod(ncom,10)+max0(0,nc10-100))
cFL      if(ncase.eq.0) call levy !FL hier bildschirmausgabe
      call levy !jeden zeitschritt ausgeben
  100 continue
c
*****************************************************************
*                                                               *
*     END OF TIME LOOP . . .                                    *
*                                                               *
*****************************************************************
      print 1234,ncom-1,nend
 1234 format(' ncom, nend ',2i6)
      STOP '++++++ THE END ++++++'
      end




