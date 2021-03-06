c      program MUAM_48  
	subroutine pyrun()
c.....COMMON-BLOCKS
      include 'com_main.fc'
      include 'com_namelist.fc' !FL
c-----------------------------------------------------------------------
      real rmyfin1(igit,nb,kgit)
c      COMMON/adams/funm1(nb,kgit,igit),fvnm1(nb,kgit,igit),
c     $            ftnm1(nb,kgit,igit)      
      integer pmac
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
      print*
c------------------------------------------------ first record               
       pmac=0.
       nc10=0
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
c          call fourier_filter_GW
c-------------------------------     UV, IR
        if(muvir.ge.1) call strobel_60
        if(muvir.ge.2) call ircool_60
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
c          call fourier_filter_GW
      END IF
c
c-------------------------------     UV, IR
      IF(ncom.lt.nphi) then
        if(muvir.ge.1.and.mod(ncom,2*ntime).eq.0) call strobel_60
        if(muvir.ge.2.and.mod(ncom,2*ntime).eq.0) call ircool_60
      ELSE
        if(muvir.ge.1.and.mod(ncom,2*ntime).eq.0) call strobel_60
        if(muvir.ge.2.and.mod(ncom,4*ntime).eq.0) call ircool_60
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
   90 continue
c
      nc10 = nc10+1
  200 continue

c     LEVY-COURANT-CHECK (dX/c/dt > 1 !)
      ncase=mod(ncom,nstep/4)*(mod(ncom,10)+max0(0,nc10-100))
      if(ncase.eq.0) call levy
  100 continue

*****************************************************************
*                                                               *
*     END OF TIME LOOP . . .                                    *
*                                                               *
*****************************************************************
      print 1234,ncom-1,nend
 1234 format(' ncom, nend ',2i6)
c      STOP '++++++ THE END ++++++'
      end

	include 'bi_diff_60.f'
	include 'CCO2GR.f' 
	include 'fourier_filter.f' 
	include 'fourier_filter_lb.f' 
	include 'fourier_matrix.f' 
	include 'geopot_60.f' 
	include 'grwaves_RSHU.f' 
	include 'init_60.f'
	include 'irc_60.f' 
	include 'ircool_60.f' 
	include 'legandr_philb.f' 
	include 'levy.f' 
	include 'molcon_60_1.f' 
	include 'ncool_PWs_60.f' 
	include 'O3model_Berlin_60.f' 
	include 'open_60.f' 
	include 'pcoolg.f' 
	include 'plwaves_LIM.f' 
	include 'species_60.f' 
	include 'strobel_60.f' 
	include 'tendenz_60_1.f' 
	include 'tendenz_60_a.f' 
	include 'vbew.f' 


