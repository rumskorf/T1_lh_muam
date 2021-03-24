      subroutine open_60
      include 'com_main.fc'
      include 'com_namelist.fc'
      
      character*4 mlu
      integer yr,mo,nmon,monat,y,m
      real c, rmse
c.....DATA FILES
!       open (15,status='old',form='formatted',file='Run32_tst_Alex_data')
!       read(15,'(a)') f2_old
!       read(15,'(a)') f2_new
! c--------------- a possibility to write F2 file as formatted
! c      read(15,'(a)') f2_nfm
!       read(15,'(a)') modul 
!       read(15,'(a)') an1dat
!       print*,'Data files: f2_old,f2_new,modul,an1dat'
!       print*,f2_old,f2_new,modul,an1dat
! c.....OPEN-I: 
!       open (17,status='old',form='formatted',file=modul)

      OPEN(UNIT=10,FILE='muam_mod.txt')

      READ(10,RUN)
      READ(10,INI)
      READ(10,FDX)
      READ(10,LOB)
      READ(10,FOR)
      READ(10,PWS)
      READ(10,PWX)
      READ(10,PWP)
      READ(10,PWF)
      READ(10,GRW) 
      READ(10,SOL) 

      CLOSE(10)
cfl	kommt in init nochmal      
!       print*,NPHI,NEND,NPRINT,NOUT,NTIME,KRET,NSUN,NDEK
!       print*,FOLD,FNEW,FO3
!       print*,FUVT,FPHI,FVER,FGFU,FGFV,FGAU,FGAV
!       print*,TLB,T00,T11,T22,T33,G00,G11,G22,G33
!       print*,MFORC,MTURB,MCOOL,MUVIR,MIOND,MMCON,MSURF
!       print*,I_2d,I_SK,I_KW,I_FK,I_NK,I_10,I_11,I_12,I_13,I_21,I_23
!       print*,Xm2d,XmSK,XmFK,XmKW,XmNK,Xm10,Xm11,Xm12,Xm13,Xm21,Xm23
!       print*,Per_2d,Per_SK,Per_FK,Per_KW,Per_NK,Per_10,Per_11,Per_12,
!      #       Per_13,Per_21,Per_23
!       print*,fi0_2d,fi0_SK,fi0_FK,fi0_KW,fi0_NK,fi0_10,fi0_11,fi0_12,
!      #       fi0_13,fi0_21,fi0_23
!       print*,XL,WAMPL,CPHASE,TETA

      
c.....OPEN-II: prognostic-fields
C------     DIRECT ACCESS FORMAT
      open (72,status='unknown',form='unformatted',file=FUVT,
     +         access='direct', recl=1*nb*kgit*igit*3)
      open (24,status='unknown',form='unformatted',file=FPHI, 
     +         access='direct',recl=1*nb*igit*kgit)
      open (25,status='unknown',form='unformatted',file=FVER, 
     +	       access='direct',recl=1*nb*igit*kgit)

C------     GW momentum flux(u,v) and acceleration ac(u,v) ! by Rliu

      open (26,status='unknown',form='unformatted',file=FGFU,
     +         access='direct',recl=1*nb*igit*kgit)
      open (27,status='unknown',form='unformatted',file=FGFV,
     +         access='direct',recl=1*nb*igit*kgit)
      open (28,status='unknown',form='unformatted',file=FGAU,
     +         access='direct',recl=1*nb*igit*kgit)
      open (29,status='unknown',form='unformatted',file=FGAV,
     +         access='direct',recl=1*nb*igit*kgit)
      open (30,status='unknown',form='unformatted',file=FGAT,
     +         access='direct',recl=1*nb*igit*kgit)
     
!       open (31,status='unknown',form='unformatted',file='nudging.dx',
!      +	       access='direct',recl=1*nb*igit*kgit)	
! c------     ASCII				      
!       open (70,status='unknown',form='formatted',
!      $                         file='zon_tst.dat')
! c      open (71,status='unknown',form='formatted',
! c     $                         file='zon_srf.dat')
! c      open (73,status='unknown',form='formatted',
! c     $                         file='tmp_srf.dat')
c------     STANDARD BINARIES
      open ( 12,status='old',form='unformatted',file=FOLD)
      open (  2,status='replace',form='unformatted',file=FNEW)
c      open (113,status='new',form='formatted  ',file=f2_nfm)
C???????????????????????????????? ------ what is the reason of defaul setting 
C.....MODULE SWITCHES (default settings)
cfl 	werden nicht gebraucht	
!       mref    =10
!       mhemi   =0
!       mfilt   =1
!       mforc   =1 
!       mturb   =1
!       mcool   =2
!       muvir   =1
!       miond   =2
!       mmcon   =1
!       msurf   =0
c---------------------------     read new settings from module file -->mod_name
!       read (17,1140) mforc,mturb, mcool, muvir, miond, mmcon, msurf
!  1140 format(/,7(/,18x, i4))
cfl      print*,mforc,mturb, mcool, muvir, miond, mmcon, msurf

c.....INITIALIZE DYNAMIC FIELDS -->fort.2
      read (12) an0
      read (12) an1
      read (12) an2
      
c      read (12) diff3
      read (12) philb0
      read (12) philb1
c      read (12) difflb
      read (12) fphi_0
c-------------------------------------------     control parameters
      read (12) ken1 
c      print*, FNEW,' OK!'
      do 273 iii=1,64
c      write(73,173) iii,ken1(iii)
      print 173, iii,ken1(iii)
  173 format(1x,' ken',i2,'=',i8)
  273 continue
c--------------------------------------------     overwrite defaults 
cfl in muam_mod
!       nstart  =ncom
!       nend    =7680
!       nprint  =5
!       nout    =192
!       ntime   =8
!      kret    =0  ! -1,0,+1 kein sichern, letzten schritt, jeden tag
      ken1(15)=igit
      ken1(16)=kgit
      ken1(17)=nb
! c-----------------     read new settings from module file -->mod_name
!       read (17,1150) 
!      +  nphi,nend, nprint,nout, ntime, kret, nsun, ndek
!  1150 format( 8(/,18x, i7))

cfl	print*,nphi,nend, nprint,nout, ntime, kret, nsun, ndek

c-----------------------     CHECK OF CONTROL PARAMETERS AND SWITCHES
cfl kommt nochmal in init_60
!       print 1040, mforc,mturb,mcool,muvir,miond,mmcon,msurf
!  1040 format(//'  processes:'/20x,
!      +     ' mforc =',i6,/20x,
!      +     ' mturb =',i6,/20x,
!      +     ' mcool =',i6,/20x,
!      +     ' muvir =',i6,/20x,
!      +     ' miond =',i6,/20x,
!      +     ' mmcon =',i6,/20x,
!      +     ' msurf =',i6,/)
!       print 1042,
!      $    ncom,nphi,nend,nprint,nout,ntime,nsec,kret,nsun,ndek,nsrc
!  1042 format(//'  controls:'/20x,
!      +     ' nstart (ncom)  =',i9,/20x,
!      +     ' nphi  =',i9,/20x,
!      +     ' nend  =',i9,/20x,
!      +     ' nprint=',i9,/20x,
!      +     ' nout  =',i9,/20x,
!      +     ' ntime =',i9,/20x,
!      +     ' nsec  =',i9,/20x,
!      +     ' kret  =',i9,/20x,
!      +     ' nsun  =',i9,/20x,
!      +     ' ndek  =',i9,/20x,
!      +     ' nsrc  =',i9,/)
C.....ABSORBER CONCENTRATIONS
      if(NDEK.le.31)                  nmon=1
      if(NDEK.ge. 32.and.NDEK.le. 59) nmon=2
      if(NDEK.ge. 60.and.NDEK.le. 90) nmon=3
      if(NDEK.ge. 91.and.NDEK.le.120) nmon=4
      if(NDEK.ge.121.and.NDEK.le.151) nmon=5
      if(NDEK.ge.152.and.NDEK.le.181) nmon=6
      if(NDEK.ge.182.and.NDEK.le.212) nmon=7
      if(NDEK.ge.213.and.NDEK.le.243) nmon=8
      if(NDEK.ge.244.and.NDEK.le.273) nmon=9
      if(NDEK.ge.274.and.NDEK.le.304) nmon=10
      if(NDEK.ge.305.and.NDEK.le.334) nmon=11
      if(NDEK.ge.335.and.NDEK.le.365) nmon=12
      
      monat=nmon

c	included by FL 2014 for exact co2 data for a certain year
      open(1,file=FCO2)
      read(1,'(70(/))')
      do y=1974,2013
	do m=1,12
	  read(1,'(A4,I4,X,I2,2(F11.3))') mlu,yr,mo,c,rmse	  
	  if (yr.eq.YEAR.and.mo.eq.nmon) then
	    print*,yr,YEAR,mo,nmon
	    co2yr=c
	    print*,co2yr
	    if (co2.lt.-999.) print*,'No CO2 data!'
	  endif	  
	enddo
      enddo      
      close(1)
      
      CALL O3model_Berlin_60(nmon)
      do k=1,kgit
	print*,'o3=',so3(23,k),'  o3ppmv=',ozppmv(23,k) 
      enddo
C.....OTHER PROFILES
c
      call species_60
       return
      end

c      include 'species_60.f'
c      include 'O3model_Berlin_60.f'
c      include 'CCO2GR.f'


