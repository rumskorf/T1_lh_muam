      subroutine heating_LH
      include 'com_main.fc'
      include 'com_namelist.fc'

      vrAD=PI/180.
      PI2=2.*PI

c##################################### calculation of the UT
c---------------------------- Alex Pogoreltsev February 2016
        rtime=float(ncom)/384.
        time=(rtime-INT(rtime))*24.+12.
        if(time.GE.24.) time=time-24.

c------------ without the diurnal variability of LH heating
c             until ncom < nphi  
c----------------------------------------------------------------------
      nfor=max0(ncom-nphi,0) 
      xsec   = 3600./float(ntime)          ! seconds per step
      tdfor  = float(nfor)*xsec/30./86400. 
c--------------------- efficiency of the generation of tides after nphi
      t_tides= 1.- exp(-tdfor) 
        DO 110 k=1,10
        do 111 j=1,nb
        do 112 i=1,igit
        rlon=float(i-1)*5.625
c------------------------------- LH ----------------------------------
	Heat_LH(i,j,k)=Heat_LH0(i,j,k)
     *  +t_tides*(
     *   Alhm0_24(j,k)*cos(PI2*time/24.-Flhm0_24(j,k)*VRAD)
     *  +Alhm0_12(j,k)*cos(PI2*time/12.-Flhm0_12(j,k)*VRAD)
c------------------------------- m = 1 --------------------------------------
     *  +Alhm1w24(j,k)*cos(PI2*time/24.+(rlon-Flhm1w24(j,k))*VRAD)
     *  +Alhm1w12(j,k)*cos(PI2*time/12.+(rlon-Flhm1w12(j,k))*VRAD)
     *  +Alhm1e24(j,k)*cos(-PI2*time/24.+(rlon-Flhm1e24(j,k))*VRAD)
     *  +Alhm1e12(j,k)*cos(-PI2*time/12.+(rlon-Flhm1e12(j,k))*VRAD)
c------------------------------- m = 2 --------------------------------------
     *  +Alhm2w24(j,k)*cos(PI2*time/24.+2.*(rlon-Flhm2w24(j,k))*VRAD)
     *  +Alhm2w12(j,k)*cos(PI2*time/12.+2.*(rlon-Flhm2w12(j,k))*VRAD)
     *  +Alhm2e24(j,k)*cos(-PI2*time/24.+2.*(rlon-Flhm2e24(j,k))*VRAD)
     *  +Alhm2e12(j,k)*cos(-PI2*time/12.+2.*(rlon-Flhm2e12(j,k))*VRAD)
c------------------------------- m = 3 --------------------------------------
     *  +Alhm3w24(j,k)*cos(PI2*time/24.+3.*(rlon-Flhm3w24(j,k))*VRAD)
     *  +Alhm3w12(j,k)*cos(PI2*time/12.+3.*(rlon-Flhm3w12(j,k))*VRAD)
     *  +Alhm3e24(j,k)*cos(-PI2*time/24.+3.*(rlon-Flhm3e24(j,k))*VRAD)
     *  +Alhm3e12(j,k)*cos(-PI2*time/12.+3.*(rlon-Flhm3e12(j,k))*VRAD)
c------------------------------- m = 4 --------------------------------------
     *  +Alhm4w24(j,k)*cos(PI2*time/24.+4.*(rlon-Flhm4w24(j,k))*VRAD)
     *  +Alhm4w12(j,k)*cos(PI2*time/12.+4.*(rlon-Flhm4w12(j,k))*VRAD)
     *  +Alhm4e24(j,k)*cos(-PI2*time/24.+4.*(rlon-Flhm4e24(j,k))*VRAD)
     *  +Alhm4e12(j,k)*cos(-PI2*time/12.+4.*(rlon-Flhm4e12(j,k))*VRAD))
      IF(Heat_LH(i,j,k).le.0.) Heat_LH(i,j,k) = 0.
c--------------------------------------------------------------------
 112  CONTINUE
 111  CONTINUE
 110  CONTINUE
      RETURN 
      END 
