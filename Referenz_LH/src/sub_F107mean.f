
C       subroutine sub_F107meam(jahr,monat,mittelwert)
C       f2py --fcompiler=gfortran -m -c sub_F107meam sub_F107meam.f

      subroutine sub_F107mean(JAHR,NDEK,f107mean)
      implicit none
!      character*3 im
      integer i,j,indexx,JAHR,datum,datumtest,mm,ij,NDEK
      real  f107,summe,f107mean

!Cf2py intent(in) ij
!Cf2py intent(in) im
!Cf2py intent(out) mittelwert
      
      open(10,File="DAILYPLT.OBS")
      open(20,File="f107mean.txt")
      
110   format (i8,4x,f6.1)
120   format (i3)

      
      if(NDEK.le. 31)                 nmon=1
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

!       if     (im.eq.'Jan') then
!        mm=1
!       elseif (im.eq.'Feb') then
!        mm=2
!       elseif (im.eq.'Mar') then
!        mm=3
!       elseif (im.eq.'Apr') then
!        mm=4
!       elseif (im.eq.'May') then
!        mm=5
!       elseif (im.eq.'Jun') then
!        mm=6
!       elseif (im.eq.'Jul') then
!        mm=7
!       elseif (im.eq.'Aug') then
!        mm=8
!       elseif (im.eq.'Sep') then
!        mm=9
!       elseif (im.eq.'Oct') then
!        mm=10
!       elseif (im.eq.'Nov') then
!        mm=11
!       elseif (im.eq.'Dec') then
!        mm=12
!       endif
      vergleich=ij*100+nmon
      summe=0.
      j=0
      do i=1,200000
         read(10,110)Datum,f107
         datumtest =datum/100
         if (datumtest.eq.vergleich) then
            if (f107.eq.0.)then
               indexx=0
            else
               indexx=1
            endif
            j=j+indexx
            summe=summe+f107
            f107mean=summe/j
            rewind(20)
            write(20,120)f107mean
         elseif (datumtest.gt.vergleich) then
            return
         endif
      enddo
      close(10)
      close(20)
      END

