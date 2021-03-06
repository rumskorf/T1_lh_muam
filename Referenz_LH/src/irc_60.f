      subroutine irc_60  

      parameter (NGIT=118,k_48=48)
 
c    interpolation fehlender werte fuer KGIT-level version
c
c------------------------------------------- KGIT --> K_48!!!!!!!!!!!
c                                            A.Pogoreltsev, Feb. 2006  

      include 'com_main.fc'
      include 'com_namelist.fc' !FL

c     ---------------- neue CO2-Parametrisierung -------------
      include 'com_newCO2.fc'
c      include 'data_newCO2.f'
c     --------------------------------------------------------
      real co2pmv,co2fak
c      dimension x(kgit+1),y(kgit+1),y2(kgit+1),xi(kgit+1)
      dimension xi(k_48+1)
c???????????????
c- N will be used in spline for rnn(kgit,j)???????????????? - remove???
c      N = KGIT + 1
C???????????????

csl      co2pmv=360.+360.      ! Verdopplung der CO2-Konzentration
csl   normalkonditionen CO2: co2pmv=360.

cfl      co2pmv=360. 
cfl      co2pmv=360.*CO2
      co2pmv=co2yr
      co2fak=co2pmv/360.
      h     =7000.
      Rh    =287.04/h  ! R/H
      Ch    =.2857/h   ! kappa/h
      R     =287.04
      dich0 =1.273e-3
      rgit  =float(igit)
      z_48 =135000.  !!!fest vorgegeben
      dzi   = z_48/(float(ngit)-0.5)   !!! ngit
      zunteni=dzi*0.5

c.....zi,roi,dichti (KGIT-->118) -->ircool_118
      do k=1,ngit   ! 0.57 < zi < 135.
         zi(k)=zunteni/1000.+float(k-1)*dzi/1000. !log-press.hoehe/km
         roi(k)=exp(-zi(k)*1000./h)		       
         dichti(k)=dich0*roi(k)
      enddo
c-------------------------------------     xi fuer h2o, cox
      do k=1,k_48+1   ! -1.44 < xi < 135.
         xi(k)=zunten/1000.+float(k-2)*dz/1000. !log-press.hoehe/km
      enddo

C     ....h2o,O3
         zh2o(1)=0.
         zint(1)=0.  
      do k=1,72  ! 0.0 < zh2o < 40.47,  0.0 < zint < 124.25
         zh2o(k)=zh2o(1)+float(k-1)*zi(1)
         zint(k)=zint(1)+float(k-1)*0.25*skh
      enddo
       print*,'step1'
c?????????????????????? - remove bellow?
c.....rni -->ircool_old
c     N=points, y1=1.abl bei x1, yN=1.abl bei xN, y2=2.abl
c      do j=1,3
c         yN=rnn(kgit,j)-rnn(kgit-1,j)
c         n0=0
c         do i=1,kgit+1
c            i1=i-1
c            if(i1.eq.0) i1=1    !level 0 ergaenzt
c            x(i)=zunten+float(i-2)*dz !log-press.-hoehe/m
c            y(i)=rnn(i1,j)
c            if(y(i).eq.0.) n0=n0+1
c         enddo
c         do i=1,n0
c            y(i)=y(n0+1)
c         enddo
c         y1=y(2)-y(1)
c         call spline(x,y,N,y1/dz,yN/dz,y2)c
c
c------------------------     altes feld=x,y; interpol werte=aa,bb
c         klo=1
c         khi=2
c      do k=1,ngit
c         aa=zi(k)*1000.
c------     Verschieben des Intervalls, so dass x(klo) <= a < x2(khi)
c         if (aa.gt.x(khi)) then
c            klo=khi
c            khi=khi+1
c         endif
c         call splint(x,y,y2,N,aa,bb,klo,khi)
c         rni(k,j)=bb
c      enddo
c      enddo
c?????????????????????????????? - remove above? 
       print*,'step2'
c.....index-arrays 118-->KGIT
         klo=1
         khi=2
      do 20 k=1,k_48
         aa=z(k)
c----     Verschieben des Intervalls, so dass zi(klo) < aa <= zi(khi)
         do kk=1,ngit
            if (aa.gt.zi(khi)) then
               klo=khi
               khi=khi+1
            else
               goto 500
            endif
         enddo
 500     continue
         ipVar118(k)=klo
 20      continue

c----     rueckwaerts KGIT-->118
         klo=1
         khi=2
      do 40 k=1,ngit
         aa=zi(k)
c----     Verschieben des Intervalls, so dass z(klo) < aa <= z(khi)
         do kk=1,k_48
            if (aa.gt.xi(khi)) then
               klo=khi
               khi=khi+1
            else
               goto 505
            endif
         enddo
 505     continue
         ip118var(k)=klo
 40   continue

c     +++++++++++++ I N D E X ++++++++++++++++++++++

c.....index-arrays 72(O3,h2o)-->KGIT
         klo=1
         khi=2
         kmx=int(124.25*1000./dz+.5)
         do 60 k=1,kmx
         aa=z(k)
c--------     Verschieben des Intervalls, so dass zi(klo) < aa <= zi(khi)
         do kk=1,70
            if (aa.gt.zint(khi)) then
               klo=khi
               khi=khi+1
            else
               goto 501
            endif
         enddo
 501     continue
         ipCOXvar(k)=klo
 60      continue
c---------------------------------------------     rueckwaerts KGIT-->72
         klo=1
         khi=2
      do 80 k=1,72
         aa=zint(k)
c---------     Verschieben des Intervalls, so dass z(klo) < aa <= z(khi)
         do kk=1,k_48
            if (aa.gt.xi(khi)) then
               klo=khi
               khi=khi+1
            else
               goto 605
            endif
         enddo
 605     continue
         ipolCOX(k)=klo
 80   continue
         klo=1
         khi=2
      do 90 k=1,72
         aa=zH2O(k)
c-----------     Verschieben des Intervalls, so dass z(klo) < aa <= z(khi)
         do kk=1,k_48
            if (aa.gt.xi(khi)) then
               klo=khi
               khi=khi+1
            else
               goto 705
            endif
         enddo
 705     continue
         ipolH2O(k)=klo
   90 continue
      print*,'step3'
c
c ---------------- Initialisation newCO2-param. (Ogibalov 99)  ----------
c
c     CO2 Volumenmischungsverhaeltnis an aktuelle Konzentration anpassen
      do i=1,73
         HCO2VMR(i)=HCO2VMR(i)*co2fak
      enddo

      do ii=43,59
        i=ii-42 ! i=1: x=12.5, -> ii=1: x=2
        pxnlte(i)=HLOGPRES(ii)        
        PCO2(i)=HCO2VMR(ii)
        PAM(i)=HAVWEI(ii)
      enddo

c *****  CO2 concentration below x=12.5                     *****
      r360CO2=0.3600E-03
      rCO2=co2pmv*1.e-6   ! 0.3600E-03  variablenname angepasst (lange oct 99)
                          ! modern value of CO2 vmr below x=12.5
                          ! should be changed if necessary

c *****  CO2 column amount (cm-2) above the level of x=16.5 *****
      utop360=1.052353E+14
      UTOP=rCO2/r360CO2*utop360

c ***************************************************************
C
c to  determine the coefficients in both LTE and NLTE regions for a certain 
c CO2 concentration.
c Should be called if CO2 concentration between x=0-16.5 has been changed
C
      print*,'step4'
      if(IREP.eq.0) then
      call cco2gr(PCO2,PAM,UTOP,rCO2)
      IREP = 1
      end if
      print*,'step5'
C
c------------------------------------------------------------------------------
      return
      end
C------------------------------------------------------------------------------
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      PARAMETER (NMAX=120)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y,KLO,KHI)
      DIMENSION XA(N),YA(N),Y2A(N)
c      KLO=1
c      KHI=N
c1     IF (KHI-KLO.GT.1) THEN
c        K=(KHI+KLO)/2
c        IF(XA(K).GT.X)THEN
c          KHI=K
c        ELSE
c          KLO=K
c        ENDIF
c      GOTO 1
c      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) PAUSE 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END








