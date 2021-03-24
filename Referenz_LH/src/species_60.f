************************ -*- Mode: Fortran -*- **********************
** species.f ---  generates N2, O2, and O related profiles for COMMA
** Author          : Alex Pogoreltsev
** Created On      : Fri Mar  8 09:59:05 2002
** Last Modified By: Alex Pogoreltsev
** Last Modified On: Fri Mar  8 10:53:58 2002
** Update Count    : 1
** Status          : Unknown, Use with caution!
*********************************************************************
      subroutine species_60
c------ nb and kgit are arbitrary
      include 'com_main.fc'
      include 'com_namelist.fc' !FL
      include 'com_newCO2.fc'
      include 'data_newCO2.f'
      REAL ZKM(kgit)
      REAL HLOGPRF(101),HN2VMRF(101),HO2VMRF(101),
     $     HOVMRF(101)
      dzkm   = 135./(48.-0.5)
      zuntkm = 0.5*dzkm    
      SKH = 7.
         DO 10 K=1,KGIT
         ZKM(K)=ZUNTKM+FLOAT(K-1)*DZKM
   10 CONTINUE
C--  Determination the full ARRAYS
      do k=1,8
        HLOGPRF(k)=0.25*float(k-1)
        HN2VMRF(k)=HN2VMR(1)
        HO2VMRF(k)=HO2VMR(1)
        HOVMRF(k) =HOVMR(1)
c        HAVWEIF(k)=HAVWEI(1)
      enddo
      do k=9,101
         ii=k-8
        HLOGPRF(k)=HLOGPRES(ii)
        HN2VMRF(k)=HN2VMR(ii)
        HO2VMRF(k)=HO2VMR(ii)
        HOVMRF(k) =HOVMR(ii)
c        HAVWEIF(k)=HAVWEI(ii)
      enddo
      do k=1,kgit
        xx=zkm(k)/skh
          RN2VMR(k)=A18LIN(xx,HLOGPRF,HN2VMRF,1,101)
          RO2VMR(k)=A18LIN(xx,HLOGPRF,HO2VMRF,1,101)
          ROXVMR(k)=A18LIN(xx,HLOGPRF,HOVMRF, 1,101)
          RM(k)    =28.*RN2VMR(k)+32.*RO2VMR(k)+16.*ROXVMR(k)
c          RM(k)    =A18LIN(xx,HLOGPRF,HAVWEIF,1,81)
       Print*,xx, rm(k)
      end do
c      do k=49,kgit
c        xx=zkm(k)/skh
c          RN2VMR(k)=RN2VMR(48)*EXP(-(zkm(k)-zkm(48))/20.)
c          RO2VMR(k)=RO2VMR(48)*EXP(-(zkm(k)-zkm(48))/17.5)
c          ROXVMR(k)=1.-RN2VMR(k)-RO2VMR(k)
c          RM(k)    =28.*RN2VMR(k)+32.*RO2VMR(k)+16.*ROXVMR(k)
c       Print*,xx, rm(k)
c      end do
      RETURN
      END

