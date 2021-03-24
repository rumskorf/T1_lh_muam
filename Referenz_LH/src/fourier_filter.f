      subroutine fourier_filter
      include 'com_main.fc'
      include 'com_namelist.fc' !FL
!FL       common/fourier/pf1(igit,igit),pf3(igit,igit),pf5 (igit,igit),
!FL      $               pf7(igit,igit),pf9(igit,igit),pf11(igit,igit),
!FL      $                              pf2(igit,igit),pf12(igit,igit)
      real xu1(igit),xv1(igit),xt1(igit)
      double precision pf1(igit,igit),pf3(igit,igit),pf5 (igit,igit),
     $                 pf7(igit,igit),pf9(igit,igit),pf11(igit,igit),
     $                 pf2(igit,igit),pf12(igit,igit), pfleg(nb,nb)
c------------------------- 1 --------------------------
      do 100 k=1,kgit
        do i=1,igit
               xu1(i)=an1(1,k,i,1)
               xv1(i)=an1(1,k,i,2)
               xt1(i)=an1(1,k,i,3)
        end do
        do i=1,igit
          sumu  = 0.
          sumv  = 0.
          sumt  = 0.
            do ii = 1, igit
c             if(nsec.le.30*86400) then 
c              sumu  = sumu + pf1(i,ii) * xu1(ii)
c              sumv  = sumv + pf1(i,ii) * xv1(ii)
c              sumt  = sumt + pf1(i,ii) * xt1(ii)
c             else    
              sumu  = sumu + pf2(i,ii) * xu1(ii)
              sumv  = sumv + pf2(i,ii) * xv1(ii)
              sumt  = sumt + pf2(i,ii) * xt1(ii)
c             end if
            end do
          an1(1,k,i,1)=sumu
          an1(1,k,i,2)=sumv
          an1(1,k,i,3)=sumt
        end do
c------------------------- 2 --------------------------
        do i=1,igit
               xu1(i)=an1(2,k,i,1)
               xv1(i)=an1(2,k,i,2)
               xt1(i)=an1(2,k,i,3)
        end do
        do i=1,igit
          sumu  = 0.
          sumv  = 0.
          sumt  = 0.
            do ii = 1, igit
c            if(nsec.le.30*86400) then 
c              sumu  = sumu + pf2(i,ii) * xu1(ii)
c              sumv  = sumv + pf2(i,ii) * xv1(ii)
c              sumt  = sumt + pf2(i,ii) * xt1(ii)
c             else    
              sumu  = sumu + pf3(i,ii) * xu1(ii)
              sumv  = sumv + pf3(i,ii) * xv1(ii)
              sumt  = sumt + pf3(i,ii) * xt1(ii)
c             end if
            end do
          an1(2,k,i,1)=sumu
          an1(2,k,i,2)=sumv
          an1(2,k,i,3)=sumt
        end do
c------------------------- 3 --------------------------
        do i=1,igit
               xu1(i)=an1(3,k,i,1)
               xv1(i)=an1(3,k,i,2)
               xt1(i)=an1(3,k,i,3)
        end do
        do i=1,igit
          sumu  = 0.
          sumv  = 0.
          sumt  = 0.
            do ii = 1, igit
              sumu  = sumu + pf5(i,ii) * xu1(ii)
              sumv  = sumv + pf5(i,ii) * xv1(ii)
              sumt  = sumt + pf5(i,ii) * xt1(ii)
            end do
          an1(3,k,i,1)=sumu
          an1(3,k,i,2)=sumv
          an1(3,k,i,3)=sumt
        end do
c------------------------- 4 --------------------------
        do i=1,igit
               xu1(i)=an1(4,k,i,1)
               xv1(i)=an1(4,k,i,2)
               xt1(i)=an1(4,k,i,3)
        end do
        do i=1,igit
          sumu  = 0.
          sumv  = 0.
          sumt  = 0.
            do ii = 1, igit
              sumu  = sumu + pf7(i,ii) * xu1(ii)
              sumv  = sumv + pf7(i,ii) * xv1(ii)
              sumt  = sumt + pf7(i,ii) * xt1(ii)
            end do
          an1(4,k,i,1)=sumu
          an1(4,k,i,2)=sumv
          an1(4,k,i,3)=sumt
        end do
c------------------------- 5 --------------------------
        do i=1,igit
               xu1(i)=an1(5,k,i,1)
               xv1(i)=an1(5,k,i,2)
               xt1(i)=an1(5,k,i,3)
        end do
        do i=1,igit
          sumu  = 0.
          sumv  = 0.
          sumt  = 0.
            do ii = 1, igit
              sumu  = sumu + pf9(i,ii) * xu1(ii)
              sumv  = sumv + pf9(i,ii) * xv1(ii)
              sumt  = sumt + pf9(i,ii) * xt1(ii)
            end do
          an1(5,k,i,1)=sumu
          an1(5,k,i,2)=sumv
          an1(5,k,i,3)=sumt
        end do
c------------------------- 6 --------------------------
        do i=1,igit
               xu1(i)=an1(6,k,i,1)
               xv1(i)=an1(6,k,i,2)
               xt1(i)=an1(6,k,i,3)
        end do
        do i=1,igit
          sumu  = 0.
          sumv  = 0.
          sumt  = 0.
            do ii = 1, igit
              sumu  = sumu + pf11(i,ii) * xu1(ii)
              sumv  = sumv + pf11(i,ii) * xv1(ii)
              sumt  = sumt + pf11(i,ii) * xt1(ii)
            end do
          an1(6,k,i,1)=sumu
          an1(6,k,i,2)=sumv
          an1(6,k,i,3)=sumt
        end do
c------------------ 12 harmonics --------------------------
      do j=7,nb-6
        do i=1,igit
               xu1(i)=an1(j,k,i,1)
               xv1(i)=an1(j,k,i,2)
               xt1(i)=an1(j,k,i,3)
        end do
        do i=1,igit
          sumu  = 0.
          sumv  = 0.
          sumt  = 0.
            do ii = 1, igit
              sumu  = sumu + pf12(i,ii) * xu1(ii)
              sumv  = sumv + pf12(i,ii) * xv1(ii)
              sumt  = sumt + pf12(i,ii) * xt1(ii)
            end do
          an1(j,k,i,1)=sumu
          an1(j,k,i,2)=sumv
          an1(j,k,i,3)=sumt
        end do
      end do
c------------------------- nb-5 --------------------------
        do i=1,igit
               xu1(i)=an1(nb-5,k,i,1)
               xv1(i)=an1(nb-5,k,i,2)
               xt1(i)=an1(nb-5,k,i,3)
        end do
        do i=1,igit
          sumu  = 0.
          sumv  = 0.
          sumt  = 0.
            do ii = 1, igit
              sumu  = sumu + pf11(i,ii) * xu1(ii)
              sumv  = sumv + pf11(i,ii) * xv1(ii)
              sumt  = sumt + pf11(i,ii) * xt1(ii)
            end do
          an1(nb-5,k,i,1)=sumu
          an1(nb-5,k,i,2)=sumv
          an1(nb-5,k,i,3)=sumt
        end do
c------------------------- nb-4 --------------------------
        do i=1,igit
               xu1(i)=an1(nb-4,k,i,1)
               xv1(i)=an1(nb-4,k,i,2)
               xt1(i)=an1(nb-4,k,i,3)
        end do
        do i=1,igit
          sumu  = 0.
          sumv  = 0.
          sumt  = 0.
            do ii = 1, igit
              sumu  = sumu + pf9(i,ii) * xu1(ii)
              sumv  = sumv + pf9(i,ii) * xv1(ii)
              sumt  = sumt + pf9(i,ii) * xt1(ii)
            end do
          an1(nb-4,k,i,1)=sumu
          an1(nb-4,k,i,2)=sumv
          an1(nb-4,k,i,3)=sumt
        end do
c------------------------- nb-3 --------------------------
        do i=1,igit
               xu1(i)=an1(nb-3,k,i,1)
               xv1(i)=an1(nb-3,k,i,2)
               xt1(i)=an1(nb-3,k,i,3)
        end do
        do i=1,igit
          sumu  = 0.
          sumv  = 0.
          sumt  = 0.
            do ii = 1, igit
              sumu  = sumu + pf7(i,ii) * xu1(ii)
              sumv  = sumv + pf7(i,ii) * xv1(ii)
              sumt  = sumt + pf7(i,ii) * xt1(ii)
            end do
          an1(nb-3,k,i,1)=sumu
          an1(nb-3,k,i,2)=sumv
          an1(nb-3,k,i,3)=sumt
        end do
c------------------------- nb-2 --------------------------
        do i=1,igit
               xu1(i)=an1(nb-2,k,i,1)
               xv1(i)=an1(nb-2,k,i,2)
               xt1(i)=an1(nb-2,k,i,3)
        end do
        do i=1,igit
          sumu  = 0.
          sumv  = 0.
          sumt  = 0.
            do ii = 1, igit
              sumu  = sumu + pf5(i,ii) * xu1(ii)
              sumv  = sumv + pf5(i,ii) * xv1(ii)
              sumt  = sumt + pf5(i,ii) * xt1(ii)
            end do
          an1(nb-2,k,i,1)=sumu
          an1(nb-2,k,i,2)=sumv
          an1(nb-2,k,i,3)=sumt
        end do
c------------------------- nb-1 --------------------------
        do i=1,igit
               xu1(i)=an1(nb-1,k,i,1)
               xv1(i)=an1(nb-1,k,i,2)
               xt1(i)=an1(nb-1,k,i,3)
        end do
        do i=1,igit
          sumu  = 0.
          sumv  = 0.
          sumt  = 0.
            do ii = 1, igit
c            if(nsec.le.30*86400) then 
c              sumu  = sumu + pf2(i,ii) * xu1(ii)
c              sumv  = sumv + pf2(i,ii) * xv1(ii)
c              sumt  = sumt + pf2(i,ii) * xt1(ii)
c             else    
              sumu  = sumu + pf3(i,ii) * xu1(ii)
              sumv  = sumv + pf3(i,ii) * xv1(ii)
              sumt  = sumt + pf3(i,ii) * xt1(ii)
c             end if 
            end do
          an1(nb-1,k,i,1)=sumu
          an1(nb-1,k,i,2)=sumv
          an1(nb-1,k,i,3)=sumt
        end do
c------------------------- nb --------------------------
        do i=1,igit
               xu1(i)=an1(nb,k,i,1)
               xv1(i)=an1(nb,k,i,2)
               xt1(i)=an1(nb,k,i,3)
        end do
        do i=1,igit
          sumu  = 0.
          sumv  = 0.
          sumt  = 0.
            do ii = 1, igit
c             if(nsec.le.30*86400) then 
c              sumu  = sumu + pf1(i,ii) * xu1(ii)
c              sumv  = sumv + pf1(i,ii) * xv1(ii)
c              sumt  = sumt + pf1(i,ii) * xt1(ii)
c             else    
              sumu  = sumu + pf2(i,ii) * xu1(ii)
              sumv  = sumv + pf2(i,ii) * xv1(ii)
              sumt  = sumt + pf2(i,ii) * xt1(ii)
c             end if
            end do
          an1(nb,k,i,1)=sumu
          an1(nb,k,i,2)=sumv
          an1(nb,k,i,3)=sumt
      end do
  100 continue
c----------------------------------------------------------
      return
      end

