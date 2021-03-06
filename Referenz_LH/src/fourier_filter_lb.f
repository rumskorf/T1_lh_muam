      subroutine fourier_filter_lb
      include 'com_main.fc'
      include 'com_namelist.fc' !FL
c      common/fourier/pf1(igit,igit),pf3(igit,igit),pf5 (igit,igit),
c     $               pf7(igit,igit),pf9(igit,igit),pf11(igit,igit),
c     $                              pf2(igit,igit),pf12(igit,igit)
      double precision pf1(igit,igit),pf3(igit,igit),pf5 (igit,igit),
     $                 pf7(igit,igit),pf9(igit,igit),pf11(igit,igit),
     $                 pf2(igit,igit),pf12(igit,igit), pfleg(nb,nb)
      dimension xw1(igit),xw2(igit)
c------------------------- 1 --------------------------
        do i=1,igit
               xw1(i)=G_1000(1,i)
               xw2(i)=T_1000(1,i)
        end do
        do i=1,igit
          sumw1  = 0.
          sumw2  = 0.
            do ii = 1, igit
              sumw1  = sumw1 + pf1(i,ii) * xw1(ii)
              sumw2  = sumw2 + pf1(i,ii) * xw2(ii)
            end do
          G_1000(1,i)=sumw1
          T_1000(1,i)=sumw2
        end do
c------------------------- 2 --------------------------
        do i=1,igit
               xw1(i)=G_1000(2,i)
               xw2(i)=T_1000(2,i)
        end do
        do i=1,igit
          sumw1  = 0.
          sumw2  = 0.
            do ii = 1, igit
              sumw1  = sumw1 + pf2(i,ii) * xw1(ii)
              sumw2  = sumw2 + pf2(i,ii) * xw2(ii)
            end do
          G_1000(2,i)=sumw1
          T_1000(2,i)=sumw2
        end do
c------------------------- 3 --------------------------
        do i=1,igit
               xw1(i)=G_1000(3,i)
               xw2(i)=T_1000(3,i)
        end do
        do i=1,igit
          sumw1  = 0.
          sumw2  = 0.
            do ii = 1, igit
              sumw1  = sumw1 + pf3(i,ii) * xw1(ii)
              sumw2  = sumw2 + pf3(i,ii) * xw2(ii)
            end do
          G_1000(3,i)=sumw1
          T_1000(3,i)=sumw2
        end do
c------------------------- 4 --------------------------
        do i=1,igit
               xw1(i)=G_1000(4,i)
               xw2(i)=T_1000(4,i)
        end do
        do i=1,igit
          sumw1  = 0.
          sumw2  = 0.
            do ii = 1, igit
              sumw1  = sumw1 + pf5(i,ii) * xw1(ii)
              sumw2  = sumw2 + pf5(i,ii) * xw2(ii)
            end do
          G_1000(4,i)=sumw1
          T_1000(4,i)=sumw2
        end do
c------------------------- 5 --------------------------
        do i=1,igit
               xw1(i)=G_1000(5,i)
               xw2(i)=T_1000(5,i)
        end do
        do i=1,igit
          sumw1  = 0.
          sumw2  = 0.
            do ii = 1, igit
              sumw1  = sumw1 + pf7(i,ii) * xw1(ii)
              sumw2  = sumw2 + pf7(i,ii) * xw2(ii)
            end do
          G_1000(5,i)=sumw1
          T_1000(5,i)=sumw2
        end do
c------------------------- 6 --------------------------
        do i=1,igit
               xw1(i)=G_1000(6,i)
               xw2(i)=T_1000(6,i)
        end do
        do i=1,igit
          sumw1  = 0.
          sumw2  = 0.
            do ii = 1, igit
              sumw1  = sumw1 + pf9(i,ii) * xw1(ii)
              sumw2  = sumw2 + pf9(i,ii) * xw2(ii)
            end do
          G_1000(6,i)=sumw1
          T_1000(6,i)=sumw2
        end do
c------------------ 12 harmonics --------------------------
      do j=7,nb-6
        do i=1,igit
               xw1(i)=G_1000(j,i)
               xw2(i)=T_1000(j,i)
        end do
        do i=1,igit
          sumw1  = 0.
          sumw2  = 0.
            do ii = 1, igit
              sumw1  = sumw1 + pf11(i,ii) * xw1(ii)
              sumw2  = sumw2 + pf11(i,ii) * xw2(ii)
            end do
          G_1000(j,i)=sumw1
          T_1000(j,i)=sumw2
        end do
      end do
c------------------------- nb-5 --------------------------
        do i=1,igit
               xw1(i)=G_1000(nb-5,i)
               xw2(i)=T_1000(nb-5,i)
        end do
        do i=1,igit
          sumw1  = 0.
          sumw2  = 0.
            do ii = 1, igit
              sumw1  = sumw1 + pf9(i,ii) * xw1(ii)
              sumw2  = sumw2 + pf9(i,ii) * xw2(ii)
            end do
          G_1000(nb-5,i)=sumw1
          T_1000(nb-5,i)=sumw2
        end do
c------------------------- nb-4 --------------------------
        do i=1,igit
               xw1(i)=G_1000(nb-4,i)
               xw2(i)=T_1000(nb-4,i)
        end do
        do i=1,igit
          sumw1  = 0.
          sumw2  = 0.
            do ii = 1, igit
              sumw1  = sumw1 + pf7(i,ii) * xw1(ii)
              sumw2  = sumw2 + pf7(i,ii) * xw2(ii)
            end do
          G_1000(nb-4,i)=sumw1
          T_1000(nb-4,i)=sumw2
        end do
c------------------------- nb-3 --------------------------
        do i=1,igit
               xw1(i)=G_1000(nb-3,i)
               xw2(i)=T_1000(nb-3,i)
        end do
        do i=1,igit
          sumw1  = 0.
          sumw2  = 0.
            do ii = 1, igit
              sumw1  = sumw1 + pf5(i,ii) * xw1(ii)
              sumw2  = sumw2 + pf5(i,ii) * xw2(ii)
            end do
          G_1000(nb-3,i)=sumw1
          T_1000(nb-3,i)=sumw2
        end do
c------------------------- nb-2 --------------------------
        do i=1,igit
               xw1(i)=G_1000(nb-2,i)
               xw2(i)=T_1000(nb-2,i)
        end do
        do i=1,igit
          sumw1  = 0.
          sumw2  = 0.
            do ii = 1, igit
              sumw1  = sumw1 + pf3(i,ii) * xw1(ii)
              sumw2  = sumw2 + pf3(i,ii) * xw2(ii)
            end do
          G_1000(nb-2,i)=sumw1
          T_1000(nb-2,i)=sumw2
        end do
c------------------------- nb-1 --------------------------
        do i=1,igit
               xw1(i)=G_1000(nb-1,i)
               xw2(i)=T_1000(nb-1,i)
        end do
        do i=1,igit
          sumw1  = 0.
          sumw2  = 0.
            do ii = 1, igit
              sumw1  = sumw1 + pf2(i,ii) * xw1(ii)
              sumw2  = sumw2 + pf2(i,ii) * xw2(ii)
            end do
          G_1000(nb-1,i)=sumw1
          T_1000(nb-1,i)=sumw2
        end do
c------------------------- nb --------------------------
        do i=1,igit
               xw1(i)=G_1000(nb,i)
               xw2(i)=T_1000(nb,i)
        end do
        do i=1,igit
          sumw1  = 0.
          sumw2  = 0.
            do ii = 1, igit
              sumw1  = sumw1 + pf1(i,ii) * xw1(ii)
              sumw2  = sumw2 + pf1(i,ii) * xw2(ii)
            end do
          G_1000(nb,i)=sumw1
          T_1000(nb,i)=sumw2
        end do
c----------------------------------------------------------
      return
      end

