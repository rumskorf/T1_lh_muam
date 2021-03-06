      subroutine legandr_philb
      include 'com_main.fc'
      include 'com_namelist.fc' !FL
c      common/legandr/pfleg(nb,nb)
      real philb_av(nb),philb_1(nb,igit)
      double precision pfleg(nb,nb)
c------------------------------ perturbation
      rgit  = float(igit)
      do j=1,nb
          philb_av(j)=0.     
          do i=1,igit
            philb_av(j)=philb_av(j)+philb1(j,i)/rgit
          end do
           do i=1,igit
            philb_1(j,i)=philb1(j,i)-philb_av(j)
           end do
      end do
c--------------------------------- smoothing       
        do i=1,igit
          do j=1,nb
            sumphi  = 0.
            do jj = 1,nb
              sumphi  = sumphi + pfleg(j,jj)*philb_1(jj,i)
            end do
            philb1(j,i) = sumphi
          end do
        end do
c-------------------------------------------
      return
      end

