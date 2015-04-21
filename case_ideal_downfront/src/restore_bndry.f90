      subroutine restore_bndry(n,step,dtim) 
!     ---------------------------------------------                     
!     restores the density below  MLD at the N-S boundaries             
!     so as to counter the effect of wind stress curl -                 
!     added on 4/21/11                                                  
      USE header

      implicit none 

      INTEGER :: i,j,k,n,nw1,nw2,ycenter,step 
      REAL(kind=rc_kind) :: rhonorth(NK),rhosouth(NK),rhomean,zdepth,   &
          dtim,restoretime(nj),trestore,fac,bndry,smean,tmean
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                        
!     Wind stress curl is in the region 1 to nw1, and nw2 to NJ         

!      nw1=0.1*nj 
!      nw2=0.9*nj 
! Restore 25 km from coast
!      bndry=25000.d0
!      nw1= int1(bndry/dx)
!      nw2= nj-nw1
! FOR FLAT-WIND ASIRI RUNS
      ycenter = 0.5*(yc(NJ)+yc(1)) 
      nw1 =   ycenter/4
      nw2 = 2*ycenter-nw1

!      write(6,*) 'nw1=',nw1,'nw2',nw2

      if (pickup_step > 0 .and. step .eq. pickup_step+1) then
!        open(unit=222,file='/home/sr71d/PSOM/github/psom-gotm/psom_v05/c&
!     ode/asiri_13N_downfront_288x576_rain2/inc/SNwall_ST.txt')
        do k=1,NK
          salt_south(k)=s(NI/2,1,k,0)
          salt_north(k)=s(NI/2,NJ,k,0)
          tsouth(k)    =T(NI/2,1,k,0)
          tnorth(k)    =T(NI/2,NJ,k,0)
!          read(222,*) salt_south(k),tsouth(k),salt_north(k),tnorth(k) 
!          write(6,*) salt_south(k),tsouth(k),salt_north(k),tnorth(k)
        end do
!      close(222)
      end if 
                                                                        
!     Restore at southern boundary, below MLD                           
!     restore time = 3 days = 3*86400/TL                                

      trestore=2.d0*86400.d0/TL 

!      do j=1,NJ
!         restoretime(j)= 10000.d0
!      end do
!      do j=1,nw1
!         restoretime(j)= float(nw1-i)/float(nw1)
!!         fac(j) = dtim/(float(nw1-i)/float(nw1))
!      end do
!      do j=nw2,nj
!         restoretime(j)= float(nj-j)/bndry
!      end do
!      do j=1,NJ
!         write(6,*) j,  restoretime(j)
!      end do
!      stop

      do k=1,NK 
         do j=1,nw1 
            smean=0.d0 
            tmean=0.d0
            do i=1,NI 
               smean= s(i,j,k,n)+ smean 
               tmean= T(i,j,k,n)+ tmean 
            end do 
            smean= smean/dble(NI) 
            tmean= tmean/dble(NI) 

            fac= dtim/trestore
            do i=1,NI 
!               zdepth=-zf(i,j,k)*DL                                    
!               if (zdepth.gt.mld(i,j)) then   !resotre the mean to rhos
                  s(i,j,k,n)= s(i,j,k,0)- fac*(smean-salt_south(k)) 
                  T(i,j,k,n)= T(i,j,k,0)- fac*(tmean-tsouth(k)) 
!               end if                                                  
            end do 
         end do 
                                                                        
         do j=nw2+1,NJ
            smean=0.d0 
            tmean=0.d0 
            do i=1,NI 
               smean= s(i,j,k,n)+ smean 
               tmean= T(i,j,k,n)+ tmean 
            end do 
            smean= smean/dble(NI) 
            tmean= tmean/dble(NI) 
            fac= dtim/trestore
            do i=1,NI 
!               zdepth=-zf(i,j,k)*DL                                    
!               if (zdepth.gt.mld(i,j)) then   !resotre the mean to rhon
                  s(i,j,k,n)= s(i,j,k,0)- fac*(smean-salt_north(k)) 
                  T(i,j,k,n)= T(i,j,k,0)- fac*(tmean-tnorth(k)) 
!               end if                                                  
            end do 
         end do 
                                                                        
      end do 
                                                                        
      return 
      END                                           
