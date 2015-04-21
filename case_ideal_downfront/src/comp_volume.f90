       subroutine comp_volume(volume,jmin,jmax,kmin,kmax)

        USE header, only:NI,NJ,NK,DL,dx,dy,zf
        implicit none

        integer jmin,jmax,kmin,kmax
        double precision volume

        integer i,j,k 
        double precision dz 


          volume=0.d0 
          do i=1,NI
           do j=jmin,jmax
            do k=kmin,kmax 
              dz = DL*(zf(i,j,k)-zf(i,j,k-1))  
              volume = volume + dx*dy*dz 
            end do
           end do
          end do 
!         write(6,*) volume 

       return
       END  
