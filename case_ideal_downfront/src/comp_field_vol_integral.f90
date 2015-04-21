       subroutine comp_field_vol_integral(field,fieldaver,jmin,jmax,    &
                                          kmin,kmax,eddyvol)

        USE header, only: NI,NJ,NK,DL,zf,dx,dy

        implicit none

        integer jmin,jmax,kmin,kmax
        double precision fieldaver
        double precision, dimension(NI,NJ,NK):: field

        integer i,j,k
        double precision dz,eddyvol   

          fieldaver=0.d0 
          do i=1,NI
           do j=jmin,jmax
            do k=kmin,kmax 
              dz = DL*(zf(i,j,k)-zf(i,j,k-1))  
              fieldaver = fieldaver + field(i,j,k)*dx*dy*dz 
            end do
           end do
          end do 


       call comp_volume(eddyvol,jmin,jmax,kmin,kmax)

       fieldaver=fieldaver/eddyvol 


       return
       END 
