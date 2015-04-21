subroutine ini_st 
  !     --------------------                                              
  USE header
!     initializes s as pot.density with a single vertical profile       
  IMPLICIT NONE 
  integer  i,j,k  
  integer NI2,NJ2,NK2
      parameter ( NI2=NI+2,NJ2=NJ+2,NK2=NK+2)

  INTEGER :: idTi, idvti, idvsi, rcode

#include "netcdf.inc"

!----- This part is for netcdf variables---------!
integer start(3), count(3)
!------------------------------------------------!

DATA start /1, 96, 1/
DATA count /NI2, NJ2, NK2/

!  print*, 'ini_grids'
!  CALL RANDOM_SEED()

! ! assign the density or temperature and salinity by either analytic functions or
! ! any particular hydrographic section.
!  CALL ini_grids()

!!-------- Add the latmix data to the code ---------!!
idTi = ncopn('tsideal.cdf', NCNOWRIT,rcode)
idvti = ncvid(idTi, 'Tideal', rcode)
idvsi = ncvid(idTi, 'Sideal', rcode)
call ncvgt(idTi, idvti, start, count, T, rcode)
call ncvgt(idTi, idvsi, start, count, s, rcode)

!  do j=(NJ/2-25),(NJ/2+25)  
!    do i=1,NI  
!      do k=(NK-8),NK  
!        T(i,j,k,0) = T(i,j,k,0) + ( 0.005d0 * sin( (i-NI/2)*3.14159d0/(NI/2)  ) );  
!      end do  
!    end do  
!  end do  
                                                            
  do k=0,NK+1 
    do i=1,NI 
      s(i,0,k,0)= s(i,1,k,0) 
      s(i,NJ+1,k,0)= s(i,NJ,k,0) 
      T(i,0,k,0)= T(i,1,k,0) 
      T(i,NJ+1,k,0)= T(i,NJ,k,0) 
    end do 
    do j=0,NJ+1 
      s(0,j,k,0)= s(NI,j,k,0) 
      s(NI+1,j,k,0)= s(1,j,k,0) 
      T(0,j,k,0)= T(NI,j,k,0) 
      T(NI+1,j,k,0)= T(1,j,k,0) 
    end do 
  end do 
                                                                        
                                                                        
return 
END                                           
