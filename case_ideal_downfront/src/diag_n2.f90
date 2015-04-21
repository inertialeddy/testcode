subroutine diag_n2(rhotest) 
 !     ---------------------------                                       
 USE header
 !     Calculate N2  (units: per second)                                 
                                                                       
 implicit none 
 integer i,j,k 
 REAL(kind=rc_kind) :: rdy,rdz,rpdz,depth,dz 
 REAL(kind=rc_kind) :: avezf(NK),aveN2(NK) 
 REAL(kind=rc_kind) :: rhotest(0:NI+1,0:NJ+1,0:NK+1)
                                                                       
 !  evalrho was called from rpevalgrad                                
  DLinv= 1.d0/DL 
  do k=1,NK 
    do j=1,NJ 
      do i=1,NI 

        rdy= 0.5*(  ( rhotest(i,j+1,k) - SUM(rhotest(:,j+1,k))/NI ) - ( rhotest(i,j-1,k) - SUM(rhotest(:,j-1,k))/NI ) ) * vy(i,j)/LEN
       !  PRINT*,"P2 ",LEN,vy(i,j),rdy                             
         if (k.eq.NK) then 
           rdz= (rhotest(i,j,k) -rhotest(i,j,k-1))*wz(i,j,k)*DLinv 
           rpdz= (rhotest(i,j,k) - SUM(rhotest(:,j,k))/NI -( rhotest(i,j,k-1)- SUM(rhotest(:,j,k-1))/NI))*wz(i,j,k)*DLinv 
          else if (k.eq.1) then 
           rdz= (rhotest(i,j,k+1) -rhotest(i,j,k))*wz(i,j,k)*DLinv 
           rpdz= (rhotest(i,j,k+1) - SUM(rhotest(:,j,k+1))/NI -( rhotest(i,j,k)- SUM(rhotest(:,j,k))/NI))*wz(i,j,k)*DLinv 
          else 
           rdz= 0.5*(rhotest(i,j,k+1) -rhotest(i,j,k-1))*wz(i,j,k)*DLinv 
           rpdz= 0.5* (rhotest(i,j,k+1) - SUM(rhotest(:,j,k-1))/NI -( rhotest(i,j,k-1)- SUM(rhotest(:,j,k-1))/NI))*wz(i,j,k)*DLinv 
         end if 
                                                   
         freqby(i,j,k)=(-gpr*10.d0/R0)*rdy                 
         freqbz(i,j,k)=(-gpr*10.d0/R0)*rpdz                 
         freqN2(i,j,k)=(-gpr*10.d0/R0)*rdz 
                                                                   
      end do 
    end do 
  end do 
                                                                    
  return 
END                                           

