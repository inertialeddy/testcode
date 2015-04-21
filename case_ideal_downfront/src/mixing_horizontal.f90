subroutine mixing_horizontal(var,vardif) 
  !     ---------------------------------------------                     
  USE header
  !     use level n                                                       
  !     computes d2u/dx2 + d2u/dy2 at the cell centers.                   
  !                                                                       
  IMPLICIT NONE 
  integer i,j,k,n,nj2
  REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: var
  REAL(kind=rc_kind) :: vardif(NI,NJ,NK)
  REAL(kind=rc_kind) :: dvardxfc(0:NI), dvardyfc(0:NJ)
  REAL(kind=rc_kind) :: vyj,uxi,fac,ubylsq,wbylsq 
  REAL(kind=rc_kind) :: Kx_l(NJ),Ky_l(NJ),e1inv,yfrac,y0,y0inv,yend,yy,KxPassTr(NJ)                                                     
  REAL(kind=rc_kind) :: sponge_dist, K_spg 

!!  integer :: selvar
  INTEGER sponge_n, sponge_s 
  LOGICAL sponge_logic 

  INTEGER OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,CHUNK,NPROC

   sponge_logic =.false. 

 if (sponge_logic) then 
  sponge_dist = 25.d0 ! (in km) 

  sponge_n=sponge_dist*1.d3/dx 
  sponge_s=sponge_n 

  K_spg=10.d0  

  write(6,*) 'sponge_n ',sponge_n, sponge_s 

 do j=1,NJ

  if (j .le. sponge_s+1) then 
    Kx_l(j) = Kx + (K_spg-Kx)*0.5*( 1.d0+dcos(PI*j/(sponge_s+1)) ) 
    Ky_l(j) = Ky + (K_spg-Ky)*0.5*( 1.d0+dcos(PI*j/(sponge_s+1)) ) 
  elseif (j .ge. NJ-sponge_n) then
    Kx_l(j) = Kx + (K_spg-Kx)*0.5*( 1.d0-dcos(PI*(j-(NJ-sponge_n))/(sponge_n)) )
    Ky_l(j) = Ky + (K_spg-Ky)*0.5*( 1.d0-dcos(PI*(j-(NJ-sponge_n))/(sponge_n)) )
  else
    Kx_l(j) = Kx 
    Ky_l(j) = Ky 
  endif 

 end do

   KxPassTr = Kx_l 

 else

  Kx_l(1:NJ) = Kx; KxPassTr= Kx_l; 
                                                                        
  Ky_l(1:NJ) = Ky 

 endif


!$OMP PARALLEL DO PRIVATE(k,i,j,dvardyfc) SHARED(vardif) COLLAPSE(2)                                                                       
  do k=1,NK 
    do i=1,NI 
      do j=1,NJ-1 
        dvardyfc(j)= 0.5*(vy(i,j)+vy(i,j+1))*(var(i,j+1,k)-var(i,j,k)) 
      end do
      dvardyfc(0)= 0.d0 
      dvardyfc(NJ)= 0.d0 
      do j=1,NJ 
        vardif(i,j,k)= Ky_l(j)*vy(i,j)*(dvardyfc(j)-dvardyfc(j-1)) 
      end do
    end do
  end do
!$OMP END PARALLEL DO



!$OMP PARALLEL DO PRIVATE(k,j,i,dvardxfc) SHARED(vardif) COLLAPSE(2)                                                                       
  do k=1,NK 
    do j=1,NJ 
      do i=1,NI-1 
        dvardxfc(i)= 0.5*(ux(i,j)+ux(i+1,j))*(var(i+1,j,k)-var(i,j,k)) 
      end do
      !     periodic-ew boundaries                                            
      dvardxfc(0)= 0.5*(ux(NI,j)+ux(1,j))*(var(1,j,k)-var(NI,j,k)) 
      dvardxfc(NI)= 0.5*(ux(NI,j)+ux(1,j))*(var(1,j,k)-var(NI,j,k)) 
      do i=1,NI 
        vardif(i,j,k)= vardif(i,j,k)+ Kx_l(j)*ux(i,j)*(dvardxfc(i)-dvardxfc(i-1))                
      end do
    end do
  end do
!$OMP END PARALLEL DO                                                                       
  fac= 1.0/(UL*LEN) 

  vardif(1:NI,1:NJ,1:NK)=fac*Jac(1:NI,1:NJ,1:NK)*vardif(1:NI,1:NJ,1:NK)

                                                                        
  return 
END subroutine
