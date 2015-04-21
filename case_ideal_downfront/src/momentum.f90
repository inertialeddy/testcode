
SUBROUTINE momentum(pcorr,step) 
  !----------------------------------------------------                   
  USE header
!  use relaxation
  INTEGER step,i,j,k 
  REAL(kind=rc_kind)                 ::  dtim,fdiv,cfcdiv,ctrdiv,edt,pcorr(maxout),temp 
  REAL(kind=rc_kind)                 :: tarray(2),tim 
!  REAL(kind=rc_kind),   dimension(  0:NI+1,0:NJ+1, 0:NK+1  )  :: rho2

!  integer ivb,ivs,ivf

  oldh= h

solverlogic = 2;   
if (solverlogic==1) then   
  write(6,*) 'mixing explicit'        
else if (solverlogic==2) then          
  write(6,*) 'mixing implicit'          
end if               

! rho_old saves the rho at the current time step, which is
! to be used by gotm only; rest everywhere rho is used.
rho_old = rho

do ivb=1,3

   if(ivb==1) then; dtim=dtf/3.d0;  tsp = dtim*1.d05; ivs=0; ivf=1; endif; 
   if(ivb==2) then; dtim=0.5d0*dtf; tsp = dtim*1.d05; ivs=1; ivf=1; endif; 
   if(ivb==3) then; dtim=dtf;       tsp = dtim*1.d05; ivs=1; ivf=0; endif; 

  CALL findzall 
  CALL sigma 
  CALL advection_and_mixing(ivs,ivf,dtim,step) 
  CALL intpol   
 if(.NOT.(use_Shchepetkin)) then 
    CALL rpevalgrad_Song(ivf) 
   else
    CALL rpevalgrad_Song(ivf) 
    CALL rpevalgrad_Sche(ivf);
    drpx(1:NI,1:NJ,1:NK)=ru4_Sche(1:NI,1:NJ,1:NK);  drpy(1:NI,1:NJ,1:NK)=rv4_Sche(1:NI,1:NJ,1:NK); 
    grpifc(0:NI,1:NJ,1:NK)=ru2_Sche(0:NI,1:NJ,1:NK); grpjfc(1:NI,0:NJ,1:NK)=rv2_Sche(1:NI,0:NJ,1:NK);
 endif
  CALL coriolis(ivs) 
  CALL srcface(ivs,step)      
  CALL hsolve(h,oldh,hdt,dtim)  
  CALL calcskfc                 
  CALL vhydro(dtim)             
  CALL cfdiv(cfcdiv)   
  CALL newcor(dtim,ivs)           
  CALL newsrc                   
  edt= EPS/dtim 
  CALL mgrid(pcorr,dtim,edt,cfcdiv)
  CALL vface(pcorr,dtim)           
  CALL vcenter(pcorr,dtim,ivf)       
  IF (fnhhy.NE.0) CALL pcorrect(pcorr)
  CALL facediv(dtim,fdiv)             
  CALL cdiv(dtim,ctrdiv,ivf)          
!  CALL evalrho(rho,ivf)   

enddo

  CALL evalrho(rho,ivf)  
!  CALL conadjust(step,0) 

  CALL diag_n2(rho) 
!  CALL sponge(0,dtim) 
 
END SUBROUTINE momentum

