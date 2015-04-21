
subroutine fluxes(step)
  !     ---------------------------------------------                     
  USE header

  !     Wind Stress can be specified here                                 
  !     Kz*du/dz = -Tx/rho,   Kz*dv/dz = -Ty/rho                          
  !     use level m                                                       
  !     computes d2s/dz2 at the cell centers.                             
  !     fricu and fricv are u,v diffusion terms (in header)               
  !     that need to be saved (if m.eq.0) for n2budget)                   

  implicit none 

  INTEGER ::i,j,k,m,step,ycenter,ymax,rst_step,step_1f,rst_tau,stress_len
  REAL(kind=rc_kind)  ::  stressmax_x, stressmax_y,timday,time_ramp,diff,Tf,f_BoB 
! REAL(kind=rc_kind)  ::  J_lambda1, J_lambda2, J_A
  logical             ::  use_winds,latewinds
!---------------------------------------------------------

!----IF YOU WANT TO USE WINDS----
    use_winds = .true. 

!---IF YOU WANT TO RAMP-UP WINDS GRADUALLY OVER ONE INERTIAL PERIOD------------------------------
!---SETTING FALSE BELOW WILL TURN ON THE WINDS TO THEIR MAX. AMPLITUDE INSTANTANEOUSLY-----------
    latewinds = .true. 
!------------------------------

   f_BoB=2*OMEGA*dsin(phi0deg*PI/180.d0) 

!  timday= dble(step)*dtf*TL/86400.d0 

   Tf=2.d0*3.14/f_BoB
   step_1f=int(dble(Tf)/dble(dtf*TL))


   if (pickup_step > 0.5) then
     rst_step = pickup_step
   else
     rst_step = 0 
   endif

 if (use_winds) then 

    time_ramp=dble(step-rst_step)*dtf*TL/Tf

    stressmax_x = 0.1d0 
    stressmax_y = 0.0d0 

      if (rst_step == 0) then 
      
         if ((step-rst_step) < step_1f) then
           stressmax_x = time_ramp*stressmax_x
           stressmax_y = time_ramp*stressmax_y
         endif 
      
      else
      
         if (latewinds) then
      
           if  ((step-rst_step) < step_1f) then
              stressmax_x = time_ramp*stressmax_x
              stressmax_y = time_ramp*stressmax_y
           endif
      
         endif  !(latewinds) 
      
      endif  !(rst_step ==0)  

 else

    stressmax_x=0.d0
    stressmax_y=0.d0 

 endif  ! (use_winds) 


  write(6,*) 'stressmax_x ',stressmax_x
  write(6,*) 'stressmax_y ',stressmax_y

    ycenter = 0.5*(yc(NJ)+yc(1)) 
    ymax = yc(NJ-10)-ycenter 
    do j=1,NJ
       diff = (yc(j)-ycenter) 
       if (dabs(diff).gt.ymax) then 
          stressx(j)= 0.d0 
          stressy(j)= 0.d0 
       else 
          stressx(j)= stressmax_x*dcos((diff/ymax)*0.5*PI) 
          stressy(j)= stressmax_y*dcos((diff/ymax)*0.5*PI) 
       end if 
    end do 


!--This part deals with the water type coefficients
  J_lambda1 = 0.6d0     
  J_lambda2 = 20.d0     
  J_A = 0.62d0      
!-----------------------------------

!--This part deals with short wave radiation (swr), qloss = Qnet - shortwave, qlatent = Latent heat loss, rainrate = precipitation (mm/hr) 
  swr = 0.d0
  qloss = 0.d0
  qlatent=0.d0
  rainrate=0.d0 
!----------------------------------

                   
return

END
