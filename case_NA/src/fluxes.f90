
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

  INTEGER :: i,j,k,m,step,ycenter,ymax,rst_step,step_1f,rst_tau,stress_len
  REAL(kind=rc_kind)  ::  stressmax_x, stressmax_y,timday,time_ramp,diff,Tf,f_BoB 
! REAL(kind=rc_kind)  ::  J_lambda1, J_lambda2, J_A
  REAL(kind=rc_kind)  ::  alpha,taux_store   
  logical             ::  use_winds,latewinds
!---------------------------------------------------------

    use_winds = .true. 
    latewinds = .false. 

     stress_len = 42752 
!     rst_tau = 2*(2000+546) 
      rst_tau = 1500-250+2000
!      rst_tau=251

!---ANGLE FOR ROTATING WINDS------
! alpha is the angle between old and new x-axes

     alpha = -36.12*3.14d0/180.d0  ! transect C 
!      alpha = 0.d0 

!------------------------------- 

!  stressx= 0.0d0 
!  stressy= 0.0d0 


!f_BoB = 3.3212d-5  ! Coriolis parameter at 13.2N
f_BoB = 4.d-5  ! Coriolis parameter at 16N

!  timday= dble(step)*dtf*TL/86400.d0 

   Tf=2.d0*3.14/f_BoB
   step_1f=int(dble(Tf)/dble(dtf*TL))


   if (pickup_step > 0.5) then
     rst_step = pickup_step
   else
     rst_step = 0 
   endif

 if (use_winds) then 

   if (step == rst_step) then
    stressx_timeseries=0.d0
    stressy_timeseries=0.d0
   endif

    time_ramp=dble(step-rst_step)*dtf*TL/Tf

! This part deals with the wind stress
!-----------------------------------

!  do j=1,NJ
!    stressx(j) =  ( ( 0.13) * sin( j * 3.14159d0 / (NJ)  ) )  - 0.03d0 ; 
!    stressx(j) = max( stressx(j), 0.d0 );
!  end do

  stressmax_y = 0.d0 

  if (step == rst_step) then 
  
  write(6,*) 'READING IN WIND STRESS' 
  open(unit=40,file='/home/sr71d/PSOM/github/psom-gotm/psom_v05/code/asiri_forced_C/inc/taux_20min_C_leg1_dt108.txt')
  open(unit=41,file='/home/sr71d/PSOM/github/psom-gotm/psom_v05/code/asiri_forced_C/inc/tauy_20min_C_leg1_dt108.txt')
  do i=1,stress_len   
     read(40,*) stressx_timeseries(i) 
     read(41,*) stressy_timeseries(i) 
  end do 
  close(40)
  close(41) 

    write(6,*) stressx_timeseries(1:2), stressy_timeseries(1:2) 

    do i=1,stress_len
    taux_store=stressx_timeseries(i) 
    stressx_timeseries(i) =  stressx_timeseries(i)*dcos(alpha) + stressy_timeseries(i)*dsin(alpha) 
    stressy_timeseries(i) =            -taux_store*dsin(alpha) + stressy_timeseries(i)*dcos(alpha) 
!     stressy_timeseries(i) = 0.d0 
    end do

!    write(6,*) 'new stress ', stressx_timeseries(1:2), stressy_timeseries(1:2) 
  
  endif  ! (if step == rst_step)  

if (rst_step == 0) then 

   if ((step-rst_step) < step_1f) then
     stressmax_x = time_ramp*stressx_timeseries(1)
     stressmax_y = time_ramp*stressy_timeseries(1)
   else
  
!     stressmax_x = stressx_timeseries( (step-rst_step-step_1f)+1) 
!     stressmax_y = stressy_timeseries( (step-rst_step-step_1f)+1) 

     stressmax_x = stressx_timeseries( 2*(step-rst_step-step_1f)+1) 
     stressmax_y = stressy_timeseries( 2*(step-rst_step-step_1f)+1) 
  
   endif 

else

   if (latewinds) then

  
     if  ((step-rst_step) < step_1f) then
        stressmax_x = time_ramp*stressx_timeseries(1)
        stressmax_y = time_ramp*stressy_timeseries(1)
     else
        stressmax_x = stressx_timeseries( 2*(step-rst_step-step_1f)+1) 
        stressmax_y = stressy_timeseries( 2*(step-rst_step-step_1f)+1) 
     endif

   else
     stressmax_x = stressx_timeseries(rst_tau + 2*(step-rst_step)) 
     stressmax_y = stressy_timeseries(rst_tau + 2*(step-rst_step)) 

!     stressmax_x = stressx_timeseries(step-rst_step+rst_tau)
!     stressmax_y = stressy_timeseries(step-rst_step+rst_tau)

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

!--This part deals with short wave radiation
  swr = 0.d0
  qloss = 0.d0
!----------------------------------

                   
return

END
