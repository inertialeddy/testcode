subroutine mixing_vertical(var,vardif,m,step,iv_compute_kzl) 
  !     ---------------------------------------------                     
  USE header
  !     Wind Stress can be specified here                                 
  !     Kz*du/dz = -Tx/rho,   Kz*dv/dz = -Ty/rho                          
  !     use level m                                                       
  !     computes d2s/dz2 at the cell centers.                             
  !     fricu and fricv are u,v diffusion terms (in header)            
  !     that need to be saved (if m.eq.0) for n2budget)                   

  implicit none 
  REAL(kind=rc_kind) :: dSin,posneg,facreverse, swrtemp
  integer i,j,k,m,step,nday, jj  
  integer iv_compute_kzl 
  REAL(kind=rc_kind) :: zdep,yw,yset,ymid,stress,ustar,Ekdepth,ff  , swrmax                                    
  REAL(kind=rc_kind), dimension(0:NI+1,0:NJ+1, 0:NK+1) :: var 
  REAL(kind=rc_kind) :: vardif(NI,NJ,NK)
  REAL(kind=rc_kind) :: dvardzfc(0:NK)
  REAL(kind=rc_kind) :: Kdudzb,Kdvdzb,Kdudzt,Kdvdzt,Kdfluxdzt,Kdfluxdzb, Kdfluxdzz(1:NK), wzkth,fact,Cp         
  REAL(kind=rc_kind) :: fac,facb,dfactor,Kzmax,facy,rhoinv,diff                  
  REAL(kind=rc_kind) :: wgt,day,KzmaxTr,ztransit,zextent,thy , Tsghost=30.0125d0, Tbghost=24.7443      
  REAL(kind=rc_kind) :: Lv ! Latent heat of vaporization 
  REAL(kind=rc_kind) :: evap_flux, precip_flux 
  PARAMETER (Kzmax= 1.d-3, KzmaxTr=1.d-3) !     The viscosity is Kz*Kzmax                                         

  INTEGER OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,CHUNK,NPROC
!-------------------------------------------------------------

  Lv=2.5d6

                                                    
  facb = RR*DL      ! Linear Drag                                                       
  ! facb = RR*DL*UL ! Quadratic drag                                                    

  fac= 1.d0/(UL*DL*delta) 
  fact = DL/UL 


!if (step>=72) then
!  swrmax = 220.00;
!else
!  swrmax = 220.00*step/72;
!end if

!do jj=1,NJ
!    if (jj <= NJ/2) then
!        swr(jj) = swrmax*(NJ/2-jj)/(NJ/2) ;
!        qloss(jj) = 0.d0;
!    else
!        swr(jj) = 0;
!        qloss(jj) = swrmax*(NJ/2-jj)/(NJ/2) ;
!    end if
!end do

!-------------------
! COMPUTATION OF KZ
!-------------------
 if(iv_compute_kzl==1) then

! write(6,*) 'Ekman diffusivity'
 do j=1,NJ 
   do i=1,NI 
          stress= sqrt(stressx(j)*stressx(j)+ stressy(j)*stressy(j) ) 
          ustar = sqrt( stress /R0 ) 
          ff= ffc(i,j)*FPAR 
          Ekdepth= 0.4d0*ustar/ff 
          zextent= 0.5d0*Ekdepth !     Set zextent
          ztransit = -Ekdepth 
          do k=NK,0,-1    
            Kz(i,j,k)= 1.d0   
            zdep = zf(i,j,k)*DL    
            thy = (1.d0 +tanh(((zdep-ztransit)/zextent)*PI))*0.5d0                                              
            Kz(i,j,k)= max(0.01d0,thy) * KzmaxTr                      
            Kz(i,j,k)= max(Kz(i,j,k), 1.e-05)

!            if (Kz(i,j,k) .le. 1.e-05) then
!              Kz(i,j,k) = 1.e-05
!            end if

          end do             
      enddo                 
    enddo                                                                 

 else if(iv_compute_kzl==2) then      
    do j=1,NJ         
      do i=1,NI           

       if (ivb==1) then        
         call shearn2(i,j,step)        
       end if        
         call couple_gotm(i,j)       

      enddo       
    enddo       
 endif     


!$OMP PARALLEL DO PRIVATE(j,i,k,dvardzfc,dfactor) SHARED(vardif) COLLAPSE(2)
  do j=1,NJ 
    do i=1,NI 

!!--- This part is for GOTM only. --------------------
     if (gotmlogic==1) then 
         if (selvar==1 .OR. selvar==2 ) then 
           Kz(i,j,1:NK) = KzTr(i,j,1:NK) 
         end if 
         if (selvar==3 .OR. selvar==4 .OR. selvar==5 ) then 
           Kz(i,j,1:NK) = KzMom(i,j,1:NK) 
         end if 
     end if 
!! ---------------------------------------------------

      !-------------------
      ! ADDITIONAL COMPUTATIONS
      !-------------------
                                                                      
      !     Quadratic drag                                                    
      !=            Kdudzb= facb*u(i,j,1,m)*abs(u(i,j,1,m))                  
      !=            Kdvdzb= facb*v(i,j,1,m)*abs(v(i,j,1,m))                  
      !     Linear Drag                                                       
      !Kdudzb= facb*u(i,j,1,m) 
      !Kdvdzb= facb*v(i,j,1,m) 
      !                                                              
!      rhoinv = 1.d0/rho(i,j,NK) 
      rhoinv = 1.d0/R0
!      Kdudzt= stressx(j)*rhoinv*fact 
!      Kdvdzt= stressy(j)*rhoinv*fact 


!! Add sinusoidal heat flux
!swrtemp = 0d0*sin( (2*3.14159d0/(24*3600))*step*dtf*(1.d05)    )
!if (swrtemp.lt.(0.d0)) then
!swrtemp = 0.d0
!end if

if (selvar==1) then ! Temperature
  Kdfluxdzt = DL*( swr(j) + qloss(j) )/(R0*4187.d0)
  do k=1,NK   
    swrd(k) = swr(j)*( J_A*exp(zf(NI/2,NJ/2,k)*DL/J_lambda1) + (1 - J_A)*exp(zf(NI/2,NJ/2,k)*DL/J_lambda2)   )   
    Kdfluxdzz(k) = DL*( swrd(k) )/(R0*4187.d0)    
  end do   
  Kdfluxdzb = 0.d0   
else if (selvar==2) then ! Salinity
  evap_flux  = var(i,j,NK)*dabs(qlatent(j))/(Lv*R0) 
! the rainrate is in mm/hr, converted below to m/s
  precip_flux= var(i,j,NK)*rainrate(j)*1d-3/3600.d0 

!  write(6,*) 'evap_flux ',evap_flux
!  write(6,*) 'precip_flux ',precip_flux

!  Kdfluxdzt = 0.d0
  Kdfluxdzt = DL*(evap_flux-precip_flux)
  Kdfluxdzb = 0.d0
  Kdfluxdzz = 0.d0
else if (selvar==3) then ! u velocity
  Kdfluxdzt = stressx(j)*rhoinv*fact 
  Kdfluxdzb = 0.d0
  Kdfluxdzz = 0.d0
else if (selvar==4) then ! v velocity	
  Kdfluxdzt = stressy(j)*rhoinv*fact 
  Kdfluxdzb = 0.d0
  Kdfluxdzz = 0.d0
else if (selvar==5) then ! w velocity
  Kdfluxdzt = 0.d0
  Kdfluxdzb = 0.d0
  Kdfluxdzz = 0.d0
end if

      !------------------------------------------------------
      ! COMPUTATION OF THE VARIATIONS OF THE VARIABLE ALONG Z
      !------------------------------------------------------
      do k=1,NK-1 
        dvardzfc(k)= wz(i,j,k)*(var(i,j,k+1)-var(i,j,k)) 
      end do 
      dvardzfc(0)= 0.d0 
      dvardzfc(NK)= 0.0 
      
      !-----------------------------------
      ! COMPUTATION OF THE FLUX DIVERGENCE
      !-----------------------------------
      k=1 
!     ---
      dfactor=  fac*Jac(i,j,k)*wz(i,j,k)           
      if (solverlogic==1) then
        vardif(i,j,k)= dfactor*(dvardzfc(k)*Kz(i,j,k)- dvardzfc(k-1)*Kz(i,j,k-1) )    ! Neumann   
!       vardif(i,j,k)= dfactor*(dvardzfc(k)*Kz(i,j,k)-  wz(i,j,k-1)*(var(i,j,k) - varbghost )*Kz(i,j,k) )     ! Dirichlet    
      else if (solverlogic==2) then   
!       mat_B(i,j,k,selvar) = (-1.d0)*dfactor*(Kz(i,j,k)*wz(i,j,k) + Kz(i,j,k)*wz(i,j,k-1) )  !! Dirichlet  
!       mat_D(i,j,k,selvar) = dfactor*Kz(i,j,k)*wz(i,j,k-1)*Tbghost                           !! Dirichlet  
        mat_B(i,j,k,selvar) = (-1.d0)*dfactor*(Kz(i,j,k)*wz(i,j,k) )            !! Neumann           
        mat_D(i,j,k,selvar) = (-1.d0)*dfactor*Kdfluxdzb                         !! Neumann           
        mat_C(i,j,k,selvar) = dfactor*Kz(i,j,k)*wz(i,j,k)          
      end if

      do k=2,NK-1 
!     -----------
      dfactor=  fac*Jac(i,j,k)*wz(i,j,k) 
      if (solverlogic==1) then
        vardif(i,j,k)= dfactor*(dvardzfc(k)*Kz(i,j,k)- dvardzfc(k-1)*Kz(i,j,k-1) + Kdfluxdzz(k) - Kdfluxdzz(k-1)  )        
      else if (solverlogic==2) then        
        mat_A(i,j,k,selvar) = dfactor*kz(i,j,k-1)*wz(i,j,k-1)           
        mat_B(i,j,k,selvar) = (-1.d0)*dfactor*( Kz(i,j,k)*wz(i,j,k) + Kz(i,j,k-1)*wz(i,j,k-1) )           
        mat_C(i,j,k,selvar) = dfactor*Kz(i,j,k)*wz(i,j,k)         
        mat_D(i,j,k,selvar) = dfactor*( Kdfluxdzz(k) - Kdfluxdzz(k-1) )            
      end if        
      end do 
                                                                                                                    
      k=NK 
!     ----                                                        
      dfactor=  fac*Jac(i,j,k)*wz(i,j,k)    
      if (solverlogic==1) then     
!       vardif(i,j,k)= dfactor*KzmaxTr*(dvardzfc(k)*Kz(i,j,k)- dvardzfc(k-1)*Kz(i,j,k-1) )     ! Previous one
        vardif(i,j,k)= dfactor*( Kdfluxdzt - dvardzfc(k-1)*Kz(i,j,k-1)  )        
      else if (solverlogic==2) then     
        mat_A(i,j,k,selvar) = dfactor*Kz(i,j,k-1)*wz(i,j,k-1)      
!       mat_B(i,j,k,selvar) = (-1.d0)*dfactor*(Kz(i,j,k-1)*wz(i,j,k-1) + Kz(i,j,k)*wz(i,j,k)  )   !! Dirichlet    
!       mat_D(i,j,k,selvar) = dfactor*wz(i,j,k)*Kz(i,j,k)*Tsghost                                 !! Dirichlet
        mat_B(i,j,k,selvar) = (-1.d0)*dfactor*Kz(i,j,k-1)*wz(i,j,k-1)                 !! Neumann    
        mat_D(i,j,k,selvar) = dfactor*Kdfluxdzt                                       !! Neumann    
      end if   

    end do ! i   
  end do ! j   
!$OMP END PARALLEL DO                                                                      
return 
END                                           
                                                                                                                                                

      subroutine viscosity(dudz,dvdz,drdz,i,j) 
!     ---------------------------------------------                     
        USE header
!     Compute a Richardson number based vertical viscosity              
!     coefficient Kz using the final velocities u,v and density         
!     with the n=0 subscript                                            
!     Kz is situated at k cell faces                                    
!     The algorithm is from Rev Geophys., p 373, 1994. Large et al.     
!                                                                       
!     Assumes that rho is evaluated just prior to this call             
                                                                        
      integer i,j,k 
                                                                        
      integer n1 
      parameter (n1=3) 
      REAL(kind=rc_kind) :: fac,bvfreq,grho,RiCr,Ri,vshear 
      REAL(kind=rc_kind) :: dudz(0:NK),dvdz(0:NK),drdz(0:NK) 
      parameter (RiCr= 0.7d0) 
                                                                        
      grho= 9.81/R0 
!     fac= DL/(UL*UL)                                                   
      DLinv = 1.0/DL 
      fac = UL*UL/(DL*DL) 
!                                                                       
!     Set kz(k=0)= 1. It is required only for the s,T equations, since  
!     for the momentum equations,  K*du/dz= ru.                         
!     Kz is not needed at k=0 and NK since stress bcs are used.         
      Kz(i,j,0)= 0.d0 
      Kz(i,j,NK)= 0.d0 
!      if  (i.eq.16) write(100,*) 'j = ', j                             
      do k=1,NK-1 
         bvfreq=  -grho*drdz(k)*DLinv 
!     BVfreq is N**2 and is in s^-2 if re-dim by DLinv                  
         vshear= (dudz(k)*dudz(k) + dvdz(k)*dvdz(k))*fac 
         if (vshear.eq.0.d0) then 
!         if (vshear.le.1.d-12) then                                    
            Ri= 100.d0 
         else 
            Ri= bvfreq/vshear 
         end if 
!     unstable density profile => Ri < 0                                
         if (Ri.le.0.d0) then 
            Kz(i,j,k)= 1.d0 
         else if (Ri.lt.RiCr) then 
            Kz(i,j,k)= (1.d0 - (Ri*Ri)/(RiCr*RiCr))**n1 
         else 
            Kz(i,j,k)= 0.d0 
         endif 
!         if ((i.eq.NI/2).and.(j.eq.NJ/2)) write(6,*) Ri ,bvfreq,vshear,
!     &        Kz(i,j,k)                                                
      end do 
                                                                        
!     The value of Kz to be used is  Kz(k) *Kzmax                       
                                                                        
      return 
      END                                           





