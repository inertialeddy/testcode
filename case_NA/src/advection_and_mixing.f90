subroutine advection_and_mixing(m,n,dtimel,step) 
#include "cppdefs.h"
!     ---------------------------------------------                     
USE header

INTEGER :: n,m,step, i,j,k

REAL(kind=rc_kind) :: dtimel

REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: var, var0             ! var=(s,T,u,v,w,Tr(it,:,:,:)                    
REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: uvarx, uvarx_advec    ! uvarx  is the source term, divergence of the advective fluxes 
REAL(kind=rc_kind), dimension(    1:NI  ,1:NJ  , 1:NK  ) :: vardif                ! vardif is the source term from diabatic processes 
REAL(kind=rc_kind), dimension(    1:NI  ,1:NJ  , 1:NK  ) :: vardif2               ! vardif is the source term from diabatic processes 

!INTEGER :: selvar
INTEGER :: iv_compute_kz, kk
INTEGER :: av_comp(5+ntr)


!Hcontent = 0.d0
!Hmixing = 0.d0

if (gotmlogic == 1) then
iv_compute_kz=2; 
end if

! By default, all variables will be subject to the following loop
av_comp=1

! If the simulation is a "rhoonly" simulation (rho is stored in s), T=0 and does not need to be updated.
#ifdef rhoonly
av_comp(1)=0
#endif

do selvar=1,5+ntr  

!--Tridiadonal solver is initialized to 0 before the loop
if (solverlogic==2) then
  mat_A = 0.d0; mat_B = 0.d0; mat_C = 0.d0; mat_D = 0.d0; mat_test = 0.d0;  
end if

!--Flux terms initialized to 0
uvarx = 0.d0; uvarx_advec = 0.d0; vardif = 0.d0  

!--var0() takes in the value of the variables at the start of the selvar loop.
  if(av_comp(selvar)==1) then
    if(selvar==1) then
      var=T(:,:,:,m); var0(:,:,:) = T(:,:,:,0);
    end if
    if(selvar==2) then
      var=s(:,:,:,m); var0(:,:,:) = s(:,:,:,0);
    end if
    if(selvar==3) then
      var=u(:,:,:,m); var0(:,:,:) = u(:,:,:,0);
    end if
    if(selvar==4) then
     var=v(:,:,:,m); var0(:,:,:) = v(:,:,:,0);
    end if
    if(selvar==5) then
     var=w(:,:,:,m); var0(:,:,:) = w(:,:,:,0);
    end if
    do it=1,ntr
      if(selvar==5+it) var=Tr(it,:,:,:,m)
    enddo

! computation of the advective fluxes, using QUICK scheme
    CALL advect(var,uvarx)

!if (selvar==2) then    
!    d_part1(1:NK) = uvarx(1,1,1:NK)
!end if

! computation of the horizontal diabatic fluxes
    if(selvar<5.5) then
      vardif=0.d0; vardif2=0.
      call mixing_horizontal(var,vardif);
     !call mixing_isopycnal(var,vardif,10.);PRINT*,"VISCOUS REDI";
     !call mixing_isopycnal_biharmonic(var,vardif,1.);PRINT*,"BIHARMONIC REDI"
      uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-vardif(1:NI,1:NJ,1:NK)
!if (selvar==2) then
!    d_part2(1:NK) = 0.d0 - vardif(1,1,1:NK)
!end if
    endif

    uvarx_advec(1:NI,1:NJ,1:NK) = uvarx(1:NI,1:NJ,1:NK)       

! computation of the vertical diabatic fluxes         
    vardif=0.d0;
    if(selvar<4.5) then        
      vardif=0.;      
      call mixing_vertical(var,vardif,m,step,iv_compute_kz);    

      if (solverlogic==1) then     
        uvarx(1:NI,1:NJ,1:NK)=uvarx_advec(1:NI,1:NJ,1:NK)-vardif(1:NI,1:NJ,1:NK)       
      end if        

      iv_compute_kz=0;
    endif

if (solverlogic==1) then
 if(selvar==3 .OR. selvar==4) then
   uvarx(1:NI,1:NJ,1)=uvarx(1:NI,1:NJ,1) + RR*(1.d0/(UL*delta) ) * ( Jac(1:NI,1:NJ,1)*wz(1:NI,1:NJ,1) ) * var(1:NI,1:NJ,1)
 end if
end if

if (solverlogic==2) then
!! summation for horizontal mixing part
if (selvar<(4.5)) then
do i=1,NI
  do j=1,NJ
    mat_A(i,j,1:NK,selvar)   = dtimel*Jacinv(i,j,1:NK)*mat_A(i,j,1:NK,selvar)      
    mat_B(i,j,1:NK,selvar)   = ( (dtimel*Jacinv(i,j,1:NK)*mat_B(i,j,1:NK,selvar)) - 1)     
    mat_C(i,j,1:NK-1,selvar) = dtimel*Jacinv(i,j,1:NK-1)*mat_C(i,j,1:NK-1,selvar)    

!!Previous statement->mat_D(i,j,1:NK,selvar)   = dtimel* (Jacinv(i,j,1:NK)*( mat_D(i,j,1:NK,selvar) - &
!!& uvarx_advec(i,j,1:NK)  ) ) + var(i,j,1:NK )      
    mat_D(i,j,1:NK,selvar) = 0.d0 - var0(i,j,1:NK ) - (dtimel*Jacinv(i,j,1:NK)*mat_D(i,j,1:NK,selvar)) + (dtimel*Jacinv(i,j,1:NK)*uvarx_advec(i,j,1:NK))   

!if (selvar==2 .AND. i==1 .AND. j==1) then    
!    d_part3(1:NK) = 0.d0 - var0(i,j,1:NK )  
!    d_part4(1:NK) = (dtimel*Jacinv(i,j,1:NK)*uvarx_advec(i,j,1:NK))
!end if

  call solve_tridiag(mat_A(i,j,1:NK,selvar), mat_B(i,j,1:NK,selvar), mat_C(i,j,1:NK,selvar), mat_D(i,j,1:NK,selvar), mat_test(i,j,1:NK,selvar), NK )    

  end do
end do
end if !! selvar<4.5
end if !! solverlogic==2

!if (selvar==2 .AND. solverlogic==1) then
!    d_part3(1:NK) = 0.d0 - var0(1,1,1:NK )
!    d_part4(1:NK) = (dtimel*Jacinv(1,1,1:NK)*uvarx_advec(1,1,1:NK))
!end if

!uvarx_x(1:NK) = uvarx_x(1:NK)*dtimel*Jacinv(32,173,1:NK)  
!uvarx_y(1:NK) = uvarx_y(1:NK)*dtimel*Jacinv(32,173,1:NK)  
!uvarx_z(1:NK) = uvarx_z(1:NK)*dtimel*Jacinv(32,173,1:NK)  

!if (selvar==2) then
!write(6,*) 'parts 1,2,3:',selvar
!  do kk=NK-10,NK
!    write(6,*) d_part1(kk), d_part2(kk), d_part4(kk)
!  end do
!write(6,*) 'Advection parts'
!  do k=NK-10,NK
!    write(6,*) uvarx_x(k),uvarx_y(k),uvarx_z(k)
!  end do
!write(6,*) 'dtimel.Jacinv.uvarx ->', selvar
!  do kk=1,NK
!    write(6,*) d_part4(kk)
!  end do
!end if

if ( selvar.lt.(4.5) ) then
if (solverlogic==2) then
    if(selvar==1) T(1:NI,1:NJ,1:NK,n)  = mat_test(1:NI,1:NJ,1:NK,selvar)
    if(selvar==2) s(1:NI,1:NJ,1:NK,n)  = mat_test(1:NI,1:NJ,1:NK,selvar)
    if(selvar==3) cx(1:NI,1:NJ,1:NK)   = mat_test(1:NI,1:NJ,1:NK,selvar)
    if(selvar==4) cy(1:NI,1:NJ,1:NK)   = mat_test(1:NI,1:NJ,1:NK,selvar)
else if (solverlogic==1) then
    if(selvar==1) T(1:NI,1:NJ,1:NK,n) = T(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
    if(selvar==2) s(1:NI,1:NJ,1:NK,n) = s(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
    if(selvar==3) cx(1:NI,1:NJ,1:NK)  = u(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
    if(selvar==4) cy(1:NI,1:NJ,1:NK)  = v(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
end if
end if

if(selvar==5) cz(1:NI,1:NJ,1:NK)  = w(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
do it=1,ntr
    if(selvar==5+it) Tr(it,1:NI,1:NJ,1:NK,n) =Tr(it,1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
enddo
endif
enddo !selvar number

!do k=NK,3,-1
!  Hcontent =  Hcontent + (    3.9d03*rho(NI/2,NJ/2,k)*T(NI/2,NJ/2,k,n)*  (   zf(NI/2,NJ/2,k+1) - zf(NI/2,NJ/2,k) )*DL   ) ;
!  Hmixing =  Hmixing - (3.93d03)*rho(NI/2,NJ/2,k)*KzTr(NI/2,NJ/2,k)*  ( T(NI/2,NJ/2,k+1,n)-T(NI/2,NJ/2,k,n))/(DL*(zc(NI/2,NJ/2,k+1)-zc(NI/2,NJ/2,k)) ) ;
!end do

!if (ivb==3) then
!write(6,*) 'TKE,EPS,L'    
!do k=NK-10,NK      
!   write(6,*) psom_tke(NI/2,NJ/2,k), psom_eps(NI/2,NJ/2,k), psom_l(NI/2,NJ/2,k)       
!end do       
!write(6,*) 'KzMom,KzTr'       
!do k=NK-10,NK      
!  write(6,*) u(NI/2,NJ/2,k,0), v(NI/2,NJ/2,k,0),shq(NI/2,NJ/2,k)      
!  write(6,*) psom_P(NI/2,NJ/2,k), psom_B(NI/2,NJ/2,k)        
!   write(6,*) KzMom(NI/2,NJ/2,k), KzTr(NI/2,NJ/2,k)           
!  write(6,*) KzMom(NI/2,NJ/2,k), psom_tke(NI/2,NJ/2,k), (zf(NI/2,NJ/2,k)*DL) ;         
!end do         
!end if       

#ifdef allow_particle
  ! PRINT*,"ALLOW PARTICLE"
#endif

return 
                                                                       
END

                                           
