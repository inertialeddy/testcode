!!! This subroutine interfaces with gotm and calculates diffusion coefficients
subroutine couple_gotm(I,J)
! Downloads copy
!!! Subroutines accessed from within couple_gotm
USE header
use turbulence,           only : eps1d => eps
use turbulence,           only : tke, L, num, nuh, nus, sig_e,sig_k
use turbulence,           only : P1d=> P, B1d=>B 
use turbulence,           only : ugradtke, vgradtke, wgradtke 
use turbulence,           only : ugradeps, vgradeps, wgradeps 
use turbulence,           only : ce1,ce2,ce3 
use turbulence,           only : init_turbulence, do_turbulence, clean_turbulence
use mtridiagonal,         only : init_tridiagonal

use turbulence,      only : wgradtke, wgradeps
use turbulence,      only : ugradtke, vgradtke, ugradeps, vgradeps

! dtf (PSOM) => dt (GOTM) ->> The time step
! NK (PSOM) => nlev (GOTM) ->> Number of levels along vertical
!!! Velocity variables
!!! u, v, w declared at an array size (0:NI+1,0:NJ+1,0:NK+1,0:1)
!!! freqN2(i,j,k) is the buoyancy frequency (mymodules.f90)
!!! u_fric(I,J) is the surface friction velocity at I,J (mymodules.f90)
!!! Note - This variable is only used for setting initial tke, eps, L

implicit none
integer, intent(in)     ::  I,J


!!! Local variables
integer                 ::  iii, ic
!!! 
integer                 ::  k, MaxItz0b 
REAL(kind=rc_kind)      ::  charnockval, rrb, clip,epsovertke 

!!!! layer thickness (resolution), shear frequency, sea grass tke(irrelevant), buoyancy frequency
!REAL(kind=rc_kind), dimension(  0:NK  )  :: thk, SS2, TKEP, NN1d

REAL(kind=rc_kind), dimension( 1:NK ) :: depsdz,dkdz 
!!! Molecular Viscosity, heat and salt diffusivities (used for gotm calculations)
REAL(kind=rc_kind)      :: avmolu , avmolt , avmols  

!!! Variables for interfacing with gotm
REAL(kind=rc_kind)      :: dpth, gotm_dt

!!! Bottom friction velocity, surface and bottom roughness parameters
REAL(kind=rc_kind)      :: u_taub, z0s_gotm, z0b_gotm

!!! Bottom roughness = h0b
!!!  Note: z0b=0.03*h0b+0.1*(molecular viscosity)/ustar,
!!!  Von-Karman constant, bottom friction coefficient,
REAL(kind=rc_kind)      :: h0b   

!---------------------------------
charnockval = 1400.d0; clip = 1.0E0; MaxItz0b = 1;
avmolu = 1.3d-6; avmolt = 1.4d-7; avmols = 1.1d-9;
h0b = 0.05;
!---------------------------------
!!! Calculate friction velocities, bottom and surface roughness lengths, shear and
!!! buoyancy frequencies, set sea grass production to 0E0

! Update the gradtke gradeps terms here
! write(6,*) 'reaching gradtke'


!============ Compute grad terms of tke and eps ===============
do k=2,(NK-1)

!  psom_wgradtke(I,J,k) = (w(I,J,k+1,0)*psom_tke(I,J,k+1)-w(I,J,k-1,0)*psom_tke(I,J,k-1))/(zc(I,J,k+1)-zc(I,J,k-1)) * (UL/1000.d0)
!  psom_wgradeps(I,J,k) = (w(I,J,k+1,0)*psom_eps(I,J,k+1)-w(I,J,k-1,0)*psom_eps(I,J,k-1))/(zc(I,J,k+1)-zc(I,J,k-1)) * (UL/1000.d0)
   psom_wgradtke(i,j,k) = 0.d0 
   psom_wgradeps(i,j,k) = 0.d0 

  if (k.eq.2) then
    psom_wgradtke(I,J,k-1) = psom_wgradtke(I,J,k)
  end if
  if (k.eq.NK-1) then
    psom_wgradeps(I,J,k+1) = psom_wgradeps(I,J,k)
  end if

!  psom_ugradtke(I,J,k) = (u(I+1,J,k,0)*psom_tke(I+1,J,k)-u(I-1,J,k,0)*psom_tke(I-1,J,k))/(xc(I+1)-xc(I-1)) * (UL/1000.d0) 
!  psom_ugradeps(I,J,k) = (u(I+1,J,k,0)*psom_eps(I+1,J,k)-u(I-1,J,k,0)*psom_eps(I-1,J,k))/(xc(I+1)-xc(I-1)) * (UL/1000.d0) 
   psom_ugradtke(i,j,k) = 0.d0
   psom_ugradeps(i,j,k) = 0.d0

!  psom_vgradtke(I,J,k) = (v(I,J+1,k,0)*psom_tke(I,J+1,k)-u(I,J-1,k,0)*psom_tke(I,J-1,k))/(yc(J+1)-yc(J-1)) * (UL/1000.d0) 
!  psom_vgradeps(I,J,k) = (u(I,J+1,k,0)*psom_eps(I,J+1,k)-u(I,J-1,k,0)*psom_eps(I,J-1,k))/(xc(J+1)-xc(J-1)) * (UL/1000.d0) 
   psom_vgradtke(i,j,k) = 0.d0
   psom_vgradeps(i,j,k) = 0.d0

end do

!!! Time Step for PSOM/GOTM is set in momentum equation
!!! according to Runge-Kutta-3rd order scheme...

!!! Calculating layer thickness - added to cell centres
DO k = 1, NK
   thk(k) = (zf(I,J,k) - zf(I,J,k-1))*DL
   lthk(I,J,k) = thk(k)
END DO
   thk(0) = thk(1)

!!! Depth
dpth = (zf(I,J,NK) - zf(I,J,0))*DL
!! Set a time step for gotm

!!! Buoyancy frequency freqN2(I,J,K) and shear squared
do k = 0, NK
   NN1d(k) = nnface(I,J,k)
   SS2(k)  = shq(I,J,k)
!!! Change the values to cell face, use subroutine nnssface
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SeaGrass Turbulent Kinetic Energy Production
TKEP(:) = 0E0

!-------------------------------------------------------

!!! Calculate surface roughness length. Surface friction velocity is known already.
!!! Initialize the wind stress profile
u_fric(I,J) = (sqrt(stressx(J)*stressx(J) + stressy(J)*stressy(J))/1027.d0)**0.5

!IF (I.eq.(NI/2) .AND. J.eq.(NJ/2))
!  write(6,*) 'Validation- friction velocity'
!  write(6,*) u_fric(I,J), u_fric(I-1,J-1), u_fric(I+1,J+1)
!ENDIF

z0s_gotm = (u_fric(I,J)**2)*charnockval/9.81

!!!!!!!!!!!! Bottom Friction and Roughness !!!!!!!!!!!!!
u_taub = 6.7e-6  ! Limiting case

do iii=1,MaxItz0b
   z0b_gotm = 0.1*avmolu/max(avmolu,u_taub)+0.03*h0b
   
!compute the factor r (version 1, with log-law)   
   rrb=k_von/(log((z0b_gotm+thk(1)/2.d0)/z0b_gotm))

!  compute the factor r (version 2, with meanvalue log-law)
!  frac=(z0b+h(1))/z0b
!  rrb=kappa/((z0b+h(1))/h(1)*log(frac)-1.)

!  compute the friction velocity at the bottom
!  u_taub = rrb * sqrt(u(I,J,1,0)**2 + v(I,J,1,0)**2)

end do
u_taub = 0.d0
gotm_dt = tsp

!!! Keeping a parameter open for gotm - putting bottom roughness 0
 z0b_gotm = 0

!! ---- write(6,*) 'gradupdate'
!!! Update GOTM variables for tke, eps, l from PSOM
do k=0, NK
   tke(k)   = psom_tke(I,J,k) 
   eps1d(k) = psom_eps(I,J,k) 
   L(k) = psom_l(I,J,k) 
   wgradtke(k) = psom_wgradtke(I,J,k)
   wgradeps(k) = psom_wgradeps(I,J,k)
   ugradtke(k) = psom_ugradtke(I,J,k)
   ugradeps(k) = psom_ugradeps(I,J,k)
   vgradtke(k) = psom_vgradtke(I,J,k)
   vgradeps(k) = psom_vgradeps(I,J,k)
end do

u_taub = 0.d0
z0b_gotm = 0.d0

!!! Now, call do_turbulence subroutine, to update diffusivities
call do_turbulence(NK, gotm_dt, dpth, u_fric(I,J), u_taub, z0s_gotm, z0b_gotm, thk, NN1d, SS2, TKEP)
!-----------------------------------
!!!WRITE OUT THE TERMS IN THE EPSILON EQUATION: (EQ. 163) IN MANUAL 

!if (ivb==3) then 
!   do k=1,NK
!     depsdz(k)=(num(k)/sig_e)*(eps1d(k)-eps1d(k-1))/thk(k)
!   end do
!   k=1 
!    D_eps(i,j,k)=(depsdz(2)-depsdz(1))/thk(k)
!   do k=2,NK-1
!    D_eps(i,j,k)=0.5*(depsdz(k+1)-depsdz(k-1))/thk(k) 
!   end do
!   k=NK
!    D_eps(i,j,k)=(depsdz(NK)-depsdz(NK-1))/thk(k) 
!
!
!    do k=1,NK
!     epsovertke=eps1d(k)/tke(k) 
!     eps_shear(i,j,k)= ce1*epsovertke*P1d(k)
!     eps_buoy(i,j,k) = ce3*epsovertke*B1d(k)
!     eps_dest(i,j,k) = ce2*epsovertke*eps1d(k) 
!    end do 
!
!   do k=1,NK
!     dkdz(k)=(num(k)/sig_k)*(tke(k)-tke(k-1))/thk(k)
!   end do
!   k=1 
!    D_tke(i,j,k)=(dkdz(k+1)-dkdz(k))/thk(k)
!   do k=2,NK-1
!    D_tke(i,j,k)=0.5*(dkdz(k+1)-dkdz(k-1))/thk(k) 
!   end do
!   k=NK
!    D_tke(i,j,k)=(dkdz(k)-dkdz(k-1))/thk(k) 
!
!   do k=1,NK
!     sgs_shear(i,j,k)=P1d(k)
!     sgs_buoy(i,j,k) =B1d(k)
!     sgs_dest(i,j,k) =eps1d(k) 
!   end do 
!
!
!
!endif 


 


!-----------------------------------
!!! Eliminate very low values/ Confine then within the limit 1.0d-05
do k=0, NK
!! Update diffusivities in PSOM
if (num(k).ge.(1.0d-05) .AND. num(k).lt.clip) then
      KzMom(I,J,k) = num(k)
    else if (num(k).ge.clip) then
      KzMom(I,J,k) = clip
    else
      KzMom(I,J,k) = 1.0d-05
endif

if (nuh(k).ge.(1.0d-005) .AND. nuh(k).le.clip) then  
     KzTr(I,J,k) = nuh(k)  
   else if (nuh(k).ge.clip) then  
     KzTr(I,J,k) = clip  
   else  
     KzTr(I,J,k)  = 1.0d-05  
endif  

!! Change tke, eps, l only when ivb=3/3rd RK step
if (ivb==3 ) then
!! Update tke, eps, L in PSOM
   psom_tke(I,J,k) = tke(k) 
   psom_eps(I,J,k) = eps1d(k) 
   psom_l(I,J,k) = L(k) 
   psom_wgradtke(I,J,k) = wgradtke(k)
   psom_wgradeps(I,J,k) = wgradeps(k)
   psom_ugradtke(I,J,k) = ugradtke(k)
   psom_ugradeps(I,J,k) = ugradeps(k)
   psom_vgradtke(I,J,k) = vgradtke(k)
   psom_vgradeps(I,J,k) = vgradeps(k)
end if

end do

return
end subroutine couple_gotm

