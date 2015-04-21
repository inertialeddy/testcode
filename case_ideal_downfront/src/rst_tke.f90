subroutine restart_tke()
! This initializes turbulence routines in GOTM

!!! Access relevant subroutines from init_tke
use header
use turbulence,           only : eps1d => eps
use turbulence,           only : tke, L, num, nuh, nus
use turbulence,           only : init_turbulence, do_turbulence, clean_turbulence
use mtridiagonal,         only: init_tridiagonal

implicit none

!!! Some local variables
integer                 :: i, j, k, kdp
!!!!!!!!!!! Initialize Turbulence Modules !!!!!!!!!!!!!!
call init_turbulence(10, 'gotmturb.nml', NK)
call init_tridiagonal(NK)
write(*,*) 'Turbulence Modules initialized'

!===== Initialize lengthscale profile =============
do i=1,NI
   do j=1,NJ

     do k=0,NK
       psom_l(i,j,k)= (2*psom_tke(i,j,k))**1.5/(B1_MY*psom_eps(i,j,k))
     end do 

!===========
  if (i.eq.NI/2 .AND. j.eq.NJ/2) then
  write(6,*) 'Length-scale values at near-surface'
  do k=NK-10,NK
    write(6,*) psom_l(i,j,k)
  end do
  end if

   end do
end do
!==================================================

write(6,*) NI, NJ, NK

write(6,*) 'Wind stress at mid point is', stx(NI/2,NJ/2)

write(6,*) 'Friction velocity at mid point is', u_fric(NI/2,NJ/2)

write(6,*) 'Turbulent Kinetic Energy (NK)', psom_tke(NI/2, NJ/2, NK)

write(6,*) 'Turbulent Kinetic Energy (NK-2)', psom_tke(NI/2, NJ/2, NK-2)

write(6,*) 'Macro Lengthscale', psom_l(NI/2, NJ/2, NK-1)

write(6,*) 'Macro Lengthscale below the surface', psom_l(NI/2, NJ/2, NK-4)

write(6,*) 'Dissipation Rate of TKEat mid point', psom_eps(NI/2, NJ/2, NK)

write(6,*) 'Dissipation Rate of TKE below the surface', psom_eps(NI/2, NJ/2, NK-4)


return
end subroutine restart_tke


