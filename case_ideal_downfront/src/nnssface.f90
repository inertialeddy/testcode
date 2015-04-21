!! This subroutine calculates shear and buoyancy frequencies at cell faces, by a linear
!! interpolation from the cell centers

subroutine nnssface(ii,jj)

use header
implicit none

!!! Local Variables
integer :: ii,jj,k

double precision :: temp

!!!! Calculate shear separately
do k=0, (NK-1)
   ush(ii,jj,k) = (u(ii,jj,k+1,1) - u(ii,jj,k,1))*UL / ((zc(ii,jj,k+1) - zc(ii,jj,k))*DL)
   vsh(ii,jj,k) = (v(ii,jj,k+1,1) - v(ii,jj,k,1))*UL / ((zc(ii,jj,k+1) - zc(ii,jj,k))*DL)
   shq(ii,jj,k) = ush(ii,jj,k)*ush(ii,jj,k) + vsh(ii,jj,k)*vsh(ii,jj,k)
end do

  ush(ii,jj,NK) = ush(ii,jj,NK-1)
  vsh(ii,jj,NK) = vsh(ii,jj,NK-1)
  shq(ii,jj,NK) = ush(ii,jj,NK)*ush(ii,jj,NK) + vsh(ii,jj,NK)*vsh(ii,jj,NK)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


CALL diag_n2()

!! Calculate buoyancy frequency by linear interpolation from cell centers
do k=1, (NK-1)
  
   temp = (freqN2(ii,jj,k)*(zc(ii,jj,k+1)-zf(ii,jj,k))) + (freqN2(ii,jj,k+1)*(zf(ii,jj,k)-zc(ii,jj,k)))
   nnface(ii,jj,k) = temp/(zc(ii,jj,k+1) - zc(ii,jj,k))
end do
    
   temp = (freqN2(ii,jj,NK)*(zf(ii,jj,NK)-zf(ii,jj,NK-1))) - (nnface(ii,jj,NK-1)*(zf(ii,jj,NK)-zc(ii,jj,NK)))
   nnface(ii,jj,NK) = temp/(zc(ii,jj,NK)-zf(ii,jj,NK-1))

   temp = (freqN2(ii,jj,1)*(zf(ii,jj,0)-zf(ii,jj,1))) - (nnface(ii,jj,1)*(zf(ii,jj,0)-zc(ii,jj,1)))
   nnface(ii,jj,0) = temp/(zc(ii,jj,1)-zf(ii,jj,1))

!! Make everything positive
do k=0,NK
  if (nnface(ii,jj,k).lt.0) then
    nnface(ii,jj,k) = nnface(ii,jj,k) * (-1)
  else
    continue
  endif
end do

!! Simply assign buoyancy frequency nnface to freqN2, thats it
!do k=0,NK
!    nnface(ii,jj,k) = freqN2(ii,jj,k)
!end do


!!! You are now calculating NN at cell faces. So no need to linearly interpolate


do k=0,NK
  Rig(ii,jj,k) = nnface(ii,jj,k)/shq(ii,jj,k)
end do

!if (ii.eq.48 .AND. jj.eq.97 .AND. cnt.eq.102) then
!write(6,*) 'Velocity and length scales for non-dimensionalization', UL, DL
!write(6,*) 'At 34th time step'
!write(6,*) 'u velocity at 2,2'
!do k=0,NK
!   write(6,*) u(ii,jj,k,1) * UL
!end do
!write(6,*) 'v velocity at 2,2'
!do k=0,NK
!   write(6,*) v(ii,jj,k,1) * UL
!end do
!write(6,*) 'cell center at 2,2'
!do k=0,NK+1
!   write(6,*) zc(ii,jj,k) * DL
!end do
!write(6,*) 'shear at 2,2'
!do k=0,NK
!   write(6,*) shq(ii,jj,k)
!end do
!write(6,*) 'All values printed'
!endif

!write(6,*) 'nnssface called correctly.. working'

return
end subroutine nnssface
