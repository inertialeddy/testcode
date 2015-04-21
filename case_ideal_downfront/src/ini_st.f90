      subroutine ini_st 
!     --------------------
!     DESCRIPTION

!     Initializes density for model using an analytic 
!     density function, derived from an N^2 profile
!     Across-front use the TANH profile with tightness
!     as a measure of frontal spread
!     Larger the factor, tighter the front and larger b_xx.

!     ____________________
!     --------------------
!     SET VARIABLES

      USE header 
!      implicit logical (a-z)
      implicit none 
!      include 'header.f'
      integer  i,j,k,trigger,m,TP
      double precision z, analyticeval, sbkgrnd, z1, z2, slfac
      double precision thy, slfacnew, amplitude, widbar
      double precision nml, n21, n2dep, wiggles, mmm, ttt
      double precision  tt, qq, slope, yinter

!     yfront is the distance of the front from the coast
      yfront = 0.5d0*dble(yc(NJ+1) +yc(0))
!     ____________________
!     --------------------
!     DECLARE VARIABLES

      tightness = 0.03d0  

!      mldepth = 65.0d0     !Mixed layer depth
!      nml = 4.0d-6         !N^2 value in mixed layer
      mldepth = 100.d0     !Mixed layer depth
      nml = 1.0d-6         !N^2 value in mixed layer
      n21 = 6.0d-5         !N^2 value at baroclinic depth
      n2dep = 1.5d-5       !N^2 value below depth of baroclinicity
      widbar = 10.d0       !Value used to determine width

!-----------when front is confined to the mixed layer----------
       z1 = mldepth     ! Depth over which frontal strength is constant
       z2 = mldepth +30.d0 ! Max depth to which front extends
!------------when front extends below the mixed layer-----------

!      z1 = 2*mldepth       !Depth in which front is constant
!      z2 = mldepth + 3.d0*mldepth  !Max depth to which front extends
!      z2= 400.d0  
!-------------------------------------------------------------
      slfac = 0.1d0    
      trigger = 0  !0 means front trigged
!      wiggles = 10.d0
!      amplitude = 0.0d0
!     ____________________
!     --------------------
!     SET BACKGROUND DENSITY USING ANALYTIC PROFILE
      write(6,*) "tightness= ", tightness
      mmm = -zc(10,10,7)
      ttt = analyticeval(mldepth,mmm,nml,n21,n2dep,widbar)
      tt = mldepth + 13*widbar-1
!      qq = mldepth + 13*widbar


      do m=0, NK+1
      z=-zc(10,10,m)*1000
      TP = m-2
      if (z.le.tt) goto 100
      end do
100   continue

      write(6,*) "TP =  ", TP 


      do j=0,NJ+1
         do i=0,NI+1
            do k=TP,NK+1
                z= -zc(i,j,k)*1000
!                write(6,*) z
!                if (k.ge.TP) then
!                  write(6,*) z
!                 write(6,*) k
                  sbkgrnd = analyticeval(mldepth,z,nml,n21,n2dep,widbar)
!                  write(6,*) sbkgrnd
                  s(i,j,k,0) =  sbkgrnd
!                else 

!                  slope = (s(1,1,TP-1,0)-s(1,1,TP-3,0))  & 
!                       / ((zc(1,1,TP-1)-zc(1,1,TP-3))*(-1000))


!                 sbkgrnd = slope * z - slope * zc(i,j,TP)*(-1000)  & 
!                + s(i,j,TP,0)
!                 s(i,j,k,0) = sbkgrnd
!               end if
            end do
         end do
      end do

!      do j=1,NJ 
      slope = (s(10,10,TP+1,0)-s(10,10,TP+2,0))  &  
             / ((zc(10,10,TP+1)-zc(10,10,TP+2))*(-1000))
!      write(6,*) "slope = ", slope
!      write(6,*) "zc1010TP= ", zc(10,10,TP)*1000
!      write(6,*) "rho1010TP= ", s(10,10,TP,0)

      yinter =  slope * zc(10,10,TP+1)*(1000) + s(10,10,TP+1,0)      
!      write(6,*) "yinter = ", yinter
      do j=0,NJ+1
         do i=0,NI+1
            do k=0,TP+1  
                z= -zc(i,j,k)*1000
                  sbkgrnd = slope * z + yinter

!- slope * zc(10,10,TP+1)*(1000) & 
!                 - s(10,10,TP+1,0)
                  s(i,j,k,0) = sbkgrnd
!                  write(6,*) sbkgrnd

            end do
         end do
      end do



!     ____________________
!     --------------------
!     ADD FRONT


!     A frontal region is introduced to the channel.
!     The front extends to a depth of z2. The front
!     is strongest to depth z1 and decreases in 
!     strength from z1 to z2

      if (trigger.eq.0) then
 
      do j=0,NJ+1
         do i=0,NI+1
 
!     WIGGLE in FRONT
            wiggles=1.d0
            amplitude= 0.1d0
            yfront= 0.5d0*dble(yc(NJ+1) +yc(0))       & 
                 + amplitude*dsin(2.d0*PI*xc(i)/  & 
                 (0.5*dble(xc(NI+1)+xc(NI))/wiggles)  )
           thy = tanh(tightness*dble(yc(j)-yfront)*PI)
           do k=1,NK
               z= DL*zc(i,j,k)
!              if (z.gt.-2*mldepth) then
                  if (z.ge.-z1) then 
                     slfacnew = slfac
                  else if (z .ge. -z2) then
!                    slfacnew = slfac*exp(4.2d0*(z+z1)/300.d0) 
                     slfacnew = slfac*(z+z2)/(z2-z1)
                  else
                    slfacnew = 0.d0 
                  end if
                  s(i,j,k,0)= slfacnew*(thy-25d-1) + s(i,NJ,k,0)
!               end if
            end do
         end do
      end do

      else 
        if (trigger.gt.0) goto 101

      end if

101   continue

!     ___________________
!     -------------------
!     SET BOUNDARY CONDITIONS

!     North-south boundaries set to solid walls.
!     East-west boundaries are periodic.

      do k=0,NK+1
         do i=1,NI
            s(i,0,k,0)= s(i,1,k,0)
            s(i,NJ+1,k,0)= s(i,NJ,k,0)
         end do
         do j=0,NJ+1
            s(0,j,k,0)= s(NI,j,k,0)
            s(NI+1,j,k,0)= s(1,j,k,0)
         end do
      end do
!     ____________________


      return
      end

