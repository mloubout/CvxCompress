subroutine computeStrainRateElaOnKite(i, j, k, idx, idy, idz, &
                   & ilo, jlo, klo, iup, jup, kup, &
                   & vx, vy, vz, &
                   & dtexx, dteyy, dtezz, dteyz2, dtexz2, dtexy2)

!  Subroutine computes the (6) centered strain rates on a staggered 
!  grid using O(16) spatial differencing. 
!
!  o Last Modified:  
!    01-25-13  Adapted from trupdate2d.f90 (file:///devl/FD/Svn/Ssg2d/trunk)
!
!  o Written by Kurt T. Nihei


implicit none

integer :: i, j, k
integer :: ilo, jlo, klo
integer :: iup, jup, kup
real :: idx, idy, idz
real, dimension(ilo:iup,jlo:jup,klo:kup) :: vx, vy, vz
real :: dtexx, dteyy, dtezz, dteyz2, dtexz2, dtexy2

real, parameter :: c0 =  41409225./33554432. ! O(16) SG FD constants
real, parameter :: c1 = -3578575./33554432. 
real, parameter :: c2 =  3864861./167772160. 
real, parameter :: c3 = -1254825./234881024.
real, parameter :: c4 =  325325./301989888.
real, parameter :: c5 = -61425./369098752.
real, parameter :: c6 =  7425./436207616.
real, parameter :: c7 = -143./167772160.
real :: dyvz, dzvy, dzvx, dxvz, dyvx, dxvy


! Compute O(16) normal strain rates:
dtexx =  (c0*(vx(i+1,j,k) - vx(i,j,k))   + &
        & c1*(vx(i+2,j,k) - vx(i-1,j,k)) + &
        & c2*(vx(i+3,j,k) - vx(i-2,j,k)) + &
        & c3*(vx(i+4,j,k) - vx(i-3,j,k)) + &
        & c4*(vx(i+5,j,k) - vx(i-4,j,k)) + &
        & c5*(vx(i+6,j,k) - vx(i-5,j,k)) + &
        & c6*(vx(i+7,j,k) - vx(i-6,j,k)) + &
        & c7*(vx(i+8,j,k) - vx(i-7,j,k)) )*idx

dteyy =  (c0*(vy(i,j,k)   - vy(i,j-1,k)) + &
        & c1*(vy(i,j+1,k) - vy(i,j-2,k)) + &
        & c2*(vy(i,j+2,k) - vy(i,j-3,k)) + &
        & c3*(vy(i,j+3,k) - vy(i,j-4,k)) + &
        & c4*(vy(i,j+4,k) - vy(i,j-5,k)) + &
        & c5*(vy(i,j+5,k) - vy(i,j-6,k)) + &
        & c6*(vy(i,j+6,k) - vy(i,j-7,k)) + &
        & c7*(vy(i,j+7,k) - vy(i,j-8,k)) )*idy

! (note: minus sign is for pos. z-axis up)
dtezz = -(c0*(vz(i,j,k)   - vz(i,j,k-1)) + &
        & c1*(vz(i,j,k+1) - vz(i,j,k-2)) + &    
        & c2*(vz(i,j,k+2) - vz(i,j,k-3)) + &
        & c3*(vz(i,j,k+3) - vz(i,j,k-4)) + &
        & c4*(vz(i,j,k+4) - vz(i,j,k-5)) + &
        & c5*(vz(i,j,k+5) - vz(i,j,k-6)) + &
        & c6*(vz(i,j,k+6) - vz(i,j,k-7)) + &
        & c7*(vz(i,j,k+7) - vz(i,j,k-8)) )*idz

! Compute O(16) shear strain rates:
! (note: minus sign is for pos. z-axis up)
dzvy = -(c0*(vy(i,j,k+1) - vy(i,j,k))   + &
       & c1*(vy(i,j,k+2) - vy(i,j,k-1)) + &    
       & c2*(vy(i,j,k+3) - vy(i,j,k-2)) + &
       & c3*(vy(i,j,k+4) - vy(i,j,k-3)) + &
       & c4*(vy(i,j,k+5) - vy(i,j,k-4)) + &
       & c5*(vy(i,j,k+6) - vy(i,j,k-5)) + &
       & c6*(vy(i,j,k+7) - vy(i,j,k-6)) + &
       & c7*(vy(i,j,k+8) - vy(i,j,k-7)) )*idz

dyvz =  (c0*(vz(i,j+1,k) - vz(i,j,k))   + &
       & c1*(vz(i,j+2,k) - vz(i,j-1,k)) + &
       & c2*(vz(i,j+3,k) - vz(i,j-2,k)) + &
       & c3*(vz(i,j+4,k) - vz(i,j-3,k)) + &
       & c4*(vz(i,j+5,k) - vz(i,j-4,k)) + &
       & c5*(vz(i,j+6,k) - vz(i,j-5,k)) + &
       & c6*(vz(i,j+7,k) - vz(i,j-6,k)) + &
       & c7*(vz(i,j+8,k) - vz(i,j-7,k)) )*idy

dteyz2 = dzvy + dyvz

! (note: minus sign is for pos. z-axis up)
dzvx = -(c0*(vx(i,j,k+1) - vx(i,j,k))   + &
       & c1*(vx(i,j,k+2) - vx(i,j,k-1)) + &    
       & c2*(vx(i,j,k+3) - vx(i,j,k-2)) + &
       & c3*(vx(i,j,k+4) - vx(i,j,k-3)) + &
       & c4*(vx(i,j,k+5) - vx(i,j,k-4)) + &
       & c5*(vx(i,j,k+6) - vx(i,j,k-5)) + &
       & c6*(vx(i,j,k+7) - vx(i,j,k-6)) + &
       & c7*(vx(i,j,k+8) - vx(i,j,k-7)) )*idz

dxvz =  (c0*(vz(i,j,k)   - vz(i-1,j,k)) + &
       & c1*(vz(i+1,j,k) - vz(i-2,j,k)) + &
       & c2*(vz(i+2,j,k) - vz(i-3,j,k)) + &
       & c3*(vz(i+3,j,k) - vz(i-4,j,k)) + &
       & c4*(vz(i+4,j,k) - vz(i-5,j,k)) + &
       & c5*(vz(i+5,j,k) - vz(i-6,j,k)) + &
       & c6*(vz(i+6,j,k) - vz(i-7,j,k)) + &
       & c7*(vz(i+7,j,k) - vz(i-8,j,k)) )*idx

dtexz2 = dzvx + dxvz
    
dyvx =  (c0*(vx(i,j+1,k) - vx(i,j,k))   + &
       & c1*(vx(i,j+2,k) - vx(i,j-1,k)) + &
       & c2*(vx(i,j+3,k) - vx(i,j-2,k)) + &
       & c3*(vx(i,j+4,k) - vx(i,j-3,k)) + &
       & c4*(vx(i,j+5,k) - vx(i,j-4,k)) + &
       & c5*(vx(i,j+6,k) - vx(i,j-5,k)) + &
       & c6*(vx(i,j+7,k) - vx(i,j-6,k)) + &
       & c7*(vx(i,j+8,k) - vx(i,j-7,k)) )*idy

dxvy =  (c0*(vy(i,j,k)   - vy(i-1,j,k)) + &
       & c1*(vy(i+1,j,k) - vy(i-2,j,k)) + &
       & c2*(vy(i+2,j,k) - vy(i-3,j,k)) + &
       & c3*(vy(i+3,j,k) - vy(i-4,j,k)) + &
       & c4*(vy(i+4,j,k) - vy(i-5,j,k)) + &
       & c5*(vy(i+5,j,k) - vy(i-6,j,k)) + &
       & c6*(vy(i+6,j,k) - vy(i-7,j,k)) + &
       & c7*(vy(i+7,j,k) - vy(i-8,j,k)) )*idx

dtexy2 = dyvx + dxvy


return
end subroutine
