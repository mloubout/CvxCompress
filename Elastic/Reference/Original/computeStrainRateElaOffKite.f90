subroutine computeStrainRateElaOffKite(i, j, k, idx, idy, idz, &
                   & ilo, jlo, klo, iup, jup, kup, &
                   & vx, vy, vz, &
                   & dtexx_tyz, dtexx_txz, dtexx_txy, &
                   & dteyy_tyz, dteyy_txz, dteyy_txy, &
                   & dtezz_tyz, dtezz_txz, dtezz_txy, &
                   & dteyz2_tkk, dteyz2_txz, dteyz2_txy, &
                   & dtexz2_tkk, dtexz2_tyz, dtexz2_txy, &
                   & dtexy2_tkk, dtexy2_tyz, dtexy2_txz)

!  Subroutine computes the "off-kite" strain rates(i.e., the 18 strain 
!  rates that multiply the cIJ's in the upper 3x3 and lower diagnol).
!  The 18 off-kite strain rates are not properly centered on the Yee
!  staggered grid, and, thus, require interpolation of the particle
!  velocities before the strain rates can be computed. Bi-linear 
!  interpolation (O(2) accuracy) is used to approximate the
!  value of the non-centered particle velocity on the line it is needed
!  for the derivative.  Derivatives are then computed using O(8) staggered
!  grid spatial differencing. 
!
!  o Last Modified:  
!    01-25-13  Adapted from te_up3d.f90 (/users/knih/NumericalCodes/Ssg3d/
!              Svn/Pre2012/Trunk)
!
!  o Written by Kurt T. Nihei


implicit none

integer :: i, j, k
integer :: ilo, jlo, klo
integer :: iup, jup, kup
real :: idx, idy, idz
real, dimension(ilo:iup,jlo:jup,klo:kup) :: vx, vy, vz
real :: dtexx_tyz, dtexx_txz, dtexx_txy 
real :: dteyy_tyz, dteyy_txz, dteyy_txy 
real :: dtezz_tyz, dtezz_txz, dtezz_txy 
real :: dteyz2_tkk, dteyz2_txz, dteyz2_txy 
real :: dtexz2_tkk, dtexz2_tyz, dtexz2_txy 
real :: dtexy2_tkk, dtexy2_tyz, dtexy2_txz

real, parameter :: c0 =  1225./1024. ! O(8) SG FD constants
real, parameter :: c1 = -245./3072.
real, parameter :: c2 =  49./5120.
real, parameter :: c3 = -5./7168.

integer :: ii, jj, kk
integer :: icnt, jcnt, kcnt
real :: dv_1, dv_2
real, dimension(8) :: v


! Compute off-kite strain rates. For a trinclinic medium, there are
! 18 off-kite strain rates, corresponding to 27 non-centered-staggered-grid
! particle velocity spatial derivatives. First, the particle velocities 
! are interpolated on to a line using 4 point averaging (bi-linear 
! interpolation). Second, the spatial derivative of the particle
! velocity are computed using O(8) staggerd grid FD.

! ..tkk: 2 deyz/dt (off-kite strain rate #1/18)
kcnt = 0
do kk = k-3, k+4
  kcnt = kcnt + 1
  v(kcnt) = vy(i,j-1,kk) + vy(i,j,kk) + vy(i,j,kk-1) + vy(i,j-1,kk-1)
enddo
dv_1 = -(c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idz  ! minus sign 

jcnt = 0
do jj = j-4, j+3
  jcnt = jcnt + 1
  v(jcnt) = vz(i,jj,k) + vz(i,jj+1,k) + vz(i,jj+1,k-1) + vz(i,jj,k-1)
enddo
dv_2 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idy

dteyz2_tkk = 0.25*(dv_1 + dv_2)
 

! ..tkk: 2 dexz/dt (off-kite strain rate #2/18)
kcnt = 0
do kk = k-3, k+4
  kcnt = kcnt + 1
  v(kcnt) = vx(i,j,kk) + vx(i+1,j,kk) + vx(i+1,j,kk-1) + vx(i,j,kk-1)
enddo
dv_1 = -(c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idz  ! minus sign 

icnt = 0
do ii = i-4, i+3
  icnt = icnt + 1
  v(icnt) = vz(ii,j,k) + vz(ii+1,j,k) + vz(ii+1,j,k-1) + vz(ii,j,k-1)
enddo
dv_2 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idx 

dtexz2_tkk = 0.25*(dv_1 + dv_2)
 

! ..tkk: 2 dexy/dt (off-kite strain rate #3/18)
jcnt = 0
do jj = j-3, j+4
  jcnt = jcnt + 1
  v(jcnt) = vx(i,jj,k) + vx(i+1,jj,k) + vx(i+1,jj-1,k) + vx(i,jj-1,k)
enddo
dv_1 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idy 
 
icnt = 0
do ii = i-4, i+3
  icnt = icnt + 1
  v(icnt) = vy(ii,j,k) + vy(ii+1,j,k) + vy(ii+1,j-1,k) + vy(ii,j-1,k)
enddo
dv_2 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idx 

dtexy2_tkk = 0.25*(dv_1 + dv_2)
 

! ..tyz: dexx/dt (off-kite strain rate #4/18)
icnt = 0
do ii = i-3, i+4
  icnt = icnt + 1
  v(icnt) = vx(ii,j,k+1) + vx(ii,j+1,k+1) + vx(ii,j+1,k) + vx(ii,j,k)
enddo
dv_1 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idx 

dtexx_tyz = 0.25*dv_1

 
! ..tyz: deyy/dt (off-kite strain rate #5/18)
jcnt = 0
do jj = j-4, j+3
  jcnt = jcnt + 1
  v(jcnt) = vy(i,jj,k+1) + vy(i,jj+1,k+1) + vy(i,jj+1,k) + vy(i,jj,k)
enddo
dv_1 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idy 

dteyy_tyz = 0.25*dv_1
 

! ..tyz: dezz/dt (off-kite strain rate #6/18)
kcnt = 0
do kk = k-3, k+4
  kcnt = kcnt + 1
  v(kcnt) = vz(i,j,kk) + vz(i,j+1,kk) + vz(i,j+1,kk-1) + vz(i,j,kk-1)
enddo
dv_1 = -(c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idz  ! minus sign 

dtezz_tyz = 0.25*dv_1
 

! ..tyz: 2 dexz/dt (off-kite strain rate #7/18)
kcnt = 0
do kk = k-3, k+4
  kcnt = kcnt + 1
  v(kcnt) = vx(i,j+1,kk) + vx(i+1,j+1,kk) + vx(i+1,j,kk) + vx(i,j,kk)
enddo
dv_1 = -(c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idz  ! minus sign 
 
icnt = 0
do ii = i-4, i+3
  icnt = icnt + 1
  v(icnt) = vz(ii,j+1,k) + vz(ii+1,j+1,k) + vz(ii+1,j,k) + vz(ii,j,k)
enddo
dv_2 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idx 
 
dtexz2_tyz = 0.25*(dv_1 + dv_2)
 

! ..tyz: 2exy (off-kite strain rate #8/18)
jcnt = 0
do jj = j-3, j+4
  jcnt = jcnt + 1
  v(jcnt) = vx(i,jj,k+1) + vx(i+1,jj,k+1) + vx(i+1,jj,k) + vx(i,jj,k)
enddo
dv_1 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idy 
 
icnt = 0
do ii = i-4, i+3
  icnt = icnt + 1
  v(icnt) = vy(ii,j,k+1) + vy(ii+1,j,k+1) + vy(ii+1,j,k) + vy(ii,j,k)
enddo
dv_2 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idx 

dtexy2_tyz = 0.25*(dv_1 + dv_2)
 

! ..txz: dexx/dt (off-kite strain rate #9/18)
icnt = 0
do ii = i-4, i+3
  icnt = icnt + 1
  v(icnt) = vx(ii,j,k+1) + vx(ii+1,j,k+1) + vx(ii+1,j,k) + vx(ii,j,k)
enddo
dv_1 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idx 

dtexx_txz = 0.25*dv_1
 

! ..txz: deyy/dt (off-kite strain rate #10/18)
jcnt = 0
do jj = j-4, j+3
  jcnt = jcnt + 1
  v(jcnt) = vy(i-1,jj,k+1) + vy(i,jj,k+1) + vy(i,jj,k) + vy(i-1,jj,k)
enddo
dv_1 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idy 

dteyy_txz = 0.25*dv_1
 

! ..txz: dezz/dt (off-kite strain rate #11/18)
kcnt = 0
do kk = k-3, k+4
  kcnt = kcnt + 1
  v(kcnt) = vz(i-1,j,kk) + vz(i,j,kk) + vz(i,j,kk-1) + vz(i-1,j,kk-1)
enddo
dv_1 = -(c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idz  ! minus sign 

dtezz_txz = 0.25*dv_1
 

! ..txz: 2 deyz/dt (off-kite strain rate #12/18)
kcnt = 0
do kk = k-3, k+4
  kcnt = kcnt + 1
  v(kcnt) = vy(i-1,j,kk) + vy(i,j,kk) + vy(i,j-1,kk) + vy(i-1,j-1,kk)
enddo
dv_1 = -(c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idz  ! minus sign 
 
jcnt = 0
do jj = j-3, j+4
  jcnt = jcnt + 1
  v(jcnt) = vz(i-1,jj,k) + vz(i,jj,k) + vz(i,jj-1,k) + vz(i-1,jj-1,k)
enddo
dv_2 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idy 
 
dteyz2_txz = 0.25*(dv_1 + dv_2)
 

! ..txz: 2 dexy/dt (off-kite strain rate #13/18)
jcnt = 0
do jj = j-4, j+3
  jcnt = jcnt + 1
  v(jcnt) = vx(i,jj,k+1) + vx(i,jj+1,k+1) + vx(i,jj+1,k) + vz(i,jj,k)
enddo
dv_1 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idy 
 
icnt = 0
do ii = i-4, i+3
  icnt = icnt + 1
  v(icnt) = vy(ii,j-1,k+1) + vy(ii,j,k+1) + vy(ii,j,k) + vy(ii,j-1,k)
enddo
dv_2 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idx 
 
dtexy2_txz = 0.25*(dv_1 + dv_2)
 

! ..txy: dexx/dt (off-kite strain rate #14/18)
icnt = 0
do ii = i-4, i+3
  icnt = icnt + 1
  v(icnt) = vx(ii,j+1,k) + vx(ii+1,j+1,k) + vx(ii+1,j,k) + vx(ii,j,k)
enddo
dv_1 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idx 

dtexx_txy = 0.25*dv_1
 

! ..txy: deyy/dt (off-kite strain rate #15/18)
jcnt = 0
do jj = j-3, j+4
  jcnt = jcnt + 1
  v(jcnt) = vy(i-1,jj,k) + vy(i,jj,k) + vy(i,jj-1,k) + vy(i-1,jj-1,k)
enddo
dv_1 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idy 

dteyy_txy = 0.25*dv_1
 
! ..txy: dezz/dt (off-kite strain rate #16/18)
kcnt = 0
do kk = k-4, k+3
  kcnt = kcnt + 1
  v(kcnt) = vz(i-1,j+1,kk) + vz(i,j+1,kk) + vz(i,j,kk) + vz(i-1,j,kk)
enddo
dv_1 = -(c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idz  ! minus sign 

dtezz_txy = 0.25*dv_1
 

! ..txy: 2 deyz/dt (off-kite strain rate #17/18)
kcnt = 0
do kk = k-3, k+4
  kcnt = kcnt + 1
  v(kcnt) = vy(i-1,j,kk) + vy(i,j,kk) + vy(i,j,kk-1) + vy(i-1,j,kk-1)
enddo
dv_1 = -(c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idz  ! minus sign 
 
jcnt = 0
do jj = j-3, j+4
  jcnt = jcnt + 1
  v(jcnt) = vz(i-1,jj,k) + vz(i,jj,k) + vz(i,jj,k-1) + vz(i-1,jj,k-1)
enddo
dv_2 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idy 
 
dteyz2_txy = 0.25*(dv_1 + dv_2)
 

! ..txy: 2 dexz/dt (off-kite strain rate #18/18)
kcnt = 0
do kk = k-3, k+4
  kcnt = kcnt + 1
  v(kcnt) = vx(i,j,kk) + vx(i,j+1,kk) + vx(i,j+1,kk-1) + vx(i,j,kk-1)
enddo
dv_1 = -(c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idz  ! minus sign 
 
icnt = 0
do ii = i-4, i+3
  icnt = icnt + 1
  v(icnt) = vz(ii,j,k) + vz(ii,j+1,k) + vz(ii,j+1,k-1) + vz(ii,j,k-1)
enddo
dv_2 =  (c0*(v(5) - v(4)) + &
       & c1*(v(6) - v(3)) + &
       & c2*(v(7) - v(2)) + &
       & c3*(v(8) - v(1)))*idx 
 
dtexz2_txy = 0.25*(dv_1 + dv_2)


return
end subroutine
