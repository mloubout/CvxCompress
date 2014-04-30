subroutine updateStressElaVOr(orderspace, dti, dx, dy, dz, &
                   & ilo, jlo, klo, iup, jup, kup, &
                   & c11, c22, c33, c44, c55, c66, c12, c13, c23, &
                   & deta, vx, vy, vz, txx, tyy, tzz, tyz, txz, txy)

!  Function for updating the stresses on a 3D staggered grid using
!  finite differencing. Constitutive law is for an elastic Vertical 
!  Orthorhombic (VOr) medium. 
!
!
!  o Last Modified:  
!    01-25-13  Adapted from trupdate2d.f90 (file:///devl/FD/Svn/Ssg2d/trunk)
!
!  o Written by Kurt T. Nihei


implicit none

integer :: orderspace
integer :: ilo, jlo, klo
integer :: iup, jup, kup
real :: dti
real :: dx, dy, dz
real, dimension(ilo:iup,jlo:jup,klo:kup) :: c11, c22, c33, c44, c55, c66, &
& c12, c13, c23, deta, vx, vy, vz, txx, tyy, tzz, tyz, txz, txy

integer :: i, j, k
integer :: imin, jmin, kmin
integer :: imax, jmax, kmax
real :: idx, idy, idz
real :: dtexx, dteyy, dtezz
real :: dteyz2, dtexz2, dtexy2
real :: dabc

imin = ilo + orderspace/2
jmin = jlo + orderspace/2
kmin = klo + orderspace/2

imax = iup - orderspace/2
jmax = jup - orderspace/2
kmax = kup - orderspace/2

idx = 1/dx
idy = 1/dy
idz = 1/dz


!$OMP PARALLEL DO PRIVATE(i, j, dtexx, dteyy, dtezz, dteyz2, dtexz2, &
!$OMP& dtexy2, dabc)
do k = kmin, kmax
  do j = jmin, jmax
    do i = imin, imax

      ! Compute strain rates using O(16) spatial FD:
      call computeStrainRateElaOnKite(i, j, k, idx, idy, idz, &
                   & ilo, jlo, klo, iup, jup, kup, vx, vy, vz, &
                   & dtexx, dteyy, dtezz, dteyz2, dtexz2, dtexy2)

      ! Absorbing boundary decay function (for Maxwell viscoelastic model):
      dabc = (1 - 0.5*deta(i,j,k)*dti)/(1 + 0.5*deta(i,j,k)*dti)

      ! Compute stresses for an elastic VOr medium:     
      txx(i,j,k) = dabc*txx(i,j,k) + dti*(c11(i,j,k)*dtexx + &
                      & c12(i,j,k)*dteyy + c13(i,j,k)*dtezz)
      tyy(i,j,k) = dabc*tyy(i,j,k) + dti*(c12(i,j,k)*dtexx + &
                      & c22(i,j,k)*dteyy + c23(i,j,k)*dtezz)
      tzz(i,j,k) = dabc*tzz(i,j,k) + dti*(c13(i,j,k)*dtexx + &
                      & c23(i,j,k)*dteyy + c33(i,j,k)*dtezz)
      tyz(i,j,k) = dabc*tyz(i,j,k) + dti*c44(i,j,k)*dteyz2
      txz(i,j,k) = dabc*txz(i,j,k) + dti*c55(i,j,k)*dtexz2
      txy(i,j,k) = dabc*txy(i,j,k) + dti*c66(i,j,k)*dtexy2
    enddo
  enddo
enddo
!$OMP END PARALLEL DO


return
end subroutine
