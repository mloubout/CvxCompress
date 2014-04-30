subroutine updateStressElaTri(orderspace, dti, dx, dy, dz, &
          & ilo, jlo, klo, iup, jup, kup, &
          & c11, c22, c33, c44, c55, c66, &
          & c12, c13, c14, c15, c16, &
          & c23, c24, c25, c26, &
          & c34, c35, c36, &
          & c45, c46, &
          & c56, &
          & deta, vx, vy, vz, txx, tyy, tzz, tyz, txz, txy)

!  Subroutine updates the stresses on a 3D staggered grid using 
!  finite differencing. Constitutive law is for a viscoelastic Triclinic
!  medium (21 cIJ's). Stresses and strains are specified in the
!  global (FD) frame.
!  
!
!
!  o Last Modified:  
!    10-04-12  Adapted from trupdate2d.f90 (file:///devl/FD/Svn/Ssg2d/trunk)
!
!  o Written by Kurt T. Nihei


implicit none

integer :: orderspace
integer :: ilo, jlo, klo
integer :: iup, jup, kup
real :: dti
real :: dx, dy, dz
real, dimension(ilo:iup,jlo:jup,klo:kup) :: c11, c22, c33, c44, c55, c66, &
& c12, c13, c14, c15, c16, &
& c23, c24, c25, c26, &
& c34, c35, c36, &
& c45, c46, &
& c56, &
& deta, vx, vy, vz, txx, tyy, tzz, tyz, txz, txy

integer :: i, j, k
integer :: imin, jmin, kmin
integer :: imax, jmax, kmax
real :: idx, idy, idz
real :: dabc
real :: dtexx, dteyy, dtezz, dteyz2, dtexz2, dtexy2
real :: dtexx_tyz, dtexx_txz, dtexx_txy 
real :: dteyy_tyz, dteyy_txz, dteyy_txy 
real :: dtezz_tyz, dtezz_txz, dtezz_txy 
real :: dteyz2_tkk, dteyz2_txz, dteyz2_txy 
real :: dtexz2_tkk, dtexz2_tyz, dtexz2_txy 
real :: dtexy2_tkk, dtexy2_tyz, dtexy2_txz


imin = ilo + orderspace/2
jmin = jlo + orderspace/2
kmin = klo + orderspace/2

imax = iup - orderspace/2
jmax = jup - orderspace/2
kmax = kup - orderspace/2

idx = 1/dx
idy = 1/dy
idz = 1/dz


!$OMP PARALLEL DO PRIVATE(i, j, dabc, &
!$OMP& dtexx, dteyy, dtezz, dteyz2, dtexz2, dtexy2, &
!$OMP& dtexx_tyz, dtexx_txz, dtexx_txy, & 
!$OMP& dteyy_tyz, dteyy_txz, dteyy_txy, & 
!$OMP& dtezz_tyz, dtezz_txz, dtezz_txy, & 
!$OMP& dteyz2_tkk, dteyz2_txz, dteyz2_txy, & 
!$OMP& dtexz2_tkk, dtexz2_tyz, dtexz2_txy, & 
!$OMP& dtexy2_tkk, dtexy2_tyz, dtexy2_txz)
do k = kmin, kmax
  do j = jmin, jmax
!DIR$ VECTOR ALWAYS
    do i = imin, imax

      ! Compute (6) centered strain rates using O(16) spatial FD:
      call computeStrainRateElaOnKite(i, j, k, idx, idy, idz, &
                   & ilo, jlo, klo, iup, jup, kup, &
                   & vx, vy, vz, &
                   & dtexx, dteyy, dtezz, dteyz2, dtexz2, dtexy2)

      ! Compute (18) non-centered strain rates using O(16) spatial FD:
      call computeStrainRateElaOffKite(i, j, k, idx, idy, idz, &
                   & ilo, jlo, klo, iup, jup, kup, &
                   & vx, vy, vz, &
                   & dtexx_tyz, dtexx_txz, dtexx_txy, & 
                   & dteyy_tyz, dteyy_txz, dteyy_txy, & 
                   & dtezz_tyz, dtezz_txz, dtezz_txy, & 
                   & dteyz2_tkk, dteyz2_txz, dteyz2_txy, & 
                   & dtexz2_tkk, dtexz2_tyz, dtexz2_txy, & 
                   & dtexy2_tkk, dtexy2_tyz, dtexy2_txz)


      ! Absorbing boundary decay function (for Maxwell viscoelastic model):
      dabc = (1 - 0.5*deta(i,j,k)*dti)/(1 + 0.5*deta(i,j,k)*dti)


      ! Compute stresses for triclinic medium:   
      txx(i,j,k) = dabc*txx(i,j,k) + dti*( &
                 & c11(i,j,k)*dtexx + &
                 & c12(i,j,k)*dteyy + &
                 & c13(i,j,k)*dtezz + &
                 & c14(i,j,k)*dteyz2_tkk + &
                 & c15(i,j,k)*dtexz2_tkk + &
                 & c16(i,j,k)*dtexy2_tkk)

      tyy(i,j,k) = dabc*tyy(i,j,k) + dti*( &
                 & c12(i,j,k)*dtexx + &
                 & c22(i,j,k)*dteyy + &
                 & c23(i,j,k)*dtezz + &
                 & c24(i,j,k)*dteyz2_tkk + &
                 & c25(i,j,k)*dtexz2_tkk + &
                 & c26(i,j,k)*dtexy2_tkk)

      tzz(i,j,k) = dabc*tzz(i,j,k) + dti*( &
                 & c13(i,j,k)*dtexx + &
                 & c23(i,j,k)*dteyy + &
                 & c33(i,j,k)*dtezz + &
                 & c34(i,j,k)*dteyz2_tkk + &
                 & c35(i,j,k)*dtexz2_tkk + &
                 & c36(i,j,k)*dtexy2_tkk)

      tyz(i,j,k) = dabc*tyz(i,j,k) + dti*( &
                 & c14(i,j,k)*dtexx_tyz + &
                 & c24(i,j,k)*dteyy_tyz + &
                 & c34(i,j,k)*dtezz_tyz + &
                 & c44(i,j,k)*dteyz2    + &
                 & c45(i,j,k)*dtexz2_tyz + &
                 & c46(i,j,k)*dtexy2_tyz)

      txz(i,j,k) = dabc*txz(i,j,k) + dti*( &
                 & c15(i,j,k)*dtexx_txz + &
                 & c25(i,j,k)*dteyy_txz + &
                 & c35(i,j,k)*dtezz_txz + &
                 & c45(i,j,k)*dteyz2_txz + &
                 & c55(i,j,k)*dtexz2     + &
                 & c56(i,j,k)*dtexy2_txz)

      txy(i,j,k) = dabc*txy(i,j,k) + dti*( &
                 & c16(i,j,k)*dtexx_txy + &
                 & c26(i,j,k)*dteyy_txy + &
                 & c36(i,j,k)*dtezz_txy + &
                 & c46(i,j,k)*dteyz2_txy + &
                 & c56(i,j,k)*dtexz2_txy + &
                 & c66(i,j,k)*dtexy2    )

    enddo
  enddo
enddo
!$OMP END PARALLEL DO


return
end subroutine
