subroutine updateStressElaVTI(orderspace, dti, dx, dy, dz, ilo, jlo, klo, &
                   & iup, jup, kup, c11, c33, c55, c66, c13, &
                   & deta, vx, vy, vz, txx, tyy, tzz, tyz, txz, txy)

!  Function for updating the stresses on a 3D staggered grid using 
!  finte differencing. Constitutive law is for an elastic Vertical 
!  Transverse Isotropic (VTI) medium.
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
real, dimension(ilo:iup,jlo:jup,klo:kup) :: c11, c33, c55, c66, c13, &
& deta, vx, vy, vz, txx, tyy, tzz, tyz, txz, txy

integer :: i, j, k
integer :: imin, jmin, kmin
integer :: imax, jmax, kmax
real :: idx, idy, idz
real :: dtexx, dteyy, dtezz
real :: dteyz2, dtexz2, dtexy2
real :: dabc
real :: c12_

real, parameter :: c0 =  41409225./33554432. ! O(16) SG FD constants
real, parameter :: c1 = -3578575./33554432.
real, parameter :: c2 =  3864861./167772160.
real, parameter :: c3 = -1254825./234881024.
real, parameter :: c4 =  325325./301989888.
real, parameter :: c5 = -61425./369098752.
real, parameter :: c6 =  7425./436207616.
real, parameter :: c7 = -143./167772160.
real :: dyvz, dzvy, dzvx, dxvz, dyvx, dxvy

integer, parameter :: iblk = 256
integer, parameter :: jblk = 10
integer :: inb, jnb, nnb
integer :: ib, i0, i1, jb, j0, j1
integer :: ii

imin = ilo + orderspace/2
jmin = jlo + orderspace/2
kmin = klo + orderspace/2

imax = iup - orderspace/2
jmax = jup - orderspace/2
kmax = kup - orderspace/2

idx = 1/dx
idy = 1/dy
idz = 1/dz

inb = (imax-imin+iblk)/iblk
jnb = (jmax-jmin+jblk)/jblk
nnb = inb * jnb

!$OMP PARALLEL DO PRIVATE(i, j, dtexx, dteyy, dtezz, dteyz2, dtexz2, &
!$OMP& dtexy2, dabc, c12_)
do ii = 1,nnb
  jb = (ii-1)/inb
  ib = ii-1-jb*inb

  j0 = jmin + jb*jblk
  j1 = j0+jblk-1
  if (j1 > jmax) j1 = jmax
  i0 = imin + ib*iblk
  i1 = i0+iblk-1
  if (i1 > imax) i1 = imax

  do k = kmin,kmax
    do j = j0,j1
      do i = i0,i1

        ! Compute strain rates using O(16) spatial FD:
        include 'csrElaOnKite.f90'

        ! Absorbing boundary decay function (for Maxwell viscoelastic model):
        dabc = (1 - 0.5*deta(i,j,k)*dti)/(1 + 0.5*deta(i,j,k)*dti)

        ! Compute stresses for a VTI medium:     
        c12_ = c11(i,j,k) - 2.*c66(i,j,k)
        txx(i,j,k) = dabc*txx(i,j,k) + dti*(c11(i,j,k)*dtexx + &
                         & c12_*dteyy + c13(i,j,k)*dtezz)
        tyy(i,j,k) = dabc*tyy(i,j,k) + dti*(c12_*dtexx + &
                         & c11(i,j,k)*dteyy + c13(i,j,k)*dtezz)
        tzz(i,j,k) = dabc*tzz(i,j,k) + dti*(c13(i,j,k)*(dtexx + dteyy) + &
                         & c33(i,j,k)*dtezz)
        tyz(i,j,k) = dabc*tyz(i,j,k) + dti*c55(i,j,k)*dteyz2
        txz(i,j,k) = dabc*txz(i,j,k) + dti*c55(i,j,k)*dtexz2
        txy(i,j,k) = dabc*txy(i,j,k) + dti*c66(i,j,k)*dtexy2

      enddo
    enddo
  enddo
enddo
!$OMP END PARALLEL DO


return
end subroutine
