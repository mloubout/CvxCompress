subroutine updateParticleVelocity(orderspace, dti, &
                   & dx, dy, dz, ilo, jlo, klo, iup, jup, kup, rho, &
                   & deta, vx, vy, vz, txx, tyy, tzz, tyz, txz, txy, &
                   & itausig, difitau, sx, sy, sz) 
 
!  Function for updating the particle velocities on a 3D staggered grid 
!  using O(16) spatial differencing. SLS viscoelasticity is applied on the
!  particle velocities.
!
!
!  o Last Modified:  
!    10-04-12  Adapted from vupdate2d.f90 (file:///devl/FD/Svn/Ssg2d/trunk)
!
!  o Written by Kurt T. Nihei


implicit none

integer :: orderspace
integer :: ilo, jlo, klo
integer :: iup, jup, kup
real :: dti
real :: dx, dy, dz
real, dimension(ilo:iup,jlo:jup,klo:kup) :: rho, deta, vx, vy, vz, &
& txx, tyy, tzz, tyz, txz, txy, itausig, difitau, sx, sy, sz

integer :: i, j, k
integer :: ii, jj, kk
integer :: ii_cnt, jj_cnt, kk_cnt
integer :: imin, jmin, kmin
integer :: imax, jmax, kmax
real, parameter :: c0 =  41409225./33554432. ! O(16) SG FD constants
real, parameter :: c1 = -3578575./33554432. 
real, parameter :: c2 =  3864861./167772160. 
real, parameter :: c3 = -1254825./234881024.
real, parameter :: c4 =  325325./301989888.
real, parameter :: c5 = -61425./369098752.
real, parameter :: c6 =  7425./436207616.
real, parameter :: c7 = -143./167772160.
real :: idx, idy, idz
real :: dxtxx, dytyy, dztzz 
real :: dytyz, dztyz 
real :: dxtxz, dztxz 
real :: dxtxy, dytxy 
real :: dabc
real :: const1, const2, const3
real, dimension(16,4) :: func  ! func vals for O(4) SG interp & O(16) SG FD
                               ! (16 lines each with 4 func vals)
integer, parameter :: iblk = 256
integer, parameter :: jblk = 10
integer :: inb, jnb, nnb
integer :: ib, i0, i1, jb, j0, j1
integer :: iib
 
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

!$OMP PARALLEL DO PRIVATE(i, j, k, dxtxx, dytyy, dztzz, &
!$OMP& dytyz, dztyz, dxtxz, dztxz, dxtxy, dytxy, &
!$OMP& const1, const2, const3, dabc, jb, ib, j0, j1, i0, i1)
do iib = 1,nnb
  jb = (iib-1)/inb
  ib = iib-1-jb*inb

  j0 = jmin + jb*jblk
  j1 = j0+jblk-1
  if (j1 > jmax) j1 = jmax
  i0 = imin + ib*iblk
  i1 = i0+iblk-1
  if (i1 > imax) i1 = imax

  do k = kmin,kmax
    do j = j0,j1
      do i = i0,i1

        ! Compute normal stress gradients (these are properly centered): 
        dxtxx  = (c0*(txx(i,j,k)   - txx(i-1,j,k)) + &
                & c1*(txx(i+1,j,k) - txx(i-2,j,k)) + &
                & c2*(txx(i+2,j,k) - txx(i-3,j,k)) + &
                & c3*(txx(i+3,j,k) - txx(i-4,j,k)) + &
                & c4*(txx(i+4,j,k) - txx(i-5,j,k)) + &
                & c5*(txx(i+5,j,k) - txx(i-6,j,k)) + &
                & c6*(txx(i+6,j,k) - txx(i-7,j,k)) + &
                & c7*(txx(i+7,j,k) - txx(i-8,j,k)) )*idx

        dytyy =  (c0*(tyy(i,j+1,k) - tyy(i,j,k))   + &
                & c1*(tyy(i,j+2,k) - tyy(i,j-1,k)) + &
                & c2*(tyy(i,j+3,k) - tyy(i,j-2,k)) + &
                & c3*(tyy(i,j+4,k) - tyy(i,j-3,k)) + &
                & c4*(tyy(i,j+5,k) - tyy(i,j-4,k)) + &
                & c5*(tyy(i,j+6,k) - tyy(i,j-5,k)) + &
                & c6*(tyy(i,j+7,k) - tyy(i,j-6,k)) + &
                & c7*(tyy(i,j+8,k) - tyy(i,j-7,k)) )*idy

        ! (note: minus sign is for pos. z-axis up)
        dztzz = -(c0*(tzz(i,j,k+1) - tzz(i,j,k))   + &
                & c1*(tzz(i,j,k+2) - tzz(i,j,k-1)) + &
                & c2*(tzz(i,j,k+3) - tzz(i,j,k-2)) + &
                & c3*(tzz(i,j,k+4) - tzz(i,j,k-3)) + &
                & c4*(tzz(i,j,k+5) - tzz(i,j,k-4)) + &
                & c5*(tzz(i,j,k+6) - tzz(i,j,k-5)) + &
                & c6*(tzz(i,j,k+7) - tzz(i,j,k-6)) + &
                & c7*(tzz(i,j,k+8) - tzz(i,j,k-7)) )*idz


        ! Compute shear stress gradients (these are properly centered
        ! for non-tilted anisotropy cases): 
        dytxy =  (c0*(txy(i,j,k)   - txy(i,j-1,k)) + &
                & c1*(txy(i,j+1,k) - txy(i,j-2,k)) + &
                & c2*(txy(i,j+2,k) - txy(i,j-3,k)) + &
                & c3*(txy(i,j+3,k) - txy(i,j-4,k)) + &
                & c4*(txy(i,j+4,k) - txy(i,j-5,k)) + &
                & c5*(txy(i,j+5,k) - txy(i,j-6,k)) + &
                & c6*(txy(i,j+6,k) - txy(i,j-7,k)) + &
                & c7*(txy(i,j+7,k) - txy(i,j-8,k)) )*idy

        ! (note: minus sign is for pos. z-axis up)
        dztxz = -(c0*(txz(i,j,k)   - txz(i,j,k-1)) + &
                & c1*(txz(i,j,k+1) - txz(i,j,k-2)) + &
                & c2*(txz(i,j,k+2) - txz(i,j,k-3)) + &
                & c3*(txz(i,j,k+3) - txz(i,j,k-4)) + &
                & c4*(txz(i,j,k+4) - txz(i,j,k-5)) + &
                & c5*(txz(i,j,k+5) - txz(i,j,k-6)) + &
                & c6*(txz(i,j,k+6) - txz(i,j,k-7)) + &
                & c7*(txz(i,j,k+7) - txz(i,j,k-8)) )*idz

        dxtxy =  (c0*(txy(i+1,j,k) - txy(i,j,k))   + &
                & c1*(txy(i+2,j,k) - txy(i-1,j,k)) + &
                & c2*(txy(i+3,j,k) - txy(i-2,j,k)) + &
                & c3*(txy(i+4,j,k) - txy(i-3,j,k)) + &
                & c4*(txy(i+5,j,k) - txy(i-4,j,k)) + &
                & c5*(txy(i+6,j,k) - txy(i-5,j,k)) + &
                & c6*(txy(i+7,j,k) - txy(i-6,j,k)) + &
                & c7*(txy(i+8,j,k) - txy(i-7,j,k)) )*idx

        ! (note: minus sign is for pos. z-axis up)
        dztyz = -(c0*(tyz(i,j,k)   - tyz(i,j,k-1)) + &
                & c1*(tyz(i,j,k+1) - tyz(i,j,k-2)) + &
                & c2*(tyz(i,j,k+2) - tyz(i,j,k-3)) + &
                & c3*(tyz(i,j,k+3) - tyz(i,j,k-4)) + &
                & c4*(tyz(i,j,k+4) - tyz(i,j,k-5)) + &
                & c5*(tyz(i,j,k+5) - tyz(i,j,k-6)) + &
                & c6*(tyz(i,j,k+6) - tyz(i,j,k-7)) + &
                & c7*(tyz(i,j,k+7) - tyz(i,j,k-8)) )*idz

        dytyz =  (c0*(tyz(i,j,k)   - tyz(i,j-1,k)) + &
                & c1*(tyz(i,j+1,k) - tyz(i,j-2,k)) + &
                & c2*(tyz(i,j+2,k) - tyz(i,j-3,k)) + &
                & c3*(tyz(i,j+3,k) - tyz(i,j-4,k)) + &
                & c4*(tyz(i,j+4,k) - tyz(i,j-5,k)) + &
                & c5*(tyz(i,j+5,k) - tyz(i,j-6,k)) + &
                & c6*(tyz(i,j+6,k) - tyz(i,j-7,k)) + &
                & c7*(tyz(i,j+7,k) - tyz(i,j-8,k)) )*idy

        dxtxz =  (c0*(txz(i+1,j,k) - txz(i,j,k))   + &
                & c1*(txz(i+2,j,k) - txz(i-1,j,k)) + &
                & c2*(txz(i+3,j,k) - txz(i-2,j,k)) + &
                & c3*(txz(i+4,j,k) - txz(i-3,j,k)) + &
                & c4*(txz(i+5,j,k) - txz(i-4,j,k)) + &
                & c5*(txz(i+6,j,k) - txz(i-5,j,k)) + &
                & c6*(txz(i+7,j,k) - txz(i-6,j,k)) + &
                & c7*(txz(i+8,j,k) - txz(i-7,j,k)) )*idx


        !   Update viscoelastic(SLS) vector field: 
        const1 = 1./(1. + 0.5*dti*itausig(i,j,k))
        const2 = (1. - 0.5*dti*itausig(i,j,k))
        const3 = dti*difitau(i,j,k)

        sx(i,j,k) = const1*(const2*sx(i,j,k) + &
                  & const3*(dxtxx + dytxy + dztxz))
        sy(i,j,k) = const1*(const2*sy(i,j,k) + &
                  & const3*(dxtxy + dytyy + dztyz))
        sz(i,j,k) = const1*(const2*sz(i,j,k) + &
                  & const3*(dxtxz + dytyz + dztzz))
                                  
        ! Absorbing boundary decay funct (for Maxwell viscoelastic model):
        dabc = (1 - 0.5*deta(i,j,k)*dti)/(1 + 0.5*deta(i,j,k)*dti)

        ! Update viscoelastic) particle velocities:
        vx(i,j,k) = dabc*vx(i,j,k) + dti/rho(i,j,k)*(dxtxx + &
                                   & dytxy + dztxz + sx(i,j,k))
        vy(i,j,k) = dabc*vy(i,j,k) + dti/rho(i,j,k)*(dxtxy + &
                                   & dytyy + dztyz + sy(i,j,k))
        vz(i,j,k) = dabc*vz(i,j,k) + dti/rho(i,j,k)*(dxtxz + &
                                   & dytyz + dztzz + sz(i,j,k))

      enddo
    enddo
  enddo
enddo
!$OMP END PARALLEL DO


return
end subroutine
