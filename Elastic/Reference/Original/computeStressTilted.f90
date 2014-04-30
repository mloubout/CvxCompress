subroutine computeStressTilted(orderspace, ilo, jlo, klo, iup, jup, kup, &
          & azim, dip, rake, txx, tyy, tzz, tyz, txz, txy, &
          & txx_loc, tyy_loc, tzz_loc, tyz_loc, txz_loc, txy_loc)

!  Subroutine computes the stresses in the local (TTI; TOr) frame
!  by collocation and Euler rotation of the global (FD) stresses. Local
!  stresses are used in FWI cIJ gradients.
!
!
!  o Last Modified:  
!    03-11-13 Created.
!
!  o Written by Kurt T. Nihei


implicit none

integer :: orderspace
integer :: ilo, jlo, klo
integer :: iup, jup, kup
real, dimension(ilo:iup,jlo:jup,klo:kup) :: azim, dip, rake, &
& txx, tyy, tzz, tyz, txz, txy, &
& txx_loc, tyy_loc, tzz_loc, tyz_loc, txz_loc, txy_loc

integer :: i, j, k
integer :: imin, jmin, kmin
integer :: imax, jmax, kmax
real :: dip_, azim_, rake_
real :: txx_, tyy_, tzz_, tyz_, txz_, txy_
real :: b(3,3)  ! Euler ZYZ rotation tensor
real :: bt(3,3)  ! Euler ZYZ rotation tensor x stress rate tensor
real :: cosa, cosd, cosr
real :: sina, sind, sinr


imin = ilo + orderspace/2
jmin = jlo + orderspace/2
kmin = klo + orderspace/2

imax = iup - orderspace/2
jmax = jup - orderspace/2
kmax = kup - orderspace/2


!$OMP PARALLEL DO PRIVATE(i, j, txx_, tyy_, tzz_, tyz_, txz_, txy_, &
!$OMP& azim_, dip_, rake_, cosa, cosd, cosr, sina, sind, sinr, b, bt)
do k = kmin, kmax
  do j = jmin, jmax
    do i = imin, imax

      ! Compute stresses in local (TTI; TOr) frame from global (FD) frame
      ! values. Use O(2) bi-linear interpolation to collocate 
      ! shear stresses (tyz; txz; txy) on to normal stress (txx; tyy; tzz)
      ! location:

      ! ..don't need any interpol for normal stresses:
      txx_ = txx(i,j,k)
      tyy_ = tyy(i,j,k)
      tzz_ = tzz(i,j,k)

      ! ..tyz (bi-linear interpol in i-plane):
      tyz_ = 0.25*(tyz(i,j,k) + tyz(i,j,k-1) + &
                 & tyz(i,j-1,k-1) + tyz(i,j-1,k))      

      ! ..txz (bi-linear interpol in j-plane):
      txz_ = 0.25*(txz(i,j,k) + txz(i+1,j,k) + &
                 & txz(i+1,j,k-1) + txz(i,j,k-1))      

      ! ..txy (bi-linear interpol in k-plane):
      txy_ = 0.25*(txy(i,j,k) + txy(i,j-1,k) + &
                 & txy(i+1,j-1,k) + txy(i+1,j,k))      


      ! Rotate stress tensor global (FD) to local (TTI; TOr) coordinates:
      azim_ = azim(i,j,k)
      dip_  = dip(i,j,k)
      rake_ = rake(i,j,k)

      cosa = cos(azim_)
      cosd = cos(dip_)
      cosr = cos(rake_)
      sina = sin(azim_)
      sind = sin(dip_)
      sinr = sin(rake_)

      ! ..Euler (ZYZ-order) rot matrix: global (FD) to local (TTI; TOr):
      b(1,1) = cosa*cosd*cosr - sina*sinr
      b(1,2) = sina*cosd*cosr + cosa*sinr
      b(1,3) = -sind*cosr
      b(2,1) = -cosa*cosd*sinr - sina*cosr
      b(2,2) = -sina*cosd*sinr + cosa*cosr
      b(2,3) = sind*sinr
      b(3,1) = cosa*sind
      b(3,2) = sina*sind
      b(3,3) = cosd

      ! ..bik x tkl product (global (FD) to local (TOr)):
      bt(1,1) = b(1,1)*txx_ + b(1,2)*txy_ + b(1,3)*txz_
      bt(1,2) = b(1,1)*txy_ + b(1,2)*tyy_ + b(1,3)*tyz_
      bt(1,3) = b(1,1)*txz_ + b(1,2)*tyz_ + b(1,3)*tzz_
      bt(2,1) = b(2,1)*txx_ + b(2,2)*txy_ + b(2,3)*txz_
      bt(2,2) = b(2,1)*txy_ + b(2,2)*tyy_ + b(2,3)*tyz_
      bt(2,3) = b(2,1)*txz_ + b(2,2)*tyz_ + b(2,3)*tzz_
      bt(3,1) = b(3,1)*txx_ + b(3,2)*txy_ + b(3,3)*txz_
      bt(3,2) = b(3,1)*txy_ + b(3,2)*tyy_ + b(3,3)*tyz_
      bt(3,3) = b(3,1)*txz_ + b(3,2)*tyz_ + b(3,3)*tzz_

      ! ..stresses in local (TTI; TOr)-aligned coordinate system:
      txx_loc(i,j,k) = bt(1,1)*b(1,1) + bt(1,2)*b(1,2) + bt(1,3)*b(1,3)
      tyy_loc(i,j,k) = bt(2,1)*b(2,1) + bt(2,2)*b(2,2) + bt(2,3)*b(2,3)
      tzz_loc(i,j,k) = bt(3,1)*b(3,1) + bt(3,2)*b(3,2) + bt(3,3)*b(3,3)
      txy_loc(i,j,k) = bt(1,1)*b(2,1) + bt(1,2)*b(2,2) + bt(1,3)*b(2,3)
      txz_loc(i,j,k) = bt(1,1)*b(3,1) + bt(1,2)*b(3,2) + bt(1,3)*b(3,3)
      tyz_loc(i,j,k) = bt(2,1)*b(3,1) + bt(2,2)*b(3,2) + bt(2,3)*b(3,3)

    enddo
  enddo
enddo
!$OMP END PARALLEL DO


return
end subroutine
