subroutine generateSincWeights(dx_frac, dy_frac, dz_frac, fsinc)

!  Computes 8x8x8 grid of weighting coefficients for use in 3D O(8) sinc
!  extrapolation of off-grid sources (fx, fy, fz, pr) and sinc interpolation
!  of off-grid receivers (vx, vy, vz, pr). Implementation assumes a 
!  staggered grid, and uses a Kaiser windowing function.
!
!
!  o Last Modified:  
!    10-03-12 Adapted from 2D to 3D. 
!
!  o Written by Kurt T. Nihei; 09-26-12


implicit none

real :: bessi0  ! function to compute I0 (modified Bessel funct of 0 order)

real :: dx_frac, dy_frac, dz_frac

integer :: ix, iy, iz
real :: x_cell, y_cell, z_cell, b_x, b_y, b_z, win_x, win_y, win_z
real :: fsinc_x, fsinc_y, fsinc_z
real, dimension(8,8,8) :: fsinc

real, parameter :: b = 4.14  ! optimal Kaiser window param for kmax = 2pi/3
real, parameter :: r = 4.    ! half-width of sinc interpolator
real :: pi
pi = acos(-1.)


! Compute sinc interpol weights:
do iz = 1, 8
  ! cells at which to sample sinc func [normalized]
  z_cell = float(iz) - r - dz_frac  

  ! compute Kaiser window:
  if (abs(z_cell) <= r) then
    b_z = b*sqrt(1 - (z_cell/r)**2)
  else
    b_z = 0.
  endif
  win_z = bessi0(b_z)/bessi0(b)

  ! compute sinc interpolation function:
  if (z_cell == 0) then
    fsinc_z = 1.  ! rec is located on a horiz grid line
  else
    fsinc_z = win_z*sin(z_cell*pi)/(z_cell*pi) 
  endif

  do iy = 1, 8
    ! cells at which to sample sinc func [normalized]
    y_cell = float(iy) - r - dy_frac  

    ! compute Kaiser window:
    if (abs(y_cell) <= r) then
      b_y = b*sqrt(1 - (y_cell/r)**2)
    else
      b_y = 0.
    endif
    win_y = bessi0(b_y)/bessi0(b)

    ! compute sinc interpolation function:
    if (y_cell == 0) then
      fsinc_y = 1.  ! rec is located on a vert grid line
    else
      fsinc_y = win_y*sin(y_cell*pi)/(y_cell*pi) 
    endif

    do ix = 1, 8
      ! cells at which to sample sinc func [normalized]
      x_cell = float(ix) - r - dx_frac  

      ! compute Kaiser window:
      if (abs(x_cell) <= r) then
        b_x = b*sqrt(1 - (x_cell/r)**2)
      else
        b_x = 0.
      endif
      win_x = bessi0(b_x)/bessi0(b)

      ! compute sinc interpolation function:
      if (x_cell == 0) then
        fsinc_x = 1.  ! rec is located on a vert grid line
      else
        fsinc_x = win_x*sin(x_cell*pi)/(x_cell*pi) 
      endif
      fsinc(ix,iy,iz) = fsinc_x*fsinc_y*fsinc_z

    enddo

  enddo

enddo

return
end subroutine
