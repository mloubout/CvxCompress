program genFracVTI
 
! /data/mod/knih/Tests2/FDElaQ_v1/genFracVTI.f90
!
! Creates flat binary model files for fractured VTI solid. Properties 
! taken from Tsvankin (1997, Geophys., 62(4), 1292-1309, Fig. 5).
!
! Last Modified: 
! 01-22-13  Created. 
!
 
implicit none 

integer, parameter :: nx = 200   
integer, parameter :: ny = 200   
integer, parameter :: nz = 200   

real, parameter :: dx = 12.5  ! grid spacings in m
real, parameter :: dy = 12.5  
real, parameter :: dz = 12.5  

integer :: i, j, k
integer :: status

real :: pi, deg2rad

real, dimension(nx,ny,nz) :: rho, &
& vp0, vs0, eps1, eps2, del1, del2, del3, gam1, gam2, azim, dip, rake, q
 
! **************************************************************************
pi = acos(-1.)
deg2rad = pi/180.

! Generate flat binary model files. 

write(*,*) 'Generating flat bin model files for fractured VTI ..'

do i = 1, nx
  do j = 1, ny
    do k = 1, nz
      rho(i,j,k)  = 2200.
      vp0(i,j,k)  = 2437.
      vs0(i,j,k)  = 1265.
      eps1(i,j,k) = 0.329
      eps2(i,j,k) = 0.258
      del1(i,j,k) =  0.083
      del2(i,j,k) = -0.078
      del3(i,j,k) = -0.106
      gam1(i,j,k) =  0.182
      gam2(i,j,k) =  0.0455
      azim(i,j,k) =  45.*deg2rad
      dip(i,j,k)  =  20.*deg2rad
      rake(i,j,k) =  30.*deg2rad
      q(i,j,k) =  80.
    enddo
  enddo
enddo

open(unit=2, file='rho.bin', status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nx*ny*nz)
write(2) (((rho(i,j,k), k=1,nz), j=1,ny), i=1,nx)    
close(unit=2)

open(unit=2, file='vp0.bin', status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nx*ny*nz)
write(2) (((vp0(i,j,k), k=1,nz), j=1,ny), i=1,nx)    
close(unit=2)

open(unit=2, file='vs0.bin', status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nx*ny*nz)
write(2) (((vs0(i,j,k), k=1,nz), j=1,ny), i=1,nx)    
close(unit=2)

open(unit=2, file='eps1.bin', status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nx*ny*nz)
write(2) (((eps1(i,j,k), k=1,nz), j=1,ny), i=1,nx)    
close(unit=2)

open(unit=2, file='eps2.bin', status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nx*ny*nz)
write(2) (((eps2(i,j,k), k=1,nz), j=1,ny), i=1,nx)    
close(unit=2)

open(unit=2, file='del1.bin', status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nx*ny*nz)
write(2) (((del1(i,j,k), k=1,nz), j=1,ny), i=1,nx)    
close(unit=2)

open(unit=2, file='del2.bin', status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nx*ny*nz)
write(2) (((del2(i,j,k), k=1,nz), j=1,ny), i=1,nx)    
close(unit=2)

open(unit=2, file='del3.bin', status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nx*ny*nz)
write(2) (((del3(i,j,k), k=1,nz), j=1,ny), i=1,nx)    
close(unit=2)

open(unit=2, file='gam1.bin', status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nx*ny*nz)
write(2) (((gam1(i,j,k), k=1,nz), j=1,ny), i=1,nx)    
close(unit=2)

open(unit=2, file='gam2.bin', status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nx*ny*nz)
write(2) (((gam2(i,j,k), k=1,nz), j=1,ny), i=1,nx)    
close(unit=2)

open(unit=2, file='azim.bin', status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nx*ny*nz)
write(2) (((azim(i,j,k), k=1,nz), j=1,ny), i=1,nx)    
close(unit=2)

open(unit=2, file='dip.bin', status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nx*ny*nz)
write(2) (((dip(i,j,k), k=1,nz), j=1,ny), i=1,nx)    
close(unit=2)

open(unit=2, file='rake.bin', status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nx*ny*nz)
write(2) (((rake(i,j,k), k=1,nz), j=1,ny), i=1,nx)    
close(unit=2)

open(unit=2, file='q.bin', status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nx*ny*nz)
write(2) (((q(i,j,k), k=1,nz), j=1,ny), i=1,nx)    
close(unit=2)

! **************************************************************************

write(*,*) 'End of program.'

stop 
end program
