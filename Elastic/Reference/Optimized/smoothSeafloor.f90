subroutine smoothseafloor(mtype, ilo, iup, klo, kup, nx, nz, dx, dz, &
                        & c11, c33, c55, c13, c15, c35, rho)
 
! Smooths seafloor using Muir-Dellinger EM averaging. Effective medium 
! elastic constants (density) are computed using area-fraction weighted 
! harmonic (geometric) averaging (Muir et al., 1992, Geophys., 59, 
! 1189-1193). Note VTI EM averaging is used for dipping seafloor as TTI
! was found to be numerically unstable.

! o Last Modified: 
!   01-10-12 Concerned about stability of the VTI cIJ's obtained from
!            Backus averaging --> outputting only isotropic BA cIJ's!
!   12-19-11 Added conditional to update appropriate cIJ's depending
!            on mtype.
!   12-19-11 Corrected error in seafloor cell location; correct 
!            location is at shifted z_raw(i) = dz*float(k).
!   12-19-11 Changed qc_seafloor.vo gocad files to Tecplot360 format
!            for Vpvert, Vsvert, rho, delta, epsilon, and changed
!            file name to qc_seafloor_props.dat.
!   12-19-11 Added Tecplot360 header to qc_seafloor.out file and 
!            changed name to qc_seafloor_profile.dat.
!   12-12-11 Extended smoother to the full-range x:[1,nx], i.e., to
!            include smoothing up the the left and right edges.
!   11-15-11 Added swap4bytes on QC binary files for Gocad viewing.
!   11-13-11 Added fit.f90 routine to generate the smoothed seafloor
!            profile (replaces gen_seabot_prof.m Matlab routine).

! o Adapted from: 
!   /data/mod/knih/Projects3/Wheatstone/2D/Wheat5/mod2d_wheat5.f90
!
 
implicit none 

integer :: ilo, iup, klo, kup
integer :: nx, nz
real :: dx, dz 
real, dimension(ilo:iup,klo:kup) :: c11, c33, c55, c13, c15, c35, rho

character*3 :: celltype
character (len=20) :: filename1
character (len=1) :: tab         ! tab=char(9)

integer :: i, k, iflag, ik, nsamp, k_qctop, k_qcbot, nz_qc
integer :: mtype, status
integer :: mwt, ifit, icen ! parameters for fit.f90
integer, parameter :: ndata = 15  ! #pts in smoother
integer, parameter :: ndata1 = 7  ! #pts on one side of the smoother

real :: xdum
real :: x, z, xx
real :: aa, bb, siga, sigb, chi2, qq  ! parameters for fit.f90
real, dimension(-ndata1+1:nx+ndata1) :: x_sea, z_raw, z_sea
real, allocatable, dimension (:) :: x_fit, z_fit, sig
real :: c11_bav, c33_bav, c55_bav, c13_bav
real :: vp, vs, f_sed, f_wat, dn_av
real :: x_lef, x_rig, z_lef, z_rig, z_top, z_bot, slope
real :: lam_sed, lam_wat, mu_sed, mu_wat
real :: vpvert, vsvert, del, eps, eta
 
!*************************************************************************!

tab = char(9)
 
! Initialize arrays to zero:
c11_bav = 0
c33_bav = 0
c55_bav = 0
c13_bav = 0

! Locate seafloor cells:
do i = 1, nx
  x_sea(i) = dx*float(i-1)
  do k = 1, nz
    if (rho(i,k) >= 1100.) z_raw(i) = dz*float(k) 
  enddo
enddo

! Define strip of thickness nz_qc for use in qc_seafloor_props.dat:
k_qctop = minval(abs(z_raw))/dz
k_qcbot = maxval(abs(z_raw))/dz + 10
if (k_qctop < 10) then
  k_qctop = 1             ! start strip from top of model
else
  k_qctop = k_qctop - 10  ! add 10 cells to top of strip
endif
nz_qc = k_qcbot - k_qctop + 1

! Pad edges with ndata1 points:
do i = -ndata1+1, 0     ! left edge
  x_sea(i) = x_sea(1)
  z_raw(i) = z_raw(1)
enddo
do i = nx+1, nx+ndata1  ! right edge
  x_sea(i) = x_sea(nx)
  z_raw(i) = z_raw(nx)
enddo

! Initialize z_sea to original, unsmoothed profile
z_sea = z_raw  

! Smooth by fitting line to ndata points:
mwt = 0
allocate(x_fit(ndata), z_fit(ndata), sig(ndata))
do i = 1, nx
  do ifit = 1, ndata 
    x_fit(ifit) = x_sea(i+ifit-1-ndata1)
    z_fit(ifit) = z_raw(i+ifit-1-ndata1)
  enddo
  call fit(x_fit, z_fit, ndata, sig, mwt, aa, bb, siga, sigb, chi2, qq)
  icen = i  ! center point of the linear regression
  z_sea(icen) = aa + bb*x_fit(1+ndata1)
enddo
deallocate(x_fit, z_fit)

! Write seafloor profile QC file for Tecplot360:
open(unit=1, file='qc_seafloor_profile.dat', status='unknown', &
   & action='write', iostat=status)
write(1,*) 'TITLE="qc_seafloor_profile"'
write(1,*) 'VARIABLES="x (m)","z (m)"'
write(1,*) 'ZONE T="before smoothing (raw)", F=POINT'
do i = 1, nx
  write(1,*) x_sea(i), z_raw(i)
enddo
write(1,*) 'ZONE T="after smoothing", F=POINT'
do i = 1, nx
  write(1,*) x_sea(i), z_sea(i)
enddo
close(unit=1)

! Write original seafloor properties to QC file for Tecplot360:
open(unit=1, file='qc_seafloor_props.dat', status='unknown', &
   & action='write', iostat=status)
write(1,*) 'TITLE="qc_seafloor_properties"'
write(1,*) 'VARIABLES="x (m)","z (m)", "Vpvert(m/s)", "Vsvert(m/s)", ' 
write(1,*) '"density(kg/m^3)", "delta", "epsilon", "eta"' 
write(1,*) 'ZONE T="before Backus ave", I=', nx,' J=', nz_qc, ', F=POINT'
do k = k_qctop, k_qcbot
    z = -dz*float(k-1)
  do i = 1, nx
    x = dx*float(i-1)
    vpvert = sqrt(c33(i,k)/rho(i,k))
    vsvert = sqrt(c55(i,k)/rho(i,k))
    del = ((c13(i,k)+c55(i,k))**2 - (c33(i,k)-c55(i,k))**2) / &
           & (2*c33(i,k)*(c33(i,k)-c55(i,k)))
    eps = (c11(i,k)-c33(i,k)) / (2*c33(i,k))
    eta = (eps - del)/(1. + 2.*del)

    write(1,*) x, z, vpvert, vsvert, rho(i,k), del, eps, eta
  enddo
enddo

! Replace seafloor nodes with Backus Averaged properties:
! NOTE!!: For this calculation, "bot" refers to closest to z = 0, and 
!         "top" to farthest from z = 0; funky logic, but this is what
!         was deemed logical at the time -- see 04-01-09 Muir-Dellinger
!         Seafloor Smoothing Notes for the "5 Cases."
do k=4, nz-2
  z_bot = float(k-1)*dz  ! "bottom" of cell is defined as closest to z=0
  z_top = float(k)*dz    ! "top" of cell is defined as farthest from z=0

  do i=1, nx
    x = float(i-1)*dx
    xx = float(i)*dx

    ! Flag any cells that contain the seafloor and compute points of 
    !  intersection of seafloor with cell walls:
    iflag = 0
    celltype = '000'

    ! seafloor intersects Left side of cell:
    if (z_sea(i) >= z_bot .and. z_sea(i) <= z_top) then   
      x_lef = x
      z_lef = z_sea(i)

      ! case Left -> Right:
      if (z_sea(i+1) >= z_bot .and. z_sea(i+1) <= z_top) then  
        iflag = 1
        celltype = 'l-r'
        x_rig = xx
        z_rig = z_sea(i+1)
        f_sed = 0.5*((z_top - z_lef) + (z_top - z_rig))*(x_rig - x_lef)/ &
              & (dx*dz)
        f_wat = 1 - f_sed
           if (f_sed > 1 .or. f_wat > 1) then
       !!!!  write(*,*) 'Left->Right: i, k =', i, k
       !!!!  write(*,*) 'f_sed =', f_sed
       !!!!  write(*,*) 'x_rig =', x_rig
       !!!!  write(*,*) 'z_rig =', z_rig
       !!!!  write(*,*) 'x_lef =', x_lef
       !!!!  write(*,*) 'z_lef =', z_lef
       !!!!  write(*,*) 'z_top =', z_top
           endif

      ! case Left -> Bottom:
      elseif (z_sea(i+1) <= z_bot .and. z_sea(i) > z_bot) then  
        iflag = 1
        celltype = 'l-b'
        z_rig = z_bot
        slope = (z_sea(i+1) - z_sea(i))/(xx - x)
        x_rig = (z_rig - z_lef)/slope + x
        f_wat = 0.5*(x_rig - x_lef)*(z_lef - z_rig)/(dx*dz)
        f_sed = 1 - f_wat
           if (f_sed > 1 .or. f_wat > 1) then
       !!!!  write(*,*) 'Left->Bot: i, k =', i, k
       !!!!  write(*,*) 'f_sed =', f_sed
       !!!!  write(*,*) 'x_rig =', x_rig
       !!!!  write(*,*) 'z_rig =', z_rig
       !!!!  write(*,*) 'x_lef =', x_lef
       !!!!  write(*,*) 'z_lef =', z_lef
       !!!!  write(*,*) 'z_bot =', z_bot
           endif

      ! case Left -> Top:
      elseif (z_sea(i+1) >= z_top .and. z_sea(i) < z_top) then  
        iflag = 1
        celltype = 'l-t'
        z_rig = z_top
        slope = (z_sea(i+1) - z_sea(i))/(xx - x)
        x_rig = (z_rig - z_lef)/slope + x
        f_sed = 0.5*(x_rig - x_lef)*(z_rig - z_lef)/(dx*dz)
        f_wat = 1 - f_sed
           if (f_sed > 1 .or. f_wat > 1) then
       !!!!  write(*,*) 'Left->Top: i, k =', i, k
       !!!!  write(*,*) 'f_sed =', f_sed
       !!!!  write(*,*) 'x_rig =', x_rig
       !!!!  write(*,*) 'z_rig =', z_rig
       !!!!  write(*,*) 'x_lef =', x_lef
       !!!!  write(*,*) 'z_lef =', z_lef
       !!!!  write(*,*) 'z_top =', z_top
           endif

      endif

    ! seafloor intersects Right side of cell:
    elseif (z_sea(i+1) >= z_bot .and. z_sea(i+1) <= z_top) then   
      x_rig = xx
      z_rig = z_sea(i+1)

      ! case Bottom -> Right:
      if (z_sea(i) <= z_bot .and. z_sea(i+1) > z_bot) then  
        iflag = 1
        celltype = 'b-r'
        z_lef = z_bot
        slope = (z_sea(i+1) - z_sea(i))/(xx - x)
        x_lef = (z_lef - z_sea(i))/slope + x
        f_wat = 0.5*(x_rig - x_lef)*(z_rig - z_lef)/(dx*dz)
        f_sed = 1 - f_wat
           if (f_sed > 1 .or. f_wat > 1) then
       !!!!  write(*,*) 'Bot->Right: i, k =', i, k
       !!!!  write(*,*) 'f_sed =', f_sed
       !!!!  write(*,*) 'x_rig =', x_rig
       !!!!  write(*,*) 'z_rig =', z_rig
       !!!!  write(*,*) 'x_lef =', x_lef
       !!!!  write(*,*) 'z_lef =', z_lef
       !!!!  write(*,*) 'z_bot =', z_bot
           endif

      ! case Top -> Right:
      elseif (z_sea(i) >= z_top .and. z_sea(i+1) < z_top) then  
        iflag = 1
        celltype = 't-r'
        z_lef = z_top
        slope = (z_sea(i+1) - z_sea(i))/(xx - x)
        x_lef = (z_lef - z_sea(i))/slope + x
        f_sed = 0.5*(x_rig - x_lef)*(z_lef - z_rig)/(dx*dz)
        f_wat = 1 - f_sed
           if (f_sed > 1 .or. f_wat > 1) then
       !!!!  write(*,*) 'Top->Right: i, k =', i, k
       !!!!  write(*,*) 'f_sed =', f_sed
       !!!!  write(*,*) 'x_rig =', x_rig
       !!!!  write(*,*) 'z_rig =', z_rig
       !!!!  write(*,*) 'x_lef =', x_lef
       !!!!  write(*,*) 'z_lef =', z_lef
       !!!!  write(*,*) 'z_top =', z_top
           endif

      endif

    endif

    ! Now do some cleaning up:
    if (iflag == 1) then
      if (rho(i,k+1) > 1100.) then  ! can't have water below seafloor
        rho(i,k+1) = rho(i,k+2)
        c33(i,k+1) = c33(i,k+2)
        c11(i,k+1) = c11(i,k+2)
        c55(i,k+1) = c55(i,k+2)
        c13(i,k+1) = c13(i,k+2)
      endif
      if (rho(i,k-1) < 1100.) then  ! can't have sediment above seafloor
        rho(i,k-1) = rho(i,k-2)
        c33(i,k-1) = c33(i,k-2)
        c11(i,k-1) = c11(i,k-2)
        c55(i,k-1) = c55(i,k-2)
        c13(i,k-1) = c13(i,k-2)
      endif
      if (rho(i,k-2) < 1100.) then  ! can't have sediment above seafloor
        rho(i,k-2) = rho(i,k-3)
        c33(i,k-2) = c33(i,k-3)
        c11(i,k-2) = c11(i,k-3)
        c55(i,k-2) = c55(i,k-3)
        c13(i,k-2) = c13(i,k-3)
      endif
      !  don't need EM seafloor smoothing:
      if (z_sea(i) == z_top .and. z_sea(i+1) == z_top) iflag = 0   
      if (z_sea(i) == z_bot .and. z_sea(i+1) == z_bot) iflag = 0   
    endif

    ! cell contains seafloor:
    if (iflag == 1) then  
      
      ! Carry-out Backus averaging:
      ! ..sediment properties:
      mu_sed = c55(i,k+2)
      lam_sed = c33(i,k+2) - 2*mu_sed

      ! ..seawater properties:
      mu_wat = 0  !!! setting mu_wat zero!
      lam_wat = c33(i,k-2) - 2*mu_wat 

      c11_bav = f_sed*4*mu_sed*(lam_sed + mu_sed)/c33(i,k+2) + &
                f_wat*4*mu_wat*(lam_wat + mu_wat)/c33(i,k-2) + &
              & 1/(f_sed/c33(i,k+2) + f_wat/c33(i,k-2)) * &
              & (f_sed*lam_sed/c33(i,k+2) + f_wat*lam_wat/c33(i,k-2))**2


      c33_bav = 1/(f_sed/c33(i,k+2) + f_wat/c33(i,k-2))

      !!!c55_bav = 1/(f_sed/mu_sed + f_wat/mu_wat)
      c55_bav = 0.   ! <--for mu_wat = 0

      c13_bav = (1/(f_sed/c33(i,k+2) + f_wat/c33(i,k-2)))* &
              & (f_sed*lam_sed/c33(i,k+2) + f_wat*lam_wat/c33(i,k-2))
      
      ! Averaged (arithmetic) density:
      dn_av = f_sed*rho(i,k+2) + f_wat*rho(i,k-2)

      ! Set seafloor to Backus averaged VTI values:
      rho(i,k) = dn_av
      c15(i,k) = 0.
      c35(i,k) = 0.
!!!!
!!!! NOT USING VTI BACKUS AVERAGING DUE TO STABILITY CONCERNS; 
!!!! (SEE COMMENT ON 01-10-12). 
!!!!  ! ..assign Backus averaged properties to appropriate material types:
!!!!  if (mtype == 3 .or. mtype == 4 .or. mtype == 13 .or. mtype == 14) then
!!!!    ! VTI, TTI, ViscVTI & ViscTTI:
!!!!    c11(i,k) = c11_bav 
!!!!    c33(i,k) = c33_bav 
!!!!    c55(i,k) = c55_bav 
!!!!    c13(i,k) = c13_bav 
!!!!  elseif (mtype == 2 .or. mtype == 12) then
!!!!    ! Iso & ViscIso:
!!!!    c33(i,k) = c33_bav 
!!!!    c55(i,k) = c55_bav 
!!!!    c11(i,k) = c33(i,k)
!!!!    c13(i,k) = c33(i,k) - 2.*c55(i,k)
!!!!  elseif (mtype == 1 .or. mtype == 11) then
!!!!    ! Aco & ViscAco:
!!!!    c33(i,k) = c33_bav 
!!!!    c11(i,k) = c33(i,k)
!!!!    c55(i,k) = 0.
!!!!    c13(i,k) = c33(i,k)
!!!!  endif
!!!!
      ! ..assign Backus averaged properties to appropriate material types:
      if (mtype == 1 .or. mtype == 11) then
        ! Aco & ViscAco:
        c33(i,k) = c33_bav 
        c11(i,k) = c33(i,k)
        c55(i,k) = 0.
        c13(i,k) = c33(i,k)
      else
        ! Iso & ViscIso:
        c33(i,k) = c33_bav 
        c55(i,k) = c55_bav 
        c11(i,k) = c33(i,k)
        c13(i,k) = c33(i,k) - 2.*c55(i,k)
      endif
    endif
  enddo
enddo


! Write Backus averaged seafloor properties to QC file for Tecplot360:
write(1,*) 'ZONE T="after Backus ave", I=', nx,' J=', nz_qc, ', F=POINT'
do k = k_qctop, k_qcbot
    z = -dz*float(k-1)
  do i = 1, nx
    x = dx*float(i-1)
    vpvert = sqrt(c33(i,k)/rho(i,k))
    vsvert = sqrt(c55(i,k)/rho(i,k))
    del = ((c13(i,k)+c55(i,k))**2 - (c33(i,k)-c55(i,k))**2) / &
           & (2*c33(i,k)*(c33(i,k)-c55(i,k)))
    eps = (c11(i,k)-c33(i,k)) / (2*c33(i,k))
    eta = (eps - del)/(1. + 2.*del)

    write(1,*) x, z, vpvert, vsvert, rho(i,k), del, eps, eta
  enddo
enddo
close(unit=1)


return
end subroutine
