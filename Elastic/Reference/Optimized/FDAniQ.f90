program FDAniQ

!--------------------------------------------------------------------------
!    FDAniQ = Finite Difference Anisotropic Modeling with Q.
!               Visco Elastic:  eIso; eVTI; eVOr; eTTI; eTOr  
!               Visco Acoustic: aTTI; aTOr
!
!             3D finite difference time domain (FDTD) code for modeling
!             elastic wave propagation in a viscoelastic medium with
!             the following types of anisotropy: 
!
!               (1) elastic isotropic (eIso; 2 cIJ's); 
!               (2) elastic vertical transverse isotropy (eVTI; 5 cIJ's); 
!               (3) elastic FD grid-aligned orthorhombic (eVOr; 9 cIJ's); 
!               (4) elastic tilted transverse isotropy (eTTI; 5 cIJ's; 
!                   2 rot angles); 
!               (5) elastic tilted orthorhombic (eTOr; 9 cIJ's; 3 
!                   rot angles); 
!               (6) acoustic tilted transverse isotropy (aTTI; 3 cIJ's; 
!                   2 rot angles); and 
!               (7) acoustic tilted orthorhombic (aTOr; 6 cIJ's; 
!                   3 rot angles).
!
!             The first-order stress-particle velocity system of 
!             PDE's (i.e.,equation of motion and Hooke's Law) are 
!             solved on a staggered grid (Yee/Vireaux) with 
!             finite differences using a leap-frog scheme. The time
!             derivatives are computed using a symplectic time integration
!             scheme with either O(2) or O(4) accuracy (Omelyan 2003; 
!             Comp. Phys. Comm., 146, 188-202). The space derivatives 
!             are computed with O(16) acuracy Taylor polynomials.
!
!             Tilted anisotropy (eTTI; eTOr; aTTI; aTOr) is modeled by
!             rotating the cIJ tensor from its local (TTI or TOr) frame
!             to the global (FD) frame via a Bond transformation with 
!             3 Euler rotations (rake -> dip -> azimuth). Stresses
!             are updated in the global (FD) frame. "Off-kite" strain
!             rates are computed using bi-linear interpolation of the
!             particle velocities, follwed by O(8) FD. "On-kite" strain
!             rates do not require interpolation and are computed with
!             O(16) FD. The particle velocity update does not require
!             interpolation of the stresses, and are computed with O(16)
!             FD.
!
!             For the acoustic titled anisotropy cases (aTTI; aTOr) the
!             full 21 cIJ tensor is used, but the values of the shear
!             moduli (c44; c55; c66) are set to a fraction of c11, c22,
!             and c33. This approach requires the same # of computations
!             as the fully elastic case, so it is included here mainly
!             for use as a QC tool for acoustic tilted anisotropy 
!             modeling. Note that this approach has no restrictions on
!             eta >= 0, as is the case with some modeling approaches
!             that are based on dispersion relations.
!
!             Attenuation is approximately modeled with a single set
!             of SLS memory variables that are applied to the particle 
!             velocities. The particle velocity memory variable 
!             formulation reduces the number of memory variable from
!             six, required by the standard stress memory variable
!             formulation, to 3. The approximate particle velocity
!             memory variable formulation assumes that Q^-1 is applied 
!             equally to  all six stresses, and that the spatial 
!             variations in Q are smooth.
!
!             This code is for a single shot composed of multiple
!             sources (i.e., no MPI shot distribution). OpenMP 
!             threading is performed on the stress update and particle
!             velocity update outer space loops.
!
!             source code loc.: file:///devl/FD/Svn/FDAniQ/trunk
!
! o Written by: Kurt Nihei
!
! o To Do List: 
!   01-09-13  Add mIJ source capability.
!
! o Last Modified (FDAniQ):
!   06-04-12  Version 1.0: Adapted from SSG2D 
!             (file:///devl/FD/Svn/Ssg2d/trunk).
!
! o Program Notes:
!
!   (1) Anisotropic stress-velocity formulation for a 3D viscoelastic
!       medium with orthotropic anisotropy using a minimal (Yee/Vireaux) 
!       staggered grid stencil:
!
!       z
!       ^          front face of unit cell:       back face of unit cell:
!       |
!       +---> x      vx >-----o txx,tyy,tzz        txy #-----* vy
!       y               |     |                        |     |
!                       |     |                        |     |
!                   txz #-----^ vz                     +-----# tyz
!       
!       Note that no cIJ averaging is required because the all six stresses
!       are computed at a single grid location (the txx,tyy,tzz location)
!       where the cIJ's are defined.  Averaging (arithmetic) of the
!       particle velocities (also defined at the txx,tyy,tzz location) is 
!       performed to center the density to the proper location required
!       for the particle velocity updates.
!
!   (2) Absorbing boundaries are achieved by ramping down the particle 
!       velocities and stresses in strip surrounding the model, and 
!       applying Maxwell viscoelastic damping, which has a similar, but
!       but much more parsimonious, form to the PML absorbing boundary.
!
!   (3) Time domain attenuation is performed using an approximate 
!       reformulation of the standard linear solid (SLS) memory 
!       variable approach from stress variables to particle velocity
!       variables to reduce the number of memory variables from 6 to 3. 
!       This approximation assumes isotropic Q^-1 with smooth spatial
!       variations.
!
!   (4) Free-surface boundary conditions are implemented using the image 
!       method (Robertsson, Geophys., 61(6), 1921-1934, 1996) with a
!       modification to the particle velocity mirror conditions. The 
!       free-surface is assumed to be flat.  The free-surface boundary 
!       condition is implemented to O(16) spatial accuracy.
!      
!
!  3-D Computational Domain with Coordinate System (note: z-axis is up, 
!  y-axis is into page, and model origin is at the top, SW corner):
!
!     (z)
!      ^                  nabc_top
!      |:::::::::::::::::::::^::::
!      |:::::::::::::::::::::|::::                  
!  j(y)+--->i(x)-------------v---+          nabc_sdy          
!      |:::               ^   :::|    +--------------------------+ 
!      v:::               |   :::|    |::::::::^.:::::::^::::::::|
!      k:::               |   :::|    |::::::::v::::::: |:::: :::|
!      |:::  Model Domain:|   :::|    |:::              |     :::|
!      |:::  [nx,ny,nz]   |   :::|    |:::             ny     :::|
!      |:::               nz  :::|    |<-------nx-------|------->|
!      |:::               |   :::|    |:::              |     :::|
!      |<::---------nx--------::>|  j(y)::              |   nabc_sdx
!      |:::               |   :::|    ^:::              |     <--> 
!      |<:> nabc_sdx      |   :::|    |:::::::::::::::: |::::::::|
!      |:::               v   :::|    |:::::::::::::::::v::::::::|
!      +-------------------------+   k+--->i(x)------------------+ 
!      ::::::::::::::::::::::^:::::                               
!      ::::::::::::::::::::::v:::::
!        	          nabc_bot
!
! o Model Domain [nx,ny,nz]: includes side absorbing boundaries, but
!                            not the top and bottom absorbing boundaries.
!
! o FD Domain [nx,ny,nabc_top+nz+nabc_bot]: model domain plus top and 
!                            bottom absorbing boundaries.
!
! o Input Files:
!   FDAniQ.in - FD parameter file (ascii); contains grid, time step, 
!               rec type/locations, sou type/locations, sou wavelet.
!   propFile  - model property files (flat binary); density, 
!               Vp0, Vs0, delta1, delta2, delta3, epsilon1,
!               epsilon2, Q, azimuth, dip, rake
!
!-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-..-.-.-.-.-.-.-.-.-.-.-.-. */


!--DEFINE VARIABLES & ARRAYS-----------------------------------------------

implicit none

! Parameters:
real, parameter :: ver = 1.0 ! current version of FDAniQ.f90

integer, parameter :: nabc_top_d = 100  ! default #grd in Maxwell ABC - top
integer, parameter :: nabc_bot_d = 60   ! "  " - bottom 
integer, parameter :: nabc_sdx_d = 60   ! "  " - xsides 
integer, parameter :: nabc_sdy_d = 60   ! "  " - ysides 

integer, parameter :: orderspace = 16   ! O(16) FD spatial accuracy

real, parameter :: courant_safe = 0.95  ! Courant# safety factor (< 1)
real, parameter :: xisym_o2  = 0.5  ! symplectic O(2)
real, parameter :: lamsym_o2 = 1.0 
real, parameter :: xisym_o4 =   0.17861789584480907522  ! symplectic O(4) 
                        ! (Omelyan et al., 2002; with T. Jin corrections)
real, parameter :: lamsym_o4 = -0.21234183106260558041  
real, parameter :: chisym_o4 = -0.06626458266981843132 

real, parameter :: c0 =  41409225./33554432. ! O(16) SG FD constants
real, parameter :: c1 = -3578575./33554432. 
real, parameter :: c2 =  3864861./167772160. 
real, parameter :: c3 = -1254825./234881024.
real, parameter :: c4 =  325325./301989888.
real, parameter :: c5 = -61425./369098752.
real, parameter :: c6 =  7425./436207616.
real, parameter :: c7 = -143./167772160.

real, parameter :: c44overc33 = 1./100. ! ratio c44/c33 for aTOr
real, parameter :: c55overc33 = 1./100. ! ratio c55/c33 for aTTI & aTOr
real, parameter :: c66overc33 = 1./100. ! ratio c66/c33 for aTTI & aTOr

real, parameter :: qbig = 1.e6


! Integers for Loops & Constants:
integer :: i, j, k, ijk, ir, kti, l, isnps, lr, llr, lro, li
integer :: ii, jj, kk
integer :: ijk_dec
integer :: sti, nsti
integer :: idum
integer :: istat   ! file read status
real :: pi         ! pi

! OpenMP Variables:
integer :: nthreads

! Grid Properties:
integer :: nabc_top, nabc_bot, nabc_sdx, nabc_sdy  ! #grd in Maxwell ABC
integer :: ngho_buf  ! ghost buffer thickness in #cells
integer :: nx, ny, nz           ! #cells in model
integer :: nx_shot, ny_shot, nz_shot  ! #cells in FD model
integer :: ilo, jlo, klo        ! min dim of FD model (incl ABCs & FD op)
integer :: iup, jup, kup        ! max dim of FD model (incl ABCs & FD op)
integer :: imin, jmin  ! min dim of FD model (minus FD op)
integer :: imax, jmax  ! max dim of FD model (minus FD op)
real :: deta_maxtop  ! max Maxwell ABC damp coeff - top
real :: deta_maxbot  ! max Maxwell ABC damp coeff - bottom 
real, allocatable, dimension(:) :: deta_maxsdx, deta_maxsdy  ! "  " - side
real :: deta_maxsdx_top  ! max Maxwell ABC damp coeff - xside-top
real :: deta_maxsdy_top  ! max Maxwell ABC damp coeff - yside-top
real :: deta_maxsdx_bot  ! max Maxwell ABC damp coeff - xside-bot
real :: deta_maxsdy_bot  ! max Maxwell ABC damp coeff - yside-bot
real :: deta_gradx   ! max Maxwell ABC damp coeff - lin gradient top to bot
real :: deta_grady   ! max Maxwell ABC damp coeff - lin gradient top to bot
real :: vpvert_avtop, vpvert_avbot  ! used to calc max Maxwell ABC damp 
real :: dx, dy, dz              ! grid spacing in x, y, z
real :: x, y, z                 ! grid locations in x, y, z
real :: dum1, dum2 
real, allocatable, dimension(:,:,:) :: deta  ! Maxwell damping coeff ABC
logical :: abcoverrideflag  ! override default ABC values flag 
logical :: freesurfflag    ! free-surface flag for top of model  

double precision :: O(3), U(3), V(3), W(3)
integer :: NU, NV, NW
double precision :: tmp1,tmp2,tmp3

! Time Properties:
integer :: nti          ! #time steps
real :: dti             ! time step in sec.
real :: ti
real :: dl_min       ! =min(dx, dz)
real :: vp_max_c11   ! =maxval(sqrt(c11/rho))
real :: vp_max_c22   ! =maxval(sqrt(c22/rho))
real :: vp_max_c33   ! =maxval(sqrt(c33/rho))
real :: vp_max       ! =max(vp_max_c11,vp_m vz kcell 
real :: dti_override  ! =hardwired dti 
real :: dti_max      ! =courant*dl_min/vp_max
real :: courant      ! Courant#
logical :: dti_overrideflag   ! override dti flag

! Symplectic Time Integration Parameters:
integer :: ordertime  ! order of time integration
real, dimension(5) :: dti_tsym
real, dimension(4) :: dti_vsym
real :: dti_sym

! Receiver Properties:
integer :: imon, jmon, kmon  ! monitor rec locs (for monitor_S1234.dat)
integer :: ndti_mon    ! #time steps between monitor_S1234.dat writes
integer :: irec_offset ! receiver number offset -- for segywrite2d.c
integer :: krec_max  ! max rec depth; used to determine ghost buffer
integer :: nrec      ! #recs
integer :: iline  ! rec iline#
integer :: xline  ! rec xline#
real :: xmon, ymon, zmon  ! monitor rec locs (for monitor_S1234.dat)
real(8) :: t0, t1, t2                               ! clock runtime info 
real(8) :: dclock                  
logical :: recglobflag                   
real, allocatable, dimension(:) :: xrec, yrec, zrec  ! rec loc's
character (len=3) :: rectype ! rec type 
logical :: ghostrecflag   ! rec ghost
external dclock

! Source Properties:
integer :: nwav_raw     ! # points in wavelet - before resampling 
integer :: nwav         ! # points in wavelet - after resampling
integer :: nsou   ! #sources 
integer :: nsouline  ! #lines in sources section
integer :: ksou_max  ! max sou depth; used to determine ghost buffer
logical :: souglobflag           
                                             ! (=p;f;v;m) 
integer :: shotnum  ! shot#
real :: dti_raw   ! time increment in source wavelet before resamp
real :: fmax      ! max freq to be preserved in sou wavelet after resamp
                  ! (-3dB corner frequency in 5-pole Butterworth low pass) 
real :: fmax_h2o  ! max freq required to minimize dispersion in water 
real :: f40db     ! freq at -40dB in 5-pole Butterworth
real :: dl_max    ! max grid size (used in water dispersion calc)
real :: vmin      ! minimum velocity to be honored in overrided fmax calc
real :: xsou_, ysou_, zsou_ ! used in SEG-Y write
real :: dum
real, allocatable, dimension(:,:) :: souamp  ! sou amplitude
real, allocatable, dimension(:) :: swav_raw  ! sou wavelet before resamp
real, allocatable, dimension(:) :: swav, swav_  ! sou wavelet & time-integ.
                                                !  sou wavelet after resamp
real, allocatable, dimension(:) :: xsou, ysou, zsou  ! sou locs
character (len=3), allocatable, dimension(:) :: soutype ! sou type
character (len=3) :: srctype ! source type holder used when reading from file
character :: cdum
logical :: ghostsouflag   ! sou ghost
logical :: fmaxoverrideflag   ! overrides fmax if needed to min disp in H2O
integer :: ierr

! Material Properties:
integer, allocatable, dimension(:) :: ismin, ismax  ! used in fq est
real :: fq, wq            ! var's used for temp. storage of SLS parameters
real :: dt_ptp                       ! used in fq estimation
real :: tt, te      ! SLS stress and strain relax times
real :: dxvx, dyvy, dzvz, c13_  ! vars used in free-surface calc.
real, allocatable, dimension(:,:,:) :: rho  ! density (kg/m^3)
                   ! moduli (dispersed moduli for viscoelastic material)
real, allocatable, dimension(:,:,:) :: c11, c22, c33, c44, c55, c66, &
                                     & c12, c13, c14, c15, c16, &
                                     & c23, c24, c25, c26, &
                                     & c34, c35, c36, &
                                     & c45, c46, &
                                     & c56
real, allocatable, dimension(:,:,:) :: dip, azim, rake
real, allocatable, dimension(:,:,:) :: itausig, difitau  ! SLS vars 
real, allocatable, dimension(:,:,:) :: ptemp1 ! temp mat prop array
real, allocatable, dimension(:,:,:) :: ptemp2 ! temp mat prop array
character (len=3) :: orderaxes  ! fast-med-slow axes: xyz, zxy, zyx
character (len=4) :: anitype  ! eIso, eVTI, eVOr, eTTI, eTOr, aTTI, aTOr 
logical :: smoothflag  ! smooth seafloor

! Field Variables:
real :: p  ! pressure (temp var.)
real, allocatable, dimension(:,:,:) :: vx, vy, vz  ! particle velocities 
real, allocatable, dimension(:,:,:) :: txx, tyy, tzz  ! normal stresses 
real, allocatable, dimension(:,:,:) :: txz, tyz, txy  ! shear stresses 
real, allocatable, dimension(:,:,:) :: sx, sy, sz  ! memory variables
real, allocatable, dimension(:,:) :: txxtop, tyytop  ! stress vars for 
                                                     !  free-surface

! I/O Variables & Parameters:
integer :: icell, jcell, kcell  ! used in rec interpolation
integer :: icellmin, jcellmin, kcellmin  ! used in rec interpolation
integer :: icellmax, jcellmax, kcellmax  ! used in rec interpolation
integer :: nsnps  ! #snapshots
integer :: ntshift  ! #time samples to shift t=0 ref (for SEG-Ys)
integer :: kti_shift ! used in SEG-Y trace storage
integer :: nx_dec, ny_dec, nz_dec  ! snapshot dims in i,j,k
integer :: nsamp    ! #samples in vector for swap4bytes conversion
integer, allocatable, dimension(:) :: itsnp  ! snapshot timestep#'s
real :: tshift  ! time shift in sec. to shift t = 0 ref (for SEG-Ys)
real, allocatable, dimension(:) :: tsnp  ! snapshot times [s]
real, allocatable, dimension(:) :: tempswap  ! for byte swapping
real, allocatable, dimension(:,:,:) :: ptemp3  ! for pr snap write
integer :: status
character (len=5) :: filenum2  ! for snapshot files
character (len=6) :: filenum3  ! for shot#
character (len=7) :: filenum4  ! for cross plots
character (len=16) :: filename1  ! name for SEG-Y files
character (len=22) :: filename2  ! name for snap files
character (len=21) :: filename3  ! name for snap gocad flat bin header file
character (len=21) :: filename4  ! name for snap gocad flat binary file
character (len=19) :: filename6  ! name for monitor_S1234.dat file
character :: vp0File*128, vs0File*128, densityFile*128, delta1File*128, & 
& delta2File*128, delta3File*128, epsilon1File*128, epsilon2File*128, &
& gamma1File*128, gamma2File*128, qFile*128, dipFile*128, &
& azimuthFile*128, rakeFile*128, propFile*128
            !Vp3, Vs0, Density, Delta, Epsilon, Dip, Azimuth, Rake, 
            !Q file names
character (len=3) :: recnum      ! used in assigning output file #'s
character (len=1) :: tab         ! tab=char(9)
character (len=1) :: snpflag  ! snap file format (=t tecplot360;=g gocad)
logical :: bigendianflag  ! swap4bytes 

! SEG-Y Stuff:
integer :: nti_segy   ! # of time steps in SEG-Y file
integer :: kti_segy   ! current SEG-Y time step 
integer :: ktir       ! loop index used in forming traces vector
integer :: it_left    ! time sample used in trace interpol
real :: tmax          ! SEG-Y trace length [s]
real :: dti_segy      ! SEG-Y time increment [s]
real :: t_segy        ! SEG-Y time [s]
real :: dt_frac       ! fractional time step used in trace interp
real, allocatable, dimension(:) :: traces  ! vector for sending traces to 
                                           !  segywrite2d.c function
real, allocatable, dimension(:,:) :: prec 
real, allocatable, dimension(:,:) :: vxrec
real, allocatable, dimension(:,:) :: vyrec
real, allocatable, dimension(:,:) :: vzrec

! Sinc Interpolation/Extrapolation:
integer :: irectype
real :: dx_frac, dy_frac, dz_frac  ! frac distances for sinc interp/extrap
real :: f_in, f_out   ! used in sinc interp/extrap
real :: data_out
real, dimension(4) :: data4_in  ! used in trace interp
real, dimension(8,8,8) :: fsinc ! 8x8x8 sinc interp/extrap function

! Debug Stuff:
integer, dimension(3) :: loc(3)

integer :: num_args, ix, iy, iz
character(len=100), dimension(:), allocatable :: args

INTEGER, PARAMETER :: SEEK_SET = 0, SEEK_CUR = 1, SEEK_END = 2
integer :: filepos, souidx
integer :: nrecline, recidx

num_args = command_argument_count()
if (num_args > 0) then
   allocate(args(num_args))
   do ix = 1, num_args
      call get_command_argument(ix, args(ix))
      write(*,*) args(ix)
   end do
else
   write(*,*) 'Error! Need name of parameter input file.'
   stop 9
endif

! Define Some Tidbits:
pi = acos(-1.)   
tab = char(9)

!--READ FD PARAMETERS------------------------------------------------------

open(unit=1, file=args(1), status='old', action='read', &
   & iostat=status)
! ..read grid info.:	
write(*,*) '..reading headers (lines 1-3)'
read(1,*)     ! first 3 lines are headers
read(1,*)
read(1,*)
write(*,*) '..reading global origin'
read(1,*) (O(i), i=1,3)
write(*,*) '..reading U vector'
read(1,*)
read(1,*) (U(i), i=1,3)
write(*,*) '..reading V vector'
read(1,*)
read(1,*) (V(i), i=1,3)
write(*,*) '..reading W vector'
read(1,*)
read(1,*) (W(i), i=1,3)
write(*,*) '..read NU, NV, NW'
read(1,*)
read(1,*) NU, NV, NW
write(*,*) '..reading orderaxes (=zyx,xyz,zxy)'
read(1,*)
read(1,*) orderaxes
if (orderaxes == 'xyz') then
   nx = NU
   ny = NV
   nz = NW
   dx = sqrt(dot_product(U,U)) / (nx-1)
   dy = sqrt(dot_product(V,V)) / (ny-1)
   dz = sqrt(dot_product(W,W)) / (nz-1)
elseif (orderaxes == 'zxy') then
   nx = NV
   ny = NW
   nz = NU
   dx = sqrt(dot_product(V,V)) / (nx-1)
   dy = sqrt(dot_product(W,W)) / (ny-1)
   dz = sqrt(dot_product(U,U)) / (nz-1)
elseif (orderaxes == 'zyx') then
   nx = NW
   ny = NV
   nz = NU
   dx = sqrt(dot_product(W,W)) / (nx-1)
   dy = sqrt(dot_product(V,V)) / (ny-1)
   dz = sqrt(dot_product(U,U)) / (nz-1)
endif
write(*, '(A,F6.2,A,F6.2,A,F6.2)') 'dx=',dx,', dy=',dy,', dz=',dz

write(*,*) '..reading freesurfflag'
read(1,*)
read(1,*) freesurfflag  ! free-surface flag at top of model

ghostrecflag = .false.
ghostsouflag = .false.
if (.not. freesurfflag) then
  write(*,*) '..reading ghostrecflag & ghostsouflag'
  read(1,*)
  read(1,*) ghostrecflag, ghostsouflag  ! sou & rec ghost flags 
endif

! ..read time info:
write(*,*) '..reading tmax, dti_segy'
read(1,*)
read(1,*) tmax, dti_segy  ! SEG-Y trace length [s]; time increment [s]

nsnps = 0
ijk_dec = 1
snpflag = 'g'  ! default is gocad format for snapshots
write(*,*) '..reading nsnps, ijk_dec, snpflag'
read(1,*)
read(1,*) nsnps, ijk_dec, snpflag  ! #snapshots; decimat fac; 
                            ! file format (=t tecplot360; g=gocad)

if (nsnps > 0) then
  allocate (tsnp(nsnps), itsnp(nsnps))
  write(*,*) '..reading tsnp(k)'
  read(1,*)
  read(1,*) (tsnp(k), k=1,nsnps)
else
  write(*,*) '..skipping tsnp(k)'
  read(1,*)
  read(1,*)
endif

! ..read receiver info.: 
ntshift = 0
write(*,*) '..reading tshift'
read(1,*)
read(1,*) tshift  ! time [s] to shift the t=0 axis

write(*,*) '..reading nrec, rectype, iline, xline'
read(1,*)
read(1,*) nrecline, rectype, iline, xline

if (rectype == 'plz') then
  nrec = nrecline * nx * ny
else
  nrec = nrecline
endif

allocate(xrec(nrec), yrec(nrec), zrec(nrec))

write(*,*) '..reading recglobflag, xrec(li), yrec(li), zrec(li)'
read(1,*)
recidx = 1
do i=1,nrecline
  read(1,*) recglobflag, tmp1, tmp2, tmp3
  if (recglobflag) call coordGlobalToLocal(O,U,V,W,orderaxes,nx,ny,nz,&
      & dx,dy,dz,tmp1,tmp2,tmp3)
  if (rectype == 'plz') then
    do iy=1,ny
      do ix=1,nx
        xrec(recidx) = (ix-1)*dx
        yrec(recidx) = (iy-1)*dy
        zrec(recidx) = tmp3
        recidx = recidx + 1
      enddo
    enddo
  elseif (rectype == 'p') then
    xrec(i) = tmp1
    yrec(i) = tmp2
    zrec(i) = tmp3
    write(*,*) xrec(i),yrec(i),zrec(i)
  endif
enddo
krec_max = maxval(zrec/dz + 1)  ! store max rec depth for ghost region

! ..read monitor rec info. (for monitor.dat output file): 
write(*,*) '..reading ndti_mon, xmon, ymon, zmon'
read(1,*)
read(1,*) ndti_mon, xmon, ymon, zmon  ! #time steps between writes, 
                                      ! monitor rec (x,y,z) loc
imon = xmon/dx + 1
jmon = ymon/dy + 1
kmon = zmon/dz + 1

! ..read source info.: 
write(*,*) '..reading nsou'
read(1,*)
! #sources (comprising the shot), shot#
read(1,*) nsou, shotnum 

allocate (xsou(nsou), ysou(nsou), zsou(nsou)) 
allocate (soutype(nsou), souamp(nsou,6))

write(*,*) '..reading souglobflag, xsou, ysou, zsou, soutype, souamp'
write(*,*) '....soutype can be either p, m, v, f or plz'
write(*,*) '....soutype == plz creates a horizontal plane wave source'
read(1,*)
do k=1,nsou
  ! sou loc., sou type [=p,f,v,m], sou amp, sou wavelet#  
  read(1,*) souglobflag, tmp1, tmp2, tmp3, srctype 
  if (souglobflag) call coordGlobalToLocal(O,U,V,W,orderaxes,&
        & nx,ny,nz,dx,dy,dz,tmp1,tmp2,tmp3)
  xsou(k) = tmp1
  ysou(k) = tmp2
  zsou(k) = tmp3
  backspace(1)    ! now go back to read in source amplitude info.
  if (srctype == 'p') then
    read(1,*)  souglobflag, dum, dum, dum, soutype(k), souamp(k,1)  ! pressure sou 
    write(*,*) soutype(k),xsou(k),ysou(k),zsou(k),souamp(k,1)
  elseif (srctype == 'plz') then
    read(1,*)  souglobflag, dum, dum, dum, soutype(k)
    souamp(k,1) = 1.0
  elseif (srctype == 'm') then
    read(1,*)  souglobflag, dum, dum, dum, soutype(k), (souamp(k,i), i=1,6)  ! mIJ sou
  elseif (srctype == 'v' .or. srctype == 'f') then
    read(1,*)  souglobflag, dum, dum, dum, soutype(k), (souamp(k,i), i=1,3)  ! v or f sou
  endif
enddo                                                        
ksou_max = maxval(zsou/dz + 1)  ! max sou depth; used in ghost buf.

read(1,*)
! sou wavelet (before resampling): #points; dt [s]; fmax [Hz]
write(*,*) '..reading nwav_raw, dti_raw, fmax'
read(1,*) nwav_raw, dti_raw, fmax 
allocate (swav_raw(nwav_raw))
allocate (ismax(nwav_raw), ismin(nwav_raw))

write(*,*) '..reading swav_raw'
read(1,*)
do k=1,nwav_raw
  read(1,*) swav_raw(k)  ! read sou wavelet:
enddo

write(*,*) '..FDAniQ anisotropy type (eIso;eVTI;eVOr;eTTI;eTOr;aTTI;aTOr)'
read(1,*)
read(1,*) anitype
if (anitype /= 'eIso' .and. anitype /= 'eVTI' .and. &
  & anitype /= 'eVOr' .and. anitype /= 'eTTI' .and. &
  & anitype /= 'eTOr' .and. anitype /= 'aTTI' .and. &
  & anitype /= 'aTOr') then 
  write(*,*) 'Invalid material anisotropy type!:'
  write(*,*) '..User specified an anitype = ', anitype 
  write(*,*) '..Valid anitype are: eIso;eVTI;eVOr;eTTI;eTOr;aTTI;aTOr'
  write(*,*) '..Terminating program execution!'
  stop
endif

write(*,*) '..reading FDAniQ material property input file names'
read(1,*)
read(1,*,iostat=istat) densityFile
read(1,*,iostat=istat) vp0File 
read(1,*,iostat=istat) vs0File
read(1,*,iostat=istat) epsilon1File
read(1,*,iostat=istat) epsilon2File
read(1,*,iostat=istat) gamma1File
read(1,*,iostat=istat) gamma2File
read(1,*,iostat=istat) delta1File
read(1,*,iostat=istat) delta2File
read(1,*,iostat=istat) delta3File
read(1,*,iostat=istat) qFile
read(1,*,iostat=istat) azimuthFile
read(1,*,iostat=istat) dipFile
read(1,*,iostat=istat) rakeFile

! ..read center frequency for SLS Q mechanism:
write(*,*) '..reading fq'
read(1,*)
read(1,*) fq

! ..read swap4bytes flag for big endian <--> little endian file conversion:
write(*,*) '..reading bigendianflag'
read(1,*)
read(1,*) bigendianflag

! ..read seafloor smoothing flag:
write(*,*) '..reading smoothflag'
read(1,*)
read(1,*) smoothflag

! ..read number of threads (i.e., NCPUS) for OpenMP:
write(*,*) '..reading nthreads'
read(1,*)
read(1,*) nthreads
write(*,*) 'OpenMP threads =', nthreads

call omp_set_num_threads(nthreads)  ! set the number of OpenMP threads

! ..read order of FD time integration:
write(*,*) '..reading ordertime'
read(1,*)
read(1,*) ordertime
write(*,*) 'FD time integration order =', ordertime
if (ordertime == 2) then
  nsti = 2  ! index for sti symplectic time integration sub-loop
else
  nsti = 5
endif

! ..read dti override; allows user to specify FD time step dti:
write(*,*) '..reading dti_overrideflag'
read(1,*)
read(1,*) dti_overrideflag

if (dti_overrideflag) then
  backspace(unit=1)
  write(*,*) '..reading dti_override'
  read(1,*) dti_overrideflag, dti_override
endif

! ..read ABC override; allows user to specify #pts in ABC:
write(*,*) '..ABC override'
read(1,*)
read(1,*) abcoverrideflag

if (abcoverrideflag) then
  backspace(unit=1)
  write(*,*) '..reading nabc_top, nabc_bot, nabc_sdx, nabc_sdy'
  read(1,*) abcoverrideflag, nabc_top, nabc_bot, nabc_sdx, nabc_sdy
else  ! set to default values
  write(*,*) '..using default nabc_top, nabc_bot, nabc_sdx, nabc_sdy'
  nabc_top = nabc_top_d
  nabc_bot = nabc_bot_d
  nabc_sdx = nabc_sdx_d
  nabc_sdy = nabc_sdy_d
endif
write(*,*) '....nabc_top = ', nabc_top
write(*,*) '....nabc_bot = ', nabc_bot
write(*,*) '....nabc_sdx = ', nabc_sdx
write(*,*) '....nabc_sdy = ', nabc_sdy

! Read fmax override flag; resets fmax to minimize disperison in water:
write(*,*) '..fmax override flag'
read(1,*)
read(1,*) fmaxoverrideflag

if (fmaxoverrideflag) then
  backspace(unit=1)
  write(*,*) '..reading vmin'
  read(1,*) fmaxoverrideflag, vmin
  dl_max = max(dx,dz)  ! max grid size
  f40db = vmin/(2.3*dl_max)  ! use 2.3 pts/lambda_min for O(16) FD
  fmax_h2o = f40db/2.5  ! 2.5 comes from -40dB in 5-pole Butterworth
  if (fmax > fmax_h2o) then 
    write(*,*) '....overriding input fmax to minimize:'
    write(*,*) '......vmin = ', vmin
    write(*,*) '......max grid size = ', dl_max
    write(*,*) '......f40dB = ', f40db
    write(*,*) '......fmax (user) = ', fmax
    write(*,*) '......fmax (H2O) = ', fmax_h2o
    fmax = fmax_h2o  ! reset fmax to minimize dispersion in H2O
  else
    write(*,*) '....input fmax is adequte for minimal dispersion in H2O:'
    write(*,*) '......fmax (user) = ', fmax
    write(*,*) '......fmax (H2O) = ', fmax_h2o
  endif
endif

close(unit=1)  ! close FDAniQ.in

!--ALLOCATE ARRAYS---------------------------------------------------------

! ..determine buffer zone required for sou & rec ghosting
if (.not. freesurfflag) then 
  if (ghostrecflag .and. ghostsouflag) then
    ngho_buf = max(ksou_max,krec_max)
  elseif (ghostrecflag .and. .not. ghostsouflag) then
    ngho_buf = krec_max
  elseif (.not. ghostrecflag .and. ghostsouflag) then
    ngho_buf = ksou_max
  else
    ngho_buf = 0
  endif
endif

! ..set array lower bounds:
ilo = -orderspace/2 + 1
jlo = -orderspace/2 + 1
if (.not. freesurfflag) then  ! top is absorb; add extra slab of nabc_top 
                              !  & ngho_buf before start of absorb region
  klo = -orderspace/2 + 1 - nabc_top - ngho_buf   
else  ! top is a free-surface
  klo = -orderspace/2 + 1              
endif

! ..set array upper bounds:
iup = nx + orderspace/2 
jup = ny + orderspace/2 
kup = nz + orderspace/2 + nabc_bot

! ..set array upper bounds (w/out FD op; equiv to 1:nx & 1:ny):
imin = ilo + orderspace/2
jmin = jlo + orderspace/2
imax = iup - orderspace/2
jmax = jup - orderspace/2

write(*,*) 'ilo, jlo, klo = ', ilo, jlo, klo
write(*,*) 'iup, jup, kup = ', iup, jup, kup

! ..allocate particle velocities:
allocate(vx(ilo:iup,jlo:jup,klo:kup)) 
allocate(vy(ilo:iup,jlo:jup,klo:kup)) 
allocate(vz(ilo:iup,jlo:jup,klo:kup)) 

! ..allocate particle stresses:
allocate(txx(ilo:iup,jlo:jup,klo:kup)) 
allocate(tyy(ilo:iup,jlo:jup,klo:kup)) 
allocate(tzz(ilo:iup,jlo:jup,klo:kup)) 
allocate(txz(ilo:iup,jlo:jup,klo:kup)) 
allocate(tyz(ilo:iup,jlo:jup,klo:kup)) 
allocate(txy(ilo:iup,jlo:jup,klo:kup)) 

! ..allocate particle velocity memory variables:
allocate(sx(ilo:iup,jlo:jup,klo:kup)) 
allocate(sy(ilo:iup,jlo:jup,klo:kup)) 
allocate(sz(ilo:iup,jlo:jup,klo:kup)) 


!--READ IN MODEL PROPERTIES & CONVERT TO CIJ'S-----------------------------

! allocate temp. material properties storage arrays:
allocate(ptemp1(ilo:iup,jlo:jup,klo:kup))
allocate(ptemp2(ilo:iup,jlo:jup,klo:kup))
ptemp1 = 0  ! initialize arrays to zero
ptemp2 = 0

!---- rho (property#1/14) -----------------------------!
! note: need density for all the cases.
allocate(rho(ilo:iup,jlo:jup,klo:kup)) ! allocate array
write(*,*) '..reading density file'
propFile = densityFile
call readModelProperties(nx, ny, nz, ilo, jlo, klo, &
              & iup, jup, kup, orderaxes, propFile, bigendianflag, &
              &  freesurfflag, ptemp1, istat)
if (istat /= 0) then
  write(*,*) '..hardwiring density = 1000. kg/m^3!'
  ptemp1 = 1000.  ! set rho = 1000. kg/m^3
endif
if (maxval(ptemp1) < 10) ptemp1 = ptemp1*1000.  ! convert g/cc to kg/m^3
rho = ptemp1     
write(*,*) 'maxval(rho) = ', maxval(rho) 

!---- c33 (property#2/14) -----------------------------!
! note: need c33 for all the cases.
! note: c33 is computed from [Vp0; rho].
allocate(c33(ilo:iup,jlo:jup,klo:kup)) ! allocate array
write(*,*) '..reading Vp0 file'
propFile = vp0File  ! Vp along local 3-axis:
call readModelProperties(nx, ny, nz, ilo, jlo, klo, &
              & iup, jup, kup, orderaxes, propFile, bigendianflag, &
              &  freesurfflag, ptemp1, istat)
if (istat /= 0) then
  write(*,*) '..hardwiring Vp0 = 1500. m/s!'
  ptemp1 = 1500.  ! set Vp0 = 1500. m/s
endif

    filename3='xz_kurt_Vp'
    open(unit=3, file=filename3, status='unknown', action='write', &
         & iostat=status)
    do k=1,nz
      do i=1,nx
        write(3,*) i,k,ptemp1(i,401,k)
      enddo
      write(3,*) ''
    enddo
    close(unit=3)

c33 = (ptemp1**2)*rho
write(*,*) 'maxval(c33) = ', maxval(c33) 

! use vp0 to define deta_max:
! compute average vpvert along top &  bottom row of model, and use this
! to define the ABC max damping coefficients;
allocate(deta_maxsdx(nz))  ! for ABC
allocate(deta_maxsdy(nz))  ! for ABC
deta_maxsdx = 0  ! zero array 
deta_maxsdy = 0  ! zero array

vpvert_avtop = 0.
vpvert_avbot = 0.
do j = 1, ny
  do i = 1, nx
    vpvert_avtop = vpvert_avtop + ptemp1(i,j,1)
    vpvert_avbot = vpvert_avbot + ptemp1(i,j,nz)
  enddo
enddo
vpvert_avtop = vpvert_avtop/float(nx*ny)
vpvert_avbot = vpvert_avbot/float(nx*ny)
write(*,*) '....vpvert_avtop = ', vpvert_avtop
write(*,*) '....vpvert_avbot = ', vpvert_avbot

! for top and bot deta_max, use average values:
deta_maxtop = 18.*vpvert_avtop/(float(nabc_top)*dz)
deta_maxbot = 18.*vpvert_avbot/(float(nabc_bot)*dz)

! for x & y sides deta_max, use linear gradient from top to bot:
deta_maxsdx_top = 18.*vpvert_avtop/(float(nabc_sdx)*dx)
deta_maxsdx_bot = 18.*vpvert_avbot/(float(nabc_sdx)*dx)
deta_maxsdy_top = 18.*vpvert_avtop/(float(nabc_sdy)*dy)
deta_maxsdy_bot = 18.*vpvert_avbot/(float(nabc_sdy)*dy)
deta_gradx = (deta_maxsdx_bot - deta_maxsdx_top)/float(nz-1)
deta_grady = (deta_maxsdy_bot - deta_maxsdy_top)/float(nz-1)
do k = 1, nz
  deta_maxsdx(k) = deta_maxsdx_top + deta_gradx*float(k-1)
  deta_maxsdy(k) = deta_maxsdy_top + deta_grady*float(k-1)
enddo

!---- c11 (property#3/14) -----------------------------!
! note: need c11 for all the anisotropic cases.
! note: c11 is computed from [epsilon2; c33].
if (anitype == 'eVTI' .or. anitype == 'eVOr' .or. &
  & anitype == 'eTTI' .or. anitype == 'eTOr' .or. &
  & anitype == 'aTTI' .or. anitype == 'aTOr') then
  allocate(c11(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..reading epsilon2 file'
  propFile = epsilon2File  ! P-wave ani btwn x and z axes
  call readModelProperties(nx, ny, nz, ilo, jlo, klo, &
                & iup, jup, kup, orderaxes, propFile, bigendianflag, &
                & freesurfflag, ptemp1, istat)
  c11 = c33*(1. + 2.*ptemp1)
  write(*,*) 'maxval(c11) = ', maxval(c11) 
endif

!---- c22 (property#4/14) -----------------------------!
! note: need c22 for the orthorhombic and TTI cases.

! note: for ortho cases, c22 is computed from [epsilon1; c33].
if (anitype == 'eVOr' .or. anitype == 'eTOr' .or. anitype == 'aTOr') then 
  allocate(c22(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..reading epsilon1 file'
  propFile = epsilon1File  ! P-wave ani btwn y and z axes
  call readModelProperties(nx, ny, nz, ilo, jlo, klo, &
                & iup, jup, kup, orderaxes, propFile, bigendianflag, &
                & freesurfflag, ptemp1, istat)
  c22 = c33*(1. + 2.*ptemp1)
  write(*,*) 'maxval(c22) = ', maxval(c22) 

! note: for TTI cases, c22 is equal to c11.
elseif (anitype == 'eTTI' .or. anitype == 'aTTI') then 
  allocate(c22(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..setting c22 = c11 for TTI case'
  c22 = c11  ! set c22 = c11 for TTI cases
  write(*,*) 'maxval(c22) = ', maxval(c22) 
endif

!---- c55 (property#5/14) -----------------------------!
! note: need c55 for all the cases.
allocate(c55(ilo:iup,jlo:jup,klo:kup)) ! allocate array

! note: for the elastic cases, c55 is computed from [Vs0; rho].
if (anitype == 'eIso' .or. anitype == 'eVTI' .or. &
  & anitype == 'eVOr' .or. anitype == 'eTTI' .or. &
  & anitype == 'eTOr') then 
  write(*,*) '..reading Vs0 file'
  propFile = vs0File  ! Vs along local 3-axis w/ polariz along 1-axis
  call readModelProperties(nx, ny, nz, ilo, jlo, klo, &
                & iup, jup, kup, orderaxes, propFile, bigendianflag, &
                & freesurfflag, ptemp2, istat)
  if (istat /= 0) then
    write(*,*) '..hardwiring Vs0 = Vp0/1.78 (Poisson ratio = 0.27)!'
    ptemp2 = ptemp1/1.78  ! set Vs0 =  Vp0/1.78
  endif

    filename3='xz_kurt_Vs'
    open(unit=3, file=filename3, status='unknown', action='write', &
         & iostat=status)
    do k=1,nz
      do i=1,nx
        write(3,*) i,k,ptemp2(i,401,k)
      enddo
      write(3,*) ''
    enddo
    close(unit=3)

  c55 = (ptemp2**2)*rho
  write(*,*) 'maxval(c55) = ', maxval(c55) 

! note: for the acoustic ani cases, c55 is computed from [c33; c55overc33].
elseif (anitype == 'aTTI' .or. anitype == 'aTOr') then
  write(*,*) '..setting c55 = c33*ratio[c55/c33] for aco ani case'
  c55 = c33*c55overc33  ! hardwire c55 using [c55/c33] ratio 
  write(*,*) 'maxval(c55) = ', maxval(c55) 
endif

!---- c66 (property#6/14) -----------------------------!
! note: need c66 for all anisotropic cases.

! note: for elastic ani cases, c66 is computed from [gamma1; c55].
if (anitype == 'eVTI' .or. anitype == 'eVOr' .or. &
  & anitype == 'eTTI' .or. anitype == 'eTOr') then 
  allocate(c66(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..reading gamma1 file'
  propFile = gamma1File  ! SH-wave ani for prop in x-z plane
  call readModelProperties(nx, ny, nz, ilo, jlo, klo, &
                & iup, jup, kup, orderaxes, propFile, bigendianflag, &
                & freesurfflag, ptemp1, istat)
  c66 = c55*(1. + 2.*ptemp1)
  write(*,*) 'maxval(c66) = ', maxval(c66) 

! note: for acoustic ani cases, c66 is computed from [c33, c66overc33].
elseif (anitype == 'aTTI' .or. anitype == 'aTOr') then
  allocate(c66(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..setting c66 = c33*ratio[c66/c33] for aco ani case'
  c66 = c33*c66overc33  ! hardwire c66 for acoustic ani cases
  write(*,*) 'maxval(c66) = ', maxval(c66) 
endif

!---- c44 (property#7/14) -----------------------------!
! note: need c44 for the orthorhombic and TTI cases.

! note: for ortho cases, c44 is computed from [gamma2; c66].
if (anitype == 'eVOr' .or. anitype == 'eTOr' .or. anitype == 'aTOr') then 
  allocate(c44(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..reading gamma2 file'
  propFile = gamma2File  ! SH-wave ani for prop in y-z plane
  call readModelProperties(nx, ny, nz, ilo, jlo, klo, &
                & iup, jup, kup, orderaxes, propFile, bigendianflag, &
                & freesurfflag, ptemp1, istat)
  c44 = c66/(1. + 2.*ptemp1)
  write(*,*) 'maxval(c44) = ', maxval(c44)

! note: for TTI cases, c44 is equal to c55.
elseif (anitype == 'eTTI' .or. anitype == 'aTTI') then
  allocate(c44(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..setting c44 = c55 for TTI case'
  c44 = c55  ! set c44 = c55 for TTI cases
  write(*,*) 'maxval(c44) = ', maxval(c44)
endif

!---- c13 (property#8/14) -----------------------------!
! note: need c13 for all the anisotropic cases.
! note: c13 is computed from [delta2; c33; c55].
if (anitype == 'eVTI' .or. anitype == 'eVOr' .or. &
  & anitype == 'eTTI' .or. anitype == 'eTOr' .or. &
  & anitype == 'aTTI' .or. anitype == 'aTOr') then
  allocate(c13(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..reading delta2 file'
  propFile = delta2File  ! delta in local x-z plane
  call readModelProperties(nx, ny, nz, ilo, jlo, klo, &
                & iup, jup, kup, orderaxes, propFile, bigendianflag, &
                & freesurfflag, ptemp1, istat)
  c13 = sqrt((c33 - c55)*(2*ptemp1*c33 + c33 - c55)) - c55
  write(*,*) 'maxval(c13) = ', maxval(c13) 
endif

!---- c12 (property#9/14) -----------------------------!
! note: need c12 for the orthorhombic and TTI cases.

! note: for ortho cases, c12 is computed from [delta3; c11; c66].
if (anitype == 'eVOr' .or. anitype == 'eTOr' .or. anitype == 'aTOr') then 
  allocate(c12(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..reading delta3 file'
  propFile = delta3File  ! delta in local x-y plane
  call readModelProperties(nx, ny, nz, ilo, jlo, klo, &
                & iup, jup, kup, orderaxes, propFile, bigendianflag, &
                & freesurfflag, ptemp1, istat)
  c12 = sqrt((c11 - c66)*(2*ptemp1*c11 + c11 - c66)) - c66
  write(*,*) 'maxval(c12) = ', maxval(c12) 

! note: for TTI cases, c12 is computed from [c11; c66].
elseif (anitype == 'eTTI' .or. anitype == 'aTTI') then
  allocate(c12(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..setting c13 = (c11 - 2c66) for TTI case'
  c12 = c11 - 2.*c66  ! c12 = c11 - 2c66 for TTI cases
endif

!---- c23 (property#10/14) ----------------------------!
! note: need c23 for the orthorhombic and TTI cases.

! note: for ortho cases, c23 is computed from [delta1; c33; c44].
if (anitype == 'eVOr' .or. anitype == 'eTOr' .or. anitype == 'aTOr') then 
  allocate(c23(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..reading delta1 file'
  propFile = delta1File  ! delta in local y-z plane
  call readModelProperties(nx, ny, nz, ilo, jlo, klo, &
                & iup, jup, kup, orderaxes, propFile, bigendianflag, &
                & freesurfflag, ptemp1, istat)
  c23 = sqrt((c33 - c44)*(2*ptemp1*c33 + c33 - c44)) - c44
  write(*,*) 'maxval(c23) = ', maxval(c23) 

! note: for TTI cases, c23 is equal to c13.
elseif (anitype == 'eTTI' .or. anitype == 'aTTI') then
  allocate(c23(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..setting c23 = c13 for TTI case'
  c23 = c13  ! set c23 = c13 for TTI cases
  write(*,*) 'maxval(c23) = ', maxval(c23) 
endif

!---- azimuth (property#11/14) ------------------------!
! note: need azim for TTI and TOr cases.
! note: azimuth is a c-clock rot around Z(0)-axis from X-axis.
if (anitype == 'eTTI' .or. anitype == 'eTOr' .or. &
  & anitype == 'aTTI' .or. anitype == 'aTOr') then 
  allocate(azim(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..reading azimuth angle file'
  propFile = azimuthFile
  call readModelProperties(nx, ny, nz, ilo, jlo, klo, &
                & iup, jup, kup, orderaxes, propFile, bigendianflag, &
                & freesurfflag, ptemp1, istat)
  if (maxval(ptemp1) > 2.*pi) ptemp1 = ptemp1/180.  ! conv from deg to rad
  azim = ptemp1
  write(*,*) 'maxval(azim) = ', maxval(azim) 
endif

!---- dip (property#12/14) ----------------------------!
! note: need dip for TTI and TOr cases.
! note: dip is a c-clock rot around Y(1)-axis.
if (anitype == 'eTTI' .or. anitype == 'eTOr' .or. &
  & anitype == 'aTTI' .or. anitype == 'aTOr') then 
  allocate(dip(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..reading dip angle file'
  propFile = dipFile
  call readModelProperties(nx, ny, nz, ilo, jlo, klo, &
                & iup, jup, kup, orderaxes, propFile, bigendianflag, &
                & freesurfflag, ptemp1, istat)
  if (maxval(ptemp1) > 0.5*pi) ptemp1 = ptemp1/180.  ! conv from deg to rad
  dip = ptemp1
  write(*,*) 'maxval(dip) = ', maxval(dip) 
endif

!---- rake (property#13/14) ---------------------------!
! note: need rake for TOr cases.
! note: rake is a c-clock rot around Z(2)-axis).
if (anitype == 'eTOr' .or. anitype == 'aTOr') then
  allocate(rake(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..reading rake angle file'
  propFile = rakeFile
  call readModelProperties(nx, ny, nz, ilo, jlo, klo, &
                & iup, jup, kup, orderaxes, propFile, bigendianflag, &
                & freesurfflag, ptemp1, istat)
  if (maxval(ptemp1) > 2.*pi) ptemp1 = ptemp1/180.  ! conv from deg to rad
  rake = ptemp1
  write(*,*) 'maxval(rake) = ', maxval(rake) 

elseif (anitype == 'eTTI' .or. anitype == 'aTTI') then
  allocate(rake(ilo:iup,jlo:jup,klo:kup)) ! allocate array
  write(*,*) '..setting rake = 0 for TTI case'
  rake = 0  ! set rake = 0 for TTI cases
  write(*,*) 'maxval(rake) = ', maxval(rake) 
endif


!---- Smooth Seafloor Option --------------------------!
! Smooth Seafloor Option:  do this before defining visco prop's!
! To Do: enable seafloor smoothing option
!!!!!if (smoothflag) then
!!!!!  call smoothSeafloor(ilo, jlo, klo, iup, jup, kup, nx, ny, nz, &
!!!!!                    & dx, dy, dz, c11, c33, c55, c66, c13, rho)
!!!!!write(*,*) '..smoothing seafloor'
!!!!!endif

!---- Q (property#14/14) ------------------------------!
! note: need Q for all cases.
! note: Standard Linear Solid (SLS) Q is applied equally to [vz; vy; vz].
allocate(itausig(ilo:iup,jlo:jup,klo:kup)) ! allocate array
allocate(difitau(ilo:iup,jlo:jup,klo:kup)) ! allocate array
itausig = 0  ! zero array
difitau = 0  ! zero array

! compute center frequency for SLS Q mechanism:
if (fq <= 0.) then  ! compute dominant freq from wavelet peak-to-peak
  write(*,*) 'User-supplied fq is non-physical; fq [Hz] = ', fq
  ismax = maxloc(swav_raw) 
  ismin = minloc(swav_raw)
  dt_ptp = 2.*dti_raw*abs(float(ismax(1) - ismin(1)))
  fq = 1./dt_ptp
  write(*,*) 'computing fq from wavelet peak-to-peak ...'
  write(*,*) 'fq = ', fq
endif

! read Q (note: assuming that Qp = Qs = Q):
write(*,*) '..reading Q file'
propFile = qFile
call readModelProperties(nx, ny, nz, ilo, jlo, klo, &
              & iup, jup, kup, orderaxes, propFile, bigendianflag, &
              & freesurfflag, ptemp1, istat)

if (istat /= 0) then
  write(*,*) 'Setting Q = ', qbig
  ptemp1 = qbig
else
  ! QC Q:
  write(*,*) 'Note: if Q <= 0 is encountered, Q will be set = ', qbig
  do k=1,nz
    do j=1,ny
      do i=1,nx
        if (ptemp1(i,j,k) <= 0) ptemp1(i,j,k) = qbig
      enddo
    enddo
  enddo
endif
write(*,*) 'maxval(Q) = ', maxval(ptemp1) 

! convert quality factors to SLS viscoelastic constants:
do k=1,nz
  do j=1,ny
    do i=1,nx
      wq = 2.*pi*fq  ! freq. corresponding to input Q
      te = (1 + sqrt(1 + ptemp1(i,j,k)**2))/(ptemp1(i,j,k)*wq)  
      tt = 1/(te*wq**2)
      itausig(i,j,k) = 1./tt
      difitau(i,j,k) = (1./te - itausig(i,j,k))
    enddo
  enddo
enddo

! End of material property definitions:
deallocate(ptemp1, ptemp2)


!--ROTATE CIJ FROM LOCAL (eTTI; eTOr) FRAME TO GLOBAL (FD) FRAME-----------
if (anitype == 'eTTI' .or. anitype == 'eTOr' .or. &
  & anitype == 'aTTI' .or. anitype == 'aTOr') then

  ! Allocate off-kite cIJ's:
  allocate(c14(ilo:iup,jlo:jup,klo:kup)) 
  allocate(c15(ilo:iup,jlo:jup,klo:kup)) 
  allocate(c16(ilo:iup,jlo:jup,klo:kup)) 
  allocate(c24(ilo:iup,jlo:jup,klo:kup)) 
  allocate(c25(ilo:iup,jlo:jup,klo:kup)) 
  allocate(c26(ilo:iup,jlo:jup,klo:kup)) 
  allocate(c34(ilo:iup,jlo:jup,klo:kup)) 
  allocate(c35(ilo:iup,jlo:jup,klo:kup)) 
  allocate(c36(ilo:iup,jlo:jup,klo:kup)) 
  allocate(c45(ilo:iup,jlo:jup,klo:kup)) 
  allocate(c46(ilo:iup,jlo:jup,klo:kup)) 
  allocate(c56(ilo:iup,jlo:jup,klo:kup)) 

  ! Bond transf of CIJ from local (eTTI;eTOr) to global (FD) frame:
  write(*,*) 'Rotating cIJ with Bond transformation'
  call rotateCIJLocalToGlobal(nx, ny, nz, ilo, jlo, klo, iup, jup, kup, &
                    & azim, dip, rake, &
                    & c11, c22, c33, c44, c55, c66, &
                    & c12, c13, c14, c15, c16, &
                    & c23, c24, c25, c26, &
                    & c34, c35, c36, &
                    & c45, c46, &
                    & c56)
  ! Don't need 3 rot angles anymore:
  deallocate(azim, dip, rake)

endif


!--EXTRAPOLATE(1D) TOP & BOT PROPS AS NEEDED FOR ABC'S---------------------

! Fill in properties on sides (x & y slabs) before z slab extrapolations:  
do k= 1, nz  
  do j= 1, ny  
    do i= ilo, 0  ! -x slab
      rho(i,j,k) = rho(1,j,k)
      itausig(i,j,k) = itausig(1,j,k)
      difitau(i,j,k) = difitau(1,j,k)
      select case (anitype)  ! eIso; eVTI; eVOr; eTTI; eTOr; aTTI; aTOr
        case ("eIso")  ! 2 cIJ's
          c33(i,j,k) = c33(1,j,k)
          c55(i,j,k) = c55(1,j,k)
        case ("eVTI")  ! 5 cIJ's
          c11(i,j,k) = c11(1,j,k)
          c33(i,j,k) = c33(1,j,k)
          c55(i,j,k) = c55(1,j,k)
          c66(i,j,k) = c66(1,j,k)
          c13(i,j,k) = c13(1,j,k)
        case ("eVOr")  ! 9 cIJ's
          c11(i,j,k) = c11(1,j,k)
          c22(i,j,k) = c22(1,j,k)
          c33(i,j,k) = c33(1,j,k)
          c44(i,j,k) = c44(1,j,k)
          c55(i,j,k) = c55(1,j,k)
          c66(i,j,k) = c66(1,j,k)
          c12(i,j,k) = c12(1,j,k)
          c13(i,j,k) = c13(1,j,k)
          c23(i,j,k) = c23(1,j,k)
        case ("eTTI", "eTOr", "aTTI", "aTOr")  ! 21 cIJ's
          c11(i,j,k) = c11(1,j,k)
          c22(i,j,k) = c22(1,j,k)
          c33(i,j,k) = c33(1,j,k)
          c44(i,j,k) = c44(1,j,k)
          c55(i,j,k) = c55(1,j,k)
          c66(i,j,k) = c66(1,j,k)
          c12(i,j,k) = c12(1,j,k)
          c13(i,j,k) = c13(1,j,k)
          c14(i,j,k) = c14(1,j,k)
          c15(i,j,k) = c15(1,j,k)
          c16(i,j,k) = c16(1,j,k)
          c23(i,j,k) = c23(1,j,k)
          c24(i,j,k) = c24(1,j,k)
          c25(i,j,k) = c25(1,j,k)
          c26(i,j,k) = c26(1,j,k)
          c34(i,j,k) = c34(1,j,k)
          c35(i,j,k) = c35(1,j,k)
          c36(i,j,k) = c36(1,j,k)
          c45(i,j,k) = c45(1,j,k)
          c46(i,j,k) = c46(1,j,k)
          c56(i,j,k) = c56(1,j,k)
      end select
    enddo
    do i= nx+1, iup  ! +x slab
      rho(i,j,k) = rho(nx,j,k)
      itausig(i,j,k) = itausig(nx,j,k)
      difitau(i,j,k) = difitau(nx,j,k)
      select case (anitype)  ! eIso; eVTI; eVOr; eTTI; eTOr; aTTI; aTOr
        case ("eIso")  ! 2 cIJ's
          c33(i,j,k) = c33(nx,j,k)
          c55(i,j,k) = c55(nx,j,k)
        case ("eVTI")  ! 5 cIJ's
          c11(i,j,k) = c11(nx,j,k)
          c33(i,j,k) = c33(nx,j,k)
          c55(i,j,k) = c55(nx,j,k)
          c66(i,j,k) = c66(nx,j,k)
          c13(i,j,k) = c13(nx,j,k)
        case ("eVOr")  ! 9 cIJ's
          c11(i,j,k) = c11(nx,j,k)
          c22(i,j,k) = c22(nx,j,k)
          c33(i,j,k) = c33(nx,j,k)
          c44(i,j,k) = c44(nx,j,k)
          c55(i,j,k) = c55(nx,j,k)
          c66(i,j,k) = c66(nx,j,k)
          c12(i,j,k) = c12(nx,j,k)
          c13(i,j,k) = c13(nx,j,k)
          c23(i,j,k) = c23(nx,j,k)
        case ("eTTI", "eTOr", "aTTI", "aTOr")  ! 21 cIJ's
          c11(i,j,k) = c11(nx,j,k)
          c22(i,j,k) = c22(nx,j,k)
          c33(i,j,k) = c33(nx,j,k)
          c44(i,j,k) = c44(nx,j,k)
          c55(i,j,k) = c55(nx,j,k)
          c66(i,j,k) = c66(nx,j,k)
          c12(i,j,k) = c12(nx,j,k)
          c13(i,j,k) = c13(nx,j,k)
          c14(i,j,k) = c14(nx,j,k)
          c15(i,j,k) = c15(nx,j,k)
          c16(i,j,k) = c16(nx,j,k)
          c23(i,j,k) = c23(nx,j,k)
          c24(i,j,k) = c24(nx,j,k)
          c25(i,j,k) = c25(nx,j,k)
          c26(i,j,k) = c26(nx,j,k)
          c34(i,j,k) = c34(nx,j,k)
          c35(i,j,k) = c35(nx,j,k)
          c36(i,j,k) = c36(nx,j,k)
          c45(i,j,k) = c45(nx,j,k)
          c46(i,j,k) = c46(nx,j,k)
          c56(i,j,k) = c56(nx,j,k)
      end select
    enddo
  enddo
enddo
  
do k= 1, nz  
  do i= ilo, iup  
    do j= jlo, 0  ! -y slab
      rho(i,j,k) = rho(i,1,k)
      itausig(i,j,k) = itausig(i,1,k)
      difitau(i,j,k) = difitau(i,1,k)
      select case (anitype)  ! eIso; eVTI; eVOr; eTTI; eTOr; aTTI; aTOr
        case ("eIso")  ! 2 cIJ's
          c33(i,j,k) = c33(i,1,k)
          c55(i,j,k) = c55(i,1,k)
        case ("eVTI")  ! 5 cIJ's
          c11(i,j,k) = c11(i,1,k)
          c33(i,j,k) = c33(i,1,k)
          c55(i,j,k) = c55(i,1,k)
          c66(i,j,k) = c66(i,1,k)
          c13(i,j,k) = c13(i,1,k)
        case ("eVOr")  ! 9 cIJ's
          c11(i,j,k) = c11(i,1,k)
          c22(i,j,k) = c22(i,1,k)
          c33(i,j,k) = c33(i,1,k)
          c44(i,j,k) = c44(i,1,k)
          c55(i,j,k) = c55(i,1,k)
          c66(i,j,k) = c66(i,1,k)
          c12(i,j,k) = c12(i,1,k)
          c13(i,j,k) = c13(i,1,k)
          c23(i,j,k) = c23(i,1,k)
        case ("eTTI", "eTOr", "aTTI", "aTOr")  ! 21 cIJ's
          c11(i,j,k) = c11(i,1,k)
          c22(i,j,k) = c22(i,1,k)
          c33(i,j,k) = c33(i,1,k)
          c44(i,j,k) = c44(i,1,k)
          c55(i,j,k) = c55(i,1,k)
          c66(i,j,k) = c66(i,1,k)
          c12(i,j,k) = c12(i,1,k)
          c13(i,j,k) = c13(i,1,k)
          c14(i,j,k) = c14(i,1,k)
          c15(i,j,k) = c15(i,1,k)
          c16(i,j,k) = c16(i,1,k)
          c23(i,j,k) = c23(i,1,k)
          c24(i,j,k) = c24(i,1,k)
          c25(i,j,k) = c25(i,1,k)
          c26(i,j,k) = c26(i,1,k)
          c34(i,j,k) = c34(i,1,k)
          c35(i,j,k) = c35(i,1,k)
          c36(i,j,k) = c36(i,1,k)
          c45(i,j,k) = c45(i,1,k)
          c46(i,j,k) = c46(i,1,k)
          c56(i,j,k) = c56(i,1,k)
      end select
    enddo
    do j= ny+1, jup  ! +y slab
      rho(i,j,k) = rho(i,ny,k)
      itausig(i,j,k) = itausig(i,ny,k)
      difitau(i,j,k) = difitau(i,ny,k)
      select case (anitype)  ! eIso; eVTI; eVOr; eTTI; eTOr; aTTI; aTOr
        case ("eIso")  ! 2 cIJ's
          c33(i,j,k) = c33(i,ny,k)
          c55(i,j,k) = c55(i,ny,k)
        case ("eVTI")  ! 5 cIJ's
          c11(i,j,k) = c11(i,ny,k)
          c33(i,j,k) = c33(i,ny,k)
          c55(i,j,k) = c55(i,ny,k)
          c66(i,j,k) = c66(i,ny,k)
          c13(i,j,k) = c13(i,ny,k)
        case ("eVOr")  ! 9 cIJ's
          c11(i,j,k) = c11(i,ny,k)
          c22(i,j,k) = c22(i,ny,k)
          c33(i,j,k) = c33(i,ny,k)
          c44(i,j,k) = c44(i,ny,k)
          c55(i,j,k) = c55(i,ny,k)
          c66(i,j,k) = c66(i,ny,k)
          c12(i,j,k) = c12(i,ny,k)
          c13(i,j,k) = c13(i,ny,k)
          c23(i,j,k) = c23(i,ny,k)
        case ("eTTI", "eTOr", "aTTI", "aTOr")  ! 21 cIJ's
          c11(i,j,k) = c11(i,ny,k)
          c22(i,j,k) = c22(i,ny,k)
          c33(i,j,k) = c33(i,ny,k)
          c44(i,j,k) = c44(i,ny,k)
          c55(i,j,k) = c55(i,ny,k)
          c66(i,j,k) = c66(i,ny,k)
          c12(i,j,k) = c12(i,ny,k)
          c13(i,j,k) = c13(i,ny,k)
          c14(i,j,k) = c14(i,ny,k)
          c15(i,j,k) = c15(i,ny,k)
          c16(i,j,k) = c16(i,ny,k)
          c23(i,j,k) = c23(i,ny,k)
          c24(i,j,k) = c24(i,ny,k)
          c25(i,j,k) = c25(i,ny,k)
          c26(i,j,k) = c26(i,ny,k)
          c34(i,j,k) = c34(i,ny,k)
          c35(i,j,k) = c35(i,ny,k)
          c36(i,j,k) = c36(i,ny,k)
          c45(i,j,k) = c45(i,ny,k)
          c46(i,j,k) = c46(i,ny,k)
          c56(i,j,k) = c56(i,ny,k)
      end select
    enddo
  enddo
enddo
  
do i= ilo, iup  
  do j= jlo, jup  
    do k= klo, 0  ! -z slab
      rho(i,j,k) = rho(i,j,1)
      itausig(i,j,k) = itausig(i,j,1)
      difitau(i,j,k) = difitau(i,j,1)
      select case (anitype)  ! eIso; eVTI; eVOr; eTTI; eTOr; aTTI; aTOr
        case ("eIso")  ! 2 cIJ's
          c33(i,j,k) = c33(i,j,1)
          c55(i,j,k) = c55(i,j,1)
        case ("eVTI")  ! 5 cIJ's
          c11(i,j,k) = c11(i,j,1)
          c33(i,j,k) = c33(i,j,1)
          c55(i,j,k) = c55(i,j,1)
          c66(i,j,k) = c66(i,j,1)
          c13(i,j,k) = c13(i,j,1)
        case ("eVOr")  ! 9 cIJ's
          c11(i,j,k) = c11(i,j,1)
          c22(i,j,k) = c22(i,j,1)
          c33(i,j,k) = c33(i,j,1)
          c44(i,j,k) = c44(i,j,1)
          c55(i,j,k) = c55(i,j,1)
          c66(i,j,k) = c66(i,j,1)
          c12(i,j,k) = c12(i,j,1)
          c13(i,j,k) = c13(i,j,1)
          c23(i,j,k) = c23(i,j,1)
        case ("eTTI", "eTOr", "aTTI", "aTOr")  ! 21 cIJ's
          c11(i,j,k) = c11(i,j,1)
          c22(i,j,k) = c22(i,j,1)
          c33(i,j,k) = c33(i,j,1)
          c44(i,j,k) = c44(i,j,1)
          c55(i,j,k) = c55(i,j,1)
          c66(i,j,k) = c66(i,j,1)
          c12(i,j,k) = c12(i,j,1)
          c13(i,j,k) = c13(i,j,1)
          c14(i,j,k) = c14(i,j,1)
          c15(i,j,k) = c15(i,j,1)
          c16(i,j,k) = c16(i,j,1)
          c23(i,j,k) = c23(i,j,1)
          c24(i,j,k) = c24(i,j,1)
          c25(i,j,k) = c25(i,j,1)
          c26(i,j,k) = c26(i,j,1)
          c34(i,j,k) = c34(i,j,1)
          c35(i,j,k) = c35(i,j,1)
          c36(i,j,k) = c36(i,j,1)
          c45(i,j,k) = c45(i,j,1)
          c46(i,j,k) = c46(i,j,1)
          c56(i,j,k) = c56(i,j,1)
      end select
    enddo
    do k= nz+1, kup  ! +z slab  
      rho(i,j,k) = rho(i,j,nz)
      itausig(i,j,k) = itausig(i,j,nz)
      difitau(i,j,k) = difitau(i,j,nz)
      select case (anitype)  ! eIso; eVTI; eVOr; eTTI; eTOr; aTTI; aTOr
        case ("eIso")  ! 2 cIJ's
          c33(i,j,k) = c33(i,j,nz)
          c55(i,j,k) = c55(i,j,nz)
        case ("eVTI")  ! 5 cIJ's
          c11(i,j,k) = c11(i,j,nz)
          c33(i,j,k) = c33(i,j,nz)
          c55(i,j,k) = c55(i,j,nz)
          c66(i,j,k) = c66(i,j,nz)
          c13(i,j,k) = c13(i,j,nz)
        case ("eVOr")  ! 9 cIJ's
          c11(i,j,k) = c11(i,j,nz)
          c22(i,j,k) = c22(i,j,nz)
          c33(i,j,k) = c33(i,j,nz)
          c44(i,j,k) = c44(i,j,nz)
          c55(i,j,k) = c55(i,j,nz)
          c66(i,j,k) = c66(i,j,nz)
          c12(i,j,k) = c12(i,j,nz)
          c13(i,j,k) = c13(i,j,nz)
          c23(i,j,k) = c23(i,j,nz)
        case ("eTTI", "eTOr", "aTTI", "aTOr")  ! 21 cIJ's
          c11(i,j,k) = c11(i,j,nz)
          c22(i,j,k) = c22(i,j,nz)
          c33(i,j,k) = c33(i,j,nz)
          c44(i,j,k) = c44(i,j,nz)
          c55(i,j,k) = c55(i,j,nz)
          c66(i,j,k) = c66(i,j,nz)
          c12(i,j,k) = c12(i,j,nz)
          c13(i,j,k) = c13(i,j,nz)
          c14(i,j,k) = c14(i,j,nz)
          c15(i,j,k) = c15(i,j,nz)
          c16(i,j,k) = c16(i,j,nz)
          c23(i,j,k) = c23(i,j,nz)
          c24(i,j,k) = c24(i,j,nz)
          c25(i,j,k) = c25(i,j,nz)
          c26(i,j,k) = c26(i,j,nz)
          c34(i,j,k) = c34(i,j,nz)
          c35(i,j,k) = c35(i,j,nz)
          c36(i,j,k) = c36(i,j,nz)
          c45(i,j,k) = c45(i,j,nz)
          c46(i,j,k) = c46(i,j,nz)
          c56(i,j,k) = c56(i,j,nz)
      end select
    enddo
  enddo
enddo


! ..QC model prop's after ABC extrapolation:
write(*,*) 'QC of model properties after Bond rot & ABC (1D) extrapol: '
write(*,*) 'minval(rho) = ', minval(rho) 
write(*,*) 'maxval(rho) = ', maxval(rho) 
select case (anitype)  ! eIso; eVTI; eVOr; eTTI; eTOr; aTTI; aTOr
  case ("eIso")  ! 2 cIJ's
    write(*,*) 'minval(c33) = ', minval(c33) 
    write(*,*) 'maxval(c33) = ', maxval(c33) 
    write(*,*) 'minval(c55) = ', minval(c55) 
    write(*,*) 'maxval(c55) = ', maxval(c55) 
  case ("eVTI")  ! 5 cIJ's
    write(*,*) 'minval(c11) = ', minval(c11) 
    write(*,*) 'maxval(c11) = ', maxval(c11) 
    write(*,*) 'minval(c33) = ', minval(c33) 
    write(*,*) 'maxval(c33) = ', maxval(c33) 
    write(*,*) 'minval(c55) = ', minval(c55) 
    write(*,*) 'maxval(c55) = ', maxval(c55) 
    write(*,*) 'minval(c66) = ', minval(c66) 
    write(*,*) 'maxval(c66) = ', maxval(c66) 
    write(*,*) 'minval(c13) = ', minval(c13) 
    write(*,*) 'maxval(c13) = ', maxval(c13) 
  case ("eVOr")  ! 9 cIJ's
    write(*,*) 'minval(c11) = ', minval(c11) 
    write(*,*) 'maxval(c11) = ', maxval(c11) 
    write(*,*) 'minval(c22) = ', minval(c22) 
    write(*,*) 'maxval(c22) = ', maxval(c22) 
    write(*,*) 'minval(c33) = ', minval(c33) 
    write(*,*) 'maxval(c33) = ', maxval(c33) 
    write(*,*) 'minval(c44) = ', minval(c44) 
    write(*,*) 'maxval(c44) = ', maxval(c44) 
    write(*,*) 'minval(c55) = ', minval(c55) 
    write(*,*) 'maxval(c55) = ', maxval(c55) 
    write(*,*) 'minval(c66) = ', minval(c66) 
    write(*,*) 'maxval(c66) = ', maxval(c66) 
    write(*,*) 'minval(c12) = ', minval(c12) 
    write(*,*) 'maxval(c12) = ', maxval(c12) 
    write(*,*) 'minval(c13) = ', minval(c13) 
    write(*,*) 'maxval(c13) = ', maxval(c13) 
    write(*,*) 'minval(c23) = ', minval(c23) 
    write(*,*) 'maxval(c23) = ', maxval(c23) 
  case ("eTTI", "eTOr", "aTTI", "aTOr")  ! 21 cIJ's
    write(*,*) 'minval(c11) = ', minval(c11) 
    write(*,*) 'maxval(c11) = ', maxval(c11) 
    write(*,*) 'minval(c22) = ', minval(c22) 
    write(*,*) 'maxval(c22) = ', maxval(c22) 
    write(*,*) 'minval(c33) = ', minval(c33) 
    write(*,*) 'maxval(c33) = ', maxval(c33) 
    write(*,*) 'minval(c44) = ', minval(c44) 
    write(*,*) 'maxval(c44) = ', maxval(c44) 
    write(*,*) 'minval(c55) = ', minval(c55) 
    write(*,*) 'maxval(c55) = ', maxval(c55) 
    write(*,*) 'minval(c66) = ', minval(c66) 
    write(*,*) 'maxval(c66) = ', maxval(c66) 
    write(*,*) 'minval(c12) = ', minval(c12) 
    write(*,*) 'maxval(c12) = ', maxval(c12) 
    write(*,*) 'minval(c13) = ', minval(c13) 
    write(*,*) 'maxval(c13) = ', maxval(c13) 
    write(*,*) 'minval(c14) = ', minval(c14) 
    write(*,*) 'maxval(c14) = ', maxval(c14) 
    write(*,*) 'minval(c15) = ', minval(c15) 
    write(*,*) 'maxval(c15) = ', maxval(c15) 
    write(*,*) 'minval(c16) = ', minval(c16) 
    write(*,*) 'maxval(c16) = ', maxval(c16) 
    write(*,*) 'minval(c23) = ', minval(c23) 
    write(*,*) 'maxval(c23) = ', maxval(c23) 
    write(*,*) 'minval(c24) = ', minval(c24) 
    write(*,*) 'maxval(c24) = ', maxval(c24) 
    write(*,*) 'minval(c25) = ', minval(c25) 
    write(*,*) 'maxval(c25) = ', maxval(c25) 
    write(*,*) 'minval(c26) = ', minval(c26) 
    write(*,*) 'maxval(c26) = ', maxval(c26) 
    write(*,*) 'minval(c34) = ', minval(c34) 
    write(*,*) 'maxval(c34) = ', maxval(c34) 
    write(*,*) 'minval(c35) = ', minval(c35) 
    write(*,*) 'maxval(c35) = ', maxval(c35) 
    write(*,*) 'minval(c36) = ', minval(c36) 
    write(*,*) 'maxval(c36) = ', maxval(c36) 
    write(*,*) 'minval(c45) = ', minval(c45) 
    write(*,*) 'maxval(c45) = ', maxval(c45) 
    write(*,*) 'minval(c46) = ', minval(c46) 
    write(*,*) 'maxval(c46) = ', maxval(c46) 
    write(*,*) 'minval(c56) = ', minval(c56) 
    write(*,*) 'maxval(c56) = ', maxval(c56) 
end select


!--SET FD TIME STEPPING PARAMETERS-----------------------------------------

! ..find smallest grid size and maximum P velocity in model:
dl_min = min(dx, dy)
dl_min = min(dl_min, dz)
vp_max_c11 = maxval(sqrt(c11/rho))
vp_max_c33 = maxval(sqrt(c33/rho))
vp_max = max(vp_max_c11, vp_max_c33)
if (anitype == 'VOr' .or. anitype == 'TOr') then 
  vp_max_c22 = maxval(sqrt(c22/rho))
  vp_max = max(vp_max, vp_max_c22)
endif

write(*,*) 'vp_max (m/s) = ', vp_max
if (dti_overrideflag) then
  dti = dti_override
  write(*,*) 'Using user-provided time step dti:'
  write(*,*) 'dti_override [s] = ', dti_override
endif

! ..Courant# for O(2) time leap frog SG FD
courant = 1./(sqrt(3.)*(c0 - c1 + c2 - c3 + c4 - c5 + c6 - c7))
courant = courant_safe*courant  ! drop Courant# by this safety factor
if (ordertime == 4) then
  courant = courant*0.8/0.54  ! O(4) time Courant# for Omelyan symplectic
endif
dti_max = courant*dl_min/vp_max
dti = dti_max  ! use time step determined by FD stability

! ..make sure that SEG-Y time sampling is greater than or equal to dti:
if (dti >= dti_segy) then
  write(*,*) 'SEG-Y time sampling is too fine!:'
  write(*,*) '..Smallest SEG-Y time sampling possible is [s]:', dti
  write(*,*) '..Terminating program execution!'
  stop
endif

nti = nint((tmax + tshift)/dti) + 1  ! add extra tshift to simul 
ntshift = nint(tshift/dti) 
nti_segy = nint(tmax/dti_segy) + 1  ! #time steps in SEG-Y's 
write(*,*) 'Max time (s) = ', tmax
write(*,*) 'Courant# = ', courant
write(*,*) 'FD time step [s] = ', dti
write(*,*) '# FD time steps = ', nti
write(*,*) 'SEG-Y time step (s) = ', dti_segy
write(*,*) '# SEG-Y time steps = ', nti_segy

if (rectype == 'p' .or. rectype == 'plz') then
  allocate(prec(nti+1,nrec))
else
  allocate(prec(nti+1,nrec), vxrec(nti+1,nrec), vyrec(nti+1,nrec), &
         & vzrec(nti+1,nrec))
endif

! ..initialize rec arrays to zero:
prec = 0  
if (rectype /= 'p' .and. rectype /= 'plz') then
  vxrec = 0  
  vyrec = 0 
  vzrec = 0
endif

! ..timesteps to output snapshots:
itsnp = tsnp/dti + 1

! ..compute symplectic time steps:
if (ordertime == 2) then  ! O(2) time standard
  dti_tsym(1) = xisym_o2*dti
  dti_tsym(2) = xisym_o2*dti

  dti_vsym(1) = lamsym_o2*dti


else  ! O(4) time Omelyan et al. (2003; Comp. Phys. Comm., 146, 188-202)
  dti_tsym(1) = xisym_o4*dti
  dti_tsym(2) = chisym_o4*dti
  dti_tsym(3) = (1 - 2*(chisym_o4 + xisym_o4))*dti
  dti_tsym(4) = chisym_o4*dti
  dti_tsym(5) = xisym_o4*dti

  dti_vsym(1) = (1 - 2*lamsym_o4)*dti*.5
  dti_vsym(2) = lamsym_o4*dti
  dti_vsym(3) = lamsym_o4*dti
  dti_vsym(4) = (1 - 2*lamsym_o4)*dti*.5
endif

!--PREP SOURCE WAVELET-----------------------------------------------------

write(*,*) '..resampling and integrating source wavelet'
nwav = nint(tmax/dti)  ! set #pts in wavelet to cover full simulation time
allocate (swav(nwav), swav_(nwav))
! ..resample sou wavelet to dti and perform time-integration:
call prepWavelet(tmax, fmax, dti_raw, dti, dti_segy, &
               & nwav_raw, nwav, swav_raw, swav, swav_)


!--SET UP (DECIMATED) SNAPSHOT OUTPUT--------------------------------------

! ..dimensions for snapshot output:
if (mod(nx,ijk_dec) == 0) then
  nx_dec = nx/ijk_dec
else
  nx_dec = nx/ijk_dec + 1
endif

if (mod(ny,ijk_dec) == 0) then
  ny_dec = ny/ijk_dec
else
  ny_dec = ny/ijk_dec + 1
endif

if (mod(nz,ijk_dec) == 0) then
  nz_dec = nz/ijk_dec
else
  nz_dec = nz/ijk_dec + 1
endif


!--SET UP BOUNDARY CONDITIONS----------------------------------------------

!--Set-up Maxwell viscoelastic absorbing boundaries:

allocate(deta(ilo:iup,jlo:jup,klo:kup))  ! for ABC
deta = 0  ! initialize array to  zero
write(*,*) '..generating absorbing decay function for Maxwell ABC'
call generateABCFunction(orderspace, freesurfflag, ilo, jlo, klo, &
        & iup, jup, kup, nz, nabc_top, nabc_bot, nabc_sdx, nabc_sdy, &
        & deta_maxtop, deta_maxbot, deta_maxsdx, deta_maxsdy, deta)

!--Allocate free-surface boundary arrays if required:

if (freesurfflag) then
  allocate (txxtop(ilo:iup,jlo:jup))
  allocate (tyytop(ilo:iup,jlo:jup))
endif


!--ZERO WAVEFIELDS---------------------------------------------------------

vx = 0.
vy = 0.
vz = 0.
txx = 0.
tyy = 0.
tzz = 0.
txz = 0.
tyz = 0.
txy = 0.
sx = 0.
sy = 0.
sz = 0.

!--OPEN MONITOR RECEIVER FILE----------------------------------------------

write(filenum3,'(i6.6)') shotnum 
filename6='monitor_S'//filenum3//'.dat'  
open(unit=7, file=filename6, status='unknown', action='write', &
   & iostat=status)
write(7,*) 'Monitor Data for FDAniQ.f90'
write(7,*) 'Version = ', ver
write(7,*) '... total number of times steps = ', nti
write(7,*) '... FD time step [s] = ', dti
write(7,*) '... Courant# = ', courant
write(7,*) '... Maxium Vp [m/s] = ', vp_max
write(7,*) '... order of FD operator in time = ', ordertime
write(7,*) '... number of threads = ', nthreads

write(7,*) 'it.# tot.time[s] int.time[s]  p(@mon.rx) &
          & vx(@mon.rx) vy(@mon.rx) vz(@mon.rx)'
t0 = dclock()  ! dclock returns elapsed time (use with Intel compiler)
t2 = t0

!--FINITE DIFFERENCE LOOP--------------------------------------------------

    filename3='xz_kurt_Vp'
    open(unit=3, file=filename3, status='unknown', action='write', &
         & iostat=status)
    do k=1,nz
      do i=1,nx
        write(3,*) i,k,-(txx(i,401,k)+tyy(i,401,k)+tzz(i,401,k))/3.0
      enddo
      write(3,*) ''
    enddo
    close(unit=3)

! I. Main FD Time Loop:
isnps = 1  ! initialize current snapshot# to 1
do kti = 1, nti    
  ti = dti*kti
  write(*,*) 'Executing time step# ', kti

  ! II. Symplectic Time Integration Subloop:
  do sti = 1, nsti

    ! II.A. Update Stress Fields:

    ! store old stresses on free-surface:
    if (freesurfflag) then
      do j=jmin, jmax
        do i=imin, imax
          txxtop(i,j) = txx(i,j,1)    ! set free-surface at k=1
          tyytop(i,j) = tyy(i,j,1)    
        enddo
      enddo
    endif
 
    dti_sym = dti_tsym(sti)    ! set symplectic time step

    ! update stresses with FD:
    select case (anitype)  ! eIso; eVTI; eVOr; eTTI; eTOr; aTTI; aTOr
      case ("eIso")
        call updateStressElaIso(orderspace, dti_sym, dx, dy, dz, &
                    & ilo, jlo, klo, iup, jup, kup, &
                    & c33, c55, &
                    & deta, vx, vy, vz, txx, tyy, tzz, tyz, txz, txy)
      case ("eVTI")
        call updateStressElaVTI(orderspace, dti_sym, dx, dy, dz, &
                    & ilo, jlo, klo, iup, jup, kup, &
                    & c11, c33, c55, c66, c13, &
                    & deta, vx, vy, vz, txx, tyy, tzz, tyz, txz, txy)
      case ("eVOr")
        call updateStressElaVOr(orderspace, dti_sym, dx, dy, dz, &
                    & ilo, jlo, klo, iup, jup, kup, &
                    & c11, c22, c33, c44, c55, c66, c12, c13, c23, &
                    & deta, vx, vy, vz, txx, tyy, tzz, tyz, txz, txy)
      case ("eTTI", "eTOr", "aTTI", "aTOr")
        call updateStressElaTri(orderspace, dti_sym, dx, dy, dz, &
                    & ilo, jlo, klo, iup, jup, kup, &
                    & c11, c22, c33, c44, c55, c66, &
                    & c12, c13, c14, c15, c16, &
                    & c23, c24, c25, c26, &
                    & c34, c35, c36, &
                    & c45, c46, &
                    & c56, &
                    & deta, vx, vy, vz, txx, tyy, tzz, tyz, txz, txy)
    end select

   !!!! DEBUG STUFF:
   !!!! write(*,*) 'loc A!'
   !!!! write(*,*) 'kti, sti = ', kti, sti
   !!!! write(*,*) 'maxval(pabs) = ', maxval(abs((txx+tyy+tzz)/3))
   !!!! loc = maxloc(abs((txx+tyy+tzz)/3))
   !!!! write(*,*) 'maxloc(pabs) = ', loc(1), loc(2), loc(3)
   !!!! write(*,*) 'rho(loc) = ', rho(loc(1),loc(2),loc(3))
   !!!! write(*,*) 'c11(loc) = ', c11(loc(1),loc(2),loc(3))
   !!!! write(*,*) 'c22(loc) = ', c22(loc(1),loc(2),loc(3))
   !!!! write(*,*) 'c33(loc) = ', c33(loc(1),loc(2),loc(3))
   !!!! write(*,*) 'c44(loc) = ', c44(loc(1),loc(2),loc(3))
   !!!! write(*,*) 'c55(loc) = ', c55(loc(1),loc(2),loc(3))
   !!!! write(*,*) 'c66(loc) = ', c66(loc(1),loc(2),loc(3))
   !!!! write(*,*) 'c12(loc) = ', c12(loc(1),loc(2),loc(3))
   !!!! write(*,*) 'c13(loc) = ', c13(loc(1),loc(2),loc(3))
   !!!! write(*,*) 'c23(loc) = ', c23(loc(1),loc(2),loc(3))
   !!!!

    ! apply pressure source:        
    write(*,*) 'kti=',kti,'nwav=',nwav,'nsou=',nsou
    if (kti <= nwav) then
      do l=1, nsou
        if (soutype(l) == 'p' .or. soutype(l) == 'plz') then  ! pressure sou
          if (soutype(l) == 'p') then
            ! pressure source: 
            ! ..nearest grid point: 
            icellmin = nint(xsou(l)/dx) + 1 ! to left of extrap pt
            jcellmin = nint(ysou(l)/dy) + 1 ! 
            kcellmin = nint(zsou(l)/dz) + 1 ! above extrap pt
            icellmax = icellmin
            jcellmax = jcellmin
            kcellmax = kcellmin
          elseif (soutype(l) == 'plz') then
            icellmin = 1
            jcellmin = 1
            kcellmin = nint(zsou(l)/dz) + 1
            icellmax = nx
            jcellmax = ny
            kcellmax = kcellmin
          endif

          ! ..fractional x, y and z distances from grid pt to sou:
          dx_frac = (xsou(l) - dx*float(icellmin - 1))/dx  
          dy_frac = (ysou(l) - dy*float(jcellmin - 1))/dy  
          dz_frac = (zsou(l) - dz*float(kcellmin - 1))/dz

          ! ..sinc extrapolation weights: 
          call generateSincWeights(dx_frac, dy_frac, dz_frac, fsinc)
       
          ! ..apply pressure sou on 8x8 grid (minus is for compression;
          !   use time-integ sou wavelet to preserve wavelet shape):
!$OMP PARALLEL DO PRIVATE(icell,jcell,kcell,i,j,k,ii,jj,kk)
          do jcell = jcellmin,jcellmax
            do icell = icellmin,icellmax
              do kcell = kcellmin,kcellmax

                do k = 1, 8
                  kk = kcell - 3 + (k-1)
                  do j = 1, 8
                    jj = jcell - 3 + (j-1)
                    do i = 1, 8
                      ii = icell - 3 + (i-1)
                      txx(ii,jj,kk) = txx(ii,jj,kk) - fsinc(i,j,k)*dti_sym* &
                                    & souamp(l,1)*swav_(kti) 
                      tyy(ii,jj,kk) = tyy(ii,jj,kk) - fsinc(i,j,k)*dti_sym* &
                                    & souamp(l,1)*swav_(kti) 
                      tzz(ii,jj,kk) = tzz(ii,jj,kk) - fsinc(i,j,k)*dti_sym* &
                                    & souamp(l,1)*swav_(kti) 
                    enddo
                  enddo
                enddo

                if (ghostsouflag) then
                  ! ghost (mirror) source has opposite sign & is located at
                  ! image location:
                  do k = 1, 8
                    kk = 5 - kcell + (k-1)  ! reset kk to above image plane
                    do j = 1, 8
                      jj = jcell - 3 + (j-1)
                      do i = 1, 8
                        ii = icell - 3 + (i-1)
                        txx(ii,jj,kk) = txx(ii,jj,kk) + fsinc(i,j,k)* &
                                      & dti_sym*souamp(l,1)*swav_(kti) 
                        tyy(ii,jj,kk) = tyy(ii,jj,kk) + fsinc(i,j,k)* &
                                      & dti_sym*souamp(l,1)*swav_(kti) 
                        tzz(ii,jj,kk) = tzz(ii,jj,kk) + fsinc(i,j,k)* &
                                      & dti_sym*souamp(l,1)*swav_(kti) 
                      enddo
                    enddo
                  enddo
                endif

                if (freesurfflag) then  ! add free-surface mirror contrib 
                  ! for sou near FS, need to (neg) mirror sinc weights
                  ! above the FS to mirror positions below the FS:
                  do j = 1, 8
                    jj = jcell - 3 + (j-1)
                    do i = 1, 8
                      ii = icell - 3 + (i-1)
                      if (kcell == 1) then
                        txx(ii,jj,2) = txx(ii,jj,2) - fsinc(i,j,3)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                        tyy(ii,jj,2) = tyy(ii,jj,2) - fsinc(i,j,3)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                        tzz(ii,jj,2) = tzz(ii,jj,2) - fsinc(i,j,3)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 

                        txx(ii,jj,3) = txx(ii,jj,3) - fsinc(i,j,2)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                        tyy(ii,jj,3) = tyy(ii,jj,3) - fsinc(i,j,2)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                        tzz(ii,jj,3) = tzz(ii,jj,3) - fsinc(i,j,2)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 

                        txx(ii,jj,4) = txx(ii,jj,4) - fsinc(i,j,1)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                        tyy(ii,jj,4) = tyy(ii,jj,4) - fsinc(i,j,1)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                        tzz(ii,jj,4) = tzz(ii,jj,4) - fsinc(i,j,1)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                      elseif (kcell == 2) then
                        txx(ii,jj,2) = txx(ii,jj,2) - fsinc(i,j,2)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                        tyy(ii,jj,2) = tyy(ii,jj,2) - fsinc(i,j,2)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                        tzz(ii,jj,2) = tzz(ii,jj,2) - fsinc(i,j,2)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 

                        txx(ii,jj,3) = txx(ii,jj,3) - fsinc(i,j,1)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                        tyy(ii,jj,3) = tyy(ii,jj,3) - fsinc(i,j,1)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                        tzz(ii,jj,3) = tzz(ii,jj,3) - fsinc(i,j,1)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                      elseif (kcell == 3) then
                        txx(ii,jj,2) = txx(ii,jj,2) - fsinc(i,j,1)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                        tyy(ii,jj,2) = tyy(ii,jj,2) - fsinc(i,j,1)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                        tzz(ii,jj,2) = tzz(ii,jj,2) - fsinc(i,j,1)*dti_sym* &
                                  & souamp(l,1)*swav_(kti) 
                      endif
                    enddo
                  enddo
                endif
              enddo
            enddo
          enddo
        endif
      enddo
    endif

    ! apply free-surface (for O(16) spatial FD & isotropic medium):
    if (freesurfflag) then
      do j=jmin, jmax
        do i=imin, imax
          tyz(i,j,-7) = -tyz(i,j,8)  ! negative mirror
          tyz(i,j,-6) = -tyz(i,j,7)  
          tyz(i,j,-5) = -tyz(i,j,6)  
          tyz(i,j,-4) = -tyz(i,j,5)  
          tyz(i,j,-3) = -tyz(i,j,4)  
          tyz(i,j,-2) = -tyz(i,j,3)
          tyz(i,j,-1) = -tyz(i,j,2)
          tyz(i,j,0)  = -tyz(i,j,1)
        
          txz(i,j,-7) = -txz(i,j,8)  ! negative mirror
          txz(i,j,-6) = -txz(i,j,7)  
          txz(i,j,-5) = -txz(i,j,6)  
          txz(i,j,-4) = -txz(i,j,5)  
          txz(i,j,-3) = -txz(i,j,4)  
          txz(i,j,-2) = -txz(i,j,3)
          txz(i,j,-1) = -txz(i,j,2)
          txz(i,j,0)  = -txz(i,j,1)

          txy(i,j,-6) = -txy(i,j,8)  ! negative mirror
          txy(i,j,-5) = -txy(i,j,7)  
          txy(i,j,-4) = -txy(i,j,6)  
          txy(i,j,-3) = -txy(i,j,5)  
          txy(i,j,-2) = -txy(i,j,4)  
          txy(i,j,-1) = -txy(i,j,3)
          txy(i,j,0)  = -txy(i,j,2)
          txy(i,j,1)  = 0.           ! zero at free-surface

          txx(i,j,-6) = txx(i,j,8)   ! positive mirror
          txx(i,j,-5) = txx(i,j,7)   
          txx(i,j,-4) = txx(i,j,6)   
          txx(i,j,-3) = txx(i,j,5)   
          txx(i,j,-2) = txx(i,j,4)   
          txx(i,j,-1) = txx(i,j,3)
          txx(i,j,0)  = txx(i,j,2)

          tzz(i,j,-6) = -tzz(i,j,8)  ! negative mirror
          tzz(i,j,-5) = -tzz(i,j,7)  
          tzz(i,j,-4) = -tzz(i,j,6)  
          tzz(i,j,-3) = -tzz(i,j,5)  
          tzz(i,j,-2) = -tzz(i,j,4)  
          tzz(i,j,-1) = -tzz(i,j,3)
          tzz(i,j,0)  = -tzz(i,j,2)
          tzz(i,j,1)  = 0.           ! zero at free-surface

          c13_ = c33(i,j,1) - 2*c55(i,j,1)

          dxvx = (c0*(vx(i+1,j,1) - vx(i,j,1))   + &
                & c1*(vx(i+2,j,1) - vx(i-1,j,1)) + &
                & c2*(vx(i+3,j,1) - vx(i-2,j,1)) + &
                & c3*(vx(i+4,j,1) - vx(i-3,j,1)) + &
                & c4*(vx(i+5,j,1) - vx(i-4,j,1)) + &
                & c5*(vx(i+6,j,1) - vx(i-5,j,1)) + &
                & c6*(vx(i+7,j,1) - vx(i-6,j,1)) + &
                & c7*(vx(i+8,j,1) - vx(i-7,j,1)) ) / dx

          dyvy = (c0*(vy(i,j,1)   - vy(i,j-1,1)) + &
                & c1*(vy(i,j+1,1) - vy(i,j-2,1)) + &
                & c2*(vy(i,j+2,1) - vy(i,j-3,1)) + &
                & c3*(vy(i,j+3,1) - vy(i,j-4,1)) + &
                & c4*(vy(i,j+4,1) - vy(i,j-5,1)) + &
                & c5*(vy(i,j+5,1) - vy(i,j-6,1)) + &
                & c6*(vy(i,j+6,1) - vy(i,j-7,1)) + &
                & c7*(vy(i,j+7,1) - vy(i,j-8,1)) ) / dy

          if (anitype =='eIso') then
             dum1 = c33(i,j,1) - (c13_**2)/c33(i,j,1)
             dum2 = c13_ - (c13_**2)/c33(i,j,1)
          else
             dum1 = c11(i,j,1) - (c13_**2)/c33(i,j,1)
             dum2 = c13_ - (c13_**2)/c33(i,j,1)
          endif
          txx(i,j,1) = txxtop(i,j) + dti_sym*(dum1*dxvx + dum2*dyvy)
          tyy(i,j,1) = tyytop(i,j) + dti_sym*(dum2*dxvx + dum1*dyvy)
        enddo
      enddo
    endif
    write(*,*) 'After freesurface mirror'


    ! Particle Velocity Update:
    if (sti < nsti) then
      dti_sym = dti_vsym(sti)  ! set symplectic time step
      write(*,*) 'Before particle velocity'
      call updateParticleVelocity(orderspace, dti_sym, &
                   & dx, dy, dz, ilo, jlo, klo, iup, jup, kup, rho, &
                   & deta, vx, vy, vz, txx, tyy, tzz, tyz, txz, txy, &
                   & itausig, difitau, sx, sy, sz)
      write(*,*) 'After particle velocity' 

      if (kti < nwav) then
        do l=1, nsou

          ! apply body force or particle velocity source:
          if (soutype(l) == 'f' .or. soutype(l) == 'v') then  

            ! fx/vx contribution:

            ! nearest grid point:
            icell = nint(xsou(l)/dx) + 1 ! to left of extrap pt
            jcell = nint(ysou(l)/dy) + 1 
            kcell = nint(zsou(l)/dz) + 1 ! above interp pt:

            ! ..fractional distances from grid pt to sou:
            ! (fx/vx sou needs to be shifted +0.5icell to colloc w/ pr)
            dx_frac = (xsou(l) + 0.5*dx - dx*float(icell - 1))/dx  
            dy_frac = (ysou(l) - dy*float(jcell - 1))/dy  
            dz_frac = (zsou(l) - dz*float(kcell - 1))/dz  
         
            ! ..sinc extrapolation weights: 
            call generateSincWeights(dx_frac, dy_frac, dz_frac, fsinc)
       
            ! ..apply fx/vx sou on 8x8 grid:
            do k = 1, 8
              kk = kcell - 3 + (k-1)
              do j = 1, 8
                jj = jcell - 3 + (j-1)
                do i = 1, 8
                  ii = icell - 3 + (i-1)
                  if (soutype(l) == 'f') then  ! fx sou
                    vx(ii,jj,kk) = vx(ii,jj,kk) + fsinc(i,j,k)* &
                              & dti_sym*souamp(l,1)*swav(kti)/rho(ii,jj,kk)  
                  else  ! vx sou
                    vx(ii,jj,kk) = vx(ii,jj,kk) + fsinc(i,j,k)* &
                              & dti_sym*souamp(l,1)*swav_(kti)  
                  endif
                enddo
              enddo
            enddo

            ! fy/vy contribution:

            ! nearest grid point:
            icell = nint(xsou(l)/dx) + 1 ! to left of extrap pt
            jcell = nint(ysou(l)/dy) + 1 
            kcell = nint(zsou(l)/dz) + 1 ! above interp pt:

            ! ..fractional distances from grid pt to sou:
            !    (fy/vy sou needs to be shifted -0.5dy to colloc w/ pr)
            dx_frac = (xsou(l) - dx*float(icell - 1))/dx  
            dy_frac = (ysou(l) - 0.5*dy - dy*float(jcell - 1))/dy  
            dz_frac = (zsou(l) - dz*float(kcell - 1))/dz  
         
            ! ..sinc extrapolation weights: 
            call generateSincWeights(dx_frac, dy_frac, dz_frac, fsinc)
       
            ! ..apply fy/vy sou on 8x8 grid:
            do k = 1, 8
              kk = kcell - 3 + (k-1)
              do j = 1, 8
                jj = jcell - 3 + (j-1)
                do i = 1, 8
                  ii = icell - 3 + (i-1)
                  if (soutype(l) == 'f') then  ! fx sou
                    vy(ii,jj,kk) = vy(ii,jj,kk) + fsinc(i,j,k)* &
                              & dti_sym*souamp(l,2)*swav(kti)/rho(ii,jj,kk)
                  else  ! vy sou
                    vy(ii,jj,kk) = vy(ii,jj,kk) + fsinc(i,j,k)* &
                              & dti_sym*souamp(l,2)*swav_(kti)  
                  endif
                enddo
              enddo
            enddo

            ! fz/vz contribution:

            ! nearest grid point:
            icell = nint(xsou(l)/dx) + 1 ! to left of extrap pt
            jcell = nint(ysou(l)/dy) + 1 
            kcell = nint(zsou(l)/dz) + 1 ! above interp pt:

            ! ..fractional distances from grid pt to sou:
            ! (fz/vz sou need to be shifted -0.5kcell to colloc w/ pr)
            dx_frac = (xsou(l) - dx*float(icell - 1))/dx  
            dy_frac = (ysou(l) - dy*float(jcell - 1))/dy  
            dz_frac = (zsou(l) - 0.5*dz - dz*float(kcell - 1))/dz  
         
            ! ..sinc extrapolation weights: 
            call generateSincWeights(dx_frac, dy_frac, dz_frac, fsinc)
       
            ! ..apply fz/vz sou on 8x8 grid:
            do k = 1, 8
              kk = kcell - 3 + (k-1)
              do j = 1, 8
                jj = jcell - 3 + (j-1)
                do i = 1, 8
                  ii = icell - 3 + (i-1)
                  if (soutype(l) == 'f') then  ! fz sou
                    vz(ii,jj,kk) = vz(ii,jj,kk) + fsinc(i,j,k)* &
                                 & dti_sym*souamp(l,3)*swav(kti)/ &
                                 & rho(ii,jj,kk)  
                  else  ! vz sou
                    vz(ii,jj,kk) = vz(ii,jj,kk) + fsinc(i,j,k)* &
                                 & dti_sym*souamp(l,3)*swav_(kti)  
                  endif
                enddo
              enddo
            enddo

            if (ghostsouflag) then

              ! fx/vx ghost (mirror) sou have same sign (for S-waves)
              ! & are located at image location:
              do k = 1, 8
                kk = 5 - kcell + (k-1)  ! reset kk to above image plane
                do j = 1, 8
                  jj = jcell - 3 + (j-1)
                  do i = 1, 8
                    ii = icell - 3 + (i-1)
                    if (soutype(l) == 'f') then ! fx sou ghost (pos mirror)
                      vx(ii,jj,kk) = vx(ii,jj,kk) + fsinc(i,j,k)* &
                                   & dti_sym*souamp(l,1)*swav(kti)/ &
                                   & rho(ii,jj,kk)  
                    else  ! vx sou ghost (pos mirror)
                      vx(ii,jj,kk) = vx(ii,jj,kk) + fsinc(i,j,k)* &
                                   & dti_sym*souamp(l,1)*swav_(kti)
                    endif
                  enddo
                enddo
              enddo

              ! fy/vy ghost (mirror) sou have same sign (for S-waves)
              ! & are located at image location:
              do k = 1, 8
                kk = 5 - kcell + (k-1)  ! reset kk to above image plane
                do j = 1, 8
                  jj = jcell - 3 + (j-1)
                  do i = 1, 8
                    ii = icell - 3 + (i-1)
                    if (soutype(l) == 'f') then ! fx sou ghost (pos mirror)
                      vy(ii,jj,kk) = vy(ii,jj,kk) + fsinc(i,j,k)* &
                                   & dti_sym*souamp(l,2)*swav(kti)/ &
                                   & rho(ii,jj,kk)  
                    else  ! vy sou ghost (pos mirror)
                      vy(ii,jj,kk) = vy(ii,jj,kk) + fsinc(i,j,k)* &
                                   & dti_sym*souamp(l,2)*swav_(kti)
                    endif
                  enddo
                enddo
              enddo

              ! fz/vz ghost (mirror) sou have same sign & are located at
              ! image location:
              do k = 1, 8
                kk = 4 - kcell + (k-1)  ! reset kk to above image plane
                do j = 1, 8
                  jj = jcell - 3 + (j-1)
                  do i = 1, 8
                    ii = icell - 3 + (i-1)
                    if (soutype(l) == 'f') then ! fz sou ghost (pos mirror)
                      vz(ii,jj,kk) = vz(ii,jj,kk) + fsinc(i,j,k)* &
                                   & dti_sym*souamp(l,3)*swav(kti)/ &
                                   & rho(ii,jj,kk)  
                    else  ! vz sou ghost (pos mirror)
                      vz(ii,jj,kk) = vz(ii,jj,kk) + fsinc(i,j,k)* &
                                  & dti_sym*souamp(l,3)*swav_(kti) 
                    endif
                  enddo
                enddo
              enddo

            endif

            if (freesurfflag) then  ! add free-surface mirror contrib 
              ! for sou near FS, need to (pos) mirror sinc weights
              ! above the FS to mirror positions below the FS:

              ! fx/vx FS mirror contribution:       
              do j = 1, 8
                jj = jcell - 3 + (j-1)
                do i = 1, 8
                  ii = icell - 3 + (i-1)
                  if (soutype(l) == 'f') then  ! fx sou FS mirror
                    if (kcell == 1) then
                      vx(ii,jj,2) = vx(ii,jj,2) + fsinc(i,j,3)* &
                              & dti_sym*souamp(l,1)*swav(kti)/rho(ii,jj,kk)
                      vx(ii,jj,3) = vx(ii,jj,3) + fsinc(i,j,2)* &
                              & dti_sym*souamp(l,1)*swav(kti)/rho(ii,jj,kk)
                      vx(ii,jj,4) = vx(ii,jj,4) + fsinc(i,j,1)* &
                              & dti_sym*souamp(l,1)*swav(kti)/rho(ii,jj,kk)
                    elseif (kcell == 2) then
                      vx(ii,jj,2) = vx(ii,jj,2) + fsinc(i,j,2)* &
                              & dti_sym*souamp(l,1)*swav(kti)/rho(ii,jj,kk)
                      vx(ii,jj,3) = vx(ii,jj,3) + fsinc(i,j,1)* &
                              & dti_sym*souamp(l,1)*swav(kti)/rho(ii,jj,kk)
                    elseif (kcell == 3) then
                      vx(ii,jj,2) = vx(ii,jj,2) + fsinc(i,j,1)* &
                              & dti_sym*souamp(l,1)*swav(kti)/rho(ii,jj,kk)
                    endif
                  else  ! vx sou FS mirror
                    if (kcell == 1) then
                      vx(ii,jj,2) = vx(ii,jj,2) + fsinc(i,j,3)* &
                                  & dti_sym*souamp(l,1)*swav_(kti)
                      vx(ii,jj,3) = vx(ii,jj,3) + fsinc(i,j,2)* &
                                  & dti_sym*souamp(l,1)*swav_(kti)
                      vx(ii,jj,4) = vx(ii,jj,4) + fsinc(i,j,1)* &
                                  & dti_sym*souamp(l,1)*swav_(kti)
                    elseif (kcell == 2) then
                      vx(ii,jj,2) = vx(ii,jj,2) + fsinc(i,j,2)* &
                                  & dti_sym*souamp(l,1)*swav_(kti)
                      vx(ii,jj,3) = vx(ii,jj,3) + fsinc(i,j,1)* &
                                  & dti_sym*souamp(l,1)*swav_(kti)
                    elseif (kcell == 3) then
                      vx(ii,jj,2) = vx(ii,jj,2) + fsinc(i,j,1)* &
                                  & dti_sym*souamp(l,1)*swav_(kti)
                    endif
                  endif
                enddo
              enddo

              ! fy/vy FS mirror contribution:       
              do j = 1, 8
                jj = jcell - 3 + (j-1)
                do i = 1, 8
                  ii = icell - 3 + (i-1)
                  if (soutype(l) == 'f') then  ! fx sou FS mirror
                    if (kcell == 1) then
                      vy(ii,jj,2) = vy(ii,jj,2) + fsinc(i,j,3)* &
                              & dti_sym*souamp(l,2)*swav(kti)/rho(ii,jj,kk)
                      vy(ii,jj,3) = vy(ii,jj,3) + fsinc(i,j,2)* &
                              & dti_sym*souamp(l,2)*swav(kti)/rho(ii,jj,kk)
                      vy(ii,jj,4) = vy(ii,jj,4) + fsinc(i,j,1)* &
                              & dti_sym*souamp(l,2)*swav(kti)/rho(ii,jj,kk)  
                    elseif (kcell == 2) then
                      vy(ii,jj,2) = vy(ii,jj,2) + fsinc(i,j,2)* &
                              & dti_sym*souamp(l,2)*swav(kti)/rho(ii,jj,kk)
                      vy(ii,jj,3) = vy(ii,jj,3) + fsinc(i,j,1)* &
                              & dti_sym*souamp(l,2)*swav(kti)/rho(ii,jj,kk)
                    elseif (kcell == 3) then
                      vy(ii,jj,2) = vy(ii,jj,2) + fsinc(i,j,1)* &
                              & dti_sym*souamp(l,2)*swav(kti)/rho(ii,jj,kk)
                    endif
                  else  ! vy sou FS mirror
                    if (kcell == 1) then
                      vy(ii,jj,2) = vy(ii,jj,2) + fsinc(i,j,3)* &
                                & dti_sym*souamp(l,2)*swav_(kti)
                      vy(ii,jj,3) = vy(ii,jj,3) + fsinc(i,j,2)* &
                                & dti_sym*souamp(l,2)*swav_(kti)
                      vy(ii,jj,4) = vy(ii,jj,4) + fsinc(i,j,1)* &
                                & dti_sym*souamp(l,2)*swav_(kti)
                    elseif (kcell == 2) then
                      vy(ii,jj,2) = vy(ii,jj,2) + fsinc(i,j,2)* &
                                & dti_sym*souamp(l,2)*swav_(kti)
                      vy(ii,jj,3) = vy(ii,jj,3) + fsinc(i,j,1)* &
                                & dti_sym*souamp(l,2)*swav_(kti)
                    elseif (kcell == 3) then
                      vy(ii,jj,2) = vy(ii,jj,2) + fsinc(i,j,1)* &
                                & dti_sym*souamp(l,2)*swav_(kti)
                    endif
                  endif
                enddo
              enddo

              ! fz/vz FS mirror contribution:
              do j = 1, 8
                jj = jcell - 3 + (j-1)
                do i = 1, 8
                  ii = icell - 3 + (i-1)
                  if (soutype(l) == 'f') then  ! fz sou FS mirror
                    if (kcell == 1) then
                      vz(ii,jj,1) = vz(ii,jj,1) + fsinc(i,j,3)* &
                              & dti_sym*souamp(l,3)*swav(kti)/rho(ii,jj,kk)
                      vz(ii,jj,2) = vz(ii,jj,2) + fsinc(i,j,2)* &
                              & dti_sym*souamp(l,3)*swav(kti)/rho(ii,jj,kk)
                      vz(ii,jj,3) = vz(ii,jj,3) + fsinc(i,j,1)* &
                              & dti_sym*souamp(l,3)*swav(kti)/rho(ii,jj,kk)
                    elseif (kcell == 2) then
                      vz(ii,jj,1) = vz(ii,jj,1) + fsinc(i,j,2)* &
                              & dti_sym*souamp(l,3)*swav(kti)/rho(ii,jj,kk)
                      vz(ii,jj,2) = vz(ii,jj,2) + fsinc(i,j,1)* &
                              & dti_sym*souamp(l,3)*swav(kti)/rho(ii,jj,kk)
                    elseif (kcell == 3) then
                      vz(ii,jj,1) = vz(ii,jj,1) + fsinc(i,j,1)* &
                              & dti_sym*souamp(l,3)*swav(kti)/rho(ii,jj,kk)
                    endif
                  else  ! vz sou FS mirror
                    if (kcell == 1) then
                      vz(ii,jj,1) = vz(ii,jj,1) + fsinc(i,j,3)* &
                                  & dti_sym*souamp(l,3)*swav_(kti)
                      vz(ii,jj,2) = vz(ii,jj,2) + fsinc(i,j,2)* &
                                  & dti_sym*souamp(l,3)*swav_(kti)
                      vz(ii,jj,3) = vz(ii,jj,3) + fsinc(i,j,1)* &
                                  & dti_sym*souamp(l,3)*swav_(kti)
                    elseif (kcell == 2) then
                      vz(ii,jj,1) = vz(ii,jj,1) + fsinc(i,j,2)* &
                                  & dti_sym*souamp(l,3)*swav_(kti)
                      vz(ii,jj,2) = vz(ii,jj,2) + fsinc(i,j,1)* &
                                  & dti_sym*souamp(l,3)*swav_(kti)
                    elseif (kcell == 3) then
                      vz(ii,jj,1) = vz(ii,jj,1) + fsinc(i,j,1)* &
                                  & dti_sym*souamp(l,3)*swav_(kti)
                    endif
                  endif
                enddo
              enddo

            endif  ! end free-surface mirror

          endif  
        enddo  
      endif  

      !..apply free-surface (for O(16) spatial FD):
      if (freesurfflag) then
        do j=jmin, jmax
          do i=imin, imax
            vx(i,j,-6) = vx(i,j,8)  ! positive mirror
            vx(i,j,-5) = vx(i,j,7)  
            vx(i,j,-4) = vx(i,j,6)  
            vx(i,j,-3) = vx(i,j,5)  
            vx(i,j,-2) = vx(i,j,4)  
            vx(i,j,-1) = vx(i,j,3)
            vx(i,j,0)  = vx(i,j,2)

            vy(i,j,-6) = vy(i,j,8)  ! positive mirror
            vy(i,j,-5) = vy(i,j,7)  
            vy(i,j,-4) = vy(i,j,6)  
            vy(i,j,-3) = vy(i,j,5)  
            vy(i,j,-2) = vy(i,j,4)  
            vy(i,j,-1) = vy(i,j,3)
            vy(i,j,0)  = vy(i,j,2)

            vz(i,j,-7) = vz(i,j,8)  ! positive mirror
            vz(i,j,-6) = vz(i,j,7)  
            vz(i,j,-5) = vz(i,j,6)  
            vz(i,j,-4) = vz(i,j,5)  
            vz(i,j,-3) = vz(i,j,4)  
            vz(i,j,-2) = vz(i,j,3)
            vz(i,j,-1) = vz(i,j,2)
            vz(i,j,0)  = vz(i,j,1)

          enddo
        enddo
      endif

    endif ! End of Velocity Update

    if (MOD(kti,24) == 1) then
       write(filenum4,'(i5.5,a,i1.1)') kti,'_',sti
       write(*,*) filenum4

       filename3='xz_kurt_'//filenum4//'_P'
       open(unit=3, file=filename3, status='unknown', action='write', &
            & iostat=status)
       do k=1,nz
         do i=1,nx
           write(3,*) i,k,-(txx(i,401,k)+tyy(i,401,k)+tzz(i,401,k))/3.0
         enddo
         write(3,*) ''
       enddo
       close(unit=3)
    endif

  enddo ! End of Symplectic Time Integration Subloop:

! Write Time Snapshots of 3-D Particle Velocities & Pressure Wavefields:
  if (isnps <= nsnps) then
    if (kti == itsnp(isnps)) then
      write(filenum2,'(i5.5)') kti
      write(filenum3,'(i6.6)') shotnum

      if (snpflag == 'g') then    ! gocad voxet flat binary format

        ! write gocad voxet flat binary file header
        filename3='snp_T'//filenum2//'_S'//filenum3//'.vo'    
        open(unit=2, file=filename3, status='unknown', action='write', &
           & iostat=status)
        write(2,*) 'GOCAD Voxet 0.01'
        write(2,*) 'HEADER {'
        filename3='snp_T'//filenum2//'_S'//filenum3   
        write(2,*) 'name: ', filename3 
        write(2,*) '}'
        write(2,*) 'AXIS_U  0.  0. -1.'   ! U is slow axis
        write(2,*) 'AXIS_V  0. -1.  0.'   ! V is intermed axis
        write(2,*) 'AXIS_W -1.  0.  0.'   ! W is fast axis
        write(2,'(a7,a1,i6,a1,i6,a1,i6)') ' AXIS_N', tab, nz_dec, &
              & tab, ny_dec, tab, nx_dec  ! number of points
        write(2,*) 'AXIS_O 0. 0. 0.'   ! coordinate origin
        write(2,'(a7,a1,f12.4,a1,f12.4,a1,f12.4)') ' AXIS_D', tab, &
              & dz*float(nz), tab, dy*float(ny), tab, &
              & dx*float(nx)  ! grid spacing
        write(2,*) 'AXIS_NAME "Z" "Y" "X"'  ! axis labels
        write(2,*) 'AXIS_UNIT "m" "m" "m"'  ! axis units

        filename4='snp_T'//filenum2//'_S'//filenum3//'_pr'  
        write(2,*) '    '
        write(2,'(a54,i5,a2,g10.4,a1)') ' PROPERTY 1 "pressure snapshot &
            & at time step#, time[s]:', kti, ';', ti,'"'
        write(2,*) 'PROP_UNIT 1 "(Pa)"'
        write(2,*) 'PROP_FILE 1 ', filename4

        filename4='snp_T'//filenum2//'_S'//filenum3//'_vx'  
        write(2,*) '    '
        write(2,'(a71,i5,a2,g10.4,a1)') ' PROPERTY 2 "horiz x particle &
            & velocity snapshot at time step#, time[s]:', kti, ';', ti,'"'
        write(2,*) 'PROP_UNIT 2 "(m/s)"'
        write(2,*) 'PROP_FILE 2 ', filename4

        filename4='snp_T'//filenum2//'_S'//filenum3//'_vy'  
        write(2,*) '    '
        write(2,'(a71,i5,a2,g10.4,a1)') ' PROPERTY 3 "horiz y particle &
            & velocity snapshot at time step#, time[s]:', kti, ';', ti,'"'
        write(2,*) 'PROP_UNIT 3 "(m/s)"'
        write(2,*) 'PROP_FILE 3 ', filename4

        filename4='snp_T'//filenum2//'_S'//filenum3//'_vz'  
        write(2,*) '    '
        write(2,'(a70,i5,a2,g10.4,a1)') ' PROPERTY 4 "vert z particle &
            & velocity snapshot at time step#, time[s]:', kti, ';', ti,'"'
        write(2,*) 'PROP_UNIT 4 "(m/s)"'
        write(2,*) 'PROP_FILE 4 ', filename4

        close(unit=2)

        ! write pr snapshot flat binary file
        allocate(ptemp3(ilo:iup,jlo:jup,klo:kup))
        ptemp3 = -(txx + tyy + tzz)/3
        filename4='snp_T'//filenum2//'_S'//filenum3//'_pr'   
        call writeGocadSnapShot(ijk_dec, nx_dec, ny_dec, nz_dec, &
                              & nx, ny, nz, ilo, jlo, klo, iup, jup, kup, &
                              & ptemp3, filename4)
        deallocate(ptemp3)

        ! write vx snapshot flat binary file
        filename4='snp_T'//filenum2//'_S'//filenum3//'_vx'   
        call writeGocadSnapShot(ijk_dec, nx_dec, ny_dec, nz_dec, &
                              & nx, ny, nz, ilo, jlo, klo, iup, jup, kup, &
                              & vx, filename4)

        ! write vy snapshot flat binary file
        filename4='snp_T'//filenum2//'_S'//filenum3//'_vy'   
        call writeGocadSnapShot(ijk_dec, nx_dec, ny_dec, nz_dec, &
                              & nx, ny, nz, ilo, jlo, klo, iup, jup, kup, &
                              & vy, filename4)

        ! write vz snapshot flat binary file
        filename4='snp_T'//filenum2//'_S'//filenum3//'_vz'   
        call writeGocadSnapShot(ijk_dec, nx_dec, ny_dec, nz_dec, &
                              & nx, ny, nz, ilo, jlo, klo, iup, jup, kup, &
                              & vz, filename4)

      else    ! tecplot360 ascii format

        filename2='snp_T'//filenum2//'_S'//filenum3//'.dat'
        open(unit=2, file=filename2, status='unknown', action='write', &
           & iostat=status)
        write(2,*) 'TITLE="time step#; time[s]',kti,';',ti,'"'
        write(2,*) 'VARIABLES="x[m]","y[m]","z[m]","p","vx","vy","vz"'
        write(2,*) 'ZONE T="snapshot",I=',nx_dec,',J=',ny_dec, &
                 & ',K=',nz_dec,',F=POINT'
        do k = 1, nz, ijk_dec
          z = dz*float(k-1)
          do j = 1, ny, ijk_dec
            y = dy*float(j-1)
            do i = 1, nx, ijk_dec
              x = dx*float(i-1)
              ! minus for compression:
              p = -(txx(i,j,k)+tyy(i,j,k)+tzz(i,j,k))/3  
              ! note: minus sign on z & vz for tecplot360 c.s.
              write(2,*) x, tab, y, tab, -z, tab, p, tab, vx(i,j,k), &
                       & tab, vy(i,j,k), tab, -vz(i,j,k)
            enddo                                                          
          enddo                                                          
        enddo
        close(unit=2)
      endif
      isnps = isnps + 1
    endif
  endif 

! Write Monitor Receiver Data to File:
  if (mod(real(kti), real(ndti_mon)) == 0) then
    t1 = t2  ! store prior elapsed time
    t2 = dclock()
    write(7,'(i5,2x,6(g10.4,2x))') kti, (t2-t0), (t2-t1), &
        & (-(txx(imon,jmon,kmon)+tyy(imon,jmon,kmon)+ &
        & tzz(imon,jmon,kmon))/3), vx(imon,jmon,kmon), &
        & vy(imon,jmon,kmon), vz(imon,jmon,kmon)
  endif

! Store Receiver Traces (for eventual write to SEG-Y files):
  if (kti > ntshift) then  ! start storing trace data 
    kti_shift = kti - ntshift + 1  ! +1 so 1st sample is at t = 0  

    ! Interpolate off-grid receivers: 
    do lr = 1, nrecline
      do irectype = 1, 4  ! loop over rec type
        if ((rectype == 'p' .or. rectype == 'plz') .and. irectype > 1) cycle

        if (rectype == 'p') then       
          ! fractional x,y,z distances from grid pt to rec:
          icellmin = nint(xrec(lr)/dx) + 1 
          jcellmin = nint(yrec(lr)/dy) + 1 
          kcellmin = nint(zrec(lr)/dz) + 1
          icellmax = icellmin
          jcellmax = jcellmin
          kcellmax = kcellmin
        elseif (rectype == 'plz') then
          icellmin = 1
          jcellmin = 1
          kcellmin = nint(zrec(llr)/dz) + 1
          icellmax = nx
          jcellmax = ny
          kcellmax = kcellmin
        endif

        ! proper fractional distances for pr rec (irectype == 1):
        dx_frac = (xrec(lr) - dx*float(icellmin - 1))/dx
        dy_frac = (yrec(lr) - dy*float(jcellmin - 1))/dy  
        dz_frac = (zrec(lr) - dz*float(kcellmin - 1))/dz  

        if (irectype == 2) then  ! vx rec
          ! locate nearest cell, with lower-bound (i,j,k) indices:
          ! note: vx rec need to be shifted  +0.5icell to colloc w/ pr
          dx_frac = (xrec(lr) + 0.5*dx - dx*float(icellmin - 1))/dx
        elseif (irectype == 3) then  ! vy rec
          ! note: vy rec need to be shifted -0.5jcell to colloc w/ pr
          dy_frac = (yrec(lr) - 0.5*dy - dy*float(jcellmin - 1))/dy  
        elseif (irectype == 4) then  ! vz rec
          ! note: vz rec need to be shifted -0.5kcell to colloc w/ pr
          dz_frac = (zrec(lr) - 0.5*dz - dz*float(kcellmin - 1))/dz  
        endif

        ! compute 8x8 sinc function interpolation weights: 
        call generateSincWeights(dx_frac, dy_frac, dz_frac, fsinc) 

        if (rectype == 'plz') then
           lro = (lr-1)*nx*ny+1;
        else
           lro = lr;
        endif

!$OMP PARALLEL DO PRIVATE(icell,jcell,kcell,i,j,k,ii,jj,kk,f_in,f_out,llr)
        do jcell = jcellmin,jcellmax
          do icell = icellmin,icellmax
            do kcell = kcellmin,kcellmax

              llr = lro + (icell-icellmin) + ((jcell-jcellmin) + (kcell-kcellmin)*(jcellmax-jcellmin+1))*(icellmax-icellmin+1)

              f_out = 0
              do k = 1, 8
                kk = kcell - 3 + (k-1)
                do j = 1, 8
                  jj = jcell - 3 + (j-1)
                  do i = 1, 8
                    ii = icell - 3 + (i-1)

                    if (irectype == 1) then
                      ! note: minus is for compression
                      f_in = -(txx(ii,jj,kk)+tyy(ii,jj,kk)+tzz(ii,jj,kk)) / 3
                    elseif (irectype == 2) then
                      f_in = vx(ii,jj,kk)
                    elseif (irectype == 3) then
                      f_in = vy(ii,jj,kk)
                    elseif (irectype == 4) then
                      f_in = vz(ii,jj,kk)
                    endif

                    ! use sinc weights to weight f_in contribution & then stack:
                    if (freesurfflag .and. kk < 1 .or. &
                      & ghostrecflag .and. kk < 1) then
                      f_out = f_out  ! don't add fields below FS/GS
                    else
                      f_out = f_out + fsinc(i,j,k)* f_in
                    endif
                  enddo
                enddo
              enddo

              if (irectype == 1) then
                prec(kti_shift,llr) = f_out
              elseif (irectype == 2) then
                vxrec(kti_shift,llr) = f_out
              elseif (irectype == 3) then
                vyrec(kti_shift,llr) = f_out
              elseif (irectype == 4) then
                vzrec(kti_shift,llr) = f_out
              endif

              ! add rec ghost (using sinc weights computed above):
              if (ghostrecflag) then  
                f_out = 0
                do k = 1, 8
                  do j = 1, 8
                    jj = jcell - 3 + (j-1)
                    do i = 1, 8
                      ii = icell - 3 + (i-1)

                      if (irectype == 1) then
                        kk = 5 - kcell + (k-1)  ! reset kk to above image plane
                        ! pressure (minus is for compression):
                        f_in = -(txx(ii,jj,kk) + tyy(ii,jj,kk) + &
                               & tzz(ii,jj,kk)) / 3
                      elseif (irectype == 2) then
                        kk = 5 - kcell + (k-1)  ! reset kk to above image plane
                        f_in = vx(ii,jj,kk)
                      elseif (irectype == 3) then
                        kk = 5 - kcell + (k-1)  ! reset kk to above image plane
                        f_in = vy(ii,jj,kk)
                      elseif (irectype == 4) then
                        kk = 4 - kcell + (k-1) ! reset kk to above image plane
                        f_in = vz(ii,jj,kk)
                      endif

                      ! use sinc weights to weigh f_in contrib & then stack:
                      if (irectype /= 4 .and. kk > 1 .or. &
                        & irectype == 4 .and. kk > 0) then
                        f_out = f_out  ! don't add fields below GS
                      else
                        f_out = f_out + fsinc(i,j,k)* f_in
                      endif
                    enddo
                  enddo
                enddo

                if (irectype == 1) then
                  ! neg mirror on pr
                  prec(kti_shift,llr) = prec(kti_shift,llr) - f_out  
                elseif (irectype == 2) then
                  ! neg mirror on vx
                  vxrec(kti_shift,llr) = vxrec(kti_shift,llr) - f_out  
                elseif (irectype == 3) then
                  ! neg mirror on vy (!!!!check this!)
                  vyrec(kti_shift,llr) = vyrec(kti_shift,llr) - f_out  
                elseif (irectype == 4) then
                  ! pos mirror on vz
                  vzrec(kti_shift,llr) = vzrec(kti_shift,llr) + f_out  
                endif
              endif

            enddo
          enddo
        enddo

      enddo
    enddo
  endif  

enddo  ! End of Main FD Time Loop

close(unit=7)  ! close Monitor receiver file

!--INTERPOLATE SHOT GATHER TRACES FROM FD TO SEG-Y TIME SAMPLE-------------
allocate(traces(nti_segy))  ! temp storage vector  

! Pressure (pr):
do ir = 1, nrec  

  do kti_segy = 1, nti_segy  
    t_segy = dti_segy*float(kti_segy - 1)  ! SEG-Y time [s]
    it_left = nint(t_segy/dti)  ! nearest sample to left of interp pt

    if (it_left >= 2 .and. it_left <= (nti - 2)) then  ! sinc interp
      data4_in(1) = prec(it_left-1,ir)
      data4_in(2) = prec(it_left,ir)
      data4_in(3) = prec(it_left+1,ir)
      data4_in(4) = prec(it_left+2,ir)
      call interpolateTrace(it_left, t_segy, dti, dti_segy, & 
                          & data4_in, data_out)
    elseif (it_left == 0) then  ! first sample in trace is for t = 0
      data_out = 0.

    else  ! linear interp
      ! normalized time from left sample to interp pt
      dt_frac = (t_segy - dti*float(it_left))/dti  

      data_out = (1 - dt_frac)*prec(it_left,ir) + &
                    &  dt_frac*prec(it_left+1,ir)
    endif
    traces(kti_segy) = data_out  ! store temporarily in traces vector

  enddo

  do kti_segy = 1, nti_segy  
    ! put interp trace back in rec array
    prec(kti_segy,ir) = traces(kti_segy)  
  enddo

enddo

if (rectype /= 'p' .and. rectype /= 'plz') then
  ! Horizontal Particle Velocity (vxr):
  do ir = 1, nrec
 
    do kti_segy = 1, nti_segy  
      t_segy = dti_segy*float(kti_segy - 1)  ! SEG-Y time [s]
      it_left = nint(t_segy/dti)  ! nearest samp to left of interp pt

      if (it_left >= 2 .and. it_left <= (nti - 2)) then ! sinc interp
        data4_in(1) = vxrec(it_left-1,ir)
        data4_in(2) = vxrec(it_left,ir)
        data4_in(3) = vxrec(it_left+1,ir)
        data4_in(4) = vxrec(it_left+2,ir)
        call interpolateTrace(it_left, t_segy, dti, dti_segy, &
                            & data4_in, data_out)

      elseif (it_left == 0) then  ! first sample in trace is for t = 0
        data_out = 0.

      else  ! linear interp
        ! normalized time from left sample to interp pt
        dt_frac = (t_segy - dti*float(it_left))/dti  

        data_out = (1 - dt_frac)*vxrec(it_left,ir) + &
                    &  dt_frac*vxrec(it_left+1,ir)
      endif
      traces(kti_segy) = data_out  ! store temp in traces vector

    enddo

    do kti_segy = 1, nti_segy  
      ! put interp trace back in rec array
      vxrec(kti_segy,ir) = traces(kti_segy)  
    enddo

  enddo

  ! Vertical Particle Velocity (vzr):
  do ir = 1, nrec
  
    do kti_segy = 1, nti_segy  
      t_segy = dti_segy*float(kti_segy - 1)  ! SEG-Y time [s]
      it_left = nint(t_segy/dti)  ! nearest samp to left of interp pt

      if (it_left >= 2 .and. it_left <= (nti - 2)) then ! sinc interp
        data4_in(1) = vzrec(it_left-1,ir)
        data4_in(2) = vzrec(it_left,ir)
        data4_in(3) = vzrec(it_left+1,ir)
        data4_in(4) = vzrec(it_left+2,ir)
        call interpolateTrace(it_left, t_segy, dti, dti_segy, &
                            & data4_in, data_out)

      elseif (it_left == 0) then  ! first sample in trace is for t = 0
        data_out = 0.

      else  ! linear interp
        ! normalized time from left sample to interp pt
        dt_frac = (t_segy - dti*float(it_left))/dti  

        data_out = (1 - dt_frac)*vzrec(it_left,ir) + &
                    &  dt_frac*vzrec(it_left+1,ir)
      endif
      traces(kti_segy) = data_out  ! store temp in traces vector

    enddo

    do kti_segy = 1, nti_segy  
      ! put interp trace back in rec array
      vzrec(kti_segy,ir) = traces(kti_segy)  
    enddo

  enddo

endif
deallocate(traces) 


! Write Traces to SEG-Y:
xsou_ = xsou(1) ! assume shot loc. coincides w/ sou#1 
ysou_ = ysou(1) 
zsou_ = zsou(1) 
irec_offset = 0

write(filenum3,'(i6.6)') shotnum

allocate(traces(nti_segy*nrec))   

! pressure (pr):
do ir = 1, nrec  ! re-package trace array as a vector for 
                 !  interlanguage call
  do kti_segy = 1, nti_segy  
    ktir = kti_segy + (ir-1)*nti_segy
    traces(ktir) = prec(kti_segy,ir)
  enddo
enddo
filename1='pr_S'//filenum3//'.segy' // char(0)
call writeSEGY(traces, iline, xline, shotnum, filename1, &
& xsou_, ysou_, zsou_, xrec, yrec, -zrec, nrec, nti_segy, &
& dti_segy, irec_offset,&
& O(1),O(2),O(3),U(1),U(2),U(3),V(1),V(2),V(3),W(1),W(2),W(3),&
& orderaxes,nx,ny,nz,dx,dy,dz)

if (rectype /= 'p' .and. rectype /= 'plz') then
  ! horizontal partical velocity (vx):
  do ir = 1, nrec  ! re-package trace array as a vector for 
                   !  interlanguage call
    do kti_segy = 1, nti_segy 
      ktir = kti_segy + (ir-1)*nti_segy
      traces(ktir) = vxrec(kti_segy,ir)
    enddo
  enddo
  filename1='vx_S'//filenum3//'.segy' // char(0)
  call writeSEGY(traces, iline, xline, shotnum, filename1, &
    & xsou_, ysou_, zsou_, xrec, yrec, -zrec, nrec, nti_segy, &
    & dti_segy, irec_offset,&
    & O(1),O(2),O(3),U(1),U(2),U(3),V(1),V(2),V(3),W(1),W(2),W(3),&
    & orderaxes,nx,ny,nz,dx,dy,dz)

  ! horizontal partical velocity (vy):
  do ir = 1, nrec  ! re-package trace array as a vector for 
                   !  interlanguage call
    do kti_segy = 1, nti_segy 
      ktir = kti_segy + (ir-1)*nti_segy
      traces(ktir) = vyrec(kti_segy,ir)
    enddo
  enddo
  filename1='vy_S'//filenum3//'.segy' // char(0)
  call writeSEGY(traces, iline, xline, shotnum, filename1, &
    & xsou_, ysou_, zsou_, xrec, yrec, -zrec, nrec, nti_segy, &
    & dti_segy, irec_offset,&
    & O(1),O(2),O(3),U(1),U(2),U(3),V(1),V(2),V(3),W(1),W(2),W(3),&
    & orderaxes,nx,ny,nz,dx,dy,dz)

  ! vertical partical velocity (vz):
  do ir = 1, nrec  ! re-package trace array as a vector for 
                   !  interlanguage call
    do kti_segy = 1, nti_segy  
      ktir = kti_segy + (ir-1)*nti_segy
      traces(ktir) = vzrec(kti_segy,ir)
    enddo
  enddo
  filename1='vz_S'//filenum3//'.segy' // char(0)
  call writeSEGY(traces, iline, xline, shotnum, filename1, &
    & xsou_, ysou_, zsou_, xrec, yrec, -zrec, nrec, nti_segy, &
    & dti_segy, irec_offset,&
    & O(1),O(2),O(3),U(1),U(2),U(3),V(1),V(2),V(3),W(1),W(2),W(3),&
    & orderaxes,nx,ny,nz,dx,dy,dz)

endif

deallocate(traces)
      

stop
end program
