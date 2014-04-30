subroutine generateABCFunction(orderspace, fsflag, ilo, jlo, klo, &
             & iup, jup, kup, nz, nabc_top, nabc_bot, nabc_sdx, nabc_sdy, &
             & deta_maxtop, deta_maxbot, deta_maxsdx, deta_maxsdy, deta)

!* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
! Computes the Maxwell Damping Coefficient Spatial Decay Function.
! Spatial decay of Maxwell damping factor is taken from the PML literature
! (Collino & Tsogka, 2001, Geophys., 66, 294-307).
!
! o Last Modified:  
!   10-08-12  Adapted from abc2d.f90.

! -.--.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

implicit none

logical :: fsflag
integer :: orderspace
integer :: nabc_top, nabc_bot, nabc_sdx, nabc_sdy
integer :: ilo, jlo, klo
integer :: iup, jup, kup
integer :: nz
real :: deta_maxtop, deta_maxbot
real, dimension(nz) :: deta_maxsdx, deta_maxsdy
real, dimension(ilo:iup,jlo:jup,klo:kup) :: deta  

integer :: imin, jmin, kmin
integer :: imax, jmax, kmax
integer :: icnt, jcnt, kcnt
integer :: i, j, k


! Set range of deta:
imin = ilo + orderspace/2
jmin = jlo + orderspace/2
kmin = klo + orderspace/2

imax = iup - orderspace/2
jmax = jup - orderspace/2
kmax = kup - orderspace/2

! Initialize Maxwell damping coefficient:
deta = 0.  ! initialize to no damping 

! Maxwell 1D Damping Function for Plates:
! (note: won't worry about beginning and ends of non-1D indices for plates 
!  since these will be overwritten subsequently by beams and cubes)

! ..top (-k) plate
if (.not. fsflag) then  
  kcnt = 0
  do k = kmin, kmin + nabc_top - 1
    kcnt = kcnt + 1
    do j = jmin, jmax  
      do i = imin, imax
        deta(i,j,k) = deta_maxtop*(float(nabc_top - kcnt)/nabc_top)**2
      enddo 
    enddo 
  enddo
endif

! ..bottom (+k) plate
kcnt = nabc_bot + 1
do k = kmax - nabc_bot + 1, kmax
  kcnt = kcnt - 1
  do j = jmin, jmax
    do i = imin, imax
      deta(i,j,k) = deta_maxbot*(float(nabc_bot - kcnt)/nabc_bot)**2
    enddo 
  enddo 
enddo

! ..west (-i) plate
do k = 1, nz
  do j = jmin, jmax
    icnt = 0
    do i = imin, imin + nabc_sdx - 1
      icnt = icnt + 1
      deta(i,j,k) = deta_maxsdx(k)*(float(nabc_sdx - icnt)/nabc_sdx)**2
    enddo 
  enddo 
enddo

! ..east (+i) plate
do k = 1, nz
  do j = jmin, jmax
    icnt = nabc_sdx + 1
    do i = imax - nabc_sdx + 1, imax
      icnt = icnt - 1
      deta(i,j,k) = deta_maxsdx(k)*(float(nabc_sdx - icnt)/nabc_sdx)**2
    enddo 
  enddo 
enddo

! ..south (-j) plate
do k = 1, nz
  jcnt = 0
  do j = jmin, jmin + nabc_sdy - 1
    jcnt = jcnt + 1
    do i = imin, imax
      deta(i,j,k) = deta_maxsdy(k)*(float(nabc_sdy - jcnt)/nabc_sdy)**2
    enddo 
  enddo 
enddo

! ..north (+j) plate
do k = 1, nz
  jcnt = nabc_sdy + 1
  do j = jmax - nabc_sdy + 1, jmax
    jcnt = jcnt - 1
    do i = imin, imax
      deta(i,j,k) = deta_maxsdy(k)*(float(nabc_sdy - jcnt)/nabc_sdy)**2
    enddo 
  enddo 
enddo


! Maxwell 2D Damping Function for Beams:
! (note: won't worry about beginning and ends of non-2D indices for beams 
!  since these will be overwritten subsequently by cubes)

! ..Top Beams
if (.not. fsflag) then  

! ....top-west beam
  kcnt = 0
  do k = kmin, kmin + nabc_top - 1
    kcnt = kcnt + 1
    do j = jmin, jmax
      icnt = 0
      do i = imin, imin + nabc_sdx - 1
        icnt = icnt + 1
        deta(i,j,k) = sqrt((deta_maxtop**2)*( &
             & (float(nabc_sdx - icnt)/nabc_sdx)**4 + &
             & (float(nabc_top - kcnt)/nabc_top)**4))
      enddo
    enddo
  enddo

! ....top-east beam
  kcnt = 0
  do k = kmin, kmin + nabc_top - 1
    kcnt = kcnt + 1
    do j = jmin, jmax
      icnt = nabc_sdx + 1
      do i = imax - nabc_sdx + 1, imax
        icnt = icnt - 1
        deta(i,j,k) = sqrt((deta_maxtop**2)*( &
             & (float(nabc_sdx - icnt)/nabc_sdx)**4 + &
             & (float(nabc_top - kcnt)/nabc_top)**4))
      enddo
    enddo
  enddo

! ....top-south beam
  kcnt = 0
  do k = kmin, kmin + nabc_top - 1
    kcnt = kcnt + 1
    jcnt = 0
    do j = jmin, jmin + nabc_sdy - 1
      jcnt = jcnt + 1
      do i = imin, imax
        deta(i,j,k) = sqrt((deta_maxtop**2)*( &
             & (float(nabc_sdy - jcnt)/nabc_sdy)**4 + &
             & (float(nabc_top - kcnt)/nabc_top)**4))
      enddo
    enddo
  enddo

! ....top-north beam
  kcnt = 0
  do k = kmin, kmin + nabc_top - 1
    kcnt = kcnt + 1
    jcnt = nabc_sdy + 1
    do j = jmax - nabc_sdy + 1, jmax
      jcnt = jcnt - 1
      do i = imin, imax
        deta(i,j,k) = sqrt((deta_maxtop**2)*( &
             & (float(nabc_sdy - jcnt)/nabc_sdy)**4 + &
             & (float(nabc_top - kcnt)/nabc_top)**4))
      enddo
    enddo
  enddo

endif

! ..Bottom Beams

! ....bot-west beam
  kcnt = nabc_bot + 1
  do k = kmax - nabc_bot + 1, kmax
    kcnt = kcnt - 1
    do j = jmin, jmax
      icnt = 0
      do i = imin, imin + nabc_sdx - 1
        icnt = icnt + 1
        deta(i,j,k) = sqrt((deta_maxbot**2)*( &
              & (float(nabc_sdx - icnt)/nabc_sdx)**4 + &
              & (float(nabc_bot - kcnt)/nabc_bot)**4))
      enddo
    enddo
  enddo

! ....bot-east beam
  kcnt = nabc_bot + 1
  do k = kmax - nabc_bot + 1, kmax
    kcnt = kcnt - 1
    do j = jmin, jmax
      icnt = nabc_sdx + 1
      do i = imax - nabc_sdx + 1, imax
        icnt = icnt - 1
        deta(i,j,k) = sqrt((deta_maxbot**2)*( &
              & (float(nabc_sdx - icnt)/nabc_sdx)**4 + &
              & (float(nabc_bot - kcnt)/nabc_bot)**4))
      enddo
    enddo
  enddo

! ....bot-south beam
  kcnt = nabc_bot + 1
  do k = kmax - nabc_bot + 1, kmax
    kcnt = kcnt - 1
    jcnt = 0
    do j = jmin, jmin + nabc_sdy - 1
      jcnt = jcnt + 1
      do i = imin, imax
        deta(i,j,k) = sqrt((deta_maxbot**2)*( &
              & (float(nabc_sdy - jcnt)/nabc_sdy)**4 + &
              & (float(nabc_bot - kcnt)/nabc_bot)**4))
      enddo
    enddo
  enddo

! ....bot-north beam
  kcnt = nabc_bot + 1
  do k = kmax - nabc_bot + 1, kmax
    kcnt = kcnt - 1
    jcnt = nabc_sdy + 1
    do j = jmax - nabc_sdy + 1, jmax
      jcnt = jcnt - 1
      do i = imin, imax
        deta(i,j,k) = sqrt((deta_maxbot**2)*( &
              & (float(nabc_sdy - jcnt)/nabc_sdy)**4 + &
              & (float(nabc_bot - kcnt)/nabc_bot)**4))
      enddo
    enddo
  enddo


! Maxwell 3D Damping Function for Cubes:

! ..Top Cubes
if (.not. fsflag) then  

! ....top-southwest cube
  kcnt = 0
  do k = kmin, kmin + nabc_top - 1
    kcnt = kcnt + 1
    jcnt = 0
    do j = jmin, jmin + nabc_sdy - 1
      jcnt = jcnt + 1
      icnt = 0
      do i = imin, imin + nabc_sdx - 1
        icnt = icnt + 1
        deta(i,j,k) = ((deta_maxtop**2 )*( &
              & (float(nabc_sdx - icnt)/nabc_sdx)**6 + &
              & (float(nabc_sdy - jcnt)/nabc_sdy)**6 + &
              & (float(nabc_top - kcnt)/nabc_top)**6))**(1./3.)
      enddo
    enddo
  enddo

! ....top-southeast cube
  kcnt = 0
  do k = kmin, kmin + nabc_top - 1
    kcnt = kcnt + 1
    jcnt = 0
    do j = jmin, jmin + nabc_sdy - 1
      jcnt = jcnt + 1
      icnt = nabc_sdx + 1
      do i = imax - nabc_sdx + 1, imax
        icnt = icnt - 1
        deta(i,j,k) = ((deta_maxtop**2 )*( &
              & (float(nabc_sdx - icnt)/nabc_sdx)**6 + &
              & (float(nabc_sdy - jcnt)/nabc_sdy)**6 + &
              & (float(nabc_top - kcnt)/nabc_top)**6))**(1./3.)
      enddo
    enddo
  enddo

! ....top-northwest cube
  kcnt = 0
  do k = kmin, kmin + nabc_top - 1
    kcnt = kcnt + 1
    jcnt = nabc_sdy + 1
    do j = jmax - nabc_sdy + 1, jmax
      jcnt = jcnt - 1
      icnt = 0
      do i = imin, imin + nabc_sdx - 1
        icnt = icnt + 1
        deta(i,j,k) = ((deta_maxtop**2)*( &
              & (float(nabc_sdx - icnt)/nabc_sdx)**6 + &
              & (float(nabc_sdy - jcnt)/nabc_sdy)**6 + &
              & (float(nabc_top - kcnt)/nabc_top)**6))**(1./3.)
      enddo
    enddo
  enddo

! ....top-northeast cube
  kcnt = 0
  do k = kmin, kmin + nabc_top - 1
    kcnt = kcnt + 1
    jcnt = nabc_sdy + 1
    do j = jmax - nabc_sdy + 1, jmax
      jcnt = jcnt - 1
      icnt = nabc_sdx + 1
      do i = imax - nabc_sdx + 1, imax
        icnt = icnt - 1
        deta(i,j,k) = ((deta_maxtop**2)*( & 
          & (float(nabc_sdx - icnt)/nabc_sdx)**6 + &
          & (float(nabc_sdy - jcnt)/nabc_sdy)**6 + &
          & (float(nabc_top - kcnt)/nabc_top)**6))**(1./3.)
      enddo
    enddo
  enddo

endif

! ..Bottom Cubes

! ....bot-southwest cube
  kcnt = nabc_bot + 1
  do k = kmax - nabc_bot + 1, kmax
    kcnt = kcnt - 1
    jcnt = 0
    do j = jmin, jmin + nabc_sdy - 1
      jcnt = jcnt + 1
      icnt = 0
      do i = imin, imin + nabc_sdx - 1
        icnt = icnt + 1
        deta(i,j,k) = ((deta_maxbot**2)*( &
              & (float(nabc_sdx - icnt)/nabc_sdx)**6 + &
              & (float(nabc_sdy - jcnt)/nabc_sdy)**6 + &
              & (float(nabc_bot - kcnt)/nabc_bot)**6))**(1./3.)
      enddo
    enddo
  enddo

! ....bot-southeast cube
  kcnt = nabc_bot + 1
  do k = kmax - nabc_bot + 1, kmax
    kcnt = kcnt - 1
    jcnt = 0
    do j = jmin, jmin + nabc_sdy - 1
      jcnt = jcnt + 1
      icnt = nabc_sdx + 1
      do i = imax - nabc_sdx + 1, imax
        icnt = icnt - 1
        deta(i,j,k) = ((deta_maxbot**2)*( &
              & (float(nabc_sdx - icnt)/nabc_sdx)**6 + &
              & (float(nabc_sdy - jcnt)/nabc_sdy)**6 + &
              & (float(nabc_bot - kcnt)/nabc_bot)**6))**(1./3.)
      enddo
    enddo
  enddo

! ....bot-northwest cube
  kcnt = nabc_bot + 1
  do k = kmax - nabc_bot + 1, kmax
    kcnt = kcnt - 1
    jcnt = nabc_sdy + 1
    do j = jmax - nabc_sdy + 1, jmax
      jcnt = jcnt - 1
      icnt = 0
      do i = imin, imin + nabc_sdx - 1
        icnt = icnt + 1
        deta(i,j,k) = ((deta_maxbot**2)*( &
              & (float(nabc_sdx - icnt)/nabc_sdx)**6 + &
              & (float(nabc_sdy - jcnt)/nabc_sdy)**6 + &
              & (float(nabc_bot - kcnt)/nabc_bot)**6))**(1./3.)
      enddo
    enddo
  enddo

! ....bot-northeast cube
  kcnt = nabc_bot + 1
  do k = kmax - nabc_bot + 1, kmax
    kcnt = kcnt - 1
    jcnt = nabc_sdy + 1
    do j = jmax - nabc_sdy + 1, jmax
      jcnt = jcnt - 1
      icnt = nabc_sdx + 1
      do i = imax - nabc_sdx + 1, imax
        icnt = icnt - 1
        deta(i,j,k) = ((deta_maxbot**2)*( & 
              & (float(nabc_sdx - icnt)/nabc_sdx)**6 + &
              & (float(nabc_sdy - jcnt)/nabc_sdy)**6 + &
              & (float(nabc_bot - kcnt)/nabc_bot)**6))**(1./3.)
      enddo
    enddo
  enddo


return
end subroutine
