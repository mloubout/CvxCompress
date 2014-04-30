subroutine readModelProperties(nx, ny, nz, ilo, jlo, klo, &
                    & iup, jup, kup, orderaxes, propFile, bigendianflag, &
                    & freesurfflag, prop, istat)

! Prepares a model property for FD modeling of a single shot: 
! (1) reads flat binary model property files; 
! (2) loads part of model to be used in FD modeling into property array; 
! (3) pads top and bottom of model, as required for absorbing boundaries.
!
! o Last Modified:
!   06-06-12  Initial edits.
!
! o Written by: Kurt Nihei; 06-06-12
!

implicit none

integer :: nx, ny, nz           ! #cells in full earth model
integer :: ilo, jlo, klo  ! array dimensions for FD model
integer :: iup, jup, kup  ! array dimensions for FD model
character (len=3) :: orderaxes  ! fast-med-slow axes ordering: xyz, zxy, xyz
character (len=128) :: propFile
logical :: bigendianflag  ! swap4bytes 
logical :: freesurfflag  ! free-surface flag for top of model

integer :: istat
integer :: i, j, k
integer :: ijk, nsamp
real, dimension (ilo:iup,jlo:jup,klo:kup) :: prop
real, allocatable, dimension(:) :: tempswap  ! for byte swapping


! Open flat binary model property file:
istat = 0  ! initial file status 
open(unit=2, file=propFile, status='old', action='read', iostat=istat, &
   & form='binary')

if (istat /= 0) then  ! check file path:
  write(*,*) 'Could not open model property binary file: ', propFile
  write(*,*) '..check filepath/filename!'
  write(*,*) '..setting this property to zero!'
  prop = 0

else

  ! Read file:
  write(*,*) 'Reading model property file: ', propFile
  if (orderaxes == 'xyz') then  ! fast-xyz
    write(*,*) 'Ordering of fast|med|slow axes is fastxyz.'
    read(2) (((prop(i,j,k), i=1,nx), j=1,ny), k=1,nz)

  elseif (orderaxes == 'zxy') then  ! fast-zxy
    write(*,*) 'Ordering of fast|med|slow axes is fastzxy.'
    read(2) (((prop(i,j,k), k=1,nz), i=1,nx), j=1,ny)

  elseif (orderaxes == 'zyx') then  ! fast-zyx
    write(*,*) 'Ordering of fast|med|slow axes is fastzyx.'
    read(2) (((prop(i,j,k), k=1,nz), j=1,ny), i=1,nx)

  else
    write(*,*) 'Non-standard ordering of fast-med-slow axis in model file.'
    write(*,*) 'Standard orderings are (fast|med|slow): xyz, zxy, zyx'
    write(*,*) 'Encountered order is: ', orderaxes
    write(*,*) 'Terminating model file read!'
    stop
  endif

  close(unit=2)

  ! Convert from big endian (e.g., gocad) to little endian (linux):
  if (bigendianflag) then
    nsamp = nx*ny*nz
    allocate(tempswap(nsamp))
    ijk = 0
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ijk = ijk + 1
          tempswap(ijk) = prop(i,j,k)  ! store in a 1D array
        enddo
      enddo
    enddo

    call swapBytes(tempswap, nsamp)  ! convert to little endian

    ijk = 0
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ijk = ijk + 1
          prop(i,j,k) = tempswap(ijk)  ! put back into 3D array
        enddo
      enddo
    enddo

    deallocate(tempswap)
  endif

  ! Extract part of full model for use in FD modeling
  ! (note: FD model retains global (i,j,k) indexing):
  do k = klo, kup   
    do j = jlo, jup 
      do i = ilo, iup 
        prop(i,j,k) = prop(i,j,k)  
      enddo
    enddo
  enddo


  ! Add bottom & top absorbing boundary condition slabs
  ! via vertical extrapolation of model properties:
  do i = ilo, iup 
    do j = jlo, jup 
      do k = nz+1, kup   
        prop(i,j,k) = prop(i,j,nz)  ! bottom abc slab
      enddo
    enddo
  enddo

  if (.not. freesurfflag) then  
    do i= ilo, iup 
      do j= jlo, jup 
        do k= klo, 0
          prop(i,j,k) = prop(i,j,1)  ! top abc slab
        enddo
      enddo
    enddo
  endif

endif

return
end subroutine
