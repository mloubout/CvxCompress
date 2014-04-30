subroutine writeGocadSnapShot(ijk_dec, nx_dec, ny_dec, nz_dec, &
                            & nx, ny, nz, ilo, jlo, klo, iup, jup, kup, &
                            & prop, filename) 

! Writes snapshot files of the wavefields in gocad or tecplot360 formats.
!
! o Last Modified:
!   06-06-12  Initial edits.
!
! o Written by: Kurt Nihei; 06-06-12
!

implicit none

integer :: ijk_dec  ! snapshot decimation factor
integer :: nx_dec, ny_dec, nz_dec  ! snapshot dimensions (after decimation)
integer :: nx, ny, nz           ! #cells in full earth model
integer :: ilo, jlo, klo
integer :: iup, jup, kup
character (len=21) :: filename

integer :: i, j, k
integer :: status
integer*8 :: ijk, nsamp
real, dimension (ilo:iup,jlo:jup,klo:kup) :: prop 
real, allocatable, dimension(:) :: tempswap  ! for byte swapping

! Allocate temporary array for conversion from little endian to 
!   big endian:
nsamp = nx_dec*ny_dec*nz_dec
allocate(tempswap(nsamp))

ijk = 0
do i = 1, nx, ijk_dec
  do j = 1, ny, ijk_dec
    do k = 1, nz, ijk_dec
      ijk = ijk + 1
      tempswap(ijk) = prop(i,j,k)  ! store in a 1D array
    enddo
  enddo
enddo
call swapBytes(tempswap, nsamp)  ! convert to big endian

open(unit=2, file=filename, status='unknown', action='write', &
   & iostat=status, form='binary', recl=4*nsamp)
write(2) (tempswap(ijk), ijk = 1, nsamp)
close(unit=2)

return
end subroutine
