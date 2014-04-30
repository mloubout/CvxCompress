subroutine interpolateTrace(it_left, t_out, dt_in, dt_out, data_in, &
                          & data_out)
 
!  Interpolates traces using 4 point sinc interpolation with a Kaiser 
!  windowing function.
!
!
!  o Last Modified:  05-08-12  
!
!  o Written by Kurt T. Nihei; 05-07-12


implicit none

real :: bessi0  ! function to compute I0 (modified Bessel funct of 0 order)

integer :: it_left
real :: t_out, dt_in, dt_out
real, dimension(4) :: data_in
real :: data_out

integer :: is
real :: dt_frac
real :: t_samp, t_win, win
real, dimension(4) :: fsinc

real, parameter :: b = 4.14  ! optimal Kaiser window param for kmax = 2pi/3
real, parameter :: r = 2.    ! half-width of sinc interpolator
real :: pi
pi = acos(-1.)


! time from left sample to interpol pt [normalized]
dt_frac = (t_out - dt_in*float(it_left))/dt_in 

do is = 1, 4
  ! times at which to sample sinc func [normalized]
  t_samp = float(is) - 2 - dt_frac  

  ! compute Kaiser window:
  if (t_samp <= r) then
    t_win = b*sqrt(1 - (t_samp/r)**2)
  else
    t_win = 0.
  endif
  win = bessi0(t_win)/bessi0(b)

  ! compute sinc interpolation function:
  if (t_samp == 0.) then
    fsinc(is) = 1.  ! no interpolation required
  else
    fsinc(is) = win*sin(t_samp*pi)/(t_samp*pi) 
!!!!
!!!!write(*,*)
!!!!write(*,*) 't_out = ', t_out
!!!!write(*,*) 'floor(t_out/dti) = ', floor(t_out/dt_in)
!!!!write(*,*) 'it_left = ', it_left
!!!!write(*,*) 't_samp = ', t_samp
!!!!write(*,*) 't_win = ', t_win
!!!!write(*,*) 'win = ', win
!!!!write(*,*) 'dt_frac = ', dt_frac
!!!!write(*,*) 'is, fsinc(is) = ', is, fsinc(is)
!!!!
  endif

enddo

! Interpolate trace with 4 point sinc:
data_out = fsinc(1)*data_in(1) + fsinc(2)*data_in(2) + &
         & fsinc(3)*data_in(3) + fsinc(4)*data_in(4)

    
return
end subroutine
