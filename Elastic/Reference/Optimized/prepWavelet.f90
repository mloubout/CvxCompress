subroutine prepwavelet(tmax, fmax, dti_raw, dti, dti_segy, &
                     & nwav_raw0, nwav, swav_raw0, swav, swav_int)

! Uses a DFT to resample an input source wavlet to specified dti and 
! integrates this resampled wavelet. Returns the resampled, time-integrated
! wavelet.
!
! tmax      = maximum time (s) in FD simulation. 
! fmax      = maximum frequency (Hz) to be preserved after resampling & 
!            integration (i.e., Butterworth filter corner freq.; -3dB 
!            down)
! dti_raw  = time sampling (s) of original (input) wavelet
! dti       = time sampling (s) of desired (output) wavelet
! dti_segy  = time sampling (s) of SEG-Y traces (used for anti-alias filter)
! nwav_raw0 = #samples in original wavelet
! nwav_raw  = #samples in original wavelet after zero-padding on trailing
!             edge
! nwav      = #samples in resampled (and integrated) wavelet
! swav_raw0 = input wavelet before resampling & integration
! swav_raw  = input wavelet before resampling & integration after 
!             zero-padding on trailing edge
! swav      = output wavelet after resampling (& filtering)
! swav_int  = output wavelet after resampling & integration (& filtering)
!
! o Last Modified:
!   11-21-13  Version 6.3: Removed 75% drop in fmax based on comparison of
!             average power above and below fmax; issuing a warning
!             instead.
!   06-14-13  Version 6.2: Changed sign of wavelet dB spectrum in
!             prepwavelet.f90.
!   03-07-13 Version 6.0: Modified prepWavelet.f90 to write frequency
!            spectra of the input and resampled/filtered wavelet. Also,
!            added check on Nyquist frequencies (of input wavelet, of
!            FD resampled wavelet, of SEG-Y resampled wavelet). The
!            smallest of these Nyquist frequencies is compared to 
!            the -40dB frequency computed from fmax. The smaller of the
!            two is used as the target for the 5-pole Butterworth filter.
!   11-02-12 Added OMP to IDFT.
!   11-02-12 Fixed bug in QC swav write.
!   11-02-12 Fixed over-specification bug in swav_raw assignment.
!   10-11-12 Replaced Gaussian filter with Butterworth filter.
!   10-10-12 Cleaned-up DFT and low-pass filter code; added OMP to DFT.
!   10-10-12 Removing any DC shift of time-domain resampled wavelet by
!            simple subtraction of the amplitude of the first sample.
!   10-09-12 For integrated wavelet, estimate value of complex wavelet at 
!            zero frequency by linear extrapolation from the RHS. 
!   10-09-12 Added zero-padding on trailing edge of input wavelet to reduce
!            DC artifacts from DFT. 
!   12-19-11 Added Tecplot360 headers to qc_wavelet files, and renamed 
!            these qc_wavelet.dat.
!   12-13-11 Removed multiplication of integrated wavelet by dti; this is
!            not necessary.
!   11-15-11 Adapted from prototype code prepwavelet.m
!
! o Created:  
!   11-16-11; KTN; /users/knih/NumericalCodes/Ssg2d/Ssg2d_O3/Ver2
!

implicit none

integer :: nwav_raw0, nwav
real :: tmax, fmax
real :: dti_raw, dti, dti_segy
real, dimension(nwav_raw0) :: swav_raw0
real, dimension(nwav) :: swav, swav_int

integer, parameter :: orderbutt = 5  ! order of Butterworth low-pass filter
integer :: ii, kk 
integer :: nwav_raw
integer :: nf_raw, nf
integer :: nfilt_start, nfilt_end, nfilt   
integer :: kk_fsafe
real, parameter :: zpad = 4.  ! zero-padding factor to reduce DC artifacts
real :: tmax_raw, df_raw, df
real :: ff, ww, tt
real :: dc_swav, dc_swav_int
real :: pi
real :: wmax, filtbutt
real :: fmax_40db, fnyq_wavelet, fnyq_fd, fnyq_segy, fsafe, ave_amp_ratio
real, allocatable, dimension(:) :: swav_raw
complex :: fswav_ave_below, fswav_ave_above
complex, allocatable, dimension(:) :: fswav, fswav_int, fswav_sav
complex, allocatable, dimension(:) :: fswav_resfil, fswav_int_resfil
complex, dimension(nwav) :: cswav, cswav_int


pi = acos(-1.)

! Perform check on important frequencies:
fmax_40db = fmax*2.5  ! converts -3dB fmax to -40dB value
fnyq_wavelet = 1./(2*dti_raw)  ! Nyquist freq of wavelet (-60dB)
fnyq_fd = 1./(2*dti)    ! Nyquist freq of FD resamp wavelet (-60dB)
fnyq_segy = 1./(2*dti_segy) ! Nyquist freq of SEG-Y resamp wavelet (-60dB)
fsafe = min(fmax_40db, fnyq_wavelet, fnyq_fd, fnyq_segy)
write(*,*) '....computing important frequencies'
write(*,*) '.....fmax_40dB = ', fmax_40db
write(*,*) '.....fNyquist_wavelet = ', fnyq_wavelet
write(*,*) '.....fNyquist_FD = ', fnyq_fd
write(*,*) '.....fNyquist_SEG-Y = ', fnyq_segy
write(*,*) '.....fsafe = ', fsafe


! Zero pad input wavelet to tmax to reduce DFT DC artifacts:
write(*,*) '....zero-padding wavelet'
nwav_raw = int(zpad*tmax/dti_raw)  ! zero-pad to reduce DC artifacts 

allocate (swav_raw(nwav_raw))
allocate(fswav(nwav_raw), fswav_int(nwav_raw), fswav_sav(nwav_raw))
allocate(fswav_resfil(nwav), fswav_int_resfil(nwav))

if (nwav_raw0 < nwav_raw) then
  do ii = 1, nwav_raw0
    swav_raw(ii) = swav_raw0(ii)
  enddo
  do ii = nwav_raw0+1, nwav_raw
    swav_raw(ii) = 0.  ! zero pad trailing end of input wavelet 
  enddo
else
  do ii = 1, nwav_raw
    swav_raw(ii) = swav_raw0(ii)
  enddo
endif


! Define properties of input wavelet (after zero-padding):
tmax_raw = dti_raw*float(nwav_raw - 1)
nf_raw = nwav_raw/2  ! freq-domain needs only 1/2 the points 
df_raw = 1/tmax_raw


! Compute DFT of input wavelet:
write(*,*) '....computing DFT of input wavelet'
!$OMP PARALLEL DO PRIVATE(ff, ww, ii, tt)
do kk = 1, nf_raw
  ff = df_raw*(kk-1)
  ww = 2*pi*ff

  fswav(kk) = cmplx(0.,0.)
  do ii = 1, nwav_raw
    tt = dti_raw*float(ii - 1)
    fswav(kk) = fswav(kk) + swav_raw(ii)*cmplx(cos(ww*tt), sin(ww*tt))
  enddo
enddo
!$OMP END PARALLEL DO
fswav_sav = fswav  ! save input wavelet freq spectrum before filtering


! Check average power above and below fsafe; reset fmax if necessary:
if (fsafe /= fmax_40db) then  ! don't reset fmax if it is the limiting freq
  fswav_ave_below = 0.
  kk_fsafe = fsafe/df_raw + 1
  do kk = 1, kk_fsafe
    fswav_ave_below = fswav_ave_below + fswav(kk)
  enddo
  fswav_ave_below = fswav_ave_below/float(kk_fsafe)

  fswav_ave_above = 0.
  do kk = kk_fsafe+1, nf_raw
    fswav_ave_above = fswav_ave_above + fswav(kk)
  enddo
  fswav_ave_above = fswav_ave_above/float(nf_raw - kk_fsafe + 1)

  write(*,*) '....checking fsafe'
  ave_amp_ratio = cabs(fswav_ave_above)/cabs(fswav_ave_below)
  write(*,*) '.....[ave_spec_above/ave_spec_below]dB = ', &
            & 20.*log10(ave_amp_ratio)
  if (ave_amp_ratio > 1.e-3) then  ! 1.e-3 is equiv to -60dB
    write(*,*) '.....Warning! Power above fsafe exceeds -60dB!'
    write(*,*) '.......Numerical dispersion may result!'
  else
    write(*,*) '.....power above fsafe is below -60dB limit'
  endif
endif


! Apply low-pass filter to transformed wavelet:
write(*,*) '....applying low-pass filter'
wmax = 2.*pi*fmax  ! Butterworth filter corner angular freq. 
do kk = 1, nf_raw
  ff = df_raw*(kk-1)
  ww = 2*pi*ff
  filtbutt = 1./sqrt(1. + (ww/wmax)**(2.*orderbutt))
  fswav(kk) = fswav(kk)*filtbutt 
enddo 


! Perform time integration by multiplication by 1/(iw):
write(*,*) '....integrating wavelet'
do kk = 2, nf_raw  ! kk=1 is singular; taken care of after do-loop
  ff = df_raw*(kk-1)
  ww = 2*pi*ff
  fswav_int(kk) = fswav(kk)/cmplx(0.,-ww)
enddo

! ..estimate first value of fswav_int via linear extrapolation:
fswav_int(1) = 2.*fswav_int(2) - fswav_int(3)


! Compute IDF at new (resampled) time locations:
write(*,*) '....computing inverse DFT of wavelet'
!$OMP PARALLEL DO PRIVATE(tt, kk, ff, ww)
do ii = 1, nwav
  tt = dti*float(ii - 1)
  cswav(ii) = cmplx(0.,0.)
  cswav_int(ii) = cmplx(0.,0.)
  do kk = 1, nf_raw
    ff = df_raw*(kk-1)
    ww = 2*pi*ff
    cswav(ii) = cswav(ii) + (1./nf_raw)*fswav(kk)* &
              & cmplx(cos(ww*tt), -sin(ww*tt))
    cswav_int(ii) = cswav_int(ii) + (1./nf_raw)* &
                  & fswav_int(kk)*cmplx(cos(ww*tt), -sin(ww*tt))
  enddo
enddo
!$OMP END PARALLEL DO


! Extract real parts:
swav = real(cswav)
swav_int = real(cswav_int)


! Remove DC shift of wavelet by subtracting amplitude of first sample:
write(*,*) '....removing DC shift'
dc_swav = swav(1)
dc_swav_int = swav_int(1)
swav = swav - dc_swav
swav_int = swav_int - dc_swav_int


! QC Wavelet Spectra: compute DFT of resampled-filtered wavelet
write(*,*) '....computing DFT of resampled-filtered wavelet'
! ..define properties of input wavelet (after zero-padding):
nf = nwav/2  ! freq-domain needs only 1/2 the points 
df = 1/(dti*float(nwav - 1))

!$OMP PARALLEL DO PRIVATE(ff, ww, ii, tt)
do kk = 1, nf
  ff = df*(kk-1)
  ww = 2*pi*ff

  fswav_resfil(kk) = cmplx(0.,0.)
  fswav_int_resfil(kk) = cmplx(0.,0.)
  do ii = 1, nwav
    tt = dti*float(ii - 1)
    fswav_resfil(kk) = fswav_resfil(kk) + &
                     & swav(ii)*cmplx(cos(ww*tt), sin(ww*tt))
    fswav_int_resfil(kk) = fswav_int_resfil(kk) + &
                     & swav_int(ii)*cmplx(cos(ww*tt), sin(ww*tt))
  enddo
enddo
!$OMP END PARALLEL DO


! Write resampled wavelet time series to file for QC:
write(*,*) '....writing wavelet QC files'
open(unit=1, file='qc_wavelet_time.dat', status='unknown')
write(1,*) 'TITLE="qc_wavelet_time"'
write(1,*) 'VARIABLES="time (s)","Amplitude"'

write(1,*) 'ZONE T="input wavelet (raw)", F=POINT'
do ii = 1, nwav_raw0
  tt = dti_raw*float(ii - 1)
  write(1,*) tt, swav_raw0(ii)
enddo

write(1,*) 'ZONE T="resampled (filtered) wavelet", F=POINT'
do ii = 1, nwav
  tt = dti*float(ii - 1)
  write(1,*) tt, swav(ii)
enddo

write(1,*) 'ZONE T="integrated-resampled (filtered) wavelet", F=POINT'
do ii = 1, nwav
  tt = dti*float(ii - 1)
  write(1,*) tt, swav_int(ii)
enddo
close(unit=1)


! Write resampled wavelet power spectra to file for QC:
open(unit=1, file='qc_wavelet_freq.dat', status='unknown')
write(1,*) 'TITLE="qc_wavelet_freq"'
write(1,*) 'VARIABLES="frequency (Hz)","Power Spectrum (dB)"'

write(1,*) 'ZONE T="input wavelet (raw)", F=POINT'
do kk = 1, nf_raw
  ff = df_raw*(kk-1)
  write(1,*) ff, 20.*log10(cabs(fswav_sav(kk)))
enddo

write(1,*) 'ZONE T="resampled (filtered) wavelet", F=POINT'
do kk = 1, nf
  ff = df*(kk-1)
  write(1,*) ff, 20.*log10(cabs(fswav_resfil(kk)))
enddo

write(1,*) 'ZONE T="integrated-resampled (filtered) wavelet", F=POINT'
do kk = 1, nf
  ff = df*(kk-1)
  write(1,*) ff, 20.*log10(cabs(fswav_int_resfil(kk)))
enddo
close(unit=1)


return
end subroutine
