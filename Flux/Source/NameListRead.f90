! NameListRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read FLUX namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListRead (INTG, NPSI, NTHETA, NNC, NTOR, QFLG, Q95, H0, ACC, ETA, DR, MMIN, MMAX, PSILIM, TIME) &
     bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  integer (kind = c_int),    intent (inout) :: NPSI
  integer (kind = c_int),    intent (inout) :: NTHETA
  integer (kind = c_int),    intent (inout) :: NNC
  integer (kind = c_int),    intent (inout) :: NTOR
  integer (kind = c_int),    intent (inout) :: QFLG
  real    (kind = c_double), intent (inout) :: Q95
  real    (kind = c_double), intent (inout) :: H0
  real    (kind = c_double), intent (inout) :: ACC
  real    (kind = c_double), intent (inout) :: ETA
  real    (kind = c_double), intent (inout) :: DR
  integer (kind = c_int),    intent (inout) :: MMIN
  integer (kind = c_int),    intent (inout) :: MMAX
  real    (kind = c_double), intent (inout) :: PSILIM
  real    (kind = c_double), intent (inout) :: TIME
  integer (kind = c_int),    intent (inout) :: INTG
  
  namelist /FLUX_CONTROL/ NPSI, NTHETA, NNC, NTOR, QFLG, Q95, H0, ACC, ETA, DR, MMIN, MMAX, PSILIM, TIME, INTG
  
  open  (unit = 100, file = 'Inputs/Flux.in', status = 'old')
  read  (unit = 100, nml  = FLUX_CONTROL)
  close (unit = 100)

endsubroutine namelistRead
