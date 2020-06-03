! NameListRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read EPEC namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListRead (IMPURITY, NEUTRAL, FREQ, PHASE, STAGE2, INTP, RESTART, OLD, TSTART, TEND, DT) &
     bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  integer (kind = c_int),    intent (inout) :: IMPURITY
  integer (kind = c_int),    intent (inout) :: NEUTRAL
  integer (kind = c_int),    intent (inout) :: FREQ
  integer (kind = c_int),    intent (inout) :: PHASE
  integer (kind = c_int),    intent (inout) :: STAGE2
  integer (kind = c_int),    intent (inout) :: INTP
  integer (kind = c_int),    intent (inout) :: RESTART
  integer (kind = c_int),    intent (inout) :: OLD
  real    (kind = c_double), intent (inout) :: TSTART
  real    (kind = c_double), intent (inout) :: TEND
  real    (kind = c_double), intent (inout) :: DT
  
  namelist /ScanInputs/ IMPURITY, NEUTRAL, FREQ, STAGE2, PHASE, INTP, RESTART, OLD, TSTART, TEND, DT
  
  open  (unit = 100, file = 'namelist.txt', status = 'old')
  read  (unit = 100, nml  = ScanInputs)
  close (unit = 100)

endsubroutine namelistRead
