! NameListRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read IslandDynamics namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListRead (FLUX_NTOR, FLUX_MMIN, FLUX_MMAX,&
     NEO_INTF, NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN,&
     PHASE, PHASE_MID, PHASE_INTN, PHASE_STAGE5, PHASE_OLD, PHASE_VER2, PHASE_FREQ, PHASE_LIN, PHASE_SCALE,&
     RESTART, TSTART, TEND, DT) &
     bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  integer (kind = c_int),    intent (inout) :: FLUX_NTOR
  integer (kind = c_int),    intent (inout) :: FLUX_MMIN
  integer (kind = c_int),    intent (inout) :: FLUX_MMAX
  integer (kind = c_int),    intent (inout) :: NEO_INTF
  integer (kind = c_int),    intent (inout) :: NEO_IMPURITY
  integer (kind = c_int),    intent (inout) :: NEO_NEUTRAL
  integer (kind = c_int),    intent (inout) :: NEO_FREQ
  integer (kind = c_int),    intent (inout) :: NEO_NTYPE
  real    (kind = c_double), intent (inout) :: NEO_NN
  real    (kind = c_double), intent (inout) :: NEO_LN
  real    (kind = c_double), intent (inout) :: NEO_YN
  integer (kind = c_int),    intent (inout) :: PHASE
  integer (kind = c_int),    intent (inout) :: PHASE_MID
  integer (kind = c_int),    intent (inout) :: PHASE_STAGE5
  integer (kind = c_int),    intent (inout) :: PHASE_INTN
  integer (kind = c_int),    intent (inout) :: PHASE_OLD
  integer (kind = c_int),    intent (inout) :: PHASE_VER2
  integer (kind = c_int),    intent (inout) :: PHASE_FREQ
  integer (kind = c_int),    intent (inout) :: PHASE_LIN
  real    (kind = c_double), intent (inout) :: PHASE_SCALE
  integer (kind = c_int),    intent (inout) :: RESTART
  real    (kind = c_double), intent (inout) :: TSTART
  real    (kind = c_double), intent (inout) :: TEND
  real    (kind = c_double), intent (inout) :: DT
  
  namelist /ISLANDDYNAMICS_CONTROL/ FLUX_NTOR, FLUX_MMIN, FLUX_MMAX,&
       NEO_INTF, NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN,&
       PHASE, PHASE_MID, PHASE_INTN, PHASE_STAGE5, PHASE_OLD, PHASE_VER2, PHASE_FREQ, PHASE_LIN, PHASE_SCALE,&
       RESTART, TSTART, TEND, DT
  
  open  (unit = 100, file = 'Inputs/Island.in', status = 'old')
  read  (unit = 100, nml  = ISLANDDYNAMICS_CONTROL)
  close (unit = 100)

endsubroutine namelistRead
