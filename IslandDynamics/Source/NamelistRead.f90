! NameListRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read IslandDynamics namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListRead (FLUX_NTOR, FLUX_MMIN, FLUX_MMAX,&
     NEO_INTP, NEO_INTF, NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN,&
     PHASE, PHASE_MID, PHASE_INTN, PHASE_INTU, PHASE_STAGE5, PHASE_OLD, PHASE_FREQ, PHASE_LIN, PHASE_SCALE,&
     PHASE_CHIR, RESTART, TSTART, TEND, DT) &
     bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  integer (kind = c_int),    intent (inout) :: FLUX_NTOR
  integer (kind = c_int),    intent (inout) :: FLUX_MMIN
  integer (kind = c_int),    intent (inout) :: FLUX_MMAX
  integer (kind = c_int),    intent (inout) :: NEO_INTP
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
  integer (kind = c_int),    intent (inout) :: PHASE_INTU
  integer (kind = c_int),    intent (inout) :: PHASE_OLD
  integer (kind = c_int),    intent (inout) :: PHASE_FREQ
  integer (kind = c_int),    intent (inout) :: PHASE_LIN
  real    (kind = c_double), intent (inout) :: PHASE_SCALE
  real    (kind = c_double), intent (inout) :: PHASE_CHIR
  integer (kind = c_int),    intent (inout) :: RESTART
  real    (kind = c_double), intent (inout) :: TSTART
  real    (kind = c_double), intent (inout) :: TEND
  real    (kind = c_double), intent (inout) :: DT
  
  namelist /ISLANDDYNAMICS_CONTROL/ FLUX_NTOR, FLUX_MMIN, FLUX_MMAX,&
       NEO_INTP, NEO_INTF, NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN,&
       PHASE, PHASE_MID, PHASE_INTN, PHASE_INTU, PHASE_STAGE5, PHASE_OLD, PHASE_FREQ, PHASE_LIN, PHASE_SCALE,&
       PHASE_CHIR, RESTART, TSTART, TEND, DT
  
  open  (unit = 100, file = 'Inputs/Island.nml', status = 'old')
  read  (unit = 100, nml  = ISLANDDYNAMICS_CONTROL)
  close (unit = 100)

endsubroutine namelistRead
