! NameListRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read IslandDynamics namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListRead (FLUX_NTOR, FLUX_MMIN, FLUX_MMAX, NEO_INTF, NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_YN,&
     PHASE_INTN, PHASE_STAGE2, PHASE_OLD, RESTART, TSTART, TEND, DT) &
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
  real    (kind = c_double), intent (inout) :: NEO_YN
  integer (kind = c_int),    intent (inout) :: PHASE_STAGE2
  integer (kind = c_int),    intent (inout) :: PHASE_INTN
  integer (kind = c_int),    intent (inout) :: PHASE_OLD
  integer (kind = c_int),    intent (inout) :: RESTART
  real    (kind = c_double), intent (inout) :: TSTART
  real    (kind = c_double), intent (inout) :: TEND
  real    (kind = c_double), intent (inout) :: DT
  
  namelist /IslandDynamicsInputs/ FLUX_NTOR, FLUX_MMIN, FLUX_MMAX, NEO_INTF, NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_YN,&
  PHASE_INTN, PHASE_STAGE2, PHASE_OLD, RESTART, TSTART, TEND, DT
  
  open  (unit = 100, file = 'namelist.txt', status = 'old')
  read  (unit = 100, nml  = IslandDynamicsInputs)
  close (unit = 100)

endsubroutine namelistRead
