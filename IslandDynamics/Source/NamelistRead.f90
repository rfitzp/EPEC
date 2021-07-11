! NameListRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read IslandDynamics namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListRead (PHASE_MID, PHASE_COPT, PHASE_CORE, PHASE_LIN,&
     PHASE_FREQ, PHASE_FFAC, PHASE_HIGH, PHASE_RATS, PHASE_NATS, PHASE_SCALE,&
     PHASE_CHIR, PHASE_CXD, PHASE_BOOT, PHASE_CURV, PHASE_POLZ,&
     RESTART, TSTART, TEND, DT)&
     bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  integer (kind = c_int),    intent (inout) :: PHASE_MID
  integer (kind = c_int),    intent (inout) :: PHASE_COPT
  integer (kind = c_int),    intent (inout) :: PHASE_LIN
  integer (kind = c_int),    intent (inout) :: PHASE_FREQ
  integer (kind = c_int),    intent (inout) :: PHASE_HIGH
  integer (kind = c_int),    intent (inout) :: PHASE_RATS
  integer (kind = c_int),    intent (inout) :: PHASE_NATS
  integer (kind = c_int),    intent (inout) :: PHASE_CXD
  integer (kind = c_int),    intent (inout) :: PHASE_BOOT
  integer (kind = c_int),    intent (inout) :: PHASE_CURV
  integer (kind = c_int),    intent (inout) :: PHASE_POLZ
  real    (kind = c_double), intent (inout) :: PHASE_FFAC
  real    (kind = c_double), intent (inout) :: PHASE_SCALE
  real    (kind = c_double), intent (inout) :: PHASE_CHIR
  real    (kind = c_double), intent (inout) :: PHASE_CORE
  integer (kind = c_int),    intent (inout) :: RESTART
  real    (kind = c_double), intent (inout) :: TSTART
  real    (kind = c_double), intent (inout) :: TEND
  real    (kind = c_double), intent (inout) :: DT
  
  namelist /ISLANDDYNAMICS_CONTROL/ PHASE_MID, PHASE_COPT, PHASE_CORE, PHASE_LIN,&
       PHASE_FREQ, PHASE_FFAC, PHASE_HIGH, PHASE_RATS, PHASE_NATS, PHASE_SCALE,&
       PHASE_CHIR, PHASE_CXD, PHASE_BOOT, PHASE_CURV, PHASE_POLZ,&
       RESTART, TSTART, TEND, DT
  
  open  (unit = 100, file = 'Inputs/Island.nml', status = 'old')
  read  (unit = 100, nml  = ISLANDDYNAMICS_CONTROL)
  close (unit = 100)

endsubroutine namelistRead
