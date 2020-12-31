! Function to read PHASE namelist

subroutine NameListRead (NFLOW, STAGE5, INTF, INTN, INTU, OLD, FREQ, LIN, MID, DT, TSTART, TEND, SCALE, PMAX,&
     NCTRL, xTCTRL, xICTRL, xPCTRL) bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  integer (kind = c_int),    intent (inout) :: STAGE5
  integer (kind = c_int),    intent (inout) :: NFLOW
  integer (kind = c_int),    intent (inout) :: INTF
  integer (kind = c_int),    intent (inout) :: INTN
  integer (kind = c_int),    intent (inout) :: INTU
  integer (kind = c_int),    intent (inout) :: OLD
  integer (kind = c_int),    intent (inout) :: FREQ
  integer (kind = c_int),    intent (inout) :: LIN
  integer (kind = c_int),    intent (inout) :: MID
  real    (kind = c_double), intent (inout) :: DT
  real    (kind = c_double), intent (inout) :: TSTART
  real    (kind = c_double), intent (inout) :: TEND
  real    (kind = c_double), intent (inout) :: SCALE
  real    (kind = c_double), intent (inout) :: PMAX
  integer (kind = c_int),    intent (inout) :: NCTRL
  
  real (kind = c_double), dimension (*), intent (inout) :: xTCTRL 
  real (kind = c_double), dimension (*), intent (inout) :: xICTRL 
  real (kind = c_double), dimension (*), intent (inout) :: xPCTRL

  real, dimension (:), allocatable :: TCTRL
  real, dimension (:), allocatable :: ICTRL
  real, dimension (:), allocatable :: PCTRL

  integer          :: i
  double precision :: pi
 
  namelist /PHASE_CONTROL/   STAGE5, NFLOW, INTF, INTN, INTU, OLD, FREQ, LIN, MID, DT, TSTART, TEND, SCALE, PMAX
  namelist /PHASE_WAVEFORM/  NCTRL
  namelist /PHASE_WAVEFORM1/ TCTRL, ICTRL, PCTRL
  
  open  (unit = 100, file = 'Inputs/Phase.nml', status = 'old')
  read  (unit = 100, nml  = PHASE_CONTROL) 
  close (unit = 100)
  
  open  (unit = 100, file = 'Inputs/Waveform.nml', status = 'old')
  read  (unit = 100, nml  = PHASE_WAVEFORM)

  allocate (TCTRL (NCTRL))
  allocate (ICTRL (NCTRL))
  allocate (PCTRL (NCTRL))

  read  (unit = 100, nml  = PHASE_WAVEFORM1)

  close (unit = 100)
  
  pi = 4.*atan(1.)
  do i = 1, NCTRL
     xTCTRL (i) = TCTRL (i)
     xICTRL (i) = ICTRL (i)
     xPCTRL (i) = PCTRL (i) * pi
  enddo
 
  deallocate (TCTRL)
  deallocate (ICTRL)
  deallocate (PCTRL)

endsubroutine namelistRead
