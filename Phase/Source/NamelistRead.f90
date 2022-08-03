! Function to read PHASE namelist and WAVEFORM namelist

subroutine NameListRead (NFLOW, STAGE5, INTF, INTN, INTU, NATS, OLD, FREQ, LIN, MID, COPT, DT,&
     TSTART, TEND, TOFF, SCALE, PMAX, CHIR, HIGH, RATS, CORE, FFAC, CXD, POEM, BOOT, CURV, POLZ,&
     WALL, TAUW, THRES, MSTOP, TYPE, NCTRL, xTCTRL, xICTRL, xPCTRL,&
     SSTART, SEND, SAMP, SPHA, SPVE, BACK, RPERIOD, RSTART, REND, RPHA)&
     bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  integer (kind = c_int),    intent (inout) :: STAGE5
  integer (kind = c_int),    intent (inout) :: NFLOW
  integer (kind = c_int),    intent (inout) :: INTF
  integer (kind = c_int),    intent (inout) :: INTN
  integer (kind = c_int),    intent (inout) :: INTU
  integer (kind = c_int),    intent (inout) :: NATS
  integer (kind = c_int),    intent (inout) :: OLD
  integer (kind = c_int),    intent (inout) :: FREQ
  integer (kind = c_int),    intent (inout) :: LIN 
  integer (kind = c_int),    intent (inout) :: MID
  integer (kind = c_int),    intent (inout) :: COPT
  integer (kind = c_int),    intent (inout) :: NCTRL
  integer (kind = c_int),    intent (inout) :: HIGH
  integer (kind = c_int),    intent (inout) :: RATS
  integer (kind = c_int),    intent (inout) :: CXD
  integer (kind = c_int),    intent (inout) :: POEM
  integer (kind = c_int),    intent (inout) :: BOOT
  integer (kind = c_int),    intent (inout) :: CURV
  integer (kind = c_int),    intent (inout) :: POLZ
  integer (kind = c_int),    intent (inout) :: TYPE
  integer (kind = c_int),    intent (inout) :: WALL
  integer (kind = c_int),    intent (inout) :: MSTOP
  real    (kind = c_double), intent (inout) :: THRES
  real    (kind = c_double), intent (inout) :: FFAC
  real    (kind = c_double), intent (inout) :: DT
  real    (kind = c_double), intent (inout) :: TSTART
  real    (kind = c_double), intent (inout) :: TEND
  real    (kind = c_double), intent (inout) :: SCALE
  real    (kind = c_double), intent (inout) :: PMAX
  real    (kind = c_double), intent (inout) :: CHIR
  real    (kind = c_double), intent (inout) :: CORE
  real    (kind = c_double), intent (inout) :: TAUW
  real    (kind = c_double), intent (inout) :: TOFF
  real    (kind = c_double), intent (inout) :: SSTART
  real    (kind = c_double), intent (inout) :: SEND
  real    (kind = c_double), intent (inout) :: SAMP
  real    (kind = c_double), intent (inout) :: SPHA
  real    (kind = c_double), intent (inout) :: SPVE
  real    (kind = c_double), intent (inout) :: BACK
  real    (kind = c_double), intent (inout) :: RPERIOD
  real    (kind = c_double), intent (inout) :: RSTART
  real    (kind = c_double), intent (inout) :: REND
  real    (kind = c_double), intent (inout) :: RPHA

  real (kind = c_double), dimension (*), intent (inout) :: xTCTRL 
  real (kind = c_double), dimension (*), intent (inout) :: xICTRL 
  real (kind = c_double), dimension (*), intent (inout) :: xPCTRL

  real, dimension (:), allocatable :: TCTRL
  real, dimension (:), allocatable :: ICTRL
  real, dimension (:), allocatable :: PCTRL

  integer          :: i
  double precision :: pi
 
  namelist /PHASE_CONTROL/ STAGE5, NFLOW, INTF, INTN, INTU, NATS, OLD, FREQ, FFAC, LIN, MID, COPT, DT,&
       TSTART, TEND, TOFF, SCALE, PMAX, CHIR, HIGH, RATS, CORE, CXD, POEM, BOOT, CURV, POLZ, WALL, TAUW,&
       THRES, MSTOP
  namelist /WAVEFORM_CONTROL/ TYPE, NCTRL
  namelist /TYPE1_WAVEFORM/ TCTRL, ICTRL, PCTRL
  namelist /TYPE2_WAVEFORM/ SSTART, SEND, SAMP, SPHA, SPVE, BACK
  namelist /TYPE3_WAVEFORM/ RPERIOD, RSTART, REND, RPHA
  
  open  (unit = 100, file = 'Inputs/Phase.nml', status = 'old')
  read  (unit = 100, nml  = PHASE_CONTROL) 
  close (unit = 100)
  
  open  (unit = 100, file = 'Inputs/Waveform.nml', status = 'old')
  read  (unit = 100, nml  = WAVEFORM_CONTROL)

  allocate (TCTRL (NCTRL))
  allocate (ICTRL (NCTRL))
  allocate (PCTRL (NCTRL))

  read  (unit = 100, nml  = TYPE1_WAVEFORM)
  read  (unit = 100, nml  = TYPE2_WAVEFORM)
  read  (unit = 100, nml  = TYPE3_WAVEFORM)

  close (unit = 100)
  
  pi = 4.*atan(1.)
  do i = 1, NCTRL
     xTCTRL (i) = TCTRL (i)
     xICTRL (i) = ICTRL (i)
     xPCTRL (i) = PCTRL (i) * pi
  enddo

  SPHA = SPHA * pi;
  RPHA = RPHA * pi;
 
  deallocate (TCTRL)
  deallocate (ICTRL)
  deallocate (PCTRL)

endsubroutine namelistRead
