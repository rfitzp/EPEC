! Function to read PHASE namelist

subroutine NameListRead (NFLOW, STAGE2, INTF, INTN, INTU, OLD, DT, TIME, NCTRL, xTCTRL, xICTRL, xPCTRL) &
     bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  integer (kind = c_int),    intent (inout) :: STAGE2
  integer (kind = c_int),    intent (inout) :: NFLOW
  integer (kind = c_int),    intent (inout) :: INTF
  integer (kind = c_int),    intent (inout) :: INTN
  integer (kind = c_int),    intent (inout) :: INTU
  integer (kind = c_int),    intent (inout) :: OLD
  real    (kind = c_double), intent (inout) :: DT
  real    (kind = c_double), intent (inout) :: TIME
  integer (kind = c_int),    intent (inout) :: NCTRL
  
  real (kind = c_double), dimension (*), intent (inout) :: xTCTRL 
  real (kind = c_double), dimension (*), intent (inout) :: xICTRL 
  real (kind = c_double), dimension (*), intent (inout) :: xPCTRL

  real, dimension (:), allocatable :: TCTRL
  real, dimension (:), allocatable :: ICTRL
  real, dimension (:), allocatable :: PCTRL

  integer          :: i
  double precision :: pi
 
  namelist /PhaseInputs/ STAGE2, NFLOW, INTF, INTN, INTU, OLD, DT, TIME, NCTRL
  namelist /WaveForm/    TCTRL, ICTRL, PCTRL
  
  open  (unit = 100, file = 'namelist.txt', status = 'old')
  read  (unit = 100, nml = PhaseInputs) 
  close (unit = 100)
  
  allocate (TCTRL (NCTRL))
  allocate (ICTRL (NCTRL))
  allocate (PCTRL (NCTRL))

  open  (unit = 100, file = 'namelist.txt', status = 'old')
  read  (unit = 100, nml = WaveForm) 
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
