! NameListRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read WindowFind namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListRead (ISTART, IEND, DI) &
     bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  real (kind = c_double), intent (inout) :: ISTART
  real (kind = c_double), intent (inout) :: IEND
  real (kind = c_double), intent (inout) :: DI
  
  namelist /WINDOWFIND_CONTROL/ ISTART, IEND, DI
  
  open  (unit = 100, file = 'Inputs/Window.nml', status = 'old')
  read  (unit = 100, nml  = WINDOWFIND_CONTROL)
  close (unit = 100)

endsubroutine namelistRead
