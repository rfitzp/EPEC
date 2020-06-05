! NameListRead.f90

! ######################################
! Function to read NEOCLASSICAL namelist
! ######################################

subroutine NameListRead (IMPURITY, NEUTRAL, FREQ, INTP, INTF, CHI, NN, LN, SVN, YN, EN, TIME, COULOMB) &
     bind (c, name = 'NameListRead')
  
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none
  
  integer (kind = c_int),    intent (inout) :: IMPURITY
  integer (kind = c_int),    intent (inout) :: NEUTRAL
  integer (kind = c_int),    intent (inout) :: FREQ
  integer (kind = c_int),    intent (inout) :: INTP
  integer (kind = c_int),    intent (inout) :: INTF
  real    (kind = c_double), intent (inout) :: CHI
  real    (kind = c_double), intent (inout) :: NN
  real    (kind = c_double), intent (inout) :: LN
  real    (kind = c_double), intent (inout) :: SVN
  real    (kind = c_double), intent (inout) :: YN
  real    (kind = c_double), intent (inout) :: EN
  real    (kind = c_double), intent (inout) :: TIME
  real    (kind = c_double), intent (inout) :: COULOMB
  
  namelist /NeoclassicalInputs/ IMPURITY, NEUTRAL, FREQ, INTP, INTF, CHI, NN, LN, SVN, YN, EN, TIME, COULOMB
  
  open  (unit = 100, file = 'namelist.txt', status = 'old')
  read  (unit = 100, nml = NeoclassicalInputs) 
  close (unit = 100)
  
endsubroutine namelistRead
