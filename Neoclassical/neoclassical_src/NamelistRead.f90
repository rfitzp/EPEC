! NameListRead.f90

! ######################################
! Function to read NEOCLASSICAL namelist
! ######################################

subroutine NameListRead (IMPURITY, NEUTRAL, FREQ, INTP, CHI, NN, LN, SVN, YN, EN, TIME, COULOMB) &
     bind (c, name = 'NameListRead')
  
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none
  
  real    (kind = c_double), intent (inout) :: CHI
  integer (kind = c_int),    intent (inout) :: IMPURITY
  integer (kind = c_int),    intent (inout) :: NEUTRAL
  integer (kind = c_int),    intent (inout) :: FREQ
  real    (kind = c_double), intent (inout) :: NN
  real    (kind = c_double), intent (inout) :: LN
  real    (kind = c_double), intent (inout) :: SVN
  real    (kind = c_double), intent (inout) :: YN
  real    (kind = c_double), intent (inout) :: EN
  real    (kind = c_double), intent (inout) :: TIME
  real    (kind = c_double), intent (inout) :: COULOMB
  integer (kind = c_int),    intent (inout) :: INTP
  
  namelist /NeoclassicalInputs/ CHI, IMPURITY, NEUTRAL, FREQ, NN, LN, SVN, YN, EN, TIME, COULOMB, INTP
  
  open  (unit = 100, file = 'namelist.txt', status = 'old')
  read  (unit = 100, nml = NeoclassicalInputs) 
  close (unit = 100)
  
endsubroutine namelistRead
