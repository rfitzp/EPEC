! NameListRead.f90

! ######################################
! Function to read NEOCLASSICAL namelist
! ######################################

subroutine NameListRead (IMPURITY, NEUTRAL, EXB, INTP, INTF, INTC, NTYPE, NN, LN, SVN, YN, EN, TIME, COULOMB, NSMOOTH, CATS) &
     bind (c, name = 'NameListRead')
  
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none
  
  integer (kind = c_int),    intent (inout) :: IMPURITY
  integer (kind = c_int),    intent (inout) :: NEUTRAL
  integer (kind = c_int),    intent (inout) :: EXB
  integer (kind = c_int),    intent (inout) :: INTP
  integer (kind = c_int),    intent (inout) :: INTF
  integer (kind = c_int),    intent (inout) :: INTC
  integer (kind = c_int),    intent (inout) :: NTYPE
  integer (kind = c_int),    intent (inout) :: NSMOOTH
  integer (kind = c_int),    intent (inout) :: CATS
  real    (kind = c_double), intent (inout) :: NN
  real    (kind = c_double), intent (inout) :: LN
  real    (kind = c_double), intent (inout) :: SVN
  real    (kind = c_double), intent (inout) :: YN
  real    (kind = c_double), intent (inout) :: EN
  real    (kind = c_double), intent (inout) :: TIME
  real    (kind = c_double), intent (inout) :: COULOMB
  
  namelist /NEOCLASSICAL_CONTROL/ IMPURITY, NEUTRAL, EXB, INTP, INTF, INTC, NTYPE, NN, LN, SVN, YN, EN, TIME, COULOMB, NSMOOTH, CATS
  
  open  (unit = 100, file = 'Inputs/Neoclassical.nml', status = 'old')
  read  (unit = 100, nml = NEOCLASSICAL_CONTROL) 
  close (unit = 100)
  
endsubroutine namelistRead
