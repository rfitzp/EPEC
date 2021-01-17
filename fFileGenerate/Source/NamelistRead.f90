! NameListRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read fFileGenerate namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListRead (FLUX_NTOR, FLUX_MMIN, FLUX_MMAX,&
     NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN,&
     TSTART, TEND, DT) &
     bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  integer (kind = c_int),    intent (inout) :: FLUX_NTOR
  integer (kind = c_int),    intent (inout) :: FLUX_MMIN
  integer (kind = c_int),    intent (inout) :: FLUX_MMAX
  integer (kind = c_int),    intent (inout) :: NEO_IMPURITY
  integer (kind = c_int),    intent (inout) :: NEO_NEUTRAL
  integer (kind = c_int),    intent (inout) :: NEO_FREQ
  integer (kind = c_int),    intent (inout) :: NEO_NTYPE
  real    (kind = c_double), intent (inout) :: NEO_NN
  real    (kind = c_double), intent (inout) :: NEO_LN
  real    (kind = c_double), intent (inout) :: NEO_YN
  real    (kind = c_double), intent (inout) :: TSTART
  real    (kind = c_double), intent (inout) :: TEND
  real    (kind = c_double), intent (inout) :: DT
  
  namelist /FFILEGENERATE_CONTROL/ FLUX_NTOR, FLUX_MMIN, FLUX_MMAX,&
       NEO_IMPURITY, NEO_NEUTRAL, NEO_FREQ, NEO_NTYPE, NEO_NN, NEO_LN, NEO_YN,&
       TSTART, TEND, DT
  
  open  (unit = 100, file = 'Inputs/fFile.nml', status = 'old')
  read  (unit = 100, nml  = FFILEGENERATE_CONTROL)
  close (unit = 100)

endsubroutine namelistRead
