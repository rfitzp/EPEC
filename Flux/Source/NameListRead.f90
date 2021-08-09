! NameListRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read FLUX namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListRead (INTG, NPSI, PACK, NTHETA, NNC, NTOR, H0, ACC, ETA, MMIN, MMAX,&
     PSILIM, TIME, PSIPED, NSMOOTH, PSIRAT, NEOANG) &
     bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  integer (kind = c_int),    intent (inout) :: INTG
  integer (kind = c_int),    intent (inout) :: NPSI
  integer (kind = c_int),    intent (inout) :: NTHETA
  integer (kind = c_int),    intent (inout) :: NNC
  integer (kind = c_int),    intent (inout) :: NTOR
  integer (kind = c_int),    intent (inout) :: MMIN
  integer (kind = c_int),    intent (inout) :: MMAX
  integer (kind = c_int),    intent (inout) :: NSMOOTH
  integer (kind = c_int),    intent (inout) :: NEOANG
  real    (kind = c_double), intent (inout) :: PACK
  real    (kind = c_double), intent (inout) :: H0
  real    (kind = c_double), intent (inout) :: ACC
  real    (kind = c_double), intent (inout) :: ETA
  real    (kind = c_double), intent (inout) :: PSILIM
  real    (kind = c_double), intent (inout) :: TIME
  real    (kind = c_double), intent (inout) :: PSIPED
  real    (kind = c_double), intent (inout) :: PSIRAT
    
  namelist /FLUX_CONTROL/ INTG, NPSI, PACK, NTHETA, NNC, NTOR, H0, ACC, ETA, MMIN, MMAX,&
       PSILIM, TIME, PSIPED, NSMOOTH, PSIRAT, NEOANG
  
  open  (unit = 100, file = 'Inputs/Flux.nml', status = 'old')
  read  (unit = 100, nml  = FLUX_CONTROL)
  close (unit = 100)

endsubroutine namelistRead
