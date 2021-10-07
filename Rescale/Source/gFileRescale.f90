! gFileRescale.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Function to perform Type 0 rescaling of gFile
!
! No rescaling
!
! Initial gFile in Inputs/gFile
! Rescaled gFile in Outputs/gFile
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine gFileRescaleType0 () bind (c, name = 'gFileRescaleType0')

  use Function_Defs_0

  implicit none

  character (len = 100) :: string
  
  integer :: i, j, i3, NRBOX, NZBOX, NBOUND, NLIM 
  
  double precision :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
  double precision :: RLEFT, RRIGHT, ZLOW, ZHIGH, RMAX, ZMAX
  
  double precision, dimension (:), allocatable :: T, P, TTp, Pp, Q, RBOUND, ZBOUND, RLIM, ZLIM
  
  double precision, dimension (:, :), allocatable :: Psi
 
  ! ................
  ! Read input gFile
  ! ................

  call ReadgFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
       B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)

  ! ..................
  ! Write output gFile
  ! ..................

  call WritegFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
       B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)
   
  ! ........
  ! Clean up
  ! ........
  
  deallocate (T)
  deallocate (P)
  deallocate (TTp)
  deallocate (Pp)
  deallocate (Q)
  deallocate (Psi)
  deallocate (RBOUND)
  deallocate (ZBOUND)
  deallocate (RLIM)
  deallocate (ZLIM)
  
end subroutine gFileRescaleType0

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Function to perform Type I rescaling of gFile
!
! Rescaling such that current increases by factor A
! at constant beta and q95
!
! Initial gFile in Inputs/gFile
! Rescaled gFile in Outputs/gFile
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine gFileRescaleTypeI (A, OPOINT, XPOINT) bind (c, name = 'gFileRescaleTypeI')

  use Function_Defs_0

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  real    (kind = c_double), intent (inout) :: A
  integer (kind = c_int),    intent (inout) :: OPOINT
  integer (kind = c_int),    intent (inout) :: XPOINT

  character (len = 100) :: string
  
  integer :: i, j, i3, NRBOX, NZBOX, NBOUND, NLIM 
  
  double precision :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
  double precision :: RLEFT, RRIGHT, ZLOW, ZHIGH, RMAX, ZMAX
  double precision :: dR, dZ, RO, ZO, RX, ZX
 
  double precision, dimension (:), allocatable :: T, P, TTp, Pp, Q, RBOUND, ZBOUND, RLIM, ZLIM, RR, ZZ
  
  double precision, dimension (:, :), allocatable :: Psi
 
  ! ................
  ! Read input gFile
  ! ................

  call ReadgFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
       B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)

  allocate (RR (NRBOX))
  allocate (ZZ (NZBOX))

  ! ..........
  ! Rescale B0
  ! ..........
  B0 = A * B0

  ! ...............
  ! Setup R, Z grid
  ! ...............
 
  RLEFT  = RBOXLFT 
  RRIGHT = RBOXLFT + RBOXLEN
  ZLOW   = - ZBOXLEN /2. + ZOFF
  ZHIGH  =   ZBOXLEN /2. + ZOFF

  RMAX = dble (NRBOX-1)
  do i = 1, NRBOX
     RR (i) = RLEFT + (RRIGHT - RLEFT) * dble (i-1) /RMAX
  enddo

  ZMAX = dble (NZBOX-1)
  do j = 1, NZBOX
     ZZ (j) =  ZLOW + (ZHIGH - ZLOW) * dble (j-1) /ZMAX
  enddo

  ! ...................
  ! Rescale equilibrium
  ! ...................
  
  do i = 1, NRBOX

     do j = 1, NZBOX
        Psi (i, j) = A * Psi (i, j)
     enddo
     
     T   (i) = A   * T   (i)
     TTp (i) = A   * TTp (i)
     P   (i) = A*A * P   (i)
     Pp  (i) = A   * Pp  (i)
     Q   (i) =       Q   (i)

  enddo

  PSIAXIS  = A * PSIAXIS
  PSIBOUND = A * PSIBOUND
  CURRENT  = A * CURRENT

  ! ............
  ! Find O-point
  ! ............

  dR = 1.e-2
  dZ = 1.e-2

  RO = RAXIS
  ZO = ZAXIS

  if (OPOINT /= 0) then
     call FindOXPoint (RO, ZO, NRBOX, NZBOX, RR, ZZ, Psi, dR, dZ, PSIAXIS)
  end if

  RAXIS = RO
  ZAXIS = ZO

  ! ...................................
  ! Find X-point (assume lower X-point)
  ! ...................................

  ZX = ZAXIS
  do i = 1, NBOUND
     if (ZBOUND (i) < ZX) then
        ZX = ZBOUND (i)
        RX = RBOUND (i)
     endif
  enddo

  if (XPOINT /= 0) then
     call FindOXPoint (RX, ZX, NRBOX, NZBOX, RR, ZZ, Psi, dR, dZ, PSIBOUND)
  end if
 
  ! ..................
  ! Write output gFile
  ! ..................

  call WritegFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
       B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)
   
  ! ........
  ! Clean up
  ! ........
  
  deallocate (T)
  deallocate (P)
  deallocate (TTp)
  deallocate (Pp)
  deallocate (Q)
  deallocate (Psi)
  deallocate (RBOUND)
  deallocate (ZBOUND)
  deallocate (RLIM)
  deallocate (ZLIM)
  deallocate (RR)
  deallocate (ZZ)
  
end subroutine gFileRescaleTypeI

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Function to perform Type II rescaling of gFile
!
! Rescaling such that size increases by factor A
! at constant beta and q95
!
! Initial gFile in Inputs/gFile
! Rescaled gFile in Outputs/gFile
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine gFileRescaleTypeII (A, OPOINT, XPOINT) bind (c, name = 'gFileRescaleTypeII')

  use Function_Defs_0

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  real    (kind = c_double), intent (inout) :: A
  integer (kind = c_int),    intent (inout) :: OPOINT
  integer (kind = c_int),    intent (inout) :: XPOINT

  character (len = 100) :: string
  
  integer :: i, j, i3, NRBOX, NZBOX, NBOUND, NLIM
  
  double precision :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
  double precision :: RLEFT, RRIGHT, ZLOW, ZHIGH, RMAX, ZMAX
  double precision :: dR, dZ, RO, ZO, RX, ZX
 
  double precision, dimension (:), allocatable :: T, P, TTp, Pp, Q, RBOUND, ZBOUND, RLIM, ZLIM, RR, ZZ
  
  double precision, dimension (:, :), allocatable :: Psi
 
  ! ................
  ! Read input gFile
  ! ................

  call ReadgFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
       B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)

  allocate (RR (NRBOX))
  allocate (ZZ (NZBOX))

  ! ...............
  ! Setup R, Z grid
  ! ...............

  RBOXLEN = A * RBOXLEN
  ZBOXLEN = A * ZBOXLEN
  R0      = A * R0
  RBOXLFT = A * RBOXLFT
  ZOFF    = A * ZOFF
  RAXIS   = A * RAXIS
  ZAXIS   = A * ZAXIS
  
  RLEFT  = RBOXLFT 
  RRIGHT = RBOXLFT + RBOXLEN
  ZLOW   = - ZBOXLEN /2. + ZOFF
  ZHIGH  =   ZBOXLEN /2. + ZOFF

  RMAX = dble (NRBOX-1)
  do i = 1, NRBOX
     RR (i) = RLEFT + (RRIGHT - RLEFT) * dble (i-1) /RMAX
  enddo

  ZMAX = dble (NZBOX-1)
  do j = 1, NZBOX
     ZZ (j) =  ZLOW + (ZHIGH - ZLOW) * dble (j-1) /ZMAX
  enddo

  do i = 1, NBOUND
     RBOUND (i) = A * RBOUND (i)
     ZBOUND (i) = A * ZBOUND (i)
  end do

  do i = 1, NLIM
     RLIM (i) = A * RLIM (i)
     ZLIM (i) = A * ZLIM (i)
  end do

   ! ...................
  ! Rescale equilibrium
  ! ...................
  
  do i = 1, NRBOX

     do j = 1, NZBOX
        Psi (i, j) = A*A * Psi (i, j)
     enddo

     T   (i) = A * T   (i)
     TTp (i) =     TTp (i) 
     P   (i) =     P   (i)
     Pp  (i) =     Pp  (i) /A/A
     Q   (i) =     Q   (i)
 
  enddo

  PSIAXIS  = A*A * PSIAXIS
  PSIBOUND = A*A * PSIBOUND
  CURRENT  = A   * CURRENT

  ! ............
  ! Find O-point
  ! ............

  dR = 1.e-2
  dZ = 1.e-2

  RO = RAXIS
  ZO = ZAXIS

  if (OPOINT /= 0) then
     call FindOXPoint (RO, ZO, NRBOX, NZBOX, RR, ZZ, Psi, dR, dZ, PSIAXIS)
  end if

  RAXIS = RO
  ZAXIS = ZO

  ! ...................................
  ! Find X-point (assume lower X-point)
  ! ...................................

  ZX = ZAXIS
  do i = 1, NBOUND
     if (ZBOUND (i) < ZX) then
        ZX = ZBOUND (i)
        RX = RBOUND (i)
     endif
  enddo

  if (XPOINT /= 0) then
     call FindOXPoint (RX, ZX, NRBOX, NZBOX, RR, ZZ, Psi, dR, dZ, PSIBOUND)
  end if
 
  ! ..................
  ! Write output gFile
  ! ..................

  call WritegFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
       B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)
   
  ! ........
  ! Clean up
  ! ........
  
  deallocate (T)
  deallocate (P)
  deallocate (TTp)
  deallocate (Pp)
  deallocate (Q)
  deallocate (Psi)
  deallocate (RBOUND)
  deallocate (ZBOUND)
  deallocate (RLIM)
  deallocate (ZLIM)
  deallocate (RR)
  deallocate (ZZ)
  
end subroutine gFileRescaleTypeII

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Function to perform Type III rescaling of gFile
!
! Rescaling such that pressure shifted by constant
! factor at constant beta and q95
!
! Initial gFile in Inputs/gFile
! Rescaled gFile in Outputs/gFile
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine gFileRescaleTypeIII (PSHIFT) bind (c, name = 'gFileRescaleTypeIII')

  use Function_Defs_0

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  real    (kind = c_double), intent (inout) :: PSHIFT

  character (len = 100) :: string
  
  integer :: i, j, i3, NRBOX, NZBOX, NBOUND, NLIM, i95
  
  double precision :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
 
  double precision, dimension (:), allocatable :: T, P, TTp, Pp, Q, RBOUND, ZBOUND, RLIM, ZLIM
  
  double precision, dimension (:, :), allocatable :: Psi
 
  ! ................
  ! Read input gFile
  ! ................

  call ReadgFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
       B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)

  ! ...................
  ! Rescale equilibrium
  ! ...................
  
  do i = 1, NRBOX

   P (i) = P (i) + PSHIFT

  enddo
 
  ! ..................
  ! Write output gFile
  ! ..................

  call WritegFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
       B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)
   
  ! ........
  ! Clean up
  ! ........
  
  deallocate (T)
  deallocate (P)
  deallocate (TTp)
  deallocate (Pp)
  deallocate (Q)
  deallocate (Psi)
  deallocate (RBOUND)
  deallocate (ZBOUND)
  deallocate (RLIM)
  deallocate (ZLIM)
  
end subroutine gFileRescaleTypeIII

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Function to perform Type 7 rescaling of gFile
!
! Rescaling such that q_95 = q_95_target but
! vacuum toroidal magnetic field remains constant
!
! Initial gFile in Inputs/gFile
! Rescaled gFile in Outputs/gFile
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine gFileRescaleType7 (Q95_NEW, OPOINT, XPOINT, q95_old, a1) bind (c, name = 'gFileRescaleType7')

  use Function_Defs_0

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  real    (kind = c_double), intent (inout) :: Q95_NEW
  integer (kind = c_int),    intent (inout) :: OPOINT
  integer (kind = c_int),    intent (inout) :: XPOINT
  real    (kind = c_double), intent (inout) :: q95_old
  real    (kind = c_double), intent (inout) :: a1

  character (len = 100) :: string
  
  integer :: i, j, i3, NRBOX, NZBOX, NBOUND, NLIM, i95
  
  double precision :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
  double precision :: RLEFT, RRIGHT, ZLOW, ZHIGH, RMAX, ZMAX
  double precision :: x95, q95, T95, Tedge, a2, Tnew, Told, dR, dZ, RO, ZO, RX, ZX
 
  double precision, dimension (:), allocatable :: T, P, TTp, Pp, Q, RBOUND, ZBOUND, RLIM, ZLIM, RR, ZZ
  
  double precision, dimension (:, :), allocatable :: Psi
 
  ! ................
  ! Read input gFile
  ! ................

  call ReadgFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
       B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)

  allocate (RR (NRBOX))
  allocate (ZZ (NZBOX))

  ! ...............
  ! Setup R, Z grid
  ! ...............
 
  RLEFT  = RBOXLFT 
  RRIGHT = RBOXLFT + RBOXLEN
  ZLOW   = - ZBOXLEN /2. + ZOFF
  ZHIGH  =   ZBOXLEN /2. + ZOFF

  RMAX = dble (NRBOX-1)
  do i = 1, NRBOX
     RR (i) = RLEFT + (RRIGHT - RLEFT) * dble (i-1) /RMAX
  enddo

  ZMAX = dble (NZBOX-1)
  do j = 1, NZBOX
     ZZ (j) =  ZLOW + (ZHIGH - ZLOW) * dble (j-1) /ZMAX
  enddo

  ! ..................
  ! Determine old q_95
  ! ..................

  x95     = 0.95 * dble (NRBOX) + 1.
  i95     = int (x95)
  q95     = (dble (i95) - x95) * Q (i95+1) + (x95 + 1. - dble (i95)) * Q (i95)
  T95     = (dble (i95) - x95) * T (i95+1) + (x95 + 1. - dble (i95)) * T (i95)
  Tedge   = T (NRBOX)
  q95_old = q95

  ! ...................
  ! Rescale equilibrium
  ! ...................

  a2 = (Q95_NEW*Q95_NEW /q95/q95 - 1.) * T95*T95
  
  do i = 1, NRBOX
     Told = T (i)
     if (Told > 0.) then
        Tnew =   sqrt (Told*Told + a2)
     else
        Tnew = - sqrt (Told*Told + a2)
     endif
     T (i) = Tnew

     Q (i) = Q (i) * Tnew /Told
   enddo

  a1 = Tedge /T (NRBOX)
 
  do i = 1, NRBOX

     do j = 1, NZBOX
        Psi (i, j) = a1 * Psi (i, j)
     enddo

     T  (i)  = a1    * T   (i)
     TTp (i) = a1    * TTp (i)
     P  (i)  = a1*a1 * P   (i) 
     Pp (i)  = a1    * Pp  (i)
     Q  (i)  =         Q   (i)

  enddo

  PSIAXIS  = a1 * PSIAXIS
  PSIBOUND = a1 * PSIBOUND
  CURRENT  = a1 * CURRENT

  ! ..................
  ! Determine new q_95
  ! ..................

  x95 = 0.95 * dble (NRBOX) + 1.
  i95 = int (x95)
  q95 = (dble (i95) - x95) * Q (i95+1) + (x95 + 1. - dble (i95)) * Q (i95)
  T95 = (dble (i95) - x95) * T (i95+1) + (x95 + 1. - dble (i95)) * T (i95)

  ! ............
  ! Find O-point
  ! ............

  dR = 1.e-2
  dZ = 1.e-2

  RO = RAXIS
  ZO = ZAXIS

  if (OPOINT /= 0) then
     call FindOXPoint (RO, ZO, NRBOX, NZBOX, RR, ZZ, Psi, dR, dZ, PSIAXIS)
  end if

  RAXIS = RO
  ZAXIS = ZO

  ! ...................................
  ! Find X-point (assume lower X-point)
  ! ...................................

  ZX = ZAXIS
  do i = 1, NBOUND
     if (ZBOUND (i) < ZX) then
        ZX = ZBOUND (i)
        RX = RBOUND (i)
     endif
  enddo

  if (XPOINT /= 0) then
     call FindOXPoint (RX, ZX, NRBOX, NZBOX, RR, ZZ, Psi, dR, dZ, PSIBOUND)
  end if
 
  ! ..................
  ! Write output gFile
  ! ..................

  call WritegFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
       B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)
   
  ! ........
  ! Clean up
  ! ........
  
  deallocate (T)
  deallocate (P)
  deallocate (TTp)
  deallocate (Pp)
  deallocate (Q)
  deallocate (Psi)
  deallocate (RBOUND)
  deallocate (ZBOUND)
  deallocate (RLIM)
  deallocate (ZLIM)
  deallocate (RR)
  deallocate (ZZ)
  
end subroutine gFileRescaleType7

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read Rescale namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListRead (TYPE, SCALE, PSHIFT, WSHIFT, Q95, OPOINT, XPOINT) bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  integer (kind = c_int),    intent (inout) :: TYPE
  integer (kind = c_int),    intent (inout) :: OPOINT
  integer (kind = c_int),    intent (inout) :: XPOINT
  real    (kind = c_double), intent (inout) :: SCALE
  real    (kind = c_double), intent (inout) :: PSHIFT
  real    (kind = c_double), intent (inout) :: WSHIFT
  real    (kind = c_double), intent (inout) :: Q95
  
  namelist /RESCALE_CONTROL/ TYPE, SCALE, PSHIFT, WSHIFT, Q95, OPOINT, XPOINT
  
  open  (unit = 100, file = 'Inputs/Rescale.nml', status = 'old')
  read  (unit = 100, nml  = RESCALE_CONTROL)
  close (unit = 100)

end subroutine namelistRead

! %%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to read gFile
! %%%%%%%%%%%%%%%%%%%%%%%%

subroutine ReadgFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
     B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)

  implicit none

  character (len = 100) :: string

  integer :: i, j, i3, NRBOX, NZBOX, NBOUND, NLIM
  
  double precision :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
 
  double precision, dimension (:), allocatable :: T, P, TTp, Pp, Q, RBOUND, ZBOUND, RLIM, ZLIM
  
  double precision, dimension (:, :), allocatable :: Psi

  open (unit = 100, file = 'Inputs/gFile', status = 'old')
  
  read (100, '(a48, 3i4)') string,  i3,      NRBOX,   NZBOX
  read (100, '(5e16.9  )') RBOXLEN, ZBOXLEN, R0,      RBOXLFT,  ZOFF
  read (100, '(5e16.9  )') RAXIS,   ZAXIS,   PSIAXIS, PSIBOUND, B0
  read (100, '(5e16.9  )') CURRENT, zero,    zero,    zero,     zero
  read (100, '(5e16.9  )') zero,    zero,    zero,    zero,     zero

  allocate (T   (NRBOX))
  allocate (P   (NRBOX))
  allocate (TTp (NRBOX))
  allocate (Pp  (NRBOX))
  allocate (Q   (NRBOX))
  allocate (Psi (NRBOX, NZBOX))
  
  read (100, '(5e16.9)') (T   (i), i = 1, NRBOX)
  read (100, '(5e16.9)') (P   (i), i = 1, NRBOX)
  read (100, '(5e16.9)') (TTp (i), i = 1, NRBOX)
  read (100, '(5e16.9)') (Pp  (i), i = 1, NRBOX)
  
  read (100, '(5e16.9)') ((Psi  (i, j), i = 1, NRBOX), j = 1, NZBOX)

  read (100, '(5e16.9)') (Q (i), i = 1, NRBOX)

  read (100, '(2i5)') NBOUND, NLIM

  allocate (RBOUND (NBOUND))
  allocate (ZBOUND (NBOUND))
  allocate (RLIM   (NLIM))
  allocate (ZLIM   (NLIM))

  read (100, '(5e16.9)') (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
  read (100, '(5e16.9)') (RLIM   (i), ZLIM   (i), i = 1, NLIM)

  close (unit = 100)

end subroutine ReadgFile

! %%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to write gFile
! %%%%%%%%%%%%%%%%%%%%%%%%%

subroutine WritegFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
     B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)

  implicit none

  character (len = 100) :: string

  integer :: i, j, i3, NRBOX, NZBOX, NBOUND, NLIM
  
  double precision :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
 
  double precision, dimension (:), allocatable :: T, P, TTp, Pp, Q, RBOUND, ZBOUND, RLIM, ZLIM
  
  double precision, dimension (:, :), allocatable :: Psi

  open (unit = 100, file = 'Outputs/gFile', status = 'unknown')

  zero = 0.
  
  write (100, '(a48, 3i4)') string,  i3,      NRBOX,    NZBOX
  write (100, '(5e16.9  )') RBOXLEN, ZBOXLEN, R0,       RBOXLFT,  ZOFF
  write (100, '(5e16.9  )') RAXIS,   ZAXIS,   PSIAXIS,  PSIBOUND, B0
  write (100, '(5e16.9  )') CURRENT, PSIAXIS, zero,     RAXIS,    zero
  write (100, '(5e16.9  )') ZAXIS,   zero,    PSIBOUND, zero,     zero
  write (100, '(5e16.9)')  (T   (i), i = 1, NRBOX)
  write (100, '(5e16.9)')  (P   (i), i = 1, NRBOX)
  write (100, '(5e16.9)')  (TTp (i), i = 1, NRBOX)
  write (100, '(5e16.9)')  (Pp  (i), i = 1, NRBOX)
  write (100, '(5e16.9)')  ((Psi  (i, j), i = 1, NRBOX), j = 1, NZBOX)
  write (100, '(5e16.9)')  (Q (i), i = 1, NRBOX)
  write (100, '(2i5)')      NBOUND, NLIM
  write (100, '(5e16.9)')  (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
  write (100, '(5e16.9)')  (RLIM   (i), ZLIM   (i), i = 1, NLIM)

  close (unit = 100)
  
end subroutine WritegFile

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to return interpolated Psi
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double precision function GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, Psi)

  implicit none

  integer :: NRBOX, NZBOX, i, j

  double precision :: R, Z, x, y

  double precision, dimension (:), allocatable :: RR, ZZ
  
  double precision, dimension (:, :), allocatable :: Psi

  i = 1 + int ((R - RR (1)) /(RR (2) - RR (1)))
  j = 1 + int ((Z - ZZ (1)) /(ZZ (2) - ZZ (1)))

  x = (R - RR (i)) /(RR (2) - RR (1))
  y = (Z - ZZ (j)) /(ZZ (2) - ZZ (1))

  GetPsi = Psi (i, j) * (1.-x)*(1.-y) + Psi (i+1, j) * x*(1.-y) + Psi (i, j+1) * (1.-x)*y + Psi(i+1, j+1) * x*y

end function GetPsi

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to return interpolated dPsi/dR
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double precision function GetPsiR (R, Z, NRBOX, NZBOX, RR, ZZ, Psi, dR)

  use Function_Defs_1

  integer :: NRBOX, NZBOX

  double precision :: R, Z, dR, R1, R2 

  double precision, dimension (:), allocatable :: RR, ZZ
  
  double precision, dimension (:, :), allocatable :: Psi

  R1 = R - dR
  R2 = R + dR

  GetPsiR = (GetPsi (R2, Z, NRBOX, NZBOX, RR, ZZ, Psi) - GetPsi (R1, Z, NRBOX, NZBOX, RR, ZZ, Psi)) /2./dR

end function GetPsiR

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to return interpolated dPsi/dZ
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double precision function GetPsiZ (R, Z, NRBOX, NZBOX, RR, ZZ, Psi, dZ)

  use Function_Defs_1
  
  implicit none

  integer :: NRBOX, NZBOX

  double precision :: R, Z, dZ, Z1, Z2 

  double precision, dimension (:), allocatable :: RR, ZZ
  
  double precision, dimension (:, :), allocatable :: Psi

  Z1 = Z - dZ
  Z2 = Z + dZ

  GetPsiZ = (GetPsi (R, Z2, NRBOX, NZBOX, RR, ZZ, Psi) - GetPsi (R, Z1, NRBOX, NZBOX, RR, ZZ, Psi)) /2./dZ

end function GetPsiZ

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to return interpolated d^2Psi/dR^2
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double precision function GetPsiRR (R, Z, NRBOX, NZBOX, RR, ZZ, Psi, dR)

  use Function_Defs_1
  
  implicit none

  integer :: NRBOX, NZBOX

  double precision :: R, Z, dR, R1, R2

  double precision, dimension (:), allocatable :: RR, ZZ
  
  double precision, dimension (:, :), allocatable :: Psi

  R1 = R - dR
  R2 = R + dR

  GetPsiRR = (GetPsi (R2, Z, NRBOX, NZBOX, RR, ZZ, Psi) - 2. * GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, Psi)&
       + GetPsi (R1, Z, NRBOX, NZBOX, RR, ZZ, Psi)) /dR/dR

end function GetPsiRR

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to return interpolated d^2Psi/dZ^2
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double precision function GetPsiZZ (R, Z, NRBOX, NZBOX, RR, ZZ, Psi, dZ)

  use Function_Defs_1
  
  implicit none
  
  integer :: NRBOX, NZBOX

  double precision :: R, Z, dZ, Z1, Z2

  double precision, dimension (:), allocatable :: RR, ZZ
  
  double precision, dimension (:, :), allocatable :: Psi

  Z1 = Z - dZ
  Z2 = Z + dZ

  GetPsiZZ = (GetPsi (R, Z2, NRBOX, NZBOX, RR, ZZ, Psi) - 2. * GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, Psi)&
       +  GetPsi (R, Z1, NRBOX, NZBOX, RR, ZZ, Psi)) /dZ/dZ

end function GetPsiZZ

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to return interpolated d^2Psi/dRdZ
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double precision function GetPsiRZ (R, Z, NRBOX, NZBOX, RR, ZZ, Psi, dR)

  use Function_Defs_1
  
  implicit none

  integer :: NRBOX, NZBOX

  double precision :: R, Z, dR, R1, R2, dZ, Z1, Z2

  double precision, dimension (:), allocatable :: RR, ZZ
  
  double precision, dimension (:, :), allocatable :: Psi

  R1 = R - dR
  R2 = R + dR
  Z1 = Z - dR
  Z2 = Z + dR

  GetPsiRZ = (GetPsi (R2, Z2, NRBOX, NZBOX, RR, ZZ, Psi) - GetPsi (R1, Z2, NRBOX, NZBOX, RR, ZZ, Psi)&
       -  GetPsi (R2, Z1, NRBOX, NZBOX, RR, ZZ, Psi) + GetPsi (R1, Z1, NRBOX, NZBOX, RR, ZZ, Psi)) /4./dR/dR

end function GetPsiRZ

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to find O- and X-points
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine FindOXPoint (R, Z, NRBOX, NZBOX, RR, ZZ, Psi, dR, dZ, p)

  use Function_Defs_1
  use Function_Defs_2
  
  implicit none

  integer :: NRBOX, NZBOX, i

  double precision :: R, Z, dR, dZ, p, pr, pz, prr, pzz, prz, det

  double precision, dimension (:), allocatable :: RR, ZZ
  
  double precision, dimension (:, :), allocatable :: Psi

  do i = 1, 50
     
     pr  = GetPsiR  (R, Z, NRBOX, NZBOX, RR, ZZ, Psi, dR)
     pz  = GetPsiZ  (R, Z, NRBOX, NZBOX, RR, ZZ, Psi, dZ)
     prr = GetPsiRR (R, Z, NRBOX, NZBOX, RR, ZZ, Psi, dR)
     pzz = GetPsiZZ (R, Z, NRBOX, NZBOX, RR, ZZ, Psi, dZ)
     prz = GetPsiRZ (R, Z, NRBOX, NZBOX, RR, ZZ, Psi, dR)

     det = prr*pzz - prz*prz

     R = R + (prz*pz - pzz*pr) /det
     Z = Z + (prz*pr - prr*pz) /det

     p = GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, Psi)
     
  enddo

end subroutine FindOXPoint

