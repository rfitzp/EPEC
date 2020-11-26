! Rescale.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Program to rescale plasma equilibrium gFile
!
! Rescaling such that q_95 = q_95_target but
! vacuum toroidal magnetic field remains constant
!
! Initial gFile in Inputs/gFile
! Rescaled gFile in Outputs/gFile
! Target q_95 in Inputs/q_95
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program Rescale

  implicit none

  interface
     subroutine ReadgFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
          B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)
       
       character (len = 100) :: string
       
       integer :: i3, NRBOX, NZBOX, NBOUND, NLIM
       
       double precision :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
       
       double precision, dimension (:), allocatable :: T, P, TTp, Pp, Q, RBOUND, ZBOUND, RLIM, ZLIM
       
       double precision, dimension (:, :), allocatable :: PSI
       
     end subroutine ReadgFile
  end interface

  interface
     subroutine WritegFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
          B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)

       character (len = 100) :: string
       
       integer :: i3, NRBOX, NZBOX, NBOUND, NLIM
       
       double precision :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
       
       double precision, dimension (:), allocatable :: T, P, TTp, Pp, Q, RBOUND, ZBOUND, RLIM, ZLIM
       
       double precision, dimension (:, :), allocatable :: PSI
    
     end subroutine WritegFile
  end interface

  interface
     subroutine FindOXPoint (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dR, dZ, p)

       integer :: NRBOX, NZBOX

       double precision :: R, Z, dR, dZ, p
       
       double precision, dimension (:), allocatable :: RR, ZZ
       
       double precision, dimension (:, :), allocatable :: PSI
    
     end subroutine FindOXPoint
  end interface

  character (len = 100) :: string
  
  integer :: i, j, i3, NRBOX, NZBOX, NBOUND, NLIM, i95
  
  double precision :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
  double precision :: RLEFT, RRIGHT, ZLOW, ZHIGH, RMAX, ZMAX
  double precision :: x95, q95, T95, Tedge, q95_new, q95_old, a2, a1, Tnew, Told, qnew, qold
  double precision :: Pold, Pnew, Psiold, Psinew, dR, dZ, RO, ZO, RX, ZX
 
  double precision, dimension (:), allocatable :: T, P, TTp, Pp, Q, RBOUND, ZBOUND, RLIM, ZLIM, RR, ZZ
  
  double precision, dimension (:, :), allocatable :: PSI
 
  print *, ''
  print *, 'Program Rescale'
  print *, ''
  
  ! ................
  ! Read input gFile
  ! ................

  call ReadgFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
       B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)

  allocate (RR (NRBOX))
  allocate (ZZ (NZBOX))

  ! ...............
  ! Setup R, Z grid
  ! ...............
 
  RLEFT  = RBOXLFT 
  RRIGHT = RBOXLFT + RBOXLEN
  ZLOW   = -ZBOXLEN /2.
  ZHIGH  =  ZBOXLEN /2.
  
  RMAX = dble (NRBOX-1)
  do i = 1, NRBOX
     RR (i) = RLEFT + (RRIGHT - RLEFT) * dble (i-1) /RMAX
  enddo

  ZMAX = dble (NZBOX-1)
  do j = 1, NZBOX
     ZZ (j) =  ZLOW + (ZHIGH - ZLOW) * dble (j-1) /ZMAX
  enddo

  ! ..........................
  ! Output initial equilibrium
  ! ..........................
   
  open  (unit = 101, file = 'Outputs/Box.txt')
  write (101, '(4e17.9)') RLEFT, ZLOW, RRIGHT, ZHIGH
  close (unit = 101)

  open  (unit = 101, file = 'Outputs/Axis.txt')
  write (101, '(4e17.9)') RAXIS, ZAXIS
  close (unit = 101)

  open  (unit = 101, file = 'Outputs/Psiaxis.txt')
  write (101, '(4e17.9)') PSIAXIS, PSIBOUND
  close (unit = 101)
 
  open  (unit = 101, file = 'Outputs/R.txt')
  do i = 1, NRBOX
     write (101, '(1e17.9)') RR (i)
  enddo
  close (unit = 101)

  open (unit = 101, file = 'Outputs/Z.txt')
  do j = 1, NZBOX
     write (101, '(1e17.9)') ZZ (j)
  enddo
  close (unit = 101)

  open  (unit = 101, file = 'Outputs/Psi.txt')
  do i = 1, NRBOX
        write (101, '(1000e17.9)') (PSI (i, j), j = 1, NZBOX)
  enddo
  close (unit = 101)

  open  (unit = 101, file = 'Outputs/Profiles.txt')
  write (101, '(7e17.9)') (dble (i-1) /RMAX, T (i), P (i), TTp (i), Pp (i), Q (i), - R0*Pp (i) - TTp (i)/R0, i = 1, NRBOX)
  close (unit = 101)

  open  (unit = 101, file = 'Outputs/Points.txt')
  write (101, '(4i5)') NRBOX, NZBOX, NBOUND, NLIM
  close (unit = 101)

  open  (unit = 101, file = 'Outputs/Boundary.txt')
  write (101, '(2e17.9)') (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
  close (unit = 101)

  open  (unit = 101, file = 'Outputs/Limiter.txt')
  write (101, '(2e17.9)') (RLIM (i), ZLIM (i), i = 1, NLIM)
  close (unit = 101)

  ! ................
  ! Read target q_95
  ! ................
  
  open (unit = 100, file = 'Inputs/q95', status = 'old')

  read (100, *) q95_new

  close (100)

  print *, 'q95_tgt = ', q95_new

  ! ..................
  ! Determine old q_95
  ! ..................

  x95   = 0.95 * dble (NRBOX) + 1.
  i95   = int (x95)
  q95   = (dble (i95) - x95) * Q (i95+1) + (x95 + 1. - dble (i95)) * Q (i95)
  T95   = (dble (i95) - x95) * T (i95+1) + (x95 + 1. - dble (i95)) * T (i95)
  Tedge = T (NRBOX)

  print *, 'q95_old = ', q95, 'Psiaxis = ', PSIAXIS, 'Psibound', PSIBOUND, 'Current = ', CURRENT

  ! ...................
  ! Rescale equilibrium
  ! ...................

  a2 = (q95_new*q95_new /q95/q95 - 1.) * T95*T95

  do i = 1, NRBOX
     Told = T (i)
     if (Told > 0.) then
        Tnew =   sqrt (Told*Told + a2)
     else
        Tnew = - sqrt (Told*Told + a2)
     endif
     T (i) = Tnew

     qold  = Q (i)
     qnew  = qold * Tnew /Told
     Q (i) = qnew
  enddo

  a1 = Tedge /T (NRBOX)

  do i = 1, NRBOX
     Told  = T (i)
     Tnew  = a1 * Told
     T (i) = Tnew

     Told    = TTp (i)
     Tnew    = a1 * Told
     TTp (i) = Tnew

     Pold  = P (i)
     Pnew  = a1*a1 * Pold
     P (i) = Pnew

     Pold   = Pp (i)
     Pnew   = a1 * Pold
     Pp (i) = Pnew

     do j = 1, NZBOX
        Psiold     = Psi (i, j)
        Psinew     = a1 * Psiold
        Psi (i, j) = Psinew
     enddo
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

  print *, 'q95_new = ', q95, 'Psiaxis = ', PSIAXIS, 'Psibound', PSIBOUND, 'Current = ', CURRENT

  ! ............
  ! Find O-point
  ! ............

  dR = 1.e-2
  dZ = 1.e-2

  RO = Raxis
  ZO = Zaxis
  
  call FindOXPoint (RO, ZO, NRBOX, NZBOX, RR, ZZ, PSI, dR, dZ, PSIAXIS)

  !print *, 'O-point at R = ', RO, ' Z = ', ZO, 'PSIAXIS = ', PSIAXIS

  Raxis = RO
  Zaxis = ZO

  ! ...................................
  ! Find X-point (assume lower X-point)
  ! ...................................

  ZX = Zaxis
  do i = 1, NBOUND
     if (ZBOUND (i) < ZX) then
        ZX = ZBOUND (i)
        RX = RBOUND (i)
     endif
  enddo
  
  call FindOXPoint (RX, ZX, NRBOX, NZBOX, RR, ZZ, PSI, dR, dZ, PSIBOUND)

  !print *, 'X-point at R = ', RX, ' Z = ', ZX, 'PSIBOUND = ', PSIBOUND

  print *, 'q95_new = ', q95, 'Psiaxis = ', PSIAXIS, 'Psibound', PSIBOUND, 'Current = ', CURRENT
  print *, ''
  print *, 'Rescale parameters: a1 = ', a1, ' a1^2 = ', a1*a1, 'a2 = ', a2
  print *, ''
  
  ! ........................
  ! Output final equilibrium
  ! ........................

  open  (unit = 101, file = 'Outputs/Psiaxis1.txt')
  write (101, '(4e17.9)') PSIAXIS, PSIBOUND
  close (unit = 101)
  
  open  (unit = 101, file = 'Outputs/Psi1.txt')
  do i = 1, NRBOX
        write (101, '(1000e17.9)') (PSI (i, j), j = 1, NZBOX)
  enddo
  close (unit = 101)

  open  (unit = 101, file = 'Outputs/Profiles1.txt')
  write (101, '(7e17.9)') (dble (i-1) /RMAX, T (i), P (i), TTp (i), Pp (i), Q (i), - R0*Pp(i) - TTp(i)/R0, i = 1, NRBOX)
  close (unit = 101)

  ! ..................
  ! Write output gFile
  ! ..................

  call WritegFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
       B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)
   
  ! ........
  ! Clean up
  ! ........
  
  deallocate (T)
  deallocate (P)
  deallocate (TTp)
  deallocate (Pp)
  deallocate (Q)
  deallocate (PSI)
  deallocate (RBOUND)
  deallocate (ZBOUND)
  deallocate (RLIM)
  deallocate (ZLIM)
  deallocate (RR)
  deallocate (ZZ)
  
end program Rescale

! %%%%%%%%%%%%%%%%%%%%%%%%
! Subroutine to read gFile
! %%%%%%%%%%%%%%%%%%%%%%%%

subroutine ReadgFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
     B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)

  implicit none

  character (len = 100) :: string

  integer :: i, j, i3, NRBOX, NZBOX, NBOUND, NLIM
  
  double precision :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
 
  double precision, dimension (:), allocatable :: T, P, TTp, Pp, Q, RBOUND, ZBOUND, RLIM, ZLIM
  
  double precision, dimension (:, :), allocatable :: PSI

  open (unit = 100, file = 'Inputs/gFile', status = 'old')
  
  read (100, '(a48, 3i4)') string,  i3,      NRBOX,   NZBOX
  read (100, '(5e16.9  )') RBOXLEN, ZBOXLEN, R0,      RBOXLFT,  zero
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
  
  read (100, '(5e16.9)') ((PSI  (i, j), i = 1, NRBOX), j = 1, NZBOX)

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

subroutine WritegFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
     B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)

  implicit none

  character (len = 100) :: string

  integer :: i, j, i3, NRBOX, NZBOX, NBOUND, NLIM
  
  double precision :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
 
  double precision, dimension (:), allocatable :: T, P, TTp, Pp, Q, RBOUND, ZBOUND, RLIM, ZLIM
  
  double precision, dimension (:, :), allocatable :: PSI

  open (unit = 100, file = 'Outputs/gFile', status = 'unknown')

  zero = 0.
  
  write (100, '(a48, 3i4)') string,  i3,      NRBOX,    NZBOX
  write (100, '(5e16.9  )') RBOXLEN, ZBOXLEN, R0,       RBOXLFT,  zero
  write (100, '(5e16.9  )') RAXIS,   ZAXIS,   PSIAXIS,  PSIBOUND, B0
  write (100, '(5e16.9  )') CURRENT, PSIAXIS, zero,     RAXIS,    zero
  write (100, '(5e16.9  )') ZAXIS,   zero,    PSIBOUND, zero,     zero
  write (100, '(5e16.9)')  (T   (i), i = 1, NRBOX)
  write (100, '(5e16.9)')  (P   (i), i = 1, NRBOX)
  write (100, '(5e16.9)')  (TTp (i), i = 1, NRBOX)
  write (100, '(5e16.9)')  (Pp  (i), i = 1, NRBOX)
  write (100, '(5e16.9)')  ((PSI  (i, j), i = 1, NRBOX), j = 1, NZBOX)
  write (100, '(5e16.9)')  (Q (i), i = 1, NRBOX)
  write (100, '(2i5)')      NBOUND, NLIM
  write (100, '(5e16.9)')  (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
  write (100, '(5e16.9)')  (RLIM   (i), ZLIM   (i), i = 1, NLIM)

  close (unit = 100)
  
end subroutine WritegFile

double precision function GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, PSI)

  implicit none

  integer :: NRBOX, NZBOX, i, j

  double precision :: R, Z, x, y

  double precision, dimension (:), allocatable :: RR, ZZ
  
  double precision, dimension (:, :), allocatable :: PSI

  i = 1 + int ((R - RR (1)) /(RR (2) - RR (1)))
  j = 1 + int ((Z - ZZ (1)) /(ZZ (2) - ZZ (1)))

  x = (R - RR (i)) /(RR (2) - RR (1))
  y = (Z - ZZ (j)) /(ZZ (2) - ZZ (1))

  GetPsi = Psi (i, j) * (1.-x)*(1.-y) + Psi (i+1, j) * x*(1.-y) + Psi (i, j+1) * (1.-x)*y + Psi(i+1, j+1) * x*y

end function GetPsi

double precision function GetPsiR (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dR)

  interface
     double precision function GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, PSI)
       
       integer :: NRBOX, NZBOX

       double precision :: R, Z
       
       double precision, dimension (:), allocatable :: RR, ZZ
  
       double precision, dimension (:, :), allocatable :: PSI

     end function GetPsi
  end interface

  integer :: NRBOX, NZBOX

  double precision :: R, Z, dR, R1, R2 

  double precision, dimension (:), allocatable :: RR, ZZ
  
  double precision, dimension (:, :), allocatable :: PSI

  R1 = R - dR
  R2 = R + dR

  GetPsiR = (GetPsi (R2, Z, NRBOX, NZBOX, RR, ZZ, PSI) - GetPsi (R1, Z, NRBOX, NZBOX, RR, ZZ, PSI)) /2./dR

end function GetPsiR

double precision function GetPsiZ (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dZ)

  implicit none

  interface
     double precision function GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, PSI)

       integer :: NRBOX, NZBOX

       double precision :: R, Z

       double precision, dimension (:), allocatable :: RR, ZZ
  
       double precision, dimension (:, :), allocatable :: PSI

     end function GetPsi
  end interface

  integer :: NRBOX, NZBOX

  double precision :: R, Z, dZ, Z1, Z2 

  double precision, dimension (:), allocatable :: RR, ZZ
  
  double precision, dimension (:, :), allocatable :: PSI

  Z1 = Z - dZ
  Z2 = Z + dZ

  GetPsiZ = (GetPsi (R, Z2, NRBOX, NZBOX, RR, ZZ, PSI) - GetPsi (R, Z1, NRBOX, NZBOX, RR, ZZ, PSI)) /2./dZ

end function GetPsiZ

double precision function GetPsiRR (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dR)

  implicit none

  interface
     double precision function GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, PSI)

       integer :: NRBOX, NZBOX

       double precision :: R, Z, x, y

       double precision, dimension (:), allocatable :: RR, ZZ
  
       double precision, dimension (:, :), allocatable :: PSI

     end function GetPsi
  end interface

  integer :: NRBOX, NZBOX

  double precision :: R, Z, dR, R1, R2

  double precision, dimension (:), allocatable :: RR, ZZ
  
  double precision, dimension (:, :), allocatable :: PSI

  R1 = R - dR
  R2 = R + dR

  GetPsiRR = (GetPsi (R2, Z, NRBOX, NZBOX, RR, ZZ, PSI) - 2. * GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, PSI)&
       + GetPsi (R1, Z, NRBOX, NZBOX, RR, ZZ, PSI)) /dR/dR

end function GetPsiRR

double precision function GetPsiZZ (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dZ)

  implicit none

  interface
     double precision function GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, PSI)

       integer :: NRBOX, NZBOX

       double precision :: R, Z

       double precision, dimension (:), allocatable :: RR, ZZ
  
       double precision, dimension (:, :), allocatable :: PSI

     end function GetPsi
  end interface

  integer :: NRBOX, NZBOX

  double precision :: R, Z, dZ, Z1, Z2

  double precision, dimension (:), allocatable :: RR, ZZ
  
  double precision, dimension (:, :), allocatable :: PSI

  Z1 = Z - dZ
  Z2 = Z + dZ

  GetPsiZZ = (GetPsi (R, Z2, NRBOX, NZBOX, RR, ZZ, PSI) - 2. * GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, PSI)&
       +  GetPsi (R, Z1, NRBOX, NZBOX, RR, ZZ, PSI)) /dZ/dZ

end function GetPsiZZ

double precision function GetPsiRZ (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dR)

  implicit none

  interface
     double precision function GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, PSI)

       integer :: NRBOX, NZBOX

       double precision :: R, Z

       double precision, dimension (:), allocatable :: RR, ZZ
  
       double precision, dimension (:, :), allocatable :: PSI

     end function GetPsi
  end interface

  integer :: NRBOX, NZBOX

  double precision :: R, Z, dR, R1, R2, dZ, Z1, Z2

  double precision, dimension (:), allocatable :: RR, ZZ
  
  double precision, dimension (:, :), allocatable :: PSI

  R1 = R - dR
  R2 = R + dR
  Z1 = Z - dR
  Z2 = Z + dR

  GetPsiRZ = (GetPsi (R2, Z2, NRBOX, NZBOX, RR, ZZ, PSI) - GetPsi (R1, Z2, NRBOX, NZBOX, RR, ZZ, PSI)&
       -  GetPsi (R2, Z1, NRBOX, NZBOX, RR, ZZ, PSI) + GetPsi (R1, Z1, NRBOX, NZBOX, RR, ZZ, PSI)) /4./dR/dR

end function GetPsiRZ

subroutine FindOXPoint (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dR, dZ, p)
  
  implicit none

  interface
     double precision function GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, PSI)

       integer :: NRBOX, NZBOX
       
       double precision :: R, Z
       
       double precision, dimension (:), allocatable :: RR, ZZ
       
       double precision, dimension (:, :), allocatable :: PSI
       
     end function GetPsi
  end interface

  interface
     double precision function GetPsiR (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dR)

       integer :: NRBOX, NZBOX

       double precision :: R, Z, dR
    
       double precision, dimension (:), allocatable :: RR, ZZ
  
       double precision, dimension (:, :), allocatable :: PSI

     end function GetPsiR
  end interface

  interface
     double precision function GetPsiZ (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dZ)

       integer :: NRBOX, NZBOX
       
       double precision :: R, Z, dZ
       
       double precision, dimension (:), allocatable :: RR, ZZ
       
       double precision, dimension (:, :), allocatable :: PSI
       
     end function GetPsiZ
  end interface

  interface
     double precision function GetPsiRR (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dR)
       
       integer :: NRBOX, NZBOX
       
       double precision :: R, Z, dR
       
       double precision, dimension (:), allocatable :: RR, ZZ
       
       double precision, dimension (:, :), allocatable :: PSI
       
     end function GetPsiRR
  end interface

  interface
     double precision function GetPsiZZ (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dZ)

       integer :: NRBOX, NZBOX

       double precision :: R, Z, dZ

       double precision, dimension (:), allocatable :: RR, ZZ
       
       double precision, dimension (:, :), allocatable :: PSI
       
     end function GetPsiZZ
  end interface

  interface
     double precision function GetPsiRZ (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dR)

       integer :: NRBOX, NZBOX

       double precision :: R, Z, dR
       
       double precision, dimension (:), allocatable :: RR, ZZ
       
       double precision, dimension (:, :), allocatable :: PSI

     end function GetPsiRZ
  end interface

  integer :: NRBOX, NZBOX, i

  double precision :: R, Z, dR, dZ, p, pr, pz, prr, pzz, prz, det

  double precision, dimension (:), allocatable :: RR, ZZ
  
  double precision, dimension (:, :), allocatable :: PSI

  do i = 1, 50
     
     pr  = GetPsiR  (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dR)
     pz  = GetPsiZ  (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dZ)
     prr = GetPsiRR (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dR)
     pzz = GetPsiZZ (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dZ)
     prz = GetPsiRZ (R, Z, NRBOX, NZBOX, RR, ZZ, PSI, dR)

     det = prr*pzz - prz*prz

     R = R + (prz*pz - pzz*pr) /det
     Z = Z + (prz*pr - prr*pz) /det

     p = GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, PSI)

     !print *, 'i = ', i, 'R = ', R, 'Z = ', Z, 'Psi = ', p
     
  enddo

end subroutine FindOXPoint

