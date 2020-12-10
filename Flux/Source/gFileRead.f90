! gFileRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read gFile and output data to Stage1 files

! Stage1/R0B0.txt ... R0(m), B0(T)

! All other quantities normalized to R0 and B0

! Stage1/Box.txt           ... Rmin, Zmin, Rmax, Zmax of bounding box
! Stage1/Axis.txt          ... R, Z coordinates of magnetic axis
! Stage1/R.txt             ... R grid
! Stage1/Z.txt             ... Z grid
! Stage1/Psi.txt           ... Psi(R, Z) in columns
! Stage1/PsiSequantial.txt ... Unnormalized Psi(R, Z) in single column
! Stage1/Profiles.txt      ... PsiN, T, P, TTP, PP, Q profiles
! Stage1/Points.txt        ... Numbers of points in R, Z, boundary and limiter arrays
! Stage1/Boundary.txt      ... R, Z coordinates of plasma boundary
! Stage1/Limiter.txt       ... R, Z coordinates of limiter
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine gFileRead () bind (c, name = 'gFileRead')

  implicit none
  
  character (len = 100) :: string
  integer               :: i,       j,        i3,      NRBOX,   NZBOX, NBOUND, NLIM
  double precision      :: RBOXLEN, ZBOXLEN,  R0,      RBOXLFT, zero,  RMAX
  double precision      :: RAXIS,   ZAXIS,    B0,      MU0,     ZMAX
  double precision      :: RLEFT,   RRIGHT,   ZLOW,    ZHIGH
  double precision      :: PSIAXIS, PSIBOUND, CURRENT

  double precision, dimension (:),    allocatable :: T,      P,      TTp,  Pp,  Q
  double precision, dimension (:),    allocatable :: RBOUND, ZBOUND, RLIM, ZLIM  
  double precision, dimension (:, :), allocatable :: PSI

  MU0 = 16. * atan(1.0) * 1.e-7

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
  
  open  (unit = 101, file = 'Outputs/Stage1/R0B0.txt')
  write (101, '(2e17.9)') R0, B0
  close (unit = 101)

  RLEFT  = RBOXLFT /R0
  RRIGHT = (RBOXLFT + RBOXLEN) /R0
  ZLOW   = -ZBOXLEN/2. /R0
  ZHIGH  =  ZBOXLEN/2. /R0

  open  (unit = 101, file = 'Outputs/Stage1/Box.txt')
  write (101, '(4e17.9)') RLEFT, ZLOW, RRIGHT, ZHIGH
  close (unit = 101)

  open  (unit = 101, file = 'Outputs/Stage1/Axis.txt')
  write (101, '(2e17.9)') RAXIS /R0, ZAXIS /R0
  close (unit = 101)
 
  RMAX = dble (NRBOX-1)

  open  (unit = 101, file = 'Outputs/Stage1/R.txt')
  do i = 1, NRBOX
     write (101, '(1e17.9)') RLEFT + (RRIGHT - RLEFT) * dble (i-1) /RMAX
  enddo
  close (unit = 101)

  ZMAX = dble (NZBOX-1)

  open  (unit = 101, file = 'Outputs/Stage1/Z.txt')
  do j = 1, NZBOX
     write (101, '(1e17.9)') ZLOW + (ZHIGH - ZLOW) * dble (j-1) /ZMAX
  enddo
  close (unit = 101)
  
  read (100, '(5e16.9)') (T   (i), i = 1, NRBOX)
  read (100, '(5e16.9)') (P   (i), i = 1, NRBOX)
  read (100, '(5e16.9)') (TTp (i), i = 1, NRBOX)
  read (100, '(5e16.9)') (Pp  (i), i = 1, NRBOX)
  
  read (100, '(5e16.9)') ((PSI  (i, j), i = 1, NRBOX), j = 1, NZBOX)

  open  (unit = 101, file = 'Outputs/Stage1/Psi.txt')
  do i = 1, NRBOX
        write (101, '(1000e17.9)') (PSI (i, j) /R0/R0/B0, j = 1, NZBOX)
  enddo
  close (unit = 101)

  open  (unit = 101, file = 'Outputs/Stage1/PsiSequential.txt')
  do i = 1, NRBOX
     do j = 1, NZBOX
         write (101, '(2i4,e17.9)') i, j, PSI (i, j)
     enddo
  enddo
  close (unit = 101)

  read (100, '(5e16.9)') (Q (i), i = 1, NRBOX)

  open  (unit = 101, file = 'Outputs/Stage1/RawProfiles.txt')
  write (101, '(7e17.9)') (dble (i-1) /RMAX, T (i), P (i), TTp (i), Pp (i), Q (i), - R0*Pp(i) - TTp(i)/R0, i = 1, NRBOX)
  close (unit = 101)

  open  (unit = 101, file = 'Outputs/Stage1/Profiles.txt')
  write (101, '(6e17.9)') (dble (i-1) /RMAX, T (i) /R0/B0, P (i) *MU0/B0/B0, TTp (i) /B0, Pp (i) *MU0*R0*R0/B0, Q (i), i = 1, NRBOX)
  close (unit = 101)

  read (100, '(2i5)') NBOUND, NLIM

  open  (unit = 101, file = 'Outputs/Stage1/Points.txt')
  write (101, '(4i5)') NRBOX, NZBOX, NBOUND, NLIM
  close (unit = 101)
  
  allocate (RBOUND (NBOUND))
  allocate (ZBOUND (NBOUND))
  allocate (RLIM   (NLIM))
  allocate (ZLIM   (NLIM))

  read (100, '(5e16.9)') (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
  read (100, '(5e16.9)') (RLIM   (i), ZLIM   (i), i = 1, NLIM)

  open  (unit = 101, file = 'Outputs/Stage1/Boundary.txt')
  write (101, '(2e17.9)') (RBOUND (i)/R0, ZBOUND (i)/R0, i = 1, NBOUND)
  close (unit = 101)

  open  (unit = 101, file = 'Outputs/Stage1/Limiter.txt')
  write (101, '(2e17.9)') (RLIM (i)/R0, ZLIM (i)/R0, i = 1, NLIM)
  close (unit = 101)

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
  
  close (unit = 100)

endsubroutine gFileRead
