! gFileInterpolate.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Functions to read multiple gFiles and write interpolated gFile
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine gFileInterpolateQuadratic () bind (c, name = 'gFileInterpolateQuadratic')

  implicit none

  character (len = 100) :: gFile1, gFile2, gFile
  double precision      :: time1,  time2,  weight1, weight2, time
  
  character (len = 100) :: string
  integer               :: i, j, i3
  double precision      :: zero

  integer                                         :: NRBOX1,   NZBOX1,   NBOUND1, NLIM1
  double precision                                :: RBOXLEN1, ZBOXLEN1, RBOXLFT1
  double precision                                :: RAXIS1,   ZAXIS1,   B01,     R01
  double precision, dimension (:),    allocatable :: T1,       P1,       TTp1,    Pp1,  Q1
  double precision, dimension (:),    allocatable :: RBOUND1,  ZBOUND1,  RLIM1,   ZLIM1  
  double precision, dimension (:, :), allocatable :: PSI1

  integer                                         :: NRBOX2,   NZBOX2,   NBOUND2, NLIM2
  double precision                                :: RBOXLEN2, ZBOXLEN2, RBOXLFT2
  double precision                                :: RAXIS2,   ZAXIS2,   B02,     R02
  double precision, dimension (:),    allocatable :: T2,       P2,       TTp2,    Pp2,  Q2
  double precision, dimension (:),    allocatable :: RBOUND2,  ZBOUND2,  RLIM2,   ZLIM2
  double precision, dimension (:, :), allocatable :: PSI2
  
  integer                                         :: NRBOX,   NZBOX,   NBOUND, NLIM
  double precision                                :: RBOXLEN, ZBOXLEN, RBOXLFT
  double precision                                :: RAXIS,   ZAXIS,   B0,     R0
  double precision, dimension (:),    allocatable :: T,       P,       TTp,    Pp,  Q
  double precision, dimension (:),    allocatable :: RBOUND,  ZBOUND,  RLIM,   ZLIM    
  double precision, dimension (:, :), allocatable :: PSI

  ! -----------------------
  ! Read file and time data
  ! -----------------------
  open (unit = 100, file = 'Interface.txt', status = 'old')

  read (100, '(a100)')  gFile1
  read (100, '(e16.9)') time1
  read (100, '(a100)')  gFile2
  read (100, '(e16.9)') time2
  read (100, '(a100)')  gFile
  read (100, '(e16.9)') time

  close (unit = 100)
  
  ! ----------------
  ! Read first gFile
  ! ----------------
  open (unit = 100, file = gFile1, status = 'old')
  
  read (100, '(a48, 3i4)') string,   i3,       NRBOX1, NZBOX1
  read (100, '(5e16.9  )') RBOXLEN1, ZBOXLEN1, R01,    RBOXLFT1,  zero
  read (100, '(5e16.9  )') RAXIS1,   ZAXIS1,   zero,   zero,      B01
  read (100, '(5e16.9  )') zero,     zero,     zero,   zero,      zero
  read (100, '(5e16.9  )') zero,     zero,     zero,   zero,      zero

  allocate (T1   (NRBOX1))
  allocate (P1   (NRBOX1))
  allocate (TTp1 (NRBOX1))
  allocate (Pp1  (NRBOX1))
  allocate (Q1   (NRBOX1))
  allocate (Psi1 (NRBOX1, NZBOX1))
  
  read (100, '(5e16.9)') (T1    (i),    i = 1, NRBOX1)
  read (100, '(5e16.9)') (P1    (i),    i = 1, NRBOX1)
  read (100, '(5e16.9)') (TTp1  (i),    i = 1, NRBOX1)
  read (100, '(5e16.9)') (Pp1   (i),    i = 1, NRBOX1)
  read (100, '(5e16.9)') ((PSI1 (i, j), i = 1, NRBOX1), j = 1, NZBOX1)
  read (100, '(5e16.9)') (Q1    (i),    i = 1, NRBOX1)

  read (100, '(2i5)') NBOUND1, NLIM1
  
  allocate (RBOUND1 (NBOUND1))
  allocate (ZBOUND1 (NBOUND1))
  allocate (RLIM1   (NLIM1))
  allocate (ZLIM1   (NLIM1))

  read (100, '(5e16.9)') (RBOUND1 (i), ZBOUND1 (i), i = 1, NBOUND1)
  read (100, '(5e16.9)') (RLIM1   (i), ZLIM1   (i), i = 1, NLIM1)

  close (unit = 100)

  ! -----------------
  ! Read second gFile
  ! -----------------
  open (unit = 100, file = gFile2, status = 'old')
  
  read (100, '(a48, 3i4)') string,   i3,       NRBOX2,  NZBOX2
  read (100, '(5e16.9  )') RBOXLEN2, ZBOXLEN2, R02,     RBOXLFT2, zero
  read (100, '(5e16.9  )') RAXIS2,   ZAXIS2,   zero,    zero,     B02
  read (100, '(5e16.9  )') zero,     zero,     zero,    zero,     zero
  read (100, '(5e16.9  )') zero,     zero,     zero,    zero,     zero

  allocate (T2   (NRBOX2))
  allocate (P2   (NRBOX2))
  allocate (TTp2 (NRBOX2))
  allocate (Pp2  (NRBOX2))
  allocate (Q2   (NRBOX2))
  allocate (Psi2 (NRBOX2, NZBOX2))
  
  read (100, '(5e16.9)') (T2    (i),    i = 1, NRBOX2)
  read (100, '(5e16.9)') (P2    (i),    i = 1, NRBOX2)
  read (100, '(5e16.9)') (TTp2  (i),    i = 1, NRBOX2)
  read (100, '(5e16.9)') (Pp2   (i),    i = 1, NRBOX2)
  read (100, '(5e16.9)') ((PSI2 (i, j), i = 1, NRBOX2), j = 1, NZBOX2)
  read (100, '(5e16.9)') (Q2    (i),    i = 1, NRBOX2)

  read (100, '(2i5)') NBOUND2, NLIM2
  
  allocate (RBOUND2 (NBOUND2))
  allocate (ZBOUND2 (NBOUND2))
  allocate (RLIM2   (NLIM2))
  allocate (ZLIM2   (NLIM2))

  read (100, '(5e16.9)') (RBOUND2 (i), ZBOUND2 (i), i = 1, NBOUND2)
  read (100, '(5e16.9)') (RLIM2   (i), ZLIM2   (i), i = 1, NLIM2)

  close (unit = 100)

  ! ......................
  ! Interpolate gFile data
  ! ......................
  if (NRBOX1 == NRBOX2) then
     NRBOX = NRBOX1
  else
     write (*,*) "Error - NRBOX mismatch:", NRBOX1, NRBOX2
     stop
  endif
  if (NZBOX1 == NZBOX2) then
     NZBOX = NZBOX1
  else
     write (*,*) "Error - NZBOX mismatch:", NZBOX1, NZBOX2
     stop
  endif
  if (NBOUND1 == NBOUND2) then
     NBOUND = NBOUND1
  else
     write (*,*) "Warning - NBOUND mismatch:", NBOUND1, NBOUND2

     if (NBOUND1 < NBOUND2) then
        NBOUND = NBOUND1
     else
        NBOUND = NBOUND2
     end if
  endif
  if (NLIM1 == NLIM2) then
     NLIM = NLIM1
  else
     write (*,*) "Warning - NLIM mismatch:", NLIM1, NLIM2

     if (NLIM1 < NLIM2) then
        NLIM = NLIM1
     else
        NLIM = NLIM2
     end if
  endif
  
  allocate (T      (NRBOX))
  allocate (P      (NRBOX))
  allocate (TTp    (NRBOX))
  allocate (Pp     (NRBOX))
  allocate (Q      (NRBOX))
  allocate (Psi    (NRBOX, NZBOX))
  allocate (RBOUND (NBOUND))
  allocate (ZBOUND (NBOUND))
  allocate (RLIM   (NLIM))
  allocate (ZLIM   (NLIM))

  weight1 = (time - time2) /(time1 - time2);
  weight2 = (time - time1) /(time2 - time1);
 
  R0      = R01      * weight1 + R02      * weight2
  B0      = B01      * weight1 + B02      * weight2
  RAXIS   = RAXIS1   * weight1 + RAXIS2   * weight2
  ZAXIS   = ZAXIS1   * weight1 + ZAXIS2   * weight2
  RBOXLEN = RBOXLEN1 * weight1 + RBOXLEN2 * weight2
  ZBOXLEN = ZBOXLEN1 * weight1 + ZBOXLEN2 * weight2
  RBOXLFT = RBOXLFT1 * weight1 + RBOXLFT2 * weight2
 
  do i = 1, NRBOX
     T   (i) = T1   (i) * weight1 + T2   (i) * weight2
     P   (i) = P1   (i) * weight1 + P2   (i) * weight2
     TTp (i) = TTp1 (i) * weight1 + TTp2 (i) * weight2
     Pp  (i) = Pp1  (i) * weight1 + Pp2  (i) * weight2
     Q   (i) = Q1   (i) * weight1 + Q2   (i) * weight2
  end do

  do i = 1, NRBOX
     do j = 1, NZBOX
        Psi (i, j) = Psi1 (i, j) * weight1 + Psi2 (i, j) * weight2 
     end do
  end do

  do i = 1, NBOUND
     RBOUND (i) = RBOUND1 (i) * weight1 + RBOUND2 (i) * weight2
     ZBOUND (i) = ZBOUND1 (i) * weight1 + ZBOUND2 (i) * weight2
  end do

  do i = 1, NLIM
     RLIM (i) = RLIM1 (i) * weight1 + RLIM2 (i) * weight2
     ZLIM (i) = ZLIM1 (i) * weight1 + ZLIM2 (i) * weight2
  end do

  ! .............................
  ! Write interpolated gFile data
  ! .............................
  open (unit = 100, file = gFile, status = 'replace')

  write (100, '(a48, 3i4)') string,  i3,      NRBOX, NZBOX
  write (100, '(5e16.9)'  ) RBOXLEN, ZBOXLEN, R0,    RBOXLFT, zero
  write (100, '(5e16.9)'  ) RAXIS,   ZAXIS,   zero,  zero,    B0
  write (100, '(5e16.9)'  ) zero,    zero,    zero,  RAXIS,   zero
  write (100, '(5e16.9)'  ) ZAXIS,   zero,    zero,  zero,    zero

  write (100, '(5e16.9)') (T    (i),    i = 1, NRBOX)
  write (100, '(5e16.9)') (P    (i),    i = 1, NRBOX)
  write (100, '(5e16.9)') (TTp  (i),    i = 1, NRBOX)
  write (100, '(5e16.9)') (Pp   (i),    i = 1, NRBOX)
  write (100, '(5e16.9)') ((PSI (i, j), i = 1, NRBOX), j = 1, NZBOX)
  write (100, '(5e16.9)') (Q    (i),    i = 1, NRBOX)

  write (100, '(2i5)'   ) NBOUND, NLIM
  write (100, '(5e16.9)') (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
  write (100, '(5e16.9)') (RLIM   (i), ZLIM   (i), i = 1, NLIM)

  close (unit = 100)

  write (*, *) "gFile Interpolation:"
  write (*, "(A, E11.4)") gFile1, weight1
  write (*, "(A, E11.4)") gFile2, weight2
  
  ! ........
  ! Clean up
  ! ........
  deallocate (T1)
  deallocate (P1)
  deallocate (TTp1)
  deallocate (Pp1)
  deallocate (Q1)
  deallocate (PSI1)
  deallocate (RBOUND1)
  deallocate (ZBOUND1)
  deallocate (RLIM1)
  deallocate (ZLIM1)
  
  deallocate (T2)
  deallocate (P2)
  deallocate (TTp2)
  deallocate (Pp2)
  deallocate (Q2)
  deallocate (PSI2)
  deallocate (RBOUND2)
  deallocate (ZBOUND2)
  deallocate (RLIM2)
  deallocate (ZLIM2)

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
  
endsubroutine gFileInterpolateQuadratic

subroutine gFileInterpolateCubic () bind (c, name = 'gFileInterpolateCubic')

  implicit none

  character (len = 100) :: gFile1, gFile2, gFile3, gFile
  double precision      :: time1,  time2,  time3,  weight1, weight2, weight3, time
  
  character (len = 100) :: string
  integer               :: i, j, i3
  double precision      :: zero

  integer                                         :: NRBOX1,   NZBOX1,   NBOUND1, NLIM1
  double precision                                :: RBOXLEN1, ZBOXLEN1, RBOXLFT1
  double precision                                :: RAXIS1,   ZAXIS1,   B01,     R01
  double precision, dimension (:),    allocatable :: T1,       P1,       TTp1,    Pp1,  Q1
  double precision, dimension (:),    allocatable :: RBOUND1,  ZBOUND1,  RLIM1,   ZLIM1  
  double precision, dimension (:, :), allocatable :: PSI1

  integer                                         :: NRBOX2,   NZBOX2,   NBOUND2, NLIM2
  double precision                                :: RBOXLEN2, ZBOXLEN2, RBOXLFT2
  double precision                                :: RAXIS2,   ZAXIS2,   B02,     R02
  double precision, dimension (:),    allocatable :: T2,       P2,       TTp2,    Pp2,  Q2
  double precision, dimension (:),    allocatable :: RBOUND2,  ZBOUND2,  RLIM2,   ZLIM2
  double precision, dimension (:, :), allocatable :: PSI2

  integer                                         :: NRBOX3,   NZBOX3,   NBOUND3, NLIM3
  double precision                                :: RBOXLEN3, ZBOXLEN3, RBOXLFT3
  double precision                                :: RAXIS3,   ZAXIS3,   B03,     R03
  double precision, dimension (:),    allocatable :: T3,       P3,       TTp3,    Pp3,  Q3
  double precision, dimension (:),    allocatable :: RBOUND3,  ZBOUND3,  RLIM3,   ZLIM3
  double precision, dimension (:, :), allocatable :: PSI3
  
  integer                                         :: NRBOX,   NZBOX,   NBOUND, NLIM
  double precision                                :: RBOXLEN, ZBOXLEN, RBOXLFT
  double precision                                :: RAXIS,   ZAXIS,   B0,     R0
  double precision, dimension (:),    allocatable :: T,       P,       TTp,    Pp,  Q
  double precision, dimension (:),    allocatable :: RBOUND,  ZBOUND,  RLIM,   ZLIM    
  double precision, dimension (:, :), allocatable :: PSI

  ! -----------------------
  ! Read file and time data
  ! -----------------------
  open (unit = 100, file = 'Interface.txt', status = 'old')

  read (100, '(a100)')  gFile1
  read (100, '(e16.9)') time1
  read (100, '(a100)')  gFile2
  read (100, '(e16.9)') time2
  read (100, '(a100)')  gFile3
  read (100, '(e16.9)') time3
  read (100, '(a100)')  gFile
  read (100, '(e16.9)') time

  close (unit = 100)
  
  ! ----------------
  ! Read first gFile
  ! ----------------
  open (unit = 100, file = gFile1, status = 'old')
  
  read (100, '(a48, 3i4)') string,   i3,       NRBOX1, NZBOX1
  read (100, '(5e16.9  )') RBOXLEN1, ZBOXLEN1, R01,    RBOXLFT1,  zero
  read (100, '(5e16.9  )') RAXIS1,   ZAXIS1,   zero,   zero,      B01
  read (100, '(5e16.9  )') zero,     zero,     zero,   zero,      zero
  read (100, '(5e16.9  )') zero,     zero,     zero,   zero,      zero

  allocate (T1   (NRBOX1))
  allocate (P1   (NRBOX1))
  allocate (TTp1 (NRBOX1))
  allocate (Pp1  (NRBOX1))
  allocate (Q1   (NRBOX1))
  allocate (Psi1 (NRBOX1, NZBOX1))
  
  read (100, '(5e16.9)') (T1    (i),    i = 1, NRBOX1)
  read (100, '(5e16.9)') (P1    (i),    i = 1, NRBOX1)
  read (100, '(5e16.9)') (TTp1  (i),    i = 1, NRBOX1)
  read (100, '(5e16.9)') (Pp1   (i),    i = 1, NRBOX1)
  read (100, '(5e16.9)') ((PSI1 (i, j), i = 1, NRBOX1), j = 1, NZBOX1)
  read (100, '(5e16.9)') (Q1    (i),    i = 1, NRBOX1)

  read (100, '(2i5)') NBOUND1, NLIM1
  
  allocate (RBOUND1 (NBOUND1))
  allocate (ZBOUND1 (NBOUND1))
  allocate (RLIM1   (NLIM1))
  allocate (ZLIM1   (NLIM1))

  read (100, '(5e16.9)') (RBOUND1 (i), ZBOUND1 (i), i = 1, NBOUND1)
  read (100, '(5e16.9)') (RLIM1   (i), ZLIM1   (i), i = 1, NLIM1)

  close (unit = 100)

  ! -----------------
  ! Read second gFile
  ! -----------------
  open (unit = 100, file = gFile2, status = 'old')
  
  read (100, '(a48, 3i4)') string,   i3,       NRBOX2,  NZBOX2
  read (100, '(5e16.9  )') RBOXLEN2, ZBOXLEN2, R02,     RBOXLFT2, zero
  read (100, '(5e16.9  )') RAXIS2,   ZAXIS2,   zero,    zero,     B02
  read (100, '(5e16.9  )') zero,     zero,     zero,    zero,     zero
  read (100, '(5e16.9  )') zero,     zero,     zero,    zero,     zero

  allocate (T2   (NRBOX2))
  allocate (P2   (NRBOX2))
  allocate (TTp2 (NRBOX2))
  allocate (Pp2  (NRBOX2))
  allocate (Q2   (NRBOX2))
  allocate (Psi2 (NRBOX2, NZBOX2))
  
  read (100, '(5e16.9)') (T2    (i),    i = 1, NRBOX2)
  read (100, '(5e16.9)') (P2    (i),    i = 1, NRBOX2)
  read (100, '(5e16.9)') (TTp2  (i),    i = 1, NRBOX2)
  read (100, '(5e16.9)') (Pp2   (i),    i = 1, NRBOX2)
  read (100, '(5e16.9)') ((PSI2 (i, j), i = 1, NRBOX2), j = 1, NZBOX2)
  read (100, '(5e16.9)') (Q2    (i),    i = 1, NRBOX2)

  read (100, '(2i5)') NBOUND2, NLIM2
  
  allocate (RBOUND2 (NBOUND2))
  allocate (ZBOUND2 (NBOUND2))
  allocate (RLIM2   (NLIM2))
  allocate (ZLIM2   (NLIM2))

  read (100, '(5e16.9)') (RBOUND2 (i), ZBOUND2 (i), i = 1, NBOUND2)
  read (100, '(5e16.9)') (RLIM2   (i), ZLIM2   (i), i = 1, NLIM2)

  close (unit = 100)

  ! ----------------
  ! Read third gFile
  ! ----------------
  open (unit = 100, file = gFile3, status = 'old')
  
  read (100, '(a48, 3i4)') string,   i3,       NRBOX3,  NZBOX3
  read (100, '(5e16.9  )') RBOXLEN3, ZBOXLEN3, R03,     RBOXLFT3, zero
  read (100, '(5e16.9  )') RAXIS3,   ZAXIS3,   zero,    zero,     B03
  read (100, '(5e16.9  )') zero,     zero,     zero,    zero,     zero
  read (100, '(5e16.9  )') zero,     zero,     zero,    zero,     zero

  allocate (T3   (NRBOX3))
  allocate (P3   (NRBOX3))
  allocate (TTp3 (NRBOX3))
  allocate (Pp3  (NRBOX3))
  allocate (Q3   (NRBOX3))
  allocate (Psi3 (NRBOX3, NZBOX3))
  
  read (100, '(5e16.9)') (T3    (i),    i = 1, NRBOX3)
  read (100, '(5e16.9)') (P3    (i),    i = 1, NRBOX3)
  read (100, '(5e16.9)') (TTp3  (i),    i = 1, NRBOX3)
  read (100, '(5e16.9)') (Pp3   (i),    i = 1, NRBOX3)
  read (100, '(5e16.9)') ((PSI3 (i, j), i = 1, NRBOX3), j = 1, NZBOX3)
  read (100, '(5e16.9)') (Q3    (i),    i = 1, NRBOX3)

  read (100, '(3i5)') NBOUND3, NLIM3
  
  allocate (RBOUND3 (NBOUND3))
  allocate (ZBOUND3 (NBOUND3))
  allocate (RLIM3   (NLIM3))
  allocate (ZLIM3   (NLIM3))

  read (100, '(5e16.9)') (RBOUND3 (i), ZBOUND3 (i), i = 1, NBOUND3)
  read (100, '(5e16.9)') (RLIM3   (i), ZLIM3   (i), i = 1, NLIM3)

  close (unit = 100)

  ! ......................
  ! Interpolate gFile data
  ! ......................
  if (NRBOX1 == NRBOX2 .and. NRBOX2 == NRBOX3) then
     NRBOX = NRBOX1
  else
     write (*,*) "Error - NRBOX mismatch:", NRBOX1, NRBOX2, NRBOX3
     stop
  endif
  if (NZBOX1 == NZBOX2 .and. NZBOX2 == NZBOX3) then
     NZBOX = NZBOX1
  else
     write (*,*) "Error - NZBOX mismatch:", NZBOX1, NZBOX2, NZBOX3
     stop
  endif
  if (NBOUND1 == NBOUND2 .and. NBOUND2 == NBOUND3) then
     NBOUND = NBOUND1
  else
     write (*,*) "Warning - NBOUND mismatch:", NBOUND1, NBOUND2, NBOUND3

     if (NBOUND1 < NBOUND2) then
        NBOUND = NBOUND1
     else
        NBOUND = NBOUND2
     end if

     if (NBOUND3 < NBOUND) then
        NBOUND = NBOUND3
     end if
  endif
  if (NLIM1 == NLIM2 .and. NLIM2 == NLIM3) then
     NLIM = NLIM1
  else
     write (*,*) "Warning - NLIM mismatch:", NLIM1, NLIM2, NLIM3

     if (NLIM1 < NLIM2) then
        NLIM = NLIM1
     else
        NLIM = NLIM2
     end if

     if (NLIM3 < NLIM) then
        NLIM = NLIM3
     end if
  endif
  
  allocate (T      (NRBOX))
  allocate (P      (NRBOX))
  allocate (TTp    (NRBOX))
  allocate (Pp     (NRBOX))
  allocate (Q      (NRBOX))
  allocate (Psi    (NRBOX, NZBOX))
  allocate (RBOUND (NBOUND))
  allocate (ZBOUND (NBOUND))
  allocate (RLIM   (NLIM))
  allocate (ZLIM   (NLIM))

  weight1 = (time - time2) * (time - time3) /(time1 - time2) /(time1 - time3)
  weight2 = (time - time1) * (time - time3) /(time2 - time1) /(time2 - time3)
  weight3 = (time - time1) * (time - time2) /(time3 - time1) /(time3 - time2)

  R0      = R01      * weight1 + R02      * weight2 + R03      * weight3
  B0      = B01      * weight1 + B02      * weight2 + B03      * weight3
  RAXIS   = RAXIS1   * weight1 + RAXIS2   * weight2 + RAXIS3   * weight3
  ZAXIS   = ZAXIS1   * weight1 + ZAXIS2   * weight2 + ZAXIS3   * weight3
  RBOXLEN = RBOXLEN1 * weight1 + RBOXLEN2 * weight2 + RBOXLEN3 * weight3
  ZBOXLEN = ZBOXLEN1 * weight1 + ZBOXLEN2 * weight2 + ZBOXLEN3 * weight3
  RBOXLFT = RBOXLFT1 * weight1 + RBOXLFT2 * weight2 + RBOXLFT3 * weight3
 
  do i = 1, NRBOX
     T   (i) = T1   (i) * weight1 + T2   (i) * weight2 + T3   (i) * weight3
     P   (i) = P1   (i) * weight1 + P2   (i) * weight2 + P3   (i) * weight3
     TTp (i) = TTp1 (i) * weight1 + TTp2 (i) * weight2 + TTp3 (i) * weight3
     Pp  (i) = Pp1  (i) * weight1 + Pp2  (i) * weight2 + Pp3  (i) * weight3
     Q   (i) = Q1   (i) * weight1 + Q2   (i) * weight2 + Q3   (i) * weight3
  end do

  do i = 1, NRBOX
     do j = 1, NZBOX
        Psi (i, j) = Psi1 (i, j) * weight1 + Psi2 (i, j) * weight2 + Psi3 (i, j) * weight3
     end do
  end do

  do i = 1, NBOUND
     RBOUND (i) = RBOUND1 (i) * weight1 + RBOUND2 (i) * weight2 + RBOUND3 (i) * weight3
     ZBOUND (i) = ZBOUND1 (i) * weight1 + ZBOUND2 (i) * weight2 + ZBOUND3 (i) * weight3
  end do

  do i = 1, NLIM
     RLIM (i) = RLIM1 (i) * weight1 + RLIM2 (i) * weight2 + RLIM3 (i) * weight3
     ZLIM (i) = ZLIM1 (i) * weight1 + ZLIM2 (i) * weight2 + ZLIM3 (i) * weight3
  end do

  ! .............................
  ! Write interpolated gFile data
  ! .............................
  open (unit = 100, file = gFile, status = 'replace')

  write (100, '(a48, 3i4)') string,  i3,      NRBOX, NZBOX
  write (100, '(5e16.9)'  ) RBOXLEN, ZBOXLEN, R0,    RBOXLFT, zero
  write (100, '(5e16.9)'  ) RAXIS,   ZAXIS,   zero,  zero,    B0
  write (100, '(5e16.9)'  ) zero,    zero,    zero,  RAXIS,   zero
  write (100, '(5e16.9)'  ) ZAXIS,   zero,    zero,  zero,    zero

  write (100, '(5e16.9)') (T    (i),    i = 1, NRBOX)
  write (100, '(5e16.9)') (P    (i),    i = 1, NRBOX)
  write (100, '(5e16.9)') (TTp  (i),    i = 1, NRBOX)
  write (100, '(5e16.9)') (Pp   (i),    i = 1, NRBOX)
  write (100, '(5e16.9)') ((PSI (i, j), i = 1, NRBOX), j = 1, NZBOX)
  write (100, '(5e16.9)') (Q    (i),    i = 1, NRBOX)

  write (100, '(2i5)'   ) NBOUND, NLIM
  write (100, '(5e16.9)') (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
  write (100, '(5e16.9)') (RLIM   (i), ZLIM   (i), i = 1, NLIM)

  close (unit = 100)

  write (*, *) "gFile Interpolation:"
  write (*, "(A, E11.4)") gFile1, weight1
  write (*, "(A, E11.4)") gFile2, weight2
  write (*, "(A, E11.4)") gFile3, weight3
  
  ! ........
  ! Clean up
  ! ........
  deallocate (T1)
  deallocate (P1)
  deallocate (TTp1)
  deallocate (Pp1)
  deallocate (Q1)
  deallocate (PSI1)
  deallocate (RBOUND1)
  deallocate (ZBOUND1)
  deallocate (RLIM1)
  deallocate (ZLIM1)
  
  deallocate (T2)
  deallocate (P2)
  deallocate (TTp2)
  deallocate (Pp2)
  deallocate (Q2)
  deallocate (PSI2)
  deallocate (RBOUND2)
  deallocate (ZBOUND2)
  deallocate (RLIM2)
  deallocate (ZLIM2)

  deallocate (T3)
  deallocate (P3)
  deallocate (TTp3)
  deallocate (Pp3)
  deallocate (Q3)
  deallocate (PSI3)
  deallocate (RBOUND3)
  deallocate (ZBOUND3)
  deallocate (RLIM3)
  deallocate (ZLIM3)

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
  
endsubroutine gFileInterpolateCubic

subroutine gFileInterpolateQuartic () bind (c, name = 'gFileInterpolateQuartic')

  implicit none

  character (len = 100) :: gFile1, gFile2, gFile3, gFile4, gFile
  double precision      :: time1,  time2,  time3,  time4, weight1, weight2, weight3, weight4, time
  
  character (len = 100) :: string
  integer               :: i, j, i3
  double precision      :: zero

  integer                                         :: NRBOX1,   NZBOX1,   NBOUND1, NLIM1
  double precision                                :: RBOXLEN1, ZBOXLEN1, RBOXLFT1
  double precision                                :: RAXIS1,   ZAXIS1,   B01,     R01
  double precision, dimension (:),    allocatable :: T1,       P1,       TTp1,    Pp1,  Q1
  double precision, dimension (:),    allocatable :: RBOUND1,  ZBOUND1,  RLIM1,   ZLIM1  
  double precision, dimension (:, :), allocatable :: PSI1

  integer                                         :: NRBOX2,   NZBOX2,   NBOUND2, NLIM2
  double precision                                :: RBOXLEN2, ZBOXLEN2, RBOXLFT2
  double precision                                :: RAXIS2,   ZAXIS2,   B02,     R02
  double precision, dimension (:),    allocatable :: T2,       P2,       TTp2,    Pp2,  Q2
  double precision, dimension (:),    allocatable :: RBOUND2,  ZBOUND2,  RLIM2,   ZLIM2
  double precision, dimension (:, :), allocatable :: PSI2

  integer                                         :: NRBOX3,   NZBOX3,   NBOUND3, NLIM3
  double precision                                :: RBOXLEN3, ZBOXLEN3, RBOXLFT3
  double precision                                :: RAXIS3,   ZAXIS3,   B03,     R03
  double precision, dimension (:),    allocatable :: T3,       P3,       TTp3,    Pp3,  Q3
  double precision, dimension (:),    allocatable :: RBOUND3,  ZBOUND3,  RLIM3,   ZLIM3
  double precision, dimension (:, :), allocatable :: PSI3

  integer                                         :: NRBOX4,   NZBOX4,   NBOUND4, NLIM4
  double precision                                :: RBOXLEN4, ZBOXLEN4, RBOXLFT4
  double precision                                :: RAXIS4,   ZAXIS4,   B04,     R04
  double precision, dimension (:),    allocatable :: T4,       P4,       TTp4,    Pp4,  Q4
  double precision, dimension (:),    allocatable :: RBOUND4,  ZBOUND4,  RLIM4,   ZLIM4
  double precision, dimension (:, :), allocatable :: PSI4
  
  integer                                         :: NRBOX,   NZBOX,   NBOUND, NLIM
  double precision                                :: RBOXLEN, ZBOXLEN, RBOXLFT
  double precision                                :: RAXIS,   ZAXIS,   B0,     R0
  double precision, dimension (:),    allocatable :: T,       P,       TTp,    Pp,  Q
  double precision, dimension (:),    allocatable :: RBOUND,  ZBOUND,  RLIM,   ZLIM    
  double precision, dimension (:, :), allocatable :: PSI

  ! -----------------------
  ! Read file and time data
  ! -----------------------
  open (unit = 100, file = 'Interface.txt', status = 'old')

  read (100, '(a100)')  gFile1
  read (100, '(e16.9)') time1
  read (100, '(a100)')  gFile2
  read (100, '(e16.9)') time2
  read (100, '(a100)')  gFile3
  read (100, '(e16.9)') time3
  read (100, '(a100)')  gFile4
  read (100, '(e16.9)') time4
  read (100, '(a100)')  gFile
  read (100, '(e16.9)') time

  close (unit = 100)
  
  ! ----------------
  ! Read first gFile
  ! ----------------
  open (unit = 100, file = gFile1, status = 'old')
  
  read (100, '(a48, 3i4)') string,   i3,       NRBOX1, NZBOX1
  read (100, '(5e16.9  )') RBOXLEN1, ZBOXLEN1, R01,    RBOXLFT1,  zero
  read (100, '(5e16.9  )') RAXIS1,   ZAXIS1,   zero,   zero,      B01
  read (100, '(5e16.9  )') zero,     zero,     zero,   zero,      zero
  read (100, '(5e16.9  )') zero,     zero,     zero,   zero,      zero

  allocate (T1   (NRBOX1))
  allocate (P1   (NRBOX1))
  allocate (TTp1 (NRBOX1))
  allocate (Pp1  (NRBOX1))
  allocate (Q1   (NRBOX1))
  allocate (Psi1 (NRBOX1, NZBOX1))
  
  read (100, '(5e16.9)') (T1    (i),    i = 1, NRBOX1)
  read (100, '(5e16.9)') (P1    (i),    i = 1, NRBOX1)
  read (100, '(5e16.9)') (TTp1  (i),    i = 1, NRBOX1)
  read (100, '(5e16.9)') (Pp1   (i),    i = 1, NRBOX1)
  read (100, '(5e16.9)') ((PSI1 (i, j), i = 1, NRBOX1), j = 1, NZBOX1)
  read (100, '(5e16.9)') (Q1    (i),    i = 1, NRBOX1)

  read (100, '(2i5)') NBOUND1, NLIM1
  
  allocate (RBOUND1 (NBOUND1))
  allocate (ZBOUND1 (NBOUND1))
  allocate (RLIM1   (NLIM1))
  allocate (ZLIM1   (NLIM1))

  read (100, '(5e16.9)') (RBOUND1 (i), ZBOUND1 (i), i = 1, NBOUND1)
  read (100, '(5e16.9)') (RLIM1   (i), ZLIM1   (i), i = 1, NLIM1)

  close (unit = 100)

  ! -----------------
  ! Read second gFile
  ! -----------------
  open (unit = 100, file = gFile2, status = 'old')
  
  read (100, '(a48, 3i4)') string,   i3,       NRBOX2,  NZBOX2
  read (100, '(5e16.9  )') RBOXLEN2, ZBOXLEN2, R02,     RBOXLFT2, zero
  read (100, '(5e16.9  )') RAXIS2,   ZAXIS2,   zero,    zero,     B02
  read (100, '(5e16.9  )') zero,     zero,     zero,    zero,     zero
  read (100, '(5e16.9  )') zero,     zero,     zero,    zero,     zero

  allocate (T2   (NRBOX2))
  allocate (P2   (NRBOX2))
  allocate (TTp2 (NRBOX2))
  allocate (Pp2  (NRBOX2))
  allocate (Q2   (NRBOX2))
  allocate (Psi2 (NRBOX2, NZBOX2))
  
  read (100, '(5e16.9)') (T2    (i),    i = 1, NRBOX2)
  read (100, '(5e16.9)') (P2    (i),    i = 1, NRBOX2)
  read (100, '(5e16.9)') (TTp2  (i),    i = 1, NRBOX2)
  read (100, '(5e16.9)') (Pp2   (i),    i = 1, NRBOX2)
  read (100, '(5e16.9)') ((PSI2 (i, j), i = 1, NRBOX2), j = 1, NZBOX2)
  read (100, '(5e16.9)') (Q2    (i),    i = 1, NRBOX2)

  read (100, '(2i5)') NBOUND2, NLIM2
  
  allocate (RBOUND2 (NBOUND2))
  allocate (ZBOUND2 (NBOUND2))
  allocate (RLIM2   (NLIM2))
  allocate (ZLIM2   (NLIM2))

  read (100, '(5e16.9)') (RBOUND2 (i), ZBOUND2 (i), i = 1, NBOUND2)
  read (100, '(5e16.9)') (RLIM2   (i), ZLIM2   (i), i = 1, NLIM2)

  close (unit = 100)

  ! ----------------
  ! Read third gFile
  ! ----------------
  open (unit = 100, file = gFile3, status = 'old')
  
  read (100, '(a48, 3i4)') string,   i3,       NRBOX3,  NZBOX3
  read (100, '(5e16.9  )') RBOXLEN3, ZBOXLEN3, R03,     RBOXLFT3, zero
  read (100, '(5e16.9  )') RAXIS3,   ZAXIS3,   zero,    zero,     B03
  read (100, '(5e16.9  )') zero,     zero,     zero,    zero,     zero
  read (100, '(5e16.9  )') zero,     zero,     zero,    zero,     zero

  allocate (T3   (NRBOX3))
  allocate (P3   (NRBOX3))
  allocate (TTp3 (NRBOX3))
  allocate (Pp3  (NRBOX3))
  allocate (Q3   (NRBOX3))
  allocate (Psi3 (NRBOX3, NZBOX3))
  
  read (100, '(5e16.9)') (T3    (i),    i = 1, NRBOX3)
  read (100, '(5e16.9)') (P3    (i),    i = 1, NRBOX3)
  read (100, '(5e16.9)') (TTp3  (i),    i = 1, NRBOX3)
  read (100, '(5e16.9)') (Pp3   (i),    i = 1, NRBOX3)
  read (100, '(5e16.9)') ((PSI3 (i, j), i = 1, NRBOX3), j = 1, NZBOX3)
  read (100, '(5e16.9)') (Q3    (i),    i = 1, NRBOX3)

  read (100, '(3i5)') NBOUND3, NLIM3
  
  allocate (RBOUND3 (NBOUND3))
  allocate (ZBOUND3 (NBOUND3))
  allocate (RLIM3   (NLIM3))
  allocate (ZLIM3   (NLIM3))

  read (100, '(5e16.9)') (RBOUND3 (i), ZBOUND3 (i), i = 1, NBOUND3)
  read (100, '(5e16.9)') (RLIM3   (i), ZLIM3   (i), i = 1, NLIM3)

  close (unit = 100)

  ! -----------------
  ! Read fourth gFile
  ! -----------------
  open (unit = 100, file = gFile4, status = 'old')
  
  read (100, '(a48, 3i4)') string,   i3,       NRBOX4,  NZBOX4
  read (100, '(5e16.9  )') RBOXLEN4, ZBOXLEN4, R04,     RBOXLFT4, zero
  read (100, '(5e16.9  )') RAXIS4,   ZAXIS4,   zero,    zero,     B04
  read (100, '(5e16.9  )') zero,     zero,     zero,    zero,     zero
  read (100, '(5e16.9  )') zero,     zero,     zero,    zero,     zero

  allocate (T4   (NRBOX4))
  allocate (P4   (NRBOX4))
  allocate (TTp4 (NRBOX4))
  allocate (Pp4  (NRBOX4))
  allocate (Q4   (NRBOX4))
  allocate (Psi4 (NRBOX4, NZBOX4))
  
  read (100, '(5e16.9)') (T4    (i),    i = 1, NRBOX4)
  read (100, '(5e16.9)') (P4    (i),    i = 1, NRBOX4)
  read (100, '(5e16.9)') (TTp4  (i),    i = 1, NRBOX4)
  read (100, '(5e16.9)') (Pp4   (i),    i = 1, NRBOX4)
  read (100, '(5e16.9)') ((PSI4 (i, j), i = 1, NRBOX4), j = 1, NZBOX4)
  read (100, '(5e16.9)') (Q4    (i),    i = 1, NRBOX4)

  read (100, '(3i5)') NBOUND4, NLIM4
  
  allocate (RBOUND4 (NBOUND4))
  allocate (ZBOUND4 (NBOUND4))
  allocate (RLIM4   (NLIM4))
  allocate (ZLIM4   (NLIM4))

  read (100, '(5e16.9)') (RBOUND4 (i), ZBOUND4 (i), i = 1, NBOUND4)
  read (100, '(5e16.9)') (RLIM4   (i), ZLIM4   (i), i = 1, NLIM3)

  close (unit = 100)

  ! ......................
  ! Interpolate gFile data
  ! ......................
  if (NRBOX1 == NRBOX2 .and. NRBOX2 == NRBOX3 .and. NRBOX3 == NRBOX4) then
     NRBOX = NRBOX1
  else
     write (*,*) "Error - NRBOX mismatch:", NRBOX1, NRBOX2, NRBOX3, NRBOX4
     stop
  endif
  if (NZBOX1 == NZBOX2 .and. NZBOX2 == NZBOX3 .and. NZBOX3 == NZBOX4) then
     NZBOX = NZBOX1
  else
     write (*,*) "Error - NZBOX mismatch:", NZBOX1, NZBOX2, NZBOX3, NZBOX4
     stop
  endif
  if (NBOUND1 == NBOUND2 .and. NBOUND2 == NBOUND3 .and. NBOUND3 == NBOUND4) then
     NBOUND = NBOUND1
  else
     write (*,*) "Warning - NBOUND mismatch:", NBOUND1, NBOUND2, NBOUND3, NBOUND4

     if (NBOUND1 < NBOUND2) then
        NBOUND = NBOUND1
     else
        NBOUND = NBOUND2
     end if

     if (NBOUND3 < NBOUND) then
        NBOUND = NBOUND3
     end if

     if (NBOUND4 < NBOUND) then
        NBOUND = NBOUND4
     end if
  endif
  if (NLIM1 == NLIM2 .and. NLIM2 == NLIM3 .and. NLIM3 == NLIM4) then
     NLIM = NLIM1
  else
     write (*,*) "Warning - NLIM mismatch:", NLIM1, NLIM2, NLIM3, NLIM4

     if (NLIM1 < NLIM2) then
        NLIM = NLIM1
     else
        NLIM = NLIM2
     end if

     if (NLIM3 < NLIM) then
        NLIM = NLIM3
     end if

     if (NLIM4 < NLIM) then
        NLIM = NLIM4
     end if
  endif
  
  allocate (T      (NRBOX))
  allocate (P      (NRBOX))
  allocate (TTp    (NRBOX))
  allocate (Pp     (NRBOX))
  allocate (Q      (NRBOX))
  allocate (Psi    (NRBOX, NZBOX))
  allocate (RBOUND (NBOUND))
  allocate (ZBOUND (NBOUND))
  allocate (RLIM   (NLIM))
  allocate (ZLIM   (NLIM))
  
  weight1 = (time - time2) * (time - time3) * (time - time4) /(time1 - time2) /(time1 - time3) /(time1 - time4)
  weight2 = (time - time1) * (time - time3) * (time - time4) /(time2 - time1) /(time2 - time3) /(time2 - time4)
  weight3 = (time - time1) * (time - time2) * (time - time4) /(time3 - time1) /(time3 - time2) /(time3 - time4)
  weight4 = (time - time1) * (time - time2) * (time - time3) /(time4 - time1) /(time4 - time2) /(time4 - time3)

  R0      = R01      * weight1 + R02      * weight2 + R03      * weight3 + R04      * weight4 
  B0      = B01      * weight1 + B02      * weight2 + B03      * weight3 + B04      * weight4
  RAXIS   = RAXIS1   * weight1 + RAXIS2   * weight2 + RAXIS3   * weight3 + RAXIS4   * weight4
  ZAXIS   = ZAXIS1   * weight1 + ZAXIS2   * weight2 + ZAXIS3   * weight3 + ZAXIS4   * weight4
  RBOXLEN = RBOXLEN1 * weight1 + RBOXLEN2 * weight2 + RBOXLEN3 * weight3 + RBOXLEN4 * weight4
  ZBOXLEN = ZBOXLEN1 * weight1 + ZBOXLEN2 * weight2 + ZBOXLEN3 * weight3 + ZBOXLEN4 * weight4
  RBOXLFT = RBOXLFT1 * weight1 + RBOXLFT2 * weight2 + RBOXLFT3 * weight3 + RBOXLFT4 * weight4
 
  do i = 1, NRBOX
     T   (i) = T1   (i) * weight1 + T2   (i) * weight2 + T3   (i) * weight3 + T4   (i) * weight4
     P   (i) = P1   (i) * weight1 + P2   (i) * weight2 + P3   (i) * weight3 + P4   (i) * weight4
     TTp (i) = TTp1 (i) * weight1 + TTp2 (i) * weight2 + TTp3 (i) * weight3 + TTp4 (i) * weight4
     Pp  (i) = Pp1  (i) * weight1 + Pp2  (i) * weight2 + Pp3  (i) * weight3 + Pp4  (i) * weight4
     Q   (i) = Q1   (i) * weight1 + Q2   (i) * weight2 + Q3   (i) * weight3 + Q4   (i) * weight4
  end do

  do i = 1, NRBOX
     do j = 1, NZBOX
        Psi (i, j) = Psi1 (i, j) * weight1 + Psi2 (i, j) * weight2 + Psi3 (i, j) * weight3 + Psi4 (i, j) * weight4
     end do
  end do

  do i = 1, NBOUND
     RBOUND (i) = RBOUND1 (i) * weight1 + RBOUND2 (i) * weight2 + RBOUND3 (i) * weight3 + RBOUND4 (i) * weight4
     ZBOUND (i) = ZBOUND1 (i) * weight1 + ZBOUND2 (i) * weight2 + ZBOUND3 (i) * weight3 + ZBOUND4 (i) * weight4
  end do

  do i = 1, NLIM
     RLIM (i) = RLIM1 (i) * weight1 + RLIM2 (i) * weight2 + RLIM3 (i) * weight3 + RLIM4 (i) * weight4
     ZLIM (i) = ZLIM1 (i) * weight1 + ZLIM2 (i) * weight2 + ZLIM3 (i) * weight3 + ZLIM4 (i) * weight4
  end do

  ! .............................
  ! Write interpolated gFile data
  ! .............................
  open (unit = 100, file = gFile, status = 'replace')

  write (100, '(a48, 3i4)') string,  i3,      NRBOX, NZBOX
  write (100, '(5e16.9)'  ) RBOXLEN, ZBOXLEN, R0,    RBOXLFT, zero
  write (100, '(5e16.9)'  ) RAXIS,   ZAXIS,   zero,  zero,    B0
  write (100, '(5e16.9)'  ) zero,    zero,    zero,  RAXIS,   zero
  write (100, '(5e16.9)'  ) ZAXIS,   zero,    zero,  zero,    zero

  write (100, '(5e16.9)') (T    (i),    i = 1, NRBOX)
  write (100, '(5e16.9)') (P    (i),    i = 1, NRBOX)
  write (100, '(5e16.9)') (TTp  (i),    i = 1, NRBOX)
  write (100, '(5e16.9)') (Pp   (i),    i = 1, NRBOX)
  write (100, '(5e16.9)') ((PSI (i, j), i = 1, NRBOX), j = 1, NZBOX)
  write (100, '(5e16.9)') (Q    (i),    i = 1, NRBOX)

  write (100, '(2i5)'   ) NBOUND, NLIM
  write (100, '(5e16.9)') (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
  write (100, '(5e16.9)') (RLIM   (i), ZLIM   (i), i = 1, NLIM)

  close (unit = 100)

  write (*, *) "gFile Interpolation:"
  write (*, "(A, E11.4)") gFile1, weight1
  write (*, "(A, E11.4)") gFile2, weight2
  write (*, "(A, E11.4)") gFile3, weight3
  write (*, "(A, E11.4)") gFile4, weight4
  
  ! ........
  ! Clean up
  ! ........
  deallocate (T1)
  deallocate (P1)
  deallocate (TTp1)
  deallocate (Pp1)
  deallocate (Q1)
  deallocate (PSI1)
  deallocate (RBOUND1)
  deallocate (ZBOUND1)
  deallocate (RLIM1)
  deallocate (ZLIM1)
  
  deallocate (T2)
  deallocate (P2)
  deallocate (TTp2)
  deallocate (Pp2)
  deallocate (Q2)
  deallocate (PSI2)
  deallocate (RBOUND2)
  deallocate (ZBOUND2)
  deallocate (RLIM2)
  deallocate (ZLIM2)

  deallocate (T3)
  deallocate (P3)
  deallocate (TTp3)
  deallocate (Pp3)
  deallocate (Q3)
  deallocate (PSI3)
  deallocate (RBOUND3)
  deallocate (ZBOUND3)
  deallocate (RLIM3)
  deallocate (ZLIM3)

  deallocate (T4)
  deallocate (P4)
  deallocate (TTp4)
  deallocate (Pp4)
  deallocate (Q4)
  deallocate (PSI4)
  deallocate (RBOUND4)
  deallocate (ZBOUND4)
  deallocate (RLIM4)
  deallocate (ZLIM4)

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
  
endsubroutine gFileInterpolateQuartic
