Module Function_Defs_2

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
  
End Module Function_Defs_2
