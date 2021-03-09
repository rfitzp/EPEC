Module Function_Defs_1

  interface
     double precision function GetPsi (R, Z, NRBOX, NZBOX, RR, ZZ, PSI)

       integer :: NRBOX, NZBOX
       
       double precision :: R, Z
       
       double precision, dimension (:), allocatable :: RR, ZZ
       
       double precision, dimension (:, :), allocatable :: PSI
       
     end function GetPsi
  end interface
  
End Module Function_Defs_1
