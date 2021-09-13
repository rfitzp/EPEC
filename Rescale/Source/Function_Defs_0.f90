Module Function_Defs_0

  interface
     subroutine NameListRead (q_95_new, OPOINT, XPOINT) 

       integer          :: OPOINT, XPOINT
       double precision :: q_95_new
  
 
     end subroutine namelistRead
  end interface   

  interface
     subroutine ReadgFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
          B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)
       
       character (len = 100) :: string
       
       integer :: i3, NRBOX, NZBOX, NBOUND, NLIM
       
       double precision :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
       
       double precision, dimension (:), allocatable :: T, P, TTp, Pp, Q, RBOUND, ZBOUND, RLIM, ZLIM
       
       double precision, dimension (:, :), allocatable :: PSI
       
     end subroutine ReadgFile
  end interface

  interface
     subroutine WritegFile (string, i3, NRBOX, NZBOX, RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, RAXIS, ZAXIS, PSIAXIS, PSIBOUND,&
          B0, CURRENT, T, P, TTp, Pp, Q, Psi, NBOUND, NLIM, RBOUND, ZBOUND, RLIM, ZLIM)

       character (len = 100) :: string
       
       integer :: i3, NRBOX, NZBOX, NBOUND, NLIM
       
       double precision :: RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZOFF, zero, RAXIS, ZAXIS, B0, PSIAXIS, PSIBOUND, CURRENT
       
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

End Module Function_Defs_0
