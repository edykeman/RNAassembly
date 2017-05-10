! ==============================================================================
! Subroutine: TINT12 (LIST,EB)
! 
! Purpose: Performs a table lookup for the internal energy for an asymmetric
!          mismatch pair between two basepairs in a helix.
!
! Method: Uses the MFOLD 3.0 energy function table for RNA @ T=37.
!
! Arguments:
!
!          LIST - Array of length 4 containing the nucleotides in
!                 numerical code (A=1,C=2,G=3,U=4) for the
!                 following locations:
!
!                         (3)           
!                 5' (1) A .    X (6) 3'
!                 3' (2) U .  . Y (7) 5'
!                         (4)(5)             
!
!                 where LIST(1) = letter code for position 1 etc.
!
!            EB - (OUTPUT) MFOLD 3.0 internal loop energy of the sequence
!                 provided in LIST.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            10/01/2014   Original Code
!
! Dependencies:
!
! Modules -
! Functions -
! Subroutines -
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      SUBROUTINE TINT12 (LIST,EB)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: list(7)
        REAL, INTENT(INOUT) :: eb

        !=== VARIABLES ===!

        INTEGER :: i,j,k,i1,i2,i3,i4,i5,i6,i7
        REAL :: au(16,24),cg(16,24),gc(16,24)
        REAL :: ua(16,24),gu(16,24),ug(16,24)

        DATA (au(1,i),i=1,24) / 3.90, 3.70, 3.10, 5.50, 3.20, 3.00, &
            & 2.40, 4.80, 3.20, 3.00, 2.40, 4.80, 3.90, 3.70, 3.10, &
            & 5.50, 3.90, 3.70, 3.10, 5.50, 3.90, 3.70, 3.10, 5.50 /
        DATA (au(2,i),i=1,24) / 3.80, 3.70, 5.50, 3.70, 3.10, 3.00, &
            & 4.80, 3.00, 3.10, 3.00, 4.80, 3.00, 3.80, 3.70, 5.50, &
            & 3.70, 3.80, 3.70, 5.50, 3.70, 3.80, 3.70, 5.50, 3.70 /
        DATA (au(3,i),i=1,24) / 3.20, 5.50, 2.30, 5.50, 2.50, 4.80, &
            & 1.60, 4.80, 2.50, 4.80, 1.60, 4.80, 3.20, 5.50, 2.30, &
            & 5.50, 3.20, 5.50, 2.30, 5.50, 3.20, 5.50, 2.30, 5.50 /
        DATA (au(4,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
            & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
            & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (au(5,i),i=1,24) / 3.60, 3.20, 3.10, 5.50, 2.90, 2.50, &
            & 2.40, 4.80, 2.90, 2.50, 2.40, 4.80, 3.60, 3.20, 3.10, &
            & 5.50, 3.60, 3.20, 3.10, 5.50, 3.60, 3.20, 3.10, 5.50 /
        DATA (au(6,i),i=1,24) / 3.70, 4.00, 5.50, 3.70, 3.00, 3.30, &
            & 4.80, 3.00, 3.00, 3.30, 4.80, 3.00, 3.70, 4.00, 5.50, &
            & 3.70, 3.70, 4.00, 5.50, 3.70, 3.70, 4.00, 5.50, 3.70 /
        DATA (au(7,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
            & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
            & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (au(8,i),i=1,24) / 5.50, 3.70, 5.50, 2.80, 4.80, 3.00, &
            & 4.80, 2.10, 4.80, 3.00, 4.80, 2.10, 5.50, 3.70, 5.50, &
            & 2.80, 5.50, 3.70, 5.50, 2.80, 5.50, 3.70, 5.50, 2.80 /
        DATA (au(9,i),i=1,24) / 2.50, 2.10, 1.90, 5.50, 1.80, 1.40, &
            & 1.20, 4.80, 1.80, 1.40, 1.20, 4.80, 2.50, 2.10, 1.90, &
            & 5.50, 2.50, 2.10, 1.90, 5.50, 2.50, 2.10, 1.90, 5.50 /
        DATA (au(10,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (au(11,i),i=1,24) / 2.30, 5.50, 3.70, 5.50, 1.60, 4.80, &
             & 3.00, 4.80, 1.60, 4.80, 3.00, 4.80, 2.30, 5.50, 3.70, &
             & 5.50, 2.30, 5.50, 3.70, 5.50, 2.30, 5.50, 3.70, 5.50 /
        DATA (au(12,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (au(13,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (au(14,i),i=1,24) / 4.00, 3.40, 5.50, 3.70, 3.30, 2.70, &
             & 4.80, 3.00, 3.30, 2.70, 4.80, 3.00, 4.00, 3.40, 5.50, &
             & 3.70, 4.00, 3.40, 5.50, 3.70, 4.00, 3.40, 5.50, 3.70 /
        DATA (au(15,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (au(16,i),i=1,24) / 5.50, 3.20, 5.50, 2.70, 4.80, 2.50, &
             & 4.80, 2.00, 4.80, 2.50, 4.80, 2.00, 5.50, 3.20, 5.50, &
             & 2.70, 5.50, 3.20, 5.50, 2.70, 5.50, 3.20, 5.50, 2.70 /

        DATA (cg(1,i),i=1,24) / 3.20, 3.00, 2.40, 4.80, 2.30, 2.20, &
            & 1.10, 4.00, 2.40, 2.20, 1.60, 4.00, 3.20, 3.00, 2.40, &
            & 4.80, 3.20, 3.00, 2.40, 4.80, 3.20, 3.00, 2.40, 4.80 /
        DATA (cg(2,i),i=1,24) / 3.10, 3.00, 4.80, 3.00, 2.30, 2.20, &
            & 4.00, 2.20, 2.30, 2.20, 4.00, 2.20, 3.10, 3.00, 4.80, &
            & 3.00, 3.10, 3.00, 4.80, 3.00, 3.10, 3.00, 4.80, 3.00 /
        DATA (cg(3,i),i=1,24) / 2.50, 4.80, 1.60, 4.80, 1.70, 4.00, &
            & 0.80, 4.00, 1.70, 4.00, 0.80, 4.00, 2.50, 4.80, 1.60, &
            & 4.80, 2.50, 4.80, 1.60, 4.80, 2.50, 4.80, 1.60, 4.80 /
        DATA (cg(4,i),i=1,24) / 4.80, 4.80, 4.80, 4.80, 4.00, 4.00, &
            & 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.80, 4.80, 4.80, &
            & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80 /
        DATA (cg(5,i),i=1,24) / 2.90, 2.50, 2.40, 4.80, 2.10, 1.70, &
            & 1.60, 4.00, 2.10, 1.70, 1.60, 4.00, 2.90, 2.50, 2.40, &
            & 4.80, 2.90, 2.50, 2.40, 4.80, 2.90, 2.50, 2.40, 4.80 /
        DATA (cg(6,i),i=1,24) / 3.00, 3.30, 4.80, 3.00, 2.20, 2.50, &
            & 4.00, 2.20, 2.20, 2.50, 4.00, 2.20, 3.00, 3.30, 4.80, &
            & 3.00, 3.00, 3.30, 4.80, 3.00, 3.00, 3.30, 4.80, 3.00 /
        DATA (cg(7,i),i=1,24) / 4.80, 4.80, 4.80, 4.80, 4.00, 4.00, &
            & 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.80, 4.80, 4.80, &
            & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80 /
        DATA (cg(8,i),i=1,24) / 4.80, 3.00, 4.80, 2.10, 4.00, 2.20, &
            & 4.00, 1.50, 4.00, 2.20, 4.00, 1.30, 4.80, 3.00, 4.80, &
            & 2.10, 4.80, 3.00, 4.80, 2.10, 4.80, 3.00, 4.80, 2.10 /
        DATA (cg(9,i),i=1,24) / 1.80, 1.40, 1.20, 4.80, 0.80, 0.60, &
            & 0.40, 4.00, 1.00, 0.60, 0.40, 4.00, 1.80, 1.40, 1.20, &
            & 4.80, 1.80, 1.40, 1.20, 4.80, 1.80, 1.40, 1.20, 4.80 /
        DATA (cg(10,i),i=1,24) / 4.80, 4.80, 4.80, 4.80, 4.00, 4.00, &
             & 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.80, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80 /
        DATA (cg(11,i),i=1,24) / 1.60, 4.80, 3.00, 4.80, 0.80, 4.00, &
             & 2.20, 4.00, 0.80, 4.00, 2.20, 4.00, 1.60, 4.80, 3.00, &
             & 4.80, 1.60, 4.80, 3.00, 4.80, 1.60, 4.80, 3.00, 4.80 /
        DATA (cg(12,i),i=1,24) / 4.80, 4.80, 4.80, 4.80, 4.00, 4.00, &
             & 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.80, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80 /
        DATA (cg(13,i),i=1,24) / 4.80, 4.80, 4.80, 4.80, 4.00, 4.00, &
             & 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.80, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80 /
        DATA (cg(14,i),i=1,24) / 3.30, 2.70, 4.80, 3.00, 2.50, 1.90, &
             & 4.00, 2.20, 2.50, 1.90, 4.00, 2.20, 3.30, 2.70, 4.80, &
             & 3.00, 3.30, 2.70, 4.80, 3.00, 3.30, 2.70, 4.80, 3.00 /
        DATA (cg(15,i),i=1,24) / 4.80, 4.80, 4.80, 4.80, 4.00, 4.00, &
             & 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.80, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80 /
        DATA (cg(16,i),i=1,24) / 4.80, 2.50, 4.80, 2.00, 4.00, 1.70, &
             & 4.00, 1.20, 4.00, 1.70, 4.00, 1.20, 4.80, 2.50, 4.80, &
             & 2.00, 4.80, 2.50, 4.80, 2.00, 4.80, 2.50, 4.80, 2.00 /

        DATA (gc(1,i),i=1,24) / 3.20, 3.00, 2.40, 4.80, 2.40, 2.20, &
            & 1.60, 4.00, 2.50, 2.20, 2.10, 4.00, 3.20, 3.00, 2.40, &
            & 4.80, 3.20, 3.00, 2.40, 4.80, 3.20, 3.00, 2.40, 4.80 /
        DATA (gc(2,i),i=1,24) / 3.10, 3.00, 4.80, 3.00, 2.30, 2.20, &
            & 4.00, 2.20, 2.30, 2.20, 4.00, 2.20, 3.10, 3.00, 4.80, &
            & 3.00, 3.10, 3.00, 4.80, 3.00, 3.10, 3.00, 4.80, 3.00 /
        DATA (gc(3,i),i=1,24) / 2.50, 4.80, 1.60, 4.80, 1.70, 4.00, &
            & 0.80, 4.00, 1.70, 4.00, 0.80, 4.00, 2.50, 4.80, 1.60, &
            & 4.80, 2.50, 4.80, 1.60, 4.80, 2.50, 4.80, 1.60, 4.80 /
        DATA (gc(4,i),i=1,24) / 4.80, 4.80, 4.80, 4.80, 4.00, 4.00, &
            & 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.80, 4.80, 4.80, &
            & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80 /
        DATA (gc(5,i),i=1,24) / 2.90, 2.50, 2.40, 4.80, 2.10, 1.70, &
            & 1.60, 4.00, 2.10, 1.70, 1.60, 4.00, 2.90, 2.50, 2.40, &
            & 4.80, 2.90, 2.50, 2.40, 4.80, 2.90, 2.50, 2.40, 4.80 /
        DATA (gc(6,i),i=1,24) / 3.00, 3.30, 4.80, 3.00, 2.20, 2.50, &
            & 4.00, 2.20, 2.20, 2.50, 4.00, 2.20, 3.00, 3.30, 4.80, &
            & 3.00, 3.00, 3.30, 4.80, 3.00, 3.00, 3.30, 4.80, 3.00 /
        DATA (gc(7,i),i=1,24) / 4.80, 4.80, 4.80, 4.80, 4.00, 4.00, &
            & 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.80, 4.80, 4.80, &
            & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80 /
        DATA (gc(8,i),i=1,24) / 4.80, 3.00, 4.80, 2.10, 4.00, 2.20, &
            & 4.00, 1.30, 4.00, 2.20, 4.00, 1.20, 4.80, 3.00, 4.80, &
            & 2.10, 4.80, 3.00, 4.80, 2.10, 4.80, 3.00, 4.80, 2.10 /
        DATA (gc(9,i),i=1,24) / 1.80, 1.40, 1.20, 4.80, 1.00, 0.60, &
            & 0.40, 4.00, 1.20, 0.60, 0.40, 4.00, 1.80, 1.40, 1.20, &
            & 4.80, 1.80, 1.40, 1.20, 4.80, 1.80, 1.40, 1.20, 4.80 /
        DATA (gc(10,i),i=1,24) / 4.80, 4.80, 4.80, 4.80, 4.00, 4.00, &
             & 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.80, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80 /
        DATA (gc(11,i),i=1,24) / 1.60, 4.80, 3.00, 4.80, 0.80, 4.00, &
             & 2.20, 4.00, 0.80, 4.00, 2.20, 4.00, 1.60, 4.80, 3.00, &
             & 4.80, 1.60, 4.80, 3.00, 4.80, 1.60, 4.80, 3.00, 4.80 /
        DATA (gc(12,i),i=1,24) / 4.80, 4.80, 4.80, 4.80, 4.00, 4.00, &
             & 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.80, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80 /
        DATA (gc(13,i),i=1,24) / 4.80, 4.80, 4.80, 4.80, 4.00, 4.00, &
             & 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.80, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80 /
        DATA (gc(14,i),i=1,24) / 3.30, 2.70, 4.80, 3.00, 2.50, 1.90, &
             & 4.00, 2.20, 2.50, 1.90, 4.00, 2.20, 3.30, 2.70, 4.80, &
             & 3.00, 3.30, 2.70, 4.80, 3.00, 3.30, 2.70, 4.80, 3.00 /
        DATA (gc(15,i),i=1,24) / 4.80, 4.80, 4.80, 4.80, 4.00, 4.00, &
             & 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.80, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 4.80 /
        DATA (gc(16,i),i=1,24) / 4.80, 2.50, 4.80, 2.00, 4.00, 1.70, &
             & 4.00, 1.20, 4.00, 1.70, 4.00, 1.20, 4.80, 2.50, 4.80, &
             & 2.00, 4.80, 2.50, 4.80, 2.00, 4.80, 2.50, 4.80, 2.00 /

        DATA (ua(1,i),i=1,24) / 3.90, 3.70, 3.10, 5.50, 3.20, 3.00, &
            & 2.40, 4.80, 3.20, 3.00, 2.40, 4.80, 3.90, 3.70, 3.10, &
            & 5.50, 3.90, 3.70, 3.10, 5.50, 3.90, 3.70, 3.10, 5.50 /
        DATA (ua(2,i),i=1,24) / 3.80, 3.70, 5.50, 3.70, 3.10, 3.00, &
            & 4.80, 3.00, 3.10, 3.00, 4.80, 3.00, 3.80, 3.70, 5.50, &
            & 3.70, 3.80, 3.70, 5.50, 3.70, 3.80, 3.70, 5.50, 3.70 /
        DATA (ua(3,i),i=1,24) / 3.20, 5.50, 2.30, 5.50, 2.50, 4.80, &
            & 1.60, 4.80, 2.50, 4.80, 1.60, 4.80, 3.20, 5.50, 2.30, &
            & 5.50, 3.20, 5.50, 2.30, 5.50, 3.20, 5.50, 2.30, 5.50 /
        DATA (ua(4,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
            & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
            & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (ua(5,i),i=1,24) / 3.60, 3.20, 3.10, 5.50, 2.90, 2.50, &
            & 2.40, 4.80, 2.90, 2.50, 2.40, 4.80, 3.60, 3.20, 3.10, &
            & 5.50, 3.60, 3.20, 3.10, 5.50, 3.60, 3.20, 3.10, 5.50 /
        DATA (ua(6,i),i=1,24) / 3.70, 4.00, 5.50, 3.70, 3.00, 3.30, &
            & 4.80, 3.00, 3.00, 3.30, 4.80, 3.00, 3.70, 4.00, 5.50, &
            & 3.70, 3.70, 4.00, 5.50, 3.70, 3.70, 4.00, 5.50, 3.70 /
        DATA (ua(7,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
            & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
            & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (ua(8,i),i=1,24) / 5.50, 3.70, 5.50, 2.80, 4.80, 3.00, &
            & 4.80, 2.10, 4.80, 3.00, 4.80, 2.10, 5.50, 3.70, 5.50, &
            & 2.80, 5.50, 3.70, 5.50, 2.80, 5.50, 3.70, 5.50, 2.80 /
        DATA (ua(9,i),i=1,24) / 2.50, 2.10, 1.90, 5.50, 1.80, 1.40, &
            & 1.20, 4.80, 1.80, 1.40, 1.20, 4.80, 2.50, 2.10, 1.90, &
            & 5.50, 2.50, 2.10, 1.90, 5.50, 2.50, 2.10, 1.90, 5.50 /
        DATA (ua(10,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (ua(11,i),i=1,24) / 2.30, 5.50, 3.70, 5.50, 1.60, 4.80, &
             & 3.00, 4.80, 1.60, 4.80, 3.00, 4.80, 2.30, 5.50, 3.70, &
             & 5.50, 2.30, 5.50, 3.70, 5.50, 2.30, 5.50, 3.70, 5.50 /
        DATA (ua(12,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (ua(13,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (ua(14,i),i=1,24) / 4.00, 3.40, 5.50, 3.70, 3.30, 2.70, &
             & 4.80, 3.00, 3.30, 2.70, 4.80, 3.00, 4.00, 3.40, 5.50, &
             & 3.70, 4.00, 3.40, 5.50, 3.70, 4.00, 3.40, 5.50, 3.70 /
        DATA (ua(15,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (ua(16,i),i=1,24) / 5.50, 3.20, 5.50, 2.70, 4.80, 2.50, &
             & 4.80, 2.00, 4.80, 2.50, 4.80, 2.00, 5.50, 3.20, 5.50, &
             & 2.70, 5.50, 3.20, 5.50, 2.70, 5.50, 3.20, 5.50, 2.70 /

        DATA (gu(1,i),i=1,24) / 3.90, 3.70, 3.10, 5.50, 3.20, 3.00, &
            & 2.40, 4.80, 3.20, 3.00, 2.40, 4.80, 3.90, 3.70, 3.10, &
            & 5.50, 3.90, 3.70, 3.10, 5.50, 3.90, 3.70, 3.10, 5.50 /
        DATA (gu(2,i),i=1,24) / 3.80, 3.70, 5.50, 3.70, 3.10, 3.00, &
            & 4.80, 3.00, 3.10, 3.00, 4.80, 3.00, 3.80, 3.70, 5.50, &
            & 3.70, 3.80, 3.70, 5.50, 3.70, 3.80, 3.70, 5.50, 3.70 /
        DATA (gu(3,i),i=1,24) / 3.20, 5.50, 2.30, 5.50, 2.50, 4.80, &
            & 1.60, 4.80, 2.50, 4.80, 1.60, 4.80, 3.20, 5.50, 2.30, &
            & 5.50, 3.20, 5.50, 2.30, 5.50, 3.20, 5.50, 2.30, 5.50 /
        DATA (gu(4,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
            & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
            & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (gu(5,i),i=1,24) / 3.60, 3.20, 3.10, 5.50, 2.90, 2.50, &
            & 2.40, 4.80, 2.90, 2.50, 2.40, 4.80, 3.60, 3.20, 3.10, &
            & 5.50, 3.60, 3.20, 3.10, 5.50, 3.60, 3.20, 3.10, 5.50 /
        DATA (gu(6,i),i=1,24) / 3.70, 4.00, 5.50, 3.70, 3.00, 3.30, &
            & 4.80, 3.00, 3.00, 3.30, 4.80, 3.00, 3.70, 4.00, 5.50, &
            & 3.70, 3.70, 4.00, 5.50, 3.70, 3.70, 4.00, 5.50, 3.70 /
        DATA (gu(7,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
            & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
            & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (gu(8,i),i=1,24) / 5.50, 3.70, 5.50, 2.80, 4.80, 3.00, &
            & 4.80, 2.10, 4.80, 3.00, 4.80, 2.10, 5.50, 3.70, 5.50, &
            & 2.80, 5.50, 3.70, 5.50, 2.80, 5.50, 3.70, 5.50, 2.80 /
        DATA (gu(9,i),i=1,24) / 2.50, 2.10, 1.90, 5.50, 1.80, 1.40, &
            & 1.20, 4.80, 1.80, 1.40, 1.20, 4.80, 2.50, 2.10, 1.90, &
            & 5.50, 2.50, 2.10, 1.90, 5.50, 2.50, 2.10, 1.90, 5.50 /
        DATA (gu(10,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (gu(11,i),i=1,24) / 2.30, 5.50, 3.70, 5.50, 1.60, 4.80, &
             & 3.00, 4.80, 1.60, 4.80, 3.00, 4.80, 2.30, 5.50, 3.70, &
             & 5.50, 2.30, 5.50, 3.70, 5.50, 2.30, 5.50, 3.70, 5.50 /
        DATA (gu(12,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (gu(13,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (gu(14,i),i=1,24) / 4.00, 3.40, 5.50, 3.70, 3.30, 2.70, &
             & 4.80, 3.00, 3.30, 2.70, 4.80, 3.00, 4.00, 3.40, 5.50, &
             & 3.70, 4.00, 3.40, 5.50, 3.70, 4.00, 3.40, 5.50, 3.70 /
        DATA (gu(15,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (gu(16,i),i=1,24) / 5.50, 3.20, 5.50, 2.70, 4.80, 2.50, &
             & 4.80, 2.00, 4.80, 2.50, 4.80, 2.00, 5.50, 3.20, 5.50, &
             & 2.70, 5.50, 3.20, 5.50, 2.70, 5.50, 3.20, 5.50, 2.70 /

        DATA (ug(1,i),i=1,24) / 3.90, 3.70, 3.10, 5.50, 3.20, 3.00, &
            & 2.40, 4.80, 3.20, 3.00, 2.40, 4.80, 3.90, 3.70, 3.10, &
            & 5.50, 3.90, 3.70, 3.10, 5.50, 3.90, 3.70, 3.10, 5.50 /
        DATA (ug(2,i),i=1,24) / 3.80, 3.70, 5.50, 3.70, 3.10, 3.00, &
            & 4.80, 3.00, 3.10, 3.00, 4.80, 3.00, 3.80, 3.70, 5.50, &
            & 3.70, 3.80, 3.70, 5.50, 3.70, 3.80, 3.70, 5.50, 3.70 /
        DATA (ug(3,i),i=1,24) / 3.20, 5.50, 2.30, 5.50, 2.50, 4.80, &
            & 1.60, 4.80, 2.50, 4.80, 1.60, 4.80, 3.20, 5.50, 2.30, &
            & 5.50, 3.20, 5.50, 2.30, 5.50, 3.20, 5.50, 2.30, 5.50 /
        DATA (ug(4,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
            & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
            & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (ug(5,i),i=1,24) / 3.60, 3.20, 3.10, 5.50, 2.90, 2.50, &
            & 2.40, 4.80, 2.90, 2.50, 2.40, 4.80, 3.60, 3.20, 3.10, &
            & 5.50, 3.60, 3.20, 3.10, 5.50, 3.60, 3.20, 3.10, 5.50 /
        DATA (ug(6,i),i=1,24) / 3.70, 4.00, 5.50, 3.70, 3.00, 3.30, &
            & 4.80, 3.00, 3.00, 3.30, 4.80, 3.00, 3.70, 4.00, 5.50, &
            & 3.70, 3.70, 4.00, 5.50, 3.70, 3.70, 4.00, 5.50, 3.70 /
        DATA (ug(7,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
            & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
            & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (ug(8,i),i=1,24) / 5.50, 3.70, 5.50, 2.80, 4.80, 3.00, &
            & 4.80, 2.10, 4.80, 3.00, 4.80, 2.10, 5.50, 3.70, 5.50, &
            & 2.80, 5.50, 3.70, 5.50, 2.80, 5.50, 3.70, 5.50, 2.80 /
        DATA (ug(9,i),i=1,24) / 2.50, 2.10, 1.90, 5.50, 1.80, 1.40, &
            & 1.20, 4.80, 1.80, 1.40, 1.20, 4.80, 2.50, 2.10, 1.90, &
            & 5.50, 2.50, 2.10, 1.90, 5.50, 2.50, 2.10, 1.90, 5.50 /
        DATA (ug(10,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (ug(11,i),i=1,24) / 2.30, 5.50, 3.70, 5.50, 1.60, 4.80, &
             & 3.00, 4.80, 1.60, 4.80, 3.00, 4.80, 2.30, 5.50, 3.70, &
             & 5.50, 2.30, 5.50, 3.70, 5.50, 2.30, 5.50, 3.70, 5.50 /
        DATA (ug(12,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (ug(13,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (ug(14,i),i=1,24) / 4.00, 3.40, 5.50, 3.70, 3.30, 2.70, &
             & 4.80, 3.00, 3.30, 2.70, 4.80, 3.00, 4.00, 3.40, 5.50, &
             & 3.70, 4.00, 3.40, 5.50, 3.70, 4.00, 3.40, 5.50, 3.70 /
        DATA (ug(15,i),i=1,24) / 5.50, 5.50, 5.50, 5.50, 4.80, 4.80, &
             & 4.80, 4.80, 4.80, 4.80, 4.80, 4.80, 5.50, 5.50, 5.50, &
             & 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50, 5.50 /
        DATA (ug(16,i),i=1,24) / 5.50, 3.20, 5.50, 2.70, 4.80, 2.50, &
             & 4.80, 2.00, 4.80, 2.50, 4.80, 2.00, 5.50, 3.20, 5.50, &
             & 2.70, 5.50, 3.20, 5.50, 2.70, 5.50, 3.20, 5.50, 2.70 /


        !         (3)            !
        ! 5' (1) A .    X (6) 3' !
        ! 3' (2) U .  . Y (7) 5' !
        !         (4)(5)         !

        i1 = list(1)
        i2 = list(2)
        i3 = list(3)
        i4 = list(4)
        i5 = list(5)
        i6 = list(6)
        i7 = list(7)

        IF ( i5 == 1 ) j = 0 + i3
        IF ( i5 == 2 ) j = 4 + i3
        IF ( i5 == 3 ) j = 8 + i3
        IF ( i5 == 4 ) j =12 + i3

        IF ( i6 == 1 .and. i7 == 4 ) k = 0 + i4
        IF ( i6 == 2 .and. i7 == 3 ) k = 4 + i4
        IF ( i6 == 3 .and. i7 == 2 ) k = 8 + i4
        IF ( i6 == 4 .and. i7 == 1 ) k =12 + i4
        IF ( i6 == 3 .and. i7 == 4 ) k =16 + i4
        IF ( i6 == 4 .and. i7 == 3 ) k =20 + i4

        IF ( i1 == 1 .and. i2 == 4 ) eb = eb + au(j,k)
        IF ( i1 == 2 .and. i2 == 3 ) eb = eb + cg(j,k)
        IF ( i1 == 3 .and. i2 == 2 ) eb = eb + gc(j,k)
        IF ( i1 == 4 .and. i2 == 1 ) eb = eb + ua(j,k)
        IF ( i1 == 3 .and. i2 == 4 ) eb = eb + gu(j,k)
        IF ( i1 == 4 .and. i2 == 3 ) eb = eb + ug(j,k)

        RETURN

      END SUBROUTINE TINT12
