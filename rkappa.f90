! ==============================================================================
! Function: RKAPPA (ICP,JCP,IDIS)
! 
! Purpose: Calculates the rate of CP+SL binding to the growing capsid shell.
!          This depends on the RNA distance between the capsid shell and next
!          CP+SL.
!
! Method:
!
! Arguments:
!
!           ICP  - First Capsid Protein number in the shell
!           JCP  - Second Capsid Protein number in the shell.
!                  The CP+SL will attach at this position in the shell.
!           iDIS - Number of nucleotides between the capsid shell and CP+SL.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            05/01/2016   Original Code
!
! Dependancies:
!
! Modules - SystemVar
! Functions -
! Subroutines -
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      DOUBLE PRECISION FUNCTION RKAPPA (ICP,JCP,IDIS)

        USE SystemVar, ONLY : rate_ab,rkap_ab

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: icp,jcp,idis

        !=== VARIABLES ===!

        INTEGER :: i,j
        DOUBLE PRECISION :: x

        INTEGER, PARAMETER :: maxd = 20


        x = rkap_ab

        !=== Dodec Version ===!

        IF ( idis > maxd ) x = 1.0d-20

        !=== MS2 Version ===!
!
!        i = INT(icp/5) + 1
!        j = INT(jcp/5) + 1
!
!        IF ( i == j ) THEN
!
!          IF ( idis > maxd ) x = 0.0d0
!
!        ELSE
!
!          IF ( idis > maxd ) x = 0.0d0
!          IF ( idis < 30 ) x = 0.0d0
!
!        ENDIF

         RKAPPA = x

        RETURN

      END FUNCTION RKAPPA
