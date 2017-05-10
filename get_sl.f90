! ==============================================================================
! Subroutine: GET_SL (VR,JOB,ISL,IC)
! 
! Purpose: Given a stem-loop number (ISL), identifies the closest upstream
!          stem-loop present in the RNA strand using a binary tree.
!
! Method: 
!
! Arguments:
!
!            VR - Class structure containing information on the
!                 RNA secondary structure and possible reactions.
!           JOB - IF JOB = 'P', find closest present SL
!                 IF JOB = 'B', find closest CP bound SL
!           ISL - Stem loop number to find the closest sl upstream.
!            IC - Stem-loop number closest to ISL.
!                 IF: IC = 0 - There are no other SLs upstream present.
!                     IC = n - SL #n is the closest SL upstream to stem-loop ISL.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            05/01/2016   Original Code
!
! Dependancies:
!
! Modules -
! Functions -
! Subroutines -
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      SUBROUTINE GET_SL (VR,JOB,ISL,IC)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(VIRAL_RNA), INTENT(IN) :: vr

        CHARACTER, INTENT(IN) :: job

        INTEGER, INTENT(IN) :: isl
        INTEGER, INTENT(OUT) :: ic

        !=== VARIABLES ===!

        INTEGER :: i,j,k,n,n1,n2
        INTEGER :: i1,i2,i5,indx


        ic = 0

        IF ( job == 'P' ) indx = 1
        IF ( job == 'B' ) indx = 2

        i5 = vr% lstem(1,isl)

        i = i5
        IF ( MOD(i,2) == 0 ) i = i - 1


        !=== Find Position of Closest SL ===!

        n = 1
        n1= 2
        n2= 4

        DO WHILE ( ic == 0 .and. n1 < vr% nn )

          i = INT(i/n2) * n2 + n1

          j = i - n
          k = i + n

          i1 = vr% itree(indx,j)
          i2 = vr% itree(indx,k)

          IF ( i2 < i5 .and. i2 /= 0 ) THEN
            ic = vr% ipres(i2)
          ELSEIF ( i1 < i5 .and. i1 /= 0 ) THEN
            ic = vr% ipres(i1)
          ENDIF

          n = n1
          n1 = n2
          n2 = 2 * n2

        ENDDO

        RETURN

      END SUBROUTINE GET_SL
