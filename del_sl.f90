! ==============================================================================
! Subroutine: DEL_SL (VR,JOB,ISL)
! 
! Purpose: Deletes a stem-loop from the list of present stem-loops using
!          a binary tree and linked list.
!
! Method: 
!
! Arguments:
!
!            VR - Class structure containing information on the
!                 RNA secondary structure and possible reactions.
!           JOB - IF JOB = 'P', add present SL to list
!                 IF JOB = 'B', add CP bound SL to list
!           ISL - Stem loop number to delete from the list of present sl.
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

      SUBROUTINE DEL_SL (VR,JOB,ISL)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(VIRAL_RNA), INTENT(INOUT) :: vr

        CHARACTER, INTENT(IN) :: job

        INTEGER, INTENT(IN) :: isl

        !=== VARIABLES ===!

        INTEGER :: i,j,k,n,n1,n2
        INTEGER :: i1,i2,ic,indx


        i = vr% lstem(2,isl)

        IF ( job == 'P' ) indx = 1
        IF ( job == 'B' ) indx = 2

        IF ( job == 'P' ) vr% ipres(i) = 0

        IF ( MOD(i,2) == 0 ) i = i - 1

        vr% itree(indx,i) = 0


        !=== Re-Make Tree ===!

        n = 1
        n1= 2
        n2= 4

        DO WHILE ( n1 < vr% nn )

          i = INT(i/n2) * n2 + n1

          j = i - n
          k = i + n

          i1 = vr% itree(indx,j)
          i2 = vr% itree(indx,k)

          vr% itree(indx,i) = MAX(i1,i2)

          n  = n1
          n1 = n2
          n2 = 2 * n2

        ENDDO


        !=== Adjust Linked List ===!

        IF ( job == 'P' ) THEN

          IF ( vr% ip == isl ) THEN

            vr% ip = vr% linkp(isl)
            vr% linkp(isl) = 0

          ELSE

            CALL GET_SL (vr,job,isl,ic)

            IF ( ic /= 0 ) THEN
              vr% linkp(ic) = vr% linkp(isl)
              vr% linkp(isl) = 0
            ENDIF

          ENDIF

        ENDIF

        IF ( job == 'B' ) THEN

          IF ( vr% ib == isl ) THEN

            vr% ib = vr% linkb(isl)
            vr% linkb(isl) = 0

          ELSE

            CALL GET_SL (vr,job,isl,ic)

            IF ( ic /= 0 ) THEN
              vr% linkb(ic) = vr% linkb(isl)
              vr% linkb(isl) = 0
            ENDIF

          ENDIF

        ENDIF

        RETURN

      END SUBROUTINE DEL_SL
