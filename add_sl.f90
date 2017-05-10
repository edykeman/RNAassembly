! ==============================================================================
! Subroutine: ADD_SL (VR,JOB,ISL)
! 
! Purpose: Adds a stem-loop to the list of present stem-loops using
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
!           ISL - Stem loop number to add to the list of present sl.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            05/01/2016   Original Code
!
! Dependancies:
!
! Modules - CLASS_VRNA
! Functions -
! Subroutines - GET_SL
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      SUBROUTINE ADD_SL (VR,JOB,ISL)

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

        IF ( job == 'P' ) vr% ipres(i) = isl

        IF ( MOD(i,2) == 0 ) i = i - 1

        vr% itree(indx,i) = vr% lstem(2,isl)


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

        CALL GET_SL (vr,job,isl,ic)

        IF ( ic == 0 ) THEN

          IF ( job == 'P' ) THEN
            vr% linkp(isl) = vr% ip
            vr% ip = isl
          ENDIF

          IF ( job == 'B' ) THEN
            vr% linkb(isl) = vr% ib
            vr% ib = isl
          ENDIF

        ELSE

          IF ( job == 'P' ) THEN
            vr% linkp(isl) = vr% linkp(ic)
            vr% linkp(ic)  = isl
          ENDIF

          IF ( job == 'B' ) THEN
            vr% linkb(isl) = vr% linkb(ic)
            vr% linkb(ic)  = isl
          ENDIF

        ENDIF

        RETURN

      END SUBROUTINE ADD_SL
