! ==============================================================================
! Subroutine: REQUEUE_SL (VR,ISL)
! 
! Purpose: RE-Queues the SL priority queue with new reactions times for
!          a given stem-loop.
!
! Method: 
!
! Arguments:
!
!            VR - Class structure containing information on the
!                 RNA secondary structure and possible reactions.
!           ISL - Stem loop number to re-queue reaction times for.
!
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

      SUBROUTINE REQUEUE_SL (VR,ISL)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(VIRAL_RNA), INTENT(INOUT) :: vr

        INTEGER, INTENT(IN) :: isl


        !=== VARIABLES ===!

        INTEGER :: i,j,k,n,n1,n2
        INTEGER :: i1,i2

        DOUBLE PRECISION :: t1,t2


        !=== Re-Queue ===!

        n = 1
        n1= 2
        n2= 4

        i = isl
        IF ( MOD(i,2) == 0 ) i = i - 1

        !=== RNA Folding Queue ===!

        t1 = vr% tnxt_sl(0,i+0)
        t2 = vr% tnxt_sl(0,i+1)

        IF ( t1 <= t2 ) vr% ique_sl(0,i) = i
        IF ( t2 <  t1 ) vr% ique_sl(0,i) = i+1

        !=== RNA Binding Queue ===!

        t1 = vr% tnxt_sl(1,i+0)
        t2 = vr% tnxt_sl(1,i+1)

        IF ( t1 <= t2 ) vr% ique_sl(1,i) = i
        IF ( t2 <  t1 ) vr% ique_sl(1,i) = i+1

        DO WHILE ( n1 < vr% ns )

          i = INT(i/n2) * n2 + n1

          j = i - n
          k = i + n

          !=== RNA Folding Queue ===!

          i1 = vr% ique_sl(0,j)
          i2 = vr% ique_sl(0,k)

          t1 = vr% tnxt_sl(0,i1)
          t2 = vr% tnxt_sl(0,i2)

          IF ( t1 <= t2 ) vr% ique_sl(0,i) = i1
          IF ( t2 <  t1 ) vr% ique_sl(0,i) = i2

          !=== RNA Binding Queue ===!

          i1 = vr% ique_sl(1,j)
          i2 = vr% ique_sl(1,k)

          t1 = vr% tnxt_sl(1,i1)
          t2 = vr% tnxt_sl(1,i2)

          IF ( t1 <= t2 ) vr% ique_sl(1,i) = i1
          IF ( t2 <  t1 ) vr% ique_sl(1,i) = i2

          n  = n1
          n1 = n2
          n2 = 2 * n2

        ENDDO

        RETURN

      END SUBROUTINE REQUEUE_SL
