! ==============================================================================
! Subroutine: FOLD_REAC (VR,ISL,INXT)
! 
! Purpose: Calculates the possible folding reactions of a stem loop
!          and estimates a time in which this reaction will occur.
!
! Method: 
!
! Arguments:
!
!            VR - Class structure containing information on the
!                 RNA secondary structure and possible reactions.
!           ISL - Stem loop number to calculate reactions for.
!          INXT - The next state that the SL will transition to. 
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
! Functions - RANDOM
! Subroutines - GET_SL, ADD_SL, DEL_SL, CAPSID_NUCL, REQUEUE_SL
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      SUBROUTINE FOLD_REAC (VR,ISL,INXT)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(VIRAL_RNA), INTENT(INOUT) :: vr
        INTEGER, INTENT(IN) :: isl,inxt

        !=== VARIABLES ===!

        INTEGER :: icur,jsl
        DOUBLE PRECISION :: r,t,tau,random


        icur = vr% info_sl(isl)

        vr% info_sl(isl) = inxt

        IF ( icur == 0 .and. inxt == 1 ) THEN

          !=== Current State is SS ===!
          !=== Next State is Folded ==!

          !r = 0.904792018d0 * RANDOM(iseed) + 0.0000454d0
          r = RANDOM(iseed)

          tau = -DLOG(r) / vr% rfold(2,isl)

          vr% tfold(isl) = tau
          vr% tnxt_sl(0,isl) = time + tau
          vr% inxt_sl(isl) = 0

          r = RANDOM(iseed)

          tau = -DLOG(r) / vr% rbind(1,isl)

          vr% tbind(isl) = tau
          vr% tnxt_sl(1,isl) = tint + tau

          CALL ADD_SL (vr,'P',isl)

        ENDIF

        IF ( icur == 1 .and. inxt == 0 ) THEN

          !=== Current State is Folded ===!
          !=== NEXT State is SS ===!

          !r = 0.904792018d0 * RANDOM(iseed) + 0.0000454d0
          r = RANDOM(iseed)

          tau = -DLOG(r) / vr% rfold(1,isl)

          vr% tfold(isl)  = tau
          vr% tnxt_sl(0,isl)= time + tau
          vr% tnxt_sl(1,isl)= tinf

          vr% inxt_sl(isl) = 1

          CALL DEL_SL (vr,'P',isl)

        ENDIF

        IF ( icur == 1 .and. inxt == 2 ) THEN

          !=== Current State is Folded ===!
          !=== Next State is CP Bound ===!

          r = RANDOM(iseed)

          tau = -DLOG(r) / vr% rbind(2,isl)

          vr% tbind(isl) = tau

          vr% tnxt_sl(0,isl) = time + tau
          vr% tnxt_sl(1,isl) = tinf

          vr% inxt_sl(isl) = 1

          npro = npro - 1

          CALL ADD_SL (vr,'B',isl)

          !=== Compute Capsid Nucleation ===!

          CALL CAPSID_NUCL (vr,isl)

          IF ( vr% nt == 0 ) THEN

            CALL GET_SL (vr,'B',isl,jsl)

            IF ( jsl /= 0 ) THEN
            CALL CAPSID_NUCL (vr,jsl)
            ENDIF

          ENDIF

        ENDIF

        IF ( icur == 2 .and. inxt == 1 ) THEN

          !=== Current State is CP Bound ===!
          !=== Next State is Folded ===!

          tau = vr% tfold(isl)

          vr% tnxt_sl(0,isl) = time + tau
          vr% inxt_sl(isl) = 0

          r = RANDOM(iseed)

          tau = -DLOG(r) / vr% rbind(1,isl)

          vr% tbind(isl) = tau
          vr% tnxt_sl(1,isl) = tint + tau

          npro = npro + 1

          IF ( vr% nt > 0 ) THEN

            CALL CAPSID_NUCL (vr,isl)

          ENDIF

          CALL DEL_SL (vr,'B',isl)

          !=== Compute Capsid Nucleation ===!

          IF ( vr% nt == 0 ) THEN
 
            CALL GET_SL (vr,'B',isl,jsl)

            IF ( jsl /= 0 ) THEN
            CALL CAPSID_NUCL (vr,jsl)
            ENDIF

          ENDIF

        ENDIF


        !=== Requeue ISL ===!

        CALL REQUEUE_SL (vr,isl)

        RETURN

      END SUBROUTINE FOLD_REAC
