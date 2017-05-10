! ==============================================================================
! Subroutine: CAPSID_FIRE (VR)
! 
! Purpose: Fires a capsid reaction for assembly around a viral RNA then
!          updates the reactions and the priority queue.
!
! Method: 
!
! Arguments:
!
!            VR - Class structure containing information on the
!                 RNA secondary structure and possible reactions.
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
! Subroutines - GET_SL, REQUEUE_CP, CAPSID_FIRE, CAPSID_NUCL
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      SUBROUTINE CAPSID_FIRE (VR)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(VIRAL_RNA), INTENT(INOUT) :: vr

        !=== VARIABLES ===!

        INTEGER :: i,j,n,icp,jcp,nt
        INTEGER :: i5,i3,inxt,isl,iend

        DOUBLE PRECISION :: t,tmin,xp


        n = vr% nc / 2
        nt= vr% nt

        xp = DBLE(npro)


        !=== Find Minimum Reaction Time ===!

        icp  = vr% ique_cp(0,n)
        tmin = vr% tnxt_cp(0,icp)
        inxt = vr% inxt_cp(icp)

        IF ( npro /= 0 ) THEN

          jcp = vr% ique_cp(1,n)

          t = time + ( vr% tnxt_cp(1,jcp) - tint ) / xp

          IF ( t < tmin ) THEN

            icp  = jcp
            tmin = t
            inxt = 2

          ENDIF

        ENDIF


        !=== CHECKS ===!

        IF ( inxt == 0 ) THEN

          !=== Removal of CP is Not Allowed ===!

          IF ( vr% info_cp(icp) /= 0 ) THEN

            !=== Find Min Time for Next Removal ===!

            tmin = tinf

            DO i=1,ncp
            IF ( vr% inxt_cp(i) == 1 ) THEN
              t = vr% tnxt_cp(0,i)
              tmin = MIN(t,tmin)
            ENDIF
            ENDDO

            t = tmin + vr% tcoat(icp)

            vr% tnxt_cp(0,icp) = t
            vr% tnxt_cp(1,icp) = tinf

            CALL REQUEUE_CP (vr,icp)

            RETURN

          ENDIF

        ELSEIF ( inxt == 1 ) THEN

          !=== Check CP was Removable for Correct Time ===!

          t = vr% trm(icp) + vr% tcoat(icp)

          IF ( tmin < t ) THEN

            vr% tnxt_cp(0,icp) = t
            vr% tnxt_cp(1,icp) = tinf

            CALL REQUEUE_CP (vr,icp)

            RETURN

          ENDIF

        ENDIF


        !=== Fire Reaction ===!

        !=== Update Capsid State ===!

        IF ( inxt == 1 ) THEN

          !=== Remove CP ===!

          isl = 0
          iend = 2

          i5 = vr% i5
          i3 = vr% i3

          i5 = vr% info_sl(i5)
          i3 = vr% info_sl(i3)

          vr% inxt_cp(icp) = 0

          IF ( icp == IABS(i5) ) THEN

            isl = vr% i5

            vr% i5 = vr% linkb(isl)

            iend = 5

            !=== Reset CP Binding Times ===!

            i = vr% linkp(isl)

            DO WHILE ( i /= vr% i5 )

              vr% tnxt_sl(1,i) = tint + vr% tbind(i)

              CALL REQUEUE_SL (vr,i)

              i = vr% linkp(i)

            ENDDO

          ENDIF

          IF ( icp == IABS(i3) ) THEN

            isl = vr% i3

            CALL GET_SL (vr,'B',isl,i)

            vr% i3 = i

            iend = 3

            !=== Reset CP Binding Times ===!

            i = vr% linkp(i)

            DO WHILE ( i /= isl )

              vr% tnxt_sl(1,i) = tint + vr% tbind(i)

              CALL REQUEUE_SL (vr,i)

              i = vr% linkp(i)

            ENDDO

          ENDIF

          IF ( isl == 0 ) THEN

            !=== Autosteric (or C/C) Removal ===!

            npro = npro + 1

            vr% info_cp(icp) = 0
            vr% nt = nt - 1

          ELSE

            !=== CP+SL Removal ===!

            vr% info_cp(icp) = 0
            vr% info_sl(isl) = 2
            vr% nt = nt - 1

            vr% tub(isl) = tmin

            !=== Check For Un-nucleation ===!

            IF ( nt == 2 ) THEN

              isl = vr% i5
              icp = vr% info_sl(isl)

              icp = IABS(icp)

              vr% info_cp(icp) = 0
              vr% info_sl(isl) = 2

              vr% tub(isl) = tmin

              vr% i5 = 0
              vr% i3 = 0
              vr% nt = 0

            ENDIF

          ENDIF

        ELSEIF ( inxt == 2 ) THEN

          !=== Autosteric (or C/C) Addition ===!

          iend = 2

          npro = npro - 1

          vr% info_cp(icp) = 1
          vr% inxt_cp(icp) = 1

          vr% nt = nt + 1

        ELSEIF ( inxt >  2 ) THEN

          iend = inxt

          !=== CP+SL Addition ===!

          IF ( inxt == 5 ) THEN
            i5 = vr% i5
            CALL GET_SL (vr,'B',i5,isl)
            vr% i5 = isl
          ENDIF

          IF ( inxt == 3 ) THEN
            i3 = vr% i3
            isl = vr% linkb(i3)
            vr% i3 = isl
          ENDIF

          vr% info_cp(icp) = -isl
          vr% info_sl(isl) = -icp

          vr% inxt_cp(icp) = 1

          vr% nt = nt + 1

        ENDIF

        !=== Update Time ===!

        tint = tint + xp * ( tmin - time )
 
        time = tmin

        !=== Compute New Capsid Reactions ===!

        IF ( vr% nt == 0 ) THEN

          vr% inxt_cp(:) = 0
          vr% tnxt_cp(:,:) = tinf

          CALL CAPSID_NUCL (vr,0)

        ELSE

          CALL CAPSID_REAC (vr,icp,iend)

        ENDIF

        RETURN

      END SUBROUTINE CAPSID_FIRE
