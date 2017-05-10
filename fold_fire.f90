!==============================================================================
! Subroutine: FOLD_FIRE (VR)
! 
! Purpose: Fires an RNA stem-loop reaction for a viral RNA then updates
!          the reactions and the priority queue.
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
! Subroutines - GET_SL, REQUEUE_SL, FOLD_REAC, CAPSID_NUCL
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      SUBROUTINE FOLD_FIRE (VR)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(VIRAL_RNA), INTENT(INOUT) :: vr

        !=== VARIABLES ===!

        INTEGER :: i,j,k,n,n1,n2,nq
        INTEGER :: i1,i2,i5,i3,inxt
        INTEGER :: isl,jsl,icp

        DOUBLE PRECISION :: t,t1,t2,tmin,xp

        LOGICAL :: ok


        n = vr% ns / 2
        xp = DBLE(npro)


        !=== Find SL with Minimum Reaction Time ===!

        isl  = vr% ique_sl(0,n)
        tmin = vr% tnxt_sl(0,isl)
        inxt = vr% inxt_sl(isl)

        IF ( npro /= 0 ) THEN

          jsl = vr% ique_sl(1,n)

          t = time + ( vr% tnxt_sl(1,jsl) - tint ) / xp

          IF ( t < tmin ) THEN

            isl  = jsl
            tmin = t
            inxt = 2

          ENDIF

        ENDIF


        !=== CHECKS ===!

        IF ( inxt == 1 ) THEN

          !=== Check if SL Folding is Allowed ===!

          IF ( vr% info_sl(isl) == 0 ) THEN

            ok = .true.

            !=== Check for Overlap With Present SLs ===!

            IF ( vr% ip /= 0 ) THEN

              t = 0.0d0

              CALL GET_SL (vr,'P',isl,jsl)

              IF ( jsl == 0 ) jsl = vr% ip

              i1 = 0
              i5 = vr% lstem(1,isl)
              i3 = vr% lstem(2,isl)

              DO WHILE ( i1 <= i3 .and. jsl /= 0 )

                i1 = vr% lstem(1,jsl)
                i2 = vr% lstem(2,jsl)

                IF ( i1 <= i3 .and. i5 <= i2 ) THEN

                  !=== SL OVERLAPS ===!

                  ok = .false.

                  i = vr% info_sl(jsl)
                  j = vr% inxt_sl(jsl)

                  icp = IABS(i)

                  IF ( i < 0 ) THEN

                    t1 = vr% tnxt_cp(0,icp)
                    t1 = t1 + vr% tbind(jsl)
                    t1 = t1 + vr% tfold(jsl)

                  ELSE


                    IF ( j <  0 ) THEN
                      t1 = vr% tsave(jsl)
                    ELSE
                      t1 = vr% tnxt_sl(0,jsl)
                    ENDIF

                    IF ( i /= 1 ) t1 = t1 + vr% tfold(jsl)

                    !t1 = vr% tnxt_sl(0,jsl)

                    !IF ( j <  0 ) t1 = t1 + vr% tbind(jsl)
                    !IF ( i /= 1 ) t1 = t1 + vr% tfold(jsl)

                  ENDIF

                  t = MAX(t,t1)

                ENDIF

                jsl = vr% linkp(jsl)

              ENDDO

              t = t + vr% tfold(isl)

            ENDIF

            !=== Check RNA was SS for Correct Time ===!

            IF ( ok ) THEN

              t = 0.0d0

              i5 = vr% lstem(1,isl)
              i3 = vr% lstem(2,isl)

              DO i=i5,i3
                t1 = vr% tss(i)
                t = MAX(t,t1)
              ENDDO

              t = t + vr% tfold(isl)

              IF ( tmin < t ) ok = .false.

            ENDIF

            IF ( .not. ok ) THEN

              vr% tnxt_sl(0,isl) = t

              CALL REQUEUE_SL (vr,isl)

              RETURN

            ENDIF

          ENDIF

          !=== Check if Protein Unbinding is Allowed ===!

          IF ( vr% info_sl(isl) < 0 ) THEN

            i = vr% info_sl(isl)
            icp = IABS(i)

            t = vr% tnxt_cp(0,icp) + vr% tbind(isl)

            vr% tnxt_sl(0,isl) = t

            CALL REQUEUE_SL (vr,isl)

            RETURN

          ENDIF

          !=== Check CP+SL was Unbound From Capsid for Correct Time ===!

          IF ( vr% info_sl(isl) == 2 ) THEN

            t = vr% tub(isl) + vr% tbind(isl)

            IF ( tmin < t ) THEN

              vr% tnxt_sl(0,isl) = t

              CALL REQUEUE_SL (vr,isl)

              RETURN

            ENDIF

          ENDIF

        ELSEIF ( inxt == 2 ) THEN

          !=== Check if Protein Binding is Allowed ===!

          i5 = vr% i5
          i3 = vr% i3

          i1 = vr% lstem(1,isl)
          i2 = vr% lstem(2,isl)

          IF ( i5 /= 0 ) i5 = vr% lstem(1,i5)
          IF ( i3 /= 0 ) i3 = vr% lstem(2,i3)

          IF ( i5 <= i1 .and. i3 >= i2 ) THEN

            vr% tnxt_sl(1,isl) = tinf
            !vr% tnxt_sl(1,isl) = tint + vr% tbind(isl)

            CALL REQUEUE_SL (vr,isl)

            RETURN

          ENDIF

        ENDIF


        !=== Fire Reaction ===!

        !=== Update Time ===!

        tint = tint + xp * ( tmin - time )
 
        time = tmin

        !=== Compute New Reactions ===!

        IF ( inxt < 0 ) THEN

          i5 = isl
          i3 = vr% linkb(i5)

          vr% info_cp(1) = -i5
          vr% info_cp(2) = -i3
          vr% info_sl(i5)= -1
          vr% info_sl(i3)= -2

          vr% i5 = i5
          vr% i3 = i3
          vr% nt = 2

          CALL CAPSID_NUCL (vr,0)

        ELSE

          CALL FOLD_REAC (vr,isl,inxt)

        ENDIF

        !=== Update SS Info ===!

        IF ( inxt == 0 ) THEN

          i5 = vr% lstem(1,isl)
          i3 = vr% lstem(2,isl)

          DO i=i5,i3
          vr% tss(i) = tmin
          ENDDO

        ENDIF

        RETURN

      END SUBROUTINE FOLD_FIRE
