! ==============================================================================
! Subroutine: CAPSID_NUCL (VR,ISL)
! 
! Purpose: Performs four tasks regarding nucleating CP+RNA complexes
!          onto a capsid shell.
!
!          IF no SL number is given (ISL == 0):
!
!          (1) Computes reaction times for neigboring pairs of CP+RNA
!          complexes to nucleate into the initial nucleus of a capsid shell.
!          (2) If the initial nucleus is formed, it computes the initial
!          elongation reactions and turns off all other nucleation reactions.
!
!          IF a specific SL number is given (ISL /= 0):
!
!          (3) Computes the time to form the initial nucleus with its
!          CP bound SL neighbor.
!          (4) If the initial nucleus is formed, computes the elongation
!          reactions of the SL attaching to the capsid.
!
! Method: 
!
! Arguments:
!
!            VR - Class structure containing information on the
!                 RNA secondary structure and possible reactions.
!           ISL - Stem-loop number in which to compute times to nucleate
!                 onto the growing capsid shell.
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
! Functions - RANDOM, RKAPPA
! Subroutines - GET_SL, REQUEUE_SL, REQUEUE_CP
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      SUBROUTINE CAPSID_NUCL (VR,ISL)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(VIRAL_RNA), INTENT(INOUT) :: vr
        INTEGER, INTENT(IN) :: isl

        !=== VARIABLES ===!

        INTEGER :: i,j,i1,i2,i5,i3,nl
        INTEGER :: idis,icp,jcp
        INTEGER :: list(2,mxcp)

        DOUBLE PRECISION :: t,t1,t2,tau,tmin
        DOUBLE PRECISION :: r,x,dg,random,rkappa


        !=== Set SL Reactions To Unbinding ===!

        IF ( isl == 0 ) THEN

          i1 = vr% ib

          DO WHILE ( i1 /= 0 )

            IF ( vr% inxt_sl(i1) < 0 ) THEN
              vr% inxt_sl(i1) = 1
              vr% tnxt_sl(0,i1) = vr% tsave(i1)
            ENDIF

            IF ( vr% nt /= 0 ) THEN

              CALL REQUEUE_SL (vr,i1)

            ENDIF

            i1 = vr% linkb(i1)

          ENDDO

        ELSE

          IF ( vr% inxt_sl(isl) < 0 ) THEN

            vr% inxt_sl(isl) = 1
            vr% tnxt_sl(0,isl) = vr% tsave(isl)

            IF ( vr% nt /= 0 ) THEN

              CALL REQUEUE_SL (vr,isl)

            ENDIF

          ENDIF

        ENDIF


        !=== CASE 1: Initial Nucleus Not Formed ===!
        !=== Calculate Nucleation Times for ALL CP+SL ===!

        IF ( isl == 0 .and. vr% nt == 0 ) THEN 

          i1 = vr% ib

          DO WHILE ( i1 /= 0 )

            i2 = vr% linkb(i1)

            IF ( i2 /= 0 ) THEN

              idis = vr% lstem(1,i2) - vr% lstem(2,i1)

              r = RANDOM(iseed)

              tau = -DLOG(r) / rkappa(1,2,idis)

              t = time + tau

              t1 = vr% tnxt_sl(0,i1)
              t2 = vr% tnxt_sl(0,i2)

              IF ( t < t1 .and. t < t2 ) THEN

                vr% inxt_sl(i1) = -1

                vr% tsave(i1) = vr% tnxt_sl(0,i1)
                vr% tnxt_sl(0,i1) = t

              ENDIF

            ENDIF

            !=== Requeue ===!

            CALL REQUEUE_SL(vr,i1)

            i1 = i2

          ENDDO

        ENDIF


        !=== CASE 2: Initial Nucleus is Formed ===!
        !=== Calculate Initial Capsid Reactions ===!

        IF ( isl == 0 .and. vr% nt > 0 ) THEN

          i5 = vr% i5
          i3 = vr% i3

          list(:,:) = 0

          nl = 2

          list(1,1) = 1
          list(1,2) = 2

          list(2,1) = 1
          list(2,2) = 1

          DO i=1,nncp(1)

            icp = mapnn(i,1)

            IF ( list(2,icp) == 0 ) THEN
              nl = nl + 1
              list(1,nl) = icp
              list(2,icp) = 2
            ENDIF

          ENDDO

          DO i=1,nncp(2)

            icp = mapnn(i,2)

            IF ( list(2,icp) == 0 ) THEN
              nl = nl + 1
              list(1,nl) = icp
              list(2,icp) = 2
            ENDIF

          ENDDO

          !=== Calculate Reactions ===!

          DO i=1,nl

            jcp = list(1,i)

            !=== Calculate Removal Time for jCP ===!

            IF ( vr% info_cp(jcp) /= 0 ) THEN

              dg = 0.0d0

              DO j=1,nncp(jcp)

                icp = mapnn(j,jcp)

                IF ( vr% info_cp(icp) == 0 ) CYCLE

                dg = dg + gb(j,jcp)

              ENDDO

              !=== Compute Off-Rate ===!

              x = dg * beta
              x = DEXP(x)

              IF ( vr% info_cp(jcp) < 0 ) THEN
                x = rate_ab / x
              ELSEIF ( jcp > nps ) THEN
                x = rate_cc / x
              ELSE
                !x = rate_au / x
                x = rate_cc / x
              ENDIF

              r = RANDOM(iseed)

              tau = -DLOG(r) / x

              vr% tcoat(jcp) = tau

              vr% inxt_cp(jcp) = 1

              vr% tnxt_cp(0,jcp) = time + tau
              vr% tnxt_cp(1,jcp) = tinf

              CALL REQUEUE_CP (vr,jcp)

            ENDIF

            !=== Calculate Addition Time for jCP ===!

            IF ( vr% info_cp(jcp) == 0 ) THEN

              !=== Add Autosteric Reaction ===!

              IF ( jcp >  nps ) x = rkap_cc
              IF ( jcp <= nps ) x = rkap_au

              r = RANDOM(iseed)

              tau = -DLOG(r) / x

              vr% inxt_cp(jcp) = 2
              vr% tnxt_cp(0,jcp) = tinf
              vr% tnxt_cp(1,jcp) = tint + tau

              !=== Add 5'(3') Reaction ===!

              !=== 5' End ===!

              icp = vr% info_sl(i5)
              icp = IABS(icp)

              CALL GET_SL (vr,'B',i5,i1)

              IF ( i1 /= 0 ) THEN

                idis = vr% lstem(1,i5) - vr% lstem(2,i1)

                DO j=1,nnps
                IF ( maphp(j,jcp) == icp ) THEN

                  r = RANDOM(iseed)

                  tau = -DLOG(r) / rkappa(jcp,icp,idis)

                  vr% inxt_cp(jcp) = 5
                  vr% tnxt_cp(0,jcp) = time + tau
                  vr% tcoat(jcp) = tinf

                ENDIF
                ENDDO

              ENDIF

              !=== 3' End ===!

              icp = vr% info_sl(i3)
              icp = IABS(icp)

              i2 = vr% linkb(i3)

              IF ( i2 /= 0 ) THEN

                idis = vr% lstem(1,i2) - vr% lstem(2,i3)

                DO j=1,nnps
                IF ( maphp(j,jcp) == icp ) THEN

                  r = RANDOM(iseed)

                  tau = -DLOG(r) / rkappa(jcp,icp,idis)

                  t = time + tau

                  IF ( t < vr% tnxt_cp(0,jcp) ) THEN
                    vr% inxt_cp(jcp) = 3
                    vr% tcoat(jcp) = vr% tnxt_cp(0,jcp)
                    vr% tnxt_cp(0,jcp) = t
                  ELSE
                    vr% tcoat(jcp) = t
                  ENDIF

                ENDIF
                ENDDO

              ENDIF

              CALL REQUEUE_CP (vr,jcp)

            ENDIF

          ENDDO

        ENDIF


        !=== CASE 3: Initial Nucleus Not Formed ===!
        !=== Calculate Nucleation Reaction for ISL ===!

        IF ( isl /= 0 .and. vr% nt == 0 ) THEN

          i1 = isl
          i2 = vr% linkb(i1)

          IF ( i2 /= 0 ) THEN

            idis = vr% lstem(1,i2) - vr% lstem(2,i1)

            r = RANDOM(iseed)

            tau = -DLOG(r) / rkappa(1,2,idis)

            t = time + tau

            t1 = vr% tnxt_sl(0,i1)
            t2 = vr% tnxt_sl(0,i2)

            IF ( t < t1 .and. t < t2 ) THEN

              vr% inxt_sl(i1) = -1

              vr% tsave(i1) = vr% tnxt_sl(0,i1)
              vr% tnxt_sl(0,i1) = t

            ENDIF

          ENDIF

          !=== Requeue ===!

          CALL REQUEUE_SL(vr,i1)

        ENDIF


        !=== CASE 4: Initial Nucleus is Formed ===!
        !=== Calculate Nucleation Reaction for ISL ===!

        IF ( isl /= 0 .and. vr% nt > 0 ) THEN

          !=== 5' or 3' end? ===!

          i1 = 0
          i2 = 0
          nl = 0

          i5 = vr% i5
          i3 = vr% i3

          IF ( i5  == vr% linkb(isl) ) THEN

            nl = 5

            i1 = isl
            i2 = i5

            IF ( vr% info_sl(isl) == 1 ) THEN
              CALL GET_SL (vr,'B',isl,i1)
            ENDIF

            icp = vr% info_sl(i5)
            icp = IABS(icp)

          ENDIF

          IF ( isl == vr% linkb(i3)  ) THEN

            nl = 3

            i1 = i3
            i2 = isl

            IF ( vr% info_sl(isl) == 1 ) THEN
              i2 = vr% linkb(isl)
            ENDIF

            icp = vr% info_sl(i3)
            icp = IABS(icp)

          ENDIF

          IF ( i1 /= 0 .and. i2 /= 0 ) THEN

            idis = vr% lstem(1,i2) - vr% lstem(2,i1)

            !=== Calculate Time for all Neighbors ===!

            DO i=1,nnps

              jcp = maphp(i,icp)

              IF ( vr% info_cp(jcp) /= 0 ) CYCLE

              r = RANDOM(iseed)

              tau = -DLOG(r) / rkappa(icp,jcp,idis)

              t = time + tau

              IF ( vr% inxt_cp(jcp) == nl ) THEN

                IF ( t < vr% tcoat(jcp) ) THEN

                  vr% tnxt_cp(0,jcp) = t

                ELSE

                  IF ( nl == 5 ) vr% inxt_cp(jcp) = 3
                  IF ( nl == 3 ) vr% inxt_cp(jcp) = 5

                  vr% tnxt_cp(0,jcp) = vr %tcoat(jcp)

                ENDIF

              ELSE

                IF ( t < vr% tnxt_cp(0,jcp) ) THEN
                  vr% inxt_cp(jcp) = nl
                  vr% tcoat(jcp) = vr% tnxt_cp(0,jcp)
                  vr% tnxt_cp(0,jcp) = t
                ELSE
                  vr% tcoat(jcp) = t
                ENDIF

              ENDIF

              !=== Requeue ===!

              CALL REQUEUE_CP(vr,jcp)

            ENDDO

          ELSEIF ( nl /= 0 ) THEN

            DO i=1,nnps

              jcp = maphp(i,icp)

              IF ( vr% info_cp(jcp) /= 0 ) CYCLE

              IF ( vr% inxt_cp(jcp) == nl ) THEN

                IF ( nl == 5 ) vr% inxt_cp(jcp) = 3
                IF ( nl == 3 ) vr% inxt_cp(jcp) = 5

                vr% tnxt_cp(0,jcp) = vr %tcoat(jcp)

              ENDIF

              !=== Requeue ===!

              CALL REQUEUE_CP(vr,jcp)

            ENDDO

          ENDIF

        ENDIF

        RETURN

      END SUBROUTINE CAPSID_NUCL
