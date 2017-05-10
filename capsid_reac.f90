! ==============================================================================
! Subroutine: CAPSID_REAC (VR,ICP,IEND)
! 
! Purpose: Calculates and updates queue times for capsid assembly reactions.
!          The given capsid protein number (ICP) is the last addition/deletion
!          of CP made to the capsid shell.
!
! Method: 
!
! Arguments:
!
!            VR - Class structure containing information on the
!                 RNA secondary structure and possible reactions.
!           ICP - Capsid protein number which was last added/subtracted
!                 from the capsid shell.
!          IEND - Which end of the RNA (5' or 3') the CP+SL complex was
!                 added or subtracted from the capsid shell.
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
! Subroutines - GET_SL, REQUEUE_CP
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      SUBROUTINE CAPSID_REAC (VR,ICP,IEND)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(VIRAL_RNA), INTENT(INOUT) :: vr
        INTEGER, INTENT(IN) :: icp,iend

        !=== VARIABLES ===!

        INTEGER :: i,j,i5,i3,nt,nc,i0,i1,i2
        INTEGER :: is,isl,jcp,kcp,inew,iold
        INTEGER :: num(ncp),low(ncp),inum
        INTEGER :: ls(2*nbond),lt(2*nbond)
        INTEGER :: ireac(mxcp),idis

        DOUBLE PRECISION :: x,r,t,dg,tau
        DOUBLE PRECISION :: random,rkappa


        i5 = vr% i5
        i3 = vr% i3
        nt = vr% nt

        !=== Find CPs That Can be Added/Removed ===!

        is = 1
        nc = 0
        inum = 0

        ls(is) = 0
        lt(is) = 0
        num(:) = 0
        low(:) = 0

        ireac(:) = 0

        DO i=1,ncp
        IF ( vr% info_cp(i) /= 0 ) THEN

          ireac(i) = 1

          IF ( ls(1) == 0 ) ls(1) = i

        ENDIF
        ENDDO

        !=== Depth First Search ===!

        DO WHILE ( is /= 0 )

          i1 = ls(is)
          i0 = lt(is)

          IF ( num(i1) == 0 ) THEN

            inum = inum + 1
            num(i1) = inum
            low(i1) = inum

            IF ( i0 /= 0 ) THEN
            IF ( num(i0) == 1 ) nc = nc + 1
            ENDIF

            !=== Add to stack ===!

            DO i=1,nncp(i1)

              i2 = mapnn(i,i1)
              j = ireac(i2)

              IF ( j == 0 ) ireac(i2) = 2

              IF ( j == 1 ) THEN
              IF ( num(i2) == 0 ) THEN

                is = is + 1

                ls(is) = i2
                lt(is) = i1

              ELSEIF ( i2 /= i0 ) THEN

                IF ( num(i2) < num(i1) ) THEN
                low(i1) = MIN(low(i1),num(i2))
                ENDIF

              ENDIF
              ENDIF

            ENDDO

          ELSE

            !=== Back Track LOW ===!

            IF ( i0 /= 0 ) THEN

              low(i0) = MIN(low(i0),low(i1))

              IF ( num(i0) == 1 ) THEN

                IF ( nc > 1 ) ireac(i0) = 0

              ELSE

                IF ( low(i1) >= num(i0) ) ireac(i0) = 0

              ENDIF

              !=== Update CP Remove Reactions ===!

              IF ( ireac(i0) == 1 ) THEN

                IF ( vr% inxt_cp(i0) == 0 ) THEN
                  vr% trm(i0) = time
                  vr% inxt_cp(i0) = 1
                ENDIF

              ELSEIF ( ireac(i0) == 0 ) THEN

                vr% inxt_cp(i0) = 0

              ENDIF

            ENDIF

            !=== Clear Stack ===!

            is = is - 1

          ENDIF

        ENDDO

        !=== Insure That We Only Remove 5' 3' HPath Ends ===!

        IF ( vr% linkb(i5) == i3 .and. nt > 2 ) THEN

          jcp = vr% info_sl(i5)
          kcp = vr% info_sl(i3)

          jcp = IABS(jcp)
          kcp = IABS(kcp)

          vr% inxt_cp(jcp) = 0
          vr% inxt_cp(kcp) = 0

        ELSEIF ( vr% linkb(i5) /= i3 ) THEN

          isl = vr% linkb(i5)

          DO WHILE ( isl /= i3 )

            jcp = vr% info_sl(isl)
            jcp = IABS(jcp)

            vr% inxt_cp(jcp) = 0

            isl = vr% linkb(isl)

          ENDDO

        ENDIF


        !=== COMPUTE NEW CP REACTIONS ===!

        isl = 0
        jcp = 0

        !=== CASE 1: iCP was Just Added ===!

        IF ( vr% info_cp(icp) /= 0 ) THEN

          IF ( iend == 5 ) THEN

            isl = vr% linkb(i5)

            jcp = vr% info_sl(isl)
            jcp = IABS(jcp)

            CALL GET_SL (vr,'B',i5,isl)

            IF ( isl /= 0 ) THEN
              idis = vr% lstem(1,i5) - vr% lstem(2,isl)
            ENDIF

          ENDIF

          IF ( iend == 3 ) THEN

            CALL GET_SL (vr,'B',i3,isl)

            jcp = vr% info_sl(isl)
            jcp = IABS(jcp)

            isl = vr% linkb(i3)

            IF ( isl /= 0 ) THEN
              idis = vr% lstem(1,isl) - vr% lstem(2,i3)
            ENDIF

          ENDIF

          inew = icp
          iold = jcp

        ENDIF

        !=== CASE 2: iCP was Just Removed ===!

        IF ( vr% info_cp(icp) == 0 ) THEN

          IF ( iend == 5 ) THEN

            jcp = vr% info_sl(i5)
            jcp = IABS(jcp)

            CALL GET_SL (vr,'B',i5,isl)

            IF ( isl /= 0 ) THEN
              idis = vr% lstem(1,i5) - vr% lstem(2,isl)
            ENDIF

          ENDIF

          IF ( iend == 3 ) THEN

            jcp = vr% info_sl(i3)
            jcp = IABS(jcp)

            isl = vr% linkb(i3)

            IF ( isl /= 0 ) THEN
              idis = vr% lstem(1,isl) - vr% lstem(2,i3)
            ENDIF

          ENDIF

          inew = jcp
          iold = icp

        ENDIF


        !=== Reaction Adjustments for All Neighbors of iCP ===!

        DO i=0,nncp(icp)

          IF ( i == 0 ) THEN
            jcp = icp
          ELSE
            jcp = mapnn(i,icp)
          ENDIF

          !=== Recalculate Removal Time for Neighbor jCP ===!

          IF ( vr% info_cp(jcp) /= 0 ) THEN

            dg = 0.0d0

            DO j=1,nncp(jcp)

              kcp = mapnn(j,jcp)

              IF ( vr% info_cp(kcp) == 0 ) CYCLE

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
            vr% tnxt_cp(0,jcp) = time + tau
            vr% tnxt_cp(1,jcp) = tinf

            CALL REQUEUE_CP (vr,jcp)

          ENDIF

          !=== Recalculate Addition Time for Neighbor jCP ===!

          IF ( vr% info_cp(jcp) == 0 ) THEN

            !=== Remove Reaction from Queue ===!

            IF ( ireac(jcp) /= 2 ) THEN

              vr% inxt_cp(jcp) = 0
              vr% tnxt_cp(:,jcp) = tinf

            ELSEIF ( vr% inxt_cp(jcp) == 0 ) THEN

              !=== Add Autosteric Reaction ===!

              IF ( jcp >  nps ) x = rkap_cc
              IF ( jcp <= nps ) x = rkap_au

              r = RANDOM(iseed)

              tau = -DLOG(r) / x

              vr% inxt_cp(jcp) = 2
              vr% tnxt_cp(0,jcp) = tinf
              vr% tnxt_cp(1,jcp) = tint + tau

              !=== Add 5'(3') Reaction ===!

              IF ( iend == 2 ) THEN

                !=== 5' End ===!

                kcp = vr% info_sl(i5)
                kcp = IABS(kcp)

                CALL GET_SL (vr,'B',i5,isl)

                IF ( isl /= 0 ) THEN

                  idis = vr% lstem(1,i5) - vr% lstem(2,isl)

                  DO j=1,nnps
                  IF ( maphp(j,jcp) == kcp ) THEN

                    r = RANDOM(iseed)

                    tau = -DLOG(r) / rkappa(jcp,kcp,idis)

                    vr% inxt_cp(jcp) = 5
                    vr% tnxt_cp(0,jcp) = time + tau
                    vr% tcoat(jcp) = tinf

                  ENDIF
                  ENDDO

                ENDIF

                !=== 3' End ===!

                kcp = vr% info_sl(i3)
                kcp = IABS(kcp)

                isl = vr% linkb(i3)

                IF ( isl /= 0 ) THEN

                  idis = vr% lstem(1,isl) - vr% lstem(2,i3)

                  DO j=1,nnps
                  IF ( maphp(j,jcp) == kcp ) THEN

                    r = RANDOM(iseed)

                    tau = -DLOG(r) / rkappa(jcp,kcp,idis)

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

              ENDIF

            ENDIF

            CALL REQUEUE_CP (vr,jcp)

          ENDIF

        ENDDO


        !=== Reaction Adjustments for 5' (3') RNA Ends ===!

        IF ( iend /= 2 ) THEN

          !=== Loop Over Ham Paths of OLD 5'(3') End ===!

          DO i=1,nnps

            jcp = maphp(i,iold)

            IF ( vr% info_cp(jcp) /= 0 ) CYCLE

            !=== Remove 5'(3') Reaction ===!

            IF ( vr% inxt_cp(jcp) == iend ) THEN

              vr% tnxt_cp(0,jcp) = vr% tcoat(jcp)

              IF ( iend == 5 ) vr% inxt_cp(jcp) = 3
              IF ( iend == 3 ) vr% inxt_cp(jcp) = 5

              CALL REQUEUE_CP (vr,jcp)

            ENDIF

            vr% tcoat(jcp) = tinf

          ENDDO

          !=== Loop Over Ham Paths of NEW 5'(3') End ===!

          IF ( isl /= 0 ) THEN

            DO i=1,nnps

              jcp = maphp(i,inew)

              IF ( vr% info_cp(jcp) /= 0 ) CYCLE

              !=== Add 5'(3') Reaction ===!

              r = RANDOM(iseed)

              tau = -DLOG(r) / rkappa(inew,jcp,idis)

              t = time + tau

              IF ( vr% inxt_cp(jcp) /= iend ) THEN

                IF ( t < vr% tnxt_cp(0,jcp) ) THEN
                  vr% inxt_cp(jcp) = iend
                  vr% tcoat(jcp) = vr% tnxt_cp(0,jcp)
                  vr% tnxt_cp(0,jcp) = t
                ELSE
                  vr% tcoat(jcp) = t
                ENDIF

              ELSE

                vr% inxt_cp(jcp) = iend
                vr% tnxt_cp(0,jcp) = t

              ENDIF

              CALL REQUEUE_CP (vr,jcp)

            ENDDO

          ENDIF

        ENDIF

        RETURN

      END SUBROUTINE CAPSID_REAC
