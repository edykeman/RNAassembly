! ==============================================================================
! Subroutine: SNREACTION (VRNA,TOUT)
! 
! Purpose: Performs a stochastic next reaction step
!
! Method: Anderson, D.F. "A modified next reaction method for simulating
!         chemical systems with time dependant propensities and delays"
!         J. Chem. Phys 127, 214107 (2007).
!
! Arguments:
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            05/01/2016   Original Code
!
! Dependancies:
!
! Modules - SYSTEMVAR, CLASS_VRNA
! Functions -
! Subroutines - FOLD_FIRE, CAPSID_FIRE
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      SUBROUTINE SNREACTION (VRNA,TOUT)

        USE SystemVar, ONLY : nrna,npro,nque,ique,time,tint
        USE Class_VRNA

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE (VIRAL_RNA), INTENT(INOUT) :: vrna(nrna)

        DOUBLE PRECISION, INTENT(IN) :: tout

        !=== VARIABLES ===!

        INTEGER :: i,j,k,n,irna,jrna
        INTEGER :: i1,i2,n1,n2,indx,ITOT(12)

        DOUBLE PRECISION :: t,t1,t2,tmin,xp


        xp = DBLE(npro)

        !=== Viral RNA Reaction Times ===!

        n = nque / 2

        irna = ique(0,n)
        tmin = vrna(irna)% t(0)
        indx = 0

        IF ( npro /= 0 ) THEN

          i = ique(1,n)

          t = time + ( vrna(i)% t(1) - tint ) / xp

          IF ( t < tmin ) THEN
            irna = i
            tmin = t
            indx = 1
          ENDIF

        ENDIF


        !=== Output Data ===!

        IF ( tmin > tout ) THEN

          itot(:) = 0
          do i=1,nrna
          j = vrna(i)% nt 
          if ( j /= 0 ) itot(j) = itot(j) + 1
          enddo

        ENDIF

        n1 = vrna(irna)% ns / 2
        n2 = vrna(irna)% nc / 2

        i1 = vrna(irna)% ique_sl(indx,n1)
        i2 = vrna(irna)% ique_cp(indx,n2)

        t1 = vrna(irna)% tnxt_sl(indx,i1)
        t2 = vrna(irna)% tnxt_cp(indx,i2)

        !=== Fire Reaction ===!

        IF ( t1 < t2 ) THEN

          CALL FOLD_FIRE (vrna(irna))

        ELSE

          CALL CAPSID_FIRE (vrna(irna))

        ENDIF

        !=== Output ===!

        IF ( time >= tmin .and. tmin > tout ) THEN

          write(*,*)tout,itot(1),itot(2),itot(3),itot(4),itot(5),&
         & itot(6),itot(7),itot(8),itot(9),itot(10),itot(11),itot(12)

        ENDIF

        !=== Find Min Times Between Queues ===!

        i1 = vrna(irna)% ique_sl(0,n1)
        i2 = vrna(irna)% ique_cp(0,n2)

        t1 = vrna(irna)% tnxt_sl(0,i1)
        t2 = vrna(irna)% tnxt_cp(0,i2)

        vrna(irna)% t(0) = MIN(t1,t2)

        i1 = vrna(irna)% ique_sl(1,n1)
        i2 = vrna(irna)% ique_cp(1,n2)

        t1 = vrna(irna)% tnxt_sl(1,i1)
        t2 = vrna(irna)% tnxt_cp(1,i2)

        vrna(irna)% t(1) = MIN(t1,t2)


        !=== Re-queue Table ===!

        n = 1
        n1= 2
        n2= 4

        i = irna
        IF ( MOD(i,2) == 0 ) i = i - 1

        !=== RNA Folding Queue ===!

        t1 = vrna(i+0)% t(0)
        t2 = vrna(i+1)% t(0)

        IF ( t1 <= t2 ) ique(0,i) = i
        IF ( t2 <  t1 ) ique(0,i) = i+1

        !=== RNA Binding Queue ===!

        t1 = vrna(i+0)% t(1)
        t2 = vrna(i+1)% t(1)

        IF ( t1 <= t2 ) ique(1,i) = i
        IF ( t2 <  t1 ) ique(1,i) = i+1


        DO WHILE ( n1 < nque )

          i = INT(i/n2) * n2 + n1

          j = i - n
          k = i + n

          !=== RNA Folding Queue ===!

          irna = ique(0,j)
          jrna = ique(0,k)

          t1 = vrna(irna)% t(0)
          t2 = vrna(jrna)% t(0)

          IF ( t1 <= t2 ) ique(0,i) = irna
          IF ( t2 <  t1 ) ique(0,i) = jrna

          !=== RNA Binding Queue ===!

          irna = ique(1,j)
          jrna = ique(1,k)

          t1 = vrna(irna)% t(1)
          t2 = vrna(jrna)% t(1)

          IF ( t1 <= t2 ) ique(1,i) = irna
          IF ( t2 <  t1 ) ique(1,i) = jrna

          n  = n1
          n1 = n2
          n2 = 2 * n2

        ENDDO

        RETURN

      END SUBROUTINE SNREACTION
