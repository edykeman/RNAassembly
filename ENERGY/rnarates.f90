! ==============================================================================
! Subroutine: RNARATES (ISEQ,IBSP,N,I5,RF,RU,RON,ROFF,LMIN)
! 
! Purpose: Computes the energy of an RNA Stem-Loop using the emperical
!          MFOLD 3.0 energy function along with CP binding rates.
!
! Method: Uses the MFOLD 3.0 energy function for RNA @ T=37 given by:
!            
!            E = SUM E_loop + E_stack
!
! Arguments:
!
!          ISEQ - Array of length N containing the sequence
!                 in numerical code (A=1,C=2,G=3,U=4)
!          IBSP - Array of dimension (N) containing the information
!                 on base pairs in the RNA fold.
!                 IBSP(i) = j [i base pairs with j]
!                 IBSP(i) = 0 [i is single stranded]
!             N - Number of nucleotides in the sequence.
!            I5 - Location of the 5' nucleotide of the closing BP.
!            RF - Rate of RNA folding for the sequence.
!            RU - Rate of RNA unfolding for the sequence.
!           RON - Coat Protein binding rate to the Stem-loop.
!          ROFF - Coat protein unbinding rate.
!          LMIN - Is stem-loop a local minimum in energy.
!                 TRUE = Yes, FALSE = No
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            04/01/2016   Original Code
!
! Dependencies:
!
! Modules - 
! Functions -
! Subroutines - ESTACK, EHAIR, EBULGE
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      SUBROUTINE RNARATES (ISEQ,IBSP,N,I5,RF,RU,RON,ROFF,LMIN)

        USE SystemVar, ONLY : ratef,beta,vol,mxbp

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: n,i5
        INTEGER, INTENT(IN) :: iseq(n),ibsp(n)

        LOGICAL, INTENT(OUT) :: lmin
        DOUBLE PRECISION, INTENT(OUT) :: rf,ru,ron,roff

        !=== VARIABLES ===!

        INTEGER :: i,j,ip,jp,ic,iloop
        INTEGER :: i3,nbp,n5p,n3p,nlp

        DOUBLE PRECISION :: xf,xu,xb
        REAL :: e,el,es,et,emax,emin

        INTEGER :: iwc(4,4)
        INTEGER, PARAMETER :: loop = 3

        DOUBLE PRECISION, PARAMETER :: rcf  = 0.11000000000000000000d+8
        DOUBLE PRECISION, PARAMETER :: cfac = 0.16605389210321896964d-8

        DATA (iwc(1,i),i=1,4) / 0,0,0,1 /
        DATA (iwc(2,i),i=1,4) / 0,0,1,0 /
        DATA (iwc(3,i),i=1,4) / 0,1,0,1 /
        DATA (iwc(4,i),i=1,4) / 1,0,1,0 /


        e = 0.0e0
        emin = 1.0e10
        emax =-1.0e10

        lmin = .true.

        i3 = ibsp(i5)


        !=== Find Hairpin Loop ===!

        iloop = 0

        DO i=i5,i3
        IF ( ibsp(i) < i .and. ibsp(i) > 0 ) THEN
          iloop = ibsp(i)
          EXIT
        ENDIF 
        ENDDO

        !=== Compute Folding Barrier Energy ===!

        i = iloop
        j = ibsp(i)

        nbp = 0
        n5p = 0
        n3p = 0
        nlp = j - i - 1

        !=== Hairpin Energy ===!

        CALL EHAIR (iseq,i,j,n,el)

        e = e + el

        IF ( e > emax ) emax = e
        IF ( e < emin ) emin = e

        !=== Check for Local Min ===!

        IF ( nlp > loop+1 ) THEN
        IF ( iwc(iseq(i+1),iseq(j-1)) == 1 ) THEN

          CALL EHAIR (iseq,i+1,j-1,n,et)
          CALL ESTACK (iseq,i,j,i+1,j-1,n,es)

          et = et + es - el

          IF ( et < 0.0e0 ) lmin = .false.

        ENDIF
        ENDIF

        CALL EHAIR (iseq,i-1,j+1,n,et)
        CALL ESTACK (iseq,i-1,j+1,i,j,n,es)

        et = et - es - el

        IF ( et < 0.0e0 ) lmin = .false.


        !=== Stacking + Bulge Energy ===!

        DO WHILE ( i > i5 )

          !=== Stacking Energy ===!

          ic = 1

          ip = i - 1
          jp = j + 1

          IF ( ip < i5 ) EXIT
          IF ( jp > i3 ) EXIT

          DO WHILE ( ibsp(ip) == jp )

            ic = ic + 1

            CALL ESTACK (iseq,ip,jp,i,j,n,es)

            e = e + es

            IF ( e > emax ) emax = e
            IF ( e < emin ) emin = e

            i = ip
            j = jp

            ip = ip - 1
            jp = jp + 1

            IF ( ip < i5 ) EXIT
            IF ( jp > i3 ) EXIT

          ENDDO

          IF ( nbp == 0 ) nbp = ic

          IF ( ip < i5 ) EXIT
          IF ( jp > i3 ) EXIT

          DO WHILE ( ibsp(ip) == 0 )
            ip = ip - 1
          ENDDO

          jp = ibsp(ip)

          IF ( n5p == 0 .and. n3p == 0 ) THEN
            n5p = i - ip - 1
            n3p = jp - j - 1
          ENDIF

          !=== Bulge Energy ===!

          CALL EBULGE (iseq,ip,jp,i,j,n,el)

          e = e + el

          IF ( e > emax ) emax = e
          IF ( e < emin ) emin = e

          !=== Check for Local Min ===!
          !=== Open/Close Top Bulge ===!

          IF ( ibsp(i-1) == 0 .and. ibsp(j+1) == 0 ) THEN
          IF ( iwc(iseq(i-1),iseq(j+1)) == 1 ) THEN

            IF ( i-ip == 2 .and. jp-j == 2 ) THEN
              CALL ESTACK (iseq,ip,jp,ip+1,jp-1,n,et)
            ELSE 
              CALL EBULGE (iseq,ip,jp,i-1,j+1,n,et)
            ENDIF

            CALL ESTACK (iseq,i-1,j+1,i,j,n,es)

            et = et + es - el

            IF ( et < 0.0e0 ) lmin = .false.

          ENDIF
          ENDIF

          CALL EBULGE (iseq,ip,jp,i+1,j-1,n,et)
          CALL ESTACK (iseq,i,j,i+1,j-1,n,es)

          et = et - es - el

          IF ( et < 0.0e0 ) lmin = .false.

          !=== Check for Local Min ===!
          !=== Open/Close Bottom Bulge ===!

          IF ( ibsp(ip+1) == 0 .and. ibsp(jp-1) == 0 ) THEN
          IF ( iwc(iseq(ip+1),iseq(jp-1)) == 1 ) THEN

            IF ( i-ip == 2 .and. jp-j == 2 ) THEN
              CALL ESTACK (iseq,i-1,j+1,i,j,n,et)
            ELSE 
              CALL EBULGE (iseq,ip+1,jp-1,i,j,n,et)
            ENDIF

            CALL ESTACK (iseq,ip,jp,ip+1,jp-1,n,es)

            et = et + es - el

            IF ( et < 0.0e0 ) lmin = .false.

          ENDIF
          ENDIF

          CALL EBULGE (iseq,ip-1,jp+1,i,j,n,et)
          CALL ESTACK (iseq,ip-1,jp+1,ip,jp,n,es)

          et = et - es - el

          IF ( et < 0.0e0 ) lmin = .false.

          i = ip
          j = jp

        ENDDO

        !=== A-U Penalty ===!

        el = 0.0e0

        IF ( iseq(i5) == 4 ) el = 0.50e0
        IF ( iseq(i3) == 4 ) el = 0.50e0

        e = e + el

        IF ( e > emax ) emax = e
        IF ( e < emin ) emin = e

        !=== Check for Local Min ===!
        !=== Open/Close Last BP ===!

        IF ( i5 > 1 .and. i3 < n .and. i3-i5+2 <= mxbp ) THEN
        IF ( ibsp(i5-1) == 0 .and. ibsp(i3+1) == 0 ) THEN
        IF ( iwc(iseq(i5-1),iseq(i3+1)) == 1 ) THEN

          et = 0.0e0

          IF ( iseq(i5-1) == 4 ) et = 0.50e0
          IF ( iseq(i3+1) == 4 ) et = 0.50e0

          CALL ESTACK (iseq,i5-1,i3+1,i5,i3,n,es)

          et = et + es - el

          IF ( et < 0.0e0 ) lmin = .false.

        ENDIF
        ENDIF
        ENDIF

        et = 0.0e0

        IF ( iseq(i5+1) == 4 ) et = 0.50e0
        IF ( iseq(i3-1) == 4 ) et = 0.50e0

        CALL ESTACK (iseq,i5,i3,i5+1,i3-1,n,es)

        et = et - es - el

        IF ( et < 0.0e0 ) lmin = .false.

        xf = emax
        xu = emax - e

        rf = ratef * DEXP(-beta*xf)
        ru = ratef * DEXP(-beta*xu)


        !=== Compute Affinity for CP ===!

        xf = 1.0d0
        xb = 1.0d0

        IF ( nlp == 4 ) THEN

          !=== Loop ===!

          i = iloop + 3
          j = iloop + 4

          IF ( iseq(i) == 1 ) xb = xb * 1.0d-2
          IF ( iseq(i) == 2 ) xb = xb * 6.0d0
          IF ( iseq(i) == 3 ) xb = xb * 1.0d-2

          IF ( iseq(j) /= 1 ) xb = xb * 1.0d-3

          !=== Bulge ===!

          i = iloop - nbp

          IF ( nbp == 2 ) THEN
            IF ( n5p == 0 ) xf = xf * 1.0d-1
            IF ( n5p == 1 ) THEN
            IF ( iseq(i) /= 1 ) xf = xf * 0.50d0
            ENDIF
          ELSEIF ( nbp == 3 ) THEN
            IF ( n5p == 0 ) xf = xf * 1.0d-1
            IF ( n5p == 1 ) xf = xf * 1.0d-2
          ELSEIF ( nbp == 4 ) THEN
            IF ( n5p == 0 ) xf = xf * 1.0d-1
            IF ( n5p == 1 ) xf = xf * 1.0d-3
          ELSE
            xf = xf * 1.0d-1
          ENDIF

          ron = rcf * cfac / vol
          roff = rcf * DEXP(-beta*12.0d0)

          ron = ron * xf
          roff = roff / xb

          IF ( n5p > 1 ) THEN
            ron = 1.0d-30
            roff = 1.0d10
          ENDIF

          IF ( xf * xb  < 1.0d-6 ) THEN
            ron = 1.0d-30
            roff = 1.0d10
          ENDIF

        ELSEIF ( nlp == 3 ) THEN

          xf = xf / 12.0d0

          !=== Loop ===!

          i = iloop + 2
          j = iloop + 3

          IF ( iseq(i) == 1 ) xb = xb * 1.0d-2
          IF ( iseq(i) == 2 ) xb = xb * 6.0d0
          IF ( iseq(i) == 3 ) xb = xb * 1.0d-2

          IF ( iseq(j) /= 1 ) xb = xb * 1.0d-3

          !=== Bulge ===!

          i = iloop - nbp

          IF ( nbp == 2 ) THEN
            IF ( n5p == 0 ) xf = xf * 1.0d-1
            IF ( n5p == 1 ) xf = xf * 1.0d-2
          ELSEIF ( nbp == 3 ) THEN
            IF ( n5p == 0 ) xf = xf * 1.0d-1
            IF ( n5p == 1 ) THEN
            IF ( iseq(i) /= 1 ) xf = xf * 0.50d0
            ENDIF
          ELSEIF ( nbp == 4 ) THEN
            IF ( n5p == 0 ) xf = xf * 1.0d-1
            IF ( n5p == 1 ) xf = xf * 1.0d-3
          ELSE
            xf = xf * 1.0d-1
          ENDIF

          ron = rcf * cfac / vol
          roff = rcf * DEXP(-beta*12.0d0)

          ron = ron * xf
          roff = roff / xb

          IF ( n5p > 1 ) THEN
            ron = 1.0d-30
            roff = 1.0d10
          ENDIF

          IF ( xf * xb  < 1.0d-6 ) THEN
            ron = 1.0d-30
            roff = 1.0d10
          ENDIF

        ELSE

          !=== Default Affinity ===!

          ron = 1.0d-30
          roff = 1.0d10

        ENDIF
        if ( xu <= 4.0d0 ) lmin = .false.
        !if ( lmin ) write(321,*)rf,ru,emax,e
        RETURN

      END SUBROUTINE RNARATES
