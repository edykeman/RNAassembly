! ==============================================================================
! Subroutine: READDATA (VRNA)
! 
! Purpose: Reads in the list of "bond" interactions between capsid proteins
!          in the capsid, the RNA binding rates, and the hamiltonian path map.
!
! Method:
!
! Arguments:
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            04/01/2016   Original Code
!
! Dependancies:
!
! Modules - SYSTEMVAR, CLASS_VRNA
! Functions -
! Subroutines -
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      SUBROUTINE READDATA (VRNA)

        USE SystemVar
        USE Class_VRNA

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE (VIRAL_RNA), INTENT(INOUT) :: vrna(nrna)

        !=== VARIABLES ===!

        INTEGER :: i,j,k,n,jj,kk,nn
        CHARACTER (LEN=70) :: fmat
        CHARACTER (LEN=1024) :: fasta

        DOUBLE PRECISION :: x,rcf,rcb

        DOUBLE PRECISION, PARAMETER :: cfac = 0.16605389210321896964d-8


        !=== Set Beta ===!

        beta = 1.0d0 / ( temp * gcons )


        !=== Set rkap_ab rkap_cc and rkap_au ===!
        !=== convert from 1/M*s to 1/s ===!

        rkap_ab = rate_ab

        rkap_cc = rate_cc * cfac / vol

        rkap_au = rate_au * cfac / vol


        !=== Read in Parameter File ===!

        READ(1,*)ncp,nbond

        ALLOCATE (lbond(2,nbond),gbond(nbond))

        READ(1,*)(lbond(1,i),lbond(2,i),gbond(i),i=1,nbond)


        !=== Read In Hamiltonian Map ===!

        READ(2,*)nps,nnps

        ALLOCATE (maphp(nnps,nps))

        fmat = '(15I5)'

        DO i=1,nps

          READ(2,fmat)(maphp(j,i),j=1,nnps)

        ENDDO


        !=== Form Neighbormap ===!

        ALLOCATE (nncp(ncp))

        nnmax = 0
        nncp(:) = 0

        DO i=1,nbond

          j = lbond(1,i)
          k = lbond(2,i)

          nncp(j) = nncp(j) + 1
          nncp(k) = nncp(k) + 1

        ENDDO

        DO i=1,ncp
        IF ( nncp(i) > nnmax ) THEN
          nnmax = nncp(i)
        ENDIF
        ENDDO

        ALLOCATE (mapnn(nnmax,ncp))
        ALLOCATE (gb(nnmax,ncp))

        nncp(:) = 0

        DO i=1,nbond

          j = lbond(1,i)
          k = lbond(2,i)

          nncp(j) = nncp(j) + 1
          nncp(k) = nncp(k) + 1

          jj = nncp(j)
          kk = nncp(k)

          mapnn(jj,j) = k
          mapnn(kk,k) = j

          gb(jj,j) = gbond(i)
          gb(kk,k) = gbond(i)

        ENDDO


        !=== Read in RNA Sequences ===!

        fmat = '(1024A1)'

        READ(3,*)fasta
        n = LEN_TRIM(fasta)

        DO i=1,nrna

!          READ(3,*)fasta
!
!          n = LEN_TRIM(fasta)

          READ(fasta,fmat)(vrna(i)% seq(j),j=1,n)

          nn = 1
          DO WHILE ( nn < n )
          nn = 2 * nn
          ENDDO

          vrna(i)% n = n
          vrna(i)% nn = nn

          DO j=1,nn
          IF ( vrna(i)% seq(j) == 'A' ) vrna(i)% iseq(j) = 1
          IF ( vrna(i)% seq(j) == 'a' ) vrna(i)% iseq(j) = 1
          IF ( vrna(i)% seq(j) == 'C' ) vrna(i)% iseq(j) = 2
          IF ( vrna(i)% seq(j) == 'c' ) vrna(i)% iseq(j) = 2
          IF ( vrna(i)% seq(j) == 'G' ) vrna(i)% iseq(j) = 3
          IF ( vrna(i)% seq(j) == 'g' ) vrna(i)% iseq(j) = 3
          IF ( vrna(i)% seq(j) == 'T' ) vrna(i)% iseq(j) = 4
          IF ( vrna(i)% seq(j) == 't' ) vrna(i)% iseq(j) = 4
          IF ( vrna(i)% seq(j) == 'U' ) vrna(i)% iseq(j) = 4
          IF ( vrna(i)% seq(j) == 'u' ) vrna(i)% iseq(j) = 4
          ENDDO

        ENDDO

        RETURN

      END SUBROUTINE READDATA
