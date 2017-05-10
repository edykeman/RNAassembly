! ==============================================================================
! Program: VASSEMBLY
!
! Discription: Performs a "black box" integration of a system of
!              chemical reactions which describe virus assembly.
!
! NOTES:
!        FILE TREE
!
!        Unit = 1  --- Parameter File
!        Unit = 2  --- Hamiltonian Path File
!        Unit = 3  --- RNA Data File
!        Unit = 4  --- Option File
!        Unit = 5  --- Standard Out (NOT USED)
!        Unit = 6  --- Standard In (NOT USED)
!        Unit = 7  --- Output File
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
! Subroutines - READDATA
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      PROGRAM VASSEMBLY

        USE SystemVar
        USE Class_VRNA

        IMPLICIT NONE

        !=== Variable Declaration ===!

        TYPE(VIRAL_RNA), DIMENSION(:), ALLOCATABLE :: vrna

        INTEGER :: i,j,k,n,io,i1,i2
        INTEGER :: iargc,narg

        DOUBLE PRECISION :: dt,t1,t2,tout

        CHARACTER(LEN=80) :: optfile,parmfile,rnafile
        CHARACTER(LEN=80) :: hpfile,outfile,arg


        !=== WELCOME ===!

        WRITE(*,*)'               WELCOME TO VASSEMBLY                 '
        WRITE(*,*)' '
        WRITE(*,*)'This program integrates a set of chemical reactions '
        WRITE(*,*)'between capsid intermediates as a function of time.'
        WRITE(*,*)' '

        !=== Default names ===!

        optfile = 'options.in'

        parmfile = 'capsid.parm'
        rnafile  = 'capsid.rna'
        hpfile   = 'capsid.hmap'
        outfile  = 'capsid.out'


        narg = IARGC ()

        DO i=1,narg,2

          CALL GETARG (i,arg)

          SELECT CASE (arg)

            CASE ('-o')
              CALL GETARG (i+1,optfile)
            CASE DEFAULT

              WRITE(*,*)arg,'Invalid Line Argument'

              STOP

          END SELECT

        ENDDO


        !=== SECTION 0 ===!

        !=== Read in Options ===!

        OPEN (UNIT = 4,FILE = optfile,STATUS='Unknown')

        READ(4,*)parmfile
        READ(4,*)hpfile
        READ(4,*)rnafile
        READ(4,*)outfile

        READ(4,*)nrna
        READ(4,*)npro
        READ(4,*)vol
        READ(4,*)temp
        READ(4,*)rate_ab
        READ(4,*)rate_cc
        READ(4,*)rate_au
        READ(4,*)ratef
        READ(4,*)time
        READ(4,*)tfinal
        READ(4,*)iseed

        !=== Set Internal Time ===!

        tint = time

        CLOSE (UNIT = 4)

        ALLOCATE (vrna(nrna))


        !=== SECTION 1 ===!

        !=== Read In Parameters ===!

        OPEN (UNIT = 1,FILE = parmfile,STATUS='Unknown')
        OPEN (UNIT = 2,FILE =   hpfile,STATUS='Unknown')
        OPEN (UNIT = 3,FILE =  rnafile,STATUS='Unknown')

        WRITE(*,*)'Reading in system data.'

        CALL READDATA (vrna)

        CLOSE (UNIT = 1)
        CLOSE (UNIT = 2)
        CLOSE (UNIT = 3)


        !=== SECTION 2 ===!

        io = 1
        dt = 1.0d-6
        tout = dt

        IF ( time >= dt ) THEN

          tout = DLOG(time) / DLOG(10.0d0)

          i = FLOOR(tout)

          dt = 10.0d0 ** i

          io = NINT(time/dt)

          io = io + 1
          tout = DBLE(io) * dt

        ENDIF


        !=== Initialize RNAs ===!

        DO i=1,nrna
        CALL RNA_INIT (vrna(i))
        ENDDO

        !=== Construct Queue ===!

        nque = 2
        DO WHILE ( nque < nrna )
        nque = 2 * nque
        ENDDO

        ALLOCATE (ique(0:1,nque))

        ique(:,:) = 0

        DO i=1,nque,2

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

        ENDDO

        k = 4
        n = 1

        DO WHILE ( k <= nque )

          j = k / 2

          DO i=j,nque,k

            !=== RNA Folding Queue ===!

            i1 = ique(0,i-n)
            i2 = ique(0,i+n)

            t1 = vrna(i1)% t(0)
            t2 = vrna(i2)% t(0)

            IF ( t1 <= t2 ) ique(0,i) = i1
            IF ( t2 <  t1 ) ique(0,i) = i2

            !=== RNA Binding Queue ===!

            i1 = ique(1,i-n)
            i2 = ique(1,i+n)

            t1 = vrna(i1)% t(1)
            t2 = vrna(i2)% t(1)

            IF ( t1 <= t2 ) ique(1,i) = i1
            IF ( t2 <  t1 ) ique(1,i) = i2

          ENDDO

          n = n * 2
          k = k * 2

        ENDDO


        !=== BEGIN STOCHASTIC SIMULATION ===!

        WRITE(*,*)'Starting SSA ...'

        DO WHILE ( time < tfinal )

          CALL SNREACTION (vrna,tout)

          !=== Increment tout ===!

          IF ( time > tout ) THEN

            tout = tout + dt

            io = io + 1

            IF ( io > 9 ) THEN
              io = 1
              dt = dt * 10.0d0
            ENDIF

          ENDIF

        ENDDO

        WRITE(*,*)'Stochastic Simulation complete.'

        !=== SECTION 4 ===!

        DEALLOCATE (vrna,ique)

      END PROGRAM VASSEMBLY
