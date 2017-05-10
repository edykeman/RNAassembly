! ==============================================================================
! Subroutine: RNA_INIT (VR)
! 
! Purpose: Initializes a viral RNA by identifying RNA stem-loops and
!          their folding and unfolding rates and binding affinites to
!          coat protein. Constructs the inital priority queue for the
!          folding times. Assumes that the RNA is single stranded.
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
! Modules -
! Functions -
! Subroutines - RANDOM, RNARATES
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      SUBROUTINE RNA_INIT (VR)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(VIRAL_RNA), INTENT(INOUT) :: vr

        !=== VARIABLES ===!

        DOUBLE PRECISION :: tau,t1,t2,r,random
        DOUBLE PRECISION :: ron,roff,rf,ru

        INTEGER :: i,j,k,l,n,ip,jp,is
        INTEGER :: ic,jc,i5,i3,nc,ns,nsl

        INTEGER :: iseq(mxnt),istk(mxnt)
        INTEGER :: jstk(mxnt),ibsp(mxnt)
        INTEGER :: icnt(10*mxnt),ifbp(10*mxnt)
        INTEGER :: iwrk(mxnt,10*mxnt)
        INTEGER :: jwrk(mxnt,10*mxnt)
        INTEGER :: kwrk(mxnt,10*mxnt)

        LOGICAL :: ipush,ilone,lmin

        INTEGER :: iwc(4,4)
        INTEGER, PARAMETER :: loop = 3
        character :: fld(mxnt)

        DATA (iwc(1,i),i=1,4) / 0,0,0,1 /
        DATA (iwc(2,i),i=1,4) / 0,0,1,0 /
        DATA (iwc(3,i),i=1,4) / 0,1,0,1 /
        DATA (iwc(4,i),i=1,4) / 1,0,1,0 /


        n = vr% n
        nsl = 0

        vr% nsl = 0

        vr% i5 = 0
        vr% i3 = 0
        vr% ip = 0
        vr% ib = 0
        vr% nt = 0

        vr% ipres(:) = 0
        vr% linkp(:) = 0
        vr% linkb(:) = 0

        vr% itree(:,:) = 0
        vr% lstem(:,:) = 0

        vr% info_sl(:) = 0
        vr% info_cp(:) = 0
        vr% inxt_sl(:) = 0
        vr% inxt_cp(:) = 0

        vr% ique_sl(:,:) = 0
        vr% ique_cp(:,:) = 0

        vr% rfold(:,:) = 0.0d0
        vr% rbind(:,:) = 0.0d0

        vr% tfold(:) = tinf
        vr% tbind(:) = tinf
        vr% tsave(:) = tinf
        vr% tcoat(:) = tinf

        vr% tnxt_sl(:,:) = tinf
        vr% tnxt_cp(:,:) = tinf

        vr% tss(:) = time
        vr% trm(:) = time
        vr% tub(:) = time
        vr% t(:) = tinf


        !=== Find All Stem-Loops ===!

        is = 1

        icnt(is) = 1
        ifbp(is) = 0

        iwrk(1,is) = 1
        jwrk(1,is) = n
        kwrk(:,is) = 0


        !=== Pull Structure From Stack ===!

 1      IF ( is == 0 ) GOTO 3

        ipush = .false.

        ic = icnt(is)
        ip = ifbp(is)

        istk(:) = iwrk(:,is)
        jstk(:) = jwrk(:,is)
        ibsp(:) = kwrk(:,is)

        is = is - 1 

        !=== Found Stem-Loop ===!

        IF ( ic == 0 ) THEN

          !=== Don't Output SS Structure ===!

          IF ( ip == 0 ) GOTO 1

          !=== Store Stem-Loop ===!

          jp = ibsp(ip)

          CALL RNARATES (vr%iseq,ibsp,n,ip,rf,ru,ron,roff,lmin)

          IF ( lmin ) THEN

            nsl = nsl + 1

            vr% lstem(1,nsl) = ip
            vr% lstem(2,nsl) = jp

            vr% rbind(1,nsl) = ron
            vr% rbind(2,nsl) = roff

            vr% rfold(1,nsl) = rf
            vr% rfold(2,nsl) = ru

          ENDIF

          GOTO 1

        ENDIF

        i = istk(ic)
        j = jstk(ic)

        ic = ic - 1

        !=== Check for Isolated BP ===!

        ilone = .false.

        IF ( i /= 1 .and. j /= n ) THEN

          i5 = i - 1
          i3 = ibsp(i5)

          ilone = .false.

          IF ( i3 == j+1 ) THEN

            ilone = .true.

            IF ( i5 > 1 ) THEN
            IF ( ibsp(i5-1) == i3+1 ) ilone = .false.
            ENDIF

          ENDIF

        ENDIF

        IF ( ilone ) THEN
        IF ( j-i <= loop ) GOTO 3
        ELSE
        IF ( j-i <= loop ) GOTO 2
        ENDIF

        !=== CASE 1: i and j basepair ===!

        IF ( j-i <= mxbp ) THEN
        IF ( iwc(vr%iseq(i),vr%iseq(j)) == 1 ) THEN

          is = is + 1

          icnt(is) = ic + 1

          IF ( ip == 0 ) ifbp(is) = i
          IF ( ip /= 0 ) ifbp(is) = ip

          iwrk(:,is) = istk(:)
          jwrk(:,is) = jstk(:)
          kwrk(:,is) = ibsp(:)

          iwrk(ic+1,is) = i+1
          jwrk(ic+1,is) = j-1

          kwrk(i,is) = j
          kwrk(j,is) = i

          ipush = .true.

        ENDIF
        ENDIF

        IF ( ilone ) GOTO 3

        !=== CASE 2: j single stranded ===!

        is = is + 1

        icnt(is) = ic + 1
        ifbp(is) = ip

        iwrk(:,is) = istk(:)
        jwrk(:,is) = jstk(:)
        kwrk(:,is) = ibsp(:)

        iwrk(ic+1,is) = i
        jwrk(ic+1,is) = j-1

        ipush = .true.

        !=== CASE 3: Internal Loop ===!

        DO k=i+1,j-1

          IF ( j-k >  mxbp ) CYCLE
          IF ( j-k <= loop ) CYCLE
          IF ( iwc(vr%iseq(j),vr%iseq(k)) == 0 ) CYCLE

          is = is + 1

          icnt(is) = ic + 1

          IF ( ip == 0 ) ifbp(is) = k
          IF ( ip /= 0 ) ifbp(is) = ip

          iwrk(:,is) = istk(:)
          jwrk(:,is) = jstk(:)
          kwrk(:,is) = ibsp(:)

          iwrk(ic+1,is) = k+1
          jwrk(ic+1,is) = j-1

          kwrk(k,is) = j
          kwrk(j,is) = k

          ipush = .true.

        ENDDO

 2      IF ( .not. ipush ) THEN

          is = is + 1

          icnt(is) = ic
          ifbp(is) = ip

          iwrk(:,is) = istk(:)
          jwrk(:,is) = jstk(:)
          kwrk(:,is) = ibsp(:)

        ENDIF

 3      IF ( is /= 0 ) GOTO 1


        !=== Compute Folding Times ===!

        DO i=1,nsl

          !r = 0.904792018d0 * RANDOM(iseed) + 0.0000454d0
          r = RANDOM(iseed)

          rf = vr% rfold(1,i)

          tau = -DLOG(r) / rf

          vr% tfold(i) = tau

          vr% inxt_sl(i) = 1
          vr% tnxt_sl(0,i) = time + tau
          vr% tnxt_sl(1,i) = tinf

        ENDDO

        !=== Construct SL Queue ===!

        ns = 2
        DO WHILE ( ns < nsl )
        ns = 2 * ns
        ENDDO

        vr% nsl = nsl
        vr% ns = ns

        DO i=1,ns,2

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

        ENDDO

        k = 4
        l = 1

        DO WHILE ( k <= ns )

          j = k / 2

          DO i=j,ns,k

            !=== RNA Folding Queue ===!

            ic = vr% ique_sl(0,i-l)
            jc = vr% ique_sl(0,i+l)

            t1 = vr% tnxt_sl(0,ic)
            t2 = vr% tnxt_sl(0,jc)

            IF ( t1 <= t2 ) vr% ique_sl(0,i) = ic
            IF ( t2 <  t1 ) vr% ique_sl(0,i) = jc

            !=== RNA Binding Queue ===!

            ic = vr% ique_sl(1,i-l)
            jc = vr% ique_sl(1,i+l)

            t1 = vr% tnxt_sl(1,ic)
            t2 = vr% tnxt_sl(1,jc)

            IF ( t1 <= t2 ) vr% ique_sl(1,i) = ic
            IF ( t2 <  t1 ) vr% ique_sl(1,i) = jc

          ENDDO

          l = l * 2
          k = k * 2

        ENDDO


        !=== Construct CP Queue ===!

        nc = 2
        DO WHILE ( nc < ncp )
        nc = 2 * nc
        ENDDO

        vr% nc = nc

        DO i=1,nc,2

          !=== CP Removal Queue ===!

          t1 = vr% tnxt_cp(0,i+0)
          t2 = vr% tnxt_cp(0,i+1)

          IF ( t1 <= t2 ) vr% ique_cp(0,i) = i
          IF ( t2 <  t1 ) vr% ique_cp(0,i) = i+1

          !=== CP Binding Queue ===!

          t1 = vr% tnxt_cp(1,i+0)
          t2 = vr% tnxt_cp(1,i+1)

          IF ( t1 <= t2 ) vr% ique_cp(1,i) = i
          IF ( t2 <  t1 ) vr% ique_cp(1,i) = i+1

        ENDDO

        k = 4
        l = 1

        DO WHILE ( k <= nc )

          j = k / 2

          DO i=j,nc,k

            !=== CP Removal Queue ===!

            ic = vr% ique_cp(0,i-l)
            jc = vr% ique_cp(0,i+l)

            t1 = vr% tnxt_cp(0,ic)
            t2 = vr% tnxt_cp(0,jc)

            IF ( t1 <= t2 ) vr% ique_cp(0,i) = ic
            IF ( t2 <  t1 ) vr% ique_cp(0,i) = jc

            !=== CP Binding Queue ===!

            ic = vr% ique_cp(1,i-l)
            jc = vr% ique_cp(1,i+l)

            t1 = vr% tnxt_cp(1,ic)
            t2 = vr% tnxt_cp(1,jc)

            IF ( t1 <= t2 ) vr% ique_cp(1,i) = ic
            IF ( t2 <  t1 ) vr% ique_cp(1,i) = jc

          ENDDO

          l = l * 2
          k = k * 2

        ENDDO

        !=== Store Min Queue Times ===!

        k = ns / 2

        ic = vr% ique_sl(0,k)
        jc = vr% ique_sl(1,k)

        vr% t(0) = vr% tnxt_sl(0,ic)
        vr% t(1) = vr% tnxt_sl(1,jc)

        RETURN

      END SUBROUTINE RNA_INIT
