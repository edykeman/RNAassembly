! ==============================================================================
! Module: CLASS_VRNA
! 
! Purpose:
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            05/01/2016   Original Code
!
! Contains:
!
! Modules -
! Functions -
! Subroutines -
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      MODULE CLASS_VRNA

        USE SystemVar, ONLY : mapnn,nncp,ncp,npro,maphp,nbond,nnps,&
                            & nps,mxnt,mxsl,mxcp,mxbp,iseed,gb,beta,&
                            & rkap_au,rkap_ab,rkap_cc,rate_au,rate_cc,&
                            & rate_ab,time,tint,tinf

        IMPLICIT NONE

        PRIVATE

        PUBLIC :: RNA_INIT, FOLD_FIRE, CAPSID_FIRE

        TYPE, PUBLIC :: VIRAL_RNA

          CHARACTER :: seq(mxnt)

          INTEGER :: iseq(mxnt)

          INTEGER :: n   !number of nucleotides
          INTEGER :: nsl !number of sl
          INTEGER :: nn  !number of entries in sl tree (power of 2 with nn >= n)
          INTEGER :: ns  !number of queued stem-lp reactions (power of 2 with ns >= nsl)
          INTEGER :: nc  !number of queued capsid reactions (power of 2 with nc >= ncp)

          INTEGER :: i5,i3
          INTEGER :: ip,ib
          INTEGER :: nt

          INTEGER :: ipres(mxnt) !list of present sl
          INTEGER :: linkp(mxsl) !linked list between present sl
          INTEGER :: linkb(mxsl) !linked list between CPbound sl

          INTEGER :: itree(2,mxnt) !binary tree of present (or CP bound) sl
          INTEGER :: lstem(2,mxsl) !list of stem loop start and stop

          INTEGER :: info_sl(mxsl) !state information for the sl
          INTEGER :: info_cp(mxcp) !state information for the cp
          INTEGER :: inxt_sl(mxsl) !next reaction type for the sl
          INTEGER :: inxt_cp(mxcp) !next reaction type for the cp

          INTEGER :: ique_sl(0:1,mxsl)
          INTEGER :: ique_cp(0:1,mxcp)

          DOUBLE PRECISION :: rfold(2,mxsl)
          DOUBLE PRECISION :: rbind(2,mxsl)

          DOUBLE PRECISION :: tfold(mxsl)
          DOUBLE PRECISION :: tbind(mxsl)
          DOUBLE PRECISION :: tsave(mxsl) !array to save unbinding reaction times
          DOUBLE PRECISION :: tcoat(mxcp)

          DOUBLE PRECISION :: tnxt_sl(0:1,mxsl)
          DOUBLE PRECISION :: tnxt_cp(0:1,mxcp)

          DOUBLE PRECISION :: tss(mxnt),trm(mxcp),tub(mxsl)
          DOUBLE PRECISION :: t(0:1)

        END TYPE VIRAL_RNA

        CONTAINS

        INCLUDE 'rna_init.f90'

        INCLUDE 'fold_fire.f90'
        INCLUDE 'fold_reac.f90'

        INCLUDE 'capsid_fire.f90'
        INCLUDE 'capsid_reac.f90'
        INCLUDE 'capsid_nucl.f90'

        INCLUDE 'add_sl.f90'
        INCLUDE 'del_sl.f90'
        INCLUDE 'get_sl.f90'

        INCLUDE 'requeue_sl.f90'
        INCLUDE 'requeue_cp.f90'

      END MODULE CLASS_VRNA
