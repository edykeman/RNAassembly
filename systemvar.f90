! ==============================================================================
! Module: SYSTEMVAR
! 
! Purpose: Contains the global variables needed in Vassembly.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            04/01/2016   Original Code
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

      MODULE SYSTEMVAR

        IMPLICIT NONE


        !=== ALLOCATABLE ARRAYS ===!

        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: gb

        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: gbond

        INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: ique
        INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: lbond
        INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: mapnn
        INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: maphp

        INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nncp


        !=== VARIABLES ===!

        DOUBLE PRECISION, SAVE :: rate_ab,rate_cc,rate_au,ratef
        DOUBLE PRECISION, SAVE :: rkap_ab,rkap_cc,rkap_au,rkapf
        DOUBLE PRECISION, SAVE :: vol,beta,temp,time,tint,tfinal,TADD

        INTEGER, SAVE :: nbond,ncp,nps,nnps,nnmax
        INTEGER, SAVE :: npro,nrna,nque,iseed,IADD


        !=== PARAMETERS ===!

        !=== Gas Constant in kcal / (mol*K) ===!
        DOUBLE PRECISION, PARAMETER :: gcons = 1.987206d-3
        DOUBLE PRECISION, PARAMETER :: tinf  = 1.000000d60

        REAL, PARAMETER :: em  = 3.40e0
        REAL, PARAMETER :: eh  = 0.40e0
        REAL, PARAMETER :: es  = 0.00e0
        REAL, PARAMETER :: eau = 0.50e0
        REAL, PARAMETER :: einf= 1.0e10

        INTEGER, PARAMETER :: mxbp = 20
        INTEGER, PARAMETER :: mxcp = 64    !(Power of 2)
        INTEGER, PARAMETER :: mxnt = 512   !(Power of 2)
        INTEGER, PARAMETER :: mxsl = 2048  !(Power of 2)

      END MODULE
