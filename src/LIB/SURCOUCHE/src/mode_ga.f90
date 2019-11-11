!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Author: J.Escobar
!
! Modifications:
!  P.Wautelet: 14/12/2018: split from fmwrit_ll.f90
!-----------------------------------------------------------------
#ifdef MNH_GA
MODULE MODE_GA
#include "mafdecls.fh"
#include "global.fh"
    !
    !  Global Array Variables
    !
    INTEGER, PARAMETER                              :: jpix=1 , jpiy = 2 , jpiz = 3
    !
    INTEGER                                         :: NIMAX_ll,NJMAX_ll, IIU_ll,IJU_ll,IKU_ll
    integer                                         :: heap=5*10**6, stack
    logical                                         :: gstatus_ga
    INTEGER, PARAMETER                              :: ndim_GA = 3
    INTEGER, DIMENSION(ndim_GA)                     :: dims_GA , chunk_GA
    INTEGER,PARAMETER                               :: CI=1 ,CJ=-1 ,CK=-1
    INTEGER                                         :: g_a
    integer, DIMENSION(ndim_GA)                     :: lo_col, hi_col , ld_col
    integer, DIMENSION(ndim_GA)                     :: lo_zplan , hi_zplan , ld_zplan
    INTEGER                                         :: NIXO_L,NIXE_L,NIYO_L,NIYE_L
    INTEGER                                         :: NIXO_G,NIXE_G,NIYO_G,NIYE_G

    LOGICAL,SAVE                                    :: GFIRST_GA  = .TRUE.
    INTEGER                                         :: IIU_ll_MAX = -1, IJU_ll_MAX = -1, IKU_ll_MAX = -1

  CONTAINS

    SUBROUTINE MNH_INIT_GA(MY_NI,MY_NJ,MY_NK,HRECFM,HRW_MODE)

!
!  Modification
!  J.Escobar 5/02/2015 : use JPHEXT from MODD_PARAMETERS_ll

      USE MODE_TOOLS_ll,       ONLY : GET_GLOBALDIMS_ll
      USE MODD_PARAMETERS_ll,  ONLY : JPHEXT
      USE MODD_IO_ll,          ONLY : ISP
      USE MODE_GATHER_ll,      ONLY : GET_DOMWRITE_ll
      USE MODE_SCATTER_ll,     ONLY : GET_DOMREAD_ll

      IMPLICIT NONE

      INTEGER,          INTENT(IN) :: MY_NI,MY_NJ,MY_NK
      CHARACTER(LEN=*), INTENT(IN) :: HRECFM   ! name of the article to write
      CHARACTER(LEN=*), INTENT(IN) :: HRW_MODE

      IF ( GFIRST_GA ) THEN
         GFIRST_GA = .FALSE.
         !
         !   Allocate memory for GA library
         !
         stack = heap
         !gstatus_ga = ma_init(MT_F_DBL, stack/ISNPROC, heap/ISNPROC)
         gstatus_ga = ma_init(MT_F_DBL, stack, heap)
         if ( .not. gstatus_ga ) STOP " MA_INIT FAILED "
         !
         !   Initialize GA library
         !
         !call ga_initialize_ltd(100000000)
         call ga_initialize()
      END IF

      CALL GET_GLOBALDIMS_ll (NIMAX_ll,NJMAX_ll)
      IIU_ll = NIMAX_ll + 2*JPHEXT
      IJU_ll = NJMAX_ll + 2*JPHEXT
      IKU_ll = MY_NK
      !
      !   configure Global array dimensions
      !
      dims_GA(JPIX) = IIU_ll
      dims_GA(JPIY) = IJU_ll
      dims_GA(JPIZ) = IKU_ll
      chunk_GA(JPIX)   = CI
      chunk_GA(JPIY)   = CJ
      chunk_GA(JPIZ)   = CK
      IF ( CI .EQ. 1 ) chunk_GA(JPIX)   = dims_GA(JPIX) ! 1 block in X direction
      IF ( CJ .EQ. 1 ) chunk_GA(JPIY)   = dims_GA(JPIY) ! 1 block in Y direction
      IF ( CK .EQ. 1 ) chunk_GA(JPIZ)   = dims_GA(JPIZ) ! 1 block in Z direction
      !
      !   (re)create global array g_a ( if to small create it ... )
      !
      IF ( ( IIU_ll .GT. IIU_ll_MAX ) .OR. ( IJU_ll .GT. IJU_ll_MAX ) .OR. ( IKU_ll .GT. IKU_ll_MAX ) ) THEN
         !
         ! reallocate the g_a , if need with bigger Z size
         !
         IF ( IKU_ll_MAX .NE. -1 ) gstatus_ga =  ga_destroy(g_a)
         IIU_ll_MAX = IIU_ll
         IJU_ll_MAX = IJU_ll
         IKU_ll_MAX = IKU_ll
         gstatus_ga = nga_create(MT_F_DBL, ndim_GA, dims_GA, HRECFM ,chunk_GA, g_a)
         call ga_sync()
      END IF
      !----------------------------------------------------------------------!
      !                                                                      !
      ! Define/describe local column data owned by this process to write     !
      !                                                                      !
      !----------------------------------------------------------------------!
      IF ( HRW_MODE .EQ. "WRITE" ) THEN
      CALL GET_DOMWRITE_ll(ISP,'local',NIXO_L,NIXE_L,NIYO_L,NIYE_L)
      CALL GET_DOMWRITE_ll(ISP,'global',NIXO_G,NIXE_G,NIYO_G,NIYE_G)
      ELSE
      CALL GET_DOMREAD_ll(ISP,NIXO_L,NIXE_L,NIYO_L,NIYE_L)
      CALL GET_DOMREAD_ll(ISP,NIXO_G,NIXE_G,NIYO_G,NIYE_G)
      END IF
      !
      ! portion of data to write/put | read/get by this proc
      !
      lo_col(JPIX) = NIXO_G
      hi_col(JPIX) = NIXE_G

      lo_col(JPIY) = NIYO_G
      hi_col(JPIY) = NIYE_G

      lo_col(JPIZ) = 1
      hi_col(JPIZ) = IKU_ll
      !
      ! declaration size of this local input column array
      !
      ld_col(JPIX) = MY_NI
      ld_col(JPIY) = MY_NJ
      ld_col(JPIZ) = MY_NK
      !
      !-----------------------------------------------------!
      !                                                     !
      !  Size of local ZSLICE_ll Write buffer on I/O proc   !
      !                                                     !
      !-----------------------------------------------------!
      !
      ! declared dimension
      !
      ld_zplan(JPIX) = IIU_ll
      ld_zplan(JPIY) = IJU_ll
      ld_zplan(JPIZ) = 1
      !
      ! write data by Z slide by I/O proc
      !
      lo_zplan(JPIX:JPIY) = 1
      hi_zplan(JPIX) = IIU_ll
      hi_zplan(JPIY) = IJU_ll
      !call ga_sync()
      !
    END SUBROUTINE MNH_INIT_GA

END MODULE MODE_GA

#endif
