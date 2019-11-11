!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 surfex 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_MNHPUT_ZS_n
!     #######################
INTERFACE
      SUBROUTINE MNHPUT_ZS_n
!
!
END SUBROUTINE MNHPUT_ZS_n
!
END INTERFACE
END MODULE MODI_MNHPUT_ZS_n
!
!     ########################################
      SUBROUTINE MNHPUT_ZS_n
!     ########################################
!
!!****  *MNHPUT_ZS_n* - routine to modify surface oropgraphy using atmospheric
!                    model orography
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!  06/2016     (G.Delautier) phasage surfex 8
!!  01/2018      (G.Delautier) SURFEX 8.1
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll
!
USE MODD_CONF,       ONLY : CPROGRAM
USE MODD_PARAMETERS, ONLY : JPHEXT
USE MODD_DIM_n,      ONLY : NIMAX, NJMAX
USE MODD_GRID_n,     ONLY : XZS
!
USE MODI_PUT_ZS_N
USE MODD_MNH_SURFEX_n
!
USE MODI_GET_LUOUT
!
USE MODE_MPPDB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER                         :: IIB, IIE, IJB, IJE
INTEGER                         :: IL
REAL, DIMENSION(:), ALLOCATABLE :: ZZS
INTEGER :: ILUOUT
!-------------------------------------------------------------------------------
!
CALL GET_LUOUT(CPROGRAM,ILUOUT)
SELECT CASE(CPROGRAM)
  CASE ('NESPGD')
    IIB = JPHEXT + 1
    IIE = JPHEXT + NIMAX
    IJB = JPHEXT + 1
    IJE = JPHEXT + NJMAX
  CASE DEFAULT
    CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
END SELECT
!
IL = (IIE-IIB+1)*(IJE-IJB+1)
!
ALLOCATE(ZZS(IL))
ZZS(:) = RESHAPE (XZS(IIB:IIE,IJB:IJE), (/ IL /) )
!
CALL PUT_ZS_n(YSURF_CUR%FM%F,YSURF_CUR%IM%S,YSURF_CUR%SM%S,&
              YSURF_CUR%U,YSURF_CUR%TM%TOP,YSURF_CUR%WM%W,'MESONH',IL,ZZS(:))
!
CALL MPPDB_CHECK_SURFEX2D(ZZS,"mnhput_zs_n:ZZS",PRECISION,ILUOUT)
CALL MPPDB_CHECK2D(XZS,"mnhput_zs_n:MODD_GRID_n::XZS",PRECISION)
!
DEALLOCATE(ZZS)
!
!==============================================================================
!
END SUBROUTINE MNHPUT_ZS_n
