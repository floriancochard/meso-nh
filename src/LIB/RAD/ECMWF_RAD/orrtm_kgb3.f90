!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD1_MODIF_IDRIS 2003/02/19 13:36:32
!-----------------------------------------------------------------
!OCL SCALAR
SUBROUTINE ORRTM_KGB3

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)
!     Reformatted for F90 by JJMorcrette, ECMWF

!     ------------------------------------------------------------------

!*    kindef: define default KIND macros
! --------------------------------------


USE PARKIND1, ONLY :&
 &JPIT, JPIS, JPIM, JPIB,&
 &JPRT, JPRS, JPRM, JPRB





! --------------------------------------


USE OYOERRTO3 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
           &FRACREFBO  ,FORREFO    ,ABSN2OAO   ,ABSN2OBO
USE OYOERRTA3 , ONLY : ABSN2OA ,ABSN2OB  ,ETAREF    ,H2OREF     ,&
           &N2OREF     ,CO2REF     ,STRRAT

!     ------------------------------------------------------------------
!USE MODI_RRTM_KGB3_A
!USE MODI_RRTM_KGB3_B
!USE MODI_RRTM_KGB3_C
!

IMPLICIT NONE
!
!! JDJDJD Decoupage en raison de pbs de compilation IDRIS
CALL RRTM_KGB3_A
CALL RRTM_KGB3_B
CALL RRTM_KGB3_C
!! JDJDJD Decoupage en raison de pbs de compilation IDRIS

!     -----------------------------------------------------------------
RETURN
END SUBROUTINE ORRTM_KGB3
