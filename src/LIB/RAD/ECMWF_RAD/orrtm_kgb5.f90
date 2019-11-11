!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD1_MODIF_IDRIS 2003/02/19 13:36:33
!-----------------------------------------------------------------
!OCL SCALAR
SUBROUTINE ORRTM_KGB5

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)
!     Reformatted for F90 by JJMorcrette, ECMWF

!     ------------------------------------------------------------------

!*    kindef: define default KIND macros
! --------------------------------------


USE PARKIND1, ONLY :&
 &JPIT, JPIS, JPIM, JPIB,&
 &JPRT, JPRS, JPRM, JPRB







! --------------------------------------


USE OYOERRTO5 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,&
           &FRACREFBO, CCL4O
USE OYOERRTA5 , ONLY : STRRAT1   ,STRRAT2

!     ------------------------------------------------------------------
!USE MODI_RRTM_KGB5_A
!USE MODI_RRTM_KGB5_B
!USE MODI_RRTM_KGB5_C
!

IMPLICIT NONE

!!!! JDJD Decoupage en 3 parties en raison de pbs de compilation IDRIS
CALL RRTM_KGB5_A
CALL RRTM_KGB5_B
CALL RRTM_KGB5_C
!!!! JDJD Decoupage en 3 parties en raison de pbs de compilation IDRIS
!     -----------------------------------------------------------------
RETURN
END SUBROUTINE ORRTM_KGB5
