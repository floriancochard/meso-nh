!!
!!    #####################
      MODULE MODN_BLOWSNW
!!    #####################
!!
!!*** *MODN_SNW*
!!
!!    PURPOSE
!!    -------
!       Namelist for snow drift scheme
!!
!!**  AUTHOR
!!    ------
!!    V. Vionnet   *CNRM*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 03/2010
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_BLOWSNW_SURF, ONLY : LSNOW_SALT, XCONC_SALT, XSNOW_PRECIP, LSNOW_PRECIP, &
                           LSNOW_WIND, XSNOW_ROUGHNESS, XSNOW_BUOYANCY,        &
                           XEMIRADIUS_SNW,XEMIALPHA_SNW,CSNOW_SALT,   &
                           CSNOW_SEDIM,LBLOWSNW_CANOSUBL,LBLOWSNW_CANODIAG,    &
                           XRSNOW_SBL,LBLOWSNW_ADV

!!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
SAVE
NAMELIST /NAM_SURF_BLOWSNW/  &
     LSNOW_SALT, XCONC_SALT, XSNOW_PRECIP, LSNOW_PRECIP, LSNOW_WIND,  &   !Parameterization type
     XSNOW_ROUGHNESS, XSNOW_BUOYANCY,XEMIRADIUS_SNW,XEMIALPHA_SNW ,   &
     CSNOW_SALT, CSNOW_SEDIM,LBLOWSNW_CANOSUBL,LBLOWSNW_CANODIAG,  &
     XRSNOW_SBL,LBLOWSNW_ADV

!
END MODULE MODN_BLOWSNW
