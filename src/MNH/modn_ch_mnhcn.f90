!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/cvsroot/MNH-VX-Y-Z/src/MNH/modn_ch_mnhcn.f90,v $ $Revision: 1.2.4.1.2.1.12.2 $ $Date: 2014/01/09 15:01:56 $
!-----------------------------------------------------------------
!!    #####################
      MODULE MODN_CH_MNHC_n
!!    #####################
!!
!!*** *MODN_CH_MNHC$n*
!!
!!    PURPOSE
!!    -------
!       Namelist for parameters that control the chemical part of MesoNH
!!
!!**  AUTHOR
!!    ------
!!    K. Suhre      *Laboratoire d'Aerologie*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 24/05/95
!!    27/07/96 (K. Suhre) restructured
!!    30/11/99 (K. Suhre) add new parameters
!!    13/07/03 (J.-P. Pinty) add flag for lightning production of NOx
!!    15/05/07 (M. Leriche) add logical for aqueous phase chemistry
!!    30/05/07 (M. Leriche) add logical and real for pH calculation
!!    25/04/08 (M. Leriche) add threshold for aqueous phase chemistry
!!    16/09/10 (M. Leriche) add logical for ice phase chemistry
!!    13/01/11 (M. Leriche) add logical for retention in ice
!!    24/03/16 (M. Leriche) remove surface option -> manage them in SURFEX
!!    01/10/16 (F. Brosse) add production/destruction terms computation
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_MNHC_n, ONLY: &
         LUSECHEM_n => LUSECHEM, &
         LUSECHAQ_n => LUSECHAQ, &
         LUSECHIC_n => LUSECHIC, &
         LCH_INIT_FIELD_n => LCH_INIT_FIELD, &
         LCH_CONV_SCAV_n => LCH_CONV_SCAV, &
         LCH_CONV_LINOX_n => LCH_CONV_LINOX, &
         LCH_PH_n => LCH_PH, &
         LCH_RET_ICE_n=>LCH_RET_ICE, &
         XCH_PHINIT_n => XCH_PHINIT, &
         XRTMIN_AQ_n => XRTMIN_AQ,  &
         CCHEM_INPUT_FILE_n => CCHEM_INPUT_FILE, &
         CCH_TDISCRETIZATION_n => CCH_TDISCRETIZATION, &
         NCH_SUBSTEPS_n => NCH_SUBSTEPS, &
         LCH_TUV_ONLINE_n => LCH_TUV_ONLINE, &
         CCH_TUV_LOOKUP_n => CCH_TUV_LOOKUP, &
         CCH_TUV_CLOUDS_n => CCH_TUV_CLOUDS, &
         XCH_TUV_ALBNEW_n => XCH_TUV_ALBNEW, &
         XCH_TUV_DOBNEW_n => XCH_TUV_DOBNEW, &
         XCH_TUV_TUPDATE_n => XCH_TUV_TUPDATE, &
         CCH_VEC_METHOD_n => CCH_VEC_METHOD, &
         NCH_VEC_LENGTH_n => NCH_VEC_LENGTH, &
         XCH_TS1D_TSTEP_n => XCH_TS1D_TSTEP, &
         CCH_TS1D_COMMENT_n => CCH_TS1D_COMMENT, &
         CCH_TS1D_FILENAME_n => CCH_TS1D_FILENAME, &
         CSPEC_PRODLOSS_n => CSPEC_PRODLOSS,  & 
         CSPEC_BUDGET_n => CSPEC_BUDGET
!
IMPLICIT NONE
!
LOGICAL  :: LUSECHEM
LOGICAL  :: LUSECHAQ
LOGICAL  :: LUSECHIC
LOGICAL  :: LCH_INIT_FIELD
LOGICAL  :: LCH_CONV_SCAV
LOGICAL  :: LCH_CONV_LINOX
LOGICAL  :: LCH_PH
lOGICAL  :: LCH_RET_ICE
REAL  :: XCH_PHINIT
REAL  :: XRTMIN_AQ
CHARACTER(LEN=80)  :: CCHEM_INPUT_FILE
CHARACTER(LEN=10)  :: CCH_TDISCRETIZATION
INTEGER  :: NCH_SUBSTEPS
LOGICAL  :: LCH_TUV_ONLINE
CHARACTER*80  :: CCH_TUV_LOOKUP
CHARACTER*4  :: CCH_TUV_CLOUDS
REAL  :: XCH_TUV_ALBNEW
REAL  :: XCH_TUV_DOBNEW
REAL  :: XCH_TUV_TUPDATE
CHARACTER*3  :: CCH_VEC_METHOD
INTEGER  :: NCH_VEC_LENGTH
REAL  :: XCH_TS1D_TSTEP
CHARACTER*80  :: CCH_TS1D_COMMENT
CHARACTER*80  :: CCH_TS1D_FILENAME
CHARACTER(LEN=1024)   :: CSPEC_PRODLOSS 
CHARACTER(LEN=1024)   :: CSPEC_BUDGET 
!
NAMELIST/NAM_CH_MNHCn/LUSECHEM,LUSECHAQ,LUSECHIC,LCH_INIT_FIELD,LCH_CONV_SCAV,&
                      LCH_CONV_LINOX,LCH_PH,LCH_RET_ICE,XCH_PHINIT,XRTMIN_AQ, &
                      CCHEM_INPUT_FILE,CCH_TDISCRETIZATION,NCH_SUBSTEPS,      &
                      LCH_TUV_ONLINE,CCH_TUV_LOOKUP,CCH_TUV_CLOUDS,           &
                      XCH_TUV_ALBNEW,XCH_TUV_DOBNEW,XCH_TUV_TUPDATE,          &
                      CCH_VEC_METHOD,NCH_VEC_LENGTH,XCH_TS1D_TSTEP,           &
                      CCH_TS1D_COMMENT,CCH_TS1D_FILENAME,CSPEC_PRODLOSS,CSPEC_BUDGET
!
CONTAINS
!
SUBROUTINE INIT_NAM_CH_MNHCn
  LUSECHEM = LUSECHEM_n
  LUSECHAQ = LUSECHAQ_n
  LUSECHIC = LUSECHIC_n
  LCH_INIT_FIELD = LCH_INIT_FIELD_n
  LCH_CONV_SCAV = LCH_CONV_SCAV_n
  LCH_CONV_LINOX = LCH_CONV_LINOX_n
  LCH_PH = LCH_PH_n
  LCH_RET_ICE = LCH_RET_ICE_n
  XCH_PHINIT = XCH_PHINIT_n
  XRTMIN_AQ = XRTMIN_AQ_n
  CCHEM_INPUT_FILE = CCHEM_INPUT_FILE_n
  CCH_TDISCRETIZATION = CCH_TDISCRETIZATION_n
  NCH_SUBSTEPS = NCH_SUBSTEPS_n
  LCH_TUV_ONLINE = LCH_TUV_ONLINE_n
  CCH_TUV_LOOKUP = CCH_TUV_LOOKUP_n
  CCH_TUV_CLOUDS = CCH_TUV_CLOUDS_n
  XCH_TUV_ALBNEW = XCH_TUV_ALBNEW_n
  XCH_TUV_DOBNEW = XCH_TUV_DOBNEW_n
  XCH_TUV_TUPDATE = XCH_TUV_TUPDATE_n
  CCH_VEC_METHOD = CCH_VEC_METHOD_n
  NCH_VEC_LENGTH = NCH_VEC_LENGTH_n
  XCH_TS1D_TSTEP = XCH_TS1D_TSTEP_n
  CCH_TS1D_COMMENT = CCH_TS1D_COMMENT_n
  CCH_TS1D_FILENAME = CCH_TS1D_FILENAME_n
  CSPEC_PRODLOSS = CSPEC_PRODLOSS_n  
  CSPEC_BUDGET = CSPEC_BUDGET_n
END SUBROUTINE INIT_NAM_CH_MNHCn

SUBROUTINE UPDATE_NAM_CH_MNHCn
  LUSECHEM_n = LUSECHEM
  LUSECHAQ_n = LUSECHAQ
  LUSECHIC_n = LUSECHIC
  LCH_INIT_FIELD_n = LCH_INIT_FIELD
  LCH_PH_n = LCH_PH
  LCH_RET_ICE_n = LCH_RET_ICE
  XCH_PHINIT_n = XCH_PHINIT
  XRTMIN_AQ_n = XRTMIN_AQ
  LCH_CONV_SCAV_n = LCH_CONV_SCAV
  LCH_CONV_LINOX_n = LCH_CONV_LINOX
  CCHEM_INPUT_FILE_n = CCHEM_INPUT_FILE
  CCH_TDISCRETIZATION_n = CCH_TDISCRETIZATION
  NCH_SUBSTEPS_n = NCH_SUBSTEPS
  LCH_TUV_ONLINE_n = LCH_TUV_ONLINE
  CCH_TUV_LOOKUP_n = CCH_TUV_LOOKUP
  CCH_TUV_CLOUDS_n = CCH_TUV_CLOUDS
  XCH_TUV_ALBNEW_n = XCH_TUV_ALBNEW
  XCH_TUV_DOBNEW_n = XCH_TUV_DOBNEW
  XCH_TUV_TUPDATE_n = XCH_TUV_TUPDATE
  CCH_VEC_METHOD_n = CCH_VEC_METHOD
  NCH_VEC_LENGTH_n = NCH_VEC_LENGTH
  XCH_TS1D_TSTEP_n = XCH_TS1D_TSTEP
  CCH_TS1D_COMMENT_n = CCH_TS1D_COMMENT
  CCH_TS1D_FILENAME_n = CCH_TS1D_FILENAME
  CSPEC_PRODLOSS_n = CSPEC_PRODLOSS 
  CSPEC_BUDGET_n = CSPEC_BUDGET
END SUBROUTINE UPDATE_NAM_CH_MNHCn

END MODULE MODN_CH_MNHC_n
