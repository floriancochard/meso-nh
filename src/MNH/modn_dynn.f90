!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     #################
      MODULE MODN_DYN_n
!     #################
!
!!****  *MODN_DYN$n* - declaration of namelist NAM_DYNn
!!
!!    PURPOSE
!!    -------
!       The purpose of this  module is to specify the namelist  NAM_DYNn
!     which concerns the dynamic control variables of one nested model.         
!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_DYN$n : contains declaration of dynamic control variables
!!
!!         XTSTEP   : Time step
!!         LSTEADY_DMASS: Switch to specify if the total dry mass Md is STEADY
!!         CPRESOPT : Choice of the pressure solver
!!         NITR     : Number of iterations for the solver pressure       
!!         XRELAX   : relaxation coefficient for the Richardson's method
!!         LRES     : set tolerence for residual
!!         XRES     : tolerence residual
!!
!!    REFERENCE
!!    ---------
!!      Book2 of Meso-NH documentation (module MODD_DYNn)
!!      Asencio N. et al., 1994, "Le projet de modele non-hydrostatique 
!!    commun CNRM-LA, specifications techniques", Note CNRM/GMME, 26, 139p,
!!    (chapters 2 and 3)
!!          
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/06/94                      
!!      Modification 16/11/94 (Lafore+Pinty+Stein) ABS_LAYER and NUM_DIFF input
!!                                                 parameters
!!      Modification  12/12/94  (J.Stein) add the relaxation parameters
!!      Modifications 06/01/95  (Lafore)  For LSTEADY_DMASS             
!!      Modifications 28/07/97  (Masson)  Supress LSTEADY_DMASS
!!      Modification  10/03/98  (J.Stein) add LHORELAX for each variables
!!      Modifications 22/01/01  (Gazen)   add LHORELAX_SVC2R2 _SVCHEM _SVLG
!!      Modifications 10/01/06  (Masson)  add LITRADJ
!!      Modifications 04/05/07  (Lac)  Separation of num diffusion 
!!      Modifications 16/07/10  (Leriche) add LHORELAX_SVCHIC 
!!      Modifications 09/06/11  (Barthe)  add LHORELAX_SVELEC in namelist
!!      Modifications 15/06/11  (Lac)     add LHORELAX for conditional sampling
!!      Modifications 12/02/12  (Pialat/Tulet) add LHORELAX_SVFF for ForeFire scalar variables
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY : JPSVMAX
USE MODD_DYN_n, ONLY : &
         XTSTEP_n => XTSTEP, &
         CPRESOPT_n => CPRESOPT, &
         NITR_n => NITR, &
         LITRADJ_n => LITRADJ, &
         LRES_n  => LRES, &
         XRES_n => XRES, &
         XRELAX_n => XRELAX, &
         LHORELAX_UVWTH_n => LHORELAX_UVWTH, &
         LHORELAX_RV_n => LHORELAX_RV, &
         LHORELAX_RC_n => LHORELAX_RC, &
         LHORELAX_RR_n => LHORELAX_RR, &
         LHORELAX_RI_n => LHORELAX_RI, &
         LHORELAX_RS_n => LHORELAX_RS, &
         LHORELAX_RG_n => LHORELAX_RG, &
         LHORELAX_RH_n => LHORELAX_RH, &
         LHORELAX_TKE_n => LHORELAX_TKE, &
         LHORELAX_SVC2R2_n => LHORELAX_SVC2R2, &
         LHORELAX_SVC1R3_n => LHORELAX_SVC1R3, &
         LHORELAX_SVELEC_n => LHORELAX_SVELEC, &
         LHORELAX_SVCHEM_n => LHORELAX_SVCHEM, &
         LHORELAX_SVCHIC_n => LHORELAX_SVCHIC, &
         LHORELAX_SVLG_n => LHORELAX_SVLG, &
         LHORELAX_SVPP_n => LHORELAX_SVPP, &
#ifdef MNH_FOREFIRE
         LHORELAX_SVFF_n => LHORELAX_SVFF, &
#endif
         LHORELAX_SVCS_n => LHORELAX_SVCS, &
         LHORELAX_SVDST_n => LHORELAX_SVDST, &
         LHORELAX_SVSLT_n => LHORELAX_SVSLT, &
         LHORELAX_SVAER_n => LHORELAX_SVAER, &
         LHORELAX_SVSNW_n => LHORELAX_SVSNW, &
         LHORELAX_SV_n => LHORELAX_SV, &
         LVE_RELAX_n => LVE_RELAX, &
         LVE_RELAX_GRD_n => LVE_RELAX_GRD, &
         NRIMX_n => NRIMX, &
         NRIMY_n => NRIMY, &
         XRIMKMAX_n => XRIMKMAX, &
         XT4DIFU_n => XT4DIFU, &
         XT4DIFTH_n => XT4DIFTH, &
         XT4DIFSV_n => XT4DIFSV
!
IMPLICIT NONE
!
REAL ,SAVE  :: XTSTEP
CHARACTER(LEN=5),SAVE  :: CPRESOPT
INTEGER ,SAVE  :: NITR
LOGICAL ,SAVE  :: LITRADJ
LOGICAL ,SAVE  :: LRES 
REAL ,SAVE     :: XRES
REAL ,SAVE     :: XRELAX
LOGICAL, SAVE  :: LHORELAX_UVWTH
LOGICAL, SAVE  :: LHORELAX_RV
LOGICAL, SAVE  :: LHORELAX_RC
LOGICAL, SAVE  :: LHORELAX_RR
LOGICAL, SAVE  :: LHORELAX_RI
LOGICAL, SAVE  :: LHORELAX_RS
LOGICAL, SAVE  :: LHORELAX_RG
LOGICAL, SAVE  :: LHORELAX_RH
LOGICAL, SAVE  :: LHORELAX_TKE
LOGICAL, SAVE  :: LHORELAX_SVC2R2
LOGICAL, SAVE  :: LHORELAX_SVC1R3
LOGICAL, SAVE  :: LHORELAX_SVELEC
LOGICAL, SAVE  :: LHORELAX_SVCHEM
LOGICAL, SAVE  :: LHORELAX_SVCHIC
LOGICAL, SAVE  :: LHORELAX_SVLG
LOGICAL, SAVE  :: LHORELAX_SVPP
#ifdef MNH_FOREFIRE
LOGICAL, SAVE  :: LHORELAX_SVFF
#endif
LOGICAL, SAVE  :: LHORELAX_SVCS
LOGICAL, SAVE  :: LHORELAX_SVDST
LOGICAL, SAVE  :: LHORELAX_SVSLT
LOGICAL, SAVE  :: LHORELAX_SVAER
LOGICAL, SAVE  :: LHORELAX_SVSNW
LOGICAL, SAVE, DIMENSION(JPSVMAX)  :: LHORELAX_SV
LOGICAL,SAVE  :: LVE_RELAX
LOGICAL,SAVE  :: LVE_RELAX_GRD
INTEGER,SAVE  :: NRIMX
INTEGER,SAVE  :: NRIMY
REAL,SAVE  :: XRIMKMAX
REAL,SAVE  :: XT4DIFU
REAL,SAVE  :: XT4DIFTH
REAL,SAVE  :: XT4DIFSV
!
NAMELIST/NAM_DYNn/XTSTEP,CPRESOPT,NITR,LITRADJ,LRES,XRES,XRELAX,LHORELAX_UVWTH, &
                  LHORELAX_RV,LHORELAX_RC,LHORELAX_RR,LHORELAX_RI, &
                  LHORELAX_RS,LHORELAX_RG,LHORELAX_RH,LHORELAX_TKE,&
                  LHORELAX_SVC2R2,LHORELAX_SVC1R3, LHORELAX_SVELEC,&
                  LHORELAX_SVCHEM, LHORELAX_SVCHIC, LHORELAX_SVLG, LHORELAX_SVDST, &
                  LHORELAX_SVSLT, LHORELAX_SVAER, LHORELAX_SVPP,LHORELAX_SVSNW, &
                  LHORELAX_SVCS, LHORELAX_SV,LVE_RELAX,LVE_RELAX_GRD,&
#ifdef MNH_FOREFIRE
                  LHORELAX_SVFF, &
#endif
                  NRIMX,NRIMY,XRIMKMAX,XT4DIFU, &
                  XT4DIFTH,XT4DIFSV
!
CONTAINS
!
SUBROUTINE INIT_NAM_DYNn
  XTSTEP = XTSTEP_n
  CPRESOPT = CPRESOPT_n
  NITR = NITR_n
  LITRADJ = LITRADJ_n
  LRES = LRES_n 
  XRES = XRES_n
  XRELAX = XRELAX_n
  LHORELAX_UVWTH = LHORELAX_UVWTH_n
  LHORELAX_RV = LHORELAX_RV_n
  LHORELAX_RC = LHORELAX_RC_n
  LHORELAX_RR = LHORELAX_RR_n
  LHORELAX_RI = LHORELAX_RI_n
  LHORELAX_RS = LHORELAX_RS_n
  LHORELAX_RG = LHORELAX_RG_n
  LHORELAX_RH = LHORELAX_RH_n
  LHORELAX_TKE = LHORELAX_TKE_n
  LHORELAX_SVC2R2 = LHORELAX_SVC2R2_n
  LHORELAX_SVC1R3 = LHORELAX_SVC1R3_n
  LHORELAX_SVCHEM = LHORELAX_SVCHEM_n
  LHORELAX_SVCHIC = LHORELAX_SVCHIC_n
  LHORELAX_SVELEC = LHORELAX_SVELEC_n
  LHORELAX_SVLG = LHORELAX_SVLG_n
  LHORELAX_SVPP = LHORELAX_SVPP_n
#ifdef MNH_FOREFIRE
  LHORELAX_SVFF = LHORELAX_SVFF_n
#endif
  LHORELAX_SVCS = LHORELAX_SVCS_n
  LHORELAX_SVDST = LHORELAX_SVDST_n
  LHORELAX_SVSLT = LHORELAX_SVSLT_n
  LHORELAX_SVAER = LHORELAX_SVAER_n
  LHORELAX_SVSNW = LHORELAX_SVSNW_n    
  LHORELAX_SV = LHORELAX_SV_n
  LVE_RELAX = LVE_RELAX_n
  LVE_RELAX_GRD = LVE_RELAX_GRD_n
  NRIMX = NRIMX_n
  NRIMY = NRIMY_n
  XRIMKMAX = XRIMKMAX_n
  XT4DIFU = XT4DIFU_n
  XT4DIFTH = XT4DIFTH_n
  XT4DIFSV = XT4DIFSV_n
END SUBROUTINE INIT_NAM_DYNn

SUBROUTINE UPDATE_NAM_DYNn
  XTSTEP_n = XTSTEP
  CPRESOPT_n = CPRESOPT
  NITR_n = NITR
  LITRADJ_n = LITRADJ
  LRES_n = LRES
  XRES_n = XRES
  XRELAX_n = XRELAX
  LHORELAX_UVWTH_n = LHORELAX_UVWTH
  LHORELAX_RV_n = LHORELAX_RV
  LHORELAX_RC_n = LHORELAX_RC
  LHORELAX_RR_n = LHORELAX_RR
  LHORELAX_RI_n = LHORELAX_RI
  LHORELAX_RS_n = LHORELAX_RS
  LHORELAX_RG_n = LHORELAX_RG
  LHORELAX_RH_n = LHORELAX_RH
  LHORELAX_TKE_n = LHORELAX_TKE
  LHORELAX_SVC2R2_n = LHORELAX_SVC2R2
  LHORELAX_SVC1R3_n = LHORELAX_SVC1R3
  LHORELAX_SVCHEM_n = LHORELAX_SVCHEM
  LHORELAX_SVCHIC_n = LHORELAX_SVCHIC
  LHORELAX_SVELEC_n = LHORELAX_SVELEC
  LHORELAX_SVLG_n = LHORELAX_SVLG
  LHORELAX_SVPP_n = LHORELAX_SVPP
#ifdef MNH_FOREFIRE
  LHORELAX_SVFF_n = LHORELAX_SVFF
#endif
  LHORELAX_SVCS_n = LHORELAX_SVCS
  LHORELAX_SVDST_n = LHORELAX_SVDST
  LHORELAX_SVSLT_n = LHORELAX_SVSLT
  LHORELAX_SVAER_n = LHORELAX_SVAER
  LHORELAX_SVSNW_n = LHORELAX_SVSNW    
  LHORELAX_SV_n = LHORELAX_SV
  LVE_RELAX_n = LVE_RELAX
  LVE_RELAX_GRD_n = LVE_RELAX_GRD
  NRIMX_n = NRIMX
  NRIMY_n = NRIMY
  XRIMKMAX_n = XRIMKMAX
  XT4DIFU_n = XT4DIFU
  XT4DIFTH_n = XT4DIFTH
  XT4DIFSV_n = XT4DIFSV
END SUBROUTINE UPDATE_NAM_DYNn

END MODULE MODN_DYN_n
