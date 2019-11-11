!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      %Z% Lib:%F%, Version:%I%, Date:%D%, Last modified:%E%
!-----------------------------------------------------------------
!      	#################
        MODULE  MODN_ELEC
!      	#################
!
!!
!!*** *MODN_ELEC*
!!
!!    PURPOSE
!!    -------
!       Namelist for the electrical scheme
!!
!!**  AUTHOR
!!    ------
!!    C. Barthe      *LACy*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 09/11/09
!!    M. Chong 26/01/10 Option for fair weather field from Helsdon-Farley
!!    M. Chong 26/08/13 Option for "Beard" effect and LLMA storage
!!    J.-P. Pinty 25/10/13 Option for "Latham" effect
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------

USE MODD_ELEC_DESCR
USE MODD_ELEC_PARAM
USE MODD_LMA_SIMULATOR
!
IMPLICIT NONE
!
NAMELIST /NAM_ELEC/ LOCG, LELEC_FIELD, LFLASH_GEOM,            &
                    LFW_HELFA, LCOSMIC_APPROX, LION_ATTACH,    &  
                    CDRIFT, LRELAX2FW_ION,                     &
                    LINDUCTIVE, LSAVE_COORD, LLNOX_EXPLICIT,   & 
                    LSERIES_ELEC, NTSAVE_SERIES, NFLASH_WRITE, &
                    CNI_CHARGING, XQTC,                        &
                    XLIM_NI_IS, XLIM_NI_IG, XLIM_NI_SG,        &
                    CLSOL, NLAPITR_ELEC, XRELAX_ELEC,          &
                    XETRIG, XEBALANCE, XEPROP, XQEXCES, XQNEUT,&
                    XDFRAC_ECLAIR, XDFRAC_L,                   &
                    XWANG_A, XWANG_B,                          &
                    LSEDIM_BEARD, LIAGGS_LATHAM, LLMA
!
END MODULE MODN_ELEC
