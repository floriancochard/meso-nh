!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/cvsroot/MNH-VX-Y-Z/src/MNH/modn_advn.f90,v $ $Revision: 1.2.4.1.18.2 $
! MASDEV4_7 modn 2007/03/27 09:58:10
!-----------------------------------------------------------------
!     #################
      MODULE MODN_ADV_n
!     #################
!
!!****  *MODN_ADV$n* - declaration of namelist NAM_ADVn
!!
!!    PURPOSE
!!    -------     
!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_ADV$n : contains declaration of scalar advection schemes
!!      parameters
!!
!!    REFERENCE
!!    ---------
!!          
!!    AUTHOR
!!    ------
!!	Vila, Lafore   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     23/10/95  (Vila, lafore) For new scalar advection schemes
!!      C.Lac        24/04/06  Introduction of CUVW_ADV_SCHEME and
!!                             removal of CFV_ADV_SCHEME
!!      J.-P. Pinty  20/03/10  Add NWENO_ORDER
!!      C.Lac, V.Masson        Add CTEMP_SCHEME and time splitting
!!                  C.LAC 10/2016 : Add OSPLIT_WENO
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_ADV_n, ONLY: &
         CUVW_ADV_SCHEME_n => CUVW_ADV_SCHEME, &
         CMET_ADV_SCHEME_n => CMET_ADV_SCHEME, &
         CSV_ADV_SCHEME_n => CSV_ADV_SCHEME, &
         CTEMP_SCHEME_n => CTEMP_SCHEME, &
         NWENO_ORDER_n => NWENO_ORDER, &
         LSPLIT_CFL_n => LSPLIT_CFL, &
         LSPLIT_WENO_n => LSPLIT_WENO, &
         LCFL_WRIT_n => LCFL_WRIT, &
         XSPLIT_CFL_n => XSPLIT_CFL
!
IMPLICIT NONE
!
CHARACTER(LEN=6)  :: CUVW_ADV_SCHEME
CHARACTER(LEN=6)  :: CMET_ADV_SCHEME
CHARACTER(LEN=6)  :: CSV_ADV_SCHEME
CHARACTER(LEN=4)  :: CTEMP_SCHEME
INTEGER  :: NWENO_ORDER
LOGICAL  :: LSPLIT_CFL     
LOGICAL  :: LSPLIT_WENO    
LOGICAL  :: LCFL_WRIT     
REAL     :: XSPLIT_CFL     
!
NAMELIST/NAM_ADVn/CUVW_ADV_SCHEME,CMET_ADV_SCHEME,CSV_ADV_SCHEME,CTEMP_SCHEME, &
                  NWENO_ORDER,LSPLIT_CFL,LSPLIT_WENO,XSPLIT_CFL,LCFL_WRIT             
!
CONTAINS
!
SUBROUTINE INIT_NAM_ADVn
  CUVW_ADV_SCHEME = CUVW_ADV_SCHEME_n
  CMET_ADV_SCHEME = CMET_ADV_SCHEME_n
  CSV_ADV_SCHEME = CSV_ADV_SCHEME_n
  CTEMP_SCHEME = CTEMP_SCHEME_n
  NWENO_ORDER = NWENO_ORDER_n
  LSPLIT_CFL = LSPLIT_CFL_n
  LSPLIT_WENO = LSPLIT_WENO_n
  LCFL_WRIT = LCFL_WRIT_n
  XSPLIT_CFL = XSPLIT_CFL_n
END SUBROUTINE INIT_NAM_ADVn

SUBROUTINE UPDATE_NAM_ADVn
  CUVW_ADV_SCHEME_n = CUVW_ADV_SCHEME
  CMET_ADV_SCHEME_n = CMET_ADV_SCHEME
  CSV_ADV_SCHEME_n = CSV_ADV_SCHEME
  CTEMP_SCHEME_n = CTEMP_SCHEME
  NWENO_ORDER_n = NWENO_ORDER
  LSPLIT_CFL_n = LSPLIT_CFL      
  LSPLIT_WENO_n = LSPLIT_WENO     
  LCFL_WRIT_n = LCFL_WRIT      
  XSPLIT_CFL_n = XSPLIT_CFL      
END SUBROUTINE UPDATE_NAM_ADVn

END MODULE MODN_ADV_n
