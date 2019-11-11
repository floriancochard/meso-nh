!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/05/18 13:07:25
!-----------------------------------------------------------------
!!   ########################
     MODULE MODI_CH_AER_INIT_SOA
!!   ########################
!!
INTERFACE
!!
SUBROUTINE CH_AER_INIT_SOA(KOUT,KVERB)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: KOUT, KVERB ! stdout output, verbosity level
END SUBROUTINE CH_AER_INIT_SOA
!!
END INTERFACE
!!
END MODULE MODI_CH_AER_INIT_SOA
!!
!!
!!   ############################################################
     SUBROUTINE CH_AER_INIT_SOA(KOUT,KVERB)
!!   ############################################################
!!
!!    PURPOSE
!!    -------
!!    Realise l'equilibre entre les moments via la masse contenue 
!!    dans les aerosols, les diametres moyens et la dispersion.
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    24/24/14 M. Leriche add ReLACS3
!!
!!    EXTERNAL
!!    --------
!!    None
!!
USE MODD_CH_AEROSOL
USE MODD_CH_M9_n,   ONLY : CNAMES, NEQ
USE MODD_CH_MNHC_n, ONLY : CCH_SCHEME
USE MODD_NSV,       ONLY : NSV_CHEM
!!
IMPLICIT NONE
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
INTEGER, INTENT(IN)  :: KOUT, KVERB ! stdout output, verbosity level
!
!*      0.2    declarations local variables
!
INTEGER  :: JN
!-------------------------------------------------------------------------------
!
!
!*       1.     TRANSFER FROM GAS TO AEROSOL MODULE
!               ------------------------------------
!        1.1    initialisation 
!
! Definition of gas chemical scheme
CCH_SCHEME = "NONE"
DO JN=1, SIZE(CNAMES)
  IF (TRIM(CNAMES(JN)) .EQ. "CO")   JP_CH_CO   = JN
  IF (TRIM(CNAMES(JN)) .EQ. "ALKA") CCH_SCHEME = "RELACS"
  IF (TRIM(CNAMES(JN)) .EQ. "HC3")  CCH_SCHEME = "RACM"
  IF (TRIM(CNAMES(JN)) .EQ. "URG1") CCH_SCHEME = "RELACS2"
  IF (TRIM(CNAMES(JN)) .EQ. "GLY")  CCH_SCHEME = "RELACS3"
  IF (TRIM(CNAMES(JN)) .EQ. "UR29") CCH_SCHEME = "CACM"
ENDDO
IF (((TRIM(CORGANIC)=="MPMPO").OR.(TRIM(CORGANIC)=="PUN")).AND.&
    ((CCH_SCHEME == "RELACS2" .OR. CCH_SCHEME == "CACM" &
      .OR. CCH_SCHEME == "RELACS3"))) THEN
    NSOA = 10
ELSE 
    NSOA = 0 ! No SOA formation
END IF
!
IF (ALLOCATED(CAERONAMES)) DEALLOCATE(CAERONAMES)
NM6_AER = 0
IF (LVARSIGI) NM6_AER = 1
IF (LVARSIGJ) NM6_AER = NM6_AER + 1
ALLOCATE(CAERONAMES((NSP+NSOA+NCARB+1)*JPMODE+NM6_AER))
!
! Index gas scheme <=> Index Orilam
JP_CH_SO4I = 1  
JP_CH_SO4J = 2  
JP_CH_NO3I = 3  
JP_CH_NO3J = 4  
JP_CH_NH3I = 5  
JP_CH_NH3J = 6  
JP_CH_H2OI = 7  
JP_CH_H2OJ = 8  
JP_CH_OCI  = 9  
JP_CH_OCJ  = 10  
JP_CH_BCI  = 11  
JP_CH_BCJ  = 12 
JP_CH_DSTI  = 13  
JP_CH_DSTJ  = 14

JP_CH_M0I = (NCARB + NSP + NSOA)*JPMODE +1
JP_CH_M0J = (NCARB + NSP + NSOA)*JPMODE +2
IF (LVARSIGI) JP_CH_M6I = (NCARB + NSP + NSOA)*JPMODE + 3
IF ((LVARSIGI).AND.(LVARSIGJ)) THEN
   JP_CH_M6J = (NCARB + NSP + NSOA)*JPMODE + 4
ELSE IF ((LVARSIGJ).AND. .NOT.(LVARSIGI)) THEN
   JP_CH_M6J = (NCARB + NSP + NSOA)*JPMODE + 3
END IF


    ! Name of aerosol component
CAERONAMES(JP_CH_SO4I) = "SO4I"
CAERONAMES(JP_CH_SO4J) = "SO4J"
CAERONAMES(JP_CH_NH3I) = "NH3I"
CAERONAMES(JP_CH_NH3J) = "NH3J"
CAERONAMES(JP_CH_H2OI) = "H2OI"
CAERONAMES(JP_CH_H2OJ) = "H2OJ"
CAERONAMES(JP_CH_NO3I) = "NO3I"
CAERONAMES(JP_CH_NO3J) = "NO3J"
CAERONAMES(JP_CH_BCI) = "BCI"
CAERONAMES(JP_CH_BCJ) = "BCJ"
CAERONAMES(JP_CH_OCI) = "OCI"
CAERONAMES(JP_CH_OCJ) = "OCJ"
CAERONAMES(JP_CH_DSTI) = "DSTI"
CAERONAMES(JP_CH_DSTJ) = "DSTJ"

IF (NSOA .EQ. 10) THEN
JP_AER_SOA1 = NCARB + NSP + 1
JP_AER_SOA2 = NCARB + NSP + 2
JP_AER_SOA3 = NCARB + NSP + 3
JP_AER_SOA4 = NCARB + NSP + 4
JP_AER_SOA5 = NCARB + NSP + 5
JP_AER_SOA6 = NCARB + NSP + 6
JP_AER_SOA7 = NCARB + NSP + 7
JP_AER_SOA8 = NCARB + NSP + 8
JP_AER_SOA9 = NCARB + NSP + 9
JP_AER_SOA10= NCARB + NSP + 10 

JP_CH_SOA1I = 15 
JP_CH_SOA1J = 16 
JP_CH_SOA2I = 17 
JP_CH_SOA2J = 18 
JP_CH_SOA3I = 19 
JP_CH_SOA3J = 20 
JP_CH_SOA4I = 21 
JP_CH_SOA4J = 22 
JP_CH_SOA5I = 23 
JP_CH_SOA5J = 24 
JP_CH_SOA6I = 25 
JP_CH_SOA6J = 26 
JP_CH_SOA7I = 27 
JP_CH_SOA7J = 28 
JP_CH_SOA8I = 29 
JP_CH_SOA8J = 30 
JP_CH_SOA9I = 31 
JP_CH_SOA9J = 32 
JP_CH_SOA10I = 33 
JP_CH_SOA10J = 34 

CAERONAMES(JP_CH_SOA1I) = "SOA1I"
CAERONAMES(JP_CH_SOA1J) = "SOA1J"
CAERONAMES(JP_CH_SOA2I) = "SOA2I"
CAERONAMES(JP_CH_SOA2J) = "SOA2J"
CAERONAMES(JP_CH_SOA3I) = "SOA3I"
CAERONAMES(JP_CH_SOA3J) = "SOA3J"
CAERONAMES(JP_CH_SOA4I) = "SOA4I"
CAERONAMES(JP_CH_SOA4J) = "SOA4J"
CAERONAMES(JP_CH_SOA5I) = "SOA5I"
CAERONAMES(JP_CH_SOA5J) = "SOA5J"
CAERONAMES(JP_CH_SOA6I) = "SOA6I"
CAERONAMES(JP_CH_SOA6J) = "SOA6J"
CAERONAMES(JP_CH_SOA7I) = "SOA7I"
CAERONAMES(JP_CH_SOA7J) = "SOA7J"
CAERONAMES(JP_CH_SOA8I) = "SOA8I"
CAERONAMES(JP_CH_SOA8J) = "SOA8J"
CAERONAMES(JP_CH_SOA9I) = "SOA9I"
CAERONAMES(JP_CH_SOA9J) = "SOA9J"
CAERONAMES(JP_CH_SOA10I) = "SOA10I"
CAERONAMES(JP_CH_SOA10J) = "SOA10J"
END IF

CAERONAMES(JP_CH_M0I) = "M0I"
CAERONAMES(JP_CH_M0J) = "M0J"
IF (LVARSIGI) CAERONAMES(JP_CH_M6I) = "M6I"
IF (LVARSIGJ) CAERONAMES(JP_CH_M6J) = "M6J"

END SUBROUTINE CH_AER_INIT_SOA
