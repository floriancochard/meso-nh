!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
MODULE MODI_COMPARE_WITH_PGD_DOMAIN
!##################################
INTERFACE
      SUBROUTINE COMPARE_WITH_PGD_DOMAIN(PLAT0,PLON0,PRPK,PBETA,                    &
                                 PLATOR,PLONOR,PXHAT,PYHAT)
!
REAL             , INTENT(IN) :: PLAT0    ! reference latitude
REAL             , INTENT(IN) :: PLON0    ! reference longitude
REAL             ,  INTENT(IN) :: PRPK     ! projecton factor
REAL             , INTENT(IN) :: PBETA    ! angle of rotation of the domain
REAL             , INTENT(IN) :: PLATOR   ! origine latitude
REAL             , INTENT(IN) :: PLONOR   ! origine longitude
REAL,DIMENSION(:), INTENT(IN) :: PXHAT    ! x coordinate of flux points
REAL,DIMENSION(:), INTENT(IN) :: PYHAT    ! y coordinate of flux points
!
END SUBROUTINE COMPARE_WITH_PGD_DOMAIN
END INTERFACE
END MODULE MODI_COMPARE_WITH_PGD_DOMAIN
!     ######spl
      SUBROUTINE COMPARE_WITH_PGD_DOMAIN(PLAT0,PLON0,PRPK,PBETA,  &
                                 PLATOR,PLONOR,PXHAT,PYHAT)
!     #####################################################
!
!!****  *COMPARE_WITH_PGD_DOMAIN* - check the equality of one domain with one in
!!                          MODD_GRID and MODD_PGDGRID
!! 
!!    PURPOSE
!!    -------
!!
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB : verbosity level for output-listing
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         TLUOUT0 : output-listing file
!!      Module MODD_GRID
!!         XBETA   : rotation of the domain
!!         XRPK    : parameter of projection
!!         XLAT0   : latitude reference for the projection
!!         XLON0   : longitude reference for the projection
!!      Module MODD_PGDGRID   : contains grid definition of PGD grid
!!         XPGDLATOR  : latitude of origine point
!!         XPGDLONOR  : longitude of origine point
!!         XPGDXHAT   : position x in the conformal plane
!!         XPGDYHAT   : position y in the conformal plane
!!         
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!
!!    AUTHOR
!!    ------
!!	
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/02/96
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_FM
USE MODE_IO_ll
USE MODE_MSG
!
USE MODD_CONF        ! declaration modules
USE MODD_LUNIT
USE MODD_GRID
USE MODD_PGDGRID
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
REAL             , INTENT(IN) :: PLAT0    ! reference latitude
REAL             , INTENT(IN) :: PLON0    ! reference longitude
REAL             , INTENT(IN) :: PRPK     ! projecton factor
REAL             , INTENT(IN) :: PBETA    ! angle of rotation of the domain
REAL             , INTENT(IN) :: PLATOR   ! origine latitude
REAL             , INTENT(IN) :: PLONOR   ! origine longitude
REAL,DIMENSION(:), INTENT(IN) :: PXHAT    ! x coordinate of flux points
REAL,DIMENSION(:), INTENT(IN) :: PYHAT    ! y coordinate of flux points
!
!*       0.2   Declaration of local variables
!              ------------------------------
REAL               :: ZEPS       ! a little number
INTEGER            :: ILUOUT0    ! logical number for listing file
INTEGER            :: IIU,    IJU
INTEGER            :: IPGDIU, IPGDJU
REAL               :: ZLON0, ZLONOR
REAL, DIMENSION(SIZE(XPGDXHAT)) :: ZDIFFX ! difference between the two x coordinates
REAL, DIMENSION(SIZE(XPGDYHAT)) :: ZDIFFY ! difference between the two y coordinates
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZATIONS
!              ---------------
!
ILUOUT0 = TLUOUT0%NLU
ZEPS=1.E-7
!
!-------------------------------------------------------------------------------
!
!*       2.    CALCULATIONS
!              ------------
!
IPGDIU=SIZE(XPGDXHAT)
IPGDJU=SIZE(XPGDYHAT)
IIU=SIZE(PXHAT)
IJU=SIZE(PYHAT)
!
IF ( (IIU == IPGDIU) .AND. (IJU ==IPGDJU) ) THEN
  ZDIFFX(:)=(XPGDXHAT(:)-0.5*(XPGDXHAT(1)+XPGDXHAT(2)) ) &
           -(PXHAT(:)   -0.5*(  PXHAT(1) + PXHAT(2)  ) )
  ZDIFFY(:)=(XPGDYHAT(:)-0.5*(XPGDYHAT(1)+XPGDYHAT(2)) ) &
           -(PYHAT(:)   -0.5*(  PYHAT(1) + PYHAT(2)  ) )
ELSE
  ZDIFFX(:)=-99.
  ZDIFFY(:)=-99.
END IF
!
!-------------------------------------------------------------------------------
!
!*       3.    TEST ON THE EQUALITY OF THE HORIZONTAL DOMAINS 
!              ----------------------------------------------
!
WRITE(ILUOUT0,*)
WRITE(ILUOUT0,*) '        input file     physiographic data'
WRITE(ILUOUT0,1) 'LAT0  ',PLAT0, ' ',XLAT0
WRITE(ILUOUT0,1) 'LON0  ',PLON0, ' ',XLON0
WRITE(ILUOUT0,1) 'RPK   ',PRPK,  ' ',XRPK
WRITE(ILUOUT0,1) 'BETA  ',PBETA, ' ',XBETA
WRITE(ILUOUT0,1) 'LATOR ',PLATOR,' ',XPGDLATOR
WRITE(ILUOUT0,1) 'LONOR ',PLONOR,' ',XPGDLONOR
WRITE(ILUOUT0,1) 'DX(1) ',(PXHAT(2)-PXHAT(1)),' ',(XPGDXHAT(2)-XPGDXHAT(1))
WRITE(ILUOUT0,1) 'DY(1) ',(PYHAT(2)-PYHAT(1)),' ',(XPGDYHAT(2)-XPGDYHAT(1))
WRITE(ILUOUT0,2) 'IU    ',IIU, ' ',IPGDIU
WRITE(ILUOUT0,2) 'JU    ',IJU, ' ',IPGDJU
WRITE(ILUOUT0,*)
!
1 FORMAT(A6,F16.9,A1,F16.9,A1,F16.9)
2 FORMAT(A6,I6,A11,I6,A11,I6)
!
ZLON0 =PLON0 +NINT((XLON0    -PLON0 )/360.)*360.
ZLONOR=PLONOR+NINT((XPGDLONOR-PLONOR)/360.)*360.
!
IF (     (ABS(PLAT0-XLAT0)>ZEPS*MAX(1.,ABS(XLAT0)))                      &
   .OR.  (ABS(ZLON0-XLON0)>ZEPS*MAX(1.,ABS(XLON0)))                      &  
   .OR.  (ABS(ABS(PRPK)-ABS(XRPK))>ZEPS*MAX(1.,ABS(XRPK)))               &  
   .OR.  (ABS(PBETA-XBETA)>ZEPS*MAX(1.,ABS(XBETA)))                      &  
   .OR.  (ABS(PLATOR-XPGDLATOR)>ZEPS*MAX(1.,ABS(XPGDLATOR)))             &  
   .OR.  (ABS(ZLONOR-XPGDLONOR)>ZEPS*MAX(1.,ABS(XPGDLONOR)))             &   
   .OR.  (ABS(IIU-IPGDIU)>0)                                             & 
   .OR.  (ABS(IJU-IPGDJU)>0)                                             & 
   .OR.  (ANY((ABS(ZDIFFX))>ZEPS*MAX(1.,ABS(XPGDXHAT(2)-XPGDXHAT(1)))))  &
   .OR.  (ANY((ABS(ZDIFFY))>ZEPS*MAX(1.,ABS(XPGDYHAT(2)-XPGDYHAT(1)))))  ) THEN
  WRITE(ILUOUT0,*)
  WRITE(ILUOUT0,*) '  +----------------------------------------------------------+'
  WRITE(ILUOUT0,*) '  | INPUT FILE AND PHYSIOGRAPHIC DATA DOMAINS ARE DIFFERENTS |'
  WRITE(ILUOUT0,*) '  +----------------------------------------------------------+'
  WRITE(ILUOUT0,*)
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','COMPARE_WITH_PGD_DOMAIN','')
ENDIF
!
!-------------------------------------------------------------------------------
!
WRITE(ILUOUT0,*) 'Routine COMPARE_WITH_PGD_DOMAIN completed'
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE COMPARE_WITH_PGD_DOMAIN
