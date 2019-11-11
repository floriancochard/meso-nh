!MNH_LIC Copyright 2000-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_RESET_EXSEG
!     #######################
!
INTERFACE
!
      SUBROUTINE RESET_EXSEG(HLUOUT)
!
CHARACTER (LEN=*),  INTENT(IN) :: HLUOUT ! Name  for output listing
!
END SUBROUTINE RESET_EXSEG
!
END INTERFACE
!
END MODULE MODI_RESET_EXSEG
!
!     ##############################
      SUBROUTINE RESET_EXSEG(HLUOUT)
!     ##############################
!
!!****  *RESET_EXSEG* - routine used to mofify the EXSEG1.nam informations
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to modify the informations read in
!     the DESFM file and corrected according to the EXSEG1.nam file.
!     For the DIAG program, we use the informations read in the DIAG1.nam
!     file to correct or reset these informations before the allocations
!     which are performed in ini_modeln
!
!!**  METHOD
!!    ------
!!
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
!!    AUTHOR
!!    ------
!!      J. Stein       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    14/09/00
!!      Modifications  04/06/02  (P Jabouille) reset radiation and convective options
!!                   02/2018 Q.Libois ECRAD
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODE_FM, ONLY : IO_FILE_OPEN_ll,IO_FILE_CLOSE_ll
USE MODE_FMREAD
USE MODE_IO_ll
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_FIND_BYNAME
USE MODE_POS
!
USE MODD_DIAG_FLAG
USE MODD_CH_MNHC_n, ONLY: LUSECHEM
USE MODD_CONF_n, ONLY: LUSERV
USE MODD_GET_n
USE MODD_IO_ll,   ONLY: TFILEDATA
USE MODD_PARAM_n, ONLY: CDCONV, CRAD
USE MODN_PARAM_KAFR_n
USE MODN_PARAM_RAD_n
USE MODN_PARAM_ECRAD_n
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
CHARACTER (LEN=*),  INTENT(IN) :: HLUOUT ! Name for output listing
!
!*       0.2   declarations of local variables
!
INTEGER :: IRESP,ILUNAM        ! return code and logical unit number
LOGICAL :: GFOUND              ! Return code when searching namelist
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
INTEGER           :: IGRID          ! IGRID : grid indicator
INTEGER           :: ILENCH         ! ILENCH : length of comment string
TYPE(TFILEDATA),POINTER :: TZNMLFILE! Namelist file
!
!-------------------------------------------------------------------------------
!
!*       1.    OPENING NAMELIST FILE
!              ---------------------
!
TZNMLFILE  => NULL()
!
CALL IO_FILE_FIND_BYNAME('DIAG1.nam',TZNMLFILE,IRESP)
CALL IO_FILE_OPEN_ll(TZNMLFILE)
ILUNAM = TZNMLFILE%NLU
!
!-------------------------------------------------------------------------------
!
!*       2.    CONVECTION INITIALIZATION CORRECTION
!              ------------------------------------
!
! if we ask to compute the convection diagnostics then the fields used in
! Meso-NH to store these diagnostics must be allocated by the ini_model1 subroutine
!
IF (NCONV_KF>=0) THEN
  CALL POSNAM(ILUNAM,'NAM_PARAM_KAFRN',GFOUND)
  IF (GFOUND) THEN
    CALL INIT_NAM_PARAM_KAFRn
    READ(UNIT=ILUNAM,NML=NAM_PARAM_KAFRN)
    PRINT*, '  namelist NAM_PARAM_KAFRN read'
  END IF
  IF (LUSERV) THEN
    LDIAGCONV=.TRUE.
  ELSE
    NCONV_KF=-1
  END IF
  CALL UPDATE_NAM_PARAM_KAFRn  !because of LDIAGCONV
  IF (CDCONV=='NONE' ) THEN
    CDCONV='KAFR'
    CGETCONV='INIT'
  END IF
END IF
!
PRINT*,'RESET_EXSEG OUTPUT: NCONV_KF=',NCONV_KF,' CDCONV=',CDCONV,' CGETCONV=',CGETCONV
!
!-------------------------------------------------------------------------------
!
!*       3.    RADIATION INITIALIZATION CORRECTION
!              -----------------------------------
!
CGETRAD='READ'
!
IF (CRAD=='NONE') THEN
  CRAD='ECMW'
  LCLEAR_SKY=.FALSE.
  CGETRAD='INIT'
END IF
!
IF(NRAD_3D>=1) THEN
  CALL POSNAM(ILUNAM,'NAM_PARAM_RADN',GFOUND)
  IF (GFOUND) THEN
    CALL INIT_NAM_PARAM_RADn
    READ(UNIT=ILUNAM,NML=NAM_PARAM_RADN)
    CALL UPDATE_NAM_PARAM_RADn
    PRINT*, '  namelist NAM_PARAM_RADN read'
  END IF
#ifdef MNH_ECRAD
  CALL POSNAM(ILUNAM,'NAM_PARAM_ECRADN',GFOUND)
  IF (GFOUND) THEN
    CALL INIT_NAM_PARAM_EcRADn
    READ(UNIT=ILUNAM,NML=NAM_PARAM_ECRADN)
    CALL UPDATE_NAM_PARAM_ECRADn
    PRINT*, '  namelist NAM_PARAM_ECRADN read'    
  END IF
#endif
ENDIF
!

IF ( NRAD_3D>=1 ) THEN
  CRAD='ECMW'
  CGETRAD='INIT'
END IF
!
IF(LEN_TRIM(CRAD_SAT) /= 0) THEN
  CRAD='ECMW'
END IF
!
PRINT*,'RESET_EXSEG OUTPUT: NRAD_3D =',NRAD_3D,' CRAD =',CRAD,' CGETRAD =',CGETRAD
!
!-------------------------------------------------------------------------------
!
!*       4.    CHEMISTRY INITIALIZATION CORRECTION
!              -----------------------------------
!
IF (LUSECHEM .AND. .NOT.LCHEMDIAG) LUSECHEM =.FALSE. 
!
PRINT*,'RESET_EXSEG OUTPUT: LUSECHEM =',LUSECHEM,' LCHEMDIAG =',LCHEMDIAG
PRINT*,' '
!
!-------------------------------------------------------------------------------
!
CALL IO_FILE_CLOSE_ll(TZNMLFILE)
!
END SUBROUTINE RESET_EXSEG
