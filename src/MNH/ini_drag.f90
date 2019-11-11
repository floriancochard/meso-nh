!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######################
      MODULE MODI_INI_DRAG
!     ######################
!
INTERFACE 
!
      SUBROUTINE INI_DRAG(OMOUNT,PZS,PHSTART,KSTART,PDRAG)
REAL, DIMENSION(:,:), INTENT(IN) :: PZS ! Orographie
LOGICAL,               INTENT(IN) :: OMOUNT !Logical switch to activate the no-slip condition on a moutain only
!
INTEGER, INTENT(IN) :: KSTART ! Index from which the no-slip condition is applied, when LMOUNT = .FALSE.  
REAL, INTENT(IN) :: PHSTART ! height above which the no-slip condition is applied, when LMOUNT = .TRUE.  
REAL, DIMENSION(:,:),          INTENT(OUT) :: PDRAG    ! Array defining where the no-slip condition is applied (-1/1)
END SUBROUTINE INI_DRAG
!
END INTERFACE
!
END MODULE MODI_INI_DRAG
!
!
!     ############################################################
        SUBROUTINE INI_DRAG(OMOUNT,PZS,PHSTART,KSTART,PDRAG)
!     ############################################################
!
!!****  *INI_DRAG* - routine to initialize the XDRAG array defining where the no-slip condition is applied
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to set the XDRAG array to -1 where the no-slip is applied and to 1 elsewhere. This values are then used in MODI_VISCi to define the no-slip condition
!!
!!**  METHOD
!!    ------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_CONF   : NVERB
!!
!!
!!    AUTHOR
!!    ------
!!  	J. Colin         * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        27/03/12
!!
!! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
USE MODE_FM
USE MODE_ll
USE MODD_PARAMETERS
!
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN) :: PZS ! Orography 
LOGICAL,               INTENT(IN) :: OMOUNT !Logical switch to activate the no-slip condition on a moutain only
!
INTEGER, INTENT(IN) :: KSTART ! Index from which the no-slip condition is applied, when LMOUNT = .FALSE.  
REAL, INTENT(IN) :: PHSTART ! height above which the no-slip condition is applied, when LMOUNT = .TRUE.
!
REAL, DIMENSION(:,:),          INTENT(OUT) :: PDRAG    ! Array defining where the no-slip condition is applied (-1/1)
! where PDRAG = 1 => Free slip
! where PDRAG = -1 => No slip
!
!*       0.2   declarations of local variables
!
INTEGER             :: IIU     ! Upper dimension in x direction (local)
INTEGER             :: IJU     ! Upper dimension in y direction (local)
INTEGER             :: IIU_ll  ! Upper dimension in x direction (global)
INTEGER             :: IJU_ll  ! Upper dimension in y direction (global)
INTEGER             :: IXORI     ! 
INTEGER             :: IXENDI     ! 
INTEGER             :: IYORI    ! 
INTEGER             :: IYENDI    ! indice intersection
INTEGER             :: IINFO    ! 
INTEGER             :: JI,JJ   ! Loop index 
!
!-------------------------------------------------------------------------------
!
PDRAG = 1.
! 
CALL GET_DIM_EXT_ll('B',IIU,IJU)
!
 IF (OMOUNT) THEN
      DO JI=1,IIU
         DO JJ=1,IJU
             IF ( PZS(JI,JJ) >= PHSTART ) THEN
          PDRAG(JI,JJ) = -1.
             ENDIF
          ENDDO
      ENDDO
    ELSE
      CALL GET_GLOBALDIMS_ll(IIU_ll,IJU_ll)
      IIU_ll = IIU_ll + 2*JPHEXT
      IJU_ll = IJU_ll + 2*JPHEXT
      CALL GET_INTERSECTION_ll(KSTART,1,IIU_ll,IJU_ll, &
            IXORI, IYORI, IXENDI, IYENDI, "EXTE",IINFO)   
        IF ((IINFO /= 1).AND.(IINFO/=-1)) THEN
        PDRAG(IXORI:IXENDI,IYORI:IYENDI) = -1
        ENDIF
    ENDIF
!
END SUBROUTINE INI_DRAG
