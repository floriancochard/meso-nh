!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!########################
MODULE MODI_ION_BOUND4DRIFT
!########################
!
INTERFACE
!
      SUBROUTINE ION_BOUND4DRIFT (HLBCX,HLBCY,PEFIELDU,PEFIELDV,PSVT)

CHARACTER(LEN=4), DIMENSION(2), INTENT(IN)  :: HLBCX,HLBCY  
REAL, DIMENSION(:,:,:,:),       INTENT(INOUT) :: PSVT
REAL, DIMENSION(:,:,:),         INTENT(IN)  :: PEFIELDU,PEFIELDV
!
END SUBROUTINE ION_BOUND4DRIFT
!
END INTERFACE
!
END MODULE MODI_ION_BOUND4DRIFT
!
!
!     ####################################################################
      SUBROUTINE ION_BOUND4DRIFT (HLBCX,HLBCY,PEFIELDU,PEFIELDV,PSVT)
!     ####################################################################
!
!!****  *ION_BOUND4DRIFT* - routine to force the Lateral Boundary Conditions for
!!                 ion variables by the fair weather concentrations
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!      Only for 'OPEN' case  Boundary Condition type
!!   
!!    EXTERNAL 
!!    --------  
!!    GET_INDICE_ll  : get physical sub-domain bounds
!!    LWEAST_ll,LEAST_ll,LNORTH_ll,LSOUTH_ll : position functions
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------  
!!      Module MODD_PARAMETERS : 
!!        JPHEXT ,JPVEXT 
!!
!!      Module MODD_CONF :
!!        CCONF
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!          M. Chong     12/2010
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!	
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!         
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_NSV,  ONLY: NSV_ELECBEG, NSV_ELECEND, XSVMIN
USE MODD_ELEC_n, ONLY : XCION_POS_FW, XCION_NEG_FW
!
USE MODE_ll
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
!
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN)  :: HLBCX,HLBCY  
REAL, DIMENSION(:,:,:,:),       INTENT(INOUT) :: PSVT
REAL, DIMENSION(:,:,:),         INTENT(IN)  :: PEFIELDU,PEFIELDV
!
!
!*       0.2   declarations of local variables
!
INTEGER :: IIB, IIE  ! index of first and last inner mass points along x
INTEGER :: IJB, IJE  ! index of first and last inner mass points along y
!
!
!------------------------------------------------------------------------------
!
!*       1.    COMPUTE DIMENSIONS OF ARRAYS AND OTHER INDICES:
!              ----------------------------------------------
!
! beginning and end indexes of the physical subdomain
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
!-------------------------------------------------------------------------------
!
!*       2.    LBC FILLING IN THE X DIRECTION (LEFT WEST SIDE):   
!              ------------------------------------------------
!
IF (LWEST_ll( ) .AND. HLBCX(1)=='OPEN') THEN
!
  WHERE ( PEFIELDU(IIB,:,:) <= 0. )         !  OUT(IN)FLOW for POS(NEG) IONS
    PSVT(IIB-1,:,:,NSV_ELECBEG) = MAX( 2.*PSVT(IIB,:,:,NSV_ELECBEG) -          &
                             PSVT(IIB+1,:,:,NSV_ELECBEG), XSVMIN(NSV_ELECBEG) )
    PSVT(IIB-1,:,:,NSV_ELECEND) = XCION_NEG_FW(IIB,:,:)  ! Nb/kg
  ELSEWHERE                            !  IN(OUT)FLOW for NEG(POS) IONS
    PSVT(IIB-1,:,:,NSV_ELECBEG) = XCION_POS_FW(IIB,:,:)  ! Nb/kg
    PSVT(IIB-1,:,:,NSV_ELECEND) = MAX( 2.*PSVT(IIB,:,:,NSV_ELECEND) -          &
                             PSVT(IIB+1,:,:,NSV_ELECEND), XSVMIN(NSV_ELECEND) )
  ENDWHERE
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       3.    LBC FILLING IN THE X DIRECTION (RIGHT EAST SIDE): 
!              -------------------------------------------------
!
IF (LEAST_ll( ) .AND. HLBCX(2)=='OPEN') THEN
! 
  WHERE ( PEFIELDU(IIE+1,:,:) >= 0. )         !  OUT(IN)FLOW for POS(NEG) IONS
    PSVT(IIE+1,:,:,NSV_ELECBEG) = MAX( 2.*PSVT(IIE,:,:,NSV_ELECBEG) -          &
                             PSVT(IIE-1,:,:,NSV_ELECBEG), XSVMIN(NSV_ELECBEG) )
    PSVT(IIE+1,:,:,NSV_ELECEND) = XCION_NEG_FW(IIE,:,:)  ! Nb/kg
  ELSEWHERE                              !  IN(OUT)FLOW for NEG(POS) IONS
    PSVT(IIE+1,:,:,NSV_ELECBEG) = XCION_POS_FW(IIE,:,:)  ! Nb/kg
    PSVT(IIE+1,:,:,NSV_ELECEND) = MAX( 2.*PSVT(IIE,:,:,NSV_ELECEND) -          &
                             PSVT(IIE-1,:,:,NSV_ELECEND), XSVMIN(NSV_ELECEND) )
  ENDWHERE
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       4.    LBC FILLING IN THE Y DIRECTION (BOTTOM SOUTH SIDE): 
!              ---------------------------------------------------
!
IF (LSOUTH_ll( ) .AND. HLBCY(1)=='OPEN') THEN
!
  WHERE ( PEFIELDV(:,IJB,:) <= 0. )         !  OUT(IN)FLOW for POS(NEG) IONS
    PSVT(:,IJB-1,:,NSV_ELECBEG) = MAX( 2.*PSVT(:,IJB,:,NSV_ELECBEG) -          &
                             PSVT(:,IJB+1,:,NSV_ELECBEG), XSVMIN(NSV_ELECBEG) )
    PSVT(:,IJB-1,:,NSV_ELECEND) = XCION_NEG_FW(:,IJB,:)  ! Nb/kg
  ELSEWHERE                            !  IN(OUT)FLOW for NEG(POS) IONS
    PSVT(:,IJB-1,:,NSV_ELECBEG) = XCION_POS_FW(:,IJB,:)  ! Nb/kg
    PSVT(:,IJB-1,:,NSV_ELECEND) = MAX( 2.*PSVT(:,IJB,:,NSV_ELECEND) -          &
                             PSVT(:,IJB+1,:,NSV_ELECEND), XSVMIN(NSV_ELECEND) )
  ENDWHERE
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       5.    LBC FILLING IN THE Y DIRECTION (TOP NORTH SIDE): 
!              ------------------------------------------------
!
IF (LNORTH_ll( ) .AND. HLBCY(2)=='OPEN') THEN
! 
  WHERE ( PEFIELDV(:,IJE+1,:) >= 0. )         !  OUT(IN)FLOW for POS(NEG) IONS
    PSVT(:,IJE+1,:,NSV_ELECBEG) = MAX( 2.*PSVT(:,IJE,:,NSV_ELECBEG) -          &
                             PSVT(:,IJE-1,:,NSV_ELECBEG), XSVMIN(NSV_ELECBEG) )
    PSVT(:,IJE+1,:,NSV_ELECEND) = XCION_NEG_FW(:,IJE,:)  ! Nb/kg
  ELSEWHERE                              !  IN(OUT)FLOW for NEG(POS) IONS
    PSVT(:,IJE+1,:,NSV_ELECBEG) = XCION_POS_FW(:,IJE,:)  ! Nb/kg
    PSVT(:,IJE+1,:,NSV_ELECEND) = MAX( 2.*PSVT(:,IJE,:,NSV_ELECEND) -          &
                             PSVT(:,IJE-1,:,NSV_ELECEND), XSVMIN(NSV_ELECEND) )
  ENDWHERE
END IF
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ION_BOUND4DRIFT
