!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 adiab 2006/12/12 15:06:20
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_ANTI_DIFF 
!     ##########################
!
INTERFACE
!
      SUBROUTINE ANTI_DIFF      (HLBCX, HLBCY,PTSTEP,PRHODJ,PAS,          &
                                 PRAUCT,PRAVCT,PRAWCT)           
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
REAL,                     INTENT(IN)    :: PTSTEP        !  Time step
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ        ! (Rho) dry * Jacobian
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PAS           ! variable at mass point
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRAUCT! Antidiffusive contravariant
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRAVCT! component
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRAWCT! of momentum
!
END SUBROUTINE ANTI_DIFF 
!
END INTERFACE
!
END MODULE MODI_ANTI_DIFF  
!
!
!
!     ######################################################################
      SUBROUTINE ANTI_DIFF      (HLBCX, HLBCY,PTSTEP,PRHODJ,PAS,          &
                                 PRAUCT,PRAVCT,PRAWCT)           
!     ######################################################################
!
!!****  *ANTI_DIFF* - calculates the antidiffusion velocities of the
!!                    MPDATA scheme 
!!
!!    PURPOSE
!!    -------
!!
!!    The purpose of the routine is to calculate the antidiffusion
!!    velocities of the MPDATA scheme.  
!
!
!!**  METHOD
!!    ------
!!    
!!    MPDATA ia an iterative advection scheme. The number of iterations
!!    taht the scheme is applied is supplied by the user with the
!!    parameter KLITER. 
!!    For every iteration a new set of antidiffusion
!!    velocities are calculated. For KLITER=2, the antidiffusion velocities
!!    are a function of the contravariant velocities. For KLITER>2 the antidiffusion
!!    velocities are a function of the previous antidiffusion velocities.
!!    
!!    EXTERNAL
!!    --------
!!    GET_DIM_EXT_ll : get extended sub-domain sizes
!!    GET_INDICE_ll  : get physical sub-domain bounds
!!    MXM,MYM,MZM Shuman operators 
!!    DXM,DYM,DZM Shuman operators 
!!    LWEAST_ll,LEAST_ll,LNORTH_ll,LSOUTH_ll : position functions
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    None
!!
!!    REFERENCE
!!    ---------
!!    Book1 of documentation (MPDATA scheme) 
!!
!!    AUTHOR
!!    ------
!!    J. Vila-Guerau   *Meteo-France*
!!    J.-P. Lafore     *Meteo-France* 
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    original     01/09/95
!!    J. Stein     01/03/98   include the cyclic case
!!    P. Jabouille 24/09/98   parallelize the code
!------------------------------------------------------------------------------
!
!*      0.   DECLARATIONS
!            ------------
!
USE MODE_ll
!
USE MODD_PARAMETERS 
USE MODI_SHUMAN 
!
IMPLICIT NONE
!
!*      0.1  DECLARATIONS OF ARGUMENTS
!
!
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
REAL,                     INTENT(IN)    :: PTSTEP        !  Time step
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ        ! (Rho) dry * Jacobian
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PAS           ! variable at mass point
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRAUCT! Antidiffusive contravariant
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRAVCT! component
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRAWCT! of momentum
!
!
!*      0.2  DECLARATIONS OF LOCAL VARIABLES
!
INTEGER:: IIB,IJB,IKB                           ! Begining useful area  in x,y,z directions
INTEGER:: IIE,IJE,IKE                           ! End useful area in x,y,z directions
INTEGER:: IIU,IJU,IKU                           ! Array sizes in i,j,k directions
!
REAL   :: ZDBLTST                                    !  Twice the time step
REAL   :: ZEPSILON                                   !  Variable to ensure antidiffusion
                                                     !  velocities equal to zero when the
                                                     !  quantity is 0 
!
REAL, DIMENSION(SIZE(PRAUCT,1),SIZE(PRAUCT,2),SIZE(PRAUCT,3)) :: ZA,ZB,ZC
                                                     ! Auxiliar variables 
                                                     ! antidiffusion velocities 
!
!*       0.3 PROLOGUE
!
CALL GET_DIM_EXT_ll    ('B',IIU,IJU)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKU=SIZE(PRAUCT,3)
IKB=1+JPVEXT
IKE=IKU-JPVEXT
!
!------------------------------------------------------------------------------
!
!*      1.   CALCULATION ANTIDIFFUSIVE VELOCITIES 
!            ------------------------------------
!
  ZDBLTST =PTSTEP
  ZEPSILON=10.E-15
!
!       1.1 Calculation auxiliar variables A,B and C 
!           ----------------------------------------
!
  ZA(:,:,:)=PRAUCT(:,:,:)*DXM(PAS(:,:,:)/PRHODJ(:,:,:))/      &
           (MXM(PAS(:,:,:))+ZEPSILON) 
  ZB(:,:,:)=PRAVCT(:,:,:)*DYM(PAS(:,:,:)/PRHODJ(:,:,:))/      &
           (MYM(PAS(:,:,:))+ZEPSILON) 
  ZC(:,:,:)=PRAWCT(:,:,:)*DZM(1,IKU,1,PAS(:,:,:)/PRHODJ(:,:,:))/      &
           (MZM(1,IKU,1,PAS(:,:,:))+ZEPSILON) 
!
!       1.2 Calculation antidiffusion velocities  
!           ------------------------------------
!
! u-component antidiffusive velocity
!
  PRAUCT(:,:,:)=PTSTEP/2.* (                                         & 
                          ZA*(                                    &
                                     MXM(PRHODJ)*SIGN(1.,PRAUCT)  &
                                    /ZDBLTST-                     &
                                     PRAUCT                       &
                                    )-                            &
                          PRAUCT*MXM(MYF(ZB)+MZF(1,IKU,1,ZC))             &
                        )
!  
! v-component antidiffusive velocity
!
  PRAVCT(:,:,:)=PTSTEP/2.* (                                         & 
                          ZB*(                                    &
                                     MYM(PRHODJ)*SIGN(1.,PRAVCT)  &
                                    /ZDBLTST-                     &
                                     PRAVCT                       &
                                    )-                            &
                          PRAVCT*MYM(MXF(ZA)+MZF(1,IKU,1,ZC))             &
                        )
!  
!  
! w-component antidiffusive velocity
!
  PRAWCT(:,:,:)=PTSTEP/2.* (                                         & 
                          ZC*(                                    &
                                     MZM(1,IKU,1,PRHODJ)*SIGN(1.,PRAWCT)  &
                                    /ZDBLTST-                     &
                                     PRAWCT                       &
                                     )-                           &
                          PRAWCT*MZM(1,IKU,1,MXF(ZA)+MYF(ZB))             &
                        )
! 
!       1.3 Limit of the antidiffusive velocities to satisfy CFL<1
!           ------------------------------------------------------
!
 PRAUCT(:,:,:)=AMIN1(PRAUCT(:,:,:),( PRHODJ(:,:,:)/ZDBLTST))
 PRAUCT(:,:,:)=AMAX1(PRAUCT(:,:,:),(-PRHODJ(:,:,:)/ZDBLTST))
 PRAVCT(:,:,:)=AMIN1(PRAVCT(:,:,:),( PRHODJ(:,:,:)/ZDBLTST))
 PRAVCT(:,:,:)=AMAX1(PRAVCT(:,:,:),(-PRHODJ(:,:,:)/ZDBLTST))
 PRAWCT(:,:,:)=AMIN1(PRAWCT(:,:,:),( PRHODJ(:,:,:)/ZDBLTST))
 PRAWCT(:,:,:)=AMAX1(PRAWCT(:,:,:),(-PRHODJ(:,:,:)/ZDBLTST))
!
!       1.4 Boundary Conditions antidiffusion velocities  
!           --------------------------------------------
!
IF ( HLBCX(1) /= 'CYCL' ) THEN
!
! open boundary conditions : u-component antidiffusive velocity
!
  IF (LWEST_ll( ))  PRAUCT(IIB,:,:)  =0.
  IF (LEAST_ll( ))  PRAUCT(IIE+1,:,:)=0.
!
END IF
!
IF ( HLBCY(1) /= 'CYCL' ) THEN
!
! open boundary conditions :  v-component antidiffusive velocity
!
  IF (LSOUTH_ll( )) PRAVCT(:,IJB,:)  =0.
  IF (LNORTH_ll( )) PRAVCT(:,IJE+1,:)=0.
!
END IF
!
! w-component antidiffusive velocity
!
  PRAWCT(:,:,IKB)  =0.
  PRAWCT(:,:,IKE+1)=0.
!
!------------------------------------------------------------------------------
END SUBROUTINE ANTI_DIFF  
