!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
!####################
MODULE MODI_RAD_BOUND
!####################
!
INTERFACE
!
      SUBROUTINE RAD_BOUND (HLBCX,HLBCY,HTURB, PCARPKMAX,             &
                        PTSTEP,PDXHAT,PDYHAT,PZHAT,                   &
                        PUT,PVT,                                      &
                        PLBXUM,PLBYVM,PLBXUS,PLBYVS,                  &
                        PCPHASE,PCPHASE_PBL,PRHODJ,                   &
                        PTKET,PRUS,PRVS,PRWS                          )
! 
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX,HLBCY   ! X and Y-direc. LBC type
CHARACTER(LEN=4),               INTENT(IN) :: HTURB         ! Turbulence scheme
!
!
REAL,                     INTENT(INOUT) :: PCARPKMAX    ! Rayleigh damping amplitude
REAL,                     INTENT(IN) :: PTSTEP      ! time step dt 
REAL,      DIMENSION(:),  INTENT(IN) :: PDXHAT      ! X-direc. meshlength 
REAL,      DIMENSION(:),  INTENT(IN) :: PDYHAT      ! Y-direc. meshlength
REAL,      DIMENSION(:),  INTENT(IN) :: PZHAT       ! height level without orography
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PUT,PVT     ! at t
!
! Lateral Boundary fields at time t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PLBXUM,PLBYVM    
! temporal derivative of the Lateral Boundary fields 
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PLBXUS,PLBYVS
!
REAL,                     INTENT(IN) :: PCPHASE     ! prescribed phase velocity
REAL,                     INTENT(IN) :: PCPHASE_PBL ! prescribed PBL phase velocity
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ      ! Jacobian * dry density 
                                                    ! of the reference state 
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PTKET ! TKE at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT):: PRUS,PRVS   ! Horizontal and Vertical
REAL, DIMENSION(:,:,:),   INTENT(INOUT):: PRWS        ! momentum tendencies
!
END SUBROUTINE RAD_BOUND
!
END INTERFACE
!
END MODULE MODI_RAD_BOUND
!
!     #################################################################
      SUBROUTINE RAD_BOUND (HLBCX,HLBCY,HTURB, PCARPKMAX,             &
                        PTSTEP,PDXHAT,PDYHAT,PZHAT,                   &
                        PUT,PVT,                                      &
                        PLBXUM,PLBYVM,PLBXUS,PLBYVS,                  &
                        PCPHASE,PCPHASE_PBL,PRHODJ,                   &
                        PTKET,PRUS,PRVS,PRWS                          )
!     #################################################################
!
!!****  *RAD_BOUND* - routine computing the velocity components normal to
!!                 Lateral boundaries, at the boundary localization 
!!                 for time t+dt.
!!
!!    PURPOSE
!!    -------
!       This computation of wind normal component to the Lateral Boundaries,
!       is required to deduce the pressure gradient at the boundaries, used
!       as boundary conditions to the elliptic pressure equation to
!       be solved in the following step. 
!
!!**  METHOD
!!    ------
!!       The computation of the adequate Normal velocity depends on the Lateral
!!      Boundary Condition type (variables CLBCX and CLBCY).
!!      
!!       For the 'OPEN' case, a radiation Sommerfeld equation's type 
!!      is used:
!!               dUn     dUn |         [ dUn       dUn |     ] 
!!              ---- =  -----|     - C [-----  -  -----|     ]
!!               dt      dt  | LS      [ dx         dx | LS  ]
!!               
!!               
!!        to account for Large Scale (LS) forcings, as proposed by 
!!       Carpenter (1982), where C is a "magic" phase velocity.
!!      
!!         A semi-implicit scheme is adopted as in Orlanski (1976).
!!      
!!       For the 'DAVI' case, the normal velocity is taken as the one of the
!!      LS field.
!!                 Un  =  Un |
!!                           | LS
!!   
!!    EXTERNAL
!!    --------  
!!    GET_INDICE_ll  : get physical sub-domain bounds
!!    LWEAST_ll,LEAST_ll,LNORTH_ll,LSOUTH_ll : position functions
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!      Module MODD_CONF   : contains configuration variable 
!!           CCONF :  Configuration of models
!!
!!      Module MODD_PARAMETERS : contains parameters commun to all models
!!        JPHEXT : Horizontal EXTernal points number (JPHEXT=1 for this version)
!!        JPVEXT : Vertical   EXTernal points number (JPVEXT=1 for this version)
!!        
!!    REFERENCE
!!    ---------
!!      Book1 and book2 of documentation (routine RAD_BOUND)
!!      Orlanski (1976), Carpenter (1982)
!!
!!    AUTHOR
!!    ------
!!	J.-P. Lafore J. Stein     * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    17/10/94  
!!      Modification 13/02/95  (Lafore)  to account for the OPEN case and
!!                                         for the LS fields introduction
!!      Modification 16/02/95  (Mallet)  bug in computating source terms
!!                                     (depends on START or RESTA configuration)
!!      Modification 06/03/95    "       change the cphase computation
!!                                     (depends on outflow or inflow)    
!!                   16/03/95  (J.Stein) remove R from the historical variables
!!      Modification 31/05/95  (Lafore)  MASTER_DEV2.1 preparation after the
!!                                       LBC tests performed by I. Mallet
!!      Modification 03/09/96  (Lafore)  computing of LS forcing at central time t
!!                                       previously it was done at t_dt
!!                                       removing R from the LS sources
!!      Modification 22/11/96  (Lafore)  tests based on LSTEADYLS replaced by the
!!                                       size of LS sources
!!      Modification 12/12/96  (Lafore)  add a relaxation toward LS fields in the
!!                                       Carpenter equation (semi-implicit treatment)
!!      Modification 06/01/97  (Masson)  bug at J=2 frontier for relaxation in
!!                                       Carpenter equation
!!      Modification 20/10/97  (Lafore, Stein) introduction of the DAVI type of lbc 
!!      Modification 12/11/97  ( Stein)  introduction of LB FIELDS
!!      Modification 05/05/04  (Escobar) INOUT argument for PRUS,PRVS,PRWS
!!      Modification   05/06   (C.Lac)   Remove the DAVI type of lbc
!!      Modification   11/09   (C.Lac)   Impose a zero prescribed phase velocity
!!                                       in the PBL
!!      Juan 25/02/2010: BUG add ZTKEX = 0.0 
!!      Modification   08/10  (V.Masson) Bug correction and add cphase_profile
!!      Escobar     9/11/2010 : cphas_profile : array bound problem if NO Turb =>  PTKET optional
!!      Lac.C.       2011     : Adaptation to FIT temporal scheme
!!      Modification 06/13     (C.Lac)   Introduction of cphase_pbl
!!      Modification 03/14     (C.Lac)   Replacement of XRIMKMAX by XCARPKMAX 
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1                             
!!      
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF         
USE MODD_PARAMETERS
USE MODD_CTURB
!
USE MODI_CPHASE_PROFILE
!
USE MODE_ll
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
! 
! 
! 
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX,HLBCY   ! X and Y-direc. LBC type
CHARACTER(LEN=4),               INTENT(IN) :: HTURB         ! Turbulence scheme
!
!
REAL,                     INTENT(INOUT) :: PCARPKMAX    ! Rayleigh damping amplitude
REAL,                     INTENT(IN) :: PTSTEP      ! time step dt 
REAL,      DIMENSION(:),  INTENT(IN) :: PDXHAT      ! X-direc. meshlength 
REAL,      DIMENSION(:),  INTENT(IN) :: PDYHAT      ! Y-direc. meshlength
REAL,      DIMENSION(:),  INTENT(IN) :: PZHAT       ! height level without orography
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PUT,PVT     ! at t
!
! Lateral Boundary fields at time t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PLBXUM,PLBYVM    
! temporal derivative of the Lateral Boundary fields 
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PLBXUS,PLBYVS
!
REAL,                     INTENT(IN) :: PCPHASE     ! prescribed phase velocity
REAL,                     INTENT(IN) :: PCPHASE_PBL ! prescribed PBL phase velocity
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ      ! Jacobian * dry density 
                                                    ! of the reference state 
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTKET ! TKE at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT):: PRUS,PRVS   ! Horizontal and Vertical
REAL, DIMENSION(:,:,:),   INTENT(INOUT):: PRWS        ! momentum tendancies
!
!*       0.2   declarations of local variables
!
INTEGER             :: IIB       ! indice I Beginning in x direction
INTEGER             :: IJB       ! indice J Beginning in y direction
INTEGER             :: IIE       ! indice I End       in x direction 
INTEGER             :: IJE       ! indice J End       in y direction 
INTEGER             :: IKB       ! indice K Beginning in z direction 
INTEGER             :: IKE       ! indice K End       in z direction 
INTEGER             :: ILBX,ILBY ! number of points of the RIM arrays
! 
REAL                :: ZINVTSTEP ! Inverse of the applicable timestep
REAL                :: ZKTSTEP   !  Rayleigh damping by the timestep
!
REAL, DIMENSION(SIZE(PUT,2),SIZE(PUT,3)) :: ZLBGU  ! LB x-gradient 
REAL, DIMENSION(SIZE(PUT,2),SIZE(PUT,3)) :: ZLBEU  ! LB temporal Evolution
REAL, DIMENSION(SIZE(PUT,2),SIZE(PUT,3)) :: ZLBXU  ! LB field at t or t+1        
REAL, DIMENSION(SIZE(PUT,2),SIZE(PUT,3)) :: ZCPHASX! Normalized Phase velocity
!                                                 ! for U field at X-boundaries
REAL, DIMENSION(SIZE(PUT,2),SIZE(PUT,3)) :: ZPHASX!  Phase velocity
!                                                 ! for U field at X-boundaries
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,3)) :: ZLBGV  ! LB y-gradient 
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,3)) :: ZLBEV  ! LB temporal Evolution
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,3)) :: ZLBYV  ! LB field at t or t+1         
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,3)) :: ZCPHASY! Normalized Phase velocity
!                                                  ! for V field at Y-boundaries
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,3)) :: ZPHASY ! Phase velocity
!                                                  ! for V field at Y-boundaries
REAL                                     :: ZALPHA2! implicitness of the damping
!
!-------------------------------------------------------------------------------
!
!*       1.    PROLOGUE
!              --------
!
!*       1.1  Compute dimensions of arrays and other indices
! 
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PUT,3) - JPVEXT
!
!
!*       1.2  Compute the inverse of the applicable timestep
!
!
ZINVTSTEP = 1./PTSTEP
IF (PCARPKMAX == XUNDEF) PCARPKMAX = 1./ (10.*PTSTEP)
ZKTSTEP   = PCARPKMAX*PTSTEP
! ZALPHA2 = O : explicit ; ZALPHA2 = 1 : implicit ; ZALPHA2 = 0.5 SI
ZALPHA2 = 1. 
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       2.    LBC FILLING IN THE X DIRECTION (LEFT WEST SIDE):
!              ------------------------------
!       ====>  It only concernes U component 
!                                ----------- 
!
IF (LWEST_ll( )) THEN
! 
SELECT CASE ( HLBCX(1) )      
!
!*       2.1  WALL CASE:
!             ==========  
!
  CASE ('WALL')
!
       PRUS (IIB,:,:) = 0.
!
!*       2.2  OPEN CASE:
!             ========= 
!
  CASE ('OPEN')
!
    IF (HTURB /= "NONE" ) THEN
       CALL CPHASE_PROFILE(PZHAT,PCPHASE,PCPHASE_PBL,ZPHASX,PTKET(IIB,:,:))
    ELSE
       CALL CPHASE_PROFILE(PZHAT,PCPHASE,PCPHASE_PBL,ZPHASX)
    END IF

    ZCPHASX(:,:) = MAX ( 0., MIN ( 1.,                                      & 
                      (-PUT(IIB,:,:) + ZPHASX(:,:) ) * PTSTEP / PDXHAT(IIB)  )    )
                      ! notice that ZCPHASX=0. when ZPHASX  < PUT(IIB,:,:)  
!
!
    IF ( SIZE(PLBXUS,1) == 0 ) THEN
      ZLBEU (:,:) = 0.
      ZLBGU (:,:) = PLBXUM(JPHEXT+1,:,:) - PLBXUM(JPHEXT,:,:)  ! 2 - 1
      ZLBXU(:,:)  = PLBXUM(JPHEXT,:,:) ! 1
    ELSE
      ZLBEU (:,:) = PLBXUS(JPHEXT,:,:) ! 1
      ZLBGU (:,:) = PLBXUM(JPHEXT+1,:,:) - PLBXUM(JPHEXT,:,:) +  & ! 2 -  1
                      PTSTEP * (PLBXUS(JPHEXT+1,:,:) - PLBXUS(JPHEXT,:,:)) ! 2 - 1
      ZLBXU(:,:)  = PLBXUM(JPHEXT,:,:) + PTSTEP * PLBXUS(JPHEXT,:,:) ! 1  + 1
    END IF
!  
!     ============================================================
!
!     PRUS (IIB,:,:) =(PRHODJ(IIB-1,:,:) + PRHODJ(IIB,:,:)) * 0.5 *           &
!                    (  (1. - ZCPHASX(:,:) - ZKTSTEP) * PUM(IIB  ,:,:)      &
!                      + 2. * ZCPHASX(:,:)            * PUT(IIB+1,:,:)      &
!           +2.* (   ZLBEU (:,:) * ZTSTEP                                   &
!                 -  ZLBGU (:,:) * ZCPHASX(:,:)                             &
!                 +  ZKTSTEP*( ZLBXU(:,:) )       )                         &
!                    ) * ZINVTSTEP / (1.+ ZCPHASX(:,:) +ZKTSTEP)
!
      PRUS (IIB,:,:) =(PRHODJ(IIB-1,:,:) + PRHODJ(IIB,:,:)) * 0.5 *          &
                       ZINVTSTEP / (1.+ ZKTSTEP * ZALPHA2 )  *               &
            (  (1. - ZCPHASX(:,:) - ZKTSTEP * (1. - ZALPHA2)) * PUT(IIB,:,:) &
                         +  ZCPHASX(:,:)            * PUT(IIB+1  ,:,:)       &
              + (   ZLBEU (:,:) * PTSTEP                                     &
                    -  ZLBGU (:,:) * ZCPHASX(:,:)                            &
                    +  ZKTSTEP*ZLBXU(:,:)   )  )  

!
!
END SELECT
!
END IF
!-------------------------------------------------------------------------------
!
!*       3.    LBC FILLING IN THE X DIRECTION (RIGHT EAST SIDE):
!              ------------------------------
!       ====>  It only concernes U component 
!                                ----------- 
!
IF (LEAST_ll( )) THEN
!
!
SELECT CASE ( HLBCX(2) )     
!
!*       3.1  WALL CASE: 
!              ========= 
!
  CASE ('WALL')
!
       PRUS (IIE+1,:,:) = 0.
!
!*       3.2  OPEN CASE:
!             =========
!
  CASE ('OPEN')
!
    IF (HTURB /= "NONE" ) THEN
       CALL CPHASE_PROFILE(PZHAT,PCPHASE,PCPHASE_PBL,ZPHASX,PTKET(IIE,:,:))
    ELSE
       CALL CPHASE_PROFILE(PZHAT,PCPHASE,PCPHASE_PBL,ZPHASX)
    END IF
!
    ZCPHASX(:,:) = MAX ( 0., MIN ( 1.,                                  &
                    ( PUT(IIE+1,:,:) + ZPHASX(:,:) ) * PTSTEP/PDXHAT(IIE)  )      )
!
! 
    ILBX=SIZE(PLBXUM,1)
    IF (SIZE(PLBXUS,1) == 0 ) THEN
      ZLBEU (:,:) = 0.
      ZLBGU (:,:) = PLBXUM(ILBX-JPHEXT+1,:,:) - PLBXUM(ILBX-JPHEXT,:,:) ! ILBX / (ILBX-1
      ZLBXU(:,:)  = PLBXUM(ILBX-JPHEXT+1,:,:)
    ELSE
      ZLBEU (:,:) = PLBXUS(ILBX-JPHEXT+1,:,:)
      ZLBGU (:,:) = PLBXUM(ILBX-JPHEXT+1,:,:) - PLBXUM(ILBX-JPHEXT,:,:) +  &
                      PTSTEP * (PLBXUS(ILBX-JPHEXT+1,:,:) - PLBXUS(ILBX-JPHEXT,:,:))
      ZLBXU(:,:)  = PLBXUM(ILBX-JPHEXT+1,:,:) + PTSTEP * PLBXUS(ILBX-JPHEXT+1,:,:)
    END IF
!     
!     ============================================================
! 
!     PRUS (IIE+1,:,:) =(PRHODJ(IIE+1,:,:) + PRHODJ(IIE,:,:)) * 0.5 *         &
!                      (  (1. - ZCPHASX(:,:) - ZKTSTEP) * PUM(IIE+1,:,:)      &
!                        + 2. * ZCPHASX(:,:)            * PUT(IIE  ,:,:)      &
!             +2.* (   ZLBEU (:,:) * ZTSTEP                                   &
!                   +  ZLBGU (:,:) * ZCPHASX(:,:)                             &
!                   +  ZKTSTEP*ZLBXU(:,:)   )                                 &
!                      ) * ZINVTSTEP / (1.+ZCPHASX(:,:) +ZKTSTEP)
! 
      PRUS (IIE+1,:,:) =(PRHODJ(IIE+1,:,:) + PRHODJ(IIE,:,:)) * 0.5 *           &
                       ZINVTSTEP / (1.+ ZKTSTEP * ZALPHA2 )  *                  &
            (  (1. - ZCPHASX(:,:) - ZKTSTEP * (1. - ZALPHA2) ) * PUT(IIE+1,:,:) &
                         +  ZCPHASX(:,:)            * PUT(IIE  ,:,:)            &
              + (   ZLBEU (:,:) * PTSTEP                                        &
                    +  ZLBGU (:,:) * ZCPHASX(:,:)                               &
                    +  ZKTSTEP*ZLBXU(:,:)   )  )   
!
!
!
END SELECT
!
END IF
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       4.    LBC FILLING IN THE Y DIRECTION (BOTTOM SOUTH SIDE):   
!              ------------------------------
!       ====>  It only concernes V component 
!                                ----------- 
!
IF (LSOUTH_ll( )) THEN
!
SELECT CASE ( HLBCY(1) )            
!
!*       4.1  WALL CASE: 
!             ========= 
!
  CASE ('WALL')
!
       PRVS (:,IJB,:) = 0.
!
!*       4.2  OPEN CASE:
!             ========= 
!
  CASE ('OPEN')
!
    IF (HTURB /= "NONE" ) THEN
       CALL CPHASE_PROFILE(PZHAT,PCPHASE,PCPHASE_PBL,ZPHASY,PTKET(:,IJB,:))
    ELSE
       CALL CPHASE_PROFILE(PZHAT,PCPHASE,PCPHASE_PBL,ZPHASY)
    END IF    
!
    ZCPHASY(:,:) = MAX ( 0., MIN ( 1.,                                      &
                    (-PVT(:,IJB,:) + ZPHASY(:,:) ) * PTSTEP/ PDYHAT(IJB)   )      )
!
    IF ( SIZE(PLBYVS,1) == 0 ) THEN
      ZLBEV (:,:) = 0.
      ZLBGV (:,:) = PLBYVM(:,JPHEXT+1,:) - PLBYVM(:,JPHEXT,:) 
      ZLBYV(:,:)  = PLBYVM(:,JPHEXT,:)
    ELSE
      ZLBEV (:,:) = PLBYVS(:,JPHEXT,:)
      ZLBGV (:,:) = PLBYVM(:,JPHEXT+1,:) - PLBYVM(:,JPHEXT,:) +  &
                      PTSTEP * (PLBYVS(:,JPHEXT+1,:) - PLBYVS(:,JPHEXT,:))
      ZLBYV(:,:)  = PLBYVM(:,JPHEXT,:) + PTSTEP * PLBYVS(:,JPHEXT,:)
    END IF
!  
!     ============================================================
!
!     PRVS (:,IJB,:) =(PRHODJ(:,IJB-1,:) + PRHODJ(:,IJB,:)) * 0.5 *           &
!                    (  (1. - ZCPHASY(:,:) - ZKTSTEP) * PVM(:,IJB  ,:)      &
!                      + 2. * ZCPHASY(:,:)  * PVT(:,IJB+1,:)                &
!           +2.* (   ZLBEV (:,:) * ZTSTEP                                   &
!                 -  ZLBGV (:,:) * ZCPHASY(:,:)                             &
!                 +  ZKTSTEP*ZLBYV(:,:)       )                             &
!                    ) * ZINVTSTEP / (1.+ ZCPHASY(:,:) +ZKTSTEP)
      PRVS (:,IJB,:) =(PRHODJ(:,IJB-1,:) + PRHODJ(:,IJB,:)) * 0.5 *        &
                       ZINVTSTEP / (1.+ ZKTSTEP * ZALPHA2 )  *             &
          (  (1. - ZCPHASY(:,:) - ZKTSTEP * (1. - ZALPHA2) ) * PVT(:,IJB,:)&
                         +  ZCPHASY(:,:)            * PVT(:,IJB+1,:)       &
              + (   ZLBEV (:,:) * PTSTEP                                   &
                    -  ZLBGV (:,:) * ZCPHASY(:,:)                          &
                    +  ZKTSTEP*ZLBYV(:,:)   )  )   
!
!
!
!
END SELECT
!
END IF
!-------------------------------------------------------------------------------
!
!*       5.    LBC FILLING IN THE Y DIRECTION (TOP NORTH SIDE):   
!              ------------------------------
!       ====>  It only concernes V component 
!                                ----------- 
!
IF (LNORTH_ll( )) THEN
!
!
SELECT CASE ( HLBCY(2) )    
!
!*       5.1  WALL CASE:
!             ========= 
!
  CASE ('WALL')
!
       PRVS (:,IJE+1,:) = 0.
!
!*       5.2  OPEN CASE:
!             ========= 
!
  CASE ('OPEN')
!
    IF (HTURB /= "NONE" ) THEN
       CALL CPHASE_PROFILE(PZHAT,PCPHASE,PCPHASE_PBL,ZPHASY,PTKET(:,IJE,:))
    ELSE
       CALL CPHASE_PROFILE(PZHAT,PCPHASE,PCPHASE_PBL,ZPHASY)
    END IF    
!
    ZCPHASY(:,:) = MAX ( 0., MIN ( 1.,                                      &
                    ( PVT(:,IJE+1,:) + ZPHASY(:,:) ) * PTSTEP/PDYHAT(IJE)  )      )
!
    ILBY=SIZE(PLBYVM,2)
    IF ( SIZE(PLBYVS,1) == 0 ) THEN
      ZLBEV (:,:) = 0.
      ZLBGV (:,:) = PLBYVM(:,ILBY-JPHEXT+1,:) - PLBYVM(:,ILBY-JPHEXT,:) 
      ZLBYV(:,:)  = PLBYVM(:,ILBY-JPHEXT+1,:)
    ELSE
      ZLBEV (:,:) = PLBYVS(:,ILBY-JPHEXT+1,:)
      ZLBGV (:,:) = PLBYVM(:,ILBY-JPHEXT+1,:) - PLBYVM(:,ILBY-JPHEXT,:) +  &
                      PTSTEP * (PLBYVS(:,ILBY-JPHEXT+1,:) - PLBYVS(:,ILBY-JPHEXT,:))
      ZLBYV(:,:)  = PLBYVM(:,ILBY-JPHEXT+1,:) + PTSTEP * PLBYVS(:,ILBY-JPHEXT+1,:)
    END IF
!  
!     ============================================================
!
!     PRVS (:,IJE+1,:) =(PRHODJ(:,IJE+1,:) + PRHODJ(:,IJE,:)) * 0.5 *           &
!                      (  (1. - ZCPHASY(:,:) - ZKTSTEP) * PVM(:,IJE+1,:)      &
!                        + 2. * ZCPHASY(:,:)  * PVT(:,IJE  ,:)                &
!             +2.* (   ZLBEV (:,:) * ZTSTEP                                   &
!                   +  ZLBGV (:,:) * ZCPHASY(:,:)                             &
!                   +  ZKTSTEP* ZLBYV(:,:)   )                                &
!                      ) * ZINVTSTEP / (1.+ ZCPHASY(:,:) +ZKTSTEP)
!
      PRVS (:,IJE+1,:) =(PRHODJ(:,IJE+1,:) + PRHODJ(:,IJE,:)) * 0.5 *         &
                       ZINVTSTEP / (1.+ ZKTSTEP * ZALPHA2 )  *                &
           (  (1. - ZCPHASY(:,:) - ZKTSTEP * (1. - ZALPHA2) ) * PVT(:,IJE+1,:)&
                         +  ZCPHASY(:,:)            * PVT(:,IJE,:)            &
              + (   ZLBEV (:,:) * PTSTEP                                      &
                    +  ZLBGV (:,:) * ZCPHASY(:,:)                             &
                    +  ZKTSTEP*ZLBYV(:,:)   )  )   
! 
!
END SELECT
!
END IF
!-------------------------------------------------------------------------------
!
!*       6.    UPPER BOUNDARY (FLAT SURFACE): 
!              ------------------------------
!
PRWS (:,:,IKE+1) = 0.
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE RAD_BOUND
