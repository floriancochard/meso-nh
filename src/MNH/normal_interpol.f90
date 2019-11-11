!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 soil 2006/05/18 13:07:25
!-----------------------------------------------------------------
!    ###########################
     MODULE MODI_NORMAL_INTERPOL
!    ###########################
!
INTERFACE
!
      SUBROUTINE NORMAL_INTERPOL(PTHT,PRVT,PPABST,                   &
                 PDIRCOSXW, PDIRCOSYW, PDIRCOSZW,                    &
                 PCOSSLOPE,PSINSLOPE,                                &
                 PDXX,PDYY,PDZZ,                                     &
                 PTHA,PRVA,PEXNA                                     )
!
!
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHT,PRVT,PPABST ! 3D fields at time t:
                                 ! PTHT  = potential temperature
                                 ! PRVT  = vapor mixing ratio
                                 ! PPABST= pressure
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDXX, PDYY, PDZZ
                                 ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PTHA,PRVA,PEXNA ! 2D fields at time t:
                                 ! PTHA = potential temperature
                                 ! PRVA = vapor mixing ratio
                                 ! PEXNA=  Exner function
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE NORMAL_INTERPOL
!
END INTERFACE
!
END MODULE MODI_NORMAL_INTERPOL
!
!     ################################################################
      SUBROUTINE NORMAL_INTERPOL(PTHT,PRVT,PPABST,                   &
                 PDIRCOSXW, PDIRCOSYW, PDIRCOSZW,                    &
                 PCOSSLOPE,PSINSLOPE,                                &
                 PDXX,PDYY,PDZZ,                                     &
                 PTHA,PRVA,PEXNA                                     )
!     ################################################################
!
!
!!****  *NORMAL_INTERPOL* - interpolates thermodynamic fields at the intersection
!!           point between  the normal to the surface and the first mass level.
!!
!!    PURPOSE
!!    -------
!!**** 
!        The purpose of this routine is to compute thermodynamic fields at 
!     the point of intersection between the normal 
!     to the orography and the first mass-level hyper-plane at PDZZ(:,:,IKB)/2 
!        
!!**  METHOD
!!    ------
!!       The values of the interpolated fields are determined
!!    by a bilinear interpolation between the 4 nearest points in the first 
!!    mass-level hyper-plane. These points are found according to the signs of 
!!    the slopes' sinus and cosinus. 
!!      
!!        Finally, the horizontal components are set at the marginal points 
!!    according to cyclic boundary conditions because this is the only case
!!    where these points can be considered.
!!
!!    EXTERNAL
!!    --------
!!       NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!       MODD_CONF      : L2D, L1D   switch for 2D model version
!!
!!
!!    REFERENCE
!!    ---------
!!      Book 1 of documentation (Chapter: Turbulence)
!!
!!    AUTHOR
!!    ------
!!      Joel Stein              * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         14/11/95
!!      
!!       J.Stein         28/01/96 : bug correction in the 2D case
!!       J.Stein         22/06/97 : use the absolute pressure 
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_CST
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
! 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHT,PRVT,PPABST ! 3D fields at time t:
                                 ! PTHT  = potential temperature
                                 ! PRVT  = vapor mixing ratio
                                 ! PPABST= pressure
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle 
                                 ! between i and the slope vector
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle 
                                 ! between i and the slope vector
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDXX, PDYY, PDZZ
                                 ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PTHA,PRVA,PEXNA ! 2D fields at time t:
                                 ! PTHA = potential temperature
                                 ! PRVA = vapor mixing ratio
                                 ! PEXNA=  Exner function
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
INTEGER, DIMENSION(SIZE(PDIRCOSXW,1),SIZE(PDIRCOSXW,2)) :: ILOC,JLOC
              ! shift index to find the 4 nearest points in x and y directions
REAL,    DIMENSION(SIZE(PDIRCOSXW,1),SIZE(PDIRCOSXW,2)) :: ZCOEFM,     &
              ! interpolation weigths for mass locations
                                                 ZTHINT,ZRVINT,ZEXNINT 
              ! intermediate values of the 2D fields after x interp.
REAL,    DIMENSION(SIZE(PPABST,1),SIZE(PPABST,2)) :: ZEXNT ! Exner function at t
INTEGER     :: IIB,IIE,IJB,IJE,IKB
              ! index values for the Beginning or the End of the physical 
              ! domain in x,y and z directions
INTEGER     :: IIU,IJU
              ! arrays' sizes for i and j indices
INTEGER     :: JI,JJ      
              ! loop index
!      
!----------------------------------------------------------------------------
!
!*      1.    PRELIMINARIES
!             -------------
!
!
IIB = 2
IJB = 2
IIU = SIZE(PTHT,1)
IJU = SIZE(PTHT,2)
IIE = IIU - 1
IJE = IJU - 1
IKB = 1+JPVEXT
!
ZEXNT(:,:)= (PPABST(:,:,IKB)/XP00) ** (XRD/XCPD)
!
!*      2.    INTERPOLATE THE THERMODYNAMIC FIELDS
!             ------------------------------------
!
ILOC(:,:)=NINT(SIGN(1.,-PCOSSLOPE(:,:)))
JLOC(:,:)=NINT(SIGN(1.,-PSINSLOPE(:,:)))
!
! interpolation in x direction
!
DO JJ = 1,IJU
  DO JI = IIB,IIE
    ZCOEFM(JI,JJ) = 1. - 0.5 * PDZZ(JI,JJ,IKB) * ABS(PDIRCOSXW(JI,JJ))    & 
                             / PDXX(JI+(ILOC(JI,JJ)+1)/2,JJ,IKB)
    ZTHINT(JI,JJ) = ZCOEFM(JI,JJ)     * PTHT(JI,JJ,IKB)              +    &
                   (1.-ZCOEFM(JI,JJ)) * PTHT(JI+ILOC(JI,JJ),JJ,IKB)
    ZRVINT(JI,JJ) = ZCOEFM(JI,JJ)     * PRVT(JI,JJ,IKB)              +    &
                   (1.-ZCOEFM(JI,JJ)) * PRVT(JI+ILOC(JI,JJ),JJ,IKB)
    ZEXNINT(JI,JJ) = ZCOEFM(JI,JJ)    * ZEXNT(JI,JJ)                 +    &
                   (1.-ZCOEFM(JI,JJ)) * ZEXNT(JI+ILOC(JI,JJ),JJ)
  END DO
END DO
!
! interpolation in y direction
!
IF (.NOT. (L2D .OR. L1D)) THEN
!
  DO JJ = IJB,IJE
    DO JI = IIB,IIE
    ZCOEFM(JI,JJ) = 1. - 0.5 * PDZZ(JI,JJ,IKB) * ABS(PDIRCOSYW(JI,JJ))   & 
                             / PDYY(JI,JJ+(JLOC(JI,JJ)+1)/2,IKB)
    PTHA(JI,JJ) = ZCOEFM(JI,JJ)      * ZTHINT(JI,JJ)                +    &
                  (1.-ZCOEFM(JI,JJ)) * ZTHINT(JI,JJ+JLOC(JI,JJ))
    PRVA(JI,JJ) = ZCOEFM(JI,JJ)      * ZRVINT(JI,JJ)                +    &
                  (1.-ZCOEFM(JI,JJ)) * ZRVINT(JI,JJ+JLOC(JI,JJ))
    PEXNA(JI,JJ) = ZCOEFM(JI,JJ)     * ZEXNINT(JI,JJ)               +    &
                  (1.-ZCOEFM(JI,JJ)) * ZEXNINT(JI,JJ+JLOC(JI,JJ))
    END DO
  END DO
! 
ELSE
!
  PTHA(IIB:IIE,IJB:IJE)  = ZTHINT(IIB:IIE,IJB:IJE)
  PRVA(IIB:IIE,IJB:IJE)  = ZRVINT(IIB:IIE,IJB:IJE)
  PEXNA(IIB:IIE,IJB:IJE) = ZEXNINT(IIB:IIE,IJB:IJE)
!
END IF
!
!*      4.    BOUNDARY CONDITIONS
!             -------------------
!
PTHA(IIE,IJE+1)=PTHA(IIE,IJE)
PTHA(IIE,IJB-1)=PTHA(IIE,IJB)
PTHA(IIB,IJE+1)=PTHA(IIB,IJE)
PTHA(IIB,IJB-1)=PTHA(IIB,IJB)
!
PRVA(IIE,IJE+1)=PRVA(IIE,IJE)
PRVA(IIE,IJB-1)=PRVA(IIE,IJB)
PRVA(IIB,IJE+1)=PRVA(IIB,IJE)
PRVA(IIB,IJB-1)=PRVA(IIB,IJB)
!
PEXNA(IIE,IJE+1)=PEXNA(IIE,IJE)
PEXNA(IIE,IJB-1)=PEXNA(IIE,IJB)
PEXNA(IIB,IJE+1)=PEXNA(IIB,IJE)
PEXNA(IIB,IJB-1)=PEXNA(IIB,IJB)
!
PTHA(IIB-1,:)=PTHA(IIE,:)
PTHA(IIE+1,:)=PTHA(IIB,:)
PRVA(IIB-1,:)=PRVA(IIE,:)
PRVA(IIE+1,:)=PRVA(IIB,:)
PEXNA(IIB-1,:)=PEXNA(IIE,:)
PEXNA(IIE+1,:)=PEXNA(IIB,:)
!
PTHA(:,IJB-1)=PTHA(:,IJE)
PTHA(:,IJE+1)=PTHA(:,IJB)
PRVA(:,IJB-1)=PRVA(:,IJE)
PRVA(:,IJE+1)=PRVA(:,IJB)
PEXNA(:,IJB-1)=PEXNA(:,IJE)
PEXNA(:,IJE+1)=PEXNA(:,IJB)
!
!----------------------------------------------------------------------------
!
END SUBROUTINE NORMAL_INTERPOL
