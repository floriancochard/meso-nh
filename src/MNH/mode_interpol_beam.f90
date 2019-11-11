!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!      #########################
       MODULE MODE_INTERPOL_BEAM 
!      #########################
!
!
!!**** INTERPOL_BEAM * * - interpolates model fields on radar grid
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to interpolate model fields on the radar
!!    grid, i.e. along the ray path.
!!
!!**  METHOD
!!    ------
!!      A bilinear method is used. The module is subdivided in two: 
!!    MODE_INTERPOL_BEAM_S interpolates model fields at one point; 
!!    MODE_INTERPOL_BEAM_A interpolates model fields on a whole radar grid. 
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq      * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    29/03/04
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!
!-----------------------------------------------------------------
USE MODD_RADAR, ONLY: NBRAD,NBELEV,NBAZIM,NBSTEPMAX,NPTS_H,NPTS_V

TYPE PAMOD
   REAL,DIMENSION(:,:,:),POINTER :: P
END TYPE PAMOD
TYPE PARAD
   REAL,DIMENSION(:,:,:,:,:,:),POINTER :: P
END TYPE PARAD

INTERFACE INTERPOL_BEAM
   MODULE PROCEDURE INTERPOL_BEAM_S, INTERPOL_BEAM_A, INTERPOL_BEAM_P
END INTERFACE

CONTAINS
!
!     #########################################################################
       SUBROUTINE INTERPOL_BEAM_S(PA,PB,PX_RAY,PY_RAY,PZ_RAY,PXHATM,PYHATM,PZM)
!     #########################################################################
! interpolates model grid-point field PA at radar bin value PB
!-------------------------------------------------------------------------------

!
!*       0.    DECLARATIONS
!              ------------

    USE MODD_PARAMETERS
    USE MODD_GRID_n
    USE MODE_ll
!
    IMPLICIT NONE
!
    REAL, DIMENSION(:,:,:),      INTENT(IN)  :: PA
    REAL,                        INTENT(OUT) :: PB 
    REAL,                        INTENT(IN)  :: PX_RAY
    REAL,                        INTENT(IN)  :: PY_RAY
    REAL,                        INTENT(IN)  :: PZ_RAY
    REAL,DIMENSION(:) ,          INTENT(IN)  :: PXHATM
    REAL,DIMENSION(:) ,          INTENT(IN)  :: PYHATM
    REAL,DIMENSION(:,:,:),       INTENT(IN)  :: PZM

!
!   
    INTEGER :: IIB,IIE          ! Loop limits for coordinate X
    INTEGER :: IJB,IJE          ! Loop limits for coordinate Y
    INTEGER :: IKB,IKE          ! Loop limits for coordinate Z
    INTEGER :: IIU,IJU,IKU      ! Loop variables of model 
    INTEGER :: II,IJ,IK00,IK10,IK01,IK11 ! local values of the previous arrays
    INTEGER :: IPAS ! index : 0 if the point of the ray is outside the model domain
                    ! 1 if the point of the ray is between ground surface and first mass level
                    ! 2 if the point of the ray is in the domain above the first mass level
    REAL    :: ZXCOEF ! x-coefficients for the interpolation
    REAL    :: ZYCOEF ! y-coefficients for the interpolation
    REAL    :: ZZCOEF00,ZZCOEF01,ZZCOEF10,ZZCOEF11 ! z-coefficients for the interpolation
    REAL    :: ZSURF00,ZSURF01,ZSURF10,ZSURF11  ! 4 nearest points surface values
!
    IPAS=0
!
    IIU=SIZE(PZM,1)
    IJU=SIZE(PZM,2)
    IKU=SIZE(PZM,3)
    CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
    IKB = JPVEXT + 1
    IKE = IKU - JPVEXT
! 
    
    II=COUNT(PXHATM(:) <= PX_RAY) ! number of mass points x-coordinates less than x-position of current ray point
    IJ=COUNT(PYHATM(:) <= PY_RAY)
    IF ( (II  <= IIE-1) .AND. (II >= IIB) .AND. (IJ <= IJE-1) .AND. (IJ >= IJB) ) THEN
       !          WRITE(ILUOUT0,*) 'inside the horizontal domain '
      ZXCOEF=(PX_RAY-PXHATM(II))/(PXHATM(II+1)-PXHATM(II))
       ZYCOEF=(PY_RAY-PYHATM(IJ))/(PYHATM(IJ+1)-PYHATM(IJ))
       ! compute nearest vertical level below the nearest horizontal points (resp.)
       IK00=COUNT(PZM(II,IJ,:)    <= PZ_RAY)
       IK10=COUNT(PZM(II+1,IJ,:)  <= PZ_RAY)
       IK01=COUNT(PZM(II,IJ+1,:)  <= PZ_RAY)
       IK11=COUNT(PZM(II+1,IJ+1,:) <= PZ_RAY)
       !         
       IF (IK00 < IKB .OR. IK01 < IKB .OR. IK10 < IKB  .OR. IK11 < IKB ) THEN
         ! We are below the lowest mass level           
          IF ((PZ_RAY >= XZS(II,IJ)).AND.(PZ_RAY >= XZS(II+1,IJ)) &
               .AND.(PZ_RAY >= XZS(II,IJ+1)).AND.(PZ_RAY >= XZS(II+1,IJ+1))) THEN  
             ! we are somewhere between the IKB mass level and the ground surface
             ! extrapolation from the lowest mass level
            IPAS=1
            ZZCOEF00=(PZ_RAY-XZS(II,IJ))/(PZM(II,IJ,IKB)-XZS(II,IJ))
            ZZCOEF10=(PZ_RAY-XZS(II+1,IJ))/(PZM(II+1,IJ,IKB)-XZS(II+1,IJ))
            ZZCOEF01=(PZ_RAY-XZS(II,IJ+1))/(PZM(II,IJ+1,IKB)-XZS(II,IJ+1))
            ZZCOEF11=(PZ_RAY-XZS(II+1,IJ+1))/(PZM(II+1,IJ+1,IKB)-XZS(II+1,IJ+1))
                                !
            ZSURF00=PA(II,IJ,IKB) -(PZM(II,IJ,IKB)     -XZS(II,IJ)    )*(PA(II,IJ,IKB+1)    -PA(II,IJ,IKB)    ) &
                 /(PZM(II,IJ,IKB+1)    -PZM(II,IJ,IKB)    )
            ZSURF01=PA(II,IJ+1,IKB)  - (PZM(II,IJ+1,IKB)  -XZS(II,IJ+1)  )*(PA(II,IJ+1,IKB+1)  -PA(II,IJ+1,IKB)  ) &
                 /(PZM(II,IJ+1,IKB+1)  -PZM(II,IJ+1,IKB)  )
            ZSURF10=PA(II+1,IJ,IKB)  - (PZM(II+1,IJ,IKB)  -XZS(II+1,IJ)  )*(PA(II+1,IJ,IKB+1)  -PA(II+1,IJ,IKB)  ) &
                 /(PZM(II+1,IJ,IKB+1)  -PZM(II+1,IJ,IKB)  )
            ZSURF11=PA(II+1,IJ+1,IKB)- (PZM(II+1,IJ+1,IKB)-XZS(II+1,IJ+1))*(PA(II+1,IJ+1,IKB+1)-PA(II+1,IJ+1,IKB)) &
                 /(PZM(II+1,IJ+1,IKB+1)-PZM(II+1,IJ+1,IKB))
                                !
                        PB =  (1.-ZYCOEF)*(1.-ZXCOEF)*((1.-ZZCOEF00)*ZSURF00+ZZCOEF00*PA(II  ,IJ  ,IKB)) &
                 +(1.-ZYCOEF)*(   ZXCOEF)*((1.-ZZCOEF10)*ZSURF10+ZZCOEF10*PA(II+1,IJ  ,IKB)) &
                 +(   ZYCOEF)*(1.-ZXCOEF)*((1.-ZZCOEF01)*ZSURF01+ZZCOEF01*PA(II  ,IJ+1,IKB)) &
                 +(   ZYCOEF)*(   ZXCOEF)*((1.-ZZCOEF11)*ZSURF11+ZZCOEF11*PA(II+1,IJ+1,IKB))
          ENDIF
       ELSE
         IF ((IK00 <= IKE-1).AND.  (IK01 <= IKE-1) .AND. (IK10 <= IKE-1) .AND. (IK11 <= IKE-1) ) THEN
         ! We are above below the lowest mass level and below the upper mass level
           IPAS=2
           ZZCOEF00=(PZ_RAY -PZM(II,IJ,IK00))    /(PZM(II,IJ,IK00+1)-PZM(II,IJ,IK00))
           ZZCOEF10=(PZ_RAY -PZM(II+1,IJ,IK10))  /(PZM(II+1,IJ,IK10+1)-PZM(II+1,IJ,IK10))
           ZZCOEF01=(PZ_RAY -PZM(II,IJ+1,IK01))  /(PZM(II,IJ+1,IK01+1)-PZM(II,IJ+1,IK01))
           ZZCOEF11=(PZ_RAY -PZM(II+1,IJ+1,IK11))/(PZM(II+1,IJ+1,IK11+1)-PZM(II+1,IJ+1,IK11))
                                !
           PB =  (1.-ZYCOEF)*(1.-ZXCOEF)*((1.-ZZCOEF00)*PA(II  ,IJ  ,IK00)+ZZCOEF00*PA(II  ,IJ  ,IK00+1)) &
                +(1.-ZYCOEF)*    ZXCOEF *((1.-ZZCOEF10)*PA(II+1,IJ  ,IK10)+ZZCOEF10*PA(II+1,IJ  ,IK10+1)) &
                +(   ZYCOEF)*(1.-ZXCOEF)*((1.-ZZCOEF01)*PA(II  ,IJ+1,IK01)+ZZCOEF01*PA(II  ,IJ+1,IK01+1)) &
                +(   ZYCOEF)*(   ZXCOEF)*((1.-ZZCOEF11)*PA(II+1,IJ+1,IK11)+ZZCOEF11*PA(II+1,IJ+1,IK11+1))
         END IF
       END IF
     END IF

    IF(IPAS == 0) PB=-XUNDEF
    
  END SUBROUTINE INTERPOL_BEAM_S

!     #########################################################################
       SUBROUTINE INTERPOL_BEAM_A(PA,PB,PX_RAY,PY_RAY,PZ_RAY,PXHATM,PYHATM,PZM)
!     #########################################################################
! interpolates model grid-point field PA at radar bin field PB

!
!*       0.    DECLARATIONS
!              ------------

    USE MODD_PARAMETERS
    USE MODD_GRID_n
    USE MODE_ll
!
    IMPLICIT NONE
!
    REAL, DIMENSION(:,:,:),    INTENT(IN)  :: PA
    REAL, DIMENSION(:,:,:,:),  INTENT(OUT) :: PB 
    REAL, DIMENSION(:,:,:,:),    INTENT(IN)  :: PX_RAY
    REAL, DIMENSION(:,:,:,:),    INTENT(IN)  :: PY_RAY
    REAL, DIMENSION(:,:,:,:),     INTENT(IN)  :: PZ_RAY
    REAL,DIMENSION(:) ,          INTENT(IN)  :: PXHATM
    REAL,DIMENSION(:) ,          INTENT(IN)  :: PYHATM
    REAL,DIMENSION(:,:,:),       INTENT(IN)  :: PZM

!
!   
    INTEGER :: IIB,IIE          ! Loop limits for coordinate X
    INTEGER :: IJB,IJE          ! Loop limits for coordinate Y
    INTEGER :: IKB,IKE          ! Loop limits for coordinate Z
    INTEGER :: IIU,IJU,IKU      ! Loop variables of model 
    INTEGER :: II,IJ,IK00,IK10,IK01,IK11 ! local values of the previous arrays
    INTEGER :: IPAS ! index : 0 if the point of the ray is outside the model domain
                    ! 1 if the point of the ray is between ground surface and first mass level
                    ! 2 if the point of the ray is in the domain above the first mass level
    REAL    :: ZXCOEF ! x-coefficients for the interpolation
    REAL    :: ZYCOEF ! y-coefficients for the interpolation
    REAL    :: ZZCOEF00,ZZCOEF01,ZZCOEF10,ZZCOEF11 ! z-coefficients for the interpolation
    REAL    :: ZSURF00,ZSURF01,ZSURF10,ZSURF11  ! 4 nearest points surface values
    INTEGER :: INBAZIM,INBSTEP,INPTS_GH_H, INPTS_GH_V
    INTEGER :: JAZ,JH,JV,JL
!
    INBAZIM=SIZE(PB,1)
    INBSTEP=SIZE(PB,2)
    INPTS_GH_H=SIZE(PB,3)
    INPTS_GH_V=SIZE(PB,4)
!
    IPAS=0
!
    IIU=SIZE(PZM,1)
    IJU=SIZE(PZM,2)
    IKU=SIZE(PZM,3)
    CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
    IKB = JPVEXT + 1
    IKE = IKU - JPVEXT
! 
    DO JAZ=1, INBAZIM
       DO JL=1,INBSTEP
          DO JH=1,INPTS_GH_H
          DO JV=1,INPTS_GH_V
            IPAS=0
            II=COUNT(PXHATM(:) <= PX_RAY(JAZ,JL,JH,JV)) ! number of mass points x-coordinates less than x-position of current ray point
            IJ=COUNT(PYHATM(:) <= PY_RAY(JAZ,JL,JH,JV))
            IF ( (II  <= IIE-1) .AND. (II >= IIB) .AND. (IJ <= IJE-1) .AND. (IJ >= IJB) ) THEN
       !          WRITE(ILUOUT0,*) 'inside the horizontal domain '
              ZXCOEF=(PX_RAY(JAZ,JL,JH,JV)-PXHATM(II))/(PXHATM(II+1)-PXHATM(II))
              ZYCOEF=(PY_RAY(JAZ,JL,JH,JV)-PYHATM(IJ))/(PYHATM(IJ+1)-PYHATM(IJ))
       ! compute nearest vertical level below the nearest horizontal points (resp.)
              IK00=COUNT(PZM(II,IJ,:)    <= PZ_RAY(JAZ,JL,JH,JV))
              IK10=COUNT(PZM(II+1,IJ,:)  <= PZ_RAY(JAZ,JL,JH,JV))
              IK01=COUNT(PZM(II,IJ+1,:)  <= PZ_RAY(JAZ,JL,JH,JV))
              IK11=COUNT(PZM(II+1,IJ+1,:) <= PZ_RAY(JAZ,JL,JH,JV))
       !         
              IF (IK00 < IKB .OR. IK01 < IKB .OR. IK10 < IKB  .OR. IK11 < IKB ) THEN
              ! We are below the lowest mass level           
                IF ((PZ_RAY(JAZ,JL,JH,JV) >= XZS(II,IJ)).AND.(PZ_RAY(JAZ,JL,JH,JV) >= XZS(II+1,IJ)) &
               .AND.(PZ_RAY(JAZ,JL,JH,JV) >= XZS(II,IJ+1)).AND.(PZ_RAY(JAZ,JL,JH,JV) >= XZS(II+1,IJ+1))) THEN  
             ! we are somewhere between the IKB mass level and the ground surface
             ! extrapolation from the lowest mass level
                  IPAS=1
                  ZZCOEF00=(PZ_RAY(JAZ,JL,JH,JV)-XZS(II,IJ))/(PZM(II,IJ,IKB)-XZS(II,IJ))
                  ZZCOEF10=(PZ_RAY(JAZ,JL,JH,JV)-XZS(II+1,IJ))/(PZM(II+1,IJ,IKB)-XZS(II+1,IJ))
                  ZZCOEF01=(PZ_RAY(JAZ,JL,JH,JV)-XZS(II,IJ+1))/(PZM(II,IJ+1,IKB)-XZS(II,IJ+1))
                  ZZCOEF11=(PZ_RAY(JAZ,JL,JH,JV)-XZS(II+1,IJ+1))/(PZM(II+1,IJ+1,IKB)-XZS(II+1,IJ+1))
                                !
                  ZSURF00=PA(II,IJ,IKB) -(PZM(II,IJ,IKB) -XZS(II,IJ) )* &
                       (PA(II,IJ,IKB+1)    -PA(II,IJ,IKB))/(PZM(II,IJ,IKB+1) -PZM(II,IJ,IKB)    )
                  ZSURF01=PA(II,IJ+1,IKB) - (PZM(II,IJ+1,IKB)-XZS(II,IJ+1))* &
                       (PA(II,IJ+1,IKB+1)-PA(II,IJ+1,IKB)) /(PZM(II,IJ+1,IKB+1)  -PZM(II,IJ+1,IKB)  )
                  ZSURF10=PA(II+1,IJ,IKB)-(PZM(II+1,IJ,IKB)-XZS(II+1,IJ))* &
                       (PA(II+1,IJ,IKB+1)-PA(II+1,IJ,IKB))/(PZM(II+1,IJ,IKB+1)  -PZM(II+1,IJ,IKB)  )
                  ZSURF11=PA(II+1,IJ+1,IKB)- (PZM(II+1,IJ+1,IKB)-XZS(II+1,IJ+1))* &
                       (PA(II+1,IJ+1,IKB+1)-PA(II+1,IJ+1,IKB))/(PZM(II+1,IJ+1,IKB+1)-PZM(II+1,IJ+1,IKB))
                  !
                  PB(JAZ,JL,JH,JV)=(1.-ZYCOEF)*(1.-ZXCOEF)*((1.-ZZCOEF00)*ZSURF00+ZZCOEF00*PA(II,IJ,IKB)) &
                       +(1.-ZYCOEF)*(   ZXCOEF)*((1.-ZZCOEF10)*ZSURF10+ZZCOEF10*PA(II+1,IJ  ,IKB)) &
                       +(   ZYCOEF)*(1.-ZXCOEF)*((1.-ZZCOEF01)*ZSURF01+ZZCOEF01*PA(II  ,IJ+1,IKB)) &
                       +(   ZYCOEF)*(   ZXCOEF)*((1.-ZZCOEF11)*ZSURF11+ZZCOEF11*PA(II+1,IJ+1,IKB))
                ENDIF
              ELSE
                IF ((IK00 <= IKE-1).AND.  (IK01 <= IKE-1) .AND. (IK10 <= IKE-1) .AND. (IK11 <= IKE-1) ) THEN
         ! We are above below the lowest mass level and below the upper mass level
                  IPAS=2
                  ZZCOEF00=(PZ_RAY(JAZ,JL,JH,JV) -PZM(II,IJ,IK00))    /(PZM(II,IJ,IK00+1)-PZM(II,IJ,IK00))
                  ZZCOEF10=(PZ_RAY(JAZ,JL,JH,JV) -PZM(II+1,IJ,IK10))  /(PZM(II+1,IJ,IK10+1)-PZM(II+1,IJ,IK10))
                  ZZCOEF01=(PZ_RAY(JAZ,JL,JH,JV) -PZM(II,IJ+1,IK01))  /(PZM(II,IJ+1,IK01+1)-PZM(II,IJ+1,IK01))
                  ZZCOEF11=(PZ_RAY(JAZ,JL,JH,JV) -PZM(II+1,IJ+1,IK11))/(PZM(II+1,IJ+1,IK11+1)-PZM(II+1,IJ+1,IK11))

                  PB(JAZ,JL,JH,JV)=(1.-ZYCOEF)*(1.-ZXCOEF)*((1.-ZZCOEF00)*PA(II  ,IJ  ,IK00)+&
                       ZZCOEF00*PA(II  ,IJ  ,IK00+1)) &
                       +(1.-ZYCOEF)*    ZXCOEF *((1.-ZZCOEF10)*PA(II+1,IJ  ,IK10)+&
                       ZZCOEF10*PA(II+1,IJ  ,IK10+1)) &
                       +(   ZYCOEF)*(1.-ZXCOEF)*((1.-ZZCOEF01)*PA(II  ,IJ+1,IK01)+&
                       ZZCOEF01*PA(II  ,IJ+1,IK01+1)) &
                       +(   ZYCOEF)*(   ZXCOEF)*((1.-ZZCOEF11)*PA(II+1,IJ+1,IK11)+&
                       ZZCOEF11*PA(II+1,IJ+1,IK11+1))
                END IF
              END IF
            END IF

            IF(IPAS == 0) PB(JAZ,JL,JH,JV)=-XUNDEF
          END DO
        END DO
      END DO
    END DO
!
  END SUBROUTINE INTERPOL_BEAM_A

!     #########################################################################
       SUBROUTINE INTERPOL_BEAM_P(PA,PB,PX_RAY,PY_RAY,PZ_RAY,PXHATM,PYHATM,PZM)
!     #########################################################################
! interpolates model grid-point field PA at radar bin field PB

!
!*       0.    DECLARATIONS
!              ------------

    USE MODD_PARAMETERS
    USE MODD_GRID_n
    USE MODE_ll
!
    IMPLICIT NONE
!
    TYPE(PAMOD),DIMENSION(:),INTENT(IN)    :: PA
    TYPE(PARAD),DIMENSION(:),INTENT(INOUT) :: PB
    REAL, DIMENSION(:,:,:,:,:,:),INTENT(IN):: PX_RAY
    REAL, DIMENSION(:,:,:,:,:,:),INTENT(IN):: PY_RAY
    REAL, DIMENSION(:,:,:,:,:,:),INTENT(IN):: PZ_RAY
    REAL,DIMENSION(:) ,      INTENT(IN)    :: PXHATM
    REAL,DIMENSION(:) ,      INTENT(IN)    :: PYHATM
    REAL,DIMENSION(:,:,:),   INTENT(IN)    :: PZM
!   
    INTEGER :: IIB,IIE          ! Loop limits for coordinate X
    INTEGER :: IJB,IJE          ! Loop limits for coordinate Y
    INTEGER :: IKB,IKE          ! Loop limits for coordinate Z
    INTEGER :: IIU,IJU,IKU      ! Loop variables of model 
    INTEGER :: II,IJ,IK00,IK10,IK01,IK11 ! local values of the previous arrays
    INTEGER :: IPAS ! index : 0 if the point of the ray is outside the model domain
                    ! 1 if the point of the ray is between ground surface and first mass level
                    ! 2 if the point of the ray is in the domain above the first mass level
    REAL    :: ZXCOEF ! x-coefficients for the interpolation
    REAL    :: ZYCOEF ! y-coefficients for the interpolation
    REAL    :: ZZCOEF00,ZZCOEF01,ZZCOEF10,ZZCOEF11 ! z-coefficients for the interpolation
    REAL    :: ZSURF00,ZSURF01,ZSURF10,ZSURF11  ! 4 nearest points surface values
    INTEGER :: INVAR,IEL
    INTEGER :: JI,JEL,JAZ,JH,JV,JL,JN
!
    
!    NBAZIM=SIZE(PB(1)%P,1)
    INVAR=SIZE(PB,1)
!
    IPAS=0
!
    CALL GET_INDICE_ll( IIB,IJB,IIE,IJE)
    IKU=SIZE(PZM,3)
    IKB = JPVEXT + 1
    IKE = IKU - JPVEXT
! 
    
    DO JI=1,NBRAD
       IEL=NBELEV(JI)
       DO JEL=1,IEL
          DO JAZ=1, NBAZIM
             DO JL=1,NBSTEPMAX+1
                DO JH=1,NPTS_H
                   DO JV=1,NPTS_V
            IPAS=0
            II=COUNT(PXHATM(:) <= PX_RAY(JI,JEL,JAZ,JL,JH,JV)) ! number of mass points x-coordinates less than x-position of current ray point
            IJ=COUNT(PYHATM(:) <= PY_RAY(JI,JEL,JAZ,JL,JH,JV))
            IF ( (II  <= IIE) .AND. (II >= IIB) .AND. (IJ <= IJE) .AND. (IJ >= IJB) ) THEN
       !          WRITE(ILUOUT0,*) 'inside the horizontal domain '
              ZXCOEF=(PX_RAY(JI,JEL,JAZ,JL,JH,JV)-PXHATM(II))/(PXHATM(II+1)-PXHATM(II))
              ZYCOEF=(PY_RAY(JI,JEL,JAZ,JL,JH,JV)-PYHATM(IJ))/(PYHATM(IJ+1)-PYHATM(IJ))
       ! compute nearest vertical level below the nearest horizontal points (resp.)
              IK00=COUNT(PZM(II,IJ,:)    <= PZ_RAY(JI,JEL,JAZ,JL,JH,JV))
              IK10=COUNT(PZM(II+1,IJ,:)  <= PZ_RAY(JI,JEL,JAZ,JL,JH,JV))
              IK01=COUNT(PZM(II,IJ+1,:)  <= PZ_RAY(JI,JEL,JAZ,JL,JH,JV))
              IK11=COUNT(PZM(II+1,IJ+1,:) <= PZ_RAY(JI,JEL,JAZ,JL,JH,JV))
       !         
              IF (IK00 < IKB .OR. IK01 < IKB .OR. IK10 < IKB  .OR. IK11 < IKB ) THEN
              ! We are below the lowest mass level           
                IF ((PZ_RAY(JI,JEL,JAZ,JL,JH,JV) >= XZS(II,IJ)).AND.(PZ_RAY(JI,JEL,JAZ,JL,JH,JV) >= XZS(II+1,IJ)) &
               .AND.(PZ_RAY(JI,JEL,JAZ,JL,JH,JV) >= XZS(II,IJ+1)).AND.(PZ_RAY(JI,JEL,JAZ,JL,JH,JV) >= XZS(II+1,IJ+1))) THEN  
             ! we are somewhere between the IKB mass level and the ground surface
             ! extrapolation from the lowest mass level
                  IPAS=1
                  ZZCOEF00=(PZ_RAY(JI,JEL,JAZ,JL,JH,JV)-XZS(II,IJ))/(PZM(II,IJ,IKB)-XZS(II,IJ))
                  ZZCOEF10=(PZ_RAY(JI,JEL,JAZ,JL,JH,JV)-XZS(II+1,IJ))/(PZM(II+1,IJ,IKB)-XZS(II+1,IJ))
                  ZZCOEF01=(PZ_RAY(JI,JEL,JAZ,JL,JH,JV)-XZS(II,IJ+1))/(PZM(II,IJ+1,IKB)-XZS(II,IJ+1))
                  ZZCOEF11=(PZ_RAY(JI,JEL,JAZ,JL,JH,JV)-XZS(II+1,IJ+1))/(PZM(II+1,IJ+1,IKB)-XZS(II+1,IJ+1))
                                !
                  DO JN=1,INVAR
                    ZSURF00=PA(JN)%P(II,IJ,IKB) -(PZM(II,IJ,IKB) -XZS(II,IJ) )* &
                         (PA(JN)%P(II,IJ,IKB+1)    -PA(JN)%P(II,IJ,IKB))/(PZM(II,IJ,IKB+1) -PZM(II,IJ,IKB)    )
                    ZSURF01=PA(JN)%P(II,IJ+1,IKB) - (PZM(II,IJ+1,IKB)-XZS(II,IJ+1))* &
                         (PA(JN)%P(II,IJ+1,IKB+1)-PA(JN)%P(II,IJ+1,IKB)) /(PZM(II,IJ+1,IKB+1)  -PZM(II,IJ+1,IKB)  )
                    ZSURF10=PA(JN)%P(II+1,IJ,IKB)-(PZM(II+1,IJ,IKB)-XZS(II+1,IJ))* &
                         (PA(JN)%P(II+1,IJ,IKB+1)-PA(JN)%P(II+1,IJ,IKB))/(PZM(II+1,IJ,IKB+1)  -PZM(II+1,IJ,IKB)  )
                    ZSURF11=PA(JN)%P(II+1,IJ+1,IKB)- (PZM(II+1,IJ+1,IKB)-XZS(II+1,IJ+1))* &
                 (PA(JN)%P(II+1,IJ+1,IKB+1)-PA(JN)%P(II+1,IJ+1,IKB))/(PZM(II+1,IJ+1,IKB+1)-PZM(II+1,IJ+1,IKB))
                                !
                    PB(JN)%P(JI,JEL,JAZ,JL,JH,JV)=(1.-ZYCOEF)*(1.-ZXCOEF)*((1.-ZZCOEF00)*ZSURF00+ZZCOEF00*PA(JN)%P(II,IJ,IKB)) &
                         +(1.-ZYCOEF)*(   ZXCOEF)*((1.-ZZCOEF10)*ZSURF10+ZZCOEF10*PA(JN)%P(II+1,IJ  ,IKB)) &
                         +(   ZYCOEF)*(1.-ZXCOEF)*((1.-ZZCOEF01)*ZSURF01+ZZCOEF01*PA(JN)%P(II  ,IJ+1,IKB)) &
                         +(   ZYCOEF)*(   ZXCOEF)*((1.-ZZCOEF11)*ZSURF11+ZZCOEF11*PA(JN)%P(II+1,IJ+1,IKB))
                  END DO
                ENDIF
              ELSE
                IF ((IK00 <= IKE).AND.  (IK01 <= IKE) .AND. (IK10 <= IKE) .AND. (IK11 <= IKE) ) THEN
         ! We are above below the lowest mass level and below the upper mass level
                  IPAS=2
                  ZZCOEF00=(PZ_RAY(JI,JEL,JAZ,JL,JH,JV) -PZM(II,IJ,IK00))    /(PZM(II,IJ,IK00+1)-PZM(II,IJ,IK00))
                  ZZCOEF10=(PZ_RAY(JI,JEL,JAZ,JL,JH,JV) -PZM(II+1,IJ,IK10))  /(PZM(II+1,IJ,IK10+1)-PZM(II+1,IJ,IK10))
                  ZZCOEF01=(PZ_RAY(JI,JEL,JAZ,JL,JH,JV) -PZM(II,IJ+1,IK01))  /(PZM(II,IJ+1,IK01+1)-PZM(II,IJ+1,IK01))
                  ZZCOEF11=(PZ_RAY(JI,JEL,JAZ,JL,JH,JV) -PZM(II+1,IJ+1,IK11))/(PZM(II+1,IJ+1,IK11+1)-PZM(II+1,IJ+1,IK11))
                  DO JN=1,INVAR               !
                    PB(JN)%P(JI,JEL,JAZ,JL,JH,JV)=(1.-ZYCOEF)*(1.-ZXCOEF)*((1.-ZZCOEF00)*PA(JN)%P(II  ,IJ  ,IK00)+&
                         ZZCOEF00*PA(JN)%P(II  ,IJ  ,IK00+1)) &
                         +(1.-ZYCOEF)*    ZXCOEF *((1.-ZZCOEF10)*PA(JN)%P(II+1,IJ  ,IK10)+&
                         ZZCOEF10*PA(JN)%P(II+1,IJ  ,IK10+1)) &
                         +(   ZYCOEF)*(1.-ZXCOEF)*((1.-ZZCOEF01)*PA(JN)%P(II  ,IJ+1,IK01)+&
                         ZZCOEF01*PA(JN)%P(II  ,IJ+1,IK01+1)) &
                         +(   ZYCOEF)*(   ZXCOEF)*((1.-ZZCOEF11)*PA(JN)%P(II+1,IJ+1,IK11)+&
                         ZZCOEF11*PA(JN)%P(II+1,IJ+1,IK11+1))
                  END DO
                END IF
              END IF
            END IF

            IF(IPAS == 0) THEN
               DO JN=1,INVAR
                  PB(JN)%P(JI,JEL,JAZ,JL,JH,JV)=-XUNDEF
               END DO
            END IF
          END DO
        END DO
      END DO
    END DO
 END DO
END DO
!
  END SUBROUTINE INTERPOL_BEAM_P

END MODULE MODE_INTERPOL_BEAM
