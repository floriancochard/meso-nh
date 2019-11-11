!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODI_ICE4_NUCLEATION_WRAPPER
INTERFACE
SUBROUTINE ICE4_NUCLEATION_WRAPPER(KIT, KJT,KKT, LDMASK, &
                           PTHT, PPABST, PRHODREF, PEXN, PLSFACT, PT, &
                           PRVT, &
                           PCIT, PRVHENI_MR)
IMPLICIT NONE
INTEGER,                        INTENT(IN)    :: KIT, KJT, KKT
LOGICAL, DIMENSION(KIT,KJT,KKT),INTENT(IN)    :: LDMASK
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PTHT    ! Theta at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PT      ! Temperature at time t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(OUT)   :: PRVHENI_MR ! Mixing ratio change due to the heterogeneous nucleation
END SUBROUTINE ICE4_NUCLEATION_WRAPPER
END INTERFACE
END MODULE MODI_ICE4_NUCLEATION_WRAPPER
SUBROUTINE ICE4_NUCLEATION_WRAPPER(KIT, KJT, KKT, LDMASK, &
                                   PTHT, PPABST, PRHODREF, PEXN, PLSFACT, PT, &
                                   PRVT, &
                                   PCIT, PRVHENI_MR)
!!
!!**  PURPOSE
!!    -------
!!      Computes the nucleation
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the splitting of rain_ice source code (nov. 2014)
!!
!!    MODIFICATIONS
!!    -------------
!!
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST, ONLY : XTT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                        INTENT(IN)    :: KIT, KJT, KKT
LOGICAL, DIMENSION(KIT,KJT,KKT),INTENT(IN)    :: LDMASK
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PTHT    ! Theta at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PT      ! Temperature at time t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(OUT)   :: PRVHENI_MR ! Mixing ratio change due to the heterogeneous nucleation
!
!*       0.2  declaration of local variables
!
INTEGER                           :: JL       ! and PACK intrinsics
LOGICAL, DIMENSION(SIZE(PRVT,1),SIZE(PRVT,2),SIZE(PRVT,3)) :: GNEGT  ! Test where to compute the HEN process
INTEGER :: INEGT
INTEGER, DIMENSION(COUNT(PT<XTT .AND. LDMASK)) :: I1,I2,I3 ! Used to replace the COUNT
REAL, DIMENSION(COUNT(PT<XTT .AND. LDMASK))  :: ZZT,      & ! Temperature
                                   ZPRES,    & ! Pressure
                                   ZRVT,     & ! Water vapor m.r. at t
                                   ZCIT,     & ! Pristine ice conc. at t
                                   ZTHT,     & ! Theta at t
                                   ZRHODREF, &
                                   ZEXN,     &
                                   ZLSFACT,  &
                                   ZRVHENI_MR, &
                                   ZB_TH, ZB_RV, ZB_RI
!-------------------------------------------------------------------------------
!
!
!
!  optimization by looking for locations where
!  the temperature is negative only !!!
!
GNEGT(:,:,:)=PT(:,:,:)<XTT .AND. LDMASK
INEGT=0
IF(COUNT(GNEGT)/=0) INEGT=ICE4_NUCLEATION_COUNTJV(GNEGT(:,:,:), KIT, KJT, KKT, SIZE(I1), I1(:), I2(:), I3(:))
!
ZB_TH(:) = 0.
ZB_RV(:) = 0.
ZB_RI(:) = 0.
!
PRVHENI_MR(:,:,:)=0.
IF(INEGT>=1) THEN
  DO JL=1, INEGT
    ZRVT(JL)=PRVT(I1(JL), I2(JL), I3(JL))
    ZCIT(JL)=PCIT(I1(JL), I2(JL), I3(JL))
    ZZT(JL)=PT(I1(JL), I2(JL), I3(JL))
    ZPRES(JL)=PPABST(I1(JL), I2(JL), I3(JL))
    ZTHT(JL)=PTHT(I1(JL), I2(JL), I3(JL))
    ZRHODREF(JL)=PRHODREF(I1(JL), I2(JL), I3(JL))
    ZEXN(JL)=PEXN(I1(JL), I2(JL), I3(JL))
    ZLSFACT(JL)=PLSFACT(I1(JL), I2(JL), I3(JL))
  ENDDO
  CALL ICE4_NUCLEATION(INEGT, .FALSE., ZZT(:)<XTT, &
                       ZTHT, ZPRES, ZRHODREF, ZEXN, ZLSFACT, ZZT, &
                       ZRVT, &
                       ZCIT, ZRVHENI_MR, ZB_TH, ZB_RV, ZB_RI)
  PRVHENI_MR(:,:,:)=UNPACK(ZRVHENI_MR(:), MASK=GNEGT(:,:,:), FIELD=0.0)
  PCIT(:,:,:)=UNPACK(ZCIT(:), MASK=GNEGT(:,:,:), FIELD=PCIT(:,:,:))
END IF
!

CONTAINS
  FUNCTION ICE4_NUCLEATION_COUNTJV(LTAB,KIT,KJT,KKT,KSIZE,I1,I2,I3) RESULT(IC)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: KIT, KJT, KKT, KSIZE
  LOGICAL, DIMENSION(KIT,KJT,KKT), INTENT(IN) :: LTAB ! Mask
  INTEGER, DIMENSION(KSIZE), INTENT(OUT) :: I1, I2, I3 ! Used to replace the COUNT and PACK
  INTEGER :: IC
  INTEGER :: JI, JJ, JK
  IC=0
  DO JK=1, SIZE(LTAB,3)
    DO JJ=1, SIZE(LTAB,2)
      DO JI=1, SIZE(LTAB,1)
        IF(LTAB(JI,JJ,JK)) THEN
          IC=IC+1
          I1(IC)=JI
          I2(IC)=JJ
          I3(IC)=JK
        END IF
      END DO
    END DO
  END DO
  END FUNCTION ICE4_NUCLEATION_COUNTJV
  !
END SUBROUTINE ICE4_NUCLEATION_WRAPPER
