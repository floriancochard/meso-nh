!MNH_LIC Copyright 2001-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_POLAR_CALC
!     ####################
INTERFACE POLAR_CALC
!
      SUBROUTINE POLAR_CALC3D(PVARIN,KICEN,KJCEN,PVARCYL,PVAR2IN,KRADMAX)
!
REAL, DIMENSION(:,:,:), INTENT(IN)        :: PVARIN ! 3-D array of input field 
INTEGER, DIMENSION(:), INTENT(IN)     :: KICEN ! center of the vortex 
INTEGER, DIMENSION(:), INTENT(IN)     :: KJCEN !along the vertical
REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PVARCYL ! 3-D array of output field
REAL, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PVAR2IN ! 2nd input field 
INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: KRADMAX ! max. radius in
                                                        !each azimuthal direction
!
END SUBROUTINE POLAR_CALC3D
!
      SUBROUTINE POLAR_CALC2D(PVARIN,KICEN,KJCEN,PVARCYL,PVAR2IN,KRADMAX)
!
REAL, DIMENSION(:,:),  INTENT(IN)     :: PVARIN ! 2-D array of input field 
INTEGER,               INTENT(IN)     :: KICEN ! center of the vortex 
INTEGER,               INTENT(IN)     :: KJCEN !
REAL, DIMENSION(:,:),  INTENT(OUT)    :: PVARCYL ! 2-D array of output field
REAL, DIMENSION(:,:),  INTENT(IN), OPTIONAL :: PVAR2IN ! 2nd input field 
INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: KRADMAX ! max. radius in
                                                        !each azimuthal direction
!
END SUBROUTINE POLAR_CALC2D
!
END INTERFACE POLAR_CALC
END MODULE MODI_POLAR_CALC
!
!     #################################################
      SUBROUTINE POLAR_CALC3D(PVARIN,KICEN,KJCEN,PVARCYL,PVAR2IN,KRADMAX)
!     #################################################
!
!!****  *POLAR_CALC* - interpolation onto cylindrical grid
!!
!!    PURPOSE
!!    -------
!         This subroutine interpolates the input field onto a cylindrical grid
!     centered on (xcen,ycen) at level pk.
!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!  	O. Nuissier           * L.A. *
!!      R. Rogers             * NOAA/AOML/HRD *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original              01/12/01
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF,   ONLY: NVERB
USE MODD_CST,    ONLY: XPI
USE MODD_DIM_n,  ONLY: NIMAX,NJMAX
USE MODD_GRID_n, ONLY: XXHAT,XYHAT
USE MODD_LUNIT,  ONLY: TLUOUT0
!
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:,:,:),  INTENT(IN)   :: PVARIN ! 3-D array of input field 
INTEGER, DIMENSION(:), INTENT(IN)     :: KICEN ! center of the vortex 
INTEGER, DIMENSION(:), INTENT(IN)     :: KJCEN !along the vertical
REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PVARCYL ! 3-D array of output field
REAL, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PVAR2IN ! 2nd input field 
INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: KRADMAX ! max. radius in
                                                        !each azimuthal direction
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PVARCYL,1),SIZE(PVARCYL,2),SIZE(PVARCYL,3)) :: ZVARCYL1
REAL, DIMENSION(:,:,:), ALLOCATABLE       :: ZVARCYL2
REAL, DIMENSION(SIZE(PVARCYL,1))          :: ZRADIUS
REAL                                      :: ZXI0,ZYJ0,ZX00,ZY00,ZXK,ZYK
REAL                                      :: ZDELTAX,ZDELTAY,ZDELTAR,ZDPHI,ZPHI
INTEGER                                   :: ILUOUT0
INTEGER                                   :: IIX,IIY
INTEGER                                   :: IR,IPHI,IP ! arrays dimensions
INTEGER                                   :: IIB,IJB,IIE,IJE
INTEGER                                   :: JR,JPHI,JP
INTEGER,DIMENSION(SIZE(PVARCYL,2))        :: IPBL
!
!-------------------------------------------------------------------------------
!
!*	 1.     INITIALIZATIONS
!               ---------------
!
ILUOUT0 = TLUOUT0%NLU
!
IR=SIZE(PVARCYL,1)
IPHI=SIZE(PVARCYL,2)
IP=SIZE(PVARCYL,3)
!
ZDELTAX = XXHAT(3) - XXHAT(2)
ZDELTAY = XYHAT(3) - XYHAT(2)
ZDELTAR = MAX(ZDELTAX,ZDELTAY)
ZDPHI = 2.*XPI/IPHI
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IF (NVERB>=5) &
  WRITE(ILUOUT0,'(A,4I4)')'routine polar_calc: indexes of MesoNH domain ', &
                          IIB,IIE,IJB,IJE
!    
DO JR = 1, IR
  ZRADIUS(JR) = JR * ZDELTAR
END DO
IPBL(:)=0
!
IF (PRESENT(PVAR2IN)) ALLOCATE(ZVARCYL2(IR,IPHI,IP))
!-----------------------------------------------------------------------------
!
!*	 2.     INTERPOLATE ONTO A CYLINDRICAL GRID
!           (R=1 CORRESPOND TO THE CENTER OF THE VORTEX) 
!           --------------------------------------------
!
DO JP = 1, IP
  ZXI0 = XXHAT(KICEN(JP)) + (ZDELTAX / 2.)
  ZYJ0 = XYHAT(KJCEN(JP)) + (ZDELTAY / 2.)
  ZX00 = XXHAT(1) + (ZDELTAX / 2.)
  ZY00 = XYHAT(1) + (ZDELTAY / 2.)
  DO JR = 1, IR
    DO JPHI = 1, IPHI
      ZPHI = (JPHI - 1) * ZDPHI
      ZXK = ZRADIUS(JR) * COS(ZPHI) + ZXI0
      ZYK = ZRADIUS(JR) * SIN(ZPHI) + ZYJ0
      IIX = (ZXK - ZX00) / ZDELTAX + 1
      IIY = (ZYK - ZY00) / ZDELTAY + 1
      !
      IF (IIX >IIB .AND. (IIX+1) <IIE       &
          .AND. IIY >IJB .AND. (IIY+1) <IJE )   THEN
        !
        KRADMAX(JPHI)=JR
        ZXK = (ZXK - XXHAT(IIX)) / ZDELTAX - 0.5
        ZYK = (ZYK - XYHAT(IIY)) / ZDELTAY - 0.5
        ZVARCYL1(JR,JPHI,JP) = PVARIN(IIX,IIY,JP)*(1-ZXK)*(1-ZYK) +      &
                               PVARIN(IIX+1,IIY,JP)*ZXK*(1-ZYK)   +      &
                               PVARIN(IIX,IIY+1,JP)*(1-ZXK)*ZYK   +      &
                               PVARIN(IIX+1,IIY+1,JP)*ZXK*ZYK
        IF (PRESENT(PVAR2IN)) THEN
          ZVARCYL2(JR,JPHI,JP) = PVAR2IN(IIX,IIY,JP)*(1-ZXK)*(1-ZYK) +      &
                                 PVAR2IN(IIX+1,IIY,JP)*ZXK*(1-ZYK)   +      &
                                 PVAR2IN(IIX,IIY+1,JP)*(1-ZXK)*ZYK   +      &
                                 PVAR2IN(IIX+1,IIY+1,JP)*ZXK*ZYK
          ! compute tangential wind if the two components are present
          PVARCYL(JR,JPHI,JP) = ZVARCYL2(JR,JPHI,JP)*COS(ZPHI)  &
                               -ZVARCYL1(JR,JPHI,JP)*SIN(ZPHI)
          ! radial component is = ZVARCYL1(:,JPHI,JP)*COS(ZPHI)  &
          !                      +ZVARCYL2(:,JPHI,JP)*SIN(ZPHI)
        ELSE
          PVARCYL(JR,JPHI,JP) = ZVARCYL1(JR,JPHI,JP)
        END IF
      ELSE
        PVARCYL(JR,JPHI,JP) = 0.
        IF (IPBL(JPHI)==0) THEN
          WRITE(ILUOUT0,'(A,I2,A,F5.1,A)')' -> Warning in polar_calc: for the azimuthal direction ',JPHI,' (',ZPHI*180./XPI,'dg)'
          WRITE(ILUOUT0,'(A,I3,A,F6.1,A)') &
             '     radial dimension (from XRADGUESS) is too big: ',JR, &
             ' grid points (',ZRADIUS(JR)/1000.,'km)'
          WRITE(ILUOUT0,'(A,I3,A,I3)') '    I=',IIX,' J=',IIY
          IPBL(JPHI)=1
        END IF
      END IF
    END DO
  END DO
END DO
!
!-------------------------------------------------------------------------------
END SUBROUTINE POLAR_CALC3D
!
!     #################################################
      SUBROUTINE POLAR_CALC2D(PVARIN,KICEN,KJCEN,PVARCYL,PVAR2IN,KRADMAX)
!     #################################################
!
!!****  *POLAR_CALC* - interpolation onto cylindrical grid
!!
!!    PURPOSE
!!    -------
!         This subroutine interpolates the input field onto a cylindrical grid
!     centered on (xcen,ycen) at level pk.
!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!  	O. Nuissier           * L.A. *
!!      R. Rogers             * NOAA/AOML/HRD *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original              01/12/01
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
INTERFACE
      SUBROUTINE POLAR_CALC3D(PVARIN,KICEN,KJCEN,PVARCYL,PVAR2IN,KRADMAX)
!
REAL, DIMENSION(:,:,:), INTENT(IN)        :: PVARIN ! 3-D array of input field 
INTEGER, DIMENSION(:), INTENT(IN)     :: KICEN ! center of the vortex 
INTEGER, DIMENSION(:), INTENT(IN)     :: KJCEN !along the vertical
REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PVARCYL ! 3-D array of output field
REAL, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: PVAR2IN ! 2nd input field 
INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: KRADMAX ! max. radius in
                                                        !each azimuthal direction
!
      END SUBROUTINE POLAR_CALC3D
END INTERFACE
!
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:,:),  INTENT(IN)     :: PVARIN ! 2-D array of input field 
INTEGER,               INTENT(IN)     :: KICEN ! center of the vortex 
INTEGER,               INTENT(IN)     :: KJCEN !
REAL, DIMENSION(:,:),  INTENT(OUT)    :: PVARCYL ! 2-D array of output field
REAL, DIMENSION(:,:),  INTENT(IN), OPTIONAL :: PVAR2IN ! 2nd input field 
INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL :: KRADMAX ! max. radius in
                                                        !each azimuthal direction
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PVARIN,1), SIZE(PVARIN,2),1)  :: ZVARIN3D
REAL, DIMENSION(:,:,:), ALLOCATABLE                :: ZVAR2IN3D
REAL, DIMENSION(SIZE(PVARCYL,1),SIZE(PVARCYL,2),1) :: ZVARCYL3D
INTEGER, DIMENSION(1) :: IICEN3D,IJCEN3D
!
!-------------------------------------------------------------------------------
!
!*	 1.     INITIALIZATIONS
!               ---------------
!
ZVARIN3D(:,:,1)=PVARIN(:,:)
IICEN3D(1) = KICEN
IJCEN3D(1) = KJCEN
!
IF (PRESENT(PVAR2IN)) THEN
  ALLOCATE(ZVAR2IN3D(SIZE(PVARIN,1),SIZE(PVARIN,2),1))
  ZVAR2IN3D(:,:,1)=PVAR2IN(:,:)
  IF (PRESENT(KRADMAX)) THEN
    CALL POLAR_CALC3D(ZVARIN3D,IICEN3D,IJCEN3D,ZVARCYL3D, &
                      PVAR2IN=ZVAR2IN3D,KRADMAX=KRADMAX)
  ELSE
    CALL POLAR_CALC3D(ZVARIN3D,IICEN3D,IJCEN3D,ZVARCYL3D, &
                      PVAR2IN=ZVAR2IN3D                   )
  END IF
ELSE
  IF (PRESENT(KRADMAX)) THEN
    CALL POLAR_CALC3D(ZVARIN3D,IICEN3D,IJCEN3D,ZVARCYL3D,KRADMAX=KRADMAX)
  ELSE
    CALL POLAR_CALC3D(ZVARIN3D,IICEN3D,IJCEN3D,ZVARCYL3D)
  END IF
END IF
!
PVARCYL(:,:) = ZVARCYL3D(:,:,1)
!-----------------------------------------------------------------------------
!
END SUBROUTINE POLAR_CALC2D
