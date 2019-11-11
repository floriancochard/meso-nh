!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_EXTRACT_VORTEX
!     ##########################
INTERFACE
!
      SUBROUTINE EXTRACT_VORTEX(PVARIN,PR0,PVARMOYR0,PICEN,PJCEN,PVARVORT)
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PVARIN ! 3-D array of input field
REAL, DIMENSION(:,:), INTENT(IN)    :: PR0 ! 2-D array of the radius of the
                                            ! filtering domain based on the
                                            ! tangential wind criterion
REAL, DIMENSION(:), INTENT(IN)      :: PVARMOYR0 ! 1-D array of the average
                                                 !along the periphery of the
                                                 !filter at the radius r0(p)
INTEGER, DIMENSION(:), INTENT(IN)   :: PICEN ! center of the vortex 
INTEGER, DIMENSION(:), INTENT(IN)   :: PJCEN ! at level pk
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PVARVORT ! 3-D array of field extracted
!
END SUBROUTINE EXTRACT_VORTEX
!
END INTERFACE
!
END MODULE MODI_EXTRACT_VORTEX
!
!
!
!     ####################################################################
      SUBROUTINE EXTRACT_VORTEX(PVARIN,PR0,PVARMOYR0,PICEN,PJCEN,PVARVORT)
!     ####################################################################
!
!!****  *EXTRACT_VORTEX* - extraction of a tropical cyclone
!!
!!    PURPOSE
!!    -------
!         This subroutine isolates an analyzed hurricane vortex hav
!           from the disturbance field hd. This is accomplished with
!         the following cylindrical filter centered at the location
!         of the analyzed tropical cyclone:
!
!            hav(i,j,k) = hd(i,j,k) - {hd(r0,phi)E(r) + hdbar(r0)[1 - E(r)]}
!
!             where r0(phi,z) represents the radius of the filtering domain
!
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
!!      Module MODD_PARAMETERS : contains parameters
!!        JPHEXT : number of horizontal external points
!!        JPVEXT : number of vertical external points
!!      Module MODD_GRID1 : contains grid variables
!!        XXHAT :  Position x in the conformal
!!                 plane or on the cartesian plane
!!        XYHAT :  Position y in the conformal
!!                 plane or on the cartesian plane
!!        XZHAT :  Position z in the grid without orography
!!      Module MODD_CST : contains physical constants
!!        XPI   :  pi 
!!      Module MODD_CONF   : model configuration
!!        L2D   : switch for 2D configuration
!!      Module MODD_DIM1   : model dimensions
!!        NIMAX,NJMAX
!!      Module MODD_LBC1   : model lateral boudary conditions
!!        CLBCX(2),CLBCY(2)
!!
!!    REFERENCE
!!    ---------
!!
!!     **??** subroutine EXTRACT_VORTEX
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
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY: XPI
USE MODD_PARAMETERS, ONLY: XUNDEF
USE MODD_DIM_n,      ONLY: NIMAX,NJMAX
USE MODD_GRID_n,     ONLY: XXHAT,XYHAT
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PVARIN ! 3-D array of input field
REAL, DIMENSION(:,:), INTENT(IN)    :: PR0 ! 2-D array of the radius of the
                                            ! filtering domain based on the
                                            ! tangential wind criterion
REAL, DIMENSION(:), INTENT(IN)      :: PVARMOYR0 ! 1-D array of the average
                                                 !along the periphery of the
                                                 !filter at the radius r0(p)
INTEGER, DIMENSION(:), INTENT(IN)   :: PICEN ! center of the vortex 
INTEGER, DIMENSION(:), INTENT(IN)   :: PJCEN ! at level pk
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PVARVORT ! 3-D array of field extracted
!
!*       0.2   Declarations of local variables
!
REAL                                          :: ZDELTAX,ZDELTAY,ZDELTAR
REAL                                          :: ZXK,ZYK
REAL                                          :: ZR,ZPHI,ZDPHI,ZDELT
REAL                                          :: ZR0PROJ,ZVARINPROJ,ZEFUNC,ZLTERM
REAL                                          :: ZXRK,ZYRK,ZXI0,ZYJ0
REAL                                          :: ZX00,ZY00
INTEGER                                       :: JI,JJ,JK
INTEGER                                       :: IXRK,IYRK,IKPHI,IKPHI1
INTEGER                                       :: IX,IY,IP,IPHI
INTEGER                                       :: IIB,IIE,IJB,IJE
!
!-------------------------------------------------------------------------------
!
IP= SIZE(PVARIN,3)
IX= SIZE(PVARIN,1)
IY= SIZE(PVARIN,2)
IPHI=SIZE(PR0,1)
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
ZDELTAX = XXHAT(3) - XXHAT(2)
ZDELTAY = XYHAT(3) - XYHAT(2)
ZDELTAR = MAX(ZDELTAX,ZDELTAY)
ZDPHI = 2.*XPI / IPHI
!
DO JK = 1, IP
  ZXI0 = XXHAT(PICEN(JK)) + (ZDELTAX / 2.)
  ZYJ0 = XYHAT(PJCEN(JK)) + (ZDELTAY / 2.)
  ZX00 = XXHAT(1) + (ZDELTAX / 2.)
  ZY00 = XYHAT(1) + (ZDELTAY / 2.)
  DO JJ = 1, IY
    DO JI = 1, IX
      ZXK = XXHAT(JI) - XXHAT(PICEN(JK))
      ZYK = XYHAT(JJ) - XYHAT(PJCEN(JK))
      ZR = SQRT(ZXK*ZXK + ZYK*ZYK)
      IF (ZR == 0. .OR. ZR < ZDELTAR) THEN
        PVARVORT(JI,JJ,JK) = PVARIN(JI,JJ,JK) - PVARMOYR0(JK)
      ELSE
        ZPHI = ATAN2(ZYK,ZXK)
        IF (ZPHI < 0.) ZPHI = ZPHI + 2.*XPI
        IKPHI = ZPHI / ZDPHI + 1
        IF (IKPHI>IPHI) IKPHI=IKPHI-IPHI
        ZDELT = ZPHI - (IKPHI - 1) * ZDPHI
        IKPHI1=IKPHI+1
        IF (IKPHI1>IPHI) IKPHI1=IKPHI1-IPHI
        ZR0PROJ = PR0(IKPHI,JK) + (PR0(IKPHI1,JK)-PR0(IKPHI,JK))* ZDELT / ZDPHI
        IF (ZR <= ZR0PROJ) THEN
          ZXRK = ZR0PROJ * COS(ZPHI) + ZXI0
          ZYRK = ZR0PROJ * SIN(ZPHI) + ZYJ0
          IXRK = (ZXRK - ZX00) / ZDELTAX + 1
          IYRK = (ZYRK - ZY00) / ZDELTAY + 1
          IF (IXRK >= IIB .AND. (IXRK+1) <= IIE      &
              .AND. IYRK >= IJB .AND. (IYRK+1) <= IJE) THEN 
            ZXRK = (ZXRK - XXHAT(IXRK)) / ZDELTAX - 0.5
            ZYRK = (ZYRK - XYHAT(IYRK)) / ZDELTAY - 0.5
            ZVARINPROJ = PVARIN(IXRK,IYRK,JK)*(1-ZXRK)*(1-ZYRK)            &
                       + PVARIN(IXRK+1,IYRK,JK)*ZXRK*(1-ZYRK)              &
                       + PVARIN(IXRK,IYRK+1,JK)*(1-ZXRK)*ZYRK              &
                       + PVARIN(IXRK+1,IYRK+1,JK)*ZXRK*ZYRK 
            ZLTERM = 0.2 * ZR0PROJ
            ZEFUNC = (EXP(-(ZR0PROJ - ZR)**2. / ZLTERM**2.) -               &
                      EXP(-(ZR0PROJ)**2. / ZLTERM**2.)      )   /            &
                     (1. - EXP(-(ZR0PROJ)**2. / ZLTERM**2.) )
            PVARVORT(JI,JJ,JK) = PVARIN(JI,JJ,JK) - (ZVARINPROJ * ZEFUNC +   &
                                                   PVARMOYR0(JK)*(1-ZEFUNC)  )
          ELSE
            PVARVORT(JI,JJ,JK) = XUNDEF
          END IF
        ELSE
          PVARVORT(JI,JJ,JK) = XUNDEF
        END IF
      END IF
    END DO 
  END DO       
END DO
!
!-------------------------------------------------------------------------------
END SUBROUTINE EXTRACT_VORTEX
