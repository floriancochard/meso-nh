!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      ########################
       MODULE MODI_BHMIE_BHCOAT
!      ########################
!
INTERFACE
      SUBROUTINE BHMIE_BHCOAT(PSIZE_PARAM_CORE,PSIZE_PARAM_COAT, &
                              PPREFR_CORE,PPREFR_COAT,PQEXT,PQBAK)
!
REAL,                  INTENT(IN)  :: PSIZE_PARAM_CORE
REAL,                  INTENT(IN)  :: PSIZE_PARAM_COAT
COMPLEX,               INTENT(IN)  :: PPREFR_CORE
COMPLEX,               INTENT(IN)  :: PPREFR_COAT
REAL,                  INTENT(OUT) :: PQEXT,PQBAK
!
END SUBROUTINE BHMIE_BHCOAT
END INTERFACE
END MODULE MODI_BHMIE_BHCOAT 
!
!     ############################################################
      SUBROUTINE BHMIE_BHCOAT(PSIZE_PARAM_CORE,PSIZE_PARAM_COAT, &
                              PPREFR_CORE,PPREFR_COAT,PQEXT,PQBAK)
!     ############################################################
!
!!***********************************************************************
!!
!! Subroutine BHCOAT calculates Q_ext, Q_back for coated sphere.
!! All bessel functions computed by upward recurrence.
!! Input:
!!        X = 2*PI*RCORE*REFMED/WAVEL
!!        Y = 2*PI*RMANT*REFMED/WAVEL
!!        RFREL1 = REFCOR/REFMED
!!        RFREL2 = REFMAN/REFMED 
!! where  REFCOR = complex refr.index of core)
!!        REFMAN = complex refr.index of mantle)
!!        REFMED = real refr.index of medium)
!!        RCORE = radius of core
!!        RMANT = radius of mantle
!!        WAVEL = wavelength of light in ambient medium
!!
!! Routine BHCOAT is taken from Bohren & Huffman (1983)
!! Obtained from C.L.Joseph
!!
!! History:
!! 92/11/24 (BTD) Explicit declaration of all variables
!!***********************************************************************
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                  INTENT(IN)  :: PSIZE_PARAM_CORE
REAL,                  INTENT(IN)  :: PSIZE_PARAM_COAT
COMPLEX,               INTENT(IN)  :: PPREFR_CORE
COMPLEX,               INTENT(IN)  :: PPREFR_COAT
REAL,                  INTENT(OUT) :: PQEXT,PQBAK
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JJ
INTEGER :: ISTOP,IFLAG
REAL (KIND(0.0D0))    :: ZPSIY,ZCHIY
REAL (KIND(0.0D0))    :: ZPSI0Y,ZPSI1Y
REAL (KIND(0.0D0))    :: ZCHI0Y,ZCHI1Y
REAL (KIND(0.0D0))    :: ZEN,ZONE,ZSIZE_PARAM_STOP
REAL (KIND(0.0D0))    :: ZDEL=1.0E-8
COMPLEX (KIND(0.0D0)) :: ZZREFR_COAT_CORE
COMPLEX (KIND(0.0D0)) :: ZZCHI0Y2,ZZCHI1Y2,ZZCHI0X2,ZZCHI1X2
COMPLEX (KIND(0.0D0)) :: ZZCHIX2,ZZCHIPX2,ZZCHIY2,ZZCHIPY2
COMPLEX (KIND(0.0D0)) :: ZZAN,ZZBN,ZZXIY,ZZXI0Y,ZZXI1Y
COMPLEX (KIND(0.0D0)) :: ZZANCAP,ZZBNCAP
COMPLEX (KIND(0.0D0)) :: ZZX1,ZZX2,ZZY2
COMPLEX (KIND(0.0D0)) :: ZZD0X1,ZZD0X2,ZZD0Y2,ZZD1X1,ZZD1X2,ZZD1Y2
COMPLEX (KIND(0.0D0)) :: ZZDBAR,ZZGBAR,ZZBAK
COMPLEX (KIND(0.0D0)) :: ZZAMESS1,ZZAMESS2,ZZAMESS3,ZZAMESS4
COMPLEX (KIND(0.0D0)) :: ZZBRACK,ZZCRACK
!
!***********************************************************************
!
!* ZDEL is the inner sphere convergence criterion
!
ZZX1 = PSIZE_PARAM_CORE*PPREFR_CORE
ZZX2 = PSIZE_PARAM_CORE*PPREFR_COAT
ZZY2 = PSIZE_PARAM_COAT*PPREFR_COAT
ZSIZE_PARAM_STOP = PSIZE_PARAM_COAT + 4.*PSIZE_PARAM_COAT**0.3333 + 2.0
ISTOP = INT(ZSIZE_PARAM_STOP)
!
ZZREFR_COAT_CORE = PPREFR_COAT/PPREFR_CORE
!
!         -----------------------------------------------------------
!              series terminated after nstop terms
!         -----------------------------------------------------------
!
!* Suffix 1 means CORE medium property
!* Suffix 2 means COAT medium property
!
!
!* Suffix x means CORE size
!* Suffix y means COAT size
!
ZZD0X1 = COS(ZZX1)/SIN(ZZX1)
ZZD0X2 = COS(ZZX2)/SIN(ZZX2)
ZZD0Y2 = COS(ZZY2)/SIN(ZZY2)
!
ZPSI0Y = COS(PSIZE_PARAM_COAT)
ZPSI1Y = SIN(PSIZE_PARAM_COAT)
ZCHI0Y =-SIN(PSIZE_PARAM_COAT)
ZCHI1Y = COS(PSIZE_PARAM_COAT)
ZZXI0Y = CMPLX(ZPSI0Y,-ZCHI0Y)
ZZXI1Y = CMPLX(ZPSI1Y,-ZCHI1Y)
!
ZZCHI0Y2 =-SIN(ZZY2)
ZZCHI1Y2 = COS(ZZY2)
ZZCHI0X2 =-SIN(ZZX2)
ZZCHI1X2 = COS(ZZX2)
!
PQEXT = 0.0
ZZBAK = (0.0,0.0)
ZONE  = 1.0
IFLAG = 0
DO JJ = 1,ISTOP
  ZEN = FLOAT(JJ)
  ZPSIY = (2.0*ZEN-1.)*ZPSI1Y/PSIZE_PARAM_COAT - ZPSI0Y
  ZCHIY = (2.0*ZEN-1.)*ZCHI1Y/PSIZE_PARAM_COAT - ZCHI0Y
  ZZXIY = CMPLX(ZPSIY,-ZCHIY)
!
  ZZD1Y2 = 1.0/(ZEN/ZZY2-ZZD0Y2) - ZEN/ZZY2
!
  IF (IFLAG/=1) THEN
    ZZD1X1 = 1.0/(ZEN/ZZX1-ZZD0X1) - ZEN/ZZX1
    ZZD1X2 = 1.0/(ZEN/ZZX2-ZZD0X2) - ZEN/ZZX2
    ZZCHIX2 = (2.0*ZEN - 1.0)*ZZCHI1X2/ZZX2 - ZZCHI0X2
    ZZCHIY2 = (2.0*ZEN - 1.0)*ZZCHI1Y2/ZZY2 - ZZCHI0Y2
    ZZCHIPX2 = ZZCHI1X2 - ZEN*ZZCHIX2/ZZX2
    ZZCHIPY2 = ZZCHI1Y2 - ZEN*ZZCHIY2/ZZY2
    ZZANCAP = (ZZREFR_COAT_CORE*ZZD1X1         - ZZD1X2  )/ &
              (ZZREFR_COAT_CORE*ZZD1X1*ZZCHIX2 - ZZCHIPX2)/ &
                      (ZZCHIX2*ZZD1X2 - ZZCHIPX2)
    ZZBRACK = ZZANCAP*(ZZCHIY2*ZZD1Y2 - ZZCHIPY2)
    ZZBNCAP = (ZZREFR_COAT_CORE*ZZD1X2   - ZZD1X1        )/ &
              (ZZREFR_COAT_CORE*ZZCHIPX2 - ZZD1X1*ZZCHIX2)/ &
                      (ZZCHIX2*ZZD1X2 - ZZCHIPX2)
    ZZCRACK = ZZBNCAP*(ZZCHIY2*ZZD1Y2 - ZZCHIPY2)
!
    ZZAMESS1 = ZZBRACK*ZZCHIPY2
    ZZAMESS2 = ZZBRACK*ZZCHIY2
    ZZAMESS3 = ZZCRACK*ZZCHIPY2
    ZZAMESS4 = ZZCRACK*ZZCHIY2
!
    IF (ABS(ZZAMESS1)<ZDEL*ABS(ZZD1Y2) .OR. &
        ABS(ZZAMESS2)<ZDEL             .OR. &
        ABS(ZZAMESS3)<ZDEL*ABS(ZZD1Y2) .OR. &
        ABS(ZZAMESS4)<ZDEL                   ) THEN
      ZZBRACK = (0.,0.)
      ZZCRACK = (0.,0.)
      IFLAG   = 1
    END IF
  END IF
  ZZDBAR = (ZZD1Y2 - ZZBRACK*ZZCHIPY2)/(1.0-ZZBRACK*ZZCHIY2)
  ZZGBAR = (ZZD1Y2 - ZZCRACK*ZZCHIPY2)/(1.0-ZZCRACK*ZZCHIY2)
  ZZAN = ((ZZDBAR/PPREFR_COAT + ZEN/PSIZE_PARAM_COAT)*ZPSIY - ZPSI1Y) / &
         ((ZZDBAR/PPREFR_COAT + ZEN/PSIZE_PARAM_COAT)*ZZXIY - ZZXI1Y)
  ZZBN = ((PPREFR_COAT*ZZGBAR + ZEN/PSIZE_PARAM_COAT)*ZPSIY - ZPSI1Y) / &
         ((PPREFR_COAT*ZZGBAR + ZEN/PSIZE_PARAM_COAT)*ZZXIY - ZZXI1Y)
!
  ZONE  = -ZONE
  ZZBAK = ZZBAK + (2.0*ZEN + 1.0)*ZONE*(ZZAN-ZZBN)
  PQEXT = PQEXT + (2.0*ZEN + 1.0)*     (REAL(ZZAN)+REAL(ZZBN))
!
  ZPSI0Y = ZPSI1Y
  ZPSI1Y = ZPSIY
  ZCHI0Y = ZCHI1Y
  ZCHI1Y = ZCHIY
  ZZXI1Y = CMPLX(ZPSI1Y,-ZCHI1Y)
!
  ZZCHI0X2 = ZZCHI1X2
  ZZCHI1X2 = ZZCHIX2
  ZZCHI0Y2 = ZZCHI1Y2
  ZZCHI1Y2 = ZZCHIY2
!
  ZZD0X1 = ZZD1X1
  ZZD0X2 = ZZD1X2
  ZZD0Y2 = ZZD1Y2
END DO
!
PQEXT =  PQEXT*(2.0/(PSIZE_PARAM_COAT*PSIZE_PARAM_COAT))
PQBAK = (ABS(ZZBAK)/ PSIZE_PARAM_COAT)**2
RETURN
END SUBROUTINE BHMIE_BHCOAT
