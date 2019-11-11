!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/mesonh/sources/dataio/fmf90_cray/s.fmread.f90, Version:1.2.1.2, Date:98/09/16, Last modified:98/06/04
!-----------------------------------------------------------------
!##################
MODULE MODI_FMREAD
!##################
!
INTERFACE FMREAD
      SUBROUTINE FMREADX0(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, &
                           INTENT(OUT)::PFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMREADX0
!
      SUBROUTINE FMREADX1(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:), &
                           INTENT(OUT)::PFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMREADX1
!
!
      SUBROUTINE FMREADX2(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:), &
                           INTENT(OUT)::PFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMREADX2
!
!
      SUBROUTINE FMREADX3(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:,:), &
                           INTENT(OUT)::PFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMREADX3
!
!
      SUBROUTINE FMREADX4(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:,:,:), &
                           INTENT(OUT)::PFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMREADX4
!
!
      SUBROUTINE FMREADX5(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:,:,:,:), &
                           INTENT(OUT)::PFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMREADX5
!
!
      SUBROUTINE FMREADX6(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:,:,:,:,:), &
                           INTENT(OUT)::PFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMREADX6
!
      SUBROUTINE FMREADN0(HFILEM,HRECFM,HFIPRI,KLENG,KFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
INTEGER, &
                           INTENT(OUT)::KFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMREADN0
!
      SUBROUTINE FMREADN1(HFILEM,HRECFM,HFIPRI,KLENG,KFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
INTEGER, DIMENSION(:), &
                           INTENT(OUT)::KFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMREADN1
!
      SUBROUTINE FMREADN2(HFILEM,HRECFM,HFIPRI,KLENG,KFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
INTEGER, DIMENSION(:,:), &
                           INTENT(OUT)::KFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMREADN2
!
      SUBROUTINE FMREADL0(HFILEM,HRECFM,HFIPRI,KLENG,OFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
LOGICAL, &
                           INTENT(OUT)::OFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMREADL0
!
      SUBROUTINE FMREADL1(HFILEM,HRECFM,HFIPRI,KLENG,OFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
LOGICAL, DIMENSION(:), &
                           INTENT(OUT)::OFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMREADL1
!
      SUBROUTINE FMREADC0(HFILEM,HRECFM,HFIPRI,KLENG,HFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
CHARACTER(LEN=*), &
                           INTENT(OUT)::HFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMREADC0
!
      SUBROUTINE FMREADT0(HFILEM,HRECFM,HFIPRI,KLENG,TFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
USE MODD_TYPE_DATE
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
TYPE (DATE_TIME), &
                           INTENT(OUT)::TFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMREADT0
!
END INTERFACE
!
END MODULE MODI_FMREAD
!     #############################################################
      SUBROUTINE FMREADX0(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMREADX0* - routine to read a real scalar into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADX0 is to convert the real into integer(kind=8)
!     by calling FM_READ without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     06/08/97
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, &
                           INTENT(OUT)::PFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*),          INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables

CHARACTER(LEN=JPXKRK) ::YCOMMENT 
REAL(KIND=8) :: ZFIELD
!
!-------------------------------------------------------------------------------
CALL FM_READ(HFILEM,HRECFM,HFIPRI,1,ZFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
IF(KRESP==0) PFIELD = ZFIELD
IF(KRESP==0) HCOMMENT=YCOMMENT
!-------------------------------------------------------------------------------
END SUBROUTINE FMREADX0
!     #############################################################
      SUBROUTINE FMREADX1(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMREADX1* - routine to read a real 1D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADX0 is to convert the real into integer(kind=8)
!     by calling FM_READ without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     06/08/97
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:), &
                           INTENT(OUT)::PFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPXKRK) ::YCOMMENT 
INTEGER :: ILENG
REAL(KIND=8),DIMENSION(SIZE(PFIELD)) :: ZFIELD
!-------------------------------------------------------------------------------
!
ILENG=SIZE(PFIELD)
CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
IF(KRESP==0) PFIELD = ZFIELD
IF(KRESP==0) HCOMMENT = YCOMMENT
!-------------------------------------------------------------------------------
END SUBROUTINE FMREADX1
!     #############################################################
      SUBROUTINE FMREADX2(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMREADX2* - routine to read a real 2D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADX0 is to convert the real into integer(kind=8)
!     by calling FM_READ without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     06/08/97
!!      Modification 15/10/97 (V.Masson)    1D and 2D cases
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
USE MODD_CONF
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:), &
                           INTENT(OUT)::PFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPXKRK) ::YCOMMENT 
INTEGER :: ILENG
REAL(KIND=8),DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2)) :: ZFIELD
!-------------------------------------------------------------------------------
IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/9
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(2,2),KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) PFIELD(:,:)=SPREAD(SPREAD(ZFIELD(2,2),DIM=1,NCOPIES=3),DIM=2,NCOPIES=3)
ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/3
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(:,2),KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) PFIELD(:,:)=SPREAD(ZFIELD(:,2),DIM=2,NCOPIES=3)
ELSE
  ILENG=SIZE(PFIELD)
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) PFIELD = ZFIELD
END IF
IF(KRESP==0) HCOMMENT = YCOMMENT
!-------------------------------------------------------------------------------
END SUBROUTINE FMREADX2
!     #############################################################
      SUBROUTINE FMREADX3(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMREADX1* - routine to read a real 3D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADX0 is to convert the real into integer(kind=8)
!     by calling FM_READ without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     06/08/97
!!      Modification 15/10/97 (V.Masson)    1D and 2D cases
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
USE MODD_CONF
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:,:), &
                           INTENT(OUT)::PFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPXKRK) ::YCOMMENT 
INTEGER :: ILENG
REAL(KIND=8),DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2),SIZE(PFIELD,3)) :: ZFIELD
!-------------------------------------------------------------------------------
IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/9
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(2,2,:),KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) PFIELD(:,:,:)=SPREAD(SPREAD(ZFIELD(2,2,:),DIM=1,NCOPIES=3),DIM=2,NCOPIES=3)
ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/3
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(:,2,:),KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) PFIELD(:,:,:)=SPREAD(ZFIELD(:,2,:),DIM=2,NCOPIES=3)
ELSE
  ILENG=SIZE(PFIELD)
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) PFIELD = ZFIELD
END IF
IF(KRESP==0) HCOMMENT = YCOMMENT
!-------------------------------------------------------------------------------
END SUBROUTINE FMREADX3
!     #############################################################
      SUBROUTINE FMREADX4(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMREADX4* - routine to read a real 4D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADX0 is to convert the real into integer(kind=8)
!     by calling FM_READ without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     06/08/97
!!      Modification 15/10/97 (V.Masson)    1D and 2D cases
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
USE MODD_CONF
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:,:,:), &
                           INTENT(OUT)::PFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER :: ILENG
CHARACTER(LEN=JPXKRK) ::YCOMMENT 
REAL(KIND=8),DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2),       &
	      SIZE(PFIELD,3),SIZE(PFIELD,4)) :: ZFIELD
!-------------------------------------------------------------------------------
IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/9
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(2,2,:,:),KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) PFIELD(:,:,:,:)=SPREAD(SPREAD(ZFIELD(2,2,:,:),DIM=1,NCOPIES=3),DIM=2,NCOPIES=3)
ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/3
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(:,2,:,:),KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) PFIELD(:,:,:,:)=SPREAD(ZFIELD(:,2,:,:),DIM=2,NCOPIES=3)
ELSE
  ILENG=SIZE(PFIELD)
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) PFIELD = ZFIELD
END IF
IF(KRESP==0) HCOMMENT = YCOMMENT
!-------------------------------------------------------------------------------
END SUBROUTINE FMREADX4
!     #############################################################
      SUBROUTINE FMREADX5(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMREADX5* - routine to read a real 5D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADX0 is to convert the real into integer(kind=8)
!     by calling FM_READ without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     06/08/97
!!      Modification 15/10/97 (V.Masson)    1D and 2D cases
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
USE MODD_CONF
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:,:,:,:), &
                           INTENT(OUT)::PFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER :: ILENG
CHARACTER(LEN=JPXKRK) ::YCOMMENT 
REAL(KIND=8),DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2),       &
    SIZE(PFIELD,3),SIZE(PFIELD,4),SIZE(PFIELD,5)) :: ZFIELD
!-------------------------------------------------------------------------------
IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/9
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(2,2,:,:,:),KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) PFIELD(:,:,:,:,:)=SPREAD(SPREAD(ZFIELD(2,2,:,:,:),DIM=1,NCOPIES=3),DIM=2,NCOPIES=3)
ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/3
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(:,2,:,:,:),KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) PFIELD(:,:,:,:,:)=SPREAD(ZFIELD(:,2,:,:,:),DIM=2,NCOPIES=3)
ELSE
  ILENG=SIZE(PFIELD)
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) PFIELD = ZFIELD
END IF
IF(KRESP==0) HCOMMENT = YCOMMENT
!-------------------------------------------------------------------------------
END SUBROUTINE FMREADX5
!     #############################################################
      SUBROUTINE FMREADX6(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMREADX6* - routine to read a real 6D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADX0 is to convert the real into integer(kind=8)
!     by calling FM_READ without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     06/08/97
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:,:,:,:,:), &
                           INTENT(OUT)::PFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER :: ILENG
CHARACTER(LEN=JPXKRK) ::YCOMMENT 
REAL(KIND=8),DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2),                     &
    SIZE(PFIELD,3),SIZE(PFIELD,4),SIZE(PFIELD,5),SIZE(PFIELD,6)) :: ZFIELD
!-------------------------------------------------------------------------------
ILENG=SIZE(PFIELD)
CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
IF(KRESP==0) PFIELD = ZFIELD
IF(KRESP==0) HCOMMENT = YCOMMENT
!-------------------------------------------------------------------------------
END SUBROUTINE FMREADX6
!     #############################################################
      SUBROUTINE FMREADN0(HFILEM,HRECFM,HFIPRI,KLENG,KFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMREADN0* - routine to read a integer scalar into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADN0 is to convert the integer into integer(kind=8)
!     by calling FM_READ without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     06/08/97
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
INTEGER, &
                           INTENT(OUT)::KFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPXKRK) ::YCOMMENT 
INTEGER(KIND=8) :: IFIELD
!-------------------------------------------------------------------------------
CALL FM_READ(HFILEM,HRECFM,HFIPRI,1,IFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
IF(KRESP==0) KFIELD = IFIELD
IF(KRESP==0) HCOMMENT = YCOMMENT
!-------------------------------------------------------------------------------
END SUBROUTINE FMREADN0
!     #############################################################
      SUBROUTINE FMREADN1(HFILEM,HRECFM,HFIPRI,KLENG,KFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMREADN1* - routine to read a integer 1D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADN1 is to convert the integer into integer(kind=8)
!     by calling FM_READ without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     06/08/97
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
INTEGER, DIMENSION(:), &
                           INTENT(OUT)::KFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPXKRK) ::YCOMMENT 
INTEGER(KIND=8), DIMENSION(SIZE(KFIELD)) :: IFIELD
INTEGER                                  :: ILENG
!-------------------------------------------------------------------------------
ILENG=SIZE(KFIELD)
CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,IFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
IF(KRESP==0) KFIELD(:)=IFIELD(:)
IF(KRESP==0) HCOMMENT = YCOMMENT
!-------------------------------------------------------------------------------
END SUBROUTINE FMREADN1
!     #############################################################
      SUBROUTINE FMREADN2(HFILEM,HRECFM,HFIPRI,KLENG,KFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMREADN2* - routine to read a integer 2D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADN1 is to convert the integer into integer(kind=8)
!     by calling FM_READ without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     06/08/97
!!      Modification 15/10/97 (V.Masson)    1D and 2D cases
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
USE MODD_CONF
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
INTEGER, DIMENSION(:,:), &
                           INTENT(OUT)::KFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPXKRK) ::YCOMMENT 
INTEGER(KIND=8), DIMENSION(SIZE(KFIELD,1),SIZE(KFIELD,2)) :: IFIELD
INTEGER                                                   :: ILENG
!-------------------------------------------------------------------------------
IF (LPACK .AND. L1D .AND. SIZE(KFIELD,1)==3 .AND. SIZE(KFIELD,2)==3) THEN
  ILENG=SIZE(KFIELD)/9
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,IFIELD(2,2),KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) KFIELD(:,:)=SPREAD(SPREAD(IFIELD(2,2),DIM=1,NCOPIES=3),DIM=2,NCOPIES=3)
ELSE IF (LPACK .AND. L2D .AND. SIZE(KFIELD,2)==3) THEN
  ILENG=SIZE(KFIELD)/3
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,IFIELD(:,2),KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) KFIELD(:,:)=SPREAD(IFIELD(:,2),DIM=2,NCOPIES=3)
ELSE
  ILENG=SIZE(KFIELD)
  CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,IFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
  IF(KRESP==0) KFIELD(:,:)=IFIELD(:,:)
END IF
IF(KRESP==0) HCOMMENT = YCOMMENT
!-------------------------------------------------------------------------------
END SUBROUTINE FMREADN2
!     #############################################################
      SUBROUTINE FMREADL0(HFILEM,HRECFM,HFIPRI,KLENG,OFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMREADL0* - routine to read a logical scalar into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADN0 is to convert the integer into integer(kind=8)
!     by calling FM_READ without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     06/08/97
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
LOGICAL, &
                           INTENT(OUT)::OFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPXKRK) ::YCOMMENT 
INTEGER(KIND=8) :: IFIELD
!-------------------------------------------------------------------------------
!
CALL FM_READ(HFILEM,HRECFM,HFIPRI,1,IFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
!
IF(KRESP==0) THEN
  IF (IFIELD==1) THEN
    OFIELD=.TRUE.
  ELSE
    OFIELD=.FALSE.
  END IF
  HCOMMENT = YCOMMENT
END IF
!-------------------------------------------------------------------------------
END SUBROUTINE FMREADL0
!     #############################################################
      SUBROUTINE FMREADL1(HFILEM,HRECFM,HFIPRI,KLENG,OFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMREADL1* - routine to read a logical array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADN0 is to convert the integer into integer(kind=8)
!     by calling FM_READ without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     06/08/97
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
LOGICAL, DIMENSION(:), &
                           INTENT(OUT)::OFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPXKRK) ::YCOMMENT 
INTEGER(KIND=8), DIMENSION(SIZE(OFIELD)) :: IFIELD
!-------------------------------------------------------------------------------
!
CALL FM_READ(HFILEM,HRECFM,HFIPRI,SIZE(IFIELD),IFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
!
IF(KRESP==0) THEN
  WHERE (IFIELD==1)
    OFIELD=.TRUE.
  ELSEWHERE
    OFIELD=.FALSE.
  END WHERE
  HCOMMENT = YCOMMENT
END IF
!-------------------------------------------------------------------------------
END SUBROUTINE FMREADL1
!     #############################################################
      SUBROUTINE FMREADC0(HFILEM,HRECFM,HFIPRI,KLENG,HFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMREADL1* - routine to read a logical scalar into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADL0 is to convert the string into arrayr of
!      integer(kind=8) and to call FM_READ without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     06/08/97
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
CHARACTER(LEN=*), &
                           INTENT(OUT)::HFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER                                      :: JLOOP
CHARACTER(LEN=JPXKRK)                        ::YCOMMENT 
INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE   :: IFIELD
INTEGER                                      :: ILENG
!-------------------------------------------------------------------------------
!
ILENG=LEN(HFIELD)
ALLOCATE(IFIELD(ILENG))
!
CALL FM_READ(HFILEM,HRECFM,HFIPRI,ILENG,IFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
!
IF(KRESP==0) THEN
  DO JLOOP=1,ILENG
   HFIELD(JLOOP:JLOOP)=ACHAR(IFIELD(JLOOP))
  END DO
  HCOMMENT = YCOMMENT
END IF
!
DEALLOCATE(IFIELD)
!-------------------------------------------------------------------------------
END SUBROUTINE FMREADC0
!     #############################################################
      SUBROUTINE FMREADT0(HFILEM,HRECFM,HFIPRI,KLENG,TFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMREADT0* - routine to read a date_time scalar into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADT0 is to call FM_READ without interface module
!      and to retrieve the date_time information
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     18/08/97
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_FMDECLAR
USE MODD_TYPE_DATE
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
TYPE (DATE_TIME), &
                           INTENT(OUT)::TFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)          ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=16)              :: YRECFM    ! Name of the article to be read
CHARACTER(LEN=JPXKRK)          :: YCOMMENT 
INTEGER(KIND=8), DIMENSION(3)  :: ITDATE
REAL(KIND=8)                   :: ZFIELD
!-------------------------------------------------------------------------------
!
YRECFM=TRIM(HRECFM)//'%TDATE'
CALL FM_READ(HFILEM,YRECFM,HFIPRI,3,ITDATE,KGRID,KLENCH,YCOMMENT,KRESP)
TFIELD%TDATE=DATE(ITDATE(1),ITDATE(2),ITDATE(3))  
HCOMMENT = YCOMMENT
!
YRECFM=TRIM(HRECFM)//'%TIME'
CALL FM_READ(HFILEM,YRECFM,HFIPRI,1,ZFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
TFIELD%TIME=ZFIELD
HCOMMENT = YCOMMENT
!-------------------------------------------------------------------------------
END SUBROUTINE FMREADT0
