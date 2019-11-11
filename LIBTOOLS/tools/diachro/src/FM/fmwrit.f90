!##################
MODULE MODI_FMWRIT
!##################
!
INTERFACE FMWRIT
      SUBROUTINE FMWRITX0(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, &
                           INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMWRITX0
!
      SUBROUTINE FMWRITX1(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:), &
                           INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMWRITX1
!
!
      SUBROUTINE FMWRITX2(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:), &
                           INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMWRITX2
!
!
      SUBROUTINE FMWRITX3(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:,:), &
                           INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMWRITX3
!
!
      SUBROUTINE FMWRITX4(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:,:,:), &
                           INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMWRITX4
!
!
      SUBROUTINE FMWRITX5(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:,:,:,:), &
                           INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMWRITX5
!
!
      SUBROUTINE FMWRITX6(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:,:,:,:,:), &
                           INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMWRITX6
!
      SUBROUTINE FMWRITN0(HFILEM,HRECFM,HFIPRI,KLENG,KFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
INTEGER, &
                           INTENT(IN) ::KFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMWRITN0
!
      SUBROUTINE FMWRITN1(HFILEM,HRECFM,HFIPRI,KLENG,KFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
INTEGER, DIMENSION(:), &
                           INTENT(IN) ::KFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMWRITN1
!
      SUBROUTINE FMWRITN2(HFILEM,HRECFM,HFIPRI,KLENG,KFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
INTEGER, DIMENSION(:,:), &
                           INTENT(IN) ::KFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMWRITN2
!
      SUBROUTINE FMWRITL0(HFILEM,HRECFM,HFIPRI,KLENG,OFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
LOGICAL, &
                           INTENT(IN) ::OFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMWRITL0
!
      SUBROUTINE FMWRITL1(HFILEM,HRECFM,HFIPRI,KLENG,OFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
LOGICAL,DIMENSION(:),  &
                           INTENT(IN) ::OFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMWRITL1
!
      SUBROUTINE FMWRITC0(HFILEM,HRECFM,HFIPRI,KLENG,HFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
CHARACTER(LEN=*), &
                           INTENT(IN) ::HFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMWRITC0
!
      SUBROUTINE FMWRITT0(HFILEM,HRECFM,HFIPRI,KLENG,TFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
USE MODD_TYPE_DATE
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
TYPE (DATE_TIME), &
                           INTENT(IN) ::TFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
END SUBROUTINE FMWRITT0
!
END INTERFACE
!
END MODULE MODI_FMWRIT
!     #############################################################
      SUBROUTINE FMWRITX0(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMWRITX0* - routine to write a real scalar into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITX0 is to convert the real into integer(kind=8)
!     by calling FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, &
                           INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string
!
CHARACTER(LEN=*),          INTENT(IN) ::HCOMMENT ! comment string
!
INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
REAL(KIND=8) :: ZFIELD
!
!-------------------------------------------------------------------------------
ZFIELD=PFIELD
CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,1,ZFIELD,KGRID,KLENCH,HCOMMENT,KRESP)
!-------------------------------------------------------------------------------
END SUBROUTINE FMWRITX0
!
!     #############################################################
      SUBROUTINE FMWRITX1(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMWRITX1* - routine to write a real 1D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITX0 is to convert the real into integer(kind=8)
!     by calling FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:), &
                           INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER :: ILENG
REAL(KIND=8), DIMENSION(SIZE(PFIELD)) :: ZFIELD
!-------------------------------------------------------------------------------
!
ILENG=SIZE(PFIELD)
ZFIELD=PFIELD
CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD,KGRID,KLENCH,HCOMMENT,KRESP)
!-------------------------------------------------------------------------------
END SUBROUTINE FMWRITX1
!
!     #############################################################
      SUBROUTINE FMWRITX2(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMWRITX2* - routine to write a real 2D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITX0 is to convert the real into integer(kind=8)
!     by calling FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
                           INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER :: ILENG
REAL(KIND=8), DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2)) :: ZFIELD
!-------------------------------------------------------------------------------
!
ZFIELD=PFIELD
IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/9
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(2:2,2:2),KGRID,KLENCH,HCOMMENT,KRESP)
ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/3
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(:,2),KGRID,KLENCH,HCOMMENT,KRESP)
ELSE
  ILENG=SIZE(PFIELD)
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD,KGRID,KLENCH,HCOMMENT,KRESP)
END IF
!-------------------------------------------------------------------------------
END SUBROUTINE FMWRITX2
!
!     #############################################################
      SUBROUTINE FMWRITX3(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMWRITX1* - routine to write a real 3D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITX0 is to convert the real into integer(kind=8)
!     by calling FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
                           INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER :: ILENG
REAL(KIND=8), DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2),SIZE(PFIELD,3)) :: ZFIELD
!-------------------------------------------------------------------------------
!
ZFIELD=PFIELD
IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/9
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(2,2,:),KGRID,KLENCH,HCOMMENT,KRESP)
ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/3
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(:,2,:),KGRID,KLENCH,HCOMMENT,KRESP)
ELSE
  ILENG=SIZE(PFIELD)
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD,KGRID,KLENCH,HCOMMENT,KRESP)
END IF
!-------------------------------------------------------------------------------
END SUBROUTINE FMWRITX3
!
!     #############################################################
      SUBROUTINE FMWRITX4(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMWRITX4* - routine to write a real 4D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITX0 is to convert the real into integer(kind=8)
!     by calling FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
                           INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER :: ILENG
REAL(KIND=8),    &
DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2),SIZE(PFIELD,3),SIZE(PFIELD,4)) :: ZFIELD
!-------------------------------------------------------------------------------
!
ZFIELD=PFIELD
IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/9
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(2,2,:,:),KGRID,KLENCH,HCOMMENT,KRESP)
ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/3
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(:,2,:,:),KGRID,KLENCH,HCOMMENT,KRESP)
ELSE
  ILENG=SIZE(PFIELD)
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD,KGRID,KLENCH,HCOMMENT,KRESP)
END IF
!-------------------------------------------------------------------------------
END SUBROUTINE FMWRITX4
!
!     #############################################################
      SUBROUTINE FMWRITX5(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMWRITX5* - routine to write a real 5D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITX0 is to convert the real into integer(kind=8)
!     by calling FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
                           INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER :: ILENG
REAL(KIND=8),    &
DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2),SIZE(PFIELD,3),SIZE(PFIELD,4),SIZE(PFIELD,5)) :: ZFIELD
!-------------------------------------------------------------------------------
!
ZFIELD=PFIELD
IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/9
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(2,2,:,:,:),KGRID,KLENCH,HCOMMENT,KRESP)
ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
  ILENG=SIZE(PFIELD)/3
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD(:,2,:,:,:),KGRID,KLENCH,HCOMMENT,KRESP)
ELSE
  ILENG=SIZE(PFIELD)
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD,KGRID,KLENCH,HCOMMENT,KRESP)
END IF
!-------------------------------------------------------------------------------
END SUBROUTINE FMWRITX5
!
!     #############################################################
      SUBROUTINE FMWRITX6(HFILEM,HRECFM,HFIPRI,KLENG,PFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMWRITX6* - routine to write a real 6D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITX0 is to convert the real into integer(kind=8)
!     by calling FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
REAL, DIMENSION(:,:,:,:,:,:), &
                           INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER :: ILENG
REAL(KIND=8),    &
DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2),SIZE(PFIELD,3),SIZE(PFIELD,4),SIZE(PFIELD,5),SIZE(PFIELD,6)) :: ZFIELD
!-------------------------------------------------------------------------------
!
ZFIELD=PFIELD
ILENG=SIZE(PFIELD)
CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,ZFIELD,KGRID,KLENCH,HCOMMENT,KRESP)
!-------------------------------------------------------------------------------
END SUBROUTINE FMWRITX6
!
!     #############################################################
      SUBROUTINE FMWRITN0(HFILEM,HRECFM,HFIPRI,KLENG,KFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMWRITN0* - routine to write a integer scalar into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITN0 is to convert the integer into integer(kind=8)
!     by calling FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
INTEGER, &
                           INTENT(IN) ::KFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER(KIND=8) :: IFIELD
!-------------------------------------------------------------------------------
IFIELD=KFIELD
CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,1,IFIELD,KGRID,KLENCH,HCOMMENT,KRESP)
!-------------------------------------------------------------------------------
END SUBROUTINE FMWRITN0
!
!     #############################################################
      SUBROUTINE FMWRITN1(HFILEM,HRECFM,HFIPRI,KLENG,KFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMWRITN1* - routine to write a integer 1D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITN1 is to convert the integer into integer(kind=8)
!     by calling FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
INTEGER, DIMENSION(:), &
                           INTENT(IN) ::KFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER(KIND=8), DIMENSION(SIZE(KFIELD)) :: IFIELD
INTEGER                                  :: ILENG
!-------------------------------------------------------------------------------
!
ILENG=SIZE(KFIELD)
IFIELD(:)=KFIELD(:)
CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,IFIELD,KGRID,KLENCH,HCOMMENT,KRESP)
!-------------------------------------------------------------------------------
END SUBROUTINE FMWRITN1
!
!     #############################################################
      SUBROUTINE FMWRITN2(HFILEM,HRECFM,HFIPRI,KLENG,KFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMWRITN2* - routine to write a integer 2D array into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITN1 is to convert the integer into integer(kind=8)
!     by calling FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
                           INTENT(IN) ::KFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER(KIND=8), DIMENSION(SIZE(KFIELD,1),SIZE(KFIELD,2)) :: IFIELD
INTEGER                                                   :: ILENG
!-------------------------------------------------------------------------------
!
IFIELD(:,:)=KFIELD(:,:)
!
IF (LPACK .AND. L1D .AND. SIZE(KFIELD,1)==3 .AND. SIZE(KFIELD,2)==3) THEN
  ILENG=SIZE(KFIELD)/9
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,IFIELD(2,2),KGRID,KLENCH,HCOMMENT,KRESP)
ELSE IF (LPACK .AND. L2D .AND. SIZE(KFIELD,2)==3) THEN
  ILENG=SIZE(KFIELD)/3
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,IFIELD(:,2),KGRID,KLENCH,HCOMMENT,KRESP)
ELSE
  ILENG=SIZE(KFIELD)
  CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,IFIELD,KGRID,KLENCH,HCOMMENT,KRESP)
END IF
!-------------------------------------------------------------------------------
END SUBROUTINE FMWRITN2
!
!     #############################################################
      SUBROUTINE FMWRITL0(HFILEM,HRECFM,HFIPRI,KLENG,OFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMWRITL0* - routine to write a logical scalar into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITN0 is to convert the integer into integer(kind=8)
!     by calling FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
LOGICAL, &
                           INTENT(IN) ::OFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER(KIND=8) :: IFIELD
!-------------------------------------------------------------------------------
!
IF (OFIELD) THEN
  IFIELD=1
ELSE
  IFIELD=0
END IF
!
CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,1,IFIELD,KGRID,KLENCH,HCOMMENT,KRESP)
!-------------------------------------------------------------------------------
END SUBROUTINE FMWRITL0
!     #############################################################
      SUBROUTINE FMWRITL1(HFILEM,HRECFM,HFIPRI,KLENG,OFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMWRITL0* - routine to write a logical scalar into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITN0 is to convert the integer into integer(kind=8)
!     by calling FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
LOGICAL, DIMENSION(:), &
                           INTENT(IN) ::OFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER(KIND=8), DIMENSION(SIZE(OFIELD)) :: IFIELD
!-------------------------------------------------------------------------------
!
WHERE (OFIELD)
  IFIELD=1
    ELSEWHERE
  IFIELD=0
END WHERE
!
CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,SIZE(IFIELD),IFIELD,KGRID,KLENCH,HCOMMENT,KRESP)
!-------------------------------------------------------------------------------
END SUBROUTINE FMWRITL1
!     #############################################################
      SUBROUTINE FMWRITC0(HFILEM,HRECFM,HFIPRI,KLENG,HFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMWRITC0* - routine to write a string scalar into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITL0 is to convert the string into arrayr of
!      integer(kind=8) and to call FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
CHARACTER(LEN=*), &
                           INTENT(IN) ::HFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)     ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER                                               :: JLOOP
INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE            :: IFIELD
INTEGER                                               :: ILENG
!-------------------------------------------------------------------------------
!
ILENG=LEN(HFIELD)
ALLOCATE(IFIELD(ILENG))
DO JLOOP=1,ILENG
 IFIELD(JLOOP)=IACHAR(HFIELD(JLOOP:JLOOP))
END DO
!
CALL FM_WRIT(HFILEM,HRECFM,HFIPRI,ILENG,IFIELD,KGRID,KLENCH,HCOMMENT,KRESP)
!
DEALLOCATE(IFIELD)
!-------------------------------------------------------------------------------
END SUBROUTINE FMWRITC0
!     #############################################################
      SUBROUTINE FMWRITT0(HFILEM,HRECFM,HFIPRI,KLENG,TFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!     #############################################################
!
!!****  *FMWRITT0* - routine to write a date scalar into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITT0 is to split a date_time scalar
!      and to call FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
                           INTENT(IN) ::TFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)          ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER(KIND=8), DIMENSION(3)  :: ITDATE    ! date array
CHARACTER(LEN=16)              :: YRECFM    ! Name of the article to be written
CHARACTER(LEN=JPXKRK)          :: YCOMMENT  ! Comment string
!
!-------------------------------------------------------------------------------
!
YRECFM=TRIM(HRECFM)//'%TDATE'   ! array of rank 3 for date is written in file
YCOMMENT='YYYYMMDD'
ITDATE(1)=TFIELD%TDATE%YEAR
ITDATE(2)=TFIELD%TDATE%MONTH
ITDATE(3)=TFIELD%TDATE%DAY
CALL FM_WRIT(HFILEM,YRECFM,HFIPRI,3,ITDATE,0,8,YCOMMENT,KRESP)
!
YRECFM=TRIM(HRECFM)//'%TIME'
YCOMMENT='SECONDS'
CALL FM_WRIT(HFILEM,YRECFM,HFIPRI,1,TFIELD%TIME,0,7,YCOMMENT,KRESP)
!
!
!-------------------------------------------------------------------------------
END SUBROUTINE FMWRITT0


