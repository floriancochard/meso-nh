!MNH_LIC Copyright 2003-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-----------------------------------------------------------------
MODULE MODE_READ_SURF_MNH_TOOLS

IMPLICIT NONE

CONTAINS

SUBROUTINE PREPARE_METADATA_READ_SURF(HREC,HDIR,KGRID,KTYPE,KDIMS,HSUBR,TPFIELD)
!
USE MODE_FIELD, ONLY: FIND_FIELD_ID_FROM_MNHNAME, TFIELDDATA, TFIELDLIST, TYPECHAR, TYPEDATE, TYPELOG
USE MODE_MSG
!
CHARACTER(LEN=LEN_HREC),INTENT(IN)  :: HREC     ! name of the article to write
CHARACTER(LEN=2),       INTENT(IN)  :: HDIR     ! Expected type of the data field (XX,XY,--...)
INTEGER,                INTENT(IN)  :: KGRID    ! Localization on the model grid
INTEGER,                INTENT(IN)  :: KTYPE    ! Datatype
INTEGER,                INTENT(IN)  :: KDIMS    ! Number of dimensions
CHARACTER(LEN=*),       INTENT(IN)  :: HSUBR    ! name of the subroutine calling
TYPE(TFIELDDATA),       INTENT(OUT) :: TPFIELD  ! metadata of field
!
CHARACTER(LEN=32) :: YTXT
INTEGER           :: IID, IRESP
!
CALL FIND_FIELD_ID_FROM_MNHNAME(TRIM(HREC),IID,IRESP,ONOWARNING=.TRUE.)
IF (IRESP==0) THEN
  TPFIELD = TFIELDLIST(IID)
  !Modify and check CLONGNAME
  IF (TRIM(TPFIELD%CLONGNAME)/=TRIM(HREC) &
      .AND. TRIM(HREC)/='VERSION' .AND. TRIM(HREC)/='BUG') THEN
    CALL PRINT_MSG(NVERB_WARNING,'IO',TRIM(HSUBR),'CLONGNAME different ('//TRIM(TPFIELD%CLONGNAME) &
                   //'/'//TRIM(HREC)//') than expected for article '//TRIM(HREC))
    TPFIELD%CLONGNAME = TRIM(HREC)
  END IF
  !Modify and check CDIR
  IF (TPFIELD%CDIR/=HDIR) THEN
    CALL PRINT_MSG(NVERB_WARNING,'IO',TRIM(HSUBR),'CDIR different ('//TRIM(TPFIELD%CDIR) &
                   //'/'//TRIM(HDIR)//') than expected for article '//TRIM(HREC))
    TPFIELD%CDIR = HDIR
  END IF
  !Modify and check NGRID
  IF (TPFIELD%NGRID/=KGRID) THEN
    WRITE(YTXT,'( I0,"/",I0 )') TPFIELD%NGRID,KGRID
    CALL PRINT_MSG(NVERB_WARNING,'IO',TRIM(HSUBR),'NGRID different ('//TRIM(YTXT) &
                    //') than expected for article '//TRIM(HREC))
    TPFIELD%NGRID = KGRID
  END IF
  !Modify and check NTYPE
  IF (TPFIELD%NTYPE/=KTYPE) THEN
    WRITE(YTXT,'( I0,"/",I0 )') TPFIELD%NTYPE,KTYPE
    CALL PRINT_MSG(NVERB_WARNING,'IO',TRIM(HSUBR),'NTYPE different ('//TRIM(YTXT) &
                    //') than expected for article '//TRIM(HREC))
    TPFIELD%NTYPE = KTYPE
  END IF
  !Modify and check NDIMS
  IF (TPFIELD%NDIMS/=KDIMS) THEN
    WRITE(YTXT,'( I0,"/",I0 )') TPFIELD%NDIMS,KDIMS
    CALL PRINT_MSG(NVERB_WARNING,'IO',TRIM(HSUBR),'NDIMS different ('//TRIM(YTXT) &
                    //') than expected for article '//TRIM(HREC))
    TPFIELD%NDIMS = KDIMS
  END IF
ELSE
  CALL PRINT_MSG(NVERB_DEBUG,'IO',TRIM(HSUBR),TRIM(HREC)//' not found in FIELDLIST. Generating default metadata')
  TPFIELD%CMNHNAME   = TRIM(HREC)
  TPFIELD%CSTDNAME   = ''
  TPFIELD%CLONGNAME  = TRIM(HREC)
  TPFIELD%CUNITS     = ''
  TPFIELD%CDIR       = HDIR
  TPFIELD%CCOMMENT   = '' !Expected comment is not known
  TPFIELD%NGRID      = KGRID
  TPFIELD%NTYPE      = KTYPE
  TPFIELD%NDIMS      = KDIMS
#if 0
  IF (TPFIELD%NDIMS==0 .OR. TPFIELD%NTYPE==TYPECHAR .OR. TPFIELD%NTYPE==TYPEDATE .OR. TPFIELD%NTYPE==TYPELOG) THEN
    TPFIELD%LTIMEDEP   = .FALSE.
  ELSE
    TPFIELD%LTIMEDEP   = .TRUE.
  END IF
#else
  TPFIELD%LTIMEDEP   = .FALSE.
#endif

END IF
!
END SUBROUTINE PREPARE_METADATA_READ_SURF

END MODULE MODE_READ_SURF_MNH_TOOLS


!     #############################################################
      SUBROUTINE READ_SURFX0_MNH(HREC,PFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READX0* - routine to read a real scalar
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READX0 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
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
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!!      10/10/2011 J.Escobar & G.Tanguy change BUGFIX/MNH to BUG/SURFEX version control
!!  01/2018      (G.Delautier) SURFEX 8.1
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_CONF,        ONLY: CPROGRAM
USE MODD_GRID,        ONLY: XRPK,XBETA,XLAT0,XLON0
USE MODD_IO_SURF_MNH, ONLY: TOUT, TPINFILE
USE MODD_PARAMETERS,  ONLY: JPHEXT, XUNDEF
!
USE MODE_FIELD,       ONLY: TFIELDDATA,TFIELDLIST,FIND_FIELD_ID_FROM_MNHNAME,TYPEREAL
USE MODE_FM
USE MODE_FMREAD
USE MODE_GRIDPROJ
USE MODE_MSG
USE MODE_READ_SURF_MNH_TOOLS
!
USE MODI_GET_SURF_UNDEF
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=LEN_HREC),INTENT(IN)  :: HREC     ! name of the article to be read
REAL,                   INTENT(OUT) :: PFIELD   ! the real scalar to be read
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
INTEGER           :: IGRID          ! IGRID : grid indicator
INTEGER           :: ILENCH         ! ILENCH : length of comment string
INTEGER           :: ILUOUT
INTEGER           :: IID,IRESP
INTEGER           :: IIMAX,IJMAX
REAL,DIMENSION(:), ALLOCATABLE :: ZXHAT,ZYHAT
REAL              :: ZLATOR,ZLONOR,ZXHATM,ZYHATM,ZLATORI,ZLONORI
REAL              :: ZRPK, ZBETA, ZLAT0, ZLON0
TYPE(TFIELDDATA)  :: TZFIELD
!-------------------------------------------------------------------------------
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','READ_SURFX0_MNH',TRIM(TPINFILE%CNAME)//': reading '//TRIM(HREC))
!
ILUOUT = TOUT%NLU
!
IF (HREC=='LONORI' .OR. HREC=='LATORI') THEN
  IF (TPINFILE%NMNHVERSION(1)<4 .OR. (TPINFILE%NMNHVERSION(1)==4 .AND. TPINFILE%NMNHVERSION(2)<=5)) THEN
      ZLATORI = XUNDEF
      ZLONORI = XUNDEF
      !* saves projection parameters of MODD_GRID
      ZLAT0 = XLAT0
      ZLON0 = XLON0
      ZRPK  = XRPK
      ZBETA = XBETA
      !* reads projection and grid data in the file
      CALL IO_READ_FIELD(TPINFILE,'LAT0',XLAT0)
      CALL IO_READ_FIELD(TPINFILE,'LON0',XLON0)
      CALL IO_READ_FIELD(TPINFILE,'RPK', XRPK)
      CALL IO_READ_FIELD(TPINFILE,'BETA',XBETA)
      !
      CALL IO_READ_FIELD(TPINFILE,'IMAX',IIMAX)
      CALL IO_READ_FIELD(TPINFILE,'JMAX',IJMAX)
      ALLOCATE(ZXHAT(IIMAX+2*JPHEXT),ZYHAT(IJMAX+2*JPHEXT))
      CALL IO_READ_FIELD(TPINFILE,'XHAT',ZXHAT)
      CALL IO_READ_FIELD(TPINFILE,'YHAT',ZYHAT)
      !
      CALL FIND_FIELD_ID_FROM_MNHNAME('LONORI',IID,IRESP)
      TZFIELD = TFIELDLIST(IID)
      TZFIELD%CMNHNAME = 'LONOR'
      CALL IO_READ_FIELD(TPINFILE,TZFIELD,ZLONOR)
      !
      CALL FIND_FIELD_ID_FROM_MNHNAME('LATORI',IID,IRESP)
      TZFIELD = TFIELDLIST(IID)
      TZFIELD%CMNHNAME = 'LATOR'
      CALL IO_READ_FIELD(TPINFILE,TZFIELD,ZLATOR)
      !
      ZXHATM = - 0.5 * (ZXHAT(1)+ZXHAT(2))
      ZYHATM = - 0.5 * (ZYHAT(1)+ZYHAT(2))
      DEALLOCATE(ZXHAT,ZYHAT)
      !* computes origin
      CALL SM_LATLON(ZLATOR,ZLONOR,ZXHATM,ZYHATM,ZLATORI,ZLONORI)
      IF (HREC=='LONORI') PFIELD = ZLONORI
      IF (HREC=='LATORI') PFIELD = ZLATORI
      !* restores projection parameters in module MODD_GRID
      XLAT0 = ZLAT0
      XLON0 = ZLON0
      XRPK  = ZRPK
      XBETA = ZBETA
      RETURN
  END IF
END IF

!-------------------------------------------------------------------------------

IF ( HREC=='LAT0' .OR. HREC=='LON0' .OR. HREC=='RPK' .OR. HREC=='BETA'  &
                 .OR. HREC=='LATORI'.OR. HREC=='LONORI'                  ) THEN
  CALL IO_READ_FIELD(TPINFILE,HREC,PFIELD,KRESP)
ELSE
  CALL PREPARE_METADATA_READ_SURF(HREC,'--',0,TYPEREAL,0,'READ_SURFX0_MNH',TZFIELD)
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,PFIELD,KRESP)
END IF

IF (KRESP /=0) THEN
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'error when reading article ', HREC,'KRESP=',KRESP
  WRITE(ILUOUT,*) 'default value may be used, who knows???'
  WRITE(ILUOUT,*) ' '
ENDIF
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX0_MNH
!
!     #############################################################
      SUBROUTINE READ_SURFX1_MNH(HREC,KL,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX1* - routine to fill a real 1D array for the externalised surface
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READ_SURFX1 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
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
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_FIELD,       ONLY: FIND_FIELD_ID_FROM_MNHNAME,TFIELDDATA,TFIELDLIST,TYPEREAL
USE MODE_FM
USE MODE_FMREAD
USE MODE_ll
USE MODE_IO_ll
USE MODE_MSG
USE MODE_READ_SURF_MNH_TOOLS
!
USE MODD_CST,         ONLY : XPI
!
USE MODD_IO_SURF_MNH, ONLY : TOUT, TPINFILE, NMASK, &
                             NIU, NJU, NIB, NJB, NIE, NJE, &
                             NIU_ALL, NJU_ALL, NIB_ALL,    &
                             NJB_ALL, NIE_ALL, NJE_ALL,    &
                             NMASK_ALL
USE MODD_PARAMETERS, ONLY: XUNDEF
!
USE MODI_PACK_2D_1D
!
USE MODI_GET_SURF_UNDEF
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=LEN_HREC),INTENT(IN) :: HREC     ! name of the article to be read
INTEGER,                INTENT(IN) :: KL       !  number of points
REAL, DIMENSION(KL),    INTENT(OUT):: PFIELD   ! array containing the data field
INTEGER,                INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT):: HCOMMENT ! comment
CHARACTER(LEN=1),       INTENT(IN) :: HDIR     ! type of field :
!                                           ! 'H' for HOR : with hor. dim.; and  distributed.
!                                           ! 'A' for ALL : with hor. dim.; and not distributed.
!                                           ! '-' : no horizontal dim.

!
!*      0.2   Declarations of local variables
!
INTEGER           :: IGRID          ! IGRID : grid indicator
INTEGER           :: ILENCH         ! ILENCH : length of comment string
INTEGER           :: ILUOUT
INTEGER           :: JI, JJ         ! loop counters

REAL, DIMENSION(:,:), ALLOCATABLE :: ZWORK  ! work array read in the file
REAL, DIMENSION(:),   ALLOCATABLE :: ZWORK1D! work array read in the file
REAL                              :: ZW     ! work value

CHARACTER(LEN=LEN_HREC) :: YREC
CHARACTER(LEN=2)  :: YSTORAGE_TYPE
!
INTEGER           :: IID, IRESP
INTEGER           :: IIU, IJU, IIB, IJB, IIE, IJE ! dimensions of horizontal fields
INTEGER, DIMENSION(:), ALLOCATABLE :: IMASK       ! mask for packing
REAL              :: ZUNDEF         ! undefined value in SURFEX
TYPE(TFIELDDATA)  :: TZFIELD
!-------------------------------------------------------------------------------
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','READ_SURFX1_MNH',TRIM(TPINFILE%CNAME)//': reading '//TRIM(HREC))
!
KRESP = 0
ILUOUT = TOUT%NLU
!
IF (HDIR=='A'.OR.HDIR=='E') THEN
  IIU = NIU_ALL
  IJU = NJU_ALL
  IIB = NIB_ALL
  IJB = NJB_ALL
  IIE = NIE_ALL
  IJE = NJE_ALL
  ALLOCATE(IMASK(SIZE(NMASK_ALL)))
  IMASK = NMASK_ALL
ELSE
  IIU = NIU
  IJU = NJU
  IIB = NIB
  IJB = NJB
  IIE = NIE
  IJE = NJE
  ALLOCATE(IMASK(SIZE(NMASK)))
  IMASK = NMASK
END IF
!
!*       2.    On traite d'abord des cas particuliers
!
IF (HREC=='LAT') THEN

  CALL IO_READ_FIELD(TPINFILE,'LAT0',ZW,KRESP)
  PFIELD(:) = ZW

ELSE IF (HREC=='LON') THEN

  CALL IO_READ_FIELD(TPINFILE,'LON0',ZW,KRESP)
  PFIELD(:) = ZW

ELSE IF (HREC=='MESH_SIZE') THEN

  PFIELD(:) = 0.
  HCOMMENT = ' '

ELSE IF (HREC=='XX') THEN
!! reading of a 1D field along X in the file
  ALLOCATE(ZWORK1D(IIU))
  ALLOCATE(ZWORK  (IIU,IJU))
  ZWORK(:,:) = 0.
  CALL FIND_FIELD_ID_FROM_MNHNAME('XHAT',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  IF (HDIR/='A'.AND.HDIR/='E') THEN
    TZFIELD%CDIR       = 'XX'
  ELSE
    TZFIELD%CDIR       = '--'
  END IF
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,ZWORK1D,KRESP)
  DO JJ = 1,IJU
    ZWORK(IIB:IIE,JJ) = 0.5 * ZWORK1D(IIB:IIE) + 0.5 * ZWORK1D(IIB+1:IIE+1)
  END DO
  CALL PACK_2D_1D(IMASK,ZWORK(IIB:IIE,IJB:IJE),PFIELD)
  DEALLOCATE(ZWORK1D)
  DEALLOCATE(ZWORK  )
ELSE IF (HREC=='DX') THEN
!! reading of a 1D field along X in the file
  ALLOCATE(ZWORK1D(IIU))
  ALLOCATE(ZWORK  (IIU,IJU))
  ZWORK(:,:) = 0.
  CALL FIND_FIELD_ID_FROM_MNHNAME('XHAT',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  IF (HDIR/='A'.AND.HDIR/='E') THEN
    TZFIELD%CDIR       = 'XX'
  ELSE
    TZFIELD%CDIR       = '--'
  END IF
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,ZWORK1D,KRESP)
  DO JJ = 1,IJU
    ZWORK(IIB:IIE,JJ) = - ZWORK1D(IIB:IIE) + ZWORK1D(IIB+1:IIE+1)
  END DO
  CALL PACK_2D_1D(IMASK,ZWORK(IIB:IIE,IJB:IJE),PFIELD)
  DEALLOCATE(ZWORK1D)
  DEALLOCATE(ZWORK  )
ELSE IF (HREC=='YY') THEN
!! reading of a 1D field along Y in the file
  ALLOCATE(ZWORK1D(IJU))
  ALLOCATE(ZWORK  (IIU,IJU))
  ZWORK(:,:) = 0.
  CALL FIND_FIELD_ID_FROM_MNHNAME('YHAT',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  IF (HDIR/='A'.AND.HDIR/='E') THEN
    TZFIELD%CDIR       = 'YY'
  ELSE
    TZFIELD%CDIR       = '--'
  END IF
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,ZWORK1D,KRESP)
  DO JI = 1,IIU
    ZWORK(JI,IJB:IJE) = 0.5 * ZWORK1D(IJB:IJE) + 0.5 * ZWORK1D(IJB+1:IJE+1)
  END DO
  CALL PACK_2D_1D(IMASK,ZWORK(IIB:IIE,IJB:IJE),PFIELD)
  DEALLOCATE(ZWORK1D)
  DEALLOCATE(ZWORK  )
ELSE IF (HREC=='DY') THEN
!! reading of a 1D field along Y in the file
  ALLOCATE(ZWORK1D(IJU))
  ALLOCATE(ZWORK  (IIU,IJU))
  ZWORK(:,:) = 0.
  CALL FIND_FIELD_ID_FROM_MNHNAME('YHAT',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  IF (HDIR/='A'.AND.HDIR/='E') THEN
    TZFIELD%CDIR       = 'YY'
  ELSE
    TZFIELD%CDIR       = '--'
  END IF
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,ZWORK1D,KRESP)
  DO JI = 1,IIU
    ZWORK(JI,IJB:IJE) = - ZWORK1D(IJB:IJE) + ZWORK1D(IJB+1:IJE+1)
  END DO
  CALL PACK_2D_1D(IMASK,ZWORK(IIB:IIE,IJB:IJE),PFIELD)
  DEALLOCATE(ZWORK1D)
  DEALLOCATE(ZWORK  )
!
ELSE
!
!! Reading of a 2D fields, masked and packed into 1D vector
!
  YREC = ' '
  YREC(1:LEN(HREC)) = HREC
  IF (HREC(1:8)=='Q_CANYON') THEN
    IF (TPINFILE%NMNHVERSION(1)<4 .OR. (TPINFILE%NMNHVERSION(1)==4 .AND. TPINFILE%NMNHVERSION(2)<=5)) THEN
      CALL IO_READ_FIELD(TPINFILE,'STORAGE_TYPE',YSTORAGE_TYPE)
      IF (YSTORAGE_TYPE=='TT') THEN
        PFIELD = 0.
        DEALLOCATE(IMASK)
        RETURN
      ELSE
        YREC = 'R_CANYON            '
      END IF
    END IF
  END IF
  IF (HREC(1:8)=='T_CANYON') THEN
    IF (TPINFILE%NMNHVERSION(1)<4 .OR. (TPINFILE%NMNHVERSION(1)==4 .AND. TPINFILE%NMNHVERSION(2)<=5)) THEN
      CALL IO_READ_FIELD(TPINFILE,'STORAGE_TYPE',YSTORAGE_TYPE)
      IF (YSTORAGE_TYPE=='TT') YREC = 'T_ROAD1             '
    END IF
  END IF
  IF (HREC(1:7)=='SSO_DIR') THEN
    IF (TPINFILE%NMNHVERSION(1)<4 .OR. (TPINFILE%NMNHVERSION(1)==4 .AND. TPINFILE%NMNHVERSION(2)<=5)) THEN
      YREC = 'SSO_DIRECTION       '
    END IF
  END IF
!
  ALLOCATE(ZWORK(IIU,IJU))
!
  IF (HDIR=='H') THEN
    CALL PREPARE_METADATA_READ_SURF(YREC,'XY',4,TYPEREAL,2,'READ_SURFX1_MNH',TZFIELD)
    CALL IO_READ_FIELD(TPINFILE,TZFIELD,ZWORK,KRESP)
  ELSEIF (HDIR=='A'.OR.HDIR=='E') THEN
    CALL PREPARE_METADATA_READ_SURF(YREC,'--',4,TYPEREAL,2,'READ_SURFX1_MNH',TZFIELD)
    CALL IO_READ_FIELD(TPINFILE,TZFIELD,ZWORK,KRESP)
  ELSE
    CALL PREPARE_METADATA_READ_SURF(YREC,'--',4,TYPEREAL,1,'READ_SURFX1_MNH',TZFIELD)
    CALL IO_READ_FIELD(TPINFILE,TZFIELD,PFIELD,KRESP)
  END IF
!
  IF (KRESP /=0) THEN
    WRITE(ILUOUT,*) 'WARNING'
    WRITE(ILUOUT,*) '-------'
    WRITE(ILUOUT,*) 'error when reading article ', HREC,'KRESP=',KRESP
    WRITE(ILUOUT,*) 'default value may be used, who knows???'
    WRITE(ILUOUT,*) ' '
  ELSE IF (HDIR=='H' .OR. HDIR=='A' .OR. HDIR=='E') THEN
    CALL PACK_2D_1D(IMASK,ZWORK(IIB:IIE,IJB:IJE),PFIELD)
    CALL GET_SURF_UNDEF(ZUNDEF)
!================================================
! 13/03/2009 : G. TANGUY
! on supprime le test sur lesvaleurs indéfinies 
! pour l'orographie pour que l'altitude 999 m 
! soit autorisée
    IF (HREC(1:2)/='ZS') THEN
      WHERE (PFIELD==XUNDEF) PFIELD=ZUNDEF
    ENDIF
!================================================

  END IF
!
  DEALLOCATE(ZWORK)

ENDIF

DEALLOCATE(IMASK)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX1_MNH
!
!     #############################################################
      SUBROUTINE READ_SURFX2_MNH(HREC,KL1,KL2,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX2* - routine to fill a real 2D array for the externalised surface
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READ_SURFX2 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
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
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_ll
USE MODE_FIELD,       ONLY : TFIELDDATA,TYPEREAL
USE MODE_FM
USE MODE_FMREAD
USE MODE_IO_ll
USE MODE_MSG
USE MODE_READ_SURF_MNH_TOOLS
!
USE MODD_IO_SURF_MNH, ONLY : TOUT, TPINFILE, NMASK, NIU, NJU, NIB, NJB, NIE, NJE, &
                             NIU_ALL, NJU_ALL, NIB_ALL, NJB_ALL, NIE_ALL, NJE_ALL, NMASK_ALL
USE MODD_PARAMETERS,  ONLY : XUNDEF
!
USE MODI_PACK_2D_1D
!
USE MODI_GET_SURF_UNDEF
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=LEN_HREC), INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                 INTENT(IN)  :: KL1      ! number of points
INTEGER,                 INTENT(IN)  :: KL2      ! second dimension
REAL, DIMENSION(KL1,KL2),INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,                 INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),      INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),        INTENT(IN)  :: HDIR     ! type of field :
!                                                ! 'H' for HOR : with hor. dim.; and  distributed.
!                                                ! 'A' for ALL : with hor. dim.; and not distributed.
!                                                ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
INTEGER           :: IGRID          ! IGRID : grid indicator
INTEGER           :: ILENCH         ! ILENCH : length of comment string
INTEGER           :: ILUOUT
INTEGER           :: JP             ! loop index

REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWORK  ! work array read in the file
REAL              :: ZUNDEF         ! undefined value in SURFEX
TYPE(TFIELDDATA)  :: TZFIELD
!-------------------------------------------------------------------------------
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','READ_SURFX2_MNH',TRIM(TPINFILE%CNAME)//': reading '//TRIM(HREC))
!
ILUOUT = TOUT%NLU
!
!! Reading of a 3D field, masked (2 first dimensions) and with
!! 2 first dimensions packed into only 1 (results in a 2D array instead of 3D)
!
!*       1.     Dimension initializations:
!               -------------------------
!
!
!
IF (HDIR=='H') THEN
  ALLOCATE(ZWORK(NIU,NJU,SIZE(PFIELD,2)))
  CALL PREPARE_METADATA_READ_SURF(HREC,'XY',4,TYPEREAL,3,'READ_SURFX2_MNH',TZFIELD)
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,ZWORK,KRESP)
ELSEIF (HDIR=='A'.OR.HDIR=='E') THEN
  ALLOCATE(ZWORK(NIU_ALL,NJU_ALL,SIZE(PFIELD,2)))
  CALL PREPARE_METADATA_READ_SURF(HREC,'--',4,TYPEREAL,3,'READ_SURFX2_MNH',TZFIELD)
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,ZWORK,KRESP)
ELSE
  CALL PREPARE_METADATA_READ_SURF(HREC,'--',4,TYPEREAL,2,'READ_SURFX2_MNH',TZFIELD)
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,PFIELD,KRESP)
END IF
!
IF (KRESP /=0) THEN
    WRITE(ILUOUT,*) 'WARNING'
    WRITE(ILUOUT,*) '-------'
    WRITE(ILUOUT,*) 'error when reading article ', HREC,'KRESP=',KRESP
    WRITE(ILUOUT,*) 'default value may be used, who knows???'
    WRITE(ILUOUT,*) ' '
    DEALLOCATE(ZWORK)
 ELSE IF (HDIR=='H') THEN
    DO JP=1,SIZE(PFIELD,2)
       CALL PACK_2D_1D(NMASK,ZWORK(NIB:NIE,NJB:NJE,JP),PFIELD(:,JP))
    END DO
    DEALLOCATE(ZWORK)
 ELSE IF (HDIR=='A'.OR.HDIR=='E') THEN
    DO JP=1,SIZE(PFIELD,2)
       CALL PACK_2D_1D(NMASK_ALL,ZWORK(NIB_ALL:NIE_ALL,NJB_ALL:NJE_ALL,JP),PFIELD(:,JP))
    END DO
    DEALLOCATE(ZWORK)
 END IF
 CALL GET_SURF_UNDEF(ZUNDEF)
 WHERE (PFIELD==XUNDEF) PFIELD=ZUNDEF
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX2_MNH
!
!     #############################################################
      SUBROUTINE READ_SURFX2COV_MNH(HREC,KL1,KL2,PFIELD,OFLAG,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX1* - routine to fill a real 2D array for the externalised surface
!!                 with Logical mask by level
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READ_SURFX1 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
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
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!!  06/2016     (G.Delautier) phasage surfex 8
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_ll
USE MODE_FIELD,          ONLY: TFIELDDATA,TYPELOG,TYPEREAL
USE MODE_FM
USE MODE_FMREAD
USE MODE_IO_ll
USE MODE_MSG
USE MODE_READ_SURF_MNH_TOOLS
!
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
USE MODD_CST,            ONLY : XPI
!
USE MODD_IO_SURF_MNH,    ONLY : TOUT, TPINFILE, NMASK, &
                                NIU, NJU, NIB, NJB, NIE, NJE, &
                                NIU_ALL, NJU_ALL, NIB_ALL,    &
                                NJB_ALL, NIE_ALL, NJE_ALL,    &
                                NMASK_ALL
!
USE MODI_PACK_2D_1D
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=LEN_HREC),   INTENT(IN) :: HREC     ! name of the article to be read
INTEGER,                   INTENT(IN) :: KL1,KL2  !  number of points
REAL, DIMENSION(KL1,KL2),  INTENT(OUT):: PFIELD   ! array containing the data field
LOGICAL,DIMENSION(JPCOVER),INTENT(IN) :: OFLAG    ! mask for array filling
INTEGER,                   INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),        INTENT(OUT):: HCOMMENT ! comment
CHARACTER(LEN=1),          INTENT(IN) :: HDIR     ! type of field :
!                                           ! 'H' for HOR : with hor. dim.; and  distributed.
!                                           ! 'A' for ALL : with hor. dim.; and not distributed.
!                                           ! '-' : no horizontal dim.

!
!*      0.2   Declarations of local variables
!
INTEGER           :: IGRID          ! IGRID : grid indicator
INTEGER           :: ILENCH         ! ILENCH : length of comment string
INTEGER           :: ILUOUT
!
CHARACTER(LEN=LEN_HREC) :: YREC
CHARACTER(LEN=2)  :: YDIR
CHARACTER(LEN=2)  :: YSTORAGE_TYPE
!
INTEGER           :: IIU, IJU, IIB, IJB, IIE, IJE ! dimensions of horizontal fields
INTEGER, DIMENSION(:), ALLOCATABLE :: IMASK       ! mask for packing
!JUANZ
INTEGER           :: NCOVER,ICOVER,JL2
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZWORK3D
!JUANZ
INTEGER  :: IRESP
INTEGER  :: IVERSION, IBUGFIX
LOGICAL  :: GCOVER_PACKED ! .T. if COVER are all packed into one field
TYPE(TFIELDDATA) :: TZFIELD
!-------------------------------------------------------------------------------
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','READ_SURFX2COV_MNH',TRIM(TPINFILE%CNAME)//': reading '//TRIM(HREC))
!
KRESP = 0
IRESP = 0
ILUOUT = TOUT%NLU
!
IF (HDIR=='A'.OR.HDIR=='E') THEN
  YDIR="--"
  IIU = NIU_ALL
  IJU = NJU_ALL
  IIB = NIB_ALL
  IJB = NJB_ALL
  IIE = NIE_ALL
  IJE = NJE_ALL
  ALLOCATE(IMASK(SIZE(NMASK_ALL)))
  IMASK = NMASK_ALL
ELSE
  YDIR="XY"
  IIU = NIU
  IJU = NJU
  IIB = NIB
  IJB = NJB
  IIE = NIE
  IJE = NJE
  ALLOCATE(IMASK(SIZE(NMASK)))
  IMASK = NMASK
END IF
!
!! Reading of a 2D fields, masked and packed into 1D vector
!
!
NCOVER=COUNT(OFLAG)
ALLOCATE (ZWORK3D(IIU,IJU,NCOVER))
ZWORK3D(:,:,:) =  0.0
!
CALL IO_READ_FIELD(TPINFILE,'VERSION',IVERSION)
CALL IO_READ_FIELD(TPINFILE,'BUG',    IBUGFIX)

IF (IVERSION<7 .OR. (IVERSION==7 .AND. IBUGFIX==0)) THEN
  GCOVER_PACKED = .FALSE.
ELSE
  TZFIELD%CMNHNAME   = 'COVER_PACKED'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'COVER_PACKED'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = '--'
  TZFIELD%CCOMMENT   = ''
  TZFIELD%NGRID      = 0
  TZFIELD%NTYPE      = TYPELOG
  TZFIELD%NDIMS      = 0
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,GCOVER_PACKED)
END IF
!
IF (.NOT. GCOVER_PACKED) THEN
  ICOVER=0
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = ''
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .FALSE.
  DO JL2=1,SIZE(OFLAG)
    WRITE(YREC,'(A5,I3.3)') 'COVER',JL2
    TZFIELD%CMNHNAME   = TRIM(YREC)
    TZFIELD%CLONGNAME  = TRIM(YREC)
    TZFIELD%CCOMMENT   = 'X_Y_'//TRIM(YREC)
    TZFIELD%CDIR       = YDIR
    IF (OFLAG(JL2)) THEN
      ICOVER=ICOVER+1
      CALL IO_READ_FIELD(TPINFILE,TZFIELD,ZWORK3D(:,:,ICOVER),IRESP)
    END IF
    IF (IRESP/=0) KRESP=IRESP
  END DO
ELSE
  CALL PREPARE_METADATA_READ_SURF(HREC,YDIR,4,TYPEREAL,3,'READ_SURFX2COV_MNH',TZFIELD)
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,ZWORK3D(:,:,:),KRESP)
END IF
!
IF (KRESP /=0) THEN
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'error when reading article ', HREC,'KRESP=',KRESP
  WRITE(ILUOUT,*) ' '
ELSE IF (HDIR=='H' .OR. HDIR=='A' .OR. HDIR=='E') THEN
   ICOVER=0
   DO JL2=1,NCOVER
     CALL PACK_2D_1D(IMASK,ZWORK3D(IIB:IIE,IJB:IJE,JL2),PFIELD(:,JL2))
   END DO
END IF
!
DEALLOCATE(ZWORK3D)
!
DEALLOCATE(IMASK)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX2COV_MNH
!
!     #############################################################
      SUBROUTINE READ_SURFX2COV_1COV_MNH(HREC,KL1,KCOVER,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX1* - routine to fill a real 2D array for the externalised surface
!!                 with Logical mask on one specified vertical level
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READ_SURFX1 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
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
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_FIELD,       ONLY : TFIELDDATA,TYPELOG,TYPEREAL
USE MODE_FM
USE MODE_FMREAD
USE MODE_ll
USE MODE_IO_ll
USE MODE_MSG
!
USE MODD_CST,         ONLY : XPI
!
USE MODD_IO_SURF_MNH, ONLY : TOUT, TPINFILE, NMASK, &
                             NIU, NJU, NIB, NJB, NIE, NJE, &
                             NIU_ALL, NJU_ALL, NIB_ALL,    &
                             NJB_ALL, NIE_ALL, NJE_ALL,    &
                             NMASK_ALL
!
USE MODI_PACK_2D_1D
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=LEN_HREC),INTENT(IN) :: HREC     ! name of the article to be read
INTEGER,                INTENT(IN) :: KL1      !  number of points
INTEGER,                INTENT(IN) :: KCOVER   ! index of the vertical level, it should be a index such that LCOVER(KCOVER)=.TRUE.
REAL, DIMENSION(KL1),   INTENT(OUT):: PFIELD   ! array containing the data field
INTEGER,                INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT):: HCOMMENT ! comment
CHARACTER(LEN=1),       INTENT(IN) :: HDIR     ! type of field :
!                                           ! 'H' for HOR : with hor. dim.; and  distributed.
!                                           ! 'A' for ALL : with hor. dim.; and not distributed.
!                                           ! '-' : no horizontal dim.

!
!*      0.2   Declarations of local variables
!
INTEGER           :: IGRID          ! IGRID : grid indicator
INTEGER           :: ILENCH         ! ILENCH : length of comment string
INTEGER           :: ILUOUT
!
CHARACTER(LEN=LEN_HREC) :: YREC
CHARACTER(LEN=2)  :: YDIR
CHARACTER(LEN=2)  :: YSTORAGE_TYPE
!
INTEGER           :: IIU, IJU, IIB, IJB, IIE, IJE ! dimensions of horizontal fields
INTEGER, DIMENSION(:), ALLOCATABLE :: IMASK       ! mask for packing
!JUANZ
INTEGER           :: NCOVER,ICOVER,JL2
REAL,DIMENSION(:,:), ALLOCATABLE :: ZWORK2D
!JUANZ
INTEGER  :: IVERSION, IBUGFIX
LOGICAL  :: GCOVER_PACKED ! .T. if COVER are all packed into one field
CHARACTER(LEN=1) :: YDIR1
TYPE(TFIELDDATA) :: TZFIELD
!-------------------------------------------------------------------------------
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','READ_SURFX2COV_1COV_MNH',TRIM(TPINFILE%CNAME)//': reading '//TRIM(HREC))
!
KRESP = 0
ILUOUT = TOUT%NLU
!YDIR1 = 'H'
!IF (PRESENT(HDIR)) YDIR1 = HDIR
YDIR1 = HDIR
!
IF (YDIR1=='A'.OR.YDIR1=='E') THEN
  YDIR="--"
  IIU = NIU_ALL
  IJU = NJU_ALL
  IIB = NIB_ALL
  IJB = NJB_ALL
  IIE = NIE_ALL
  IJE = NJE_ALL
  ALLOCATE(IMASK(SIZE(NMASK_ALL)))
  IMASK = NMASK_ALL
ELSE
  YDIR="XY"
  IIU = NIU
  IJU = NJU
  IIB = NIB
  IJB = NJB
  IIE = NIE
  IJE = NJE
  ALLOCATE(IMASK(SIZE(NMASK)))
  IMASK = NMASK
END IF
!
!! Reading of a 2D fields, masked and packed into 1D vector
!
!
ALLOCATE (ZWORK2D(IIU,IJU))
ZWORK2D(:,:) =  0.0
!
CALL IO_READ_FIELD(TPINFILE,'VERSION',IVERSION)
CALL IO_READ_FIELD(TPINFILE,'BUG',    IBUGFIX)

IF (IVERSION<7 .OR. (IVERSION==7 .AND. IBUGFIX==0)) THEN
  GCOVER_PACKED = .FALSE.
ELSE
  TZFIELD%CMNHNAME   = 'COVER_PACKED'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'COVER_PACKED'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = '--'
  TZFIELD%CCOMMENT   = ''
  TZFIELD%NGRID      = 0
  TZFIELD%NTYPE      = TYPELOG
  TZFIELD%NDIMS      = 0
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,GCOVER_PACKED,KRESP)
END IF
!
IF (.NOT. GCOVER_PACKED) THEN
  WRITE(YREC,'(A5,I3.3)') 'COVER',KCOVER
  TZFIELD%CMNHNAME   = TRIM(YREC)
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = TRIM(YREC)
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = YDIR1
  TZFIELD%CCOMMENT   = 'X_Y_'//TRIM(YREC)
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,ZWORK2D,KRESP)
ELSE
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'error : GCOVER_PACKED = ', GCOVER_PACKED, ' and we try to read the covers one by one '
  WRITE(ILUOUT,*) ' '
  CALL ABORT
END IF
!
IF (KRESP /=0) THEN
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'error when reading article ', HREC,'KRESP=',KRESP
  WRITE(ILUOUT,*) ' '
ELSE IF (YDIR1=='H' .OR. YDIR1=='A' .OR. YDIR1=='E') THEN
   CALL PACK_2D_1D(IMASK,ZWORK2D(IIB:IIE,IJB:IJE),PFIELD(:))
END IF
!
DEALLOCATE(ZWORK2D)


DEALLOCATE(IMASK)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX2COV_1COV_MNH
!
!     #############################################################
      SUBROUTINE READ_SURFN0_MNH(HREC,KFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READN0 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
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
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_ll
USE MODE_FIELD,       ONLY: TFIELDDATA,TYPEINT
USE MODE_FM
USE MODE_FMREAD
USE MODE_MSG
USE MODE_READ_SURF_MNH_TOOLS
!
USE MODD_IO_SURF_MNH, ONLY : TOUT, TPINFILE, NMASK, &
                             NIU, NJU, NIB, NJB, NIE, NJE
USE MODD_CONF,        ONLY : CPROGRAM
!
!
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=LEN_HREC),INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                INTENT(OUT) :: KFIELD   ! the integer to be read
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
INTEGER          :: IIMAX, IJMAX
INTEGER          :: ILUOUT
TYPE(TFIELDDATA) :: TZFIELD
!
!-------------------------------------------------------------------------------
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','READ_SURFN0_MNH',TRIM(TPINFILE%CNAME)//': reading '//TRIM(HREC))
!
KRESP=0
ILUOUT = TOUT%NLU
!
IF (HREC=='DIM_FULL' .AND. ( CPROGRAM=='IDEAL ' .OR.  &
                                  CPROGRAM=='SPAWN ' .OR. CPROGRAM=='ZOOMPG' ))THEN
   CALL IO_READ_FIELD(TPINFILE,'IMAX',IIMAX)
   CALL IO_READ_FIELD(TPINFILE,'JMAX',IJMAX)
   KFIELD = IIMAX * IJMAX
ELSE
   CALL PREPARE_METADATA_READ_SURF(HREC,'--',0,TYPEINT,0,'READ_SURFN0_MNH',TZFIELD)
   CALL IO_READ_FIELD(TPINFILE,TZFIELD,KFIELD,KRESP)

    IF (KRESP /=0) THEN
      WRITE(ILUOUT,*) 'WARNING'
      WRITE(ILUOUT,*) '-------'
      WRITE(ILUOUT,*) 'error when reading article ', HREC,'KRESP=',KRESP
      WRITE(ILUOUT,*) 'default value may be used, who knows???'
      WRITE(ILUOUT,*) ' '
   ENDIF

ENDIF
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFN0_MNH
!
!     #############################################################
      SUBROUTINE READ_SURFN1_MNH(HREC,KL,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READN0 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
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
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_FIELD,       ONLY: TFIELDDATA,TYPEINT
USE MODE_FM
USE MODE_FMREAD
USE MODE_MSG
USE MODE_READ_SURF_MNH_TOOLS
!
USE MODD_IO_SURF_MNH, ONLY : TOUT, TPINFILE, NMASK, &
                             NIU, NJU, NIB, NJB, NIE, NJE
!
USE MODI_PACK_2D_1D
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=LEN_HREC),INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                INTENT(IN)  :: KL       ! number of points
INTEGER, DIMENSION(KL), INTENT(OUT) :: KFIELD   ! the integer to be read
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
!                                               ! 'H' : field with
!                                               !       horizontal spatial dim.
!                                               ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
INTEGER           :: IGRID          ! IGRID : grid indicator
INTEGER           :: ILENCH         ! ILENCH : length of comment string
INTEGER           :: ILUOUT
!
INTEGER, DIMENSION(:,:), ALLOCATABLE :: IWORK  ! work array read in the file
TYPE(TFIELDDATA) :: TZFIELD
!---------------------------------------------------------------------
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','READ_SURFN1_MNH',TRIM(TPINFILE%CNAME)//': reading '//TRIM(HREC))
!
ILUOUT = TOUT%NLU
!
IF (HDIR=='-') THEN
!
  CALL PREPARE_METADATA_READ_SURF(HREC,'--',0,TYPEINT,1,'READ_SURFN1_MNH',TZFIELD)
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,KFIELD,KRESP)
!
ELSE IF (HDIR=='H') THEN
  ALLOCATE(IWORK(NIU,NJU))
!
  CALL PREPARE_METADATA_READ_SURF(HREC,'XY',4,TYPEINT,2,'READ_SURFN1_MNH',TZFIELD)
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,IWORK,KRESP)
!
 IF (KRESP /=0) THEN
    WRITE(ILUOUT,*) 'WARNING'
    WRITE(ILUOUT,*) '-------'
    WRITE(ILUOUT,*) 'error when reading article ', HREC,'KRESP=',KRESP
    WRITE(ILUOUT,*) 'default value may be used, who knows???'
    WRITE(ILUOUT,*) ' '
 ELSE
    CALL PACK_2D_1D(NMASK,IWORK(NIB:NIE,NJB:NJE),KFIELD)
 END IF
!
DEALLOCATE(IWORK)

ENDIF
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFN1_MNH
!
!     #############################################################
      SUBROUTINE READ_SURFC0_MNH(HREC,HFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READC0* - routine to read an integer
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READC0 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
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
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_ll
USE MODE_FIELD,       ONLY : TFIELDDATA,TYPECHAR
USE MODE_FMREAD
USE MODE_MSG
USE MODE_POS
USE MODE_READ_SURF_MNH_TOOLS
!
USE MODD_IO_SURF_MNH, ONLY : TOUT, TPINFILE
USE MODD_CONF,        ONLY : LCARTESIAN, CPROGRAM
USE MODD_LUNIT,       ONLY : TPGDFILE
!
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=LEN_HREC),INTENT(IN)  :: HREC      ! name of the article to be read
CHARACTER(LEN=40),      INTENT(OUT) :: HFIELD    ! the integer to be read
INTEGER,                INTENT(OUT) :: KRESP     ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT  ! comment
!
!*      0.2   Declarations of local variables
!
INTEGER           :: IRESP          ! return code
INTEGER           :: IGRID          ! IGRID : grid indicator
INTEGER           :: ILENCH         ! ILENCH : length of comment string
INTEGER           :: ILUOUT
!
INTEGER           :: ILUDES       ! .des file logical unit
!
LOGICAL           :: GFOUND
CHARACTER(LEN=4)  :: CTURB,CRAD,CGROUND,CCLOUD,CDCONV,CELEC
CHARACTER(LEN=6)  :: CSEA_FLUX
TYPE(TFIELDDATA)  :: TZFIELD
!
NAMELIST/NAM_PARAMn/CTURB,CRAD,CGROUND,CCLOUD,CDCONV,CSEA_FLUX, CELEC
!----------------------------------------------------------------------------
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','READ_SURFC0_MNH',TRIM(TPINFILE%CNAME)//': reading '//TRIM(HREC))
!
KRESP = 0
ILUOUT = TOUT%NLU
!
IF (TPINFILE%NMNHVERSION(1)<4 .OR. (TPINFILE%NMNHVERSION(1)==4 .AND. TPINFILE%NMNHVERSION(2)<6)) THEN
  SELECT CASE(TRIM(HREC))
    CASE('SNOW_VEG_TYPE')
      HFIELD='D95'
    CASE('SNOW_ROAD_TYPE')
      HFIELD='1-L'
    CASE('SNOW_ROOF_TYPE')
      HFIELD='1-L'
    CASE('PHOTO')
      HFIELD='NON'
    CASE('ISBA')
      HFIELD = '3-L'
    CASE('GRID_TYPE')
      IF (LCARTESIAN) THEN
        HFIELD="CARTESIAN "
      ELSE
        HFIELD='CONF PROJ '
      END IF
    CASE('NATURE','SEA','WATER','TOWN')
      IF (CPROGRAM=='REAL  ' .OR. CPROGRAM=='IDEAL ') THEN
        CGROUND='ISBA'
      ELSE
        CGROUND='NONE'
        ILUDES = TPINFILE%TDESFILE%NLU
        CALL POSNAM(ILUDES,'NAM_PARAMN',GFOUND,ILUOUT)
        IF (GFOUND) READ(UNIT=ILUDES,NML=NAM_PARAMn)
      END IF
      IF (CGROUND=='NONE') THEN
        HFIELD ='NONE  '
      ELSE IF (CGROUND=='FLUX') THEN
        HFIELD ='FLUX  '
      ELSE IF (CGROUND=='ISBA') THEN
        IF(HREC=='SEA   ') HFIELD ='SEAFLX'
        IF(HREC=='WATER ') HFIELD ='WATFLX'
        IF(HREC=='NATURE') HFIELD ='ISBA  '
        IF(HREC=='TOWN  ') HFIELD ='TEB   '
      ELSE
        CALL PRINT_MSG(NVERB_FATAL,'IO','READ_SURFC0_MNH',TRIM(TPINFILE%CNAME)//': error when reading article '//TRIM(HREC)// &
                       ' with CGROUND='//TRIM(CGROUND))
      END IF
    CASE DEFAULT
      CALL PREPARE_METADATA_READ_SURF(HREC,'--',0,TYPECHAR,0,'READ_SURFC0_MNH',TZFIELD)
      CALL IO_READ_FIELD(TPINFILE,TZFIELD,HFIELD,KRESP)
      !
      IF (KRESP /=0) THEN
        CALL PRINT_MSG(NVERB_FATAL,'IO','READ_SURFC0_MNH',TRIM(TPINFILE%CNAME)//': error when reading article '//TRIM(HREC)// &
                       ' default value may be used, who knows???')
      ENDIF
  END SELECT
ELSE IF ( HREC=='GRID_TYPE'.AND. ( &
                           (CPROGRAM=='IDEAL ' .AND. .NOT.ASSOCIATED(TPGDFILE,TOUT)) .OR. &
                           (CPROGRAM=='SPAWN ' .AND. .NOT.ASSOCIATED(TPGDFILE,TOUT)) .OR. &
                            CPROGRAM=='ZOOMPG'                         )) THEN
  IF (LCARTESIAN) THEN
    HFIELD="CARTESIAN "
  ELSE
    HFIELD='CONF PROJ '
  END IF
ELSE
  CALL PREPARE_METADATA_READ_SURF(HREC,'--',0,TYPECHAR,0,'READ_SURFC0_MNH',TZFIELD)
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,HFIELD,KRESP)
  !
  IF (KRESP /=0) THEN
        CALL PRINT_MSG(NVERB_FATAL,'IO','READ_SURFC0_MNH',TRIM(TPINFILE%CNAME)//': error when reading article '//TRIM(HREC)// &
                       ' default value may be used, who knows???')
  ENDIF
ENDIF
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFC0_MNH
!
!     #############################################################
      SUBROUTINE READ_SURFL1_MNH(HREC,KL,OFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READL1* - routine to read a logical array
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READL1 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
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
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_IO_SURF_MNH, ONLY : TOUT, TPINFILE, NMASK, &
                             NIU, NJU, NIB, NJB, NIE, NJE
!
USE MODE_FIELD,       ONLY : TFIELDDATA,TYPEINT,TYPELOG
USE MODE_FM
USE MODE_FMREAD
USE MODE_MSG
USE MODE_READ_SURF_MNH_TOOLS
!
USE MODI_PACK_2D_1D
!
!
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=LEN_HREC),INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                INTENT(IN)  :: KL       ! number of points
LOGICAL, DIMENSION(KL), INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
!                                               ! 'H' : field with
!                                               !       horizontal spatial dim.
!                                               ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
INTEGER           :: IGRID          ! IGRID : grid indicator
INTEGER           :: ILENCH         ! ILENCH : length of comment string
INTEGER           :: ILUOUT
LOGICAL, DIMENSION(:,:), ALLOCATABLE :: GWORK  ! work array read in the file
INTEGER, DIMENSION(:,:), ALLOCATABLE :: IWORK  ! work array read in the file
TYPE(TFIELDDATA)  :: TZFIELD
!-------------------------------------------------------------------------------
CALL PRINT_MSG(NVERB_DEBUG,'IO','READ_SURFL1_MNH',TRIM(TPINFILE%CNAME)//': reading '//TRIM(HREC))
!
ILUOUT = TOUT%NLU
!
IF (HDIR=='-') THEN
  CALL PREPARE_METADATA_READ_SURF(HREC,'--',0,TYPELOG,1,'READ_SURFL1_MNH',TZFIELD)
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,OFIELD,KRESP)

  IF (KRESP /=0) THEN
    WRITE(ILUOUT,*) 'WARNING'
    WRITE(ILUOUT,*) '-------'
    WRITE(ILUOUT,*) 'error when reading article ', HREC,'KRESP=',KRESP
    WRITE(ILUOUT,*) 'default value may be used, who knows???'
    WRITE(ILUOUT,*) ' '
  ENDIF
ELSE IF (HDIR=='H') THEN
  ALLOCATE(GWORK(NIU,NJU))
  GWORK = .FALSE.
!
  ALLOCATE(IWORK(NIU,NJU))
  CALL PREPARE_METADATA_READ_SURF(HREC,'XY',4,TYPEINT,2,'READ_SURFL1_MNH',TZFIELD)
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,IWORK,KRESP)
  WHERE (IWORK==1) GWORK = .TRUE.
  DEALLOCATE(IWORK)
!
  IF (KRESP /=0) THEN
    WRITE(ILUOUT,*) 'WARNING'
    WRITE(ILUOUT,*) '-------'
    WRITE(ILUOUT,*) 'error when reading article ', HREC,'KRESP=',KRESP
    WRITE(ILUOUT,*) 'default value may be used, who knows???'
    WRITE(ILUOUT,*) ' '
  ELSE
    CALL PACK_2D_1D(NMASK,GWORK(NIB:NIE,NJB:NJE),OFIELD)
  END IF
!
  DEALLOCATE(GWORK)
END IF
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFL1_MNH
!
!     #############################################################
      SUBROUTINE READ_SURFL0_MNH(HREC,OFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READL0* - routine to read a logical
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
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
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_FIELD,       ONLY : TFIELDDATA,TYPELOG
USE MODE_FM
USE MODE_FMREAD
USE MODE_MSG
USE MODE_READ_SURF_MNH_TOOLS
!
USE MODD_IO_SURF_MNH, ONLY : TOUT, TPINFILE
!
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=LEN_HREC),INTENT(IN)  :: HREC     ! name of the article to be read
LOGICAL,                INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
INTEGER           :: ILUOUT
TYPE(TFIELDDATA)  :: TZFIELD
!-------------------------------------------------------------------------------
CALL PRINT_MSG(NVERB_DEBUG,'IO','READ_SURFL0_MNH',TRIM(TPINFILE%CNAME)//': reading '//TRIM(HREC))
!
ILUOUT = TOUT%NLU
!
IF (HREC(1:4)=='BUDC') THEN
  IF (TPINFILE%NMNHVERSION(1)<4 .OR. (TPINFILE%NMNHVERSION(1)==4 .AND. TPINFILE%NMNHVERSION(2)<=5)) THEN
    OFIELD = .FALSE.
    KRESP = 0
    RETURN
  END IF
END IF
!
IF (HREC=='ECOCLIMAP') THEN
  IF (TPINFILE%NMNHVERSION(1)<4 .OR. (TPINFILE%NMNHVERSION(1)==4 .AND. TPINFILE%NMNHVERSION(2)<=6)) THEN
    OFIELD = .TRUE.
    KRESP = 0
    RETURN
  END IF
END IF
!
CALL PREPARE_METADATA_READ_SURF(HREC,'--',0,TYPELOG,0,'READ_SURFL0_MNH',TZFIELD)
CALL IO_READ_FIELD(TPINFILE,TZFIELD,OFIELD,KRESP)
HCOMMENT = TZFIELD%CCOMMENT
!
IF (KRESP /=0) THEN
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'error when reading article ', HREC,'KRESP=',KRESP
  WRITE(ILUOUT,*) 'default value may be used, who knows???'
  WRITE(ILUOUT,*) ' '
ENDIF
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFL0_MNH
!
!     #############################################################
      SUBROUTINE READ_SURFT0_MNH(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT0* - routine to read a MESO-NH date_time scalar
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READT0 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
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
USE MODE_FIELD,       ONLY : TFIELDDATA,TYPECHAR
USE MODE_FM
USE MODE_FMREAD
USE MODE_MSG
!
USE MODD_IO_SURF_MNH, ONLY : TOUT, TPINFILE
USE MODD_TYPE_DATE
!
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=LEN_HREC),INTENT(IN)    :: HREC     ! name of the article to be read
INTEGER,                INTENT(OUT)   :: KYEAR    ! year
INTEGER,                INTENT(OUT)   :: KMONTH   ! month
INTEGER,                INTENT(OUT)   :: KDAY     ! day
REAL,                   INTENT(OUT)   :: PTIME    ! time
INTEGER,                INTENT(OUT)   :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT)   :: HCOMMENT ! comment

!*      0.2   Declarations of local variables
!
INTEGER           :: IGRID          ! IGRID : grid indicator
INTEGER           :: ILENCH         ! ILENCH : length of comment string
INTEGER           :: ILUOUT
!
CHARACTER(LEN=LEN_HREC)        :: YRECFM     ! Name of the article to be written
CHARACTER(LEN=40)              :: YFILETYPE40! MESONH file type
CHARACTER(LEN=2)               :: YFILETYPE2 ! MESONH file type
INTEGER, DIMENSION(3)  :: ITDATE
TYPE(TFIELDDATA)       :: TZFIELD
TYPE(DATE_TIME)        :: TZDATETIME
!-------------------------------------------------------------------------------
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','READ_SURFT0_MNH',TRIM(TPINFILE%CNAME)//': reading '//TRIM(HREC))
!
ILUOUT = TOUT%NLU
!
IF (TPINFILE%NMNHVERSION(1)<4 .OR. (TPINFILE%NMNHVERSION(1)==4 .AND. TPINFILE%NMNHVERSION(2)<6)) THEN
  CALL IO_READ_FIELD(TPINFILE,'STORAGE_TYPE',YFILETYPE2)
ELSE
  TZFIELD%CMNHNAME   = 'STORAGETYPE'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'STORAGETYPE'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = '--'
  TZFIELD%CCOMMENT   = ''
  TZFIELD%NGRID      = 0
  TZFIELD%NTYPE      = TYPECHAR
  TZFIELD%NDIMS      = 0
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,YFILETYPE40)
  YFILETYPE2 = YFILETYPE40(1:2)
END IF
IF (YFILETYPE2(1:2)=='PG') THEN
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'Date is not read in a PGD file'
  WRITE(ILUOUT,*) 'Atmospheric model value is kept'
  WRITE(ILUOUT,*) ' '
  KRESP = -2
  RETURN
END IF
!
CALL IO_READ_FIELD(TPINFILE,HREC,TZDATETIME,KRESP)
!
IF (KRESP /=0) THEN
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'error when reading article ',YRECFM,'KRESP=',KRESP
  WRITE(ILUOUT,*) 'default value may be used, who knows???'
  WRITE(ILUOUT,*) ' '
ENDIF
!
KYEAR  = TZDATETIME%TDATE%YEAR
KMONTH = TZDATETIME%TDATE%MONTH
KDAY   = TZDATETIME%TDATE%DAY
PTIME  = TZDATETIME%TIME
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFT0_MNH

!     #############################################################
      SUBROUTINE READ_SURFT1_MNH(HREC,KL1,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT0* - routine to read a MESO-NH date_time vector
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READT1 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
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
!!      G. TANGUY      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     03/2009
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_FIELD,       ONLY : TFIELDDATA,TYPECHAR,TYPEINT,TYPEREAL
USE MODE_FM
USE MODE_FMREAD
USE MODE_MSG
!
USE MODD_IO_SURF_MNH, ONLY : TOUT, TPINFILE
!
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=LEN_HREC), INTENT(IN)    :: HREC     ! name of the article to be read
INTEGER,                 INTENT(IN)    :: KL1      ! number of points

INTEGER, DIMENSION(KL1), INTENT(OUT)   :: KYEAR    ! year
INTEGER, DIMENSION(KL1), INTENT(OUT)   :: KMONTH   ! month
INTEGER, DIMENSION(KL1), INTENT(OUT)   :: KDAY     ! day
REAL,    DIMENSION(KL1), INTENT(OUT)   :: PTIME    ! time
INTEGER,                 INTENT(OUT)   :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),      INTENT(OUT)   :: HCOMMENT ! comment

!*      0.2   Declarations of local variables
!
INTEGER           :: IGRID          ! IGRID : grid indicator
INTEGER           :: ILENCH         ! ILENCH : length of comment string
INTEGER           :: ILUOUT
!
CHARACTER(LEN=LEN_HREC)        :: YRECFM     ! Name of the article to be written
CHARACTER(LEN=40)              :: YFILETYPE40! MESONH file type
CHARACTER(LEN=2)               :: YFILETYPE2 ! MESONH file type
INTEGER, DIMENSION(3,KL1)  :: ITDATE
TYPE(TFIELDDATA)       :: TZFIELD
!-------------------------------------------------------------------------------
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','READ_SURFT1_MNH',TRIM(TPINFILE%CNAME)//': reading '//TRIM(HREC))
!
ILUOUT = TOUT%NLU
!
IF (TPINFILE%NMNHVERSION(1)<4 .OR. (TPINFILE%NMNHVERSION(1)==4 .AND. TPINFILE%NMNHVERSION(2)<6)) THEN
  CALL IO_READ_FIELD(TPINFILE,'STORAGE_TYPE',YFILETYPE2)
ELSE
  TZFIELD%CMNHNAME   = 'STORAGETYPE'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'STORAGETYPE'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = '--'
  TZFIELD%CCOMMENT   = ''
  TZFIELD%NGRID      = 0
  TZFIELD%NTYPE      = TYPECHAR
  TZFIELD%NDIMS      = 0
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_READ_FIELD(TPINFILE,TZFIELD,YFILETYPE40)
  YFILETYPE2 = YFILETYPE40(1:2)
END IF
!IF (YFILETYPE2(1:2)=='PG') THEN
!  WRITE(ILUOUT,*) 'WARNING'
!  WRITE(ILUOUT,*) '-------'
!  WRITE(ILUOUT,*) 'Date is not read in a PGD file'
!  WRITE(ILUOUT,*) 'Atmospheric model value is kept'
!  WRITE(ILUOUT,*) ' '
!  KRESP = -2
!  RETURN
!END IF
!
TZFIELD%CMNHNAME   = TRIM(HREC)//'%TDATE'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
TZFIELD%CUNITS     = ''
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = TRIM(HCOMMENT)
TZFIELD%NGRID      = 0
TZFIELD%NTYPE      = TYPEINT
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .FALSE.
!
CALL IO_READ_FIELD(TPINFILE,TZFIELD,ITDATE(:,:),KRESP)
!
KYEAR(:)  = ITDATE(1,:)
KMONTH(:) = ITDATE(2,:)
KDAY(:)   = ITDATE(3,:)
!
IF (KRESP /=0) THEN
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'error when reading article ',YRECFM,'KRESP=',KRESP
  WRITE(ILUOUT,*) 'default value may be used, who knows???'
  WRITE(ILUOUT,*) ' '
ENDIF
!
TZFIELD%CMNHNAME   = TRIM(HREC)//'%TIME'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
TZFIELD%CUNITS     = ''
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = TRIM(HCOMMENT)
TZFIELD%NGRID      = 0
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 1
TZFIELD%LTIMEDEP   = .FALSE.
!
CALL IO_READ_FIELD(TPINFILE,TZFIELD,PTIME(:),KRESP)
!
IF (KRESP /=0) THEN
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'error when reading article ',YRECFM,'KRESP=',KRESP
  WRITE(ILUOUT,*) 'default value may be used, who knows???'
  WRITE(ILUOUT,*) ' '
ENDIF
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFT1_MNH
