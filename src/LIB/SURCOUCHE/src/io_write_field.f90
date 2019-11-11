!MNH_LIC Copyright 2016-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Original version:
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!
MODULE MODE_IO_WRITE_FIELD
!
USE MODD_IO_ll, ONLY: TOUTBAK
USE MODE_FIELD
USE MODE_FMWRIT
!
IMPLICIT NONE
!
CONTAINS
!
SUBROUTINE IO_WRITE_FIELDLIST(TPOUTPUT)
!
USE MODE_MODELN_HANDLER, ONLY : GET_CURRENT_MODEL_INDEX
!
IMPLICIT NONE
!
TYPE(TOUTBAK),    INTENT(IN)  :: TPOUTPUT !Output structure
!
INTEGER :: IDX
INTEGER :: IMI
INTEGER :: JI
!
IMI = GET_CURRENT_MODEL_INDEX()
!
DO JI = 1,SIZE(TPOUTPUT%NFIELDLIST)
  IDX = TPOUTPUT%NFIELDLIST(JI)
  SELECT CASE (TFIELDLIST(IDX)%NDIMS)
    !
    !0D output
    !
    CASE (0)
      SELECT CASE (TFIELDLIST(IDX)%NTYPE)
        !
        !0D real
        !
        CASE (TYPEREAL)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_X0D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X0D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_X0D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X0D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_X0D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 0D logical fields')
          END IF
        !
        !0D integer
        !
        CASE (TYPEINT)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_N0D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N0D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_N0D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N0D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_N0D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 0D integer fields')
          END IF
        !
        !0D logical
        !
        CASE (TYPELOG)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_L0D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_L0D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_L0D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_L0D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_L0D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 0D logical fields')
          END IF
        !
        !0D string
        !
        CASE (TYPECHAR)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_C0D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_C0D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_C0D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_C0D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_C0D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 0D character fields')
          END IF
        !
        !0D date/time
        !
        CASE (TYPEDATE)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_T0D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_T0D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_T0D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_T0D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_T0D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 0D date/time fields')
          END IF
        !
        !0D other types
        !
        CASE DEFAULT
          PRINT *,'FATAL: IO_WRITE_FIELDLIST: type not yet supported for 0D output of ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
          STOP
      END SELECT
    !
    !1D output
    !
    CASE (1)
      SELECT CASE (TFIELDLIST(IDX)%NTYPE)
        !
        !1D real
        !
        CASE (TYPEREAL)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_X1D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X1D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_X1D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X1D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_X1D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 1D real fields')
          END IF
!         !
!         !1D integer
!         !
!         CASE (TYPEINT)
!           IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_N1D) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N1D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_N1D(IMI)%DATA) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N1D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
!             CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_N1D(IMI)%DATA)
!           ELSE
!             CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 1D integer fields')
!           END IF
!         !
!         !1D logical
!         !
!         CASE (TYPELOG)
!           IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_L1D) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_L1D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_L1D(IMI)%DATA) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_L1D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
!             CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_L1D(IMI)%DATA)
!           ELSE
!             CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 1D logical fields')
!           END IF
!         !
!         !1D string
!         !
!         CASE (TYPECHAR)
!           IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_C1D) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_C1D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_C1D(IMI)%DATA) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_C1D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
!             CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_C1D(IMI)%DATA)
!           ELSE
!             CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 1D character fields')
!           END IF
        !
        !1D other types
        !
        CASE DEFAULT
          PRINT *,'FATAL: IO_WRITE_FIELDLIST: type not yet supported for 1D output of ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
          STOP
      END SELECT
    !
    !2D output
    !
    CASE (2)
      SELECT CASE (TFIELDLIST(IDX)%NTYPE)
        !
        !2D real
        !
        CASE (TYPEREAL)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_X2D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X2D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_X2D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X2D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_X2D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 2D real fields')
          END IF
        !
        !2D integer
        !
        CASE (TYPEINT)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_N2D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N2D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_N2D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N2D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_N2D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 2D integer fields')
          END IF
        !
        !2D other types
        !
        CASE DEFAULT
          PRINT *,'FATAL: IO_WRITE_FIELDLIST: type not yet supported for 2D output of ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
          STOP
      END SELECT
    !
    !3D output
    !
    CASE (3)
      SELECT CASE (TFIELDLIST(IDX)%NTYPE)
        !
        !3D real
        !
        CASE (TYPEREAL)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_X3D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X3D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_X3D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X3D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_X3D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not (yet) allowed for 3D real fields')
            !PW: TODO?: add missing field in TFIELDLIST?
            !CALL IO_WRITE_FIELD_LB(TPOUTPUT%TFILE,TFIELDLIST(IDX),***,TFIELDLIST(IDX)%TFIELD_X3D(IMI)%DATA)
          END IF
        !
        !3D integer
        !
        CASE (TYPEINT)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_N3D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N3D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_N3D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N3D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_N3D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not (yet) allowed for 3D integer fields')
            !PW: TODO?: add missing field in TFIELDLIST?
            !CALL IO_WRITE_FIELD_LB(TPOUTPUT%TFILE,TFIELDLIST(IDX),***,TFIELDLIST(IDX)%TFIELD_N3D(IMI)%DATA)
          END IF
        !
        !3D other types
        !
        CASE DEFAULT
          PRINT *,'FATAL: IO_WRITE_FIELDLIST: type not yet supported for 3D output of ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
          STOP
      END SELECT
    !
    !4D output
    !
    CASE (4)
      SELECT CASE (TFIELDLIST(IDX)%NTYPE)
        !
        !4D real
        !
        CASE (TYPEREAL)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_X4D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X4D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_X4D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X4D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_X4D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not (yet) allowed for 4D real fields')
            !PW: TODO?: add missing field in TFIELDLIST?
            !CALL IO_WRITE_FIELD_LB(TPOUTPUT%TFILE,TFIELDLIST(IDX),***,TFIELDLIST(IDX)%TFIELD_X4D(IMI)%DATA)
          END IF
        !
        !4D other types
        !
        CASE DEFAULT
          PRINT *,'FATAL: IO_WRITE_FIELDLIST: type not yet supported for 4D output of ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
          STOP
      END SELECT
!     !
!     !5D output
!     !
!     CASE (5)
!       SELECT CASE (TFIELDLIST(IDX)%NTYPE)
!         !
!         !5D real
!         !
!         CASE (TYPEREAL)
!           IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_X5D) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X5D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_X5D(IMI)%DATA) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X5D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
!             CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_X5D(IMI)%DATA)
!           ELSE
!             CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not (yet) allowed for 5D real fields')
!             !PW: TODO?: add missing field in TFIELDLIST?
!             !CALL IO_WRITE_FIELD_LB(TPOUTPUT%TFILE,TFIELDLIST(IDX),***,TFIELDLIST(IDX)%TFIELD_X5D(IMI)%DATA)
!           END IF
!         !
!         !5D other types
!         !
!         CASE DEFAULT
!           PRINT *,'FATAL: IO_WRITE_FIELDLIST: type not yet supported for 5D output of ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!           STOP
!       END SELECT
!     !
!     !6D output
!     !
!     CASE (6)
!       SELECT CASE (TFIELDLIST(IDX)%NTYPE)
!         !
!         !6D real
!         !
!         CASE (TYPEREAL)
!           IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_X6D) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X6D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_X6D(IMI)%DATA) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X6D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
!             CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_X6D(IMI)%DATA)
!           ELSE
!             CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not (yet) allowed for 6D real fields')
!             !PW: TODO?: add missing field in TFIELDLIST?
!             !CALL IO_WRITE_FIELD_LB(TPOUTPUT%TFILE,TFIELDLIST(IDX),***,TFIELDLIST(IDX)%TFIELD_X6D(IMI)%DATA)
!           END IF
!         !
!         !6D other types
!         !
!         CASE DEFAULT
!           PRINT *,'FATAL: IO_WRITE_FIELDLIST: type not yet supported for 4D output of ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!           STOP
!       END SELECT
    !
    !Other number of dimensions
    !
    CASE DEFAULT
      PRINT *,'FATAL: IO_WRITE_FIELDLIST: number of dimensions not yet supported for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
      STOP
  END SELECT
END DO
!
END SUBROUTINE IO_WRITE_FIELDLIST
!
!
!
SUBROUTINE IO_WRITE_FIELD_USER(TPOUTPUT)
!
#if 0
USE MODD_PARAMETERS, ONLY : JPVEXT
USE MODD_DYN_n,    ONLY: XTSTEP
USE MODD_FIELD_n,    ONLY: XUT, XVT, XRT, XTHT
USE MODD_PRECIP_n, ONLY: XINPRR
#endif
!
IMPLICIT NONE
!
TYPE(TOUTBAK),    INTENT(IN)  :: TPOUTPUT !Output structure
!
TYPE(TFIELDDATA) :: TZFIELD
!
#if 0
INTEGER          :: IKB
!
IKB=JPVEXT+1
!
TZFIELD%CMNHNAME   = 'UTLOW'
TZFIELD%CSTDNAME   = 'x_wind'
TZFIELD%CLONGNAME  = ''
TZFIELD%CUNITS     = 'm s-1'
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_Z_U component of wind at lowest physical level'
TZFIELD%NGRID      = 2
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XUT(:,:,IKB))
!
TZFIELD%CMNHNAME   = 'VTLOW'
TZFIELD%CSTDNAME   = 'y_wind'
TZFIELD%CLONGNAME  = ''
TZFIELD%CUNITS     = 'm s-1'
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_Z_V component of wind at lowest physical level'
TZFIELD%NGRID      = 3
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XVT(:,:,IKB))
!
TZFIELD%CMNHNAME   = 'THTLOW'
TZFIELD%CSTDNAME   = 'air_potential_temperature'
TZFIELD%CLONGNAME  = ''
TZFIELD%CUNITS     = 'K'
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_Z_potential temperature at lowest physical level'
TZFIELD%NGRID      = 1
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XTHT(:,:,IKB))
!
TZFIELD%CMNHNAME   = 'RVTLOW'
!TZFIELD%CSTDNAME   = 'humidity_mixing_ratio' !ratio of the mass of water vapor to the mass of dry air
TZFIELD%CSTDNAME   = 'specific_humidity'     !mass fraction of water vapor in (moist) air
TZFIELD%CLONGNAME  = ''
TZFIELD%CUNITS     = 'kg kg-1'
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_Z_Vapor mixing Ratio at lowest physical level'
TZFIELD%NGRID      = 1
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XRT(:,:,IKB,1))
!
TZFIELD%CMNHNAME   = 'ACPRRSTEP'
TZFIELD%CSTDNAME   = 'rainfall_amount'
TZFIELD%CLONGNAME  = ''
TZFIELD%CUNITS     = 'kg m-2'
TZFIELD%CDIR       = ''
TZFIELD%CCOMMENT   = 'X_Y_ACcumulated Precipitation Rain Rate during timestep'
TZFIELD%NGRID      = 1
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
!XACPRR is multiplied by 1000. to convert from m to kg m-2 (water density is assumed to be 1000 kg m-3)
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XINPRR*XTSTEP*1.0E3)
#endif
!
END SUBROUTINE IO_WRITE_FIELD_USER
!
END MODULE MODE_IO_WRITE_FIELD
