!MNH_LIC Copyright 1995-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODI_READ_DUMMY_GR_FIELD_n
!
! 
!
INTERFACE
!
      SUBROUTINE READ_DUMMY_GR_FIELD_n(TPINIFILE,KIINF,KISUP,KJINF,KJSUP,OREAD_ALL)
!
USE MODD_IO_ll, ONLY : TFILEDATA
!
!*       0.1   declarations of arguments
!
TYPE(TFILEDATA),  INTENT(IN)  :: TPINIFILE    ! Initial file
INTEGER,          INTENT(IN)  :: KIINF,KISUP  ! Lower and upper Dimensions in x direction for working window
INTEGER,          INTENT(IN)  :: KJINF,KJSUP  ! Lower and upper dimensions in y direction for working window
LOGICAL,          INTENT(IN)  :: OREAD_ALL    ! Flag to read the entire 2D fields in the file.
!
END SUBROUTINE READ_DUMMY_GR_FIELD_n
!
!
!
END INTERFACE
!
!
!
END MODULE MODI_READ_DUMMY_GR_FIELD_n
!
!
!
!     #############################################################################
      SUBROUTINE READ_DUMMY_GR_FIELD_n(TPINIFILE,KIINF,KISUP,KJINF,KJSUP,OREAD_ALL)
!     #############################################################################
!
!!****  *READ_DUMMY_GR_FIELD_n* - routine to read dummy surface fields
!!
!!    PURPOSE
!!    -------
!       Allocates and initialize dummy surface fields
!       by reading their value in initial MNH file.
!
!!**  METHOD
!!    ------
!!    
!!    
!!
!!    EXTERNAL
!!    --------
!!      FMREAD   : to read data in LFIFM file
!!       
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_CONF   : contains configuration variables
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (routine READ_DUMMY_GR_FIELD)
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        2/05/95 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_DUMMY_GR_FIELD_n
USE MODE_FIELD,      ONLY : TFIELDDATA,TYPEINT,TYPEREAL
USE MODD_GRID_n
USE MODD_IO_ll,      ONLY : TFILEDATA
USE MODD_PARAMETERS, ONLY : JPHEXT, NMNHNAMELGTMAX
!
USE MODE_FMREAD
USE MODE_MSG
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
TYPE(TFILEDATA),  INTENT(IN)  :: TPINIFILE    ! Initial file
INTEGER,          INTENT(IN)  :: KIINF,KISUP  ! Lower and upper Dimensions in x direction for working window
INTEGER,          INTENT(IN)  :: KJINF,KJSUP  ! Lower and upper dimensions in y direction for working window
LOGICAL,          INTENT(IN)  :: OREAD_ALL    ! Flag to read the entire 2D fields in the file.
!
!*       0.2   declarations of local variables
!
INTEGER             :: IRESP       ! File management
CHARACTER (LEN=NMNHNAMELGTMAX) :: YRECFM ! variables
CHARACTER (LEN=20 ) :: YSTRING20   ! string
CHARACTER (LEN=3  ) :: YSTRING03   ! string
INTEGER             :: JDUMMY      ! Loop index for cover data
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWORK  ! work array read in the file
INTEGER                           :: IIWORK ! 1st dimension of work array
INTEGER                           :: IJWORK ! 2nd dimension of work array
INTEGER                           :: IIINF  ! lower I index
INTEGER                           :: IISUP  ! upper I index
INTEGER                           :: IJINF  ! lower J index
INTEGER                           :: IJSUP  ! upper J index
TYPE(TFIELDDATA)                  :: TZFIELD
!
!-------------------------------------------------------------------------------
!
!*       1..   READ DIMENSIONS IN THE FILE
!              ---------------------------
!
IF (OREAD_ALL) THEN
  !
  !* General case: the entire field is read
  !  --------------------------------------
  !
  IIINF = 1
  IISUP = SIZE(XXHAT)
  IJINF = 1
  IJSUP = SIZE(XYHAT)
  ALLOCATE(ZWORK(IISUP,IJSUP))
  !
ELSE
  !
  !* ONE processor ONLY. A small part of the field is read
  !  -----------------------------------------------------
  !
  IIINF = KIINF
  IISUP = KISUP
  IJINF = KJINF
  IJSUP = KJSUP
  !
  CALL IO_READ_FIELD(TPINIFILE,'IMAX',IIWORK)
  CALL IO_READ_FIELD(TPINIFILE,'JMAX',IJWORK)
  !
  ALLOCATE(ZWORK(IIWORK+2*JPHEXT,IJWORK+2*JPHEXT))
END IF
!
!-------------------------------------------------------------------------------
!
!*       2..    READ DUMMY VARIABLES
!              ----------------------
!
!
IF (TPINIFILE%NMNHVERSION(1)>=4) THEN
  TZFIELD%CMNHNAME   = 'DUMMY_GR_NBR'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'DUMMY_GR_NBR'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = '--'
  TZFIELD%CCOMMENT   = 'number of dummy pgd fields chosen by user'
  TZFIELD%NGRID      = 0
  TZFIELD%NTYPE      = TYPEINT
  TZFIELD%NDIMS      = 0
  TZFIELD%LTIMEDEP   = .FALSE.
  !
  CALL IO_READ_FIELD(TPINIFILE,TZFIELD,NDUMMY_GR_NBR,IRESP)
  !
  IF (IRESP/=0) THEN
    !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'IO','READ_DUMMY_GR_FIELD_n','problem when reading DUMMY_GR_NBR')
  ENDIF
ELSE
  NDUMMY_GR_NBR = 0
END IF
!
CDUMMY_GR_NAME(:) = '                     '
CDUMMY_GR_AREA(:) = '   '
!
IF (.NOT. ASSOCIATED(XDUMMY_GR_FIELDS)) &
ALLOCATE(XDUMMY_GR_FIELDS(SIZE(XXHAT),SIZE(XYHAT),NDUMMY_GR_NBR))
!
DO JDUMMY=1,NDUMMY_GR_NBR
  WRITE(YRECFM,'(A8,I3.3)') 'DUMMY_GR',JDUMMY
  TZFIELD%CMNHNAME   = TRIM(YRECFM)
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = TRIM(YRECFM)
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = 'XY'
  ! Expected comment is not known but is in the following form:
  ! 'X_Y_'//TRIM(YRECFM)//YSTRING20//YSTRING03
  TZFIELD%CCOMMENT   = ''
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  CALL IO_READ_FIELD(TPINIFILE,TZFIELD,ZWORK(:,:),IRESP)
  !
  IF (IRESP/=0) THEN
    !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'IO','READ_DUMMY_GR_FIELD_n','problem when reading '//TRIM(TZFIELD%CMNHNAME))
  ENDIF
  XDUMMY_GR_FIELDS(:,:,JDUMMY) = ZWORK(IIINF:IISUP,IJINF:IJSUP)
  !
  YSTRING20=TZFIELD%CCOMMENT( 4+LEN_TRIM(YRECFM)+1                : 4+LEN_TRIM(YRECFM)+LEN(YSTRING20) )
  YSTRING03=TZFIELD%CCOMMENT( 4+LEN_TRIM(YRECFM)+LEN(YSTRING20)+1 : 4+LEN_TRIM(YRECFM)+LEN(YSTRING20)+LEN(YSTRING03) )
  !
  CDUMMY_GR_NAME(JDUMMY) = YSTRING20
  CDUMMY_GR_AREA(JDUMMY) = YSTRING03
END DO
!
!-------------------------------------------------------------------------------
DEALLOCATE(ZWORK)
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_DUMMY_GR_FIELD_n
