PROGRAM testlibcomp
IMPLICIT NONE 

INTEGER, PARAMETER :: IDIMX = 2
INTEGER, PARAMETER :: IDIMY = 3
INTEGER, PARAMETER :: IDIMZ = 7

REAL(KIND=8),DIMENSION(IDIMX,IDIMY,IDIMZ) :: XORGTAB,XINTAB,XOUTTAB
INTEGER :: DATASIZE ! original size of array
INTEGER :: COMPSIZE ! size of compressed array
INTEGER :: JI
INTEGER :: INBELT,ITYPECOD

! Level 1 : constant level
XORGTAB(:,:,1) = -1.5

! Level 2 : 2 values in level
XORGTAB(:,:,2) = -10.4
XORGTAB(1,3,2) = -5.3

! Level 3 : 3 values in level
XORGTAB(:,:,3) = -8.2
XORGTAB(2,2,3) = 10.3
XORGTAB(1,3,3) = -9999.99

! Level 4 : normal
XORGTAB(:,:,4) = RESHAPE((/ (-(JI/1000.),JI=1,6) /),(/ IDIMX,IDIMY /))

! Level 5 : Min exclus
XORGTAB(:,:,5) = XORGTAB(:,:,4)
XORGTAB(2,1,5) = -5.5

! Level 6 : Max exclus
XORGTAB(:,:,6) = XORGTAB(:,:,4)
XORGTAB(2,2,6) = 10.8

! Level 7 : Min et Max exclus
XORGTAB(:,:,7) = XORGTAB(:,:,4)
XORGTAB(2,1,7) = -5.5
XORGTAB(2,2,7) = 10.8

XINTAB(:,:,:) = XORGTAB(:,:,:)
DATASIZE = IDIMX * IDIMY * IDIMZ
CALL COMPRESS_FIELD(XINTAB,IDIMX,IDIMY,DATASIZE,COMPSIZE)
PRINT *,"---> org size = ",DATASIZE,", comp size = ",COMPSIZE

! Now XINTAB is compressed
CALL  GET_COMPHEADER(XINTAB,DATASIZE,INBELT,ITYPECOD)
IF (INBELT /= DATASIZE) THEN
  PRINT *, "Fatal error in testlibcomp !"
  STOP
END IF
CALL DECOMPRESS_FIELD(XOUTTAB,DATASIZE,XINTAB,COMPSIZE,ITYPECOD)
! XOUTTAB contains the uncompressed data

DO JI=1,IDIMZ
  PRINT *,"Level ",JI
  PRINT *,"  Original    : ",XORGTAB(:,:,JI)
  PRINT *,"  comp/uncomp : ",XOUTTAB(:,:,JI)
  PRINT *,"  Difference  : ",XORGTAB(:,:,JI)-XOUTTAB(:,:,JI)
END DO

END PROGRAM testlibcomp
