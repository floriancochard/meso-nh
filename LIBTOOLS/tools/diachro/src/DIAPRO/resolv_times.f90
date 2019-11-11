!     ######spl
      SUBROUTINE RESOLV_TIMES(K)
!     ##########################
!
!!****  *RESOLV_TIMES* -  Resolution des differentes dates du modele
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      None
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       16/01/95
!!      Updated   PM 
!-------------------------------------------------------------------------------
USE MODD_TIME
USE MODD_TYPE_AND_LH
USE MODD_GRID
USE MODD_CONF
USE MODD_TITLE
USE MODD_TIME1
USE MODD_ALLOC_FORDIACHRO
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1  Dummy argument
!
INTEGER          :: K
!
!*       0.1  local variables
!          

INTEGER          :: JJ
INTEGER           :: ITIM, IHOUR, IMINU, ISECD
CHARACTER(LEN=10) :: YTIM1, YTIM2
CHARACTER(LEN=LEN(CLEGEND)) :: YLEGEND
!
!-------------------------------------------------------------------------------
!
YLEGEND(1:LEN(YLEGEND))=' '
TDTEXP%TDATE%YEAR=XDATIME(1,K);  TDTEXP%TDATE%MONTH=XDATIME(2,K)
TDTEXP%TDATE%DAY=XDATIME(3,K);   TDTEXP%TIME=XDATIME(4,K)
TDTSEG%TDATE%YEAR=XDATIME(5,K);  TDTSEG%TDATE%MONTH=XDATIME(6,K)
TDTSEG%TDATE%DAY=XDATIME(7,K);   TDTSEG%TIME=XDATIME(8,K)
TDTMOD%TDATE%YEAR=XDATIME(9,K);  TDTMOD%TDATE%MONTH=XDATIME(10,K)
TDTMOD%TDATE%DAY=XDATIME(11,K);   TDTMOD%TIME=XDATIME(12,K)
TDTCUR%TDATE%YEAR=XDATIME(13,K);  TDTCUR%TDATE%MONTH=XDATIME(14,K)
TDTCUR%TDATE%DAY=XDATIME(15,K);   TDTCUR%TIME=XDATIME(16,K)

YTIM1='          '
YTIM2='          '
DO JJ=1,2
IF(JJ == 1)ITIM=TDTMOD%TIME
IF(JJ == 2)ITIM=TDTCUR%TIME
IHOUR=ITIM/3600
IMINU=(ITIM-IHOUR*3600)/60
ISECD=ITIM-(IHOUR*3600 + IMINU*60)
  IF(JJ == 1)THEN
    WRITE(YTIM1,'(I3,''H'',I2,''M'',I2,''S'')')IHOUR,IMINU,ISECD
  ELSE
    WRITE(YTIM2,'(I3,''H'',I2,''M'',I2,''S'')')IHOUR,IMINU,ISECD
  ENDIF
ENDDO

CLEGEND2(1:LEN(CLEGEND2))=' '
IF(CSTORAGE_TYPE /= 'PG')THEN
WRITE(CLEGEND2,1001)TDTMOD%TDATE%YEAR,TDTMOD%TDATE%MONTH,  &
                    TDTMOD%TDATE%DAY,YTIM1,                &
		    TDTCUR%TDATE%YEAR,TDTCUR%TDATE%MONTH,  &
		    TDTCUR%TDATE%DAY,YTIM2
ENDIF

YTIM1='          '
YTIM2='          '
DO JJ=1,2
IF(JJ == 1)ITIM=TDTEXP%TIME
IF(JJ == 2)ITIM=TDTSEG%TIME
IHOUR=ITIM/3600
IMINU=(ITIM-IHOUR*3600)/60
ISECD=ITIM-(IHOUR*3600 + IMINU*60)
  IF(JJ == 1)THEN
    WRITE(YTIM1,'(I3,''H'',I2,''M'',I2,''S'')')IHOUR,IMINU,ISECD
  ELSE
    WRITE(YTIM2,'(I3,''H'',I2,''M'',I2,''S'')')IHOUR,IMINU,ISECD
  ENDIF
ENDDO

YLEGEND=CLEGEND
CLEGEND(1:LEN(CLEGEND))=' '

SELECT CASE(CTYPE)

  CASE('CART','MASK')

    IF(CSTORAGE_TYPE /= 'PG')THEN
    IF(LCARTESIAN)WRITE(CLEGEND,1000)TDTEXP%TDATE%YEAR,TDTEXP%TDATE%MONTH,  &
    				 TDTEXP%TDATE%DAY,YTIM1,                &
                                     TDTSEG%TDATE%YEAR,TDTSEG%TDATE%MONTH,  &
    		                 TDTSEG%TDATE%DAY,YTIM2,                &
    			         'CARTESIEN                    '
    ENDIF
    IF(.NOT.LCARTESIAN)THEN
      IF(XRPK.EQ.0. .AND. CSTORAGE_TYPE /= 'PG')                              &
        WRITE(CLEGEND,1000)TDTEXP%TDATE%YEAR,TDTEXP%TDATE%MONTH,  &
                           TDTEXP%TDATE%DAY,YTIM1,                &
                           TDTSEG%TDATE%YEAR,TDTSEG%TDATE%MONTH,  &
                           TDTSEG%TDATE%DAY,YTIM2,                &
                           'MERCATOR                     '
      IF(XRPK.EQ.0. .AND. CSTORAGE_TYPE == 'PG')                              &
        WRITE(CLEGEND,'(2X,A29,79X)')'PROJECTION MERCATOR          '
      IF(ABS(XRPK).EQ.1. .AND. CSTORAGE_TYPE /= 'PG')                              &
        WRITE(CLEGEND,1000)TDTEXP%TDATE%YEAR,TDTEXP%TDATE%MONTH,  &
                           TDTEXP%TDATE%DAY,YTIM1,                &
                           TDTSEG%TDATE%YEAR,TDTSEG%TDATE%MONTH,  &
                           TDTSEG%TDATE%DAY,YTIM2,                &
    		       'STEREOG. POLAIRE             '
      IF(ABS(XRPK).EQ.1. .AND. CSTORAGE_TYPE == 'PG')                              &
        WRITE(CLEGEND,'(2X,A29,79X)')'PROJ. STEREOGRAPHIQUE POLAIRE'
      IF(ABS(XRPK).GT.0..AND.ABS(XRPK).LT.1. .AND. CSTORAGE_TYPE /= 'PG')                &
        WRITE(CLEGEND,1000)TDTEXP%TDATE%YEAR,TDTEXP%TDATE%MONTH,  &
                           TDTEXP%TDATE%DAY,YTIM1,                &
                           TDTSEG%TDATE%YEAR,TDTSEG%TDATE%MONTH,  &
                           TDTSEG%TDATE%DAY,YTIM2,                &
                           'LAMBERT                      '
      IF(ABS(XRPK).GT.0..AND.ABS(XRPK).LT.1. .AND. CSTORAGE_TYPE == 'PG')                &
        WRITE(CLEGEND,'(2X,A29,79X)')'PROJECTION LAMBERT           '
    END IF

  CASE DEFAULT

    WRITE(CLEGEND,1002)TDTEXP%TDATE%YEAR,TDTEXP%TDATE%MONTH,  &
                           TDTEXP%TDATE%DAY,YTIM1,                &
                           TDTSEG%TDATE%YEAR,TDTSEG%TDATE%MONTH,  &
                           TDTSEG%TDATE%DAY,YTIM2
END SELECT
CLEGEND(104:108)=YLEGEND(104:108)

1000 FORMAT('DATE EXP.. ',I4,2('/',I2),1X,A10,3X,         &
            'DATE SEG. ',I4,2('/',I2),1X,A10,3X,A26)
1001 FORMAT('DATE MOD. ',I4,2('/',I2),1X,A10,3X,         &
            'DATE CUR. ',I4,2('/',I2),1X,A10)
1002 FORMAT('DATE EXP.. ',I4,2('/',I2),1X,A10,3X,         &
            'DATE SEG. ',I4,2('/',I2),1X,A10)

!
!
RETURN
END SUBROUTINE  RESOLV_TIMES 
