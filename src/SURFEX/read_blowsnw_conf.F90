!     #########
      SUBROUTINE READ_BLOWSNW_CONF(HPROGRAM)
!     #######################################################
!
!!****  *READ_BLOWSNW_CONF* - routine to read the configuration for BLOWSNW
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Vionnet   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    Based on read_dst_conf.f90
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_POS_SURF
!
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
!
USE MODN_BLOWSNW
USE MODD_BLOWSNW_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6), INTENT(IN)  :: HPROGRAM ! program calling BLOWSNW

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
LOGICAL           :: GFOUND         ! Return code when searching namelist
INTEGER           :: ILUOUT         ! logical unit of output file
INTEGER           :: INAM           ! logical unit of namelist file
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!* get output listing file logical unit
!
IF (LHOOK) CALL DR_HOOK('READ_BLOWSNW_CONF',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!* open namelist file
!
 CALL OPEN_NAMELIST(HPROGRAM,INAM)
!
!* reading of namelist
!  -------------------
!
 CALL POSNAM(INAM,'NAM_SURF_BLOWSNW',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=INAM,NML=NAM_SURF_BLOWSNW)

! Check consistency between options for blowing snow scheme
IF (XEMIALPHA_SNW .NE. 3 .AND. CSNOW_SEDIM=='TABC') THEN
  WRITE(ILUOUT,*) '*****************************************'
  WRITE(ILUOUT,*) '* XEMIALPHA_SNW must be set to 3 when                '
  WRITE(ILUOUT,*) '* CSNOW_SEDIM = TABC                                 '
  WRITE(ILUOUT,*) '* Update the look-up table in BLOWSNW_SEDIM_LKT1D    '
  WRITE(ILUOUT,*) '* to use TABC with a different value of XEMIALPHA_SNW'
  WRITE(ILUOUT,*) '*****************************************'
  CALL ABOR1_SFX('BLOWSNW_CONF: inconsistency between XEMIALPHA_SNW and CSNOW_SEDIM')
END IF
!
!
!* close namelist file
!
 CALL CLOSE_NAMELIST(HPROGRAM,INAM)
IF (LHOOK) CALL DR_HOOK('READ_BLOWSNW_CONF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_BLOWSNW_CONF
