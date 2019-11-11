!      ############
       PROGRAM FMMORE
!      ############
!
!!****  *FMMORE* - routine to list the content of a LFI file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMMORE is to list the content of a LFI file
!
!!**  METHOD
!!    ------
!!
!!      The FM and LFI routines are used to open, list and close the LFI file
!!    This routine is embedded in a Unix shell script to mimic the "more"
!!    function.
!!
!!    EXTERNAL
!!    --------
!!
!!      FMOPEN, FMLOOK, LFINAF, LFILAF, FMCLOS
!!
!!    calls: READUNTOUCH containing FMREAD
!!
!!    REFERENCE
!!    ---------
!!
!!      The structure and content of the Meso-NH files (C. Fischer)
!!
!!    AUTHOR
!!    ------
!!
!!      C. FISCHER      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                        03/95
!!      new I/O      (Mallet)                                           03/02
!!
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
#ifdef NAGf95
  USE F90_UNIX
#endif
!
! en attendant une Surcouche officielle...
!USE MODE_FM  
!
IMPLICIT NONE
!
!*      0.2   Declarations of local variables
!
INTEGER :: krep
INTEGER :: KNPRAR, KFTYPE,KVERB,KNINAR,KNUMBR
INTEGER :: KNALDO, KNTROU, KNARES, KNAMAX
LOGICAL :: LDTOUT
CHARACTER(LEN=32) :: CLUOUT,YLFINAME
CHARACTER(LEN=28) :: CFNAME
! reading of filename as input argument
#ifndef NAGf95
INTEGER :: IARGC
! CRAY specific
INTEGER :: arglen
!!!!!!!!!!!!!!!!!
#endif
INTEGER :: inarg,iresp
CHARACTER(LEN=50) :: yexe
!
!*      1.    INITIALIZATION
!             --------------
!
KFTYPE=2      ! pas de transfert dans fmclos
KVERB=0
!
CLUOUT='output_listing'
!
knaldo=0 ; kntrou=0 ; knares=0 ; knamax=0
LDTOUT=.TRUE.
!
!*      2.    READING FILENAME
!             ----------------
!READ(5,FMT='(A28)') CFNAME
INARG = IARGC()

#if defined(F90HP)
#define HPINCR 1
#else
#define HPINCR 0
#endif

#if defined(FUJI) || defined(NAGf95) || defined(NEC) || defined(HP) || defined(pgf) || defined(G95) || defined(GFORTRAN)
CALL GETARG(0+HPINCR,yexe)
IF (LEN_TRIM(yexe) == 0) THEN
  PRINT *, 'FATAL ERROR : Recompiler avec la macro -DF90HP'
  STOP
END IF
#else
CALL PXFGETARG(0,yexe,arglen,iresp)
#endif
!  PRINT *,yexe, ' avec ',INARG,' arguments.'
IF (INARG == 1) THEN 
#if defined(FUJI) || defined(NAGf95) || defined(NEC) || defined(HP) || defined(pgf) || defined(G95) || defined(GFORTRAN)
     CALL GETARG(1+HPINCR,CFNAME)
#else
     CALL PXFGETARG(1,CFNAME,arglen,iresp)
#endif
ELSE 
  PRINT *,'Usage : ', TRIM(yexe), ' [fichier fm]'
  STOP
END IF
!
!*      3.    OPENING FILE
!             ------------

! en attendant une Surcouche officielle...
!CALL FMOPEN_ll(CFNAME,'READ',CLUOUT,KNPRAR,KFTYPE,KVERB,&
CALL FMOPEN(CFNAME,'OLD',CLUOUT,KNPRAR,KFTYPE,KVERB,&
                        KNINAR,krep)
IF (krep.NE.0) GOTO 1000
!
!*      4.    
!
YLFINAME=ADJUSTL(ADJUSTR(CFNAME)//'.lfi')
! en attendant une Surcouche officielle...
!CALL FMLOOK_ll(YLFINAME,CLUOUT,knumbr,krep)
CALL FMLOOK(YLFINAME,CLUOUT,knumbr,krep)
IF (krep.NE.0) GOTO 1000
CALL LFINAF(krep,knumbr,knaldo,kntrou,knares,knamax)
IF (krep.NE.0) GOTO 1000
!WRITE(6,*) knaldo,kntrou,knares,knamax
IF (krep.NE.0) GOTO 1000
CALL LFILAF(krep,knumbr,LDTOUT)
!
CALL READUNTOUCH(CFNAME,CLUOUT)
!
! en attendant une Surcouche officielle...
!CALL FMCLOS_ll(CFNAME,'KEEP',CLUOUT,krep)
CALL FMCLOS(CFNAME,'KEEP',CLUOUT,krep)
IF (krep.NE.0) THEN
  GOTO 1000
ELSE
  GOTO 1010
ENDIF
!
1000   WRITE (0,*) ' exit in FMMORE with :',krep
1010   CONTINUE
!
END PROGRAM
