!ORILAM_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!!   ############################
     MODULE MODI_CH_AER_MOD_INIT
!!   ############################
!!
INTERFACE
SUBROUTINE CH_AER_MOD_INIT
END SUBROUTINE CH_AER_MOD_INIT
END INTERFACE
!!
END MODULE MODI_CH_AER_MOD_INIT
!!
!!
!!   ####################################
     SUBROUTINE CH_AER_MOD_INIT
!!   ####################################
!!
!!    PURPOSE
!!    -------
!!     initialize the aerosol module (to be called only once)
!!
!!    METHOD
!!    ------
!!
!!      allocate all arrays and initialize the basic variables (i.e. densities
!!    and molar weights)
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Vincent Crassier (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    20/03/03   P . Tulet (CNRM/GMEI)   add  initialization tabulation
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_AEROSOL
USE MODD_GLO
USE MODD_IO_ll,            ONLY: TFILEDATA
USE MODD_UNIFACPARAM
!
USE MODE_FM,               ONLY: IO_FILE_OPEN_ll,IO_FILE_CLOSE_ll
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_ADD2LIST
USE MODE_UNIFAC
!
!
!*       0.     DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
!*       0.2   Declarations of local variables
!
INTEGER, PARAMETER :: nc=22, nh=16, nt=11 ! inorganic interpolation
INTEGER            :: JI, JJ, JK, JL, JM ! loop counter
INTEGER :: IRESP                   ! return code in FM routines
INTEGER :: ILU                     ! logical unit
TYPE(TFILEDATA),POINTER :: TZFILE
!
!---------------------------------------------------------------------------
!
!
!
!        1.1    initialisation
!
TZFILE => NULL()
!
! Initialize the mineral tabulation 
IF (CMINERAL == 'NARES') THEN
!       .. the file ares.w contains the weights of the model
        CALL IO_FILE_ADD2LIST(TZFILE,'ares1A.w','CHEMTAB','READ')
        CALL IO_FILE_OPEN_ll(TZFILE)
        ILU = TZFILE%NLU
        READ(ILU,*) I1IA,J1JA,K1KA 
        DO JI=1,I1IA 
          READ(ILU,*) X1MAXA(1,JI),X1MINA(1,JI),X1MODA(1,JI) 
        ENDDO
        DO JI=1,K1KA 
          READ(ILU,*) X1MAXA(2,JI),X1MINA(2,JI),X1MODA(2,JI) 
        ENDDO
        DO JI=1,I1IA+1 
          READ(ILU,*) (W1IJA(JI,JJ),JJ=1,J1JA) 
        ENDDO
        DO JJ=1,J1JA+1 
          READ(ILU,*) (W1JKA(JJ,JK),JK=1,K1KA) 
        ENDDO
        CALL IO_FILE_CLOSE_ll(TZFILE)
        TZFILE => NULL()
        !
        CALL IO_FILE_ADD2LIST(TZFILE,'ares1C.w','CHEMTAB','READ')
        CALL IO_FILE_OPEN_ll(TZFILE)
        ILU = TZFILE%NLU
        READ(ILU,*) I1IC,J1JC,K1KC 
        DO JI=1,I1IC 
          READ(ILU,*) X1MAXC(1,JI),X1MINC(1,JI),X1MODC(1,JI) 
        ENDDO
        DO JI=1,K1KC 
          READ(ILU,*) X1MAXC(2,JI),X1MINC(2,JI),X1MODC(2,JI) 
        ENDDO
        DO JI=1,I1IC+1 
          READ(ILU,*) (W1IJC(JI,JJ),JJ=1,J1JC) 
        ENDDO
        DO JJ=1,J1JC+1 
          READ(ILU,*) (W1JKC(JJ,JK),JK=1,K1KC) 
        ENDDO
        CALL IO_FILE_CLOSE_ll(TZFILE)
        TZFILE => NULL()
        !
        CALL IO_FILE_ADD2LIST(TZFILE,'ares2A.w','CHEMTAB','READ')
        CALL IO_FILE_OPEN_ll(TZFILE)
        ILU = TZFILE%NLU
        READ(ILU,*) I2IA,J2JA,K2KA 
        DO JI=1,I2IA 
          READ(ILU,*) X2MAXA(1,JI),X2MINA(1,JI),X2MODA(1,JI) 
        ENDDO
        DO JI=1,K2KA 
          READ(ILU,*) X2MAXA(2,JI),X2MINA(2,JI),X2MODA(2,JI) 
        ENDDO
        DO JI=1,I2IA+1 
          READ(ILU,*) (W2IJA(JI,JJ),JJ=1,J2JA) 
        ENDDO
        DO JJ=1,J2JA+1 
          READ(ILU,*) (W2JKA(JJ,JK),JK=1,K2KA) 
        ENDDO
        CALL IO_FILE_CLOSE_ll(TZFILE)
        TZFILE => NULL()
        !
        CALL IO_FILE_ADD2LIST(TZFILE,'ares2B.w','CHEMTAB','READ')
        CALL IO_FILE_OPEN_ll(TZFILE)
        ILU = TZFILE%NLU
        READ(ILU,*) I2IB,J2JB,K2KB 
        DO JI=1,I2IB 
          READ(ILU,*) X2MAXB(1,JI),X2MINB(1,JI),X2MODB(1,JI) 
        ENDDO
        DO JI=1,K2KB 
          READ(ILU,*) X2MAXB(2,JI),X2MINB(2,JI),X2MODB(2,JI) 
        ENDDO
        DO JI=1,I2IB+1 
          READ(ILU,*) (W2IJB(JI,JJ),JJ=1,J2JB) 
        ENDDO
        DO JJ=1,J2JB+1 
          READ(ILU,*) (W2JKB(JJ,JK),JK=1,K2KB) 
        ENDDO
        CALL IO_FILE_CLOSE_ll(TZFILE)
        TZFILE => NULL()
        !
        CALL IO_FILE_ADD2LIST(TZFILE,'ares2C.w','CHEMTAB','READ')
        CALL IO_FILE_OPEN_ll(TZFILE)
        ILU = TZFILE%NLU
        READ(ILU,*) I2IC,J2JC,K2KC 
        DO JI=1,I2IC 
          READ(ILU,*) X2MAXC(1,JI),X2MINC(1,JI),X2MODC(1,JI) 
        ENDDO
        DO JI=1,K2KC 
          READ(ILU,*) X2MAXC(2,JI),X2MINC(2,JI),X2MODC(2,JI) 
        ENDDO
        DO JI=1,I2IC+1 
          READ(ILU,*) (W2IJC(JI,JJ),JJ=1,J2JC) 
        ENDDO
        DO JJ=1,J2JC+1 
          READ(ILU,*) (W2JKC(JJ,JK),JK=1,K2KC) 
        ENDDO
        CALL IO_FILE_CLOSE_ll(TZFILE)
        TZFILE => NULL()
        !
END IF
!
IF (CMINERAL == 'TABUL') THEN
  IF(.NOT.ALLOCATED(rhi)) ALLOCATE(rhi(16))
  IF(.NOT.ALLOCATED(tempi)) ALLOCATE(tempi(11))
  IF(.NOT.ALLOCATED(zsu)) ALLOCATE(zsu(22))
  IF(.NOT.ALLOCATED(znh)) ALLOCATE(znh(22))
  IF(.NOT.ALLOCATED(zni)) ALLOCATE(zni(22))
  IF(.NOT.ALLOCATED(zf)) ALLOCATE(zf(16,11,22,22,22,3))
  CALL IO_FILE_ADD2LIST(TZFILE,'AEROMIN_NEW','CHEMTAB','READ')
  CALL IO_FILE_OPEN_ll(TZFILE)
  ILU = TZFILE%NLU

  WRITE(*,*) 'LOADING MINERAL AEROSOL DATA ...'
  DO JI=1,nh
    READ(ILU,*) rhi(JI)
  ENDDO
  DO JI=1,nt
    READ(ILU,*) tempi(JI)
  ENDDO
  DO JI=1,nc
    READ(ILU,*) zsu(JI)
  ENDDO
  DO JI=1,nc
    READ(ILU,*) znh(JI)
  ENDDO
  DO JI=1,nc
    READ(ILU,*) zni(JI)
  ENDDO
  DO JI=1,nh
  DO JJ=1,nt
  DO JK=1,nc
  DO JL=1,nc
  DO JM=1,nc
    READ (ILU,*) zf(JI,JJ,JK,JL,JM,1:3)
  ENDDO
  ENDDO
  ENDDO
  ENDDO
  ENDDO
  WRITE(*,*) 'END LOADING'
  CALL IO_FILE_CLOSE_ll(TZFILE)
  TZFILE => NULL()
ENDIF

IF(TRIM(CORGANIC).eq."MPMPO")THEN
  !Set unifac coefficients for group a
  CALL AQ_UNIFAC_INI()
 
  !Set unifac coefficients for group b
  CALL ORG_UNIFAC_INI()
 
  !Calculate non time varying unifac stuff for aquous phase
  CALL UNIFAC_INI(   &
       QG_AQ         & !I [m2] surface of functional groups
       ,RG_AQ        & !I [m3] volume of functional groups
       ,NU_AQ        & !I [nbr] number of functional groups in molec
       ,THTAGP_AQ    & !O [frc] surface fraction of group (j) in molecule (i)
       ,Q_AQ         & !O [m2] surface of molecule
       ,R_AQ         & !O [m3] volume of molecule
       ,L_AQ         & !O [?] UNIFAC parameter for molecule
       ,NMOL_AQ      & !I [nbr] number of molecules used
       ,NFUNC_AQ     & !I [nbr] number of functional groups used
       )

  !Calculate non time varying unifac stuff for group organic phase
  CALL UNIFAC_INI(   &
       QG_ORG        & !I [m2] surface of functional groups
       ,RG_ORG       & !I [m3] volume of functional groups
       ,NU_ORG       & !I [nbr] number of functional groups in molec
       ,THTAGP_ORG   & !O [frc] surface fraction of group (j) in molecule (i)
       ,Q_ORG        & !O [m2] surface of molecule
       ,R_ORG        & !O [m3] volume of molecule
       ,L_ORG        & !O [?] UNIFAC parameter for molecule
       ,NMOL_ORG     & !I [nbr] number of molecules used
       ,NFUNC_ORG    & !I [nbr] number of functional groups used
       )
  
  !Set molality of solvent in binary mix with water at several RH
  CALL ZSR_INI_MPMPO()

ELSEIF(TRIM(CORGANIC).eq."PUN")THEN
   
   !Set Unifac coefficients for Pun's group A
   CALL AUNIFAC_INI
   
   !Set Unifac coefficients for Pun's group B
   CALL BUNIFAC_INI
   
   !Calculate non time varying unifac stuff for aquous phase
   CALL UNIFAC_INI(  &
        QG_A         & !I [m2] surface of functional groups
        ,RG_A        & !I [m3] volume of functional groups
        ,NU_A        & !I [nbr] number of functional groups in molec 
        ,THTAGP_A    & !O [frc] surface fraction of group (j) in molecule (i)
        ,Q_A         & !O [m2] surface of molecule
        ,R_A         & !O [m3] volume of molecule
        ,L_A         & !O [?] UNIFAC parameter for molecule
        ,NMOL_A      & !I [nbr] number of molecules used
        ,NFUNC_A     & !I [nbr] number of functional groups used
        )
   
   !Calculate non time varying unifac stuff for group organic phase
   CALL UNIFAC_INI(  &
        QG_B         & !I [m2] surface of functional groups
          ,RG_B      & !I [m3] volume of functional groups
          ,NU_B      & !I [nbr] number of functional groups in molec 
          ,THTAGP_B  & !O [frc] surface fraction of group (j) in molecule (i)
          ,Q_B       & !O [m2] surface of molecule
          ,R_B       & !O [m3] volume of molecule
          ,L_B       & !O [?] UNIFAC parameter for molecule
          ,NMOL_B    & !I [nbr] number of molecules used
          ,NFUNC_B   & !I [nbr] number of functional groups used
          )

   !Get zsr coefficients for pun's code
   CALL ZSR_INI_PUN()


ENDIF
!
!
END SUBROUTINE CH_AER_MOD_INIT
