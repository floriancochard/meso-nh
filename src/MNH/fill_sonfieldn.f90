!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_FILL_SONFIELD_n
!     ##########################
!
INTERFACE 
!
      SUBROUTINE FILL_SONFIELD_n(KMI,YFIELD,PNESTFIELD,KLSON)
!
INTEGER ,                 INTENT(IN)     :: KMI    ! son model number
CHARACTER(LEN=6),         INTENT(IN)     :: YFIELD ! name of the field to nest
REAL, DIMENSION(:,:,:,:), INTENT(INOUT)  :: PNESTFIELD
INTEGER,                  INTENT(IN)     :: KLSON  ! rank of son model in PNESTFIELD
!
END SUBROUTINE FILL_SONFIELD_n
END INTERFACE
!
END MODULE MODI_FILL_SONFIELD_n
!
!
!
!     ##################################################
      SUBROUTINE FILL_SONFIELD_n(KMI,YFIELD,PNESTFIELD,KLSON)
!     ##################################################
!
!!****  *FILL_SONFIELD_n* - fill the working array for nesting of pgd files
!!                          with        son model index= _n
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!       
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        27/09/96
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!        M.Moge        01/2016 bug fix for parallel execution
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_GRID_n
USE MODD_NESTING
USE MODD_PARAMETERS
USE MODE_SPLITTING_ll, ONLY : SPLIT2, DEF_SPLITTING2
USE MODD_VAR_ll, ONLY : NPROC, IP, YSPLITTING, NMNH_COMM_WORLD
USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
!
USE MODE_MODELN_HANDLER
!
!USE MODE_TOOLS_ll, ONLY : GET_OR_ll
!USE MODE_LS_ll
!USE MODD_LSFIELD_n, ONLY : SET_LSFIELD_1WAY_ll
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
INTEGER ,                 INTENT(IN)     :: KMI    ! son model number
CHARACTER(LEN=6),         INTENT(IN)     :: YFIELD ! name of the field to nest
REAL, DIMENSION(:,:,:,:), INTENT(INOUT)  :: PNESTFIELD
INTEGER,                  INTENT(IN)     :: KLSON  ! rank of son model in PNESTFIELD
!
!
!*       0.2   declarations of local variables
!
INTEGER :: IIB1,IIE1,IJB1,IJE1 ! limits of physical domain of KDAD model
INTEGER :: JI1,JJ1             ! loop counters   in domain of KDAD model
!
INTEGER :: JI2INF, JI2SUP      ! limits of a grid mesh of domain of KDAD model
INTEGER :: JJ2INF,JJ2SUP       ! relatively to son domain
INTEGER :: IMI                 ! current model index
INTEGER :: JLAYER              ! loop counter
INTEGER :: IINFO_ll
INTEGER :: IXSIZE, IYSIZE  ! sizes of global son domain in father grid
INTEGER :: IXSIZE_F, IYSIZE_F  ! sizes of global father domain
TYPE(ZONE_ll), DIMENSION(:), ALLOCATABLE :: TZSPLITTING
INTEGER :: IXOR, IYOR  ! origin of local subdomain
INTEGER :: IXOR_C, IYOR_C, IXEND_C, IYEND_C  ! origin and end of local physical son subdomain in father grid
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZSUM
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZSUM_C
INTEGER :: IDIMX_C, IDIMY_C ! size of extended local son subdomain in father grid obtained with GET_CHILD_DIM_ll
INTEGER :: IXDOMAINS, IYDOMAINS               ! number of subdomains in X and Y directions
LOGICAL :: GPREM                              ! needed for DEF_SPLITTING2, true if NPROC is a prime number
!-------------------------------------------------------------------------------
!
!*       1.    initializations
!              ---------------
!
IMI = GET_CURRENT_MODEL_INDEX()
CALL GET_OR_ll( YSPLITTING, IXOR, IYOR )
CALL GOTO_MODEL(KMI)
CALL GO_TOMODEL_ll(KMI, IINFO_ll)
!
IF (KLSON/=1) THEN
  ! get sizes of global son domain in father grid
  IXSIZE = NXEND_ALL(KMI) - NXOR_ALL (KMI) + 1 - 2*JPHEXT ! - 1
  IYSIZE = NYEND_ALL(KMI) - NYOR_ALL (KMI) + 1 - 2*JPHEXT ! - 1
  ! get splitting of current model KMI in father grid
  IXSIZE_F = NXEND_ALL(NDAD(KMI)) - NXOR_ALL (NDAD(KMI))  + 1 - 2*JPHEXT
  IYSIZE_F = NYEND_ALL(NDAD(KMI)) - NYOR_ALL (NDAD(KMI))  + 1 - 2*JPHEXT
  ALLOCATE(TZSPLITTING(NPROC))
! we want the same domain partitioning for the child domain and for the father domain
  CALL DEF_SPLITTING2(IXDOMAINS,IYDOMAINS,IXSIZE_F,IYSIZE_F,NPROC,GPREM)
  CALL SPLIT2 ( IXSIZE, IYSIZE, 1, NPROC, TZSPLITTING, YSPLITTING, IXDOMAINS, IYDOMAINS )
  IIB1 = JPHEXT + 1
  IIE1 = TZSPLITTING(IP)%NXEND - TZSPLITTING(IP)%NXOR + JPHEXT + 1
  IJB1 = JPHEXT + 1
  IJE1 = TZSPLITTING(IP)%NYEND - TZSPLITTING(IP)%NYOR + JPHEXT + 1
!  IIB1 = NXOR_ALL(KMI) + TZSPLITTING(IP)%NXOR - JPHEXT
!  IIE1 = NXOR_ALL(KMI) + TZSPLITTING(IP)%NXEND - JPHEXT
!  IJB1 = NYOR_ALL(KMI) + TZSPLITTING(IP)%NYOR - JPHEXT
!  IJE1 = NYOR_ALL(KMI) + TZSPLITTING(IP)%NYEND - JPHEXT
ENDIF
!
!* correct only if JPHEXT = 1
!
!JUAN A REVOIR TODO_JPHEXT !!!
! <<<<<<< fill_sonfieldn.f90
!IIB1 = NXOR_ALL (KMI)+1
!IIE1 = NXEND_ALL(KMI)-1
!IJB1 = NYOR_ALL (KMI)+1
!IJE1 = NYEND_ALL(KMI)-1
! =======
!IIB1 = NXOR_ALL (KMI)+JPHEXT
!IIE1 = NXEND_ALL(KMI)-JPHEXT
!IJB1 = NYOR_ALL (KMI)+JPHEXT
!IJE1 = NYEND_ALL(KMI)-JPHEXT
! >>>>>>> 1.2.4.1.18.2.2.1
!
DO JLAYER=1,SIZE(PNESTFIELD,4)
  PNESTFIELD(:,:,KLSON,JLAYER) = XUNDEF
END DO
!
!-------------------------------------------------------------------------------
IF (KLSON==1) THEN
!
!*       2.    case KLSON=1 : father itself
!              ----------------------------
!
      SELECT CASE(YFIELD)
        CASE ('ZS    ')
          PNESTFIELD(:,:,KLSON,1) = XZS(:,:)
         CASE ('ZSMT  ')   ! smooth topography for SLEVE coordinate
          PNESTFIELD(:,:,KLSON,1) = XZSMT(:,:)
        CASE DEFAULT
          CALL GOTO_MODEL(IMI)
          CALL GO_TOMODEL_ll(IMI, IINFO_ll)
      END SELECT
!
!-------------------------------------------------------------------------------
ELSE
!
!*       3.    case KLSON>1 : one son
!              ----------------------
!
!  ALLOCATE( ZSUM(SIZE(PNESTFIELD,1), SIZE(PNESTFIELD,2)) )
  ALLOCATE( ZSUM(SIZE(XZS,1), SIZE(XZS,2)) )
  !
  CALL GOTO_MODEL( NDAD(KMI) )
  CALL GO_TOMODEL_ll( NDAD(KMI), IINFO_ll )
  CALL GET_CHILD_DIM_ll(KMI, IDIMX_C, IDIMY_C, IINFO_ll)
  CALL GOTO_MODEL( KMI )
  CALL GO_TOMODEL_ll( KMI, IINFO_ll )
  ALLOCATE( ZSUM_C(IDIMX_C, IDIMY_C) )
  !
  DO JI1 = IIB1,IIE1
    DO JJ1 = IJB1,IJE1
      JI2INF= (JI1-IIB1)  *NDXRATIO_ALL(KMI)+1+JPHEXT
      JI2SUP= (JI1-IIB1+1)*NDXRATIO_ALL(KMI)  +JPHEXT
      JJ2INF= (JJ1-IJB1)  *NDYRATIO_ALL(KMI)+1+JPHEXT
      JJ2SUP= (JJ1-IJB1+1)*NDYRATIO_ALL(KMI)  +JPHEXT

      SELECT CASE(YFIELD)
         CASE ('ZS    ')
           ZSUM_C(1+JPHEXT+(JI1-IIB1+1),1+JPHEXT+(JJ1-IJB1+1)) = SUM ( XZS(JI2INF:JI2SUP,JJ2INF:JJ2SUP ) )&
                                     / ( NDXRATIO_ALL(KMI)*NDYRATIO_ALL(KMI) )
         CASE ('ZSMT  ')  ! smooth topography for SLEVE coordinate
           ZSUM_C(1+JPHEXT+(JI1-IIB1+1),1+JPHEXT+(JJ1-IJB1+1)) = SUM ( XZSMT(JI2INF:JI2SUP,JJ2INF:JJ2SUP ) )&
                                     / ( NDXRATIO_ALL(KMI)*NDYRATIO_ALL(KMI) )
        CASE DEFAULT
          CALL GOTO_MODEL(IMI)
          CALL GO_TOMODEL_ll(IMI, IINFO_ll)
          RETURN
      END SELECT

    END DO
  END DO
  !switch to father model to set the LSFIELD and do the communications with LS_FEEDBACK_ll
CALL GET_FEEDBACK_COORD_ll(IXOR_C,IYOR_C,IXEND_C,IYEND_C,IINFO_ll) ! physical domain's origin and end
  CALL SET_LSFIELD_2WAY_ll(PNESTFIELD(IXOR_C:IXEND_C,IYOR_C:IYEND_C,KLSON,1), ZSUM_C)
  CALL LS_FEEDBACK_ll(IINFO_ll)
  CALL UNSET_LSFIELD_1WAY_ll()
!
!-------------------------------------------------------------------------------
END IF
!
CALL GOTO_MODEL(IMI)
CALL GO_TOMODEL_ll(IMI, IINFO_ll)
!-------------------------------------------------------------------------------
!
END SUBROUTINE FILL_SONFIELD_n
