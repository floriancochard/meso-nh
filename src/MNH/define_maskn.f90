!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ########################
      MODULE MODI_DEFINE_MASK_n
!     ########################
!
!
INTERFACE 
!
      SUBROUTINE DEFINE_MASK_n
!
END SUBROUTINE DEFINE_MASK_n
!
END INTERFACE
!
END MODULE MODI_DEFINE_MASK_n
!
!
!     #########################
      SUBROUTINE DEFINE_MASK_n
!     #########################
!
!!****  *DEFINE_MASK_n* - allocates arrays for nesting of pgds
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
!!      Original        26/09/96
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!   P. Wautelet : 14/04/2016 : better way to get father coordinates
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_DIM_n
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_LUNIT
USE MODD_NESTING
USE MODD_NEST_PGD_n
!
USE MODE_FM
USE MODE_IO_ll
USE MODE_MSG
USE MODE_MODELN_HANDLER
!
USE MODD_VAR_ll, ONLY : YSPLITTING
USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
USE MODE_TOOLS_ll, ONLY : INTERSECTION
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
!*       0.2   declarations of local variables
!
INTEGER :: ILUOUT0
INTEGER :: IRESP
INTEGER :: ISON
INTEGER :: JLOOP
INTEGER :: IMI
INTEGER     :: IXOR_F, IYOR_F    ! origin of local father subdomain (global coord)
INTEGER     :: IXEND_F, IYEND_F    ! end of local father subdomain (global coord)
INTEGER     :: IXOR_C, IYOR_C    ! origin of intersection between son model and local father subdomain (global coord)
INTEGER     :: IXEND_C, IYEND_C    ! end of intersection between son model and local father subdomain (global coord)
TYPE(ZONE_ll), DIMENSION(1) :: TZSPLITTING
TYPE(ZONE_ll) :: TZCOARSESONGLB ! global son domain in father grid
TYPE(ZONE_ll), DIMENSION(1) :: TZCOARSESONLCL ! intersection of global son domain and local father subdomain
!-------------------------------------------------------------------------------
!
ILUOUT0 = TLUOUT0%NLU
IMI=GET_CURRENT_MODEL_INDEX()
!
ALLOCATE ( NNESTMASK (NIMAX+2*JPHEXT,NJMAX+2*JPHEXT,1+COUNT(NDAD(:)==IMI)))
ALLOCATE ( NSON      (                              1+COUNT(NDAD(:)==IMI)))
!
! get splitting of father model
!
CALL GET_OR_ll(YSPLITTING(1:1),IXOR_F,IYOR_F)
IXEND_F = IXOR_F+NIMAX-1
IYEND_F = IYOR_F+NJMAX-1
!
TZSPLITTING(1)%NXOR = IXOR_F+JPHEXT
TZSPLITTING(1)%NYOR = IYOR_F+JPHEXT
TZSPLITTING(1)%NXEND = IXEND_F+JPHEXT
TZSPLITTING(1)%NYEND = IYEND_F+JPHEXT
!
TZCOARSESONGLB%NZOR = TZSPLITTING(1)%NZOR    ! there is no splitting in Z direction
TZCOARSESONGLB%NZEND = TZSPLITTING(1)%NZEND  ! there is no splitting in Z direction
TZCOARSESONGLB%NUMBER = TZSPLITTING(1)%NUMBER
!
!
NNESTMASK(:,:,:) = 0
NSON(1) = IMI
!
ISON=1
DO JLOOP=1,NMODEL
  IF (NDAD(JLOOP)/=IMI) CYCLE
  ISON=ISON+1
  NSON(ISON)=JLOOP
  !
  ! init global son zone in father grid coords
  !
  TZCOARSESONGLB%NXOR = NXOR_ALL(JLOOP)+JPHEXT
  TZCOARSESONGLB%NYOR = NYOR_ALL(JLOOP)+JPHEXT
  TZCOARSESONGLB%NXEND = NXEND_ALL(JLOOP)-JPHEXT
  TZCOARSESONGLB%NYEND = NYEND_ALL(JLOOP)-JPHEXT
  ! get the intersection  with local father subdomain -> TZCOARSESONLCL
  CALL INTERSECTION( TZSPLITTING, 1, TZCOARSESONGLB, TZCOARSESONLCL)
  IXOR_C = TZCOARSESONLCL(1)%NXOR
  IXEND_C = TZCOARSESONLCL(1)%NXEND
  IYOR_C = TZCOARSESONLCL(1)%NYOR
  IYEND_C = TZCOARSESONLCL(1)%NYEND
  IF ( IXEND_C/=0 .AND. IYEND_C/=0 ) THEN
    ! the intersection is non empty
    NNESTMASK( (IXOR_C-IXOR_F+1):(IXEND_C-IXOR_F+1), (IYOR_C-IYOR_F+1):(IYEND_C-IYOR_F+1), ISON) = 1
  ENDIF
END DO
!
IF (ANY (SUM(NNESTMASK(:,:,:),DIM=3)>1) ) THEN
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','DEFINE_MASK_n','two nested models with the same father are overlapping each other')
END IF
!
NNESTMASK(:,:,1) = 1.-SUM(NNESTMASK(:,:,:),DIM=3)
!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFINE_MASK_n
