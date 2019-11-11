!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      MODULE MODE_GRIDCART
!     ####################
!
!!****  *MODE_GRIDCART* -  module routine SM_GRIDCART 
!!
!!    PURPOSE
!!    -------
!       The purpose of this executive module  is to package 
!     the routine SM_GRIDCART 
!    
!      
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!       NONE          
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/05/94 
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!--------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!-------------------------------------------------------------------------------
!
CONTAINS
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       1.   ROUTINE SM_GRIDCART
!             -------------------
!-------------------------------------------------------------------------------
!     #########################################################################
      SUBROUTINE SM_GRIDCART(PXHAT,PYHAT,PZHAT,PZS,OSLEVE,PLEN1,PLEN2,PZSMT,PDXHAT,PDYHAT,PZZ,PJ)
!     #########################################################################
!
!!****  *SM_GRIDCART * - routine to compute J 
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute the Jacobian (J) in the case
!     of a cartesian geometry 
!      
!
!!**  METHOD
!!    ------
!!       The height z is first determined, and then J is computed 
!!     
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : contains array border depths
!! 
!!        JPHEXT,JPVEXT : Arrays border zone depth
!!
!!      Module MODD_CONF       : contains  configuration variables for 
!!                               all models
!
!!        NVERB        : Listing verbosity
!!
!!    REFERENCE
!!    ---------
!!      Technical Specifications Report of the Meso-NH project (chapters 2 and 3)
!!
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/05/94 
!!      updated                 V. Ducrocq  *Meteo France*   27/06/94 
!!      Updated                 P.M.        *LA*             22/07/94
!!      Updated                 V. Ducrocq  *Meteo France*   23/08/94 
!!      Updated                 J. Escobar  *LA*             20/11/01
!!                                 rewrite "EOSHIFT" for portability (IBM-SP3 bug)
!!      Sleve coordinate        G. Zangler  *LA*             nov 2005
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST1D_ll
USE MODD_LUNIT_n,     ONLY : TLUOUT
!
USE MODD_PARAMETERS       
USE MODD_CONF
!
USE MODI_VERT_COORD
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:),     INTENT(IN)  :: PXHAT,PYHAT,PZHAT ! positions x,y,z in 
                                                         ! the cartesian plane
REAL, DIMENSION(:,:),   INTENT(IN)  :: PZS               ! orography
LOGICAL,                INTENT(IN)  :: OSLEVE            ! flag for SLEVE coordinate
REAL,                   INTENT(IN)  :: PLEN1             ! Decay scale for smooth topography
REAL,                   INTENT(IN)  :: PLEN2             ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:,:),   INTENT(IN)  :: PZSMT             ! smooth orography
!
REAL, DIMENSION(:),     INTENT(OUT) :: PDXHAT            ! meshlength in x 
                                                         ! direction
REAL, DIMENSION(:),     INTENT(OUT) :: PDYHAT            ! meshlength in y 
                                                         ! direction 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PZZ               ! Height z
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PJ                ! Jacobian of the
                                                         ! GCS transformation
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PXHAT,1),SIZE(PYHAT,1),SIZE(PZHAT,1)) :: ZDZ ! meshlength in
                                                                  ! z direction 
REAL, DIMENSION(SIZE(PZS,1),SIZE(PZS,2)) :: ZBOUNDZ          ! Extrapolated
REAL                                     :: ZBOUNDX,ZBOUNDY  ! value for the 
                                                             ! upper bounds in 
                                                             ! z,x,y directions  
!
TYPE(LIST1D_ll), POINTER :: TZHALO1_ll ! ! pointer for the list  of 1D fields
                                    ! that have  to be communicated immediatly
INTEGER :: IINFO_ll ! return code of the // routine
INTEGER      :: IIB,IJB,IKB      ! beginning of useful area of PXHAT,PYHAT,PZHAT  
INTEGER      :: IIE,IJE,IKE      ! end of useful area of PXHAT,PYHAT,PZHAT  
INTEGER      :: IIU,IJU,IKU      ! upper bounds of PXHAT,PYHAT,PZHAT  
INTEGER      :: IKLOOP           ! index for prints
INTEGER      :: ILUOUT,IRESP     ! logical unit number for prints, error code
!
INTEGER      :: JI,JJ,JK         ! loop index                         
!-------------------------------------------------------------------------------
!
!*       1    RETRIEVE LOGICAL UNIT NUMBERFOR OUTPUT-LISTING AND  DIMENSIONS 
!             --------------------------------------------------------------
!
ILUOUT = TLUOUT%NLU
!
IIU = UBOUND(PXHAT,1)         
IJU = UBOUND(PYHAT,1)        
IKU = UBOUND(PZHAT,1)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKE = IKU-JPVEXT
IKB = 1+JPVEXT
NULLIFY(TZHALO1_ll)
!
IF(NVERB >= 10) THEN                         ! Parameter checking
  WRITE(ILUOUT,*) 'SM_GRIDCART: IIU,IJU,IKU=',IIU,IJU,IKU
  WRITE(ILUOUT,*) 'SM_GRIDCART: IIE,IJE,IKE=',IIE,IJE,IKE
  WRITE(ILUOUT,*) 'SM_GRIDCART: IIB,IJB,IKB=',IIB,IJB,IKB
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTE Z
!              ---------
!
CALL VERT_COORD(OSLEVE,PZS,PZSMT,PLEN1,PLEN2,PZHAT,PZZ)
!
IF(NVERB >= 10) THEN                               !Value control
  WRITE(ILUOUT,*) 'SM_GRIDCART: Some PZS values:'
  WRITE(ILUOUT,*)  PZS(1,1),PZS(IIU/2,IJU/2),PZS(IIU,IJU)  
  WRITE(ILUOUT,*) 'SM_GRIDCART: Some PZZ values:'
  DO IKLOOP=1,IKU
    WRITE(ILUOUT,*) PZZ(1,1,IKLOOP),PZZ(IIU/2,IJU/2,IKLOOP), &
                    PZZ(IIU,IJU,IKLOOP)  
  END DO
ENDIF
!-------------------------------------------------------------------------------
!
!
!*       3.    COMPUTE J
!              ---------
!
ZBOUNDX      = 2.*PXHAT(IIU)   - PXHAT(IIU-1)
ZBOUNDY      = 2.*PYHAT(IJU)   - PYHAT(IJU-1)
ZBOUNDZ(:,:) = 2.*PZZ(:,:,IKU) - PZZ(:,:,IKU-1)
PDXHAT(:)  = EOSHIFT(PXHAT(:) ,1,ZBOUNDX)      - PXHAT(:)
PDYHAT(:)  = EOSHIFT(PYHAT(:) ,1,ZBOUNDY)      - PYHAT(:)
! update PDXHAT and PDYHAT 
CALL ADD1DFIELD_ll("X",TZHALO1_ll,PDXHAT)
CALL ADD1DFIELD_ll("Y",TZHALO1_ll,PDYHAT)
CALL UPDATE_1DHALO_ll(TZHALO1_ll,IINFO_ll)
DEALLOCATE(TZHALO1_ll)
!
ZDZ(:,:,1:IKU-1) =  PZZ(:,:,2:IKU) - PZZ(:,:,1:IKU-1)
ZDZ(:,:,IKU)     =  ZBOUNDZ(:,:)  - PZZ(:,:,IKU)
DO JK=1,IKU ; DO JJ=1,IJU ; DO JI=1,IIU
   PJ(JI,JJ,JK)  = PDXHAT(JI) * PDYHAT(JJ) * ZDZ(JI,JJ,JK) 
ENDDO ; ENDDO ; ENDDO
!
IF(NVERB >= 10) THEN                               !Value control
  WRITE(ILUOUT,*) 'Some PJ values:'
  DO IKLOOP=1,IKU
    WRITE(ILUOUT,*) PJ(1,1,IKLOOP),PJ(IIU/2,IJU/2,IKLOOP),  &
                    PJ(IIU,IJU,IKLOOP)  
  END DO
ENDIF
! 
!-------------------------------------------------------------------------------
!
END SUBROUTINE SM_GRIDCART
!-------------------------------------------------------------------------------
END MODULE MODE_GRIDCART
