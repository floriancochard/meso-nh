!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_GET_SIZEX_LB
!     ##########################
!
INTERFACE
!
      SUBROUTINE GET_SIZEX_LB(KIMAX_ll,KJMAX_ll,KRIMX,               &
                              KISIZEXF,KJSIZEXF,KISIZEXFU,KJSIZEXFU, &
                              KISIZEX4,KJSIZEX4,KISIZEX2,KJSIZEX2    )
!
INTEGER,               INTENT(IN)   :: KIMAX_ll  !  Dimensions  in x direction 
                                                 ! of the physical domain,
INTEGER,               INTENT(IN)   :: KJMAX_ll  !  Dimensions  in y direction 
                                                 ! of the physical domain,
INTEGER, INTENT(IN)  :: KRIMX ! Size of the riming zone 
INTEGER, INTENT(OUT) :: KISIZEXF,KJSIZEXF ! size of the full LB ZONE (x direction)
INTEGER, INTENT(OUT) :: KISIZEXFU,KJSIZEXFU ! size of the full LB ZONE 
                                            ! for u-grid points (x direction)
INTEGER, INTENT(OUT) ::KISIZEX4,KJSIZEX4 ! size of the 4 grid points LB zone
INTEGER, INTENT(OUT) ::KISIZEX2,KJSIZEX2 ! size of the 2 grid points LB zone
END SUBROUTINE GET_SIZEX_LB
!
END INTERFACE
!
END MODULE MODI_GET_SIZEX_LB
!
!
!
!
!
!     #################################################################
      SUBROUTINE GET_SIZEX_LB(KIMAX_ll,KJMAX_ll,KRIMX,               &
                              KISIZEXF,KJSIZEXF,KISIZEXFU,KJSIZEXFU, &
                              KISIZEX4,KJSIZEX4,KISIZEX2,KJSIZEX2    )
!     ################################################################
!
!!****  *GET_SIZEX_LB* - routine to get the size of the West-East LB arrays
!!                       for each sub-domain
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to define the size of the West-East 
!!    LB arrays for each sub-domain. All the possible LB configurations are
!     considered: full LB for mass and u points, 4 points LB and 2 points LB  
!
!!**  METHOD
!!    ------
!!      Sequentially, the following sizes are initialized :
!!     - KISIZEXF,KJSIZEXF for full LB for mass points 
!!     and   KISIZEXFU,KJSIZEXFU for full LB for u points 
!!     - KISIZEX4,KJSIZEX4 for 4 points  LB 
!!     -  KISIZEXF2,KJSIZEXF2 for 2 points  LB 
!!     For one sub-domain, the   size of the eastern side is added
!!     to the size of  the western side; as the sub-domains are  rectangles, 
!!    the sizes of along North/South are the same. 
!!  
!!    EXTERNAL
!!    --------   
!!      Module MODE_TOOLS_ll : 
!!             GET_INTERSECTION_ll : to get the global indices of the LB 
!!     arrays for each sub-domain.
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!    Module MODD_PARAMETERS :  JPHEXT
!!
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
!!      Original    23/09/98
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!   J.Escobar 28/03/2019: for very small domain , force N/S/E/W check on getting LB bounds
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------ 
USE MODD_PARAMETERS
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
INTEGER,               INTENT(IN)   :: KIMAX_ll  !  Dimensions  in x direction 
                                                 ! of the physical domain,
INTEGER,               INTENT(IN)   :: KJMAX_ll  !  Dimensions  in y direction 
                                                 ! of the physical domain,
INTEGER, INTENT(IN)  :: KRIMX ! Size of the riming zone 
INTEGER, INTENT(OUT) :: KISIZEXF,KJSIZEXF ! size of the full LB ZONE (x direction)
INTEGER, INTENT(OUT) :: KISIZEXFU,KJSIZEXFU ! size of the full LB ZONE 
                                            ! for u-grid points (x direction)
INTEGER, INTENT(OUT) ::KISIZEX4,KJSIZEX4 ! size of the 4 grid points LB zone
INTEGER, INTENT(OUT) ::KISIZEX2,KJSIZEX2 ! size of the 2 grid points LB zone
!
!*       0.2   declarations of local variables
!
INTEGER :: IINFO ! return code of // routines
INTEGER :: IXOR, IYOR, IXEND, IYEND ! Coordinates of the LB zone
                                    ! (in global coordinates)
INTEGER :: IXORI, IYORI, IXENDI, IYENDI ! Coordinates of
                                        ! the intersection in local coordinates
!-------------------------------------------------------------------------------
!
!*       1.   INITIALIZATIONS
!              -------------- 
!   
KISIZEXF=0
KJSIZEXF=0
KISIZEXFU=0
KJSIZEXFU=0
KISIZEX4=0
KJSIZEX4=0
KISIZEX2=0
KJSIZEX2=0
!------------------------------------------------------------------------------
!
!*       2.    COMPUTE THE SIZES
!              -------------------------- 
!
!*       2.1 Full West-East LB zone 
!
IF (KRIMX /=0) THEN 
!   Western side, mass points   : 1:NRIMX+1,1:NJMAX_ll+ 2 * JPHEXT
  IXOR=1
  IXEND=KRIMX+JPHEXT ! +1
  IYOR=1 
  IYEND=KJMAX_ll+ 2 * JPHEXT
  CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,"EXTE",IINFO)
  IF (IINFO/=1 .AND. LWEST_ll() ) THEN  ! no empty intersection
    KISIZEXF=KISIZEXF + (IXENDI - IXORI +1)
    KJSIZEXF= IYENDI - IYORI +1 
  ENDIF
!   Eastern side , mass  and u-grid points:
!            NIMAX_ll+JPHEXT-NRIMX +1:NIMAX_ll+ 2 *JPHEXT,1:NJMAX_ll+2 *JPHEXT
  IXOR=KIMAX_ll + 2 * JPHEXT-KRIMX-JPHEXT+1 ! -KRIMX
  IXEND=KIMAX_ll+ 2 * JPHEXT
  IYOR=1
  IYEND=KJMAX_ll+ 2 * JPHEXT
  CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,"EXTE",IINFO)
  IF (IINFO/=1 .AND. LEAST_ll() ) THEN ! no empty intersection
    KISIZEXF=KISIZEXF + (IXENDI - IXORI +1) ! added to the western side 
    KJSIZEXF= IYENDI - IYORI +1 
    KISIZEXFU=KISIZEXFU + ( IXENDI - IXORI +1)
    KJSIZEXFU= IYENDI - IYORI +1 
  ENDIF
! Western side, u-grid points : 2:NRIMX+2,1:NJMAX_ll+ 2 * JPHEXT
  IXOR=2
  IXEND=KRIMX+JPHEXT+1 ! +2
  IYOR=1 
  IYEND=KJMAX_ll+ 2 * JPHEXT
  CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,"EXTE",IINFO)
  IF (IINFO/=1 .AND. LWEST_ll() ) THEN ! no empty intersection
    KISIZEXFU=KISIZEXFU + (IXENDI - IXORI +1)
    KJSIZEXFU= IYENDI - IYORI +1 
  ENDIF
ENDIF
!
!*       2.2  West-East LB zone with only  2 points at each side 
!
!   Western side : 2:3,1:NJMAX_ll+ 2 * JPHEXT
IXOR=2 ! 2
IXEND=JPHEXT+2 ! 3
IYOR=1
IYEND=KJMAX_ll+ 2 * JPHEXT
CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,"EXTE",IINFO)
IF (IINFO /=1 .AND. LWEST_ll() ) THEN ! no empty intersection
  KISIZEX4=KISIZEX4 + ( IXENDI - IXORI +1)
  KJSIZEX4= IYENDI - IYORI +1 
ENDIF
!   Eastern side,  : NIMAX_ll+JPHEXT:NIMAX_ll+ 2 *JPHEXT,1:NJMAX_ll+2 *JPHEXT
IXOR=KIMAX_ll + 2 * JPHEXT - JPHEXT      ! + JPHEXT 
IXEND=KIMAX_ll+ 2 * JPHEXT - JPHEXT + JPHEXT  ! + 2*JPHEXT
IYOR=1
IYEND=KJMAX_ll+ 2 * JPHEXT
CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,"EXTE",IINFO)
IF (IINFO/=1 .AND. LEAST_ll() ) THEN ! no empty intersection
  KISIZEX4=KISIZEX4 + (IXENDI - IXORI +1)
  KJSIZEX4= IYENDI - IYORI +1 
ENDIF
!
!*       2.3  West-East LB zone with only  1 point at each side 
!
!   Western side : 1:1,1:NJMAX_ll+ 2 * JPHEXT
IXOR=1  ! 1
IXEND=JPHEXT ! 1
IYOR=1
IYEND=KJMAX_ll+ 2 * JPHEXT
CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,"EXTE",IINFO)
IF (IINFO/=1 .AND. LWEST_ll() ) THEN   ! no empty intersection
  KISIZEX2=KISIZEX2 + ( IXENDI - IXORI +1)
  KJSIZEX2= IYENDI - IYORI +1 
ENDIF
!   East boundary, 1 point LB zone : NIMAX_ll+ 2 * JPHEXT:NIMAX_ll+ 2 *JPHEXT,1:NJMAX_ll+2 *JPHEXT
IXOR=KIMAX_ll  + 2 * JPHEXT - JPHEXT + 1 !  + 2 * JPHEXT
IXEND=KIMAX_ll + 2 * JPHEXT - JPHEXT + JPHEXT !  + 2 * JPHEXT
IYOR=1
IYEND=KJMAX_ll+ 2 * JPHEXT
CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,"EXTE",IINFO)
IF (IINFO/=1 .AND. LEAST_ll() ) THEN    ! no empty intersection
  KISIZEX2=KISIZEX2 + ( IXENDI - IXORI +1)
  KJSIZEX2= IYENDI - IYORI +1 
ENDIF
!
END SUBROUTINE GET_SIZEX_LB
   
