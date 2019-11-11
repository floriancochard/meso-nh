!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_GET_SIZEY_LB
!     ##########################
!
INTERFACE
!
      SUBROUTINE GET_SIZEY_LB(KIMAX_ll,KJMAX_ll,KRIMY,               &
                              KISIZEYF,KJSIZEYF,KISIZEYFV,KJSIZEYFV, &
                              KISIZEY4,KJSIZEY4,KISIZEY2,KJSIZEY2    )
!
INTEGER,               INTENT(IN)   :: KIMAX_ll  !  Dimensions  in x direction 
                                                 ! of the physical domain,
INTEGER,               INTENT(IN)   :: KJMAX_ll  !  Dimensions  in y direction 
                                                 ! of the physical domain,
INTEGER, INTENT(IN)  :: KRIMY ! Size of the riming zone 
INTEGER, INTENT(OUT) :: KISIZEYF,KJSIZEYF ! size of the full LB ZONE (x direction)
INTEGER, INTENT(OUT) :: KISIZEYFV,KJSIZEYFV ! size of the full LB ZONE 
                                            ! for u-grid points (x direction)
INTEGER, INTENT(OUT) ::KISIZEY4,KJSIZEY4 ! size of the 4 grid points LB zone
INTEGER, INTENT(OUT) ::KISIZEY2,KJSIZEY2 ! size of the 2 grid points LB zone
END SUBROUTINE GET_SIZEY_LB
!
END INTERFACE
!
END MODULE MODI_GET_SIZEY_LB
!
!
!
!
!
!     ################################################################
      SUBROUTINE GET_SIZEY_LB(KIMAX_ll,KJMAX_ll,KRIMY,               &
                              KISIZEYF,KJSIZEYF,KISIZEYFV,KJSIZEYFV, &
                              KISIZEY4,KJSIZEY4,KISIZEY2,KJSIZEY2    )
!     ################################################################
!
!!****  *GET_SIZEY_LB* - routine to get the size of the South-North LB arrays
!!                       for each sub-domain
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to define the size of the  South-North
!!    LB arrays for each sub-domain. All the possible LB configurations are
!     considered: full LB for mass and v points, 4 points LB and 2 points LB  
!
!!**  METHOD
!!    ------
!!      Sequentially, the following sizes are initialized :
!!     - KISIZEYF,KJSIZEYF for full LB for mass points 
!!     and   KISIZEYFV,KJSIZEYFV for full LB for v points 
!!     - KISIZEY4,KJSIZEY4 for 4 points  LB 
!!     -  KISIZEYF2,KJSIZEYF2 for 2 points  LB 
!!     For one sub-domain, the   size of the northern side is added
!!     to the size of  the southern  side; as the sub-domains are  rectangles, 
!!    the sizes of arrays along West/East are the same. 
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
INTEGER, INTENT(IN)  :: KRIMY ! Size of the riming zone 
INTEGER, INTENT(OUT) :: KISIZEYF,KJSIZEYF ! size of the full LB ZONE (y direction)
INTEGER, INTENT(OUT) :: KISIZEYFV,KJSIZEYFV ! size of the full LB ZONE 
                                            ! for u-grid points (y direction)
INTEGER, INTENT(OUT) ::KISIZEY4,KJSIZEY4 ! size of the 4 grid points LB zone
INTEGER, INTENT(OUT) ::KISIZEY2,KJSIZEY2 ! size of the 2 grid points LB zone
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
KISIZEYF=0
KJSIZEYF=0
KISIZEYFV=0
KJSIZEYFV=0
KISIZEY4=0
KJSIZEY4=0
KISIZEY2=0
KJSIZEY2=0
!------------------------------------------------------------------------------
!
!*       2.    COMPUTE THE SIZES
!              -------------------------- 
!
!*       2.1 Full South_north LB zone 
!
IF (KRIMY /=0) THEN 
!  Southern  side, mass points   : 1:NIMAX_ll+ 2 * JPHEXT, 1:NRIMY+1,
  IXOR=1
  IXEND=KIMAX_ll+ 2 * JPHEXT
  IYOR=1 
  IYEND=KRIMY+JPHEXT ! +1
  CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,"EXTE",IINFO)
  IF (IINFO/=1 .AND. LSOUTH_ll() ) THEN  ! no empty intersection
    KISIZEYF= IXENDI - IXORI +1
    KJSIZEYF= KJSIZEYF + (IYENDI - IYORI +1)
  ENDIF
!   Northern  side , mass  and v-grid points:
!            1:NIMAX_ll+2 *JPHEXT,NJMAX_ll+JPHEXT-NRIMY +1:NJMAX_ll+ 2 *JPHEXT,
  IXOR=1
  IXEND=KIMAX_ll+ 2 * JPHEXT
  IYOR=KJMAX_ll + 2 * JPHEXT-KRIMY-JPHEXT+1  ! -KRIMY
  IYEND=KJMAX_ll+ 2 * JPHEXT
  CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,"EXTE",IINFO)
  IF (IINFO/=1 .AND. LNORTH_ll() ) THEN ! no empty intersection
    KISIZEYF=IXENDI - IXORI +1 
    KJSIZEYF= KJSIZEYF + (IYENDI - IYORI +1 )! added to the southern  side 
    KISIZEYFV= IXENDI - IXORI +1
    KJSIZEYFV= KJSIZEYFV + (IYENDI - IYORI +1 )
  ENDIF
! Southern  side, v-grid points : 1:NIMAX_ll+ 2 * JPHEXT,2:NRIMY+2
  IXOR=1 
  IXEND=KIMAX_ll+ 2 * JPHEXT
  IYOR=2
  IYEND=KRIMY+JPHEXT+1   !+2
  CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,"EXTE",IINFO)
  IF (IINFO/=1 .AND. LSOUTH_ll() ) THEN ! no empty intersection
    KISIZEYFV= IXENDI - IXORI +1
    KJSIZEYFV= KJSIZEYFV + (IYENDI - IYORI +1 )
  ENDIF
ENDIF
!
!*       2.2  South-north  LB zone with only  2 points at each side 
!
!   Southern  side : 1:NIMAX_ll+ 2 * JPHEXT,2:3
IXOR=1
IXEND=KIMAX_ll+ 2 * JPHEXT
IYOR=2  !2
IYEND=JPHEXT+2 !3
CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,"EXTE",IINFO)
IF (IINFO /=1  .AND. LSOUTH_ll() ) THEN ! no empty intersection
  KISIZEY4= IXENDI - IXORI +1
  KJSIZEY4= KJSIZEY4 + (IYENDI - IYORI +1)
ENDIF
!   Northern side,  : 1:NIMAX_ll+2 *JPHEXT,NJMAX_ll+JPHEXT:NJMAX_ll+ 2 *JPHEXT
IXOR= 1
IXEND=KIMAX_ll+ 2 * JPHEXT
IYOR=KJMAX_ll + 2 * JPHEXT - JPHEXT      ! + JPHEXT
IYEND=KJMAX_ll+ 2 * JPHEXT - JPHEXT + JPHEXT  ! + 2*JPHEXT
CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,"EXTE",IINFO)
IF (IINFO/=1 .AND. LNORTH_ll() ) THEN ! no empty intersection
  KISIZEY4=IXENDI - IXORI +1
  KJSIZEY4=KJSIZEY4 + ( IYENDI - IYORI +1 )
ENDIF
!
!*       2.3  South-north LB zone with only  1 point at each side 
!
!   Southern  side : 1:NIMAX_ll+ 2 * JPHEXT,1:1
IXOR=1
IXEND=KIMAX_ll+ 2 * JPHEXT
IYOR=1   ! 1
IYEND=JPHEXT  ! 1
CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,"EXTE",IINFO)
IF (IINFO /=1 .AND. LSOUTH_ll() ) THEN   ! no empty intersection
  KISIZEY2= IXENDI - IXORI +1
  KJSIZEY2=KJSIZEY2 + (IYENDI - IYORI +1 )
ENDIF
!   Northern boundary, 1 point LB zone :1:NIMAX_ll+2 *JPHEXT, NJMAX_ll+ 2 * JPHEXT:NJMAX_ll+ 2 *JPHEXT,
IXOR=1
IXEND=KIMAX_ll + 2 * JPHEXT
IYOR=KJMAX_ll + 2 * JPHEXT - JPHEXT + 1 !  + 2 * JPHEXT
IYEND=KJMAX_ll+ 2 * JPHEXT - JPHEXT + JPHEXT !  + 2 * JPHEXT
CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,"EXTE",IINFO)
IF (IINFO /=1 .AND. LNORTH_ll() ) THEN    ! no empty intersection
  KISIZEY2=  IXENDI - IXORI +1
  KJSIZEY2= KJSIZEY2 + (IYENDI - IYORI +1 )
ENDIF
!
END SUBROUTINE GET_SIZEY_LB
   
