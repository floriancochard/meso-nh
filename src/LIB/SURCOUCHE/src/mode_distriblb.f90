!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------
!Correction :
!  J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-----------------------------------------------------------------

!     #############################
      MODULE MODE_DISTRIB_LB
!     #############################
!
IMPLICIT NONE 

PRIVATE 

PUBLIC GET_DISTRIB_LB

CONTAINS 

SUBROUTINE GET_DISTRIB_LB(HLBTYPE,KIP,HCOORD,HMODE,KRIM,KIB,KIE,KJB,KJE)

CHARACTER(LEN=*),INTENT(IN) :: HLBTYPE ! LB type : 'LBX','LBXU','LBY','LBYV'
INTEGER,         INTENT(IN) :: KIP     ! Processor number
CHARACTER(LEN=*),INTENT(IN) :: HCOORD ! 'LOC' : local indices
                                      ! 'FM'  : FM file indices
CHARACTER(LEN=*),INTENT(IN) :: HMODE  ! 'READ' or 'WRITE' mode
INTEGER, INTENT(IN)  :: KRIM ! Size of the riming zone 
! western or south side 
INTEGER, INTENT(OUT) :: KIB,KIE,KJB,KJE
                                    ! coordinates of the LB array on each sub domain


SELECT CASE(HLBTYPE)
CASE('LBX')
  CALL GET_DISTRIBX_LB(KIP,HCOORD,HMODE,KRIM,.FALSE.,KIB,KIE,KJB,KJE)
CASE('LBXU')
  CALL GET_DISTRIBX_LB(KIP,HCOORD,HMODE,KRIM,.TRUE.,KIB,KIE,KJB,KJE)
CASE('LBY')
  CALL GET_DISTRIBY_LB(KIP,HCOORD,HMODE,KRIM,.FALSE.,KIB,KIE,KJB,KJE)
CASE('LBYV')
  CALL GET_DISTRIBY_LB(KIP,HCOORD,HMODE,KRIM,.TRUE.,KIB,KIE,KJB,KJE) 

CASE default
  WRITE(*,*) 'Error in GET_DISTRIB_LB with unknown LB type: ',HLBTYPE
END SELECT

END SUBROUTINE GET_DISTRIB_LB


!     ###################################################################
      SUBROUTINE GET_DISTRIBX_LB(KIP,HCOORD,HMODE,KRIMX,OLU,KIB,KIE,KJB,KJE)
!     ###################################################################
!
!!****  *DISTRIBX_LB* - routine to get the indices of the West-East LB arrays
!!                       for each sub-domain
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to define the indices of the West-East
!     LB arrays for each sub-domain and the corresponding indices
!      for the LB global  arrays in the FM-file.
!    
!
!!**  METHOD
!!    ------
!!     At first, the  LB areas are defined in the global extended domain.
!!    Then the intersection of these areas with each sub-domain is retrieved
!!    by returning the  local indices of the LB areas for each sub-domain.
!!    From these indices, it is determined :
!!     - the indices of the LB arrays for their western and eastern sides
!!     - the corresponding indices in  the global West-east  array in FM-file
!!  
!!    EXTERNAL
!!    --------   
!!      Module MODE_TOOLS_ll : 
!!             GET_INTERSECTION_ll : to get the global indices of the LB 
!!             GET_GLOBAL_DIMS
!!     arrays for KIP sub-domain.
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!    Module MODD_PARAMETERS_ll :  JPHEXT
!!
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
!!      Modif 
!!     J.Escobar 28/03/2019: for very small domain , force N/S/E/W check on getting LB bounds
!-------------------------------------------------------------------------------
!
USE MODD_PARAMETERS_ll,ONLY : JPHEXT
USE MODD_VAR_ll,       ONLY : TCRRT_PROCONF
USE MODD_STRUCTURE_ll,  ONLY : MODELSPLITTING_ll
USE MODE_TOOLS_ll,     ONLY : GET_INTERSECTION_ll,GET_GLOBALDIMS_ll,LWEST_ll,LEAST_ll

!*       0.    DECLARATIONS
!              ------------ 
!*       0.1   declarations of arguments
!
CHARACTER (LEN=*),     INTENT(IN)  :: HCOORD    ! 'LOC' : local indices
                                                ! 'FM'  : FM file indices
CHARACTER (LEN=*),     INTENT(IN)  :: HMODE     ! 'READ' or 'WRITE' mode
INTEGER,               INTENT(IN)   :: KIP       ! Processor number 
INTEGER, INTENT(IN)  :: KRIMX ! Size of the riming zone 
LOGICAL, INTENT(IN)  :: OLU ! switch for u-grid-points
INTEGER, INTENT(OUT) :: KIB,KIE,KJB,KJE 
                                    ! coordinates of the LB array on each sub domain

!
!*       0.2   declarations of local variables
!
INTEGER :: IINFO ! return code of // routines
INTEGER :: IIMAX_ll  !  Dimensions  in x direction of the physical domain,
INTEGER :: IJMAX_ll  !  Dimensions  in y direction of the physical domain,
INTEGER :: IXOR, IYOR, IXEND, IYEND ! Coordinates of the LB zone
                                    ! (global indices)
INTEGER :: IXORI, IYORI, IXENDI, IYENDI ! Coordinates of
                                        ! the intersection (local indices)
INTEGER :: IXOR3DX,IYOR3DX ! orignin's coordinate of the extended subdomain
CHARACTER(LEN=4) :: YMODE ! 'EXTE' or 'PHYS' mode for GET_INTERSECTION_ll routine
TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT ! 2way-splitting
INTEGER :: IL3DX

!-------------------------------------------------------------------------------
!
!*       1.   INITIALIZATIONS
!              -------------- 
!
!   
SELECT CASE(HMODE)
CASE('READ')
  YMODE = 'EXTE'
CASE('WRITE')
  YMODE = 'PHYS'
CASE default
  WRITE(*,*) 'Error in GET_DISTRIBX_LB...'
  STOP
END SELECT
!
CALL  GET_GLOBALDIMS_ll(IIMAX_ll, IJMAX_ll)
!
KIB = 0
KIE = 0
KJB = 0
KJE = 0

!------------------------------------------------------------------------------
!
!*       2.    GET THE INDICES OF THE LB ARRAYS IN MODEL AND IN FM-FILE
!              -------------------------------------------------------- 
!
!
!*       2.1   Western side
!
IF (KRIMX /=0) THEN 
  IF (OLU) THEN
! full LB zone for u grid point  : 2:NRIMX+2,1:NJMAX_ll+ 2 * JPHEXT
    IXOR=2
    IXEND=KRIMX+JPHEXT+1 ! +2
    IYOR=1 
    IYEND=IJMAX_ll+ 2 * JPHEXT
  ELSE
! full LB zone, mass point  : 1:NRIMX+1,1:NJMAX_ll+ 2 * JPHEXT
    IXOR=1
    IXEND=KRIMX+JPHEXT  ! +1
    IYOR=1 
    IYEND=IJMAX_ll+ 2 * JPHEXT
  ENDIF
ELSE
! 1 point LB zone : 1:1,1:NJMAX_ll+ 2 * JPHEXT
  IXOR=1
  IXEND=JPHEXT ! 1
  IYOR=1
  IYEND=IJMAX_ll+ 2 * JPHEXT
ENDIF
CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,YMODE,IINFO,KIP)
IF (IINFO/=1 .AND. LWEST_ll(KIP) ) THEN  ! no empty intersection
  IF (HCOORD == 'LOC') THEN
    KIB=IXORI
    KIE=IXENDI
    IF (OLU .AND. LWEST_ll(KIP)) THEN
      KIB = KIB-1
      KIE = KIE-1
    END IF
    KJB=IYORI
    KJE=IYENDI

  ELSE IF (HCOORD == 'FM') THEN

    TZSPLIT => TCRRT_PROCONF%TSPLITS_B(KIP)
    IXOR3DX = TZSPLIT%NXORE
    IYOR3DX = TZSPLIT%NYORE

    IF (OLU) THEN
      KIB=IXORI + IXOR3DX -2 ! global coordinates (LB arrays begin 
                                   ! at 1 in FM file)
      KIE=IXENDI+ IXOR3DX -2
    ELSE
      KIB=IXORI + IXOR3DX  -1 ! global coordinates
      KIE=IXENDI + IXOR3DX -1
    ENDIF
    KJB=IYORI + IYOR3DX -1
    KJE=IYENDI+ IYOR3DX- 1
  ELSE
    WRITE(*,*) 'Error in GET_DISTRIBX_LB...'
    STOP
  ENDIF
END IF

!
!*       2.1   Eastern side
!                          
!     
IF (KRIMX /=0) THEN 
!  full LB zone : NIMAX_ll+JPHEXT-NRIMX +1:NIMAX_ll+ 2 *JPHEXT,1:NJMAX_ll+2 *JPHEXT
  IXOR =IIMAX_ll+ 2 * JPHEXT - KRIMX-JPHEXT +1 ! -KRIMX
  IXEND=IIMAX_ll+ 2 * JPHEXT
  IYOR=1
  IYEND=IJMAX_ll+ 2 * JPHEXT
ELSE
!    1 point LB zone : NIMAX_ll+ 2 * JPHEXT:NIMAX_ll+ 2 *JPHEXT,1:NJMAX_ll+2 *JPHEXT
  IXOR=IIMAX_ll + 2 * JPHEXT  -JPHEXT +1      ! + 2 * JPHEXT
  IXEND=IIMAX_ll + 2 * JPHEXT -JPHEXT +JPHEXT ! + 2 * JPHEXT
  IYOR=1
  IYEND=IJMAX_ll+ 2 * JPHEXT
ENDIF
CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,YMODE,IINFO,KIP)
IF (IINFO/=1  .AND. LEAST_ll(KIP) ) THEN
  IF (HCOORD == 'LOC') THEN
    IF (KIB == 0) KIB=1+KIE
    KIE=KIE+1+IXENDI-IXORI
    KJB=IYORI
    KJE=IYENDI
  ELSE IF (HCOORD == 'FM') THEN
    TZSPLIT => TCRRT_PROCONF%TSPLITS_B(KIP)
    IXOR3DX = TZSPLIT%NXORE
    IYOR3DX = TZSPLIT%NYORE
!
    IL3DX = 2*(KRIMX+JPHEXT) ! +1
    
    IF (KIB == 0) KIB = IL3DX - (IIMAX_ll+2*JPHEXT - IXORI - IXOR3DX  +1)
    KIE=IL3DX - (IIMAX_ll + 2 *JPHEXT - IXENDI- IXOR3DX  +1)
    KJB=IYORI  + IYOR3DX -1
    KJE=IYENDI + IYOR3DX -1
  ENDIF
END IF
!
END SUBROUTINE GET_DISTRIBX_LB


!
!     ###################################################################
      SUBROUTINE GET_DISTRIBY_LB(KIP,HCOORD,HMODE,KRIMY,OLV,KIB,KIE,KJB,KJE)
!     ###################################################################
!
!!****  *DISTRIBY_LB* - routine to get the indices of the South-North LB arrays
!!                       for each sub-domain
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to define the indices of the South-North
!     LB arrays for each sub-domain and the corresponding indices
!      for the LB global  arrays in the FM-file.
!    
!
!!**  METHOD
!!    ------
!!     At first, the  LB areas are defined in the global extended domain.
!!    Then the intersection of these areas with each sub-domain is retrieved
!!    by returning the  local indices of the LB areas for each sub-domain.
!!    From these indices, we  determine :
!!     - the indices of the LB arrays for their southern  and northern sides
!!     - the corresponding indices in  the global south-north  array in FM-file
!!  
!!    EXTERNAL
!!    --------   
!!      Module MODE_TOOLS_ll : 
!!             GET_INTERSECTION_ll : to get the global indices of the LB 
!!     arrays for each sub-domain.
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!    Module MODD_PARAMETERS_ll :  JPHEXT
!!
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
!-------------------------------------------------------------------------------
!
USE MODD_PARAMETERS_ll,ONLY : JPHEXT
USE MODD_VAR_ll,       ONLY : TCRRT_PROCONF
USE MODD_STRUCTURE_ll,  ONLY : MODELSPLITTING_ll
USE MODE_TOOLS_ll,     ONLY : GET_INTERSECTION_ll,GET_GLOBALDIMS_ll,LSOUTH_ll,LNORTH_ll
!*       0.    DECLARATIONS
!              ------------ 
!*       0.1   declarations of arguments
!
CHARACTER (LEN=*),     INTENT(IN)  :: HCOORD    ! 'LOC' : local indices
                                                ! 'FM'  : FM file indices
CHARACTER (LEN=*),     INTENT(IN)  :: HMODE     ! 'READ' or 'WRITE' mode
INTEGER,               INTENT(IN)   :: KIP       ! Processor number 
INTEGER, INTENT(IN)  :: KRIMY ! Size of the riming zone 

LOGICAL, INTENT(IN)  :: OLV ! switch for v-grid-points
INTEGER, INTENT(OUT) :: KIB,KIE,KJB,KJE
!
!*       0.2   declarations of local variables
!
INTEGER :: IINFO ! return code of // routines
INTEGER :: IIMAX_ll  !  Dimensions  in x direction of the physical domain,
INTEGER :: IJMAX_ll  !  Dimensions  in y direction of the physical domain,
INTEGER :: IXOR, IYOR, IXEND, IYEND ! Coordinates of the LB zone
                                    ! (global indices)
INTEGER :: IXORI, IYORI, IXENDI, IYENDI ! Coordinates of
                                        ! the intersection (local indices)
INTEGER :: IXOR3DY,IYOR3DY ! orignin's coordinate of the extended subdomain
CHARACTER(LEN=4) :: YMODE ! 'EXTE' or 'PHYS' mode for GET_INTERSECTION_ll routine
TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT ! 2way-splitting
INTEGER      :: IL3DY ! Size of the LB array in FM-file without gap 
!-------------------------------------------------------------------------------
!
!*       1.   INITIALIZATIONS
!              -------------- 
!
!
SELECT CASE(HMODE)
CASE('READ')
  YMODE = 'EXTE'
CASE('WRITE')
  YMODE = 'PHYS'
CASE default
  WRITE(*,*) 'Error in GET_DISTRIBX_LB...'
  STOP
END SELECT
!   
CALL  GET_GLOBALDIMS_ll(IIMAX_ll, IJMAX_ll)
!
KIB=0
KIE=0
KJB=0
KJE=0
!------------------------------------------------------------------------------
!
!*       2.    GET THE INDICES OF THE LB ARRAYS IN MODEL AND IN FM-FILE
!              -------------------------------------------------------- 
!
!
!*       2.1   Southern side
!
IF (KRIMY /=0) THEN 
  IF (OLV) THEN
! full LB zone for v grid point  : 1:NIMAX_ll+ 2 * JPHEXT, 2:NRIMY+2
    IXOR=1
    IXEND=IIMAX_ll+ 2 * JPHEXT
    IYOR=2 
    IYEND=KRIMY+JPHEXT+1   !+2
  ELSE
! full LB zone, mass point  : 1:NIMAX_ll+ 2 * JPHEXT,1:NRIMY+1
    IXOR=1
    IXEND=IIMAX_ll+ 2 * JPHEXT
    IYOR=1 
    IYEND=KRIMY+JPHEXT  !+1 
  ENDIF
ELSE
! 1 point LB zone : 1:NIMAX_ll+ 2 * JPHEXT,1:1
  IXOR=1
  IXEND=IIMAX_ll+ 2 * JPHEXT
  IYOR=1
  IYEND=JPHEXT !1
ENDIF
CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,YMODE,IINFO,KIP)
IF (IINFO/=1 .AND. LSOUTH_ll(KIP) ) THEN  ! no empty intersection
  IF (HCOORD == 'LOC') THEN
    KIB=IXORI
    KIE=IXENDI
    KJB=IYORI
    KJE=IYENDI
    IF (OLV .AND. LSOUTH_ll(KIP)) THEN
      KJB=KJB-1
      KJE=KJE-1
    END IF
  ELSE IF (HCOORD == 'FM') THEN

    TZSPLIT => TCRRT_PROCONF%TSPLITS_B(KIP)
    IXOR3DY = TZSPLIT%NXORE
    IYOR3DY = TZSPLIT%NYORE

    IF (OLV) THEN
      KJB=IYORI + IYOR3DY-2 ! global coordinates (LB arrays begin 
      KJE=IYENDI+ IYOR3DY-2 ! at 1 in FM file)
    ELSE
      KJB=IYORI + IYOR3DY-1 ! global coordinates
      KJE=IYENDI+ IYOR3DY-1
    ENDIF
    KIB=IXORI   + IXOR3DY-1 ! global coordinates
    KIE=IXENDI  + IXOR3DY-1
  END IF
ENDIF
 !
!*       2.1   Northern side
!                           

!    
IF (KRIMY /=0) THEN 
!  full LB zone :1:NIMAX_ll+2 *JPHEXT, NJMAX_ll+JPHEXT-NRIMY +1:NJMAX_ll+ 2 *JPHEXT,
  IXOR =1
  IXEND=IIMAX_ll+ 2 * JPHEXT
  IYOR = IJMAX_ll+ 2 * JPHEXT -KRIMY-JPHEXT +1 ! - KRIMY
  IYEND=IJMAX_ll+ 2 * JPHEXT
ELSE
!    1 point LB zone : 1:NJMAX_ll+2 *JPHEXT,NJMAX_ll+ 2 * JPHEXT:NJMAX_ll+ 2 *JPHEXT
  IXOR=1
  IXEND=IIMAX_ll+ 2 * JPHEXT
  IYOR=IJMAX_ll + 2 * JPHEXT  - JPHEXT + 1       ! + 2 * JPHEXT 
  IYEND=IJMAX_ll + 2 * JPHEXT - JPHEXT + JPHEXT  ! + 2 * JPHEXT 
ENDIF
CALL GET_INTERSECTION_ll(IXOR,IYOR,IXEND,IYEND,IXORI,IYORI,IXENDI,IYENDI,YMODE,IINFO,KIP)
IF (IINFO/=1 .AND. LNORTH_ll(KIP) ) THEN
  IF (HCOORD == 'LOC') THEN
    KIB=IXORI
    KIE=IXENDI
    IF (KJB == 0) KJB=1 + KJE
    KJE=1+KJE + IYENDI-IYORI
  ELSE IF (HCOORD == 'FM') THEN
    TZSPLIT => TCRRT_PROCONF%TSPLITS_B(KIP)
    IXOR3DY = TZSPLIT%NXORE
    IYOR3DY = TZSPLIT%NYORE
!
    IL3DY = 2*(KRIMY+JPHEXT ) ! +1
    KIB=IXORI  + IXOR3DY -1
    KIE=IXENDI + IXOR3DY -1
    IF (KJB == 0) KJB = IL3DY - (IJMAX_ll + 2 *JPHEXT - IYORI - IYOR3DY  +1)
    KJE = IL3DY - (IJMAX_ll + 2 *JPHEXT - IYENDI- IYOR3DY  +1)
  END IF
ENDIF
!
END SUBROUTINE GET_DISTRIBY_LB
   
END MODULE MODE_DISTRIB_LB
