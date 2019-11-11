!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###############################
      MODULE MODI_ADVEC_4TH_ORDER_AUX
!     ###############################
!
INTERFACE
!
      SUBROUTINE ADVEC_4TH_ORDER_ALGO(HLBCX, HLBCY, PFIELDT, KGRID, & 
                                      PMEANX, PMEANY,TPHALO2 )
!
USE MODE_ll
!
USE MODD_CONF
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PMEANX, PMEANY ! fluxes
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PFIELDT  ! variable at t
INTEGER,                INTENT(IN)    :: KGRID    ! C grid localisation
!
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
END SUBROUTINE ADVEC_4TH_ORDER_ALGO
!
!-------------------------------------------------------------------------------
!
      FUNCTION MZF4(PA)  RESULT(PMZF4)
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at flux
                                                            !         side
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PMZF4  ! result at mass
                                                            ! localization 
!
END FUNCTION MZF4
!
!-------------------------------------------------------------------------------
!
      FUNCTION MZM4(PA)  RESULT(PMZM4)
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at mass 
                                                            ! localization
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PMZM4  ! result at flux 
                                                            ! localization 
END FUNCTION MZM4
!
!-------------------------------------------------------------------------------
!
END INTERFACE
!
END MODULE MODI_ADVEC_4TH_ORDER_AUX
!
!-------------------------------------------------------------------------------
!
!     ########################################################################
      SUBROUTINE ADVEC_4TH_ORDER_ALGO(HLBCX, HLBCY, PFIELDT, KGRID, & 
                                      PMEANX, PMEANY,TPHALO2 )
!     ########################################################################
!!
!!****  *ADVEC_4TH_ORDER_ALGO * - routine used to compute 4th order horizontal
!!                                advection fluxes of 3D prognostic variables
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute 2sd or 4th order horizontal
!!    advection fluxes of a prognostic variable.
!!
!!**  METHOD
!!    ------
!!      In case of cyclic LBCs, the routine returns the scalar component of the 
!!    advection fluxes by applying a 4th order horizontal averaging operator to
!!    the prognostic variable on each grid level. In the case of open LBCs, the
!!    averaging operator degenerates to a 2nd order one on the first ring 
!!    inside the computationnal domain.
!!      The "halo2" (or the second layer of the halo) of the prognostic
!!    variable is passed as argument.
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODULE MODD_ARGSLIST
!!         HALO2LIST_ll : type for a list of "HALO2_lls"
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine ADVEC_4TH_ORDER )
!!      User Interface for the MesoNH Parallel Package
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   25/10/05
!!      Modif
!!      J.Escobar 21/03/2013: for HALOK comment all NHALO=1 test
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll
!
USE MODD_CONF
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
USE MODE_IO_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PMEANX, PMEANY ! fluxes
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PFIELDT  ! variable at t
INTEGER,                INTENT(IN)    :: KGRID    ! C grid localisation
!
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
!*       0.2   Declarations of local variables :
!
INTEGER:: IW,IE,IS,IN,IT,IB,IWF,IEF,ISF,INF   ! Coordinate of forth order diffusion area
!
INTEGER:: IIB,IJB        ! Begining useful area in x,y directions
INTEGER:: IIE,IJE        ! End useful area in x,y directions
!
INTEGER:: ILUOUT,IRESP   ! for prints
!
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
!-------------------------------------------------------------------------------
!
!*       0.4.   INITIALIZE THE FIELDS 
!               ---------------------
!
PMEANX(:,:,:) = 0.0
PMEANY(:,:,:) = 0.0
!
!-------------------------------------------------------------------------------
!
!
!*       1.     CALCULATE THE NUMERICAL MEAN IN THE X DIRECTION
!               -----------------------------------------------
!
SELECT CASE ( HLBCX(1) ) ! X direction LBC type: (1) for left side
!
!*       1.1    CYCLIC CASE IN THE X DIRECTION:
!
CASE ('CYCL')          ! In that case one must have HLBCX(1) == HLBCX(2)
!
    IW=IIB+1
    IE=IIE
!
  IF(KGRID == 2) THEN
    IWF=IW-1
    IEF=IE-1
  ELSE
    IWF=IW
    IEF=IE
  END IF
!
!* lateral boundary conditions
  PMEANX(IWF-1,:,:) = (7.0*( PFIELDT(IW-1,:,:)+PFIELDT(IW-2,:,:) ) -  &
                           ( PFIELDT(IW,:,:)+TPHALO2%WEST(:,:) ) )/12.0
!
  PMEANX(IEF+1,:,:) = (7.0*( PFIELDT(IE+1,:,:)+PFIELDT(IE,:,:) ) -  &
                       ( TPHALO2%EAST(:,:)+PFIELDT(IE-1,:,:) ) )/12.0
!
!* inner domain
  PMEANX(IWF:IEF,:,:) = (7.0*( PFIELDT(IW:IE,:,:)+PFIELDT(IW-1:IE-1,:,:) ) -  &
                       ( PFIELDT(IW+1:IE+1,:,:)+PFIELDT(IW-2:IE-2,:,:) ) )/12.0
!
!!$!
!!$
!!$  IF(NHALO == 1) THEN
!!$    PMEANX(IWF-1,:,:) = (7.0*( PFIELDT(IW-1,:,:)+PFIELDT(IW-2,:,:) ) -  &
!!$                             ( PFIELDT(IW,:,:)+TPHALO2%WEST(:,:) ) )/12.0
!!$!
!!$    PMEANX(IEF+1,:,:) = (7.0*( PFIELDT(IE+1,:,:)+PFIELDT(IE,:,:) ) -  &
!!$                         ( TPHALO2%EAST(:,:)+PFIELDT(IE-1,:,:) ) )/12.0
!!$  ENDIF
!!$!
!!$  PMEANX(IWF:IEF,:,:) = (7.0*( PFIELDT(IW:IE,:,:)+PFIELDT(IW-1:IE-1,:,:) ) -  &
!!$                       ( PFIELDT(IW+1:IE+1,:,:)+PFIELDT(IW-2:IE-2,:,:) ) )/12.0
!!$!
!*       1.2    NON CYCLIC CASE IN THE X DIRECTION 
!
CASE ('OPEN','WALL','NEST')
!
  IF (LWEST_ll()) THEN
    IF(KGRID == 2) THEN
      IW=IIB+2          ! special case of C grid
    ELSE
      IW=IIB+1
    END IF
  ELSE
!!$    IF(NHALO == 1) THEN
      IW=IIB+1
!!$    ELSE
!!$      IW=IIB
!!$    ENDIF
  ENDIF
!!$  IF (LEAST_ll() .OR. NHALO == 1) THEN
  IF (LEAST_ll() ) THEN
! T. Maric
!    IE=IIE-1 ! original
    IE=IIE
  ELSE
    IE=IIE
  END IF  
!
  IF(KGRID == 2) THEN
    IWF=IW-1
    IEF=IE-1
  ELSE
    IWF=IW
    IEF=IE
  END IF
!
! T. Maric. 16.1.2006.
!  write(*,*)' IW, IE, IWF, IEF = ',IW, IE, IWF, IEF
!  stop 'Stopping in advec_4th_order_aux.f90'
!
!*             Use a second order scheme at the physical border
!
  IF (LWEST_ll()) THEN
    PMEANX(IWF-1,:,:) = 0.5*( PFIELDT(IW-1,:,:)+PFIELDT(IW-2,:,:) )
    ! T. Maric
    ! PMEANX(1,:,:) = PMEANX(IWF-1,:,:)
    ! extrapolate
    !PMEANX(1,:,:) = 0.5*(3.0*PFIELDT(1,:,:) - PFIELDT(2,:,:))
!!$  ELSE IF (NHALO == 1) THEN
  ELSE
    PMEANX(IWF-1,:,:) = (7.0*( PFIELDT(IW-1,:,:)+PFIELDT(IW-2,:,:) ) -  &
                             ( PFIELDT(IW,:,:)+TPHALO2%WEST(:,:) ) )/12.0
  ENDIF
!
  IF (LEAST_ll()) THEN
    PMEANX(IEF+1,:,:) = 0.5*( PFIELDT(IE+1,:,:)+PFIELDT(IE,:,:) )
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    PMEANX(IEF+1,:,:) = (7.0*( PFIELDT(IE+1,:,:)+PFIELDT(IE,:,:) ) -  &
                         ( TPHALO2%EAST(:,:)+PFIELDT(IE-1,:,:) ) )/12.0
  ENDIF
!
!*             Use a fourth order scheme elsewhere 
!
  PMEANX(IWF:IEF,:,:) = (7.0*( PFIELDT(IW:IE,:,:)+PFIELDT(IW-1:IE-1,:,:) ) -  &
                       ( PFIELDT(IW+1:IE+1,:,:)+PFIELDT(IW-2:IE-2,:,:) ) )/12.0
END SELECT
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTES THE 4TH ORDER MEAN IN THE Y DIRECTION
!               ----------------------------------------------
!
IF ( .NOT. L2D ) THEN
  SELECT CASE ( HLBCY(1) ) ! Y direction LBC type: (1) for left side
!
!*       2.1    CYCLIC CASE IN THE Y DIRECTION:
!
  CASE ('CYCL')          ! In that case one must have HLBCY(1) == HLBCY(2)
!
!
      IS=IJB+1
      IN=IJE
!
    IF(KGRID == 3) THEN
      ISF=IS-1
      INF=IN-1
    ELSE
      ISF=IS
      INF=IN
    END IF
!
!* lateral boundary conditions
    PMEANY(:,ISF-1,:) = (7.0*( PFIELDT(:,IS-1,:)+PFIELDT(:,IS-2,:) ) -  &
                        ( PFIELDT(:,IS,:)+TPHALO2%SOUTH(:,:) ) )/12.0
!
    PMEANY(:,INF+1,:) = (7.0*( PFIELDT(:,IN+1,:)+PFIELDT(:,IN,:) ) -  &
                        ( TPHALO2%NORTH(:,:)+PFIELDT(:,IN-1,:) ) )/12.0
!
!* inner domain
    PMEANY(:,ISF:INF,:) = (7.0*( PFIELDT(:,IS:IN,:)+PFIELDT(:,IS-1:IN-1,:)) -  &
                         ( PFIELDT(:,IS+1:IN+1,:)+PFIELDT(:,IS-2:IN-2,:) ))/12.0
!!$!
!!$    IF(NHALO == 1) THEN
!!$      PMEANY(:,ISF-1,:) = (7.0*( PFIELDT(:,IS,:)+PFIELDT(:,IS-1,:) ) -  &
!!$                          ( PFIELDT(:,IS+1,:)+TPHALO2%SOUTH(:,:) ) )/12.0
!!$!
!!$      PMEANY(:,ISF+1,:) = (7.0*( PFIELDT(:,IS,:)+PFIELDT(:,IS-1,:) ) -  &
!!$                          ( TPHALO2%NORTH(:,:)+PFIELDT(:,IS-2,:) ) )/12.0
!!$    ENDIF
!!$!
!!$    PMEANY(:,ISF:INF,:) = (7.0*( PFIELDT(:,IS:IN,:)+PFIELDT(:,IS-1:IN-1,:)) -  &
!!$                         ( PFIELDT(:,IS+1:IN+1,:)+PFIELDT(:,IS-2:IN-2,:) ))/12.0
!!$!
!*       2.2    NON CYCLIC CASE IN THE Y DIRECTION
!
  CASE ('OPEN','WALL','NEST')
!
    IF (LSOUTH_ll()) THEN
      IF(KGRID == 3) THEN
        IS=IJB+2          ! special case of C grid
      ELSE
        IS=IJB+1
      END IF
    ELSE
!!$      IF(NHALO == 1) THEN
        IS=IJB+1
!!$      ELSE
!!$        IS=IJB
!!$      ENDIF
    ENDIF
!!$    IF (LNORTH_ll() .OR. NHALO == 1) THEN
    IF (LNORTH_ll()) THEN
! T. Maric
!      IN=IJE-1  ! original
      IN=IJE
    ELSE
      IN=IJE
    END IF  
!
    IF(KGRID == 3) THEN
      ISF=IS-1
      INF=IN-1
    ELSE
      ISF=IS
      INF=IN
    END IF
!
!*             Use a second order scheme at the physical border
!
    IF (LSOUTH_ll()) THEN
      PMEANY(:,ISF-1,:) = 0.5*( PFIELDT(:,IS-1,:)+PFIELDT(:,IS-2,:) )
      ! T. Maric
      ! PMEANY(:,1,:) = PMEANY(:,ISF-1,:)
      ! extrapolate
      !PMEANY(:,1,:) = 0.5*(3.0*PFIELDT(:,1,:) - PFIELDT(:,2,:))
!!$    ELSEIF (NHALO == 1) THEN
    ELSE
!!$      PMEANY(:,ISF-1,:) = (7.0*( PFIELDT(:,IS,:)+PFIELDT(:,IS-1,:)) -  &
!!$                          ( PFIELDT(:,IS+1,:)+TPHALO2%SOUTH(:,:) ))/12.0
       PMEANY(:,ISF-1,:) = (7.0*( PFIELDT(:,IS-1,:)+PFIELDT(:,IS-2,:)) -  &
                           ( PFIELDT(:,IS,:)+TPHALO2%SOUTH(:,:) ))/12.0
    ENDIF
!
    IF (LNORTH_ll()) THEN
      PMEANY(:,INF+1,:) = 0.5*( PFIELDT(:,IN+1,:)+PFIELDT(:,IN,:) )
!!$    ELSEIF (NHALO == 1) THEN
    ELSE
!!$      PMEANY(:,INF+1,:) = (7.0*( PFIELDT(:,IN,:)+PFIELDT(:,IN-1,:)) -  &
!!$                          ( TPHALO2%NORTH(:,:)+PFIELDT(:,IN-2,:) ))/12.0
       PMEANY(:,INF+1,:) = (7.0*( PFIELDT(:,IN+1,:)+PFIELDT(:,IN,:)) -  &
                           ( TPHALO2%NORTH(:,:)+PFIELDT(:,IN-1,:) ))/12.0
    ENDIF
!
!*             Use a fourth order scheme elsewhere 
!
    PMEANY(:,ISF:INF,:) = (7.0*( PFIELDT(:,IS:IN,:)+PFIELDT(:,IS-1:IN-1,:)) -  &
                         ( PFIELDT(:,IS+1:IN+1,:)+PFIELDT(:,IS-2:IN-2,:) ))/12.0
!
  END SELECT
ELSE
  PMEANY(:,:,:) = 0.0
ENDIF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ADVEC_4TH_ORDER_ALGO
!
!-------------------------------------------------------------------------------
!
!     ################################
      FUNCTION MZF4(PA)  RESULT(PMZF4)
!     ################################
!
!!****  *MZF4* - 4th order Shuman operator : mean operator in z direction for a 
!!                                 variable at a flux side
!!
!!    PURPOSE
!!    -------
!!      The purpose of this function is to compute a 4th order mean value
!!    along the z direction (K index) for a field PA localized at a z-flux
!!    point (w point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------ 
!!        The result PMZF4(:,:,k) is defined by 
!!        PMZF4(:,:,k)=0.5*(PA(:,:,k)+PA(:,:,k+1)) at k=1 and size(PA,3)-1
!!        PMZF4(:,:,k)=-999.                       at k=size(PA,3)
!!        PMZF4(:,:,k)=7/12*(PA(:,:,k)+PA(:,:,k+1))
!!                    -1/12*(PA(:,:,k-1)+PA(:,:,k+2)) elsewhere
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty       * Lab Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/10/05
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at flux
                                                            !         side
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PMZF4  ! result at mass
                                                            ! localization 
!
!*       0.2   Declarations of local variables
!
!
INTEGER :: JK               ! loop index in z direction
INTEGER :: IKU              ! upper bound in z direction of PA 
!     
INTEGER :: IIU,IJU,IIJU     ! upper bounds in the x and y directions of PA
INTEGER :: JIJ,JIJK         ! running loop indexes after linearisation
INTEGER :: JIJKOR1,JIJKEND1 ! loop boundaries
INTEGER :: JIJKOR2,JIJKEND2 ! loop boundaries
INTEGER :: JIJKOR3,JIJKEND3 ! loop boundaries
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MZF4
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
IKU = SIZE(PA,3)
!
IIJU = IIU*IJU
!
JIJKOR1  = 1 + IIJU
JIJKEND1 = 2*IIJU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR1 , JIJKEND1
   PMZF4(JIJK-IIJU,1,1) = 0.5*( PA(JIJK-IIJU,1,1)+PA(JIJK,1,1) )
END DO
!
JIJKOR2  = 1 + JIJKEND1
JIJKEND2 = IIJU*IKU - IIJU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR2 , JIJKEND2
   PMZF4(JIJK-IIJU,1,1) = (7.0*( PA(JIJK,1,1)+PA(JIJK-IIJU,1,1) ) -      &
                          ( PA(JIJK+IIJU,1,1)+PA(JIJK-2*IIJU,1,1) ) )/12.0
END DO
!
JIJKOR3  = 1 + JIJKEND2
JIJKEND3 = IIJU*IKU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR3 , JIJKEND3
   PMZF4(JIJK-IIJU,1,1) = 0.5*( PA(JIJK-IIJU,1,1)+PA(JIJK,1,1) )
END DO
!
!CDIR NODEP
!OCL NOVREC
DO JIJ=1,IIJU
   PMZF4(JIJ,1,IKU)    = -999.
END DO
!
!-------------------------------------------------------------------------------
!
END FUNCTION MZF4
!
!-------------------------------------------------------------------------------
!
!     ################################
      FUNCTION MZM4(PA)  RESULT(PMZM4)
!     ################################
!
!!****  *MZM4* - 4th order Shuman operator : mean operator in z direction for a 
!!                                 mass variable 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this function is to compute a 4th order mean value
!!    along the z direction (K index) for a field PA localized at a mass
!!    point. The result is localized at a z-flux point (w point).
!!
!!**  METHOD
!!    ------ 
!!        The result PMZM4(:,:,k) is defined by 
!!        PMZM4(:,:,k)=0.5*(PA(:,:,k)+PA(:,:,k+1)) at k=2 and size(PA,3)
!!        PMZM4(:,:,k)=-999.                       at k=1
!!        PMZM4(:,:,k)=7/12*(PA(:,:,k)+PA(:,:,k+1))
!!                    -1/12*(PA(:,:,k-1)+PA(:,:,k+2)) elsewhere
!!    
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty       * Lab Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/10/05
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at mass 
                                                            ! localization
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PMZM4  ! result at flux 
                                                            ! localization 
!
!*       0.2   Declarations of local variables
!
!
INTEGER :: JK               ! loop index in z direction
INTEGER :: IKU              ! upper bound in z direction of PA 
!     
INTEGER :: IIU,IJU,IIJU     ! upper bounds in the x and y directions of PA
INTEGER :: JIJ,JIJK         ! running loop indexes after linearisation
INTEGER :: JIJKOR1,JIJKEND1 ! loop boundaries
INTEGER :: JIJKOR2,JIJKEND2 ! loop boundaries
!           
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MZM4
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
IKU = SIZE(PA,3)
!
IIJU = IIU*IJU
!
JIJKOR1  = 1 + IIJU
JIJKEND1 = JIJKOR1 + IIJU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR1 , JIJKEND1
   PMZM4(JIJK,1,1) = 0.5*( PA(JIJK,1,1)+PA(JIJK-IIJU,1,1) )
END DO
!
JIJKOR2  = 1 + JIJKEND1
JIJKEND2 = IIJU*IKU - IIJU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR2 , JIJKEND2
   PMZM4(JIJK,1,1) = (7.0*( PA(JIJK,1,1)+PA(JIJK-IIJU,1,1) ) -      &
                          ( PA(JIJK+IIJU,1,1)+PA(JIJK-2*IIJU,1,1) ) )/12.0
END DO
!
!CDIR NODEP
!OCL NOVREC
DO JIJ=1,IIJU
   PMZM4(JIJ,1,IKU) = 0.5*( PA(JIJ,1,IKU)+PA(JIJ-IIJU,1,IKU) )
END DO
!
!CDIR NODEP
!OCL NOVREC
DO JIJ=1,IIJU
   PMZM4(JIJ,1,1)    = -999.
END DO
!
!-------------------------------------------------------------------------------
!
END FUNCTION MZM4
