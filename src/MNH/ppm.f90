!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-----------------------------------------------------------------
!     ###############
      MODULE MODI_PPM
!     ###############
!
INTERFACE
!
FUNCTION PPM_01_X(HLBCX, KGRID, PSRC, PCR, PRHO, PTSTEP) &
               RESULT(PR)
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
END FUNCTION PPM_01_X
!
FUNCTION PPM_01_Y(HLBCY, KGRID, PSRC, PCR, PRHO, PTSTEP) &
               RESULT(PR)
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY  ! Y direction LBC type
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
END FUNCTION PPM_01_Y
!
FUNCTION PPM_01_Z(KGRID, PSRC, PCR, PRHO, PTSTEP) RESULT(PR)
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
END FUNCTION PPM_01_Z
!
FUNCTION PPM_S0_X(HLBCX, KGRID, PSRC, PCR, PRHO, PTSTEP) &
               RESULT(PR)
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
END FUNCTION PPM_S0_X
!
FUNCTION PPM_S0_Y(HLBCY, KGRID, PSRC, PCR, PRHO, PTSTEP) &
               RESULT(PR)
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
END FUNCTION PPM_S0_Y
!
FUNCTION PPM_S0_Z(KGRID, PSRC, PCR, PRHO, PTSTEP) &
               RESULT(PR)
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
END FUNCTION PPM_S0_Z
!
FUNCTION PPM_S1_X(HLBCX, KGRID, PSRC, PCR, PRHO, PRHOT, &
                        PTSTEP) RESULT(PR)
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHOT ! density at t+dt
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
END FUNCTION PPM_S1_X
!
FUNCTION PPM_S1_Y(HLBCY, KGRID, PSRC, PCR, PRHO, PRHOT, &
                        PTSTEP) RESULT(PR)
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY  ! X direction LBC type
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHOT ! density at t+dt
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
END FUNCTION PPM_S1_Y
!
FUNCTION PPM_S1_Z(KGRID, PSRC, PCR, PRHO, PRHOT, PTSTEP) &
                        RESULT(PR)
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHOT ! density at t+dt
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
END FUNCTION PPM_S1_Z
!
END INTERFACE
!
END MODULE MODI_PPM
!
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ########################################################################
      FUNCTION PPM_01_X(HLBCX, KGRID, PSRC, PCR, PRHO, PTSTEP) &
               RESULT(PR)
!     ########################################################################
!!
!!****  PPM_01_X - PPM_01 fully monotonic PPM advection scheme in X direction
!!                 Colella notation
!!
!!    MODIFICATIONS
!!    -------------
!!    
!!    11.5.2006.  T. Maric - original version
!!      J.Escobar 21/03/2013: for HALOK comment all NHALO=1 test
!!      J.Escobar 28/06/2018: limit computation on TAB(:,IJS:IJN,:) to avoid unneeded NaN
!!      J.Escobr  16/07/2018: still NaN pb => reintroduce initialization of temporary local array
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODE_IO_ll
USE MODI_SHUMAN
USE MODI_GET_HALO
!
USE MODD_CONF
!USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER:: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER:: IIE,IJE    ! End useful area in x,y,z directions
!
! terms used in parabolic interpolation, dmq, qL, qR, dq, q6
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZQL,ZQR
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZDQ,ZQ6
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZDMQ
!
! extra variables for the initial guess of parabolae parameters
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZQL0,ZQR0,ZQ60
!
! advection fluxes
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZFPOS, ZFNEG
!
!BEG JUAN PPM_LL
INTEGER                          :: IJS,IJN
!END JUAN PPM_LL
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IJS=IJB
IJN=IJE
!
!BEG JUAN PPM_LL
!
!*              initialise & update halo & halo2 for PSRC
!
CALL GET_HALO(PSRC)
PR=PSRC
ZQL=PSRC
ZQR=PSRC
ZDQ=PSRC
ZQ6=PSRC
ZDMQ=PSRC
ZQL0=PSRC
ZQR0=PSRC
ZQ60=PSRC
ZFPOS=PSRC
ZFNEG=PSRC
!
!-------------------------------------------------------------------------------
!
SELECT CASE ( HLBCX(1) ) ! X direction LBC type: (1) for left side
!
!        1.1   CYCLIC BOUNDARY CONDITIONS IN X DIRECTION
!              -----------------------------------------
!
CASE ('CYCL','WALL')          ! In that case one must have HLBCX(1) == HLBCX(2)
!
! calculate dmq
   ZDMQ = DIF2X(PSRC)
!
! monotonize the difference followinq eq. 5 in Lin94
!
!BEG JUAN PPM_LL01
!
!  ZDMQ(i) = Fct[ ZDMQ(i),PSRC(i),PSRC(i-1),PSRC(i+1) ]
!
   ZDMQ(IIB:IIE,IJS:IJN,:) = &
        SIGN( (MIN( ABS(ZDMQ(IIB:IIE,IJS:IJN,:)),2.0*(PSRC(IIB:IIE,IJS:IJN,:) - &
        MIN(PSRC(IIB-1:IIE-1,IJS:IJN,:),PSRC(IIB:IIE,IJS:IJN,:),PSRC(IIB+1:IIE+1,IJS:IJN,:))),    &
        2.0*(MAX(PSRC(IIB-1:IIE-1,IJS:IJN,:),PSRC(IIB:IIE,IJS:IJN,:),PSRC(IIB+1:IIE+1,IJS:IJN,:)) - &
        PSRC(IIB:IIE,IJS:IJN,:)) )), ZDMQ(IIB:IIE,IJS:IJN,:) )
!
!  WEST BOUND
!
!!$   ZDMQ(IIB-1,:,:) = & 
!!$        SIGN( (MIN( ABS(ZDMQ(IIB-1,:,:)), 2.0*(PSRC(IIB-1,:,:) - &
!!$        MIN(PSRC(IIE-1,:,:),PSRC(IIB-1,:,:),PSRC(IIB,:,:))),   &
!!$        2.0*(MAX(PSRC(IIE-1,:,:),PSRC(IIB-1,:,:),PSRC(IIB,:,:)) - &
!!$        PSRC(IIB-1,:,:)) )), ZDMQ(IIB-1,:,:) )
!
!  EAST BOUND
!
!!$   ZDMQ(IIE+1,:,:) = &
!!$        SIGN( (MIN( ABS(ZDMQ(IIE+1,:,:)), 2.0*(PSRC(IIE+1,:,:) - &
!!$        MIN(PSRC(IIE,:,:),PSRC(IIE+1,:,:),PSRC(IIB+1,:,:))),  &
!!$        2.0*(MAX(PSRC(IIE,:,:),PSRC(IIE+1,:,:),PSRC(IIB+1,:,:)) - &
!!$        PSRC(IIE+1,:,:)) )), ZDMQ(IIE+1,:,:) )
!
!  update ZDMQ HALO before next/further  utilisation 
!
   CALL  GET_HALO(ZDMQ)
!
! calculate qL and qR with the modified dmq
!
!  ZQL0(i) = Fct[ PSRC(i),PSRC(i-1),ZDMQ(i),ZDMQ(i-1) ]
!
   ZQL0(IIB:IIE+1,IJS:IJN,:) = 0.5*(PSRC(IIB:IIE+1,IJS:IJN,:) + PSRC(IIB-1:IIE,IJS:IJN,:)) - &
        (ZDMQ(IIB:IIE+1,IJS:IJN,:) - ZDMQ(IIB-1:IIE,IJS:IJN,:))/6.0
!
   CALL  GET_HALO(ZQL0)
!  
!  WEST BOUND
!
!!$   ZQL0(IIB-1,:,:) = ZQL0(IIE,:,:) JUAN PPMLL01
!
   ZQR0(IIB-1:IIE,IJS:IJN,:) = ZQL0(IIB:IIE+1,IJS:IJN,:)
!
   CALL  GET_HALO(ZQR0)
!
!  EAST BOUND
!
!!$   ZQR0(IIE+1,:,:) = ZQR0(IIB,:,:) JUAN PPMLL01
!
! determine initial coefficients of the parabolae
!
   ZDQ(:,IJS:IJN,:) = ZQR0(:,IJS:IJN,:) - ZQL0(:,IJS:IJN,:)
   ZQ60(:,IJS:IJN,:) = 6.0*(PSRC(:,IJS:IJN,:) - 0.5*(ZQL0(:,IJS:IJN,:) + ZQR0(:,IJS:IJN,:)))
!
! initialize final parabolae parameters
!
   ZQL(:,IJS:IJN,:) = ZQL0(:,IJS:IJN,:)
   ZQR(:,IJS:IJN,:) = ZQR0(:,IJS:IJN,:)
   ZQ6(:,IJS:IJN,:) = ZQ60(:,IJS:IJN,:) 
!
! eliminate over and undershoots and create qL and qR as in Lin96
!
   WHERE ( ZDMQ(:,IJS:IJN,:) == 0.0 )
      ZQL(:,IJS:IJN,:) = PSRC(:,IJS:IJN,:)
      ZQR(:,IJS:IJN,:) = PSRC(:,IJS:IJN,:)
      ZQ6(:,IJS:IJN,:) = 0.0
   ELSEWHERE ( ZQ60(:,IJS:IJN,:)*ZDQ(:,IJS:IJN,:) < -(ZDQ(:,IJS:IJN,:))**2 )
      ZQ6(:,IJS:IJN,:) = 3.0*(ZQL0(:,IJS:IJN,:) - PSRC(:,IJS:IJN,:))
      ZQR(:,IJS:IJN,:) = ZQL0(:,IJS:IJN,:) - ZQ6(:,IJS:IJN,:)
      ZQL(:,IJS:IJN,:) = ZQL0(:,IJS:IJN,:)
   ELSEWHERE ( ZQ60(:,IJS:IJN,:)*ZDQ(:,IJS:IJN,:) > (ZDQ(:,IJS:IJN,:))**2 )
      ZQ6(:,IJS:IJN,:) = 3.0*(ZQR0(:,IJS:IJN,:) - PSRC(:,IJS:IJN,:))
      ZQL(:,IJS:IJN,:) = ZQR0(:,IJS:IJN,:) - ZQ6(:,IJS:IJN,:)
      ZQR(:,IJS:IJN,:) = ZQR0(:,IJS:IJN,:)
   END WHERE
!
! recalculate coefficients of the parabolae
!
   ZDQ(:,IJS:IJN,:) = ZQR(:,IJS:IJN,:) - ZQL(:,IJS:IJN,:)
!
! and finally calculate fluxes for the advection
!
!  ZFPOS(i) = Fct[ ZQR(i-1),PCR(i),ZDQ(i-1),ZQ6(i-1) ]
!
   ZFPOS(IIB:IIE+1,IJS:IJN,:) = ZQR(IIB-1:IIE,IJS:IJN,:) - 0.5*PCR(IIB:IIE+1,IJS:IJN,:) * &            
        (ZDQ(IIB-1:IIE,IJS:IJN,:) - (1.0 - 2.0*PCR(IIB:IIE+1,IJS:IJN,:)/3.0)        &
        * ZQ6(IIB-1:IIE,IJS:IJN,:))
!
   CALL GET_HALO(ZFPOS)
!
!  WEST BOUND
!
! PPOSX(IIB-1,:,:) is not important for the calc of advection so 
! we set it to 0
!!$   ZFPOS(IIB-1,:,:) = 0.0 JUANPPMLL01
!
   ZFNEG(:,IJS:IJN,:) = ZQL(:,IJS:IJN,:) - 0.5*PCR(:,IJS:IJN,:) *      &            
        ( ZDQ(:,IJS:IJN,:) + (1.0 + 2.0*PCR(:,IJS:IJN,:)/3.0) * ZQ6(:,IJS:IJN,:) )
!
   CALL GET_HALO(ZFNEG)
!
! advect the actual field in X direction by U*dt
!
   PR = DXF( PCR*MXM(PRHO)*( ZFPOS*(0.5+SIGN(0.5,PCR)) + & 
                             ZFNEG*(0.5-SIGN(0.5,PCR)) ) )
   CALL GET_HALO(PR)   
!
!
!*       1.2    NON-CYCLIC BOUNDARY CONDITIONS IN THE X DIRECTION 
!               -------------------------------------------------
!
CASE('OPEN')
!
! calculate dmq
!
   ZDMQ = DIF2X(PSRC)
!
! overwrite the values on the boundary to get second order difference
! for qL and qR at the boundary
!
!  WEST BOUND
!
  IF (LWEST_ll()) THEN
   ZDMQ(IIB-1,IJS:IJN,:) = -ZDMQ(IIB,IJS:IJN,:)
  ENDIF
!
!  EAST BOUND
!
  IF (LEAST_ll()) THEN
   ZDMQ(IIE+1,IJS:IJN,:) = -ZDMQ(IIE,IJS:IJN,:)
  ENDIF
!
! monotonize the difference followinq eq. 5 in Lin94
!
!  ZDMQ(i) = Fct[ ZDMQ(i),PSRC(i),PSRC(i-1),PSRC(i+1) ]
!
   ZDMQ(IIB:IIE,IJS:IJN,:) = &
        SIGN( (MIN( ABS(ZDMQ(IIB:IIE,IJS:IJN,:)),2.0*(PSRC(IIB:IIE,IJS:IJN,:) - &
        MIN(PSRC(IIB-1:IIE-1,IJS:IJN,:),PSRC(IIB:IIE,IJS:IJN,:),PSRC(IIB+1:IIE+1,IJS:IJN,:))),    &
        2.0*(MAX(PSRC(IIB-1:IIE-1,IJS:IJN,:),PSRC(IIB:IIE,IJS:IJN,:),PSRC(IIB+1:IIE+1,IJS:IJN,:)) - &
        PSRC(IIB:IIE,IJS:IJN,:)) )), ZDMQ(IIB:IIE,IJS:IJN,:) )
!
!  WEST BOUND
!
!!$   ZDMQ(IIB-1,:,:) = & 
!!$        SIGN( (MIN( ABS(ZDMQ(IIB-1,:,:)), 2.0*(PSRC(IIB-1,:,:) - &
!!$        MIN(PSRC(IIE-1,:,:),PSRC(IIB-1,:,:),PSRC(IIB,:,:))),   &
!!$        2.0*(MAX(PSRC(IIE-1,:,:),PSRC(IIB-1,:,:),PSRC(IIB,:,:)) - &
!!$        PSRC(IIB-1,:,:)) )), ZDMQ(IIB-1,:,:) )
!
!  EAST BOUND
!
!!$   ZDMQ(IIE+1,:,:) = &
!!$        SIGN( (MIN( ABS(ZDMQ(IIE+1,:,:)), 2.0*(PSRC(IIE+1,:,:) - &
!!$        MIN(PSRC(IIE,:,:),PSRC(IIE+1,:,:),PSRC(IIB+1,:,:))),  &
!!$        2.0*(MAX(PSRC(IIE,:,:),PSRC(IIE+1,:,:),PSRC(IIB+1,:,:)) - &
!!$        PSRC(IIE+1,:,:)) )), ZDMQ(IIE+1,:,:) )
!
!
!  update ZDMQ HALO before next/further  utilisation 
!
   CALL  GET_HALO(ZDMQ)
!
! calculate qL and qR
!
!  ZQL0(i) = Fct[ PSRC(i),PSRC(i-1),ZDMQ(i),ZDMQ(i-1) ]
!
   ZQL0(IIB:IIE+1,IJS:IJN,:) = 0.5*(PSRC(IIB:IIE+1,IJS:IJN,:) + PSRC(IIB-1:IIE,IJS:IJN,:)) - &
        (ZDMQ(IIB:IIE+1,IJS:IJN,:) - ZDMQ(IIB-1:IIE,IJS:IJN,:))/6.0
!
   CALL  GET_HALO(ZQL0)
!  
!  WEST BOUND
!
  IF (LWEST_ll()) THEN
   ZQL0(IIB-1,IJS:IJN,:) = ZQL0(IIB,IJS:IJN,:)
  ENDIF
!
   ZQR0(IIB-1:IIE,IJS:IJN,:) = ZQL0(IIB:IIE+1,IJS:IJN,:)
!
   CALL  GET_HALO(ZQR0)
!
!  EAST BOUND
!
  IF (LEAST_ll()) THEN
   ZQR0(IIE+1,IJS:IJN,:) = ZQR0(IIE,IJS:IJN,:)
  ENDIF
!
! determine initial coefficients of the parabolae
!
   ZDQ(:,IJS:IJN,:) = ZQR0(:,IJS:IJN,:) - ZQL0(:,IJS:IJN,:)
   ZQ60(:,IJS:IJN,:) = 6.0*(PSRC(:,IJS:IJN,:) - 0.5*(ZQL0(:,IJS:IJN,:) + ZQR0(:,IJS:IJN,:)))
!
! initialize final parabolae parameters
!
   ZQL(:,IJS:IJN,:) = ZQL0(:,IJS:IJN,:)
   ZQR(:,IJS:IJN,:) = ZQR0(:,IJS:IJN,:)
   ZQ6(:,IJS:IJN,:) = ZQ60(:,IJS:IJN,:)
!
! eliminate over and undershoots and create qL and qR as in Lin96
!
   WHERE ( ZDMQ(:,IJS:IJN,:) == 0.0 )
      ZQL(:,IJS:IJN,:) = PSRC(:,IJS:IJN,:)
      ZQR(:,IJS:IJN,:) = PSRC(:,IJS:IJN,:)
      ZQ6(:,IJS:IJN,:) = 0.0
   ELSEWHERE ( ZQ60(:,IJS:IJN,:)*ZDQ(:,IJS:IJN,:) < -(ZDQ(:,IJS:IJN,:))**2 )
      ZQ6(:,IJS:IJN,:) = 3.0*(ZQL0(:,IJS:IJN,:) - PSRC(:,IJS:IJN,:))
      ZQR(:,IJS:IJN,:) = ZQL0(:,IJS:IJN,:) - ZQ6(:,IJS:IJN,:)
      ZQL(:,IJS:IJN,:) = ZQL0(:,IJS:IJN,:)
   ELSEWHERE ( ZQ60(:,IJS:IJN,:)*ZDQ(:,IJS:IJN,:) > (ZDQ(:,IJS:IJN,:))**2 )
      ZQ6(:,IJS:IJN,:) = 3.0*(ZQR0(:,IJS:IJN,:) - PSRC(:,IJS:IJN,:))
      ZQL(:,IJS:IJN,:) = ZQR0(:,IJS:IJN,:) - ZQ6(:,IJS:IJN,:)
      ZQR(:,IJS:IJN,:) = ZQR0(:,IJS:IJN,:)
   END WHERE
!
! recalculate coefficients of the parabolae
!
   ZDQ(:,IJS:IJN,:) = ZQR(:,IJS:IJN,:) - ZQL(:,IJS:IJN,:)
!
! and finally calculate fluxes for the advection
!
!
!  ZFPOS(i) = Fct[ ZQR(i-1),PCR(i),ZDQ(i-1),ZQ6(i-1) ]
!
!!$   ZFPOS(IIB+1:IIE+1,:,:) = ZQR(IIB:IIE,:,:) - 0.5*PCR(IIB+1:IIE+1,:,:) * &            
!!$        (ZDQ(IIB:IIE,:,:) - (1.0 - 2.0*PCR(IIB+1:IIE+1,:,:)/3.0)          &
!!$        * ZQ6(IIB:IIE,:,:))
   ZFPOS(IIB:IIE+1,IJS:IJN,:) = ZQR(IIB-1:IIE,IJS:IJN,:) - 0.5*PCR(IIB:IIE+1,IJS:IJN,:) * &            
        (ZDQ(IIB-1:IIE,IJS:IJN,:) - (1.0 - 2.0*PCR(IIB:IIE+1,IJS:IJN,:)/3.0)        &
        * ZQ6(IIB-1:IIE,IJS:IJN,:))
!
   CALL GET_HALO(ZFPOS)
!
!
!  WEST BOUND
!
! advection flux at open boundary when u(IIB) > 0
! 
  IF (LWEST_ll()) THEN 
   ZFPOS(IIB,IJS:IJN,:) = (PSRC(IIB-1,IJS:IJN,:) - ZQR(IIB-1,IJS:IJN,:))*PCR(IIB,IJS:IJN,:) + &
                    ZQR(IIB-1,IJS:IJN,:)
! PPOSX(IIB-1,:,:) is not important for the calc of advection so 
! we set it to 0
!!$   ZFPOS(IIB-1,:,:) = 0.0
   ENDIF
!
!!$   ZFNEG(IIB-1:IIE,:,:) = ZQL(IIB-1:IIE,:,:) - 0.5*PCR(IIB-1:IIE,:,:) *  &            
!!$        (ZDQ(IIB-1:IIE,:,:) + (1.0 + 2.0*PCR(IIB-1:IIE,:,:)/3.0)         &
!!$        * ZQ6(IIB-1:IIE,:,:))
   ZFNEG(:,IJS:IJN,:) = ZQL(:,IJS:IJN,:) - 0.5*PCR(:,IJS:IJN,:) *      &            
        ( ZDQ(:,IJS:IJN,:) + (1.0 + 2.0*PCR(:,IJS:IJN,:)/3.0) * ZQ6(:,IJS:IJN,:) )
!
   CALL GET_HALO(ZFNEG)
!
!  EAST BOUND
!
! advection flux at open boundary when u(IIE+1) < 0
  IF (LEAST_ll()) THEN
   ZFNEG(IIE+1,IJS:IJN,:) = (ZQR(IIE,IJS:IJN,:)-PSRC(IIE+1,IJS:IJN,:))*PCR(IIE+1,IJS:IJN,:) + &
                      ZQR(IIE,IJS:IJN,:)
  ENDIF
!
! advect the actual field in X direction by U*dt
!
   PR = DXF( PCR*MXM(PRHO)*( ZFPOS*(0.5+SIGN(0.5,PCR)) + & 
                             ZFNEG*(0.5-SIGN(0.5,PCR)) ) )
   CALL GET_HALO(PR)   
!
!
END SELECT
!
CONTAINS
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION DIF2X(PQ) RESULT(DQ)
!     ########################################################################
!!
!!****  DIF2X - leap-frog difference operator in X direction
!!
!!    Calculates the difference assuming periodic BC (CYCL). 
!!
!!    DQ(I) = 0.5 * (PQ(I+1) - PQ(I-1))
!!
!!    MODIFICATIONS
!!    -------------
!!    
!!    18.3.2006.  T. Maric - original version
!!    07/2010     J.Escobar : Correction for reproducility
!!    04/2017     J.Escobar : initialize realistic value in all HALO pts
!-------------------------------------------------------------------------------
!
!
USE MODE_ll
!
IMPLICIT NONE
! 
!*       0.1   Declarations of dummy arguments :
!   
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PQ
REAL, DIMENSION(SIZE(PQ,1),SIZE(PQ,2),SIZE(PQ,3)) :: DQ
!
!*       0.2   Declarations of local variables :
!   
INTEGER :: IIB,IJB        ! Begining useful area in x,y directions
INTEGER :: IIE,IJE        ! End useful area in x,y directions
!
!-------------------------------------------------------------------------------
!
!*       1.0.     COMPUTE THE DOMAIN DIMENSIONS
!                 -----------------------------
!
!!$CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IIB=2 ; IIE = SIZE(PQ,1) -1
IJB=2 ; IJE = SIZE(PQ,2) -1
!
!-------------------------------------------------------------------------------
!
!*       2.0.     COMPUTE THE DIFFERENCE
!                 ----------------------
!   
DQ(IIB:IIE,:,:) = PQ(IIB+1:IIE+1,:,:) - PQ(IIB-1:IIE-1,:,:)
DQ(IIB-1,:,:) = PQ(IIB,:,:) - PQ(IIE-1,:,:)
DQ(IIE+1,:,:) = PQ(IIB+1,:,:) - PQ(IIE,:,:)  
DQ = 0.5*DQ  
!
END FUNCTION DIF2X
!
END FUNCTION PPM_01_X
!
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION PPM_01_Y(HLBCY, KGRID, PSRC, PCR, PRHO, PTSTEP) &
               RESULT(PR)
!     ########################################################################
!!
!!****  PPM_01_Y - PPM_01 fully monotonic PPM advection scheme in Y direction
!!                 Colella notation
!!
!!    MODIFICATIONS
!!    -------------
!!    
!!    11.5.2006.  T. Maric - original version
!!      J.Escobar 21/03/2013: for HALOK comment all NHALO=1 test
!!      J.Escobar 28/06/2018: limit computation on TAB(IIW:IIA,:,:) to avoid unneeded NaN 
!!      J.Escobr  16/07/2018: still NaN pb => reintroduce initialization of temporary local array
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODE_IO_ll
USE MODI_SHUMAN
USE MODI_GET_HALO
!
USE MODD_CONF
!USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY  ! Y direction LBC type
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER:: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER:: IIE,IJE    ! End useful area in x,y,z directions
!
! terms used in parabolic interpolation, dmq, qL, qR, dq, q6
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZQL,ZQR
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZDQ,ZQ6
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZDMQ
!
! extra variables for the initial guess of parabolae parameters
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZQL0,ZQR0,ZQ60
!
! advection fluxes
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZFPOS, ZFNEG
!
!BEG JUAN PPM_LL
INTEGER                          :: IIW,IIA
!END JUAN PPM_LL
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IIW=IIB
IIA=IIE
CALL GET_HALO(PSRC)

!
!-------------------------------------------------------------------------------
!
!
PR=PSRC
ZQL=PSRC
ZQR=PSRC
ZDQ=PSRC
ZQ6=PSRC
ZDMQ=PSRC
ZQL0=PSRC
ZQR0=PSRC
ZQ60=PSRC
ZFPOS=PSRC
ZFNEG=PSRC
!
SELECT CASE ( HLBCY(1) ) ! Y direction LBC type: (1) for left side
!
!*       2.1    CYCLIC BOUNDARY CONDITIONS IN THE Y DIRECTION
!               ---------------------------------------------
!
CASE ('CYCL','WALL')          ! In that case one must have HLBCY(1) == HLBCY(2)
!
! calculate dmq
   ZDMQ = DIF2Y(PSRC)
!
! monotonize the difference followinq eq. 5 in Lin94
!BEG JUAN PPM_LL01
!
!  ZDMQ(j) = Fct[ ZDMQ(j),PSRC(j),PSRC(j-1),PSRC(j+1) ]
!
   ZDMQ(IIW:IIA,IJB:IJE,:) = &
        SIGN( (MIN( ABS(ZDMQ(IIW:IIA,IJB:IJE,:)),2.0*(PSRC(IIW:IIA,IJB:IJE,:) - &
        MIN(PSRC(IIW:IIA,IJB-1:IJE-1,:),PSRC(IIW:IIA,IJB:IJE,:),PSRC(IIW:IIA,IJB+1:IJE+1,:))),    &
        2.0*(MAX(PSRC(IIW:IIA,IJB-1:IJE-1,:),PSRC(IIW:IIA,IJB:IJE,:),PSRC(IIW:IIA,IJB+1:IJE+1,:)) - &
        PSRC(IIW:IIA,IJB:IJE,:)) )), ZDMQ(IIW:IIA,IJB:IJE,:) )
!
!  SOUTH BOUND
!
!!$   ZDMQ(:,IJB-1,:) = & 
!!$        SIGN( (MIN( ABS(ZDMQ(:,IJB-1,:)), 2.0*(PSRC(:,IJB-1,:) - &
!!$        MIN(PSRC(:,IJE-1,:),PSRC(:,IJB-1,:),PSRC(:,IJB,:))),   &
!!$        2.0*(MAX(PSRC(:,IJE-1,:),PSRC(:,IJB-1,:),PSRC(:,IJB,:)) - &
!!$        PSRC(:,IJB-1,:)) )), ZDMQ(:,IJB-1,:) )
!
!  NORTH BOUND
!
!!$   ZDMQ(:,IJE+1,:) = &
!!$        SIGN( (MIN( ABS(ZDMQ(:,IJE+1,:)), 2.0*(PSRC(:,IJE+1,:) - &
!!$        MIN(PSRC(:,IJE,:),PSRC(:,IJE+1,:),PSRC(:,IJB+1,:))),  &
!!$        2.0*(MAX(PSRC(:,IJE,:),PSRC(:,IJE+1,:),PSRC(:,IJB+1,:)) - &
!!$        PSRC(:,IJE+1,:)) )), ZDMQ(:,IJE+1,:) )   
!
!  update ZDMQ HALO before next/further  utilisation 
!
   CALL  GET_HALO(ZDMQ)
!
! calculate qL and qR with the modified dmq
!
   ZQL0(IIW:IIA,IJB:IJE+1,:) = 0.5*(PSRC(IIW:IIA,IJB:IJE+1,:) + PSRC(IIW:IIA,IJB-1:IJE,:)) - &
        (ZDMQ(IIW:IIA,IJB:IJE+1,:) - ZDMQ(IIW:IIA,IJB-1:IJE,:))/6.0
!
   CALL  GET_HALO(ZQL0)
!  
!  SOUTH BOUND
!
!!$   ZQL0(:,IJB-1,:) = ZQL0(:,IJE,:) JUAN PPMLL01
!
   ZQR0(IIW:IIA,IJB-1:IJE,:) = ZQL0(IIW:IIA,IJB:IJE+1,:)
!
   CALL  GET_HALO(ZQR0)
!
!  NORTH BOUND
!
!!$   ZQR0(:,IJE+1,:) = ZQR0(:,IJB,:) JUAN PPMLL01
!
! determine initial coefficients of the parabolae
!
   ZDQ(IIW:IIA,:,:) = ZQR0(IIW:IIA,:,:) - ZQL0(IIW:IIA,:,:)
   ZQ60(IIW:IIA,:,:) = 6.0*(PSRC(IIW:IIA,:,:) - 0.5*(ZQL0(IIW:IIA,:,:) + ZQR0(IIW:IIA,:,:)))
!
! initialize final parabolae parameters
!
   ZQL(IIW:IIA,:,:) = ZQL0(IIW:IIA,:,:)
   ZQR(IIW:IIA,:,:) = ZQR0(IIW:IIA,:,:)
   ZQ6(IIW:IIA,:,:) = ZQ60(IIW:IIA,:,:) 
!
! eliminate over and undershoots and create qL and qR as in Lin96
!
   WHERE ( ZDMQ(IIW:IIA,:,:) == 0.0 )
      ZQL(IIW:IIA,:,:) = PSRC(IIW:IIA,:,:)
      ZQR(IIW:IIA,:,:) = PSRC(IIW:IIA,:,:)
      ZQ6(IIW:IIA,:,:) = 0.0
   ELSEWHERE ( ZQ60(IIW:IIA,:,:)*ZDQ(IIW:IIA,:,:) < -(ZDQ(IIW:IIA,:,:))**2 )
      ZQ6(IIW:IIA,:,:) = 3.0*(ZQL0(IIW:IIA,:,:) - PSRC(IIW:IIA,:,:))
      ZQR(IIW:IIA,:,:) = ZQL0(IIW:IIA,:,:) - ZQ6(IIW:IIA,:,:)
      ZQL(IIW:IIA,:,:) = ZQL0(IIW:IIA,:,:)
   ELSEWHERE ( ZQ60(IIW:IIA,:,:)*ZDQ(IIW:IIA,:,:) > (ZDQ(IIW:IIA,:,:))**2 )
      ZQ6(IIW:IIA,:,:) = 3.0*(ZQR0(IIW:IIA,:,:) - PSRC(IIW:IIA,:,:))
      ZQL(IIW:IIA,:,:) = ZQR0(IIW:IIA,:,:) - ZQ6(IIW:IIA,:,:)
      ZQR(IIW:IIA,:,:) = ZQR0(IIW:IIA,:,:)
   END WHERE
!
! recalculate coefficients of the parabolae
!
   ZDQ(IIW:IIA,:,:) = ZQR(IIW:IIA,:,:) - ZQL(IIW:IIA,:,:)
!
! and finally calculate fluxes for the advection
!
!  ZFPOS(j) = Fct[ ZQR(j-1),PCR(i),ZDQ(j-1),ZQ6(j-1) ]
!
   ZFPOS(IIW:IIA,IJB:IJE+1,:) = ZQR(IIW:IIA,IJB-1:IJE,:) - 0.5*PCR(IIW:IIA,IJB:IJE+1,:) * &            
        (ZDQ(IIW:IIA,IJB-1:IJE,:) - (1.0 - 2.0*PCR(IIW:IIA,IJB:IJE+1,:)/3.0)        &
        * ZQ6(IIW:IIA,IJB-1:IJE,:))
!
   CALL GET_HALO(ZFPOS)
!
! SOUTH BOUND
!
! PPOSX(:,IJB-1,:) is not important for the calc of advection so 
! we set it to 0
!!$   ZFPOS(:,IJB-1,:) = 0.0 JUANPPMLL01
!
   ZFNEG(IIW:IIA,:,:) = ZQL(IIW:IIA,:,:) - 0.5*PCR(IIW:IIA,:,:) *      &            
        ( ZDQ(IIW:IIA,:,:) + (1.0 + 2.0*PCR(IIW:IIA,:,:)/3.0) * ZQ6(IIW:IIA,:,:) )
!
   CALL GET_HALO(ZFNEG)
!
! advect the actual field in Y direction by V*dt
!
   PR = DYF( PCR*MYM(PRHO)*( ZFPOS*(0.5+SIGN(0.5,PCR)) + & 
                             ZFNEG*(0.5-SIGN(0.5,PCR)) ) )
   CALL GET_HALO(PR) 
!
!*       2.2    NON-CYCLIC BOUNDARY CONDITIONS IN THE Y DIRECTION
!               -------------------------------------------------
!
CASE('OPEN')
!
! calculate dmq
   ZDMQ = DIF2Y(PSRC)
! overwrite the values on the boundary to get second order difference
! for qL and qR at the boundary
!
!  SOUTH BOUND
!
   IF (LSOUTH_ll()) THEN
    ZDMQ(IIW:IIA,IJB-1,:) = -ZDMQ(IIW:IIA,IJB,:)
   ENDIF
!
!  NORTH BOUND
!
   IF (LNORTH_ll()) THEN
    ZDMQ(IIW:IIA,IJE+1,:) = -ZDMQ(IIW:IIA,IJE,:)
   ENDIF
!
! monotonize the difference followinq eq. 5 in Lin94
   ZDMQ(IIW:IIA,IJB:IJE,:) = &
        SIGN( (MIN( ABS(ZDMQ(IIW:IIA,IJB:IJE,:)),2.0*(PSRC(IIW:IIA,IJB:IJE,:) - &
        MIN(PSRC(IIW:IIA,IJB-1:IJE-1,:),PSRC(IIW:IIA,IJB:IJE,:),PSRC(IIW:IIA,IJB+1:IJE+1,:))),    &
        2.0*(MAX(PSRC(IIW:IIA,IJB-1:IJE-1,:),PSRC(IIW:IIA,IJB:IJE,:),PSRC(IIW:IIA,IJB+1:IJE+1,:)) - &
        PSRC(IIW:IIA,IJB:IJE,:)) )), ZDMQ(IIW:IIA,IJB:IJE,:) )
!!$   ZDMQ(:,IJB-1,:) = & 
!!$        SIGN( (MIN( ABS(ZDMQ(:,IJB-1,:)), 2.0*(PSRC(:,IJB-1,:) - &
!!$        MIN(PSRC(:,IJE-1,:),PSRC(:,IJB-1,:),PSRC(:,IJB,:))),   &
!!$        2.0*(MAX(PSRC(:,IJE-1,:),PSRC(:,IJB-1,:),PSRC(:,IJB,:)) - &
!!$        PSRC(:,IJB-1,:)) )), ZDMQ(:,IJB-1,:) )
!!$   ZDMQ(:,IJE+1,:) = &
!!$        SIGN( (MIN( ABS(ZDMQ(:,IJE+1,:)), 2.0*(PSRC(:,IJE+1,:) - &
!!$        MIN(PSRC(:,IJE,:),PSRC(:,IJE+1,:),PSRC(:,IJB+1,:))),  &
!!$        2.0*(MAX(PSRC(:,IJE,:),PSRC(:,IJE+1,:),PSRC(:,IJB+1,:)) - &
!!$        PSRC(:,IJE+1,:)) )), ZDMQ(:,IJE+1,:) ) 
!
!  update ZDMQ HALO before next/further  utilisation 
!
   CALL  GET_HALO(ZDMQ)  
!
! calculate qL and qR with the modified dmq
!
   ZQL0(IIW:IIA,IJB:IJE+1,:) = 0.5*(PSRC(IIW:IIA,IJB:IJE+1,:) + PSRC(IIW:IIA,IJB-1:IJE,:)) - &
        (ZDMQ(IIW:IIA,IJB:IJE+1,:) - ZDMQ(IIW:IIA,IJB-1:IJE,:))/6.0
!
   CALL  GET_HALO(ZQL0)
!  
!  SOUTH BOUND
!
   IF (LSOUTH_ll()) THEN
    ZQL0(IIW:IIA,IJB-1,:) = ZQL0(IIW:IIA,IJB,:)
   ENDIF
!
   ZQR0(IIW:IIA,IJB-1:IJE,:) = ZQL0(IIW:IIA,IJB:IJE+1,:)
!
!  NORTH BOUND
!
   IF (LNORTH_ll()) THEN
    ZQR0(IIW:IIA,IJE+1,:) = ZQR0(IIW:IIA,IJE,:)
   ENDIF
!
! determine initial coefficients of the parabolae
!
   ZDQ(IIW:IIA,:,:) = ZQR0(IIW:IIA,:,:) - ZQL0(IIW:IIA,:,:)
   ZQ60(IIW:IIA,:,:) = 6.0*(PSRC(IIW:IIA,:,:) - 0.5*(ZQL0(IIW:IIA,:,:) + ZQR0(IIW:IIA,:,:)))
!
! initialize final parabolae parameters
!
   ZQL(IIW:IIA,:,:) = ZQL0(IIW:IIA,:,:)
   ZQR(IIW:IIA,:,:) = ZQR0(IIW:IIA,:,:)
   ZQ6(IIW:IIA,:,:) = ZQ60(IIW:IIA,:,:) 
!
! eliminate over and undershoots and create qL and qR as in Lin96
!
   WHERE ( ZDMQ(IIW:IIA,:,:) == 0.0 )
      ZQL(IIW:IIA,:,:) = PSRC(IIW:IIA,:,:)
      ZQR(IIW:IIA,:,:) = PSRC(IIW:IIA,:,:)
      ZQ6(IIW:IIA,:,:) = 0.0
   ELSEWHERE ( ZQ60(IIW:IIA,:,:)*ZDQ(IIW:IIA,:,:) < -(ZDQ(IIW:IIA,:,:))**2 )
      ZQ6(IIW:IIA,:,:) = 3.0*(ZQL0(IIW:IIA,:,:) - PSRC(IIW:IIA,:,:))
      ZQR(IIW:IIA,:,:) = ZQL0(IIW:IIA,:,:) - ZQ6(IIW:IIA,:,:)
      ZQL(IIW:IIA,:,:) = ZQL0(IIW:IIA,:,:)
   ELSEWHERE ( ZQ60(IIW:IIA,:,:)*ZDQ(IIW:IIA,:,:) > (ZDQ(IIW:IIA,:,:))**2 )
      ZQ6(IIW:IIA,:,:) = 3.0*(ZQR0(IIW:IIA,:,:) - PSRC(IIW:IIA,:,:))
      ZQL(IIW:IIA,:,:) = ZQR0(IIW:IIA,:,:) - ZQ6(IIW:IIA,:,:)
      ZQR(IIW:IIA,:,:) = ZQR0(IIW:IIA,:,:)
   END WHERE
!
! recalculate coefficients of the parabolae
!
   ZDQ(IIW:IIA,:,:) = ZQR(IIW:IIA,:,:) - ZQL(IIW:IIA,:,:)
!
! and finally calculate fluxes for the advection
!!$   ZFPOS(:,IJB+1:IJE+1,:) = ZQR(:,IJB:IJE,:) - 0.5*PCR(:,IJB+1:IJE+1,:) * &            
!!$        (ZDQ(:,IJB:IJE,:) - (1.0 - 2.0*PCR(:,IJB+1:IJE+1,:)/3.0)        &
!!$        * ZQ6(:,IJB:IJE,:))
   ZFPOS(IIW:IIA,IJB:IJE+1,:) = ZQR(IIW:IIA,IJB-1:IJE,:) - 0.5*PCR(IIW:IIA,IJB:IJE+1,:) * &            
        (ZDQ(IIW:IIA,IJB-1:IJE,:) - (1.0 - 2.0*PCR(IIW:IIA,IJB:IJE+1,:)/3.0)        &
        * ZQ6(IIW:IIA,IJB-1:IJE,:))
!
   CALL GET_HALO(ZFPOS)
!
!
! advection flux at open boundary when u(IJB) > 0
!  
!  SOUTH BOUND
!
   IF (LSOUTH_ll()) THEN
    ZFPOS(IIW:IIA,IJB,:) = (PSRC(IIW:IIA,IJB-1,:) - ZQR(IIW:IIA,IJB-1,:))*PCR(IIW:IIA,IJB,:) + &
                      ZQR(IIW:IIA,IJB-1,:)
   ENDIF
!
! PPOSX(:,IJB-1,:) is not important for the calc of advection so 
! we set it to 0
!!$   ZFPOS(:,IJB-1,:) = 0.0 ! JUAN PPMLL01
!
!!$   ZFNEG(:,IJB-1:IJE,:) = ZQL(:,IJB-1:IJE,:) - 0.5*PCR(:,IJB-1:IJE,:) * &            
!!$        ( ZDQ(:,IJB-1:IJE,:) + (1.0 + 2.0*PCR(:,IJB-1:IJE,:)/3.0) * &
!!$        ZQ6(:,IJB-1:IJE,:) )
   ZFNEG(IIW:IIA,:,:) = ZQL(IIW:IIA,:,:) - 0.5*PCR(IIW:IIA,:,:) *      &            
        ( ZDQ(IIW:IIA,:,:) + (1.0 + 2.0*PCR(IIW:IIA,:,:)/3.0) * ZQ6(IIW:IIA,:,:) )
!
   CALL GET_HALO(ZFNEG)
!
! advection flux at open boundary when u(IJE+1) < 0
!
!  NORTH BOUND
!
   IF (LNORTH_ll()) THEN
    ZFNEG(IIW:IIA,IJE+1,:) = (ZQR(IIW:IIA,IJE,:)-PSRC(IIW:IIA,IJE+1,:))*PCR(IIW:IIA,IJE+1,:) + &
                        ZQR(IIW:IIA,IJE,:)
   ENDIF
!
! advect the actual field in X direction by U*dt
!
   PR = DYF( PCR*MYM(PRHO)*( ZFPOS*(0.5+SIGN(0.5,PCR)) + & 
                             ZFNEG*(0.5-SIGN(0.5,PCR)) ) )
!
   CALL GET_HALO(PR)
!
!
END SELECT
!
CONTAINS
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION DIF2Y(PQ) RESULT(DQ)
!     ########################################################################
!!
!!****  DIF2Y - leap-frog difference operator in Y direction
!!
!!    Calculates the difference assuming periodic BC (CYCL). 
!!
!!    DQ(J) = 0.5 * (PQ(J+1) - PQ(J-1))
!!
!!
!!    MODIFICATIONS
!!    -------------
!!    
!!    18.3.2006.  T. Maric - original version, works only for periodic boundary
!!                           conditions and on one domain
!!    04/2017     J.Escobar : initialize realistic value in all HALO pts
!!
!-------------------------------------------------------------------------------
!
!
USE MODE_ll
!
IMPLICIT NONE
! 
!*       0.1   Declarations of dummy arguments :
!   
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PQ
REAL, DIMENSION(SIZE(PQ,1),SIZE(PQ,2),SIZE(PQ,3)) :: DQ
!
!*       0.2   Declarations of local variables :
!   
INTEGER :: IIB,IJB        ! Begining useful area in x,y directions
INTEGER :: IIE,IJE        ! End useful area in x,y directions
!
!-------------------------------------------------------------------------------
!
!*       1.0.     COMPUTE THE DOMAIN DIMENSIONS
!                 -----------------------------
!
!!$CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IIB=2 ; IIE = SIZE(PQ,1) -1
IJB=2 ; IJE = SIZE(PQ,2) -1
!
!-------------------------------------------------------------------------------
!
!*       2.0.     COMPUTE THE DIFFERENCE
!                 ----------------------
!
DQ(:,IJB:IJE,:) = PQ(:,IJB+1:IJE+1,:) - PQ(:,IJB-1:IJE-1,:)
DQ(:,IJB-1,:) = PQ(:,IJB,:) - PQ(:,IJE-1,:)
DQ(:,IJE+1,:) = PQ(:,IJB+1,:) - PQ(:,IJE,:) 
DQ = 0.5 * DQ   
!
END FUNCTION DIF2Y
!
END FUNCTION PPM_01_Y
!
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION PPM_01_Z(KGRID, PSRC, PCR, PRHO, PTSTEP) RESULT(PR)
!     ########################################################################
!!
!!****  PPM_01_Z - PPM_01 fully monotonic PPM advection scheme in Z direction
!!                 Colella notation
!!
!!    MODIFICATIONS
!!    -------------
!!    
!!    11.5.2006.  T. Maric - original version
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODI_SHUMAN
USE MODI_GET_HALO
!
USE MODD_CONF
USE MODD_PARAMETERS
!USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER:: IKB    ! Begining useful area in x,y,z directions
INTEGER:: IKE    ! End useful area in x,y,z directions
INTEGER:: IKU
!
! terms used in parabolic interpolation, dmq, qL, qR, dq, q6
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZQL,ZQR
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZDQ,ZQ6
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZDMQ
!
! extra variables for the initial guess of parabolae parameters
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZQL0,ZQR0,ZQ60
!
! advection fluxes
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZFPOS, ZFNEG
!
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
IKB = 1 + JPVEXT
IKE = SIZE(PSRC,3) - JPVEXT
IKU = SIZE(PSRC,3)
!
!-------------------------------------------------------------------------------
!
!*       3.     PPM ADVECTION IN THE Z DIRECTION
!               --------------------------------
! 
! calculate dmq
ZDMQ = DIF2Z(PSRC)
!
! monotonize the difference followinq eq. 5 in Lin94
! use the periodic BC here, it doesn't matter for vertical (hopefully) 
!
ZDMQ(:,:,IKB:IKE) = &
     SIGN( (MIN( ABS(ZDMQ(:,:,IKB:IKE)),2.0*(PSRC(:,:,IKB:IKE) - &
     MIN(PSRC(:,:,IKB-1:IKE-1),PSRC(:,:,IKB:IKE),PSRC(:,:,IKB+1:IKE+1))),    &
     2.0*(MAX(PSRC(:,:,IKB-1:IKE-1),PSRC(:,:,IKB:IKE),PSRC(:,:,IKB+1:IKE+1)) - &
     PSRC(:,:,IKB:IKE)) )), ZDMQ(:,:,IKB:IKE) )
ZDMQ(:,:,IKB-1) = & 
     SIGN( (MIN( ABS(ZDMQ(:,:,IKB-1)), 2.0*(PSRC(:,:,IKB-1) - &
     MIN(PSRC(:,:,IKE-1),PSRC(:,:,IKB-1),PSRC(:,:,IKB))),   &
     2.0*(MAX(PSRC(:,:,IKE-1),PSRC(:,:,IKB-1),PSRC(:,:,IKB)) - &
     PSRC(:,:,IKB-1)) )), ZDMQ(:,:,IKB-1) )
ZDMQ(:,:,IKE+1) = &
     SIGN( (MIN( ABS(ZDMQ(:,:,IKE+1)), 2.0*(PSRC(:,:,IKE+1) - &
     MIN(PSRC(:,:,IKE),PSRC(:,:,IKE+1),PSRC(:,:,IKB+1))),  &
     2.0*(MAX(PSRC(:,:,IKE),PSRC(:,:,IKE+1),PSRC(:,:,IKB+1)) - &
     PSRC(:,:,IKE+1)) )), ZDMQ(:,:,IKE+1) )
!
! calculate qL and qR with the modified dmq
!
ZQL0(:,:,IKB:IKE+1) = 0.5*(PSRC(:,:,IKB:IKE+1) + PSRC(:,:,IKB-1:IKE)) - &
     (ZDMQ(:,:,IKB:IKE+1) - ZDMQ(:,:,IKB-1:IKE))/6.0
ZQL0(:,:,IKB-1) = ZQL0(:,:,IKE)
!
ZQR0(:,:,IKB-1:IKE) = ZQL0(:,:,IKB:IKE+1)
ZQR0(:,:,IKE+1) = ZQR0(:,:,IKB)
!
! determine initial coefficients of the parabolae
!
ZDQ = ZQR0 - ZQL0
ZQ60 = 6.0*(PSRC - 0.5*(ZQL0 + ZQR0))
!
! initialize final parabolae parameters
!
ZQL = ZQL0
ZQR = ZQR0
ZQ6 = ZQ60 
!
! eliminate over and undershoots and create qL and qR as in Lin96
!
WHERE ( ZDMQ == 0.0 )
   ZQL = PSRC
   ZQR = PSRC
   ZQ6 = 0.0
ELSEWHERE ( ZQ60*ZDQ < -(ZDQ)**2 )
   ZQ6 = 3.0*(ZQL0 - PSRC)
   ZQR = ZQL0 - ZQ6
   ZQL = ZQL0
ELSEWHERE ( ZQ60*ZDQ > (ZDQ)**2 ) 
   ZQ6 = 3.0*(ZQR0 - PSRC)
   ZQL = ZQR0 - ZQ6
   ZQR = ZQR0
END WHERE
!
! recalculate coefficients of the parabolae
!
ZDQ = ZQR - ZQL
!
! and finally calculate fluxes for the advection
!
ZFPOS(:,:,IKB+1:IKE+1) = ZQR(:,:,IKB:IKE) - 0.5*PCR(:,:,IKB+1:IKE+1) * &            
     (ZDQ(:,:,IKB:IKE) - (1.0 - 2.0*PCR(:,:,IKB+1:IKE+1)/3.0)        &
     * ZQ6(:,:,IKB:IKE))
!
! advection flux at open boundary when u(IKB) > 0
ZFPOS(:,:,IKB) = (PSRC(:,:,IKB-1) - ZQR(:,:,IKB-1))*PCR(:,:,IKB) + &
                 ZQR(:,:,IKB-1)
!
! PPOSX(IKB-1) is not important for the calc of advection so 
! we set it to 0
ZFPOS(:,:,IKB-1) = 0.0
!
ZFNEG(:,:,IKB-1:IKE) = ZQL(:,:,IKB-1:IKE) - 0.5*PCR(:,:,IKB-1:IKE) *      &            
     ( ZDQ(:,:,IKB-1:IKE) + (1.0 + 2.0*PCR(:,:,IKB-1:IKE)/3.0) * &
       ZQ6(:,:,IKB-1:IKE) )
!
! advection flux at open boundary when u(IKE+1) < 0
ZFNEG(:,:,IKE+1) = (ZQR(:,:,IKE)-PSRC(:,:,IKE+1))*PCR(:,:,IKE+1) + &
                   ZQR(:,:,IKE)
!
! advect the actual field in Z direction by W*dt
!
PR = DZF(1,IKU,1, PCR*MZM(1,IKU,1,PRHO)*( ZFPOS*(0.5+SIGN(0.5,PCR)) + & 
                          ZFNEG*(0.5-SIGN(0.5,PCR)) ) )
CALL GET_HALO(PR)
!
CONTAINS
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION DIF2Z(PQ) RESULT(DQ)
!     ########################################################################
!!
!!****  DIF2Z - leap-frog difference operator in Z direction
!!
!!    Calculates the difference assuming periodic BC (CYCL). 
!!
!!    DQ(K) = 0.5 * (PQ(K+1) - PQ(K-1))
!!
!!
!!    MODIFICATIONS
!!    -------------
!!    
!!    18.3.2006.  T. Maric - original version
!!
!-------------------------------------------------------------------------------
!
!
USE MODE_ll
USE MODD_CONF
USE MODD_PARAMETERS
!
IMPLICIT NONE
! 
!*       0.1   Declarations of dummy arguments :
!   
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PQ
REAL, DIMENSION(SIZE(PQ,1),SIZE(PQ,2),SIZE(PQ,3)) :: DQ
!
!*       0.2   Declarations of local variables :
!   
INTEGER :: IKB    ! Begining useful area in z directions
INTEGER :: IKE    ! End useful area in z directions
!
!-------------------------------------------------------------------------------
!
!*       1.0.     COMPUTE THE DOMAIN DIMENSIONS
!                 -----------------------------
!
!CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PQ,3) - JPVEXT
!
!-------------------------------------------------------------------------------
!
!*       2.0.     COMPUTE THE DIFFERENCE
!                 ----------------------
!   
DQ(:,:,IKB:IKE) = PQ(:,:,IKB+1:IKE+1) - PQ(:,:,IKB-1:IKE-1)
DQ(:,:,IKB-1) = -DQ(:,:,IKB)
DQ(:,:,IKE+1) = -DQ(:,:,IKE)
DQ = 0.5 * DQ
!
END FUNCTION DIF2Z
!
END FUNCTION PPM_01_Z
!
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION PPM_S0_X(HLBCX, KGRID, PSRC, PCR, PRHO, PTSTEP) &
               RESULT(PR)
!     ########################################################################
!!
!!****  PPM_S0_X - PPM  advection scheme in X direction in Skamarock 2006 
!!                 notation - NO CONSTRAINTS
!!
!!    MODIFICATIONS
!!    -------------
!!    
!!    20.6.2006.  T. Maric - original version
!!      J.Escobar 21/03/2013: for HALOK comment all NHALO=1 test
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODE_IO_ll
USE MODI_SHUMAN
USE MODI_GET_HALO
!
USE MODD_CONF
!BEG JUAN PPM_LL
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!END JUAN PPM_LL
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER:: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER:: IIE,IJE    ! End useful area in x,y,z directions
!
! advection fluxes
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZFPOS, ZFNEG
!
! variable at cell edges
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZPHAT
!
!BEG JUAN PPM_LL
TYPE(HALO2LIST_ll), POINTER      :: TZ_PSRC_HALO2_ll         ! halo2 for PSRC
INTEGER                          :: IJS,IJN
!END JUAN PPM_LL
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IJS=IJB
IJN=IJE
!!$IJS=IJB-1
!!$IJN=IJE+1
!
!BEG JUAN PPM_LL
!
!*              initialise & update halo & halo2 for PSRC
!
CALL GET_HALO2(PSRC,TZ_PSRC_HALO2_ll)
ZPHAT=PSRC
ZFPOS=PSRC
ZFNEG=PSRC
PR=PSRC
!!$!
!
!END JUAN PPM_LL
!-------------------------------------------------------------------------------
!
! calculate 4th order fluxes at cell edges  
!
!BEG JUAN PPM_LL
!
!   ZPATH(i) = Fct[ PSRC(i),PSRC(i-1),PSRC(i+1),PSRC(i-2)]
!
!       inner domain       
!
ZPHAT(IIB+1:IIE,IJS:IJN,:) = ( 7.0 * &
                       ( PSRC(IIB+1:IIE  ,IJS:IJN,:) + PSRC(IIB  :IIE-1,IJS:IJN,:) ) - &
                       ( PSRC(IIB+2:IIE+1,IJS:IJN,:) + PSRC(IIB-1:IIE-2,IJS:IJN,:) ) ) / 12.0
!
SELECT CASE ( HLBCX(1) ) ! X direction LBC type: (1) for left side
CASE ('CYCL','WALL')            ! In that case one must have HLBCX(1) == HLBCX(2)
!
!!$   ZPHAT(IIB,:,:) = (7.0 * &
!!$                    (PSRC(IIB,:,:) + PSRC(IIB-1,:,:)) - &
!!$                    (PSRC(IIB+1,:,:) + PSRC(IIE-1,:,:))) / 12.0
!!$
!!$!
!!$   ZPHAT(IIE+1,:,:) = ZPHAT(IIB,:,:)
!!$   ZPHAT(IIB-1,:,:) = ZPHAT(IIE,:,:)
!
!  WEST BOUND 
!
   ZPHAT(IIB  ,IJS:IJN,:) = ( 7.0 * &
                      ( PSRC(IIB  ,IJS:IJN,:) + PSRC(IIB-1,IJS:IJN,:)                  ) - &
                      ( PSRC(IIB+1,IJS:IJN,:) + TZ_PSRC_HALO2_ll%HALO2%WEST(IJS:IJN,:) ) ) / 12.0
! <=>  WEST BOUND     ( PSRC(IIB+1,IJS:IJN,:) + PSRC(IIB-2,IJS:IJN,:)                  ) ) / 12.0
!
!  The ZPHAT(IIB-1,:,:) doesn't matter only define an realistic value
!
!!$   ZPHAT(IIB-1,:,:) = ZPHAT(IIB,:,:) ! JUANTEST1
!
! EAST BOUND
!
!  The ZPHAT(IIE+1,:,:) doesn't matter only define an realistic value
!
!!$   ZPHAT(IIE+1,:,:) = ZPHAT(IIE,:,:) ! JUANTEST1
!
!
!   update ZPHAT HALO before next/further  utilisation 
!
CALL  GET_HALO(ZPHAT)
!
   ZFPOS(IIB:IIE+1,IJS:IJN,:) = ZPHAT(IIB:IIE+1,IJS:IJN,:) - & 
        PCR(IIB:IIE+1,IJS:IJN,:)*(ZPHAT(IIB:IIE+1,IJS:IJN,:) - PSRC(IIB-1:IIE,IJS:IJN,:)) - &
        PCR(IIB:IIE+1,IJS:IJN,:)*(1.0 - PCR(IIB:IIE+1,IJS:IJN,:)) * &
        (ZPHAT(IIB-1:IIE,IJS:IJN,:) - 2.0*PSRC(IIB-1:IIE,IJS:IJN,:) + ZPHAT(IIB:IIE+1,IJS:IJN,:))
!
!!$   ZFPOS(IIB-1,:,:) = ZFPOS(IIE,:,:) !JUAN
CALL GET_HALO(ZFPOS) ! JUAN
!
   ZFNEG(IIB-1:IIE,IJS:IJN,:) = ZPHAT(IIB-1:IIE,IJS:IJN,:) + & 
        PCR(IIB-1:IIE,IJS:IJN,:)*(ZPHAT(IIB-1:IIE,IJS:IJN,:) - PSRC(IIB-1:IIE,IJS:IJN,:)) + &
        PCR(IIB-1:IIE,IJS:IJN,:)*(1.0 + PCR(IIB-1:IIE,IJS:IJN,:)) * &
        (ZPHAT(IIB-1:IIE,IJS:IJN,:) - 2.0*PSRC(IIB-1:IIE,IJS:IJN,:) + ZPHAT(IIB:IIE+1,IJS:IJN,:))
!
! define fluxes for CYCL BC outside physical domain
!!$   ZFNEG(IIE+1,:,:) = ZFNEG(IIB,:,:) !JUAN
CALL GET_HALO(ZFNEG) ! JUAN

!
! calculate the advection
!
   PR = PSRC * PRHO - &
        DXF( PCR*MXM(PRHO)*( ZFPOS*(0.5+SIGN(0.5,PCR)) + & 
                             ZFNEG*(0.5-SIGN(0.5,PCR)) ) )
   CALL GET_HALO(PR) ! JUAN
!
CASE ('OPEN')
!
!!$   ZPHAT(IIB,:,:) = 0.5*(PSRC(IIB-1,:,:) + PSRC(IIB,:,:))
!!$   ZPHAT(IIB-1,:,:) = ZPHAT(IIB,:,:)   ! not used
!!$   ZPHAT(IIE+1,:,:) = 0.5*(PSRC(IIE,:,:) + PSRC(IIE+1,:,:))
!
!  WEST BOUND 
!
  IF (.NOT. LWEST_ll()) THEN
   ZPHAT(IIB  ,IJS:IJN,:) = ( 7.0 * &
                      ( PSRC(IIB  ,IJS:IJN,:) + PSRC(IIB-1,IJS:IJN,:)                  ) - &
                      ( PSRC(IIB+1,IJS:IJN,:) + TZ_PSRC_HALO2_ll%HALO2%WEST(IJS:IJN,:) ) ) / 12.0
! <=>  WEST BOUND     ( PSRC(IIB+1,IJS:IJN,:) + PSRC(IIB-2,IJS:IJN,:)                  ) ) / 12.0
  ENDIF
!
CALL  GET_HALO(ZPHAT)
!
  IF (LWEST_ll()) THEN
   ZPHAT(IIB  ,IJS:IJN,:) = 0.5*(PSRC(IIB-1,IJS:IJN,:) + PSRC(IIB,IJS:IJN,:))
   ZPHAT(IIB-1,IJS:IJN,:) = ZPHAT(IIB,IJS:IJN,:)
  ENDIF
!
! EAST BOUND
!
  IF (LEAST_ll()) THEN
   ZPHAT(IIE+1,IJS:IJN,:) = 0.5*(PSRC(IIE,IJS:IJN,:) + PSRC(IIE+1,IJS:IJN,:))
  ENDIF
!
!   update ZPHAT HALO before next/further  utilisation 
!
!!$CALL  GET_HALO(ZPHAT)
!
!!$   ZFPOS(IIB+1:IIE+1,:,:) = ZPHAT(IIB+1:IIE+1,:,:) - & 
!!$        PCR(IIB+1:IIE+1,:,:)*(ZPHAT(IIB+1:IIE+1,:,:) - PSRC(IIB:IIE,:,:)) - &
!!$        PCR(IIB+1:IIE+1,:,:)*(1.0 - PCR(IIB+1:IIE+1,:,:)) * &
!!$        (ZPHAT(IIB:IIE,:,:) - 2.0*PSRC(IIB:IIE,:,:) + ZPHAT(IIB+1:IIE+1,:,:))
   ZFPOS(IIB:IIE+1,IJS:IJN,:) = ZPHAT(IIB:IIE+1,IJS:IJN,:) - & 
        PCR(IIB:IIE+1,IJS:IJN,:)*(ZPHAT(IIB:IIE+1,IJS:IJN,:) - PSRC(IIB-1:IIE,IJS:IJN,:)) - &
        PCR(IIB:IIE+1,IJS:IJN,:)*(1.0 - PCR(IIB:IIE+1,IJS:IJN,:)) * &
        (ZPHAT(IIB-1:IIE,IJS:IJN,:) - 2.0*PSRC(IIB-1:IIE,IJS:IJN,:) + ZPHAT(IIB:IIE+1,IJS:IJN,:))
!
CALL GET_HALO(ZFPOS) ! JUAN
!
! positive flux on the WEST boundary
  IF (LWEST_ll()) THEN
   ZFPOS(IIB,IJS:IJN,:) = (PSRC(IIB-1,IJS:IJN,:) - ZPHAT(IIB,IJS:IJN,:))*PCR(IIB,IJS:IJN,:) + &
                     ZPHAT(IIB,IJS:IJN,:) 
! this is not used
   ZFPOS(IIB-1,IJS:IJN,:) = 0.0
  ENDIF
!
! negative fluxes
!!$   ZFNEG(IIB:IIE,:,:) = ZPHAT(IIB:IIE,:,:) + & 
!!$        PCR(IIB:IIE,:,:)*(ZPHAT(IIB:IIE,:,:) - PSRC(IIB:IIE,:,:)) + &
!!$        PCR(IIB:IIE,:,:)*(1.0 + PCR(IIB:IIE,:,:)) * &
!!$        (ZPHAT(IIB:IIE,:,:) - 2.0*PSRC(IIB:IIE,:,:) + ZPHAT(IIB+1:IIE+1,:,:))
   ZFNEG(IIB-1:IIE,IJS:IJN,:) = ZPHAT(IIB-1:IIE,IJS:IJN,:) + & 
        PCR(IIB-1:IIE,IJS:IJN,:)*(ZPHAT(IIB-1:IIE,IJS:IJN,:) - PSRC(IIB-1:IIE,IJS:IJN,:)) + &
        PCR(IIB-1:IIE,IJS:IJN,:)*(1.0 + PCR(IIB-1:IIE,IJS:IJN,:)) * &
        (ZPHAT(IIB-1:IIE,IJS:IJN,:) - 2.0*PSRC(IIB-1:IIE,IJS:IJN,:) + ZPHAT(IIB:IIE+1,IJS:IJN,:))
!
   CALL GET_HALO(ZFNEG) ! JUAN
!
  IF (LEAST_ll()) THEN
!
! in OPEN case PCR(IIB-1) is not used, so we also set ZFNEG(IIB-1) = 0
!
   ZFNEG(IIB-1,IJS:IJN,:) = 0.0
!
! modified negative flux on EAST boundary. We use linear function instead of a
! parabola to represent the tracer field, so it simplifies the flux expresion
!
   ZFNEG(IIE+1,IJS:IJN,:) = (ZPHAT(IIE+1,IJS:IJN,:) - PSRC(IIE+1,IJS:IJN,:))*PCR(IIE+1,IJS:IJN,:) + &
                       ZPHAT(IIE+1,IJS:IJN,:)
  ENDIF
!
! calculate the advection
!
   PR = PSRC * PRHO - &
        DXF( PCR*MXM(PRHO)*( ZFPOS*(0.5+SIGN(0.5,PCR)) + & 
                             ZFNEG*(0.5-SIGN(0.5,PCR)) ) )
!
! in OPEN case fix boundary conditions
!
  IF (LWEST_ll()) THEN
   WHERE ( PCR(IIB,IJS:IJN,:) <= 0. ) !  OUTFLOW condition
      PR(IIB-1,IJS:IJN,:) = 2.*PR(IIB,IJS:IJN,:) - PR(IIB+1,IJS:IJN,:)
   ELSEWHERE
      PR(IIB-1,IJS:IJN,:) = PR(IIB,IJS:IJN,:)
   END WHERE
  ENDIF
!
  IF (LEAST_ll()) THEN 
   WHERE ( PCR(IIE,IJS:IJN,:) >= 0. ) !  OUTFLOW condition
      PR(IIE+1,IJS:IJN,:) = 2.*PR(IIE,IJS:IJN,:) - PR(IIE-1,IJS:IJN,:)
   ELSEWHERE
      PR(IIE+1,IJS:IJN,:) = PR(IIE,IJS:IJN,:)
   END WHERE
  ENDIF
!
!
END SELECT
!
CALL GET_HALO(PR) 
!
!-------------------------------------------------------------------------------
CALL  DEL_HALO2_ll(TZ_PSRC_HALO2_ll)
!
END FUNCTION PPM_S0_X
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION PPM_S0_Y(HLBCY, KGRID, PSRC, PCR, PRHO, PTSTEP) &
               RESULT(PR)
!     ########################################################################
!!
!!****  PPM_S0_Y - PPM  advection scheme in Y direction in Skamarock 2006 
!!                 notation - NO CONSTRAINTS
!!
!!    MODIFICATIONS
!!    -------------
!!    
!!    20.6.2006.  T. Maric - original version
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODE_IO_ll
USE MODI_SHUMAN
USE MODI_GET_HALO
!
USE MODD_CONF
!USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER:: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER:: IIE,IJE    ! End useful area in x,y,z directions
!
! advection fluxes
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZFPOS, ZFNEG
!
! variable at cell edges
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZPHAT
!
!BEG JUAN PPM_LL
TYPE(HALO2LIST_ll), POINTER      :: TZ_PSRC_HALO2_ll         ! halo2 for PSRC
TYPE(HALO2LIST_ll), POINTER      :: TZ_PHAT_HALO2_ll         ! halo2 for ZPHAT
INTEGER                          :: IIW,IIA
!END JUAN PPM_LL
!
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IIW=IIB
IIA=IIE
!!$IIW=IIB-1
!!$IIA=IIE+1
!
!-------------------------------------------------------------------------------
!
IF ( L2D ) THEN
   PR = PSRC*PRHO
   RETURN
END IF
!
CALL GET_HALO2(PSRC,TZ_PSRC_HALO2_ll)
!
! Initialize with relalistic value all work array 
!
ZPHAT=PSRC
ZFPOS=PSRC
ZFNEG=PSRC
PR=PSRC
!
!-------------------------------------------------------------------------------
!
! calculate 4th order fluxes at cell edges in the inner domain
!
ZPHAT(IIW:IIA,IJB+1:IJE,:) = (7.0 * &
                       (PSRC(IIW:IIA,IJB+1:IJE,:) + PSRC(IIW:IIA,IJB:IJE-1,:)) - &
                       (PSRC(IIW:IIA,IJB+2:IJE+1,:) + PSRC(IIW:IIA,IJB-1:IJE-2,:))) / 12.0
!
SELECT CASE ( HLBCY(1) ) ! Y direction LBC type: (1) for left side
CASE ('CYCL','WALL')            ! In that case one must have HLBCY(1) == HLBCY(2)
!
!!$   ZPHAT(:,IJB,:) = (7.0 * &
!!$                    (PSRC(:,IJB,:) + PSRC(:,IJB-1,:)) - &
!!$                    (PSRC(:,IJB+1,:) + PSRC(:,IJE-1,:))) / 12.0
!!$!
!!$   ZPHAT(:,IJE+1,:) = ZPHAT(:,IJB,:)
!!$   ZPHAT(:,IJB-1,:) = ZPHAT(:,IJE,:)
!
!  SOUTH BOUND 
!
   ZPHAT(IIW:IIA,IJB,:) = ( 7.0 * &
                    ( PSRC(IIW:IIA,IJB  ,:) + PSRC(IIW:IIA,IJB-1,:) ) - &
                    ( PSRC(IIW:IIA,IJB+1,:) + TZ_PSRC_HALO2_ll%HALO2%SOUTH(IIW:IIA,:) ) ) / 12.0
! <=> SOUTH B       ( PSRC(IIW:IIA,IJB+1,:) +    PSRC(IIW:IIA,IJB-2,:)                ) ) / 12.0
!
!  The ZPHAT(:,IJB-1,:) doesn't matter only define an realistic value
!
!!$   ZPHAT(:,IJB-1,:) = ZPHAT(:,IJB,:)
!
!  NORTH BOUND
!
!  The ZPHAT(:IJE+1,:) doesn't matter only define an realistic value
!
!!$   ZPHAT(:,IJE+1,:) = ZPHAT(:,IJE,:)
!
!   update ZPHAT HALO before next/further  utilisation 
!
CALL  GET_HALO(ZPHAT)
!
! calculate the fluxes:
!
   ZFPOS(IIW:IIA,IJB:IJE+1,:) = ZPHAT(IIW:IIA,IJB:IJE+1,:) - & 
        PCR(IIW:IIA,IJB:IJE+1,:)*(ZPHAT(IIW:IIA,IJB:IJE+1,:) - PSRC(IIW:IIA,IJB-1:IJE,:)) - &
        PCR(IIW:IIA,IJB:IJE+1,:)*(1.0 - PCR(IIW:IIA,IJB:IJE+1,:)) * &
        (ZPHAT(IIW:IIA,IJB-1:IJE,:) - 2.0*PSRC(IIW:IIA,IJB-1:IJE,:) + ZPHAT(IIW:IIA,IJB:IJE+1,:))
!
!!$   ZFPOS(:,IJB-1,:) = ZFPOS(:,IJE,:)
CALL GET_HALO(ZFPOS) ! JUAN
!
   ZFNEG(IIW:IIA,IJB-1:IJE,:) = ZPHAT(IIW:IIA,IJB-1:IJE,:) + & 
        PCR(IIW:IIA,IJB-1:IJE,:)*(ZPHAT(IIW:IIA,IJB-1:IJE,:) - PSRC(IIW:IIA,IJB-1:IJE,:)) + &
        PCR(IIW:IIA,IJB-1:IJE,:)*(1.0 + PCR(IIW:IIA,IJB-1:IJE,:)) * &
        (ZPHAT(IIW:IIA,IJB-1:IJE,:) - 2.0*PSRC(IIW:IIA,IJB-1:IJE,:) +ZPHAT(IIW:IIA,IJB:IJE+1,:))
!

!
! define fluxes for CYCL BC outside physical domain
!!$   ZFNEG(:,IJE+1,:) = ZFNEG(:,IJB,:)
CALL GET_HALO(ZFNEG) ! JUAN
!
! calculate the advection
!
   PR = PSRC * PRHO - &
        DYF( PCR*MYM(PRHO)*( ZFPOS*(0.5+SIGN(0.5,PCR)) + & 
                             ZFNEG*(0.5-SIGN(0.5,PCR)) ) )
!
CASE ('OPEN')
!
!!$   ZPHAT(:,IJB,:) = 0.5*(PSRC(:,IJB-1,:) + PSRC(:,IJB,:))
!!$   ZPHAT(:,IJB-1,:) = ZPHAT(:,IJB,:)   ! not used
!!$   ZPHAT(:,IJE+1,:) = 0.5*(PSRC(:,IJE,:) + PSRC(:,IJE+1,:))
!
!
!  SOUTH BOUND 
!
  IF ( .NOT. LSOUTH_ll()) THEN
   ZPHAT(IIW:IIA,IJB  ,:) = (7.0 * &
                      (PSRC(IIW:IIA,IJB  ,:) + PSRC(IIW:IIA,IJB-1,:)) - &
                      (PSRC(IIW:IIA,IJB+1,:) + TZ_PSRC_HALO2_ll%HALO2%SOUTH(IIW:IIA,:) )) / 12.0
! <=> SOUTH BOUND     (PSRC(IIW:IIA,IJB+1,:) + PSRC(IIW:IIA,IJB-2,:)                   )) / 12.0
  ENDIF
!
CALL  GET_HALO(ZPHAT)
!
  IF (LSOUTH_ll()) THEN
   ZPHAT(IIW:IIA,IJB  ,:) = 0.5*(PSRC(IIW:IIA,IJB-1,:) + PSRC(IIW:IIA,IJB,:))
   ZPHAT(IIW:IIA,IJB-1,:) = ZPHAT(IIW:IIA,IJB,:)
  ENDIF
!
! NORTH BOUND
!
  IF (LNORTH_ll()) THEN
   ZPHAT(IIW:IIA,IJE+1,:) =  0.5*(PSRC(IIW:IIA,IJE,:) + PSRC(IIW:IIA,IJE+1,:))
  ENDIF
!
!
!   update ZPHAT HALO before next/further  utilisation 
!
!!$CALL  GET_HALO(ZPHAT)
!
! calculate the fluxes:
! positive fluxes
!!$   ZFPOS(:,IJB+1:IJE+1,:) = ZPHAT(:,IJB+1:IJE+1,:) - & 
!!$        PCR(:,IJB+1:IJE+1,:)*(ZPHAT(:,IJB+1:IJE+1,:) - PSRC(:,IJB:IJE,:)) - &
!!$        PCR(:,IJB+1:IJE+1,:)*(1.0 - PCR(:,IJB+1:IJE+1,:)) * &
!!$        (ZPHAT(:,IJB:IJE,:) - 2.0*PSRC(:,IJB:IJE,:) + ZPHAT(:,IJB+1:IJE+1,:))
ZFPOS(IIW:IIA,IJB:IJE+1,:) = ZPHAT(IIW:IIA,IJB:IJE+1,:) - & 
        PCR(IIW:IIA,IJB:IJE+1,:)*( ZPHAT(IIW:IIA,IJB:IJE+1,:) - PSRC(IIW:IIA,IJB-1:IJE  ,:) ) - &
        PCR(IIW:IIA,IJB:IJE+1,:)*( 1.0                  -  PCR(IIW:IIA,IJB  :IJE+1,:) ) * &
        (ZPHAT(IIW:IIA,IJB-1:IJE,:) - 2.0*PSRC(IIW:IIA,IJB-1:IJE,:) + ZPHAT(IIW:IIA,IJB:IJE+1,:))
!
CALL GET_HALO(ZFPOS) ! JUAN
!
! positive flux on the SOUTH boundary
  IF (LSOUTH_ll()) THEN
   ZFPOS(IIW:IIA,IJB,:) = (PSRC(IIW:IIA,IJB-1,:) - ZPHAT(IIW:IIA,IJB,:))*PCR(IIW:IIA,IJB,:) + &
                     ZPHAT(IIW:IIA,IJB,:)
!
! this is not used
   ZFPOS(IIW:IIA,IJB-1,:) = 0.0
  ENDIF
! 
! negative fluxes
!!$   ZFNEG(:,IJB:IJE,:) = ZPHAT(:,IJB:IJE,:) + & 
!!$        PCR(:,IJB:IJE,:)*(ZPHAT(:,IJB:IJE,:) - PSRC(:,IJB:IJE,:)) + &
!!$        PCR(:,IJB:IJE,:)*(1.0 + PCR(:,IJB:IJE,:)) * &
!!$        (ZPHAT(:,IJB:IJE,:) - 2.0*PSRC(:,IJB:IJE,:) +ZPHAT(:,IJB+1:IJE+1,:))
   ZFNEG(IIW:IIA,IJB-1:IJE,:) = ZPHAT(IIW:IIA,IJB-1:IJE,:) + & 
        PCR(IIW:IIA,IJB-1:IJE,:)*(ZPHAT(IIW:IIA,IJB-1:IJE,:) - PSRC(IIW:IIA,IJB-1:IJE,:)) + &
        PCR(IIW:IIA,IJB-1:IJE,:)*(1.0 + PCR(IIW:IIA,IJB-1:IJE,:)) * &
        (ZPHAT(IIW:IIA,IJB-1:IJE,:) - 2.0*PSRC(IIW:IIA,IJB-1:IJE,:) +ZPHAT(IIW:IIA,IJB:IJE+1,:))
!
   CALL GET_HALO(ZFNEG) ! JUAN
!
  IF (LNORTH_ll()) THEN
! this is not used
   ZFNEG(IIW:IIA,IJB-1,:) = 0.0
!
! negative flux on the NORTH boundary
   ZFNEG(IIW:IIA,IJE+1,:) = (ZPHAT(IIW:IIA,IJE+1,:) - PSRC(IIW:IIA,IJE+1,:))*PCR(IIW:IIA,IJE+1,:) + &
                       ZPHAT(IIW:IIA,IJE+1,:)
  ENDIF
!
! calculate the advection
!
   PR = PSRC * PRHO - &
        DYF( PCR*MYM(PRHO)*( ZFPOS*(0.5+SIGN(0.5,PCR)) + & 
                             ZFNEG*(0.5-SIGN(0.5,PCR)) ) )
!
! in OPEN case fix boundary conditions
!
  IF (LSOUTH_ll()) THEN
   WHERE ( PCR(IIW:IIA,IJB,:) <= 0. ) !  OUTFLOW condition
      PR(IIW:IIA,IJB-1,:) = 1.0 * 2.*PR(IIW:IIA,IJB,:) - PR(IIW:IIA,IJB+1,:)
   ELSEWHERE
      PR(IIW:IIA,IJB-1,:) = PR(IIW:IIA,IJB,:) 
   END WHERE
  ENDIF
!
  IF (LNORTH_ll()) THEN
   WHERE ( PCR(IIW:IIA,IJE,:) >= 0. ) !  OUTFLOW condition
      PR(IIW:IIA,IJE+1,:) = 1.0 * 2.*PR(IIW:IIA,IJE,:) - PR(IIW:IIA,IJE-1,:)
   ELSEWHERE
      PR(IIW:IIA,IJE+1,:) = PR(IIW:IIA,IJE,:) 
   END WHERE
  ENDIF
!
! 
!
END SELECT
!
CALL GET_HALO(PR) 
!
CALL  DEL_HALO2_ll(TZ_PSRC_HALO2_ll)
!
END FUNCTION PPM_S0_Y
!
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION PPM_S0_Z(KGRID, PSRC, PCR, PRHO, PTSTEP) &
               RESULT(PR)
!     ########################################################################
!!
!!****  PPM_S0_Z - PPM  advection scheme in Z direction in Skamarock 2006 
!!                 notation - NO CONSTRAINTS
!!
!!    MODIFICATIONS
!!    -------------
!!    
!!    20.6.2006.  T. Maric - original version
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODI_SHUMAN
USE MODI_GET_HALO
!
USE MODD_CONF
USE MODD_PARAMETERS
!USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER:: IKB    ! Begining useful area in x,y,z directions
INTEGER:: IKE    ! End useful area in x,y,z directions
INTEGER:: IKU
!
! advection fluxes
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZFPOS, ZFNEG
!
! interpolated variable at cell edges
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZPHAT
!
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
IKB = 1 + JPVEXT
IKE = SIZE(PSRC,3) - JPVEXT
IKU = SIZE(PSRC,3)
!
!-------------------------------------------------------------------------------
!
! calculate 4th order fluxes at cell edges in the inner domain
!
CALL GET_HALO(PSRC)
!
!
ZPHAT(:,:,IKB+1:IKE) = (7.0 * &
                       (PSRC(:,:,IKB+1:IKE) + PSRC(:,:,IKB:IKE-1)) - &
                       (PSRC(:,:,IKB+2:IKE+1) + PSRC(:,:,IKB-1:IKE-2))) / 12.0
!
! set OPEN BC at the top and bottom
ZPHAT(:,:,IKB) = 0.5*(PSRC(:,:,IKB-1) + PSRC(:,:,IKB))
ZPHAT(:,:,IKB-1) = ZPHAT(:,:,IKB)  ! not used
ZPHAT(:,:,IKE+1) = 0.5*(PSRC(:,:,IKE) + PSRC(:,:,IKE+1))
!
!!$CALL  GET_HALO(ZPHAT)
!
! calculate fluxes through cell edges for positive and negative Courant numbers
! (for inflow or outflow situation)
!
ZFPOS(:,:,IKB+1:IKE+1) = ZPHAT(:,:,IKB+1:IKE+1) - & 
     PCR(:,:,IKB+1:IKE+1)*(ZPHAT(:,:,IKB+1:IKE+1) - PSRC(:,:,IKB:IKE)) - &
     PCR(:,:,IKB+1:IKE+1)*(1.0 - PCR(:,:,IKB+1:IKE+1)) * &
     (ZPHAT(:,:,IKB:IKE) - 2.0*PSRC(:,:,IKB:IKE) + ZPHAT(:,:,IKB+1:IKE+1))
!
!!$CALL GET_HALO(ZFPOS) ! JUAN
!
! positive flux on the BOTTOM boundary
ZFPOS(:,:,IKB) = (PSRC(:,:,IKB-1) - ZPHAT(:,:,IKB))*PCR(:,:,IKB) + &
                  ZPHAT(:,:,IKB)
!
! below bottom flux - not used
ZFPOS(:,:,IKB-1) = 0.0
!
! negative fluxes:
!
ZFNEG(:,:,IKB:IKE) = ZPHAT(:,:,IKB:IKE) + & 
     PCR(:,:,IKB:IKE)*(ZPHAT(:,:,IKB:IKE) - PSRC(:,:,IKB:IKE)) + &
     PCR(:,:,IKB:IKE)*(1.0 + PCR(:,:,IKB:IKE)) * &
     (ZPHAT(:,:,IKB:IKE) - 2.0*PSRC(:,:,IKB:IKE) +ZPHAT(:,:,IKB+1:IKE+1))
!
!!$   CALL GET_HALO(ZFNEG) ! JUAN
!
! set bottom negative flux to 0
ZFNEG(:,:,IKB-1) = 0.0 
!
! negative flux at the TOP
ZFNEG(:,:,IKE+1) = (ZPHAT(:,:,IKE+1) - PSRC(:,:,IKE+1))*PCR(:,:,IKE+1) + &
                    ZPHAT(:,:,IKE+1) 
!
! calculate the advection
!
PR = PSRC * PRHO - &
     DZF(1,IKU,1, PCR*MZM(1,IKU,1,PRHO)*( ZFPOS*(0.5+SIGN(0.5,PCR)) + & 
                          ZFNEG*(0.5-SIGN(0.5,PCR)) ) )
!
! in OPEN case fix boundary conditions
!
      PR(:,:,IKB-1) = PR(:,:,IKB)
      PR(:,:,IKE+1) = PR(:,:,IKE)
!
   CALL GET_HALO(PR) ! JUAN
!
END FUNCTION PPM_S0_Z
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION PPM_S1_X(HLBCX, KGRID, PSRC, PCR, PRHO, PRHOT, &
                        PTSTEP) RESULT(PR)
!     ########################################################################
!!
!!****  PPM_S1_X - PPM  advection scheme in X direction in Skamarock 2006 
!!                 notation - with flux limiting for monotonicity
!!
!!    MODIFICATIONS
!!    -------------
!!    
!!    23.6.2006.  T. Maric - original version
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODE_IO_ll
USE MODI_SHUMAN
!
USE MODD_CONF
USE MODD_PARAMETERS
!USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHOT ! density at t+dt
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER:: IIB,IJB,IKB    ! Begining useful area in x,y,z directions
INTEGER:: IIE,IJE,IKE    ! End useful area in x,y,z directions
!
! variable at cell edges
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZPHAT, ZRUT
!
! advection fluxes, upwind and correction
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZFUP, ZFCOR
!
! ratios for limiting the correction flux
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZRPOS, ZRNEG
!
! variables for limiting the correction flux
REAL :: ZSRCMAX, ZSRCMIN, ZFIN, ZFOUT
!
REAL, PARAMETER :: ZEPS = 1.0E-16
!
INTEGER :: II, IJ, IK
INTEGER                          :: IRESP             ! for prints
!
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PSRC,3) - JPVEXT
!
!-------------------------------------------------------------------------------
!
! Calculate contravariant component rho*u/dx
!
ZRUT = PCR/PTSTEP * MXM(PRHO)
!
! calculate 4th order fluxes at cell edges in the inner domain
!
ZPHAT(IIB+1:IIE,:,:) = (7.0 * &
                       (PSRC(IIB+1:IIE,:,:) + PSRC(IIB:IIE-1,:,:)) - &
                       (PSRC(IIB+2:IIE+1,:,:) + PSRC(IIB-1:IIE-2,:,:))) / 12.0
!
SELECT CASE ( HLBCX(1) ) ! X direction LBC type: (1) for left side
CASE ('CYCL','WALL')            ! In that case one must have HLBCX(1) == HLBCX(2)
!
   ZPHAT(IIB,:,:) = (7.0 * &
                    (PSRC(IIB,:,:) + PSRC(IIB-1,:,:)) - &
                    (PSRC(IIB+1,:,:) + PSRC(IIE-1,:,:))) / 12.0
!
   ZPHAT(IIE+1,:,:) = ZPHAT(IIB,:,:)
   ZPHAT(IIB-1,:,:) = ZPHAT(IIE,:,:)
!
CASE ('OPEN')
!
   ZPHAT(IIB,:,:) = 0.5*(PSRC(IIB-1,:,:) + PSRC(IIB,:,:))
   ZPHAT(IIB-1,:,:) = ZPHAT(IIB,:,:)
   ZPHAT(IIE+1,:,:) = 0.5*(PSRC(IIE,:,:) + PSRC(IIE+1,:,:))
!
!
END SELECT
!
! calculate upwind and correction fluxes. upwind flux is upstream value of the
! scalar variable, and correction flux is the correction to the upstream flux
! that makes it equivalent to the PPM flux
! flux_ppm = flux_up + flux_corr
!
WHERE ( PCR(IIB:IIE,:,:) .GT. 0.0 )
   ZFUP(IIB:IIE,:,:) = ZRUT(IIB:IIE,:,:) * PSRC(IIB-1:IIE-1,:,:)
   ZFCOR(IIB:IIE,:,:) = ZRUT(IIB:IIE,:,:) * &
        (1.0 - PCR(IIB:IIE,:,:)) * &
        (ZPHAT(IIB:IIE,:,:) - PSRC(IIB-1:IIE-1,:,:) - PCR(IIB:IIE,:,:) * &
        (ZPHAT(IIB-1:IIE-1,:,:) - 2.0*PSRC(IIB-1:IIE-1,:,:)+ZPHAT(IIB:IIE,:,:)))
ELSEWHERE
   ZFUP(IIB:IIE,:,:) = ZRUT(IIB:IIE,:,:) * PSRC(IIB:IIE,:,:)
   ZFCOR(IIB:IIE,:,:) = ZRUT(IIB:IIE,:,:) * &
        (1.0 + PCR(IIB:IIE,:,:)) * &
        (ZPHAT(IIB:IIE,:,:) - PSRC(IIB:IIE,:,:) + PCR(IIB:IIE,:,:) * &
        (ZPHAT(IIB:IIE,:,:) - 2.0*PSRC(IIB:IIE,:,:) + ZPHAT(IIB+1:IIE+1,:,:)))
END WHERE
!
! set boundaries to CYCL
!
WHERE ( PCR(IIB-1,:,:) .GT. 0.0 )
   ZFUP(IIB-1,:,:) = ZRUT(IIB-1,:,:) * PSRC(IIE-1,:,:)
   ZFCOR(IIB-1,:,:) =  ZRUT(IIB-1,:,:) * &
        (1.0 - PCR(IIB-1,:,:)) * &
        (ZPHAT(IIB-1,:,:) - PSRC(IIE-1,:,:) - PCR(IIB-1,:,:) * &
        (ZPHAT(IIE-1,:,:) - 2.0*PSRC(IIE-1,:,:) + ZPHAT(IIB-1,:,:)))
ELSEWHERE
   ZFUP(IIB-1,:,:) = ZRUT(IIB-1,:,:) * PSRC(IIB-1,:,:)
   ZFCOR(IIB-1,:,:) =  ZRUT(IIB-1,:,:) * &
        (1.0 + PCR(IIB-1,:,:)) * &
        (ZPHAT(IIB-1,:,:) - PSRC(IIB-1,:,:) + PCR(IIB-1,:,:) * &
        (ZPHAT(IIB-1,:,:) - 2.0*PSRC(IIB-1,:,:) + ZPHAT(IIB,:,:)))
END WHERE
!
WHERE ( PCR(IIE+1,:,:) .GT. 0.0 )
   ZFUP(IIE+1,:,:) = ZRUT(IIE+1,:,:) * PSRC(IIE,:,:)
   ZFCOR(IIE+1,:,:) =  ZRUT(IIE+1,:,:) * &
        (1.0 - PCR(IIE+1,:,:)) * &
        (ZPHAT(IIE+1,:,:) - PSRC(IIE,:,:) - PCR(IIE+1,:,:) * &
        (ZPHAT(IIE,:,:) - 2.0*PSRC(IIE,:,:) + ZPHAT(IIE+1,:,:)))
ELSEWHERE
   ZFUP(IIE+1,:,:) = ZRUT(IIE+1,:,:) * PSRC(IIE+1,:,:)
   ZFCOR(IIE+1,:,:) =  ZRUT(IIE+1,:,:) * &
        (1.0 + PCR(IIE+1,:,:)) * &
        (ZPHAT(IIE+1,:,:) - PSRC(IIE+1,:,:) + PCR(IIE+1,:,:) * &
        (ZPHAT(IIE+1,:,:) - 2.0*PSRC(IIE+1,:,:) + ZPHAT(IIB+1,:,:)))
END WHERE
!
! Perform limiting of the fluxes
!
! 1. calculate upwind tendency of the source
!
PR = PSRC*PRHO - PTSTEP*DXF(ZFUP)
!
!-------------------------------------------------------------------------------
! compute and apply the limiters
!
DO II = IIB,IIE
   DO IJ = IJB-1,IJE+1
      DO IK = IKB-1,IKE+1         
!
! 2. find local extrema in the source 
!
         ZSRCMAX = MAX( PSRC(II-1,IJ,IK),PSRC(II,IJ,IK),PSRC(II+1,IJ,IK) )
         ZSRCMIN = MIN( PSRC(II-1,IJ,IK),PSRC(II,IJ,IK),PSRC(II+1,IJ,IK) )
!
! 3. compute incoming and outgoing fluxes for this cell
!
         ZFOUT = MAX(ZEPS,MAX(0.,ZFCOR(II+1,IJ,IK)) - MIN(0.,ZFCOR(II,IJ,IK)))
         ZFIN  = MAX(ZEPS,MAX(0.,ZFCOR(II,IJ,IK)) - MIN(0.,ZFCOR(II+1,IJ,IK)))
!
! 4. calculate fraction of outgoing and incoming flux which will drive scalar
!    values outside the local extrema
!
         ZRNEG(II,IJ,IK) = MAX(0.,MIN(1., &
              (PR(II,IJ,IK) - PRHOT(II,IJ,IK)*ZSRCMIN) &
              / PTSTEP / ZFOUT))
!
         ZRPOS(II,IJ,IK) = MAX(0.,MIN(1., &
              (PRHOT(II,IJ,IK)*ZSRCMAX - PR(II,IJ,IK)) &
              / PTSTEP / ZFIN))
      END DO
   END DO
END DO
!
! set CYCL boundaries
!
DO IJ = IJB-1,IJE+1
   DO IK = IKB-1,IKE+1         
!
      ZSRCMAX = MAX( PSRC(IIE-1,IJ,IK),PSRC(IIB-1,IJ,IK),PSRC(IIB,IJ,IK) )
      ZSRCMIN = MIN( PSRC(IIE-1,IJ,IK),PSRC(IIB-1,IJ,IK),PSRC(IIB,IJ,IK) )
!
      ZFOUT = MAX(ZEPS,MAX(0.,ZFCOR(IIB,IJ,IK)) - MIN(0.,ZFCOR(IIB-1,IJ,IK)))
      ZFIN  = MAX(ZEPS,MAX(0.,ZFCOR(IIB-1,IJ,IK)) - MIN(0.,ZFCOR(IIB,IJ,IK)))
!
      ZRNEG(IIB-1,IJ,IK) = MAX(0.,MIN(1., &
           (PR(IIB-1,IJ,IK) - PRHOT(IIB-1,IJ,IK)*ZSRCMIN) &
           / PTSTEP / ZFOUT))
!
      ZRPOS(IIB-1,IJ,IK) = MAX(0.,MIN(1., &
           (PRHOT(IIB-1,IJ,IK)*ZSRCMAX - PR(IIB-1,IJ,IK)) &
           / PTSTEP / ZFIN))
!
! 
      ZSRCMAX = MAX( PSRC(IIE,IJ,IK),PSRC(IIE+1,IJ,IK),PSRC(IIB+1,IJ,IK) )
      ZSRCMIN = MIN( PSRC(IIE,IJ,IK),PSRC(IIE+1,IJ,IK),PSRC(IIB+1,IJ,IK) )
!
      ZFOUT = MAX(ZEPS,MAX(0.,ZFCOR(IIB+1,IJ,IK)) - MIN(0.,ZFCOR(IIE+1,IJ,IK)))
      ZFIN  = MAX(ZEPS,MAX(0.,ZFCOR(IIE+1,IJ,IK)) - MIN(0.,ZFCOR(IIB+1,IJ,IK)))
!
      ZRNEG(IIE+1,IJ,IK) = MAX(0.,MIN(1., &
           (PR(IIE+1,IJ,IK) - PRHOT(IIE+1,IJ,IK)*ZSRCMIN) &
           / PTSTEP / ZFOUT))
!
      ZRPOS(IIE+1,IJ,IK) = MAX(0.,MIN(1., &
           (PRHOT(IIE+1,IJ,IK)*ZSRCMAX - PR(IIE+1,IJ,IK)) &
           / PTSTEP / ZFIN))
!
   END DO
END DO
!
! 5. apply the limit to the fluxes where needed
!
ZFCOR(IIB:IIE+1,:,:) = MAX( &
     MIN(ZRNEG(IIB:IIE+1,:,:),ZRPOS(IIB-1:IIE,:,:)) * ZFCOR(IIB:IIE+1,:,:), &
     ZFCOR(IIB:IIE+1,:,:) )
ZFCOR(IIB-1,:,:) = MAX( &
     MIN(ZRNEG(IIB-1,:,:),ZRPOS(IIE-1,:,:))*ZFCOR(IIB-1,:,:),ZFCOR(IIB-1,:,:))
!ZFCOR(IIB-1,:,:) = MAX( ZRNEG(IIB-1,:,:)*ZFCOR(IIB-1,:,:), ZFCOR(IIB-1,:,:) )
!
ZFCOR(IIB:IIE+1,:,:) = MIN( &
     MIN(ZRPOS(IIB:IIE+1,:,:),ZRNEG(IIB-1:IIE,:,:)) * ZFCOR(IIB:IIE+1,:,:), &
     ZFCOR(IIB:IIE+1,:,:) )
ZFCOR(IIB-1,:,:) = MIN( &
     MIN(ZRPOS(IIB-1,:,:),ZRNEG(IIE-1,:,:))*ZFCOR(IIB-1,:,:),ZFCOR(IIB-1,:,:))
!ZFCOR(IIB-1,:,:) = MIN( ZRPOS(IIB-1,:,:)*ZFCOR(IIB-1,:,:), ZFCOR(IIB-1,:,:) )

!-------------------------------------------------------------------------------
! 6. apply the limited flux correction to scalar field
!
PR = PR - PTSTEP*DXF(ZFCOR)
!
!
END FUNCTION PPM_S1_X
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION PPM_S1_Y(HLBCY, KGRID, PSRC, PCR, PRHO, PRHOT, &
                        PTSTEP) RESULT(PR)
!     ########################################################################
!!
!!****  PPM_S1_Y - PPM  advection scheme in Y direction in Skamarock 2006 
!!                 notation - with flux limiting for monotonicity
!!
!!    MODIFICATIONS
!!    -------------
!!    
!!    23.6.2006.  T. Maric - original version
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODE_IO_ll
USE MODI_SHUMAN
!
USE MODD_CONF
USE MODD_PARAMETERS
!USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY  ! X direction LBC type
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHOT ! density at t+dt
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER:: IIB,IJB,IKB   ! Begining useful area in x,y,z directions
INTEGER:: IIE,IJE,IKE   ! End useful area in x,y,z directions
!
! variable at cell edges
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZPHAT, ZRVT
!
! advection fluxes, upwind and correction
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZFUP, ZFCOR
!
! ratios for limiting the correction flux
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZRPOS, ZRNEG
!
! variables for limiting the correction flux
REAL :: ZSRCMAX, ZSRCMIN, ZFIN, ZFOUT
!
!
REAL, PARAMETER :: ZEPS = 1.0E-16
!
INTEGER :: II, IJ, IK
INTEGER  :: IRESP   ! Return code of FM-routines
!
!-------------------------------------------------------------------------------
!
IF ( L2D ) THEN
   PR = PSRC*PRHO
   RETURN
END IF
!
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PSRC,3) - JPVEXT
!
!-------------------------------------------------------------------------------
!
ZRVT = PCR/PTSTEP * MYM(PRHO)
!
! calculate 4th order fluxes at cell edges in the inner domain !
ZPHAT(:,IJB+1:IJE,:) = (7.0 * &
                       (PSRC(:,IJB+1:IJE,:) + PSRC(:,IJB:IJE-1,:)) - &
                       (PSRC(:,IJB+2:IJE+1,:) + PSRC(:,IJB-1:IJE-2,:))) / 12.0
!
SELECT CASE ( HLBCY(1) ) ! X direction LBC type: (1) for left side
CASE ('CYCL','WALL')            ! In that case one must have HLBCY(1) == HLBCY(2)
!
   ZPHAT(:,IJB,:) = (7.0 * &
                    (PSRC(:,IJB,:) + PSRC(:,IJB-1,:)) - &
                    (PSRC(:,IJB+1,:) + PSRC(:,IJE-1,:))) / 12.0
!
   ZPHAT(:,IJE+1,:) = ZPHAT(:,IJB,:)
   ZPHAT(:,IJB-1,:) = ZPHAT(:,IJE,:)
!
CASE ('OPEN')
!
   ZPHAT(:,IJB,:) = 0.5*(PSRC(:,IJB-1,:) + PSRC(:,IJB,:))
   ZPHAT(:,IJB-1,:) = ZPHAT(:,IJB,:)
   ZPHAT(:,IJE+1,:) = 0.5*(PSRC(:,IJE,:) + PSRC(:,IJE+1,:))
!
!
END SELECT
!
! calculate upwind and correction fluxes. upwind flux is upstream value of the
! scalar variable, and correction flux is the correction to the upstream flux
! that makes it equivalent to the PPM flux
! flux_ppm = flux_up + flux_corr
!
WHERE ( PCR(:,IJB:IJE,:) .GT. 0.0 )
   ZFUP(:,IJB:IJE,:) = ZRVT(:,IJB:IJE,:) * PSRC(:,IJB-1:IJE-1,:)
   ZFCOR(:,IJB:IJE,:) = ZRVT(:,IJB:IJE,:) * &
        (1.0 - PCR(:,IJB:IJE,:)) * &
        (ZPHAT(:,IJB:IJE,:) - PSRC(:,IJB-1:IJE-1,:) - PCR(:,IJB:IJE,:) * &
        (ZPHAT(:,IJB-1:IJE-1,:) - 2.0*PSRC(:,IJB-1:IJE-1,:)+ZPHAT(:,IJB:IJE,:)))
ELSEWHERE
   ZFUP(:,IJB:IJE,:) = ZRVT(:,IJB:IJE,:) * PSRC(:,IJB:IJE,:)
   ZFCOR(:,IJB:IJE,:) = ZRVT(:,IJB:IJE,:) * &
        (1.0 + PCR(:,IJB:IJE,:)) * &
        (ZPHAT(:,IJB:IJE,:) - PSRC(:,IJB:IJE,:) + PCR(:,IJB:IJE,:) * &
        (ZPHAT(:,IJB:IJE,:) - 2.0*PSRC(:,IJB:IJE,:) + ZPHAT(:,IJB+1:IJE+1,:)))
END WHERE
!
! set boundaries to CYCL
!
WHERE ( PCR(:,IJB-1,:) .GT. 0.0 )
   ZFUP(:,IJB-1,:) = ZRVT(:,IJB-1,:) * PSRC(:,IJE-1,:)
   ZFCOR(:,IJB-1,:) =  ZRVT(:,IJB-1,:) * &
        (1.0 - PCR(:,IJB-1,:)) * &
        (ZPHAT(:,IJB-1,:) - PSRC(:,IJE-1,:) - PCR(:,IJB-1,:) * &
        (ZPHAT(:,IJE-1,:) - 2.0*PSRC(:,IJE-1,:) + ZPHAT(:,IJB-1,:)))
ELSEWHERE
   ZFUP(:,IJB-1,:) = ZRVT(:,IJB-1,:) * PSRC(:,IJB-1,:)
   ZFCOR(:,IJB-1,:) =  ZRVT(:,IJB-1,:) * &
        (1.0 + PCR(:,IJB-1,:)) * &
        (ZPHAT(:,IJB-1,:) - PSRC(:,IJB-1,:) + PCR(:,IJB-1,:) * &
        (ZPHAT(:,IJB-1,:) - 2.0*PSRC(:,IJB-1,:) + ZPHAT(:,IJB,:)))
END WHERE
!
WHERE ( PCR(:,IJE+1,:) .GT. 0.0 )
   ZFUP(:,IJE+1,:) = ZRVT(:,IJE+1,:) * PSRC(:,IJE,:)
   ZFCOR(:,IJE+1,:) =  ZRVT(:,IJE+1,:) * &
        (1.0 - PCR(:,IJE+1,:)) * &
        (ZPHAT(:,IJE+1,:) - PSRC(:,IJE,:) - PCR(:,IJE+1,:) * &
        (ZPHAT(:,IJE,:) - 2.0*PSRC(:,IJE,:) + ZPHAT(:,IJE+1,:)))
ELSEWHERE
   ZFUP(:,IJE+1,:) = ZRVT(:,IJE+1,:) * PSRC(:,IJE+1,:)
   ZFCOR(:,IJE+1,:) =  ZRVT(:,IJE+1,:) * &
        (1.0 + PCR(:,IJE+1,:)) * &
        (ZPHAT(:,IJE+1,:) - PSRC(:,IJE+1,:) + PCR(:,IJE+1,:) * &
        (ZPHAT(:,IJE+1,:) - 2.0*PSRC(:,IJE+1,:) + ZPHAT(:,IJB+1,:)))
END WHERE
!
! Perform limiting of the fluxes
!
! 1. calculate upwind tendency of the source
!
PR = PSRC*PRHO - PTSTEP*DYF(ZFUP)
!
!-------------------------------------------------------------------------------
! compute and apply the limiters
!
DO II = IIB-1,IIE+1
   DO IJ = IJB,IJE
      DO IK = IKB-1,IKE+1         
!
! 2. find local extrema in the source 
!
         ZSRCMAX = MAX( PSRC(II,IJ-1,IK),PSRC(II,IJ,IK),PSRC(II,IJ+1,IK) )
         ZSRCMIN = MIN( PSRC(II,IJ-1,IK),PSRC(II,IJ,IK),PSRC(II,IJ+1,IK) )
!
! 3. compute incoming and outgoing fluxes for this cell
!
         ZFOUT = MAX(ZEPS,MAX(0.,ZFCOR(II,IJ+1,IK)) - MIN(0.,ZFCOR(II,IJ,IK)))
         ZFIN  = MAX(ZEPS,MAX(0.,ZFCOR(II,IJ,IK)) - MIN(0.,ZFCOR(II,IJ+1,IK)))
!
! 4. calculate fraction of outgoing and incoming flux which will drive scalar
!    values outside the local extrema
!
         ZRNEG(II,IJ,IK) = MAX(0.,MIN(1., &
              (PR(II,IJ,IK) - PRHOT(II,IJ,IK)*ZSRCMIN) &
              / PTSTEP / ZFOUT))
!
         ZRPOS(II,IJ,IK) = MAX(0.,MIN(1., &
              (PRHOT(II,IJ,IK)*ZSRCMAX - PR(II,IJ,IK)) &
              / PTSTEP / ZFIN))
      END DO
   END DO
END DO
!
! set CYCL boundaries
!
DO II = IIB-1,IIE+1
   DO IK = IKB-1,IKE+1         
!
      ZSRCMAX = MAX( PSRC(II,IJE-1,IK),PSRC(II,IJB-1,IK),PSRC(II,IJB,IK) )
      ZSRCMIN = MIN( PSRC(II,IJE-1,IK),PSRC(II,IJB-1,IK),PSRC(II,IJB,IK) )
!
      ZFOUT = MAX(ZEPS,MAX(0.,ZFCOR(II,IJB,IK)) - MIN(0.,ZFCOR(II,IJB-1,IK)))
      ZFIN  = MAX(ZEPS,MAX(0.,ZFCOR(II,IJB-1,IK)) - MIN(0.,ZFCOR(II,IJB,IK)))
!
      ZRNEG(II,IJB-1,IK) = MAX(0.,MIN(1., &
           (PR(II,IJB-1,IK) - PRHOT(II,IJB-1,IK)*ZSRCMIN) &
           / PTSTEP / ZFOUT))
!
      ZRPOS(II,IJB-1,IK) = MAX(0.,MIN(1., &
           (PRHOT(II,IJB-1,IK)*ZSRCMAX - PR(II,IJB-1,IK)) &
           / PTSTEP / ZFIN))
!
! 
      ZSRCMAX = MAX( PSRC(II,IJE,IK),PSRC(II,IJE+1,IK),PSRC(II,IJB+1,IK) )
      ZSRCMIN = MIN( PSRC(II,IJE,IK),PSRC(II,IJE+1,IK),PSRC(II,IJB+1,IK) )
!
      ZFOUT = MAX(ZEPS,MAX(0.,ZFCOR(II,IJB+1,IK)) - MIN(0.,ZFCOR(II,IJE+1,IK)))
      ZFIN  = MAX(ZEPS,MAX(0.,ZFCOR(II,IJE+1,IK)) - MIN(0.,ZFCOR(II,IJB+1,IK)))
!
      ZRNEG(II,IJE+1,IK) = MAX(0.,MIN(1., &
           (PR(II,IJE+1,IK) - PRHOT(II,IJE+1,IK)*ZSRCMIN) &
           / PTSTEP / ZFOUT))
!
      ZRPOS(II,IJE+1,IK) = MAX(0.,MIN(1., &
           (PRHOT(II,IJE+1,IK)*ZSRCMAX - PR(II,IJE+1,IK)) &
           / PTSTEP / ZFIN))
!
   END DO
END DO
!
! 5. apply the limit to the fluxes where needed
!
ZFCOR(:,IJB:IJE+1,:) = MAX( &
     MIN(ZRNEG(:,IJB:IJE+1,:),ZRPOS(:,IJB-1:IJE,:)) * ZFCOR(:,IJB:IJE+1,:), &
     ZFCOR(:,IJB:IJE+1,:) )
ZFCOR(:,IJB-1,:) = MAX( &
     MIN(ZRNEG(:,IJB-1,:),ZRPOS(:,IJE-1,:))*ZFCOR(:,IJB-1,:),ZFCOR(:,IJB-1,:))
!
ZFCOR(:,IJB:IJE+1,:) = MIN( &
     MIN(ZRPOS(:,IJB:IJE+1,:),ZRNEG(:,IJB-1:IJE,:)) * ZFCOR(:,IJB:IJE+1,:), &
     ZFCOR(:,IJB:IJE+1,:) )
ZFCOR(:,IJB-1,:) = MIN( &
     MIN(ZRPOS(:,IJB-1,:),ZRNEG(:,IJE-1,:))*ZFCOR(:,IJB-1,:),ZFCOR(:,IJB-1,:))
!
!-------------------------------------------------------------------------------
! 6. apply the limited flux correction to scalar field
!
PR = PR - PTSTEP*DYF(ZFCOR)
!
!
END FUNCTION PPM_S1_Y
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION PPM_S1_Z(KGRID, PSRC, PCR, PRHO, PRHOT, PTSTEP) &
                        RESULT(PR)
!     ########################################################################
!!
!!****  PPM_S1_Z - PPM  advection scheme in Z direction in Skamarock 2006 
!!                 notation - with flux limiting for monotonicity
!!
!!    MODIFICATIONS
!!    -------------
!!    
!!    23.6.2006.  T. Maric - original version
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODI_SHUMAN
!
USE MODD_CONF
USE MODD_PARAMETERS
!USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                INTENT(IN)  :: KGRID   ! C grid localisation
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC    ! variable at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PCR     ! Courant number
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHO  ! density
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHOT ! density at t+dt
REAL,                   INTENT(IN)  :: PTSTEP  ! Time step 
!
! output source term
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER:: IIB,IJB,IKB   ! Begining useful area in x,y,z directions
INTEGER:: IIE,IJE,IKE   ! End useful area in x,y,z directions
INTEGER:: IKU
!
! variable at cell edges
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZPHAT, ZRVT
!
! advection fluxes, upwind and correction
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZFUP, ZFCOR
!
! ratios for limiting the correction flux
REAL, DIMENSION(SIZE(PCR,1),SIZE(PCR,2),SIZE(PCR,3)) :: ZRPOS, ZRNEG
!
! variables for limiting the correction flux
REAL :: ZSRCMAX, ZSRCMIN, ZFIN, ZFOUT
!
REAL, PARAMETER :: ZEPS = 1.0E-16
!
INTEGER :: II, IJ, IK
!
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PSRC,3) - JPVEXT
IKU =  SIZE(PSRC,3)
!
!-------------------------------------------------------------------------------
!
ZRVT = PCR/PTSTEP * MZM(1,IKU,1,PRHO)
!
! calculate 4th order fluxes at cell edges in the inner domain !
ZPHAT(:,:,IKB+1:IKE) = (7.0 * &
                       (PSRC(:,:,IKB+1:IKE) + PSRC(:,:,IKB:IKE-1)) - &
                       (PSRC(:,:,IKB+2:IKE+1) + PSRC(:,:,IKB-1:IKE-2))) / 12.0
!
! set BC to WALL
!
ZPHAT(:,:,IKB) = (7.0 * &
                 (PSRC(:,:,IKB) + PSRC(:,:,IKB+1)) - &
                 (PSRC(:,:,IKB+1) + PSRC(:,:,IKB+2))) / 12.0
ZPHAT(:,:,IKB-1) = ZPHAT(:,:,IKB+1)
ZPHAT(:,:,IKE+1) = (7.0 * &
                   (PSRC(:,:,IKE+1) + PSRC(:,:,IKE)) - &
                   (PSRC(:,:,IKE) + PSRC(:,:,IKE-1))) / 12.0
!
! set BC to OPEN
!
!!$ZPHAT(:,:,IKB) = 0.5*(PSRC(:,:,IKB-1) + PSRC(:,:,IKB))
!!$ZPHAT(:,:,IKB-1) = ZPHAT(:,:,IKB)
!!$ZPHAT(:,:,IKE+1) = 0.5*(PSRC(:,:,IKE) + PSRC(:,:,IKE+1))
!
! calculate upwind and correction fluxes. upwind flux is upstream value of the
! scalar variable, and correction flux is the correction to the upstream flux
! that makes it equivalent to the PPM flux
! flux_ppm = flux_up + flux_corr
!
WHERE ( PCR(:,:,IKB:IKE) .GT. 0.0 )
   ZFUP(:,:,IKB:IKE) = ZRVT(:,:,IKB:IKE) * PSRC(:,:,IKB-1:IKE-1)
   ZFCOR(:,:,IKB:IKE) = ZRVT(:,:,IKB:IKE) * &
        (1.0 - PCR(:,:,IKB:IKE)) * &
        (ZPHAT(:,:,IKB:IKE) - PSRC(:,:,IKB-1:IKE-1) - PCR(:,:,IKB:IKE) * &
        (ZPHAT(:,:,IKB-1:IKE-1) - 2.0*PSRC(:,:,IKB-1:IKE-1)+ZPHAT(:,:,IKB:IKE)))
ELSEWHERE
   ZFUP(:,:,IKB:IKE) = ZRVT(:,:,IKB:IKE) * PSRC(:,:,IKB:IKE)
   ZFCOR(:,:,IKB:IKE) = ZRVT(:,:,IKB:IKE) * &
        (1.0 + PCR(:,:,IKB:IKE)) * &
        (ZPHAT(:,:,IKB:IKE) - PSRC(:,:,IKB:IKE) + PCR(:,:,IKB:IKE) * &
        (ZPHAT(:,:,IKB:IKE) - 2.0*PSRC(:,:,IKB:IKE) + ZPHAT(:,:,IKB+1:IKE+1)))
END WHERE
!
! set BC to WALL
!
WHERE ( PCR(:,:,IKB-1) .GT. 0.0 )
   ZFUP(:,:,IKB-1) = ZRVT(:,:,IKB-1) * PSRC(:,:,IKB+2)
   ZFCOR(:,:,IKB-1) =  ZRVT(:,:,IKB-1) * &
        (1.0 - PCR(:,:,IKB-1)) * &
        (ZPHAT(:,:,IKB+1) - PSRC(:,:,IKB+2) - PCR(:,:,IKB+1) * &
        (ZPHAT(:,:,IKB+2) - 2.0*PSRC(:,:,IKB+2) + ZPHAT(:,:,IKB+1)))
ELSEWHERE
   ZFUP(:,:,IKB-1) = ZRVT(:,:,IKB-1) * PSRC(:,:,IKB+1)
   ZFCOR(:,:,IKB-1) =  ZRVT(:,:,IKB-1) * &
        (1.0 + PCR(:,:,IKB-1)) * &
        (ZPHAT(:,:,IKB+1) - PSRC(:,:,IKB+1) + PCR(:,:,IKB+1) * &
        (ZPHAT(:,:,IKB+1) - 2.0*PSRC(:,:,IKB+1) + ZPHAT(:,:,IKB)))
END WHERE
!
WHERE ( PCR(:,:,IKE+1) .GT. 0.0 )
   ZFUP(:,:,IKE+1) = ZRVT(:,:,IKE+1) * PSRC(:,:,IKE)
   ZFCOR(:,:,IKE+1) =  ZRVT(:,:,IKE+1) * &
        (1.0 - PCR(:,:,IKE+1)) * &
        (ZPHAT(:,:,IKE+1) - PSRC(:,:,IKE) - PCR(:,:,IKE+1) * &
        (ZPHAT(:,:,IKE) - 2.0*PSRC(:,:,IKE) + ZPHAT(:,:,IKE+1)))
ELSEWHERE
   ZFUP(:,:,IKE+1) = ZRVT(:,:,IKE+1) * PSRC(:,:,IKE+1)
   ZFCOR(:,:,IKE+1) =  ZRVT(:,:,IKE+1) * &
        (1.0 + PCR(:,:,IKE+1)) * &
        (ZPHAT(:,:,IKE+1) - PSRC(:,:,IKE+1) + PCR(:,:,IKE+1) * &
        (ZPHAT(:,:,IKE+1) - 2.0*PSRC(:,:,IKE+1) + ZPHAT(:,:,IKE)))
END WHERE
!
!
!!$! set boundaries to CYCL
!!$!
!!$WHERE ( PCR(:,:,IKB-1) .GT. 0.0 )
!!$   ZFUP(:,:,IKB-1) = ZRVT(:,:,IKB-1) * PSRC(:,:,IKE-1)
!!$   ZFCOR(:,:,IKB-1) =  ZRVT(:,:,IKB-1) * &
!!$        (1.0 - PCR(:,:,IKB-1)) * &
!!$        (ZPHAT(:,:,IKB-1) - PSRC(:,:,IKE-1) - PCR(:,:,IKB-1) * &
!!$        (ZPHAT(:,:,IKE-1) - 2.0*PSRC(:,:,IKE-1) + ZPHAT(:,:,IKB-1)))
!!$ELSEWHERE
!!$   ZFUP(:,:,IKB-1) = ZRVT(:,:,IKB-1) * PSRC(:,:,IKB-1)
!!$   ZFCOR(:,:,IKB-1) =  ZRVT(:,:,IKB-1) * &
!!$        (1.0 + PCR(:,:,IKB-1)) * &
!!$        (ZPHAT(:,:,IKB-1) - PSRC(:,:,IKB-1) + PCR(:,:,IKB-1) * &
!!$        (ZPHAT(:,:,IKB-1) - 2.0*PSRC(:,:,IKB-1) + ZPHAT(:,:,IKB)))
!!$END WHERE
!!$!
!!$WHERE ( PCR(:,:,IKE+1) .GT. 0.0 )
!!$   ZFUP(:,:,IKE+1) = ZRVT(:,:,IKE+1) * PSRC(:,:,IKE)
!!$   ZFCOR(:,:,IKE+1) =  ZRVT(:,:,IKE+1) * &
!!$        (1.0 - PCR(:,:,IKE+1)) * &
!!$        (ZPHAT(:,:,IKE+1) - PSRC(:,:,IKE) - PCR(:,:,IKE+1) * &
!!$        (ZPHAT(:,:,IKE) - 2.0*PSRC(:,:,IKE) + ZPHAT(:,:,IKE+1)))
!!$ELSEWHERE
!!$   ZFUP(:,:,IKE+1) = ZRVT(:,:,IKE+1) * PSRC(:,:,IKE+1)
!!$   ZFCOR(:,:,IKE+1) =  ZRVT(:,:,IKE+1) * &
!!$        (1.0 + PCR(:,:,IKE+1)) * &
!!$        (ZPHAT(:,:,IKE+1) - PSRC(:,:,IKE+1) + PCR(:,:,IKE+1) * &
!!$        (ZPHAT(:,:,IKE+1) - 2.0*PSRC(:,:,IKE+1) + ZPHAT(:,:,IKB+1)))
!!$END WHERE
!
! Perform limiting of the fluxes
!
! 1. calculate upwind tendency of the source
!
PR = PSRC*PRHO - PTSTEP*DZF(1,IKU,1,ZFUP)
!
!-------------------------------------------------------------------------------
! compute and apply the limiters
!
DO II = IIB-1,IIE+1
   DO IJ = IJB-1,IJE+1
      DO IK = IKB,IKE         
!
! 2. find local extrema in the source 
!
         ZSRCMAX = MAX( PSRC(II,IJ,IK-1),PSRC(II,IJ,IK),PSRC(II,IJ,IK+1) )
         ZSRCMIN = MIN( PSRC(II,IJ,IK-1),PSRC(II,IJ,IK),PSRC(II,IJ,IK+1) )
!
! 3. compute incoming and outgoing fluxes for this cell
!
         ZFOUT = MAX(ZEPS,MAX(0.,ZFCOR(II,IJ,IK+1)) - MIN(0.,ZFCOR(II,IJ,IK)))
         ZFIN  = MAX(ZEPS,MAX(0.,ZFCOR(II,IJ,IK)) - MIN(0.,ZFCOR(II,IJ,IK+1)))
!
! 4. calculate fraction of outgoing and incoming flux which will drive scalar
!    values outside the local extrema
!
         ZRNEG(II,IJ,IK) = MAX(0.,MIN(1., &
              (PR(II,IJ,IK) - PRHOT(II,IJ,IK)*ZSRCMIN) &
              / PTSTEP / ZFOUT))
!
         ZRPOS(II,IJ,IK) = MAX(0.,MIN(1., &
              (PRHOT(II,IJ,IK)*ZSRCMAX - PR(II,IJ,IK)) &
              / PTSTEP / ZFIN))
      END DO
   END DO
END DO
!
! set WALL boundaries
!
DO II = IIB-1,IIE+1
   DO IJ = IJB-1,IJE+1         
!
      ZSRCMAX = MAX( PSRC(II,IJ,IKB+1),PSRC(II,IJ,IKB) )
      ZSRCMIN = MIN( PSRC(II,IJ,IKB+1),PSRC(II,IJ,IKB) )
!
      ZFOUT = MAX(ZEPS,MAX(0.,ZFCOR(II,IJ,IKB)) - MIN(0.,ZFCOR(II,IJ,IKB-1)))
      ZFIN  = MAX(ZEPS,MAX(0.,ZFCOR(II,IJ,IKB-1)) - MIN(0.,ZFCOR(II,IJ,IKB)))
!
      ZRNEG(II,IJ,IKB-1) = MAX(0.,MIN(1., &
           (PR(II,IJ,IKB-1) - PRHOT(II,IJ,IKB-1)*ZSRCMIN) &
           / PTSTEP / ZFOUT))
!
      ZRPOS(II,IJ,IKB-1) = MAX(0.,MIN(1., &
           (PRHOT(II,IJ,IKB-1)*ZSRCMAX - PR(II,IJ,IKB-1)) &
           / PTSTEP / ZFIN))
!
! 
      ZSRCMAX = MAX( PSRC(II,IJ,IKE),PSRC(II,IJ,IKE+1) )
      ZSRCMIN = MIN( PSRC(II,IJ,IKE),PSRC(II,IJ,IKE+1) )
!
      ZFOUT = MAX(ZEPS,MAX(0.,ZFCOR(II,IJ,IKE)) - MIN(0.,ZFCOR(II,IJ,IKE+1)))
      ZFIN  = MAX(ZEPS,MAX(0.,ZFCOR(II,IJ,IKE+1)) - MIN(0.,ZFCOR(II,IJ,IKE)))
!
      ZRNEG(II,IJ,IKE+1) = MAX(0.,MIN(1., &
           (PR(II,IJ,IKE+1) - PRHOT(II,IJ,IKE+1)*ZSRCMIN) &
           / PTSTEP / ZFOUT))
!
      ZRPOS(II,IJ,IKE+1) = MAX(0.,MIN(1., &
           (PRHOT(II,IJ,IKE+1)*ZSRCMAX - PR(II,IJ,IKE+1)) &
           / PTSTEP / ZFIN))
!
   END DO
END DO
!
! set CYCL boundaries
!
!!$DO II = IIB-1,IIE+1
!!$   DO IJ = IJB-1,IJE+1         
!!$!
!!$      ZSRCMAX = MAX( PSRC(II,IJ,IKE-1),PSRC(II,IJ,IKB-1),PSRC(II,IJ,IKB) )
!!$      ZSRCMIN = MIN( PSRC(II,IJ,IKE-1),PSRC(II,IJ,IKB-1),PSRC(II,IJ,IKB) )
!!$!
!!$      ZFOUT = MAX(ZEPS,MAX(0.,ZFCOR(II,IJ,IKB)) - MIN(0.,ZFCOR(II,IJ,IKB-1)))
!!$      ZFIN  = MAX(ZEPS,MAX(0.,ZFCOR(II,IJ,IKB-1)) - MIN(0.,ZFCOR(II,IJ,IKB)))
!!$!
!!$      ZRNEG(II,IJ,IKB-1) = MAX(0.,MIN(1., &
!!$           (PR(II,IJ,IKB-1) - PRHOT(II,IJ,IKB-1)*ZSRCMIN) &
!!$           / PTSTEP / ZFOUT))
!!$!
!!$      ZRPOS(II,IJ,IKB-1) = MAX(0.,MIN(1., &
!!$           (PRHOT(II,IJ,IKB-1)*ZSRCMAX - PR(II,IJ,IKB-1)) &
!!$           / PTSTEP / ZFIN))
!!$!
!!$! 
!!$      ZSRCMAX = MAX( PSRC(II,IJ,IKE),PSRC(II,IJ,IKE+1),PSRC(II,IJ,IKB+1) )
!!$      ZSRCMIN = MIN( PSRC(II,IJ,IKE),PSRC(II,IJ,IKE+1),PSRC(II,IJ,IKB+1) )
!!$!
!!$      ZFOUT = MAX(ZEPS,MAX(0.,ZFCOR(II,IJ,IKB+1)) - MIN(0.,ZFCOR(II,IJ,IKE+1)))
!!$      ZFIN  = MAX(ZEPS,MAX(0.,ZFCOR(II,IJ,IKE+1)) - MIN(0.,ZFCOR(II,IJ,IKB+1)))
!!$!
!!$      ZRNEG(II,IJ,IKE+1) = MAX(0.,MIN(1., &
!!$           (PR(II,IJ,IKE+1) - PRHOT(II,IJ,IKE+1)*ZSRCMIN) &
!!$           / PTSTEP / ZFOUT))
!!$!
!!$      ZRPOS(II,IJ,IKE+1) = MAX(0.,MIN(1., &
!!$           (PRHOT(II,IJ,IKE+1)*ZSRCMAX - PR(II,IJ,IKE+1)) &
!!$           / PTSTEP / ZFIN))
!!$!
!!$   END DO
!!$END DO
!
! 5. apply the limit to the fluxes where needed
!
ZFCOR(:,:,IKB:IKE+1) = MAX( &
     MIN(ZRNEG(:,:,IKB:IKE+1),ZRPOS(:,:,IKB-1:IKE)) * ZFCOR(:,:,IKB:IKE+1), &
     ZFCOR(:,:,IKB:IKE+1) )
ZFCOR(:,:,IKB-1) = MAX( &
     MIN(ZRNEG(:,:,IKB-1),ZRPOS(:,:,IKB+2))*ZFCOR(:,:,IKB-1),ZFCOR(:,:,IKB-1))
!!$ZFCOR(:,:,IKB-1) = MAX( &
!!$     MIN(ZRNEG(:,:,IKB-1),ZRPOS(:,:,IKE-1))*ZFCOR(:,:,IKB-1),ZFCOR(:,:,IKB-1))
!
ZFCOR(:,:,IKB:IKE+1) = MIN( &
     MIN(ZRPOS(:,:,IKB:IKE+1),ZRNEG(:,:,IKB-1:IKE)) * ZFCOR(:,:,IKB:IKE+1), &
     ZFCOR(:,:,IKB:IKE+1) )
ZFCOR(:,:,IKB-1) = MIN( &
     MIN(ZRPOS(:,:,IKB-1),ZRNEG(:,:,IKB+2))*ZFCOR(:,:,IKB-1),ZFCOR(:,:,IKB-1))
!!$ZFCOR(:,:,IKB-1) = MIN( &
!!$     MIN(ZRPOS(:,:,IKB-1),ZRNEG(:,:,IKE-1))*ZFCOR(:,:,IKB-1),ZFCOR(:,:,IKB-1))
!
!-------------------------------------------------------------------------------
! 6. apply the limited flux correction to scalar field
!
PR = PR - PTSTEP*DZF(1,IKU,1,ZFCOR)
!
!
END FUNCTION PPM_S1_Z
