!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_EDDY_FLUX_n
!     #######################
!
INTERFACE
!
      SUBROUTINE EDDY_FLUX_n (KMI,KTCOUNT,PVM,PTHM,PRHODJ,PRTHS,&
                                      PVTH_FLUX_M,PWTH_FLUX_M )
!
!
INTEGER,                  INTENT(IN)    :: KMI   ! Model index
INTEGER,                  INTENT(IN)    :: KTCOUNT  ! iteration count
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ            ! dry density of
!                                 ! anelastic reference state * Jacobian
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PVM
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHM
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PVTH_FLUX_M
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PWTH_FLUX_M
!

END SUBROUTINE EDDY_FLUX_n
!
END INTERFACE
!
END MODULE MODI_EDDY_FLUX_n
!
!     ##################################################################################
      SUBROUTINE EDDY_FLUX_n (KMI,KTCOUNT,PVM,PTHM,PRHODJ,PRTHS,PVTH_FLUX_M,PWTH_FLUX_M)
!     ##################################################################################
!
!!    PURPOSE
!!    -------
!!    Baroclinic fluxes parameterization (v'T' and w'T') for 2D transect (latitude, altitude) model
!!
!!**  METHOD
!!    ------
!!   
!!    The points in the domain where baroclinic  instability is reached are first searched for using the
!!    criterion that (Coriolis parameter) * (vertical mean of meridional gradient of potential temperature) be negative. 
!!    It corresponds to vertical wind shear suitable for a baroclinic to grow and be unstable. 
!!    The eddy flux are v'T'= Kyy*(DT/dy) ; w'T'=Kyz*(DT_dy). The main purpose
!!    of the routine is to derive the expression of the Kyy coefficient from linear
!!    growth of a baroclinic wave following Zou and Gal-chen (1991). 

!!    Ultimately, the tendency is derived for the potential temperature 
!!
!!    REFERENCE
!!    ---------
!! 
!!    Zou and Gal-Chen (1991) : 
!!    Peyrillé et al. (2007) : 
!!
!!
!!    AUTHOR
!!    ------
!!	  P.Peyrille          * Meteo-France *
!!
!!   EXTERNAL: MEAN_Z is used to perform vertical weighted average
!!
!!    MODIFICATIONS
!!    -------------
!!      Original  18/02/04
!!    M.Tomasini 05/05/09 Grid-nesting
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!!
USE MODE_ll
!
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_BUDGET
!
USE MODI_SHUMAN
USE MODI_BUDGET
USE MODD_CST
!
USE MODD_DIM_n
USE MODD_CONF
USE MODD_CONF_n
USE MODD_GRID
USE MODD_GRID_n
USE MODD_METRICS_n
USE MODD_TIME
USE MODD_TIME_n
USE MODD_DYN_n
USE MODD_FIELD_n
USE MODD_CURVCOR_n
USE MODI_GRADIENT_M
USE MODI_GRADIENT_W
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_SHUMAN
USE MODE_GRIDPROJ
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_REF
USE MODD_LATZ_EDFLX
!
USE MODI_MEAN_Z

!

IMPLICIT NONE
!
INTEGER,                INTENT(IN)   :: KMI     ! Model index
INTEGER,                INTENT(IN)   :: KTCOUNT ! iteration count
!
REAL, DIMENSION(:,:,:), INTENT(IN)   :: PRHODJ  ! dry density of
!                         ! anelastic reference state * Jacobian
REAL, DIMENSION(:,:,:), INTENT(IN)   :: PVM
!
REAL, DIMENSION(:,:,:), INTENT(IN)   :: PTHM
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PRTHS
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PVTH_FLUX_M
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PWTH_FLUX_M
 
!*       0.2   Declarations of local variables :
!
INTEGER:: IIB,IJB        ! Begining useful area  in x,y directions
INTEGER:: IIE,IJE        ! End useful area in x,y directions
INTEGER:: JJ,JI,JK,IKB,IKE,IIU,IKU
INTEGER:: ILUOUT         ! Logical unit number for the output listing
INTEGER:: IRESP          ! return code in FM routines
INTEGER:: NI_INFN,NI_INFS
!
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZCORIOZ,ZTHM
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZBETA        ! f and Beta 
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZH           ! density scale height
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZND,ZWORK32  ! Brunt Vaisala frequency
REAL, DIMENSION(:,:)  , ALLOCATABLE:: ZWORK31      ! Brunt Vaisala frequency
!
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZDTHM_DY,ZDTHM_DZ ! gradients of theta
REAL:: ZDELTAZ
REAL:: ZRC ! to define ZKC according to (ZOU, Gal-Chen,99) = 1.83
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZVTH_FLUX, ZDIV_YTHFLUX
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZKYY, ZKYZ,ZGAMMAF2
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZADTDX
REAL,DIMENSION(:,:),    ALLOCATABLE:: ZDTHM_DYW,ZDTHM_DZW
REAL,DIMENSION(:,:),    ALLOCATABLE:: ZKCW,ZTHM_W
REAL,DIMENSION(:,:),    ALLOCATABLE:: ZHW, ZNDW,ZDKW,ZGAMMAW
REAL, DIMENSION(:,:)  , ALLOCATABLE:: ZADTDXW
REAL, DIMENSION(:,:)  , ALLOCATABLE:: ZGAMMAF2W,ZDW
REAL, DIMENSION(:)    , ALLOCATABLE:: ZA0
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZZ1,ZPOLI,ZWTH_FLUX,ZDIV_ZTHFLUX
REAL:: ZRF
!
! ===============================================================                         
!*       1. INITIALISATION
! ===============================================================                         

CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)

IKU=SIZE(XZHAT) ! nb points sur Z
IKB=1+JPVEXT
IKE=IKU-JPVEXT

IIU = SIZE(XXHAT)

! RECOVER THE LOGICAL UNIT NUMBER FOR THE OUTPUT PRINTS
ILUOUT = TLUOUT%NLU
!
!       1.1 ALLOCATION DES TABLEAUX
!           -----------------------
ALLOCATE(ZCORIOZ(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)))
ALLOCATE(ZTHM(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)))
ALLOCATE(ZADTDX(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)))
ALLOCATE(ZADTDXW(SIZE(PTHM,1),SIZE(PTHM,2)))


ALLOCATE(ZBETA(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3))) 
ALLOCATE(ZH(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3))) 
ALLOCATE(ZND(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)))
ALLOCATE(ZWORK32(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)))
!
ALLOCATE(ZDTHM_DY(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)))
ALLOCATE(ZDTHM_DZ(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)))
ALLOCATE(ZDTHM_DZW(SIZE(PTHM,1),SIZE(PTHM,2)))

ALLOCATE(ZHW(SIZE(PTHM,1),SIZE(PTHM,2))  )
ALLOCATE(ZWORK31(SIZE(PTHM,1),SIZE(PTHM,2))  )
ALLOCATE(ZNDW(SIZE(PTHM,1),SIZE(PTHM,2))  )
ALLOCATE(ZDKW(SIZE(PTHM,1),SIZE(PTHM,2)) ) 
ALLOCATE(ZKCW(SIZE(PTHM,1),SIZE(PTHM,2)) ) 
ALLOCATE(ZTHM_W(SIZE(PTHM,1),SIZE(PTHM,2)))
ALLOCATE(ZDTHM_DYW(SIZE(PTHM,1),SIZE(PTHM,2)) )
ALLOCATE(ZVTH_FLUX(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)))
ALLOCATE(ZDIV_YTHFLUX(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)))
ALLOCATE(ZKYY(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3))) 
ALLOCATE(ZKYZ(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3))) 
ALLOCATE(ZGAMMAF2(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3))) 
ALLOCATE(ZGAMMAF2W(SIZE(PTHM,1),SIZE(PTHM,2))) 
ALLOCATE(ZGAMMAW(SIZE(PTHM,1),SIZE(PTHM,2))) 
ALLOCATE(ZDW(SIZE(PTHM,1),SIZE(PTHM,2))) 
ALLOCATE(ZA0(SIZE(PTHM,1)))
!
ALLOCATE(ZZ1(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)))
ALLOCATE(ZPOLI(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)))
ALLOCATE(ZWTH_FLUX(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)))
ALLOCATE(ZDIV_ZTHFLUX(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)))

!       1.2  FIRST COMPUTATIONS
!          ------------------
! Intialisation and set value for coefficient in the expresion of Kyy
ZKYY(:,:,:)= 0.
ZGAMMAW(:,:)=0.
ZDW(:,:)=0.
ZDKW(:,:)=0.
ZKCW(:,:)=0.
ZA0(:)= XTH_FLX !0.74 
ZRC = 1.83
!
! assumed BL height
ZDELTAZ=1000.
! Init northern and southern point of the domain
! 
! Modif MT
ZVTH_FLUX=0.
!
NI_INFN = IIE
NI_INFS = IIB

!Coriolis parmeter and Beta
ZCORIOZ(:,:,:)= SPREAD(XCORIOZ(:,:),3,IKU)
ZBETA(:,:,:) = GX_M_U(1,IKU,1,ZCORIOZ(:,:,:),XDXX,XDZZ,XDZX)
ZCORIOZ(:,:,:)= MXM(ZCORIOZ(:,:,:))
! Dry Brunt Vaisal frequency

ZWORK32(:,:,:)=DZM(1,IKU,1,PTHM(:,:,:))/ MZM(1,IKU,1,PTHM(:,:,:))
DO JK=1,(IKE+1)
   DO JJ=1,(IJE+1)
      DO JI=1,(IIE+1)
         IF(ZWORK32(JI,JJ,JK)<0.) THEN
           ZND(JI,JJ,JK)= -1.*SQRT( ABS( XG*ZWORK32(JI,JJ,JK)/ XDZZ(JI,JJ,JK) ) )
         ELSE
           ZND(JI,JJ,JK)= SQRT( ABS( XG*ZWORK32(JI,JJ,JK)/ XDZZ(JI,JJ,JK) ) )
         ENDIF
      ENDDO
   ENDDO
ENDDO
ZND(:,:,:) = MXM(ZND(:,:,:))
!! latitudinal gradient of TH
ZDTHM_DY(:,:,:) = GX_M_U(1,IKU,1,PTHM,XDXX,XDZZ,XDZX)
ZDTHM_DZ(:,:,:) = MXM(GZ_M_M(1,IKU,1,PTHM,XDZZ))
! density scale height
ZH(:,:,:) = PTHM(:,:,:) * XRD * (XG**(-1))
ZH(:,:,:) = MXM(ZH)
!
! COMPUTE VERTICAL WEIGHTED AVERAGE 
!< ZND >
 ZNDW(:,:) = MW_Z(ZND,PRHODJ)
! < DTH/DY >
 ZDTHM_DYW(:,:) = MW_Z(ZDTHM_DY,PRHODJ)
 ZADTDX(:,:,:) = ABS(ZDTHM_DY(:,:,:))
 ZADTDXW(:,:) = MW_Z(ZADTDX,PRHODJ)
! < DTH/DZ >
 ZDTHM_DZW(:,:) = MW_Z(ZDTHM_DZ,PRHODJ)
! < TH >
 ZTHM_W(:,:) = MW_Z(PTHM,PRHODJ)
!< H >
 ZHW(:,:) = MW_Z(ZH,PRHODJ)
!
!        1.3 COMPUTaTION OF  V'T'
!           -------------------
!
DO JJ=IJB,IJE
   DO JI=IIB,IIE

  ! WE CONSIDER ONLY LATITUDES GREATER THAN 10N AND 10S        
  !
  IF ( ((ABS( 0.5*XLAT(JI,JJ)+0.5*XLAT(JI-1,JJ) ))>10.) ) THEN

   ZWORK31(JI,JJ) = (ZCORIOZ(JI,JJ,2)*ZDTHM_DYW(JI,JJ))
    ! IF (f.dtheta_dy) < 0 => eddy flux is computed
    IF ( ZWORK31(JI,JJ)< (-0.1E-10) ) THEN
      ! < GAMMA >        
      ZGAMMAW(JI,JJ)= &
      ABS( (ZBETA(JI,JJ,2)*ZHW(JI,JJ)*ZDTHM_DZW(JI,JJ))/ &
          (ZCORIOZ(JI,JJ,2)* ZDTHM_DYW(JI,JJ) ) )
      !
      !mean < Kc >
      ZKCW(JI,JJ) =SQRT( ABS( (((1+ZGAMMAW(JI,JJ))/ZRC)**2) - 0.25 ))
      ! 
      ! < dk >
      ZDKW(JI,JJ) = ZHW(JI,JJ) / &
         ( ( ((4*(ZKCW(JI,JJ)**(2))) + 1  )**(0.5)) -1)     
      !
      ! < D >
      IF (ZKCW(JI,JJ) > 1.E-10 ) THEN
       ZDW(JI,JJ) = ZHW(JI,JJ) / ZKCW(JI,JJ)
      ELSE
       ZDW(JI,JJ) = 5.
      END IF
      !
      ! 1.4 CALCUL KYY TRANSFERTS COEFFICIENTS 
      !     ----------------------------------
      DO JK=IKB,IKE
      !
        ZKYY(JI,JJ,JK) =  &
            (ZA0(JI) * XG*ZNDW(JI,JJ)/(ZTHM_W(JI,JJ)*(ZCORIOZ(JI,JJ,JK)**2)) ) &
            * (ZDW(JI,JJ)**2) * (EXP(-0.5*(XZHAT(JK)+XZHAT(JK+1))/ZDW(JI,JJ))) &
            * ZADTDXW(JI,JJ)   
      ENDDO       
      !
    ! CASE WHERE NO INSTABILITY IS DETECTED
    ELSE
      ZKYY(JI,JJ,:)=0.
    ENDIF
  ! CASE WHERE LAT BETWEEN 10S AND 10N 
  ELSE
    ZKYY(JI,JJ,:) = 0.
  ENDIF   
!
  ENDDO
ENDDO
!
!   2. COMPUTE  V'T' AND DIV(V'T')
!  ---------------------------------------------------------------

!
DO JK=IKB,IKE
   ZVTH_FLUX(:,:,JK) = - 0.5 * ZKYY(:,:,JK)*ZDTHM_DY(:,:,JK) * &
                    (1-EXP(-0.5*(XZHAT(JK)+XZHAT(JK+1))/ZDELTAZ))
END DO
!
!       2.1 Smoothing in equatorial region
!           ------------------------------

DO JI=IIB,IIE
   DO JJ=IJB,IJE
      IF ((0.5*(XLAT(JI,JJ)+XLAT(JI-1,JJ))>10.).AND.(0.5*(XLAT(JI,JJ)+XLAT(JI-1,JJ))<11.)) THEN
         NI_INFN=JI
      ENDIF
        
      IF ((0.5*(XLAT(JI,JJ)+XLAT(JI-1,JJ))<-10.).AND.(0.5*(XLAT(JI,JJ)+XLAT(JI-1,JJ))>-11.))THEN
         NI_INFS=JI
      ENDIF
   ENDDO
ENDDO
 
DO JI=IIB,IIE
   DO JJ=IJB,IJE
! 
      IF ( (0.5*(XLAT(JI,JJ)+XLAT(JI-1,JJ))>=0.).AND.&
           (0.5*(XLAT(JI,JJ)+XLAT(JI-1,JJ))<0.5*(XLAT(NI_INFN,JJ)+XLAT(NI_INFN-1,JJ))) ) THEN

         ZVTH_FLUX(JI,JJ,:) = ZVTH_FLUX(JI,JJ,:) + ZVTH_FLUX(NI_INFN,JJ,:)* &
         EXP( - (((0.5*(XLAT(NI_INFN,JJ)+XLAT(NI_INFN-1,JJ)) - 0.5*(XLAT(JI,JJ)+XLAT(JI-1,JJ)))/(5.))**(2)))

      ENDIF
        
      IF ( (0.5*(XLAT(JI,JJ)+XLAT(JI-1,JJ))<=0.).AND.&
           (0.5*(XLAT(JI,JJ)+XLAT(JI-1,JJ))>0.5*(XLAT(NI_INFS,JJ)+XLAT(NI_INFS-1,JJ))) ) THEN

        ZVTH_FLUX(JI,JJ,:) = ZVTH_FLUX(JI,JJ,:) + ZVTH_FLUX(NI_INFS,JJ,:)*&
        EXP( - (((0.5*(XLAT(NI_INFN,JJ)+XLAT(NI_INFN-1,JJ)) - 0.5*(XLAT(JI,JJ)+XLAT(JI-1,JJ)))/(5.))**(2)))

      ENDIF
   ENDDO
ENDDO   

!       2.2 BOUNDARY CONDITIONS FOR V'T'       
!           ---------------------------
ZVTH_FLUX((IIE),:,:) =  ZVTH_FLUX((IIE-1),:,:)       
ZVTH_FLUX((IIB),:,:) =  ZVTH_FLUX((IIB+1),:,:) 
ZVTH_FLUX((IIU),:,:) =  ZVTH_FLUX((IIE-1),:,:)       
ZVTH_FLUX(1,:,:)     =  ZVTH_FLUX((IIB+1),:,:)         
! 
!       3. COMPUTE  W'T'
!       ----------------

ZWTH_FLUX=0.
DO JK=IKB,IKE
   DO JI=IIB,IIE
      DO JJ=IJB,IJE
!
         IF ( (ZKYY(JI,JJ,JK)) /= 0. ) THEN
!
!            placed at a W point 
             ZZ1  (JI,JJ,JK) = XZHAT(JK)/( 0.5*ZDKW(JI,JJ)+0.5*ZDKW(JI-1,JJ) )
             ZPOLI(JI,JJ,JK) = ZZ1(JI,JJ,JK) - 0.25*(ZZ1(JI,JJ,JK)**2)
             ZWTH_FLUX(JI,JJ,JK) =  &
(-(0.5*ZDTHM_DYW(JI,JJ)+0.5*ZDTHM_DYW(JI-1,JJ))/(0.5*ZDTHM_DZW(JI,JJ)+0.5*ZDTHM_DZW(JI-1,JJ))) &
*ZPOLI(JI,JJ,JK)* 0.5 * ( 0.5*(ZVTH_FLUX(JI,JJ,JK)+ZVTH_FLUX(JI-1,JJ,JK)) + & 
                          0.5*(ZVTH_FLUX(JI,JJ,JK-1)+ZVTH_FLUX(JI-1,JJ,JK-1)) )
!
         ELSE
             ZWTH_FLUX(JI,JJ,JK) = 0.
         END IF
      ENDDO
   ENDDO
ENDDO
!
! Boundary Conditions
! -------------------
ZWTH_FLUX(:,:,IKB+1) = 0.
ZWTH_FLUX(:,:,IKB)   = 0.
ZWTH_FLUX(:,:,IKE-1) = 0.
ZWTH_FLUX(:,:,IKE)   = 0.
ZWTH_FLUX(:,:,IKB)   = 0.
ZWTH_FLUX(:,:,1  )   = 0.
ZWTH_FLUX(:,:,IKE)   = 0.
ZWTH_FLUX(:,:,IKU)   = 0.
!
!       4. Temporal smoothing and store flux at (t-dt)
!      ----------------------------------------------
!
IF ((KTCOUNT).NE.(1)) THEN
             
  ! V'T'
  !
  ZVTH_FLUX(:,:,:) =  PVTH_FLUX_M(:,:,:)*0.5 + 0.5* ZVTH_FLUX(:,:,:)
  PVTH_FLUX_M(:,:,:) = ZVTH_FLUX(:,:,:)
  !
  ! W'T'
  !
  ZWTH_FLUX(:,:,:) = PWTH_FLUX_M(:,:,:)*0.5 + 0.5* ZWTH_FLUX(:,:,:)
  PWTH_FLUX_M(:,:,:) = ZWTH_FLUX(:,:,:) 
  !
ELSE
  PVTH_FLUX_M(:,:,:) = ZVTH_FLUX(:,:,:)
  PWTH_FLUX_M(:,:,:) = ZWTH_FLUX(:,:,:) 
ENDIF
! 
!       5. DIV v'T' and w'T'
!      --------------------      
! operator GX_U_M used for gradient of v'T' (flux point) placed at a mass point        
!
ZDIV_YTHFLUX(:,:,:) = GX_U_M(1,IKU,1,ZVTH_FLUX,XDXX,XDZZ,XDZX) 
!
ZDIV_ZTHFLUX(:,:,:) = GZ_W_M(1,IKU,1,ZWTH_FLUX,XDZZ)

!
! Control test for the sign of the flux 
DO JI=IIB,IIE
  DO JK=IKB,IKE
     ZRF =  (ABS(PRTHS(JI,2,JK))- ABS(PRHODJ(JI,2,JK)*ZDIV_YTHFLUX(JI,2,JK)))
     IF ( ZRF < 0.) THEN
        WRITE(ILUOUT,FMT=*) 'PB  FLUX SUPERIEUR A THETA ',ZRF , 'JI',JI,'JK',JK
        WRITE(ILUOUT,FMT=*) 'GAMMA', ZGAMMAW(JI,2)
        WRITE(ILUOUT,FMT=*) 'KCW', ZKCW(JI,2)
        WRITE(ILUOUT,FMT=*) 'DK', ZDKW(JI,2)
        WRITE(ILUOUT,FMT=*) 'D', ZDW(JI,2)
        WRITE(ILUOUT,FMT=*) 'BETA', ZBETA(JI,2,JK)
        WRITE(ILUOUT,FMT=*) 'CORIO', ZCORIOZ(JI,2,JK)
        WRITE(ILUOUT,FMT=*) 'ND', ZNDW(JI,2)
        WRITE(ILUOUT,FMT=*) 'KYY', ZKYY(JI,2,JK)
        WRITE(ILUOUT,FMT=*) 'DTHDY', ZDTHM_DY(JI,2,JK)
        WRITE(ILUOUT,FMT=*) 'FLUX', ZVTH_FLUX(JI,2,JK)
        WRITE(ILUOUT,FMT=*) 'DIV', ZDIV_YTHFLUX(JI,2,JK)
        WRITE(ILUOUT,FMT=*) 'RHODJ', PRHODJ(JI,2,JK)
        WRITE(ILUOUT,FMT=*) 'FLUX I-1', ZVTH_FLUX(JI-1,2,JK)
        WRITE(ILUOUT,FMT=*) 'RTHS', PRTHS(JI,2,JK)
     END IF 
  ENDDO
ENDDO
! 
!       6. COMPUTE NEW TEMPERATURE source term
!       -----------------------------------------
PRTHS(:,:,:) = PRTHS(:,:,:) - PRHODJ(:,:,:)* ZDIV_YTHFLUX(:,:,:)
PRTHS(:,:,:) = PRTHS(:,:,:) - PRHODJ(:,:,:)* ZDIV_ZTHFLUX(:,:,:) 

! --------------------------------------------------------------------------

! DEALOCATE ARRAYS 
! ----------------
DEALLOCATE(ZCORIOZ)
DEALLOCATE(ZBETA)
DEALLOCATE(ZH)
DEALLOCATE(ZND)
!
DEALLOCATE(ZDTHM_DY)
DEALLOCATE(ZDTHM_DZ)
DEALLOCATE(ZHW)
DEALLOCATE(ZNDW)
DEALLOCATE(ZDKW)
DEALLOCATE(ZKCW)
DEALLOCATE(ZTHM_W)
DEALLOCATE(ZDTHM_DYW)
DEALLOCATE(ZVTH_FLUX)
DEALLOCATE(ZDIV_YTHFLUX)
DEALLOCATE(ZKYY)
DEALLOCATE(ZKYZ)
DEALLOCATE(ZGAMMAF2)
DEALLOCATE(ZGAMMAF2W)
DEALLOCATE(ZA0)
DEALLOCATE(ZGAMMAW)
DEALLOCATE(ZDW)
DEALLOCATE(ZWORK32)
DEALLOCATE(ZWORK31)
DEALLOCATE(ZTHM)
DEALLOCATE(ZADTDX)
DEALLOCATE(ZADTDXW)
DEALLOCATE(ZZ1)
DEALLOCATE(ZPOLI)
DEALLOCATE(ZWTH_FLUX)
DEALLOCATE(ZDIV_ZTHFLUX)
DEALLOCATE(ZDTHM_DZW)
!
END SUBROUTINE EDDY_FLUX_n
