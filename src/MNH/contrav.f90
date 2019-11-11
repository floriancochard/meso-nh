!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 operators 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_CONTRAV
!     ####################
!
INTERFACE
!
      SUBROUTINE CONTRAV(HLBCX,HLBCY,PRUT,PRVT,PRWT,PDXX,PDYY,PDZZ,PDZX,PDZY,  &
                         PRUCT,PRVCT,PRWCT,KADV_ORDER                       )


CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY  ! Y direction LBC type                         
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PRUT       ! Cartesian comp along x
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PRVT       ! Cartesian comp along y
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PRWT       ! Cartesian comp along z
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PDXX       ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PDYY       ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PDZZ       ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PDZX       ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PDZY       ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(OUT)   ::  PRUCT      ! Contrav comp along x-bar
REAL, DIMENSION(:,:,:),  INTENT(OUT)   ::  PRVCT      ! Contrav comp along y-bar
REAL, DIMENSION(:,:,:),  INTENT(OUT)   ::  PRWCT      ! Contrav comp along z-bar
INTEGER,                 INTENT(IN)    ::  KADV_ORDER ! Order of the advection
                                                      ! scheme
!
END SUBROUTINE CONTRAV
!
END INTERFACE
!
END MODULE MODI_CONTRAV 
!
!
!
!     ##############################################################
      SUBROUTINE CONTRAV(HLBCX,HLBCY,PRUT,PRVT,PRWT,PDXX,PDYY,PDZZ,PDZX,PDZY,  &
                         PRUCT,PRVCT,PRWCT,KADV_ORDER                        )
!     ##############################################################
!
!!****  *CONTRAV * - computes the contravariant components from the
!!       cartesian components
!!
!!    PURPOSE
!!    -------
!       This routine computes the contravariant components of vector
!     defined by its cartesian components (U,V,W) , using the following
!     formulae:
!     UC = U / DXX
!     VC = V / DYY
!               (     ----------x    ----------y )  
!               (           ---z           ---z  )
!           1   (            U              V    )
!     WC = ---  ( W - DZX * ---    - DZY * ---   )
!          DZZ  (           DXX            DYY   )
!
!  
!       In the no-topography case, WC = W / DZZ
!
!
!!**  METHOD
!!    ------
!!      We employ the Shuman operators to compute the averages. The metric
!!    coefficients PDXX, PDYY, PDZX, PDZY, PDZZ are dummy arguments
!!
!!
!!    EXTERNAL 
!!    --------
!!      MXF, MYF, MZM         : Shuman functions (mean operators)
!!       
!!      Module MODI_SHUMAN    : Interface for Shuman functions   
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF   : contains configuration variable 
!!           LFLAT : Logical for topography
!!                  = .TRUE.  if Zs = 0 (Flat terrain)
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (subroutine CONTRAV)
!!
!!
!!    AUTHOR
!!    ------
!!      J.L. Redelsperger     * CNRM *
!!	J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   27/07/94
!!      Corrections 3/08/94 (by J.P. Lafore)
!!      Corrections 17/10/94 (by J.P. Lafore) WC modified for w-advection
!!      Corrections 19/01/11 (by J.P. Pinty) WC 4th order
!!      Corrections 28/03/11 (by V.Masson) // of WC 4th order
!!      J.Escobar 21/03/2013: for HALOK comment all NHALO=1 test
!----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_CONF
USE MODD_PARAMETERS
USE MODD_GRID_n, ONLY: XZZ
!
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
USE MODE_ll
!
USE MODI_SHUMAN
USE MODI_GET_HALO
!
USE MODE_MPPDB
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments    
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
!
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PRUT     ! Cartesian comp along x
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PRVT     ! Cartesian comp along y
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PRWT     ! Cartesian comp along z
!
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PDXX     ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PDYY     ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PDZZ     ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PDZX     ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN)    ::  PDZY     ! Metric coefficients
!
!
REAL, DIMENSION(:,:,:),  INTENT(OUT)   ::  PRUCT    ! Contrav comp along x-bar
REAL, DIMENSION(:,:,:),  INTENT(OUT)   ::  PRVCT    ! Contrav comp along y-bar
REAL, DIMENSION(:,:,:),  INTENT(OUT)   ::  PRWCT    ! Contrav comp along z-bar
INTEGER,                 INTENT(IN)    ::  KADV_ORDER ! Order of the advection
                                                      ! scheme
!
!
!*       0.2   declarations of local variables
!              
REAL, DIMENSION(SIZE(PDXX,1),SIZE(PDXX,2),SIZE(PDXX,3)):: Z1,Z2
INTEGER :: IIB,IIE,IJB,IJE,IKB,IKE
INTEGER :: IIU, IJU, IKU
INTEGER:: IW,IE,IS,IN   ! Coordinate of forth order diffusion area
!
TYPE(LIST_ll),      POINTER :: TZFIELD_U, TZFIELD_V, TZFIELD_DZX, TZFIELD_DZY
TYPE(HALO2LIST_ll), POINTER :: TZHALO2_U, TZHALO2_V, TZHALO2_DZX, TZHALO2_DZY
INTEGER                     :: IINFO_ll
!JUAN
REAL          :: XPRECISION
!-----------------------------------------------------------------------
!
!*       1.    Compute the horizontal contravariant components
!              -----------------------------------------------
!
CALL MPPDB_CHECK3DM("contrav big ::PRU/V/WT",PRECISION,PRUT,PRVT,PRWT)                    
!
IIU= SIZE(PDXX,1)
IJU= SIZE(PDXX,2)
IKU= SIZE(PDXX,3)
!
CALL GET_INDICE_ll( IIB,IJB,IIE,IJE)
!
IKB=1+JPVEXT
IKE=IKU - JPVEXT
!
PRUCT(:,:,:) = PRUT(:,:,:) / PDXX(:,:,:)
PRVCT(:,:,:) = PRVT(:,:,:) / PDYY(:,:,:)
!
IF (KADV_ORDER == 4 ) THEN
 IF( .NOT. LFLAT) THEN 
  NULLIFY(TZFIELD_U)
  NULLIFY(TZFIELD_V)
  CALL ADD3DFIELD_ll(TZFIELD_U, PRUCT)
  CALL ADD3DFIELD_ll(TZFIELD_V, PRVCT)
  CALL UPDATE_HALO_ll(TZFIELD_U,IINFO_ll)
  CALL UPDATE_HALO_ll(TZFIELD_V,IINFO_ll)
!!$ IF( NHALO==1 ) THEN 
  NULLIFY(TZFIELD_DZX)
  NULLIFY(TZFIELD_DZY)
  CALL ADD3DFIELD_ll(TZFIELD_DZX, PDZX)
  CALL ADD3DFIELD_ll(TZFIELD_DZY, PDZY)
  NULLIFY(TZHALO2_U)
  NULLIFY(TZHALO2_V)
  NULLIFY(TZHALO2_DZX)
  NULLIFY(TZHALO2_DZY)
  CALL INIT_HALO2_ll(TZHALO2_U,1,IIU,IJU,IKU)
  CALL INIT_HALO2_ll(TZHALO2_V,1,IIU,IJU,IKU)
  CALL INIT_HALO2_ll(TZHALO2_DZX,1,IIU,IJU,IKU)
  CALL INIT_HALO2_ll(TZHALO2_DZY,1,IIU,IJU,IKU)
  CALL UPDATE_HALO2_ll(TZFIELD_U, TZHALO2_U, IINFO_ll)
  CALL UPDATE_HALO2_ll(TZFIELD_V, TZHALO2_V, IINFO_ll)
  CALL UPDATE_HALO2_ll(TZFIELD_DZX, TZHALO2_DZX, IINFO_ll)
  CALL UPDATE_HALO2_ll(TZFIELD_DZY, TZHALO2_DZY, IINFO_ll)
!!$ END IF
 END IF
END IF
!
!
!*       2.    Compute the vertical contravariant components (flat case)
!              ------------------------------------
!
IF (LFLAT) THEN
  PRWCT(:,:,:) = PRWT(:,:,:) / PDZZ(:,:,:)
  RETURN
END IF
!
!*       3.    Compute the vertical contravariant components (general case)
!              ------------------------------------
!
Z1 = 0.
Z2 = 0.
!
IF (KADV_ORDER == 2 ) THEN
!
  Z1(IIB:IIE,:,IKB:IKE+1)=                                             &
      (PRUCT(IIB:IIE,:,IKB:IKE+1)+PRUCT(IIB:IIE,:,IKB-1:IKE) )         &
       *PDZX(IIB:IIE,:,IKB:IKE+1) *0.25                                &
     +(PRUCT(IIB+1:IIE+1,:,IKB:IKE+1)+PRUCT(IIB+1:IIE+1,:,IKB-1:IKE) ) &
       *PDZX(IIB+1:IIE+1,:,IKB:IKE+1) *0.25   
                        
  Z2(:,IJB:IJE,IKB:IKE+1)=                                             &
      (PRVCT(:,IJB:IJE,IKB:IKE+1)+PRVCT(:,IJB:IJE,IKB-1:IKE) )         &
       *PDZY(:,IJB:IJE,IKB:IKE+1) *0.25                                &
     +(PRVCT(:,IJB+1:IJE+1,IKB:IKE+1)+PRVCT(:,IJB+1:IJE+1,IKB-1:IKE) ) &
       *PDZY(:,IJB+1:IJE+1,IKB:IKE+1) *0.25   
  PRWCT=0.             
  PRWCT(IIB:IIE,IJB:IJE,IKB:IKE+1) =                 &
      (   PRWT(IIB:IIE,IJB:IJE,IKB:IKE+1)            &
        - Z1(IIB:IIE,IJB:IJE,IKB:IKE+1)              &
        - Z2(IIB:IIE,IJB:IJE,IKB:IKE+1)              &
      ) / PDZZ(IIB:IIE,IJB:IJE,IKB:IKE+1)  
!
ELSE IF (KADV_ORDER == 4 ) THEN
!
!!$   IF (NHALO == 1) THEN
      IF ( LWEST_ll() .AND. HLBCX(1)/='CYCL' ) THEN
         IW=IIB+2 -1
      ELSE
         IW=IIB+1 -1
      END IF
      IE=IIE-1
!!$   ELSE
!!$      IF (LWEST_ll()) THEN
!!$         IW=IIB+1
!!$      ELSE
!!$         IW=IIB
!!$      END IF
!!$      IF (LEAST_ll() .AND. HLBCX(2)/='CYCL' ) THEN
!!$         IE=IIE-1
!!$      ELSE
!!$         IE=IIE
!!$      END IF
!!$   END IF
   !
!!$   IF(NHALO == 1) THEN
      IF ( LSOUTH_ll() .AND. HLBCY(1)/='CYCL' ) THEN
         IS=IJB+2 -1
      ELSE
         IS=IJB+1 -1
      END IF
      IN=IJE-1
!!$   ELSE
!!$      IF (LSOUTH_ll()) THEN
!!$         IS=IJB+1
!!$      ELSE
!!$         IS=IJB
!!$      END IF
!!$      IF (LNORTH_ll() .AND. HLBCY(2)/='CYCL' ) THEN
!!$         IN=IJE-1
!!$      ELSE
!!$         IN=IJE
!!$      END IF
!!$   END IF
   !
   !
   !*       3.1    interior of the processor subdomain
!
!
   Z1(IW:IE,:,IKB:IKE+1)=                                               &
       7.0*( (PRUCT(IW:IE,:,IKB:IKE+1)+PRUCT(IW:IE,:,IKB-1:IKE)) &
            *( 9.0*PDZX(IW:IE,:,IKB:IKE+1)-(PDZX(IW+1:IE+1,:,IKB:IKE+1) &
                  +PDZX(IW:IE,:,IKB:IKE+1)+PDZX(IW-1:IE-1,:,IKB:IKE+1))/3.0)/8.0 * 0.5 &
           +(PRUCT(IW+1:IE+1,:,IKB:IKE+1)+PRUCT(IW+1:IE+1,:,IKB-1:IKE)) &
            *( 9.0*PDZX(IW+1:IE+1,:,IKB:IKE+1)-(PDZX(IW+2:IE+2,:,IKB:IKE+1) &
                  +PDZX(IW+1:IE+1,:,IKB:IKE+1)+PDZX(IW:IE,:,IKB:IKE+1))/3.0)/8.0 * 0.5 )/12.0 &
          -( (PRUCT(IW-1:IE-1,:,IKB:IKE+1)+PRUCT(IW-1:IE-1,:,IKB-1:IKE)) &
            *PDZX(IW-1:IE-1,:,IKB:IKE+1) *0.5 &
            +(PRUCT(IW+2:IE+2,:,IKB:IKE+1)+PRUCT(IW+2:IE+2,:,IKB-1:IKE)) &
            *PDZX(IW+2:IE+2,:,IKB:IKE+1) *0.5)/12.0

!
   Z2(:,IS:IN,IKB:IKE+1)=                                               &
       7.0*( (PRVCT(:,IS:IN,IKB:IKE+1)+PRVCT(:,IS:IN,IKB-1:IKE)) &
            *( 9.0*PDZY(:,IS:IN,IKB:IKE+1)-(PDZY(:,IS+1:IN+1,IKB:IKE+1) &
                  +PDZY(:,IS:IN,IKB:IKE+1)+PDZY(:,IS-1:IN-1,IKB:IKE+1))/3.0)/8.0 * 0.5 &
           +(PRVCT(:,IS+1:IN+1,IKB:IKE+1)+PRVCT(:,IS+1:IN+1,IKB-1:IKE)) &
            *( 9.0*PDZY(:,IS+1:IN+1,IKB:IKE+1)-(PDZY(:,IS+2:IN+2,IKB:IKE+1) &
                  +PDZY(:,IS+1:IN+1,IKB:IKE+1)+PDZY(:,IS:IN,IKB:IKE+1))/3.0)/8.0 * 0.5 )/12.0 &
          -( (PRVCT(:,IS-1:IN-1,IKB:IKE+1)+PRVCT(:,IS-1:IN-1,IKB-1:IKE)) &
            *PDZY(:,IS-1:IN-1,IKB:IKE+1) *0.5 &
            +(PRVCT(:,IS+2:IN+2,IKB:IKE+1)+PRVCT(:,IS+2:IN+2,IKB-1:IKE)) &
            *PDZY(:,IS+2:IN+2,IKB:IKE+1) *0.5)/12.0
!
!*       3.2    limits of the processor subdomain (inside the whole domain or in cyclic conditions)
!
!!$  IF (NHALO==1) THEN

   Z1(IIE,:,IKB:IKE+1)=                                               &
       7.0*( (PRUCT(IIE,:,IKB:IKE+1)+PRUCT(IIE,:,IKB-1:IKE)) &
            *( 9.0*PDZX(IIE,:,IKB:IKE+1)-(PDZX(IIE+1,:,IKB:IKE+1) &
                  +PDZX(IIE,:,IKB:IKE+1)+PDZX(IIE-1,:,IKB:IKE+1))/3.0)/8.0 * 0.5 &
           +(PRUCT(IIE+1,:,IKB:IKE+1)+PRUCT(IIE+1,:,IKB-1:IKE)) &
            *( 9.0*PDZX(IIE+1,:,IKB:IKE+1)-(TZHALO2_DZX%HALO2%EAST(:,IKB:IKE+1) &
                  +PDZX(IIE+1,:,IKB:IKE+1)+PDZX(IIE,:,IKB:IKE+1))/3.0)/8.0 * 0.5 )/12.0 &
          -( (PRUCT(IIE-1,:,IKB:IKE+1)+PRUCT(IIE-1,:,IKB-1:IKE)) &
            *PDZX(IIE-1,:,IKB:IKE+1) *0.5 &
            +(TZHALO2_U%HALO2%EAST(:,IKB:IKE+1)+TZHALO2_U%HALO2%EAST(:,IKB-1:IKE)) &
            *TZHALO2_DZX%HALO2%EAST(:,IKB:IKE+1) *0.5)/12.0
!
   Z2(:,IJE,IKB:IKE+1)=                                               &
       7.0*( (PRVCT(:,IJE,IKB:IKE+1)+PRVCT(:,IJE,IKB-1:IKE)) &
            *( 9.0*PDZY(:,IJE,IKB:IKE+1)-(PDZY(:,IJE+1,IKB:IKE+1) &
                  +PDZY(:,IJE,IKB:IKE+1)+PDZY(:,IJE-1,IKB:IKE+1))/3.0)/8.0 * 0.5 &
           +(PRVCT(:,IJE+1,IKB:IKE+1)+PRVCT(:,IJE+1,IKB-1:IKE)) &
            *( 9.0*PDZY(:,IJE+1,IKB:IKE+1)-(TZHALO2_DZY%HALO2%NORTH(:,IKB:IKE+1) &
                  +PDZY(:,IJE+1,IKB:IKE+1)+PDZY(:,IJE,IKB:IKE+1))/3.0)/8.0 * 0.5 )/12.0 &
          -( (PRVCT(:,IJE-1,IKB:IKE+1)+PRVCT(:,IJE-1,IKB-1:IKE)) &
            *PDZY(:,IJE-1,IKB:IKE+1) *0.5 &
            +(TZHALO2_V%HALO2%NORTH(:,IKB:IKE+1)+TZHALO2_V%HALO2%NORTH(:,IKB-1:IKE)) &
            *TZHALO2_DZY%HALO2%NORTH(:,IKB:IKE+1) *0.5)/12.0
!!$  END IF
!
!*       3.3    non-CYCLIC CASE IN THE X DIRECTION: 2nd order case
!
  IF (HLBCX(1)/='CYCL' .AND. LWEST_ll()) THEN
!
   Z1(IIB,:,IKB:IKE+1)=                                     &
       (PRUCT(IIB,:,IKB:IKE+1)+PRUCT(IIB,:,IKB-1:IKE) )     &
        *PDZX(IIB,:,IKB:IKE+1) *0.25                        &
      +(PRUCT(IIB+1,:,IKB:IKE+1)+PRUCT(IIB+1,:,IKB-1:IKE) ) &
        *PDZX(IIB+1,:,IKB:IKE+1) *0.25   
  END IF
!
  IF (HLBCX(2)/='CYCL' .AND. LEAST_ll()) THEN
!
   Z1(IIE,:,IKB:IKE+1)=                                     &
       (PRUCT(IIE,:,IKB:IKE+1)+PRUCT(IIE,:,IKB-1:IKE) )     &
        *PDZX(IIE,:,IKB:IKE+1) *0.25                        &
      +(PRUCT(IIE+1,:,IKB:IKE+1)+PRUCT(IIE+1,:,IKB-1:IKE) ) &
        *PDZX(IIE+1,:,IKB:IKE+1) *0.25   
  END IF
!
!*       3.4    non-CYCLIC CASE IN THE Y DIRECTION: 2nd order case
!
  IF (HLBCY(1)/='CYCL' .AND. LSOUTH_ll()) THEN
!
   Z2(:,IJB,IKB:IKE+1)=                                     &
       (PRVCT(:,IJB,IKB:IKE+1)+PRVCT(:,IJB,IKB-1:IKE) )     &
        *PDZY(:,IJB,IKB:IKE+1) *0.25                        &
      +(PRVCT(:,IJB+1,IKB:IKE+1)+PRVCT(:,IJB+1,IKB-1:IKE) ) &
        *PDZY(:,IJB+1,IKB:IKE+1) *0.25   
!
  END IF
!
  IF (HLBCY(2)/='CYCL' .AND. LNORTH_ll()) THEN
!
   Z2(:,IJE,IKB:IKE+1)=                                     &
       (PRVCT(:,IJE,IKB:IKE+1)+PRVCT(:,IJE,IKB-1:IKE) )     &
        *PDZY(:,IJE,IKB:IKE+1) *0.25                        &
      +(PRVCT(:,IJE+1,IKB:IKE+1)+PRVCT(:,IJE+1,IKB-1:IKE) ) &
        *PDZY(:,IJE+1,IKB:IKE+1) *0.25   
!
  END IF
!
!*       3.5    Vertical contyravariant wind
!
!
!!$  CALL GET_HALO(Z1)
!!$  CALL GET_HALO(Z2)
!!$
!!$  CALL MPPDB_CHECK3DM("contrav ::Z1/Z2/ PDZZ",PRECISION,Z1,Z2,PDZZ)                    
  PRWCT=0.             
  PRWCT(IIB:IIE,IJB:IJE,IKB:IKE+1) =                 &
     (   PRWT(IIB:IIE,IJB:IJE,IKB:IKE+1)            &
       - Z1(IIB:IIE,IJB:IJE,IKB:IKE+1)              &
       - Z2(IIB:IIE,IJB:IJE,IKB:IKE+1)              &
     ) / PDZZ(IIB:IIE,IJB:IJE,IKB:IKE+1)  
!
END IF
!
PRWCT(:,:,1) = - PRWCT(:,:,3)     ! Mirror hypothesis
!
IF (KADV_ORDER == 4 ) THEN
  CALL CLEANLIST_ll(TZFIELD_U)
  CALL CLEANLIST_ll(TZFIELD_V)
!!$  IF (NHALO==1) THEN
    CALL CLEANLIST_ll(TZFIELD_DZX)
    CALL CLEANLIST_ll(TZFIELD_DZY)
    CALL DEL_HALO2_ll(TZHALO2_U)
    CALL DEL_HALO2_ll(TZHALO2_V)
    CALL DEL_HALO2_ll(TZHALO2_DZX)
    CALL DEL_HALO2_ll(TZHALO2_DZY)
!!$  END IF
END IF
!-----------------------------------------------------------------------
CALL MPPDB_CHECK3DM("contrav end ::PRU/V/WCT",PRECISION,PRUCT,PRVCT,PRWCT)                    
!
END SUBROUTINE CONTRAV
