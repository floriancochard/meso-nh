!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
!     ######################  
      MODULE MODI_SET_DIRCOS
!     ######################
INTERFACE
!
      SUBROUTINE SET_DIRCOS(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,                 &  
                            TPINITHALO2D_ll,PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,   &
                            PCOSSLOPE,PSINSLOPE                          )
!
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) ::  HLBCX ! X-direction LBC type at left(1)
                                             ! and right(2) boundaries  
CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) ::  HLBCY  ! Y-direction LBC type at left(1)
                                             ! and right(2) boundaries
REAL, DIMENSION(:,:,:),   INTENT(IN)         :: PDXX,PDYY,PDZX,PDZY
                                             ! metric coefficients
TYPE(LIST_ll), POINTER :: TPINITHALO2D_ll ! pointer for the list 
                                      ! of fields that have  to be communicated
REAL, DIMENSION(:,:),   INTENT(OUT)     ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW 
! Director Cosinus in x, y and z directions at the surface w-point
REAL, DIMENSION(:,:),   INTENT(OUT)     ::  PCOSSLOPE       ! cosinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(:,:),   INTENT(OUT)     ::  PSINSLOPE       ! sinus of the angle
                                 ! between i and the slope vector
!
END SUBROUTINE SET_DIRCOS
!
END INTERFACE
!
END MODULE MODI_SET_DIRCOS
!
!
!
!     ########################################################################
      SUBROUTINE SET_DIRCOS(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,                 &  
                            TPINITHALO2D_ll,PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,   &
                            PCOSSLOPE,PSINSLOPE                          )
!     ########################################################################
!
!
!!****  *SET_DIRCOS* - computes the director Cosinus at the surface
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to compute the director cosinus 
!!    at the surface
!
!!**  METHOD
!!    ------
!!      The following steps are made:
!!      1- The metric coefficients are recovered from the grid knowledge.
!!      2- The director Cosinus, located at the surface w-point, are
!!         computed according to the following formulae:
!!      
!!                ______x      _______x     _______y
!!      -   -     (d_zx)       ( d_zx )^2   ( d_zy ) ^2
!!      n . i = - (----)* (1 + (------)   + (------)      ) ^(-1/2)
!!                (___z)       ( ___z )     ( ___z )
!!                (d_xx)       ( d_xx )     ( d_yy )
!!
!!      
!!                ______y      _______x     _______y
!!      -   -     (d_zy)       ( d_zx )^2   ( d_zy ) ^2
!!      n . j = - (----)* (1 + (------)   + (------)      ) ^(-1/2)
!!                (___z)       ( ___z )     ( ___z )
!!                (d_yy)       ( d_xx )     ( d_yy )
!!
!!      
!!                             _______x     _______y
!!      -   -                  ( d_zx )^2   ( d_zy ) ^2
!!      n . k = +         (1 + (------)   + (------)      ) ^(-1/2)
!!                             ( ___z )     ( ___z )
!!                             ( d_xx )     ( d_yy )
!!                                                                -
!!     3- compute the Cosinus and the sinus of the angle between  i and the
!!     slope vector defined by its components in the base (i,j):
!!       
!!       ----->     (             )
!!       slope    = (dzs/dx,dzs/dy)
!!
!!
!!    EXTERNAL
!!    --------
!!      SUBROUTINE MXF,MYF   : computes the averages along the x and y
!!                  direction for a variable located at the Flux point
!!
!!      Module MODI_SHUMAN           : interface module for the Shuman
!!                                     operators
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!       MODD_PARAMETERS : JPVEXT  number of marginal vertical points
!!                         
!!
!!    REFERENCE
!!    ---------
!!      Book 2 of documentation (routine SET_DIRCOS)
!!      Book 1 of documentation (Chapter Turbulence)
!!
!!    AUTHOR
!!    ------
!!      Joan Cuxart             * INM and Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         14/02/95
!!      (J.Stein)        15/11/95  add the slope angle
!!      V. DUCROCQ       14/08/98  //
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
USE MODD_PARAMETERS
!
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_CST, ONLY : XMNH_TINY
!
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
!
CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) ::  HLBCX ! X-direction LBC type at left(1)
                                             ! and right(2) boundaries  
CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) ::  HLBCY  ! Y-direction LBC type at left(1)
                                             ! and right(2) boundaries
REAL, DIMENSION(:,:,:),   INTENT(IN)         :: PDXX,PDYY,PDZX,PDZY
                                             ! metric coefficients
TYPE(LIST_ll), POINTER :: TPINITHALO2D_ll ! pointer for the list 
                                      ! of fields that have  to be communicated
REAL, DIMENSION(:,:),   INTENT(OUT)     ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW 
! Director Cosinus in x, y and z directions at the surface w-point
REAL, DIMENSION(:,:),   INTENT(OUT)     ::  PCOSSLOPE       ! cosinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(:,:),   INTENT(OUT)     ::  PSINSLOPE       ! sinus of the angle
                                 ! between i and the slope vector
!
!
!       0.2  declaration of local variables

!
REAL, DIMENSION(SIZE(PDXX,1),SIZE(PDXX,2),1+JPVEXT:1+JPVEXT) :: ZWORK1, ZWORK2 
! pseudo-3D array allowing matricial computation with the Shuman operators
!
INTEGER     :: IKB,IIB,IIE,IJB,IJE
                        ! index values for the Beginning or the End of the physical
!
!----------------------------------------------------------------------------
!
!*      1. COMPUTE THE DIRECTOR COSINUS
!          ----------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT              
!
ZWORK1(:,:,IKB:IKB) =      &
   MXF( PDZX(:,:,IKB:IKB) / (0.5* (PDXX(:,:,IKB-1:IKB-1) + PDXX(:,:,IKB:IKB))) )
!
IF ( HLBCX(1) /= "CYCL" .AND. LWEST_ll()) THEN
  ZWORK1(IIB-1,:,IKB:IKB) =  ZWORK1(IIB,:,IKB:IKB)
END IF
IF ( HLBCX(2) /= "CYCL" .AND. LEAST_ll()) THEN
  ZWORK1(IIE+1,:,IKB:IKB) =  ZWORK1(IIE,:,IKB:IKB)
END IF
!
ZWORK2(:,:,IKB:IKB) =      &
   MYF( PDZY(:,:,IKB:IKB) / (0.5* (PDYY(:,:,IKB-1:IKB-1) + PDYY(:,:,IKB:IKB))) )
!
IF ( HLBCY(1) /= "CYCL" .AND. LSOUTH_ll()) THEN
  ZWORK2(:,IJB-1,IKB:IKB) = ZWORK2(:,IJB,IKB:IKB)
END IF
IF ( HLBCY(2) /= "CYCL" .AND. LNORTH_ll()) THEN
  ZWORK2(:,IJE+1,IKB:IKB) = ZWORK2(:,IJE,IKB:IKB)
END IF
!
PDIRCOSZW(:,:) = 1. / SQRT     &
 ( 1. + ZWORK1(:,:,IKB)**2 + ZWORK2(:,:,IKB)**2 )
!
PDIRCOSXW(:,:) = - ZWORK1(:,:,IKB) * PDIRCOSZW(:,:)
!
PDIRCOSYW(:,:) = - ZWORK2(:,:,IKB) * PDIRCOSZW(:,:)
!
CALL ADD2DFIELD_ll(TPINITHALO2D_ll,PDIRCOSZW)
CALL ADD2DFIELD_ll(TPINITHALO2D_ll,PDIRCOSXW)
CALL ADD2DFIELD_ll(TPINITHALO2D_ll,PDIRCOSYW)
!
!*      2. COMPUTE THE SLOPE COSINUS AND SINUS
!          -----------------------------------
!
WHERE (  ZWORK1(:,:,IKB)**2 + ZWORK2(:,:,IKB)**2 < XMNH_TINY )
  PCOSSLOPE(:,:) = 1.
  !
  PSINSLOPE(:,:) = 0.
ELSEWHERE
  PCOSSLOPE(:,:) = ZWORK1(:,:,IKB) /                &
    SQRT ( ZWORK1(:,:,IKB)**2 + ZWORK2(:,:,IKB)**2 )
  !
  PSINSLOPE(:,:) = ZWORK2(:,:,IKB) /                &
    SQRT ( ZWORK1(:,:,IKB)**2 + ZWORK2(:,:,IKB)**2 )
END WHERE
CALL ADD2DFIELD_ll(TPINITHALO2D_ll,PCOSSLOPE)
CALL ADD2DFIELD_ll(TPINITHALO2D_ll,PSINSLOPE)
!
!
!------------------------------------------------------------------------------
!
END SUBROUTINE SET_DIRCOS
