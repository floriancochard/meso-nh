!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ######spl
      MODULE MODI_VER_SHIFT
!     #####################
INTERFACE
      FUNCTION VER_SHIFT(PVAR_MX,PZS_LS,PZS) RESULT(PVAR_SH)
!
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PVAR_MX
REAL,   DIMENSION(:,:),   INTENT(IN)  :: PZS_LS
REAL,   DIMENSION(:,:),   INTENT(IN)  :: PZS
REAL,   DIMENSION(SIZE(PVAR_MX,1),SIZE(PVAR_MX,2),SIZE(PVAR_MX,3)) :: &
                                         PVAR_SH
!
END FUNCTION VER_SHIFT
END INTERFACE
END MODULE MODI_VER_SHIFT
!     ######spl
      FUNCTION VER_SHIFT(PVAR_MX,PZS_LS,PZS) RESULT(PVAR_SH)
!     ########################################################
!
!!****  *VER_SHIFT* - vertical shift function
!!
!!    PURPOSE
!!    -------
!     This routine applies the shift function f to an array of 
!!    altitudes PVAR_MX.
!!    ---------
!!
!!**  METHOD
!!    ------
!! 1  The function f is first defined as:
!!
!!    f(z)= z - (zsA-zsNH) * (1.-  (tanh ((z-zsA)/zscale) )**2  )
!!
!!    zsA   :orography of the large scaale model; variable PZS_LS
!!    zsNH  :orography of MESO-NH               ; variable PZS
!!    zscale:typical height                     ; variable ZZSCALE
!!
!!    maximum of variation for: z = argsh(1/sqrt(2))*zscale = 0.6547*zscale
!!    We take zscale=2000m
!!
!! 2  But, in case of peaks in new orography:
!!       zscale=max(2000m, 32/(3sqrt3)*(zsNH-zsA) ) to assure f'(z)>0.875
!!       because f'(z)|extremum = 1 + 4*(zsA-zsNH)/zscale / (3sqrt3)
!!    In this case, the wind acceleration is only 25%
!!
!!    In order to keep the maximum of compression at the same height above
!!    the ground, a term of vertical displacement d is added:
!!
!!    f(z)= z - (zsA-zsNH) * (       1-(tanh ((z+d-zsA)/zscale))**2 
!!                               +C*(1-(tanh ((z  -zsA)/zscale))**2)    )
!!
!!     with          d = argsh(1/sqrt2) * ( zscale - 2000m )
!!     and therefore C = (tanh (d/zscale) )**2               < 0.330 and smaller
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    Book 2
!!
!!    AUTHOR
!!    ------
!!	
!     V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/12/94
!!                  14/05/96 (V. Masson) displacement of the function in case
!!                           of high change of orography
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PVAR_MX
REAL,   DIMENSION(:,:),   INTENT(IN)  :: PZS_LS
REAL,   DIMENSION(:,:),   INTENT(IN)  :: PZS
REAL,   DIMENSION(SIZE(PVAR_MX,1),SIZE(PVAR_MX,2),SIZE(PVAR_MX,3)) :: &
                                         PVAR_SH
!
!*       0.2   Declaration of local variables
!              ------------------------------
REAL                                                 :: ZDEFAULTZSCALE
REAL, DIMENSION(SIZE(PVAR_MX,1),SIZE(PVAR_MX,2)) :: ZZSCALE
REAL, DIMENSION(SIZE(PVAR_MX,1),SIZE(PVAR_MX,2)) :: ZD
REAL, DIMENSION(SIZE(PVAR_MX,1),SIZE(PVAR_MX,2)) :: ZC
INTEGER                                              :: IIU,IJU,IKU
!
INTEGER      :: JI,JJ,JK         ! loop index
!          
!-------------------------------------------------------------------------------
!
IIU=SIZE(PVAR_MX,1)
IJU=SIZE(PVAR_MX,2)
IKU=SIZE(PVAR_MX,3)
!
ZDEFAULTZSCALE=2000.
!
ZZSCALE(:,:)=MAX(ZDEFAULTZSCALE,6.158*(PZS(:,:)-PZS_LS(:,:)))
ZD(:,:)= 0.6547 * ( ZZSCALE(:,:) - ZDEFAULTZSCALE )
ZC(:,:)= (TANH (ZD(:,:)/ZZSCALE(:,:)) )**2
!
DO JK=1,IKU ; DO JJ=1,IJU ; DO JI=1,IIU
PVAR_SH(JI,JJ,JK)= PVAR_MX(JI,JJ,JK)                        &
                 - ( PZS_LS(JI,JJ)-PZS(JI,JJ) )             &
    * ( 1.0 - (TANH                                         &
      ((PVAR_MX(JI,JJ,JK) + ZD(JI,JJ)-PZS_LS(JI,JJ) )     &
      /ZZSCALE(JI,JJ)  )                                    &
                    )**2                                    &
       + ZC(JI,JJ)                                          &
       * (1.0 - (TANH                                       &
      ((PVAR_MX(JI,JJ,JK) - PZS_LS(JI,JJ) )                 &
      /ZZSCALE(JI,JJ)  )                                    & 
                    )**2 ) )
ENDDO ; ENDDO ; ENDDO

!-------------------------------------------------------------------------------
!
END FUNCTION VER_SHIFT
