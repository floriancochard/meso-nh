!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 newsrc 2006/11/16 16:45:25
!-----------------------------------------------------------------
!!   ##############################
     MODULE MODI_SALT_FILTER
!!   ##############################
!!
INTERFACE
!
SUBROUTINE SALT_FILTER(PSV, PRHODREF)

IMPLICIT NONE

REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSV
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PRHODREF

END SUBROUTINE SALT_FILTER
!!
END INTERFACE
!!
END MODULE MODI_SALT_FILTER
!!
!!   #######################################
     SUBROUTINE SALT_FILTER(PSV, PRHODREF)
!!   #######################################
!!
!!   PURPOSE
!!   -------
!!
!!   REFERENCE
!!   ---------
!!   none
!!
!!   AUTHOR
!!    ------
!!    Pierre TULET (CNRM/GMEI) 
!!
!!   MODIFICATIONS
!!    -------------
!!   Original
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!   IMPLICIT ARGUMENTS
!
USE MODD_SALT
USE MODE_SALT_PSD

!+ Marine
USE MODI_INIT_SALT
!- Marine

!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSV
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PRHODREF
!
!*      0.2    declarations local variables
!
REAL                                  :: ZRHOI               ! [kg/m3] density of aerosol
REAL                                  :: ZMI                 ! [kg/mol] molar weight of aerosol
REAL                                  :: ZRGMIN              ! [um] minimum radius accepted
REAL                                  :: ZSIGMIN             ! minimum standard deviation accepted
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZM                  ! [aerosol units] local array which goes to output later
REAL,DIMENSION(:),       ALLOCATABLE  :: ZMMIN               ! [aerosol units] minimum values for N, sigma, M
INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM0                 ! [idx] index for Mode 0 in passed variables
INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM3                 ! [idx] indexes for Mode 3 in passed variables
INTEGER,DIMENSION(:),    ALLOCATABLE  :: NM6                 ! [idx] indexes for Mode 6 in passed variables
REAL,DIMENSION(:),       ALLOCATABLE  :: ZINIRADIUS          ! initial mean radius
REAL,DIMENSION(:),       ALLOCATABLE  :: ZINISIGMA           ! initial standard deviation
INTEGER                               :: JN,IMODEIDX,JJ      ! [idx] loop counters

!-------------------------------------------------------------------------------
!+ Marine
CALL INIT_SALT
!- Marine


ALLOCATE (NM0(NMODE_SLT))
ALLOCATE (NM3(NMODE_SLT))
ALLOCATE (NM6(NMODE_SLT))
ALLOCATE (ZM(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3), NMODE_SLT*3))
ALLOCATE (ZMMIN(NMODE_SLT*3))
ALLOCATE (ZINIRADIUS(NMODE_SLT))
ALLOCATE (ZINISIGMA(NMODE_SLT))

PSV(:,:,:,:) = MAX(PSV(:,:,:,:), XMNH_TINY)

DO JN=1,NMODE_SLT
  IMODEIDX = JPSALTORDER(JN)
  !Calculations here are for one mode only
   ZINISIGMA(JN) = XINISIG_SLT(IMODEIDX)
  IF (CRGUNITS=="MASS") THEN
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX) * EXP(-3.*(LOG(XINISIG_SLT(IMODEIDX)))**2)
  ELSE
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX)
  END IF

  !Set counter for number, M3 and M6
  NM0(JN) = 1+(JN-1)*3
  NM3(JN) = 2+(JN-1)*3
  NM6(JN) = 3+(JN-1)*3
  !Get minimum values possible
  ZMMIN(NM0(JN)) = XN0MIN_SLT(IMODEIDX)
  ZRGMIN         = ZINIRADIUS(JN)
  IF (LVARSIG_SLT) THEN
    ZSIGMIN = XSIGMIN_SLT
  ELSE
    ZSIGMIN = XINISIG_SLT(IMODEIDX)
  ENDIF
  ZMMIN(NM3(JN)) = ZMMIN(NM0(JN)) * (ZRGMIN**3)*EXP(4.5 * LOG(ZSIGMIN)**2)
  ZMMIN(NM6(JN)) = ZMMIN(NM0(JN)) * (ZRGMIN**6)*EXP(18. * LOG(ZSIGMIN)**2)
END DO
!
!Set density of aerosol, here dust (kg/m3)
ZRHOI = XDENSITY_SALT
!Set molecular weight of dust !NOTE THAT THIS IS NOW IN KG
ZMI   = XMOLARWEIGHT_SALT
!
DO JN=1,NMODE_SLT
     !calculate moment 3 from total aerosol mass (molec_{aer}/molec_{air} ==>
     IF (LVARSIG_SLT) THEN
     ZM(:,:,:,NM3(JN)) =                        &
          PSV(:,:,:,2+(JN-1)*3)                & !molec_{aer}/molec_{aer}
          * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
          * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
          * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
          * XM3TOUM3_SALT                            & !==>um3_{aer}/m3_{air}
          / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
    ELSE
       IF ((LRGFIX_SLT)) THEN
          ZM(:,:,:,NM3(JN)) =                   &
          PSV(:,:,:,JN)                        & !molec_{aer}/molec_{aer}
          * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
          * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
          * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
          * XM3TOUM3_SALT                            & !==>um3_{aer}/m3_{air}
          / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
      ELSE
          ZM(:,:,:,NM3(JN)) =                   &
          PSV(:,:,:,2+(JN-1)*2)                & !molec_{aer}/molec_{aer}
          * (ZMI/XMD)                           & !==>kg_{aer}/kg_{aer}
          * PRHODREF(:,:,:)                     & !==>kg_{aer}/m3_{air}
          * (1.d0/ZRHOI)                        & !==>m3_{aer}/m3_{air}
          * XM3TOUM3_SALT                            & !==>um3_{aer}/m3_{air}
          / (XPI * 4./3.)                         !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
      END IF
    END IF

! calculate moment 0 from dispersion and mean radius
     ZM(:,:,:,NM0(JN))=  ZM(:,:,:,NM3(JN))/&
       ((ZINIRADIUS(JN)**3)*EXP(4.5 * LOG(ZINISIGMA(JN))**2))


! calculate moment 6 from dispersion and mean radius
     ZM(:,:,:,NM6(JN)) = ZM(:,:,:,NM0(JN)) * (ZINIRADIUS(JN)**6) * &
               EXP(18 *(LOG(ZINISIGMA(JN)))**2)

     IF (LVARSIG_SLT) THEN
     WHERE ((ZM(:,:,:,NM0(JN)) .LT. ZMMIN(NM0(JN))).OR.&
            (ZM(:,:,:,NM3(JN)) .LT. ZMMIN(NM3(JN))).OR.&
            (ZM(:,:,:,NM6(JN)) .LT. ZMMIN(NM6(JN))))
     ZM(:,:,:,NM0(JN)) = ZMMIN(NM0(JN))
     ZM(:,:,:,NM3(JN)) = ZMMIN(NM3(JN))
     ZM(:,:,:,NM6(JN)) = ZMMIN(NM6(JN))
     END WHERE

     ELSE IF (.NOT.(LRGFIX_SLT)) THEN

     WHERE ((ZM(:,:,:,NM0(JN)) .LT. ZMMIN(NM0(JN))).OR.&
            (ZM(:,:,:,NM3(JN)) .LT. ZMMIN(NM3(JN))))
     ZM(:,:,:,NM0(JN)) = ZMMIN(NM0(JN))
     ZM(:,:,:,NM3(JN)) = ZMMIN(NM3(JN))
     END WHERE
     ENDIF

    ! return to concentration #/m3 =>  (#/molec_{air}
     IF (LVARSIG_SLT) THEN
     PSV(:,:,:,1+(JN-1)*3) = ZM(:,:,:,NM0(JN)) * XMD / &
                              (XAVOGADRO*PRHODREF(:,:,:))

     PSV(:,:,:,2+(JN-1)*3) = ZM(:,:,:,NM3(JN)) * XMD  * XPI * 4./3 * ZRHOI / &
                              (ZMI*PRHODREF(:,:,:)*XM3TOUM3_SALT)

     PSV(:,:,:,3+(JN-1)*3) = ZM(:,:,:,NM6(JN)) * XMD  / &
                              ( XAVOGADRO*PRHODREF(:,:,:) * 1.d-6)
     ELSE IF (LRGFIX_SLT) THEN
      PSV(:,:,:,JN)        = ZM(:,:,:,NM3(JN)) * XMD * XPI * 4./3. * ZRHOI  / &
                              (ZMI*PRHODREF(:,:,:)*XM3TOUM3_SALT)
     ELSE
     PSV(:,:,:,1+(JN-1)*2) = ZM(:,:,:,NM0(JN)) * XMD / &
                              (XAVOGADRO*PRHODREF(:,:,:))

     PSV(:,:,:,2+(JN-1)*2) = ZM(:,:,:,NM3(JN)) * XMD * XPI * 4./3. * ZRHOI  / &
                              (ZMI*PRHODREF(:,:,:)*XM3TOUM3_SALT)
     END IF

!
END DO  !Loop on modes

DEALLOCATE(ZINIRADIUS)
DEALLOCATE(ZMMIN)
DEALLOCATE(ZINISIGMA)
DEALLOCATE(ZM)
DEALLOCATE(NM6)
DEALLOCATE(NM3)
DEALLOCATE(NM0)


END SUBROUTINE SALT_FILTER
