!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 prep_real 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_TRUNC_FIELD
!     #######################

INTERFACE TRUNC_FIELD

      SUBROUTINE TRUNC_FIELD2D(KIU,KJU,KIINF,KISUP,KJINF,KJSUP,PF,PF_LS)

INTEGER,           INTENT(IN)  :: KIU        ! size of new domain in x direction
INTEGER,           INTENT(IN)  :: KJU        ! size of new domain in y direction
INTEGER,           INTENT(IN)  :: KIINF      !
INTEGER,           INTENT(IN)  :: KISUP      !
INTEGER,           INTENT(IN)  :: KJINF      ! zoom
INTEGER,           INTENT(IN)  :: KJSUP      !
REAL,DIMENSION(:,:), INTENT(IN):: PF         ! large-scale field on atm. domain
REAL,DIMENSION(:,:), INTENT(OUT)::PF_LS      ! large-scale field on new domain

END SUBROUTINE TRUNC_FIELD2D

      SUBROUTINE TRUNC_FIELD3D(KIU,KJU,KIINF,KISUP,KJINF,KJSUP,PF,PF_LS)

INTEGER,           INTENT(IN)  :: KIU        ! size of new domain in x direction
INTEGER,           INTENT(IN)  :: KJU        ! size of new domain in y direction
INTEGER,           INTENT(IN)  :: KIINF      !
INTEGER,           INTENT(IN)  :: KISUP      !
INTEGER,           INTENT(IN)  :: KJINF      ! zoom
INTEGER,           INTENT(IN)  :: KJSUP      !
REAL,DIMENSION(:,:,:), INTENT(IN):: PF       ! large-scale field on atm. domain
REAL,DIMENSION(:,:,:), INTENT(OUT)::PF_LS    ! large-scale field on new domain

END SUBROUTINE TRUNC_FIELD3D

      SUBROUTINE TRUNC_FIELD4D(KIU,KJU,KIINF,KISUP,KJINF,KJSUP,PF,PF_LS)

INTEGER,           INTENT(IN)  :: KIU        ! size of new domain in x direction
INTEGER,           INTENT(IN)  :: KJU        ! size of new domain in y direction
INTEGER,           INTENT(IN)  :: KIINF      !
INTEGER,           INTENT(IN)  :: KISUP      !
INTEGER,           INTENT(IN)  :: KJINF      ! zoom
INTEGER,           INTENT(IN)  :: KJSUP      !
REAL,DIMENSION(:,:,:,:), INTENT(IN):: PF     ! large-scale field on atm. domain
REAL,DIMENSION(:,:,:,:), INTENT(OUT)::PF_LS  ! large-scale field on new domain

END SUBROUTINE TRUNC_FIELD4D

END INTERFACE

END MODULE MODI_TRUNC_FIELD

!     ##################################################################
      SUBROUTINE TRUNC_FIELD2D(KIU,KJU,KIINF,KISUP,KJINF,KJSUP,PF,PF_LS)
!     ##################################################################

INTEGER,           INTENT(IN)  :: KIU        ! size of new domain in x direction
INTEGER,           INTENT(IN)  :: KJU        ! size of new domain in y direction
INTEGER,           INTENT(IN)  :: KIINF      !
INTEGER,           INTENT(IN)  :: KISUP      !
INTEGER,           INTENT(IN)  :: KJINF      ! zoom
INTEGER,           INTENT(IN)  :: KJSUP      !
REAL,DIMENSION(:,:), INTENT(IN):: PF         ! large-scale field on atm. domain
REAL,DIMENSION(:,:), INTENT(OUT)::PF_LS      ! large-scale field on new domain


IF (SIZE(PF)==0) THEN
  RETURN
END IF

PF_LS(:,:) = PF(KIINF:KISUP,KJINF:KJSUP)

END SUBROUTINE TRUNC_FIELD2D

!-------------------------------------------------------------------------------
!
!     ##################################################################
      SUBROUTINE TRUNC_FIELD3D(KIU,KJU,KIINF,KISUP,KJINF,KJSUP,PF,PF_LS)
!     ##################################################################

INTEGER,           INTENT(IN)  :: KIU        ! size of new domain in x direction
INTEGER,           INTENT(IN)  :: KJU        ! size of new domain in y direction
INTEGER,           INTENT(IN)  :: KIINF      !
INTEGER,           INTENT(IN)  :: KISUP      !
INTEGER,           INTENT(IN)  :: KJINF      ! zoom
INTEGER,           INTENT(IN)  :: KJSUP      !
REAL,DIMENSION(:,:,:), INTENT(IN):: PF       ! large-scale field on atm. domain
REAL,DIMENSION(:,:,:), INTENT(OUT)::PF_LS    ! large-scale field on new domain


IF (SIZE(PF)==0) THEN
  RETURN
END IF

PF_LS(:,:,:) = PF(KIINF:KISUP,KJINF:KJSUP,:)

END SUBROUTINE TRUNC_FIELD3D

!-------------------------------------------------------------------------------
!
!     ##################################################################
      SUBROUTINE TRUNC_FIELD4D(KIU,KJU,KIINF,KISUP,KJINF,KJSUP,PF,PF_LS)
!     ##################################################################

INTEGER,           INTENT(IN)  :: KIU        ! size of new domain in x direction
INTEGER,           INTENT(IN)  :: KJU        ! size of new domain in y direction
INTEGER,           INTENT(IN)  :: KIINF      !
INTEGER,           INTENT(IN)  :: KISUP      !
INTEGER,           INTENT(IN)  :: KJINF      ! zoom
INTEGER,           INTENT(IN)  :: KJSUP      !
REAL,DIMENSION(:,:,:,:), INTENT(IN):: PF     ! large-scale field on atm. domain
REAL,DIMENSION(:,:,:,:), INTENT(OUT)::PF_LS  ! large-scale field on new domain


IF (SIZE(PF)==0) THEN
  RETURN
END IF

PF_LS(:,:,:,:) = PF(KIINF:KISUP,KJINF:KJSUP,:,:)

END SUBROUTINE TRUNC_FIELD4D
