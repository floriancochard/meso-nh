!     ######spl
      MODULE  MODD_CTL_AXES_AND_STYL
!     ##############################
!
!!****  *MODD_CTL_AXES_AND_STYL* - 
!!
!!    PURPOSE
!!    -------
!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     None
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      JD    "LA"
!!
!!
!!    MODIFICATIONS
!!    -------------
!!
!!     original        01/02/96
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!

IMPLICIT NONE


! Controle des graduations majeures (labellees) et mineures sur les axes X et Y
! On donne le nb d'intervalles
!******************************************************************************
! CH Projection cartographique _K_  _Z_  _PR_  _TK_  _EV_
INTEGER:: NCHPCITVXMJ=0,NCHPCITVYMJ=0,NCHPCITVXMN=0,NCHPCITVYMN=0
! CH Cartesien _K_  _Z_  _PR_  _TK_ 
INTEGER:: NCHITVXMJ=5,NCHITVYMJ=4,NCHITVXMN=1,NCHITVYMN=1
! CV  _CV_ et _PVT_
INTEGER:: NCVITVXMJ=5,NCVITVYMJ=10,NCVITVXMN=1,NCVITVYMN=1
! PV  _PV_
INTEGER:: NPVITVXMJ=4,NPVITVYMJ=0,NPVITVXMN=1,NPVITVYMN=1
! FT  _FT_  _PVKT_
INTEGER:: NFTITVXMJ=5,NFTITVYMJ=5,NFTITVXMN=2,NFTITVYMN=2
! FT1  _FT1_
INTEGER:: NFT1ITVXMJ=5,NFT1ITVYMJ=5,NFT1ITVXMN=1,NFT1ITVYMN=1
! XY  _XY_
INTEGER:: NXYITVXMJ=5,NXYITVYMJ=5,NXYITVXMN=1,NXYITVYMN=1
! MASK  _MASK_
INTEGER:: NMASKITVXMJ=5,NMASKITVYMJ=5,NMASKITVXMN=1,NMASKITVYMN=1

! Axes labelles en latitude, longitude pour CH Proj. cart.
LOGICAL :: LGEOG=.FALSE.

! Axes labelles en indices de grilles
LOGICAL,SAVE :: LINDAX=.FALSE.

! Gestion de la taille titres en X
REAL :: XSZTITXL=0., XSZTITXM=0., XSZTITXR=0.

! Controle du type de trait avec _FT1_
LOGICAL :: LFT1STYLUSER=.FALSE.
LOGICAL :: LFTSTYLUSER=.FALSE.
LOGICAL :: LTITFTUSER=.FALSE.
!
!
END MODULE MODD_CTL_AXES_AND_STYL
