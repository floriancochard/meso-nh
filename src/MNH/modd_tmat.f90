!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     #################
      MODULE MODD_TMAT
!     #################
!!
!!    PURPOSE
!!    -------
!!
!!      Initializations of variables for mode_tmat.f90
!!      (old varibales COMMON in FORTRAN 77)
!!
!!**  METHOD
!!    ------
!!
!!
!!    REFERENCE
!!    ---------
!!      None
!!
!!    AUTHOR
!!    ------
!!      G. TANGUY  *CNRM-MESONH*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    Original 23/03/2010
!!
!-------------------------------------------------------------------------------
!

IMPLICIT NONE





 

INTEGER,PARAMETER :: NPN1=200
INTEGER,PARAMETER :: NPNG1=600
INTEGER,PARAMETER :: NPNG2=2*NPNG1
INTEGER,PARAMETER :: NPN2=2*NPN1  
INTEGER,PARAMETER :: NPN4=NPN1
INTEGER,PARAMETER :: NPN6=NPN4+1
           

 
! COMMON /TMAT/ :dimensions : (NPN6,NPN4,NPN4)
REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: XRT11,XRT12
REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: XRT21,XRT22
REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: XIT11,XIT12
REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: XIT21,XIT22


!COMMON /CT/ dimensions : (NPN2,NPN2)
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,SAVE :: XTR1,XTI1


!COMMON /CTT/ dimensions : (NPN2,NPN2)
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,SAVE :: XQR,XQI,XRGQR,XRGQI

!
!COMMON /CBESS/ dimensions (NPNG2,NPN1)
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,SAVE ::XJ,XY,XJR,XJI,XDJ
REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE,SAVE ::XDJR,XDJI,XDY
             

END MODULE MODD_TMAT      

