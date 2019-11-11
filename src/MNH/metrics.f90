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
!     ###################
      MODULE MODI_METRICS
!     ###################
INTERFACE
!
SUBROUTINE METRICS(PMAP,PDXHAT,PDYHAT,PZZ,                            &
                   PDXX,PDYY,PDZX,PDZY,PDZZ)
REAL, DIMENSION(:,:),   INTENT(IN)  :: PMAP    ! Map factor
REAL, DIMENSION(:),     INTENT(IN)  :: PDXHAT  ! Stretching in x direction 
REAL, DIMENSION(:),     INTENT(IN)  :: PDYHAT  ! Stretching in y direction 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PZZ     ! Height in z direction
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PDXX  ! metric coefficient dxx
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PDYY  ! metric coefficient dyy 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PDZX  ! metric coefficient dzx 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PDZY  ! metric coefficient dzy 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PDZZ  ! metric coefficient dzz  
!
END SUBROUTINE METRICS
!
END INTERFACE
!
END MODULE MODI_METRICS
!
!
!
!     #################################################################
      SUBROUTINE METRICS(PMAP,PDXHAT,PDYHAT,PZZ,                      &
                         PDXX,PDYY,PDZX,PDZY,PDZZ)
!     #################################################################
!
!!****  *METRICS* - routine to compute metric coefficients
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to  compute the metric coefficients
!     dxx,dyy,dzz,dzx,dzy
!
!!**  METHOD
!!    ------
!!      The horizontal coefficients dxx and dyy (PDXX and PDYY arrays)
!!    are computed according to  the thinshell or no thinshell approximation
!!    and to the cartesian or spherical geometry.  
!!   
!!    EXTERNAL
!!    --------   
!!      MXM,MYM,MZM : Shuman functions (mean operators)
!!      DXM,DYM,DZM : Shuman functions (finite differences operators)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_CST : contains physical constants 
!!
!!        XRADIUS : earth radius 
!!
!!      Module MODD_CONF : contains configuration variables
!!
!!        LTHINSHELL : Logical for thinshell approximation
!!                     .TRUE. = Thinshell approximation done
!!        LCARTESIAN : Logical for cartesian geometry
!!                     .TRUE. = Cartesian geometry used
!!        
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine METRICS)
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    12/07/94 
!!                  14/02/01 (V. Masson and J. Stein) PDZZ initialized below the surface
!!                           (influences the 3D turbulence of W) and PDXX,PDYY,PDZZ at the top
!!                  19/03/2008 (J.Escobar) remove spread !!!
!!		    2014 (M.Faivre)
!!		    25/02/2015 (M.Moge) minor bug fix with MPPDB_CHECK
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!         
USE MODD_CONF
USE MODD_CST
!
USE MODI_SHUMAN
!
!20131024
USE MODE_MPPDB
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
REAL, DIMENSION(:,:),   INTENT(IN)  :: PMAP     ! Map factor
REAL, DIMENSION(:),     INTENT(IN)  :: PDXHAT   ! Stretching in x direction 
REAL, DIMENSION(:),     INTENT(IN)  :: PDYHAT   ! Stretching in y direction 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PZZ      ! Height (z) 
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PDXX  ! metric coefficient dxx
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PDYY  ! metric coefficient dyy 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PDZX  ! metric coefficient dzx 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PDZY  ! metric coefficient dzy 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PDZZ  ! metric coefficient dzz  
!
!*       0.2   declarations of local variables
!
INTEGER             :: IIU      ! Upper dimension in x direction
INTEGER             :: IJU      ! Upper dimension in y direction
INTEGER             :: IKU      ! Upper dimension in z direction
REAL                :: ZD1      ! DELTA1 (switch 0/1) for thinshell
                                ! approximation
INTEGER :: JI,JJ,JK
REAL, DIMENSION(SIZE(PDXHAT),SIZE(PDYHAT),SIZE(PZZ,3)) :: ZDZZ
!20131024
REAL, DIMENSION(SIZE(PDXHAT),SIZE(PDYHAT)) :: TEMP2D_PDXHAT
REAL, DIMENSION(SIZE(PDXHAT),SIZE(PDYHAT)) :: TEMP2D_PDYHAT
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE DIMENSIONS OF ARRAYS :
!              ----------------------------
IIU = SIZE(PDXHAT)
IJU = SIZE(PDYHAT)
IKU = SIZE(PZZ,3)
!
!-------------------------------------------------------------------------------
!
!*       2.   COMPUTE PDXX and PDYY  : 
!            --------------------
!
!20131024
CALL MPPDB_CHECK3D(PZZ,"METRICS::PZZ",PRECISION)
IF (.NOT.LCARTESIAN) THEN
  CALL MPPDB_CHECK2D(PMAP,"METRICS::PMAP",PRECISION)
ENDIF
!20131024
DO JI=1,IIU
TEMP2D_PDXHAT(JI,:) = PDXHAT(JI)
END DO
DO JJ=1,IJU
TEMP2D_PDYHAT(:,JJ) = PDYHAT(JJ)
END DO
CALL MPPDB_CHECK2D(TEMP2D_PDXHAT,"METRICS::PDXHAT",PRECISION)
CALL MPPDB_CHECK2D(TEMP2D_PDYHAT,"METRICS::PDYHAT",PRECISION)
!
IF (LTHINSHELL) THEN
  ZD1=0.
ELSE
  ZD1=1.
END IF
IF (.NOT.LCARTESIAN) THEN
  ZDZZ(:,:,:) = MZF(1,IKU,1, 1.+ ZD1*PZZ(:,:,:)/XRADIUS)
  DO JK=1,IKU ; DO JJ=1,IJU ; DO JI=1,IIU
    PDXX(JI,JJ,JK) = ZDZZ(JI,JJ,JK) * PDXHAT(JI) /PMAP(JI,JJ)
    PDYY(JI,JJ,JK) = ZDZZ(JI,JJ,JK) * PDYHAT(JJ) /PMAP(JI,JJ)
  ENDDO ; ENDDO ; ENDDO
  !20140710
  CALL MPPDB_CHECK3D(PDXX,"METRICSbefMXM::PDXX",PRECISION)
  CALL MPPDB_CHECK3D(PDYY,"METRICSbefMYM::PDYY",PRECISION)
  PDXX(:,:,:)=MXM(PDXX(:,:,:))
  PDXX(:,:,IKU)=PDXX(:,:,IKU-1)
  PDYY(:,:,:)=MYM(PDYY(:,:,:))
  PDYY(:,:,IKU)=PDYY(:,:,IKU-1)
ELSE
  DO JK=1,IKU ; DO JJ=1,IJU ; DO JI=1,IIU
    PDXX(JI,JJ,JK) = PDXHAT(JI)
    PDYY(JI,JJ,JK) = PDYHAT(JJ)
  ENDDO ; ENDDO ; ENDDO
  PDXX(:,:,:)=MXM(PDXX(:,:,:))
  PDYY(:,:,:)=MYM(PDYY(:,:,:))
END IF
!
!20131024
CALL MPPDB_CHECK3D(PDXX,"METRICSaftMXM::PDXX",PRECISION)
CALL MPPDB_CHECK3D(PDYY,"METRICSaftMYM::PDYY",PRECISION)
!
!-------------------------------------------------------------------------------
!
!*       3.  COMPUTE PDZX AND PDZY  :
!            ----------------------
!
PDZX(:,:,:) = DXM(PZZ(:,:,:))
! 
PDZY(:,:,:) = DYM(PZZ(:,:,:))
!
!-------------------------------------------------------------------------------
!
!*       4.  COMPUTE PDZZ  :
!            -------------
!
PDZZ(:,:,:) = DZM(1,IKU,1,MZF(1,IKU,1,PZZ(:,:,:)))
PDZZ(:,:,IKU) = PZZ(:,:,IKU) - PZZ(:,:,IKU-1)  ! same delta z in IKU and IKU -1
PDZZ(:,:,1)   = PDZZ(:,:,2)                    ! same delta z in 1   and 2
!20131024
CALL MPPDB_CHECK3D(PDZZ,"METRICS::PDZZ",PRECISION)
!-----------------------------------------------------------------------------
!
END SUBROUTINE METRICS
