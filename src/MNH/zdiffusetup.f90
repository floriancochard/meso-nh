!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 newsrc 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_ZDIFFUSETUP
!     ####################
!
INTERFACE
!
!     ######################################################################
      SUBROUTINE ZDIFFUSETUP (PZZ,&
                              PZDIFFU_HALO2)
!     ######################################################################

!JUAN
USE MODE_TYPE_ZDIFFU
IMPLICIT NONE
!JUAN 
REAL, DIMENSION(:,:,:), INTENT(IN)        :: PZZ       ! Height
!JUAN
TYPE(TYPE_ZDIFFU_HALO2)                       :: PZDIFFU_HALO2
!JUAN
END SUBROUTINE ZDIFFUSETUP
!
END INTERFACE
!
END MODULE MODI_ZDIFFUSETUP
!
!
!
!     #########################################################################
      SUBROUTINE ZDIFFUSETUP (PZZ,&
                              PZDIFFU_HALO2)
!     #########################################################################
!
!!****  *ZDIFFUSETUP* - routine to calculate the interpolation coefficient needed
!!                   for the truly horizontal diffusion scheme
!!
!!    REFERENCE
!!    ---------
!!
!!      Zängl, G., 2002: An improved method for computing horizontal diffusion in a
!!                       sigma-coordinate model and its application to simulations
!!                       over mountainous topography. Mon. Wea. Rev. 130, 1423-1432.
!!
!!    AUTHOR
!!    ------
!!
!!      G. Zängl       * University of Munich*
!!      J.Escobar 7/10/2015 : remove print
!!
!
!*       0.    DECLARATIONS
!              ------------ 
USE MODD_CST
USE MODD_PARAMETERS
USE MODD_CONF
USE MODI_RELAX
USE MODE_ll
USE MODD_VAR_ll, ONLY : NMNH_COMM_WORLD
USE MODI_SHUMAN
USE MODD_ARGSLIST_ll, ONLY : LIST_ll, HALO2LIST_ll
!
!
!JUAN
USE MODE_TYPE_ZDIFFU
USE MODE_SUM_LL
!JUAN 
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!

REAL, DIMENSION(:,:,:), INTENT(IN)        :: PZZ       ! Height
!JUAN
TYPE(TYPE_ZDIFFU_HALO2)                       :: PZDIFFU_HALO2
!JUAN

! local variables

INTEGER                    :: IIB_ll,IJB_ll ! Lower bounds of the physical
                                            !  global domain in x and y directions
INTEGER                    :: IIE_ll,IJE_ll ! Upper bounds of the physical
                                            ! gloBal domain in x and y directions
                                            ! sub-domain:
INTEGER                    :: IIB,IJB,IKB   ! Lower bounds of the physical
                                            ! sub-domain in x,y and z directions
INTEGER                    :: IIE,IJE,IKE   ! Upper bounds of the physical
                                            ! sub-domain in x,y and z directions
INTEGER                    :: IIU,IJU,IKU   ! domain sizes
INTEGER                    :: JI,JJ,JK   ! Loop indexes
INTEGER                    :: IIMAX_ll,IJMAX_ll ! Number of points of
                                                ! Global physical domain
                                                ! in the x and y directions
INTEGER, DIMENSION(:,:),   ALLOCATABLE :: IKIP1,IKIP2,IKIM1,IKIM2,IKJP1,IKJP2,IKJM1,IKJM2,IKMAX

REAL,    DIMENSION(:,:),   ALLOCATABLE :: ZN4HGTI,ZN4HGTJ,ZMDHGTI,ZMDHGTJ
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZZMASS
REAL SCAL

!JUAN
INTEGER, DIMENSION(:,:),   ALLOCATABLE :: IKIP1_HALO2,IKIP2_HALO2,IKIM1_HALO2,IKIM2_HALO2
INTEGER, DIMENSION(:,:),   ALLOCATABLE :: IKJP1_HALO2,IKJP2_HALO2,IKJM1_HALO2,IKJM2_HALO2
INTEGER, DIMENSION(:,:),   ALLOCATABLE :: IKMAX_HALO2

REAL,    DIMENSION(:,:),   ALLOCATABLE :: ZN4HGTI_HALO2,ZN4HGTJ_HALO2,ZMDHGTI_HALO2,ZMDHGTJ_HALO2

TYPE(LIST_ll), POINTER  :: TZHGTMASS_ll
TYPE(HALO2LIST_ll), POINTER  :: TZHGTHALO2_ll

INTEGER :: KZDLB_ll,IERR
 
!JUAN
INTEGER:: IERROR                 ! DUMMY VARIABLE FOR ERROR MESSAGES

!
!*       1.    COMPUTE THE PHYSICAL SUBDOMAIN BOUNDS
!              ---------------------------------------
!
CALL GET_INDICE_ll( IIB,IJB,IIE,IJE)
IKE = SIZE(PZZ,3) - JPVEXT
IKB = 1 + JPVEXT
! Global physical dimensions
CALL GET_GLOBALDIMS_ll ( IIMAX_ll,IJMAX_ll)
IIU = SIZE(PZZ,1)
IJU = SIZE(PZZ,2)
IKU = SIZE(PZZ,3)

!PRINT*,'Interpolation coefficients for truly horizontal diffusion are computed'
ALLOCATE (ZZMASS(IIU,IJU,IKU))

!JUAN
ALLOCATE (IKIP1_HALO2(IIB-2:IIE+2,IJB-2:IJE+2),IKIP2_HALO2(IIB-2:IIE+2,IJB-2:IJE+2),&
          IKIM1_HALO2(IIB-2:IIE+2,IJB-2:IJE+2),IKIM2_HALO2(IIB-2:IIE+2,IJB-2:IJE+2),&
          IKJP1_HALO2(IIB-2:IIE+2,IJB-2:IJE+2),IKJP2_HALO2(IIB-2:IIE+2,IJB-2:IJE+2),&
          IKJM1_HALO2(IIB-2:IIE+2,IJB-2:IJE+2),IKJM2_HALO2(IIB-2:IIE+2,IJB-2:IJE+2),&
          IKMAX_HALO2(IIB-2:IIE+2,IJB-2:IJE+2))

ALLOCATE (ZN4HGTI_HALO2(IIB-2:IIE+2,IJB-2:IJE+2),ZN4HGTJ_HALO2(IIB-2:IIE+2,IJB-2:IJE+2),&
          ZMDHGTI_HALO2(IIB-2:IIE+2,IJB-2:IJE+2),ZMDHGTJ_HALO2(IIB-2:IIE+2,IJB-2:IJE+2))
!JUAN


NULLIFY(TZHGTMASS_ll,TZHGTHALO2_ll)

! Compute height field at mass points
ZZMASS = MZF(1,IKU,1,PZZ)

CALL INIT_HALO2_ll(TZHGTHALO2_ll,1,IIU,IJU,IKU)
CALL ADD3DFIELD_ll(TZHGTMASS_ll,ZZMASS)

CALL UPDATE_HALO2_ll(TZHGTMASS_ll,TZHGTHALO2_ll,IERROR)
!JUAN
PZDIFFU_HALO2%XZZ(IIB-2,IJB:IJE,:)   = TZHGTHALO2_ll%HALO2%WEST(IJB:IJE,:)
PZDIFFU_HALO2%XZZ(IIE+2,IJB:IJE,:)   = TZHGTHALO2_ll%HALO2%EAST(IJB:IJE,:)
PZDIFFU_HALO2%XZZ(IIB:IIE,IJB-2,:)   = TZHGTHALO2_ll%HALO2%SOUTH(IIB:IIE,:)
PZDIFFU_HALO2%XZZ(IIB:IIE,IJE+2,:)   = TZHGTHALO2_ll%HALO2%NORTH(IIB:IIE,:)
PZDIFFU_HALO2%XZZ(1:IIU,1:IJU,1:IKU) = ZZMASS
!print *,"zdiffu :: IIB-2",IIB-2
!print *,"zdiffu :: IIE+2",IIE+2
!print *,"zdiffu :: IJB-2",IJB-2
!print *,"zdiffu :: IJE+2",IJE+2

!DO JI=IIB-2,IIE+2
!print *,"zdiffu :: PZDIFFU_HALO2%X=JI=",JI,(PZDIFFU_HALO2%XZZ(JI,JJ,2),JJ=IJB-2,IJE+2)
!ENDDO



! Compute the vertical index of the remote grid points having the same height as the local
! grid point. This index is a real number composed of the index of the lower neighbouring model
! level and the weighting coefficient needed for the linear vertical interpolation between the
! model levels.

!JUAN
!JUAN
!JUAN

IKIP1_HALO2 = IKB
IKIP2_HALO2 = IKB
IKIM1_HALO2 = IKB
IKIM2_HALO2 = IKB
IKJP1_HALO2 = IKB
IKJP2_HALO2 = IKB
IKJM1_HALO2 = IKB
IKJM2_HALO2 = IKB

PZDIFFU_HALO2%XRKIP1 = IKB
PZDIFFU_HALO2%XRKIP2 = IKB
PZDIFFU_HALO2%XRKIM1 = IKB
PZDIFFU_HALO2%XRKIM2 = IKB
PZDIFFU_HALO2%XRKJP1 = IKB
PZDIFFU_HALO2%XRKJP2 = IKB
PZDIFFU_HALO2%XRKJM1 = IKB
PZDIFFU_HALO2%XRKJM2 = IKB

PZDIFFU_HALO2%XREDFACI = 1.0
PZDIFFU_HALO2%XREDFACJ = 1.0

CALL INDINT_HALO2( 1, 0,PZDIFFU_HALO2%XZZ,PZDIFFU_HALO2%XRKIP1,IKIP1_HALO2,IIB,IJB)
CALL INDINT_HALO2( 2, 0,PZDIFFU_HALO2%XZZ,PZDIFFU_HALO2%XRKIP2,IKIP2_HALO2,IIB,IJB)
CALL INDINT_HALO2(-1, 0,PZDIFFU_HALO2%XZZ,PZDIFFU_HALO2%XRKIM1,IKIM1_HALO2,IIB,IJB)
CALL INDINT_HALO2(-2, 0,PZDIFFU_HALO2%XZZ,PZDIFFU_HALO2%XRKIM2,IKIM2_HALO2,IIB,IJB)
CALL INDINT_HALO2( 0, 1,PZDIFFU_HALO2%XZZ,PZDIFFU_HALO2%XRKJP1,IKJP1_HALO2,IIB,IJB)
CALL INDINT_HALO2( 0, 2,PZDIFFU_HALO2%XZZ,PZDIFFU_HALO2%XRKJP2,IKJP2_HALO2,IIB,IJB)
CALL INDINT_HALO2( 0,-1,PZDIFFU_HALO2%XZZ,PZDIFFU_HALO2%XRKJM1,IKJM1_HALO2,IIB,IJB)
CALL INDINT_HALO2( 0,-2,PZDIFFU_HALO2%XZZ,PZDIFFU_HALO2%XRKJM2,IKJM2_HALO2,IIB,IJB)

PZDIFFU_HALO2%NZDI = MAX(IKIP1_HALO2,IKIP2_HALO2,IKIM1_HALO2,IKIM2_HALO2)
PZDIFFU_HALO2%NZDJ = MAX(IKJP1_HALO2,IKJP2_HALO2,IKJM1_HALO2,IKJM2_HALO2)

IKMAX_HALO2 = MAX(PZDIFFU_HALO2%NZDI,PZDIFFU_HALO2%NZDJ)

PZDIFFU_HALO2%NZDLB = MAXVAL(IKMAX_HALO2) ! Model level, above which a truly horizontal computation of diffusion
                              ! is possible at all grid points
!JUAN
CALL MPI_ALLREDUCE(PZDIFFU_HALO2%NZDLB ,KZDLB_ll, 1, MPI_INTEGER, MPI_MAX, NMNH_COMM_WORLD, IERR)

!print*,"zdiffusetup:: PZDIFFU_HALO2%NZDLB=",PZDIFFU_HALO2%NZDLB,KZDLB_ll
PZDIFFU_HALO2%NZDLB = KZDLB_ll
!JUAN

! Compute reduction factors for diffusion coefficient in I and J direction
! Their purpose is to decrease the diffusion coefficient when there is a large height difference
! among the grid points entering into the computation of diffusion

SCAL =  1.E9  ! scaling factor for coefficient reduction

DO JI = IIB,IIE
 DO JJ = IJB-1,IJE+1
    ZN4HGTI_HALO2(JI,JJ) =    PZDIFFU_HALO2%XZZ(JI+2,JJ,IKB)+PZDIFFU_HALO2%XZZ(JI-2,JJ,IKB)-  &
                     4*(PZDIFFU_HALO2%XZZ(JI+1,JJ,IKB)+PZDIFFU_HALO2%XZZ(JI-1,JJ,IKB))+ &
                     6* PZDIFFU_HALO2%XZZ(JI,JJ,IKB)
    ZMDHGTI_HALO2(JI,JJ) =   (PZDIFFU_HALO2%XZZ(JI+2,JJ,IKB)+PZDIFFU_HALO2%XZZ(JI-2,JJ,IKB)+  &
                        PZDIFFU_HALO2%XZZ(JI+1,JJ,IKB)+PZDIFFU_HALO2%XZZ(JI-1,JJ,IKB)-  &
                     4* PZDIFFU_HALO2%XZZ(JI,JJ,IKB))/4.
    PZDIFFU_HALO2%XREDFACI(JI,JJ) =  SCAL/(SCAL+ZN4HGTI_HALO2(JI,JJ)**4+ZMDHGTI_HALO2(JI,JJ)**4)

    IF (((ZN4HGTI_HALO2(JI,JJ).GT.100).AND.(ZMDHGTI_HALO2(JI,JJ).GT.10)).OR. &
        ((ZN4HGTI_HALO2(JI,JJ).GT.10).AND.(ZMDHGTI_HALO2(JI,JJ).GT.100))) THEN
      PZDIFFU_HALO2%XREDFACI(JI,JJ) = MIN(PZDIFFU_HALO2%XREDFACI(JI,JJ),0.1)
    ENDIF

  ENDDO
ENDDO
 DO JI = IIB-1,IIE+1
  DO JJ = IJB,IJE
    ZN4HGTJ_HALO2(JI,JJ) =    PZDIFFU_HALO2%XZZ(JI,JJ+2,IKB)+PZDIFFU_HALO2%XZZ(JI,JJ-2,IKB)-  &
                     4*(PZDIFFU_HALO2%XZZ(JI,JJ+1,IKB)+PZDIFFU_HALO2%XZZ(JI,JJ-1,IKB))+ &
                     6* PZDIFFU_HALO2%XZZ(JI,JJ,IKB)
    ZMDHGTJ_HALO2(JI,JJ) =   (PZDIFFU_HALO2%XZZ(JI,JJ+2,IKB)+PZDIFFU_HALO2%XZZ(JI,JJ-2,IKB)+  &
                        PZDIFFU_HALO2%XZZ(JI,JJ+1,IKB)+PZDIFFU_HALO2%XZZ(JI,JJ-1,IKB)-  &
                     4* PZDIFFU_HALO2%XZZ(JI,JJ,IKB))/4.
    PZDIFFU_HALO2%XREDFACJ(JI,JJ) =  SCAL/(SCAL+ZN4HGTJ_HALO2(JI,JJ)**4+ZMDHGTJ_HALO2(JI,JJ)**4)

    IF (((ZN4HGTJ_HALO2(JI,JJ).GT.100).AND.(ZMDHGTJ_HALO2(JI,JJ).GT.10)).OR. &
        ((ZN4HGTJ_HALO2(JI,JJ).GT.10).AND.(ZMDHGTJ_HALO2(JI,JJ).GT.100))) THEN
      PZDIFFU_HALO2%XREDFACJ(JI,JJ) = MIN(PZDIFFU_HALO2%XREDFACJ(JI,JJ),0.1)
    ENDIF

  ENDDO
ENDDO





CONTAINS


SUBROUTINE INDINT_HALO2(KII,KIJ,PZMASS,PKIND,KKMIN,KIB,KJB)

IMPLICIT NONE

INTEGER, INTENT(IN) :: KII,KIJ    ! Relative position of remote points
INTEGER, INTENT(IN) :: KIB,KJB    ! definition of domain begin
REAL, DIMENSION(KIB-2:,KJB-2:,:), INTENT(IN)       :: PZMASS  ! Height of mass points
REAL, DIMENSION(KIB-2:,KJB-2:,:), INTENT(INOUT)    :: PKIND  ! Real k index for vertical interpolation
INTEGER, DIMENSION(KIB-2:,KJB-2:),INTENT(INOUT)    :: KKMIN  ! Lowest model level for which truly horizontal computation of
                                                 ! the diffusion is possible without intersecting the ground

! Local variables
! Domain sizes
INTEGER IIB,IJB,IIE,IJE,IKB,IKE,II1,II2,IJ1,IJ2
! Loop indices
INTEGER JI,JJ,JK,JIR,JJR,JK2
REAL ZHGT   ! Height of the local grid point



CALL GET_INDICE_ll( IIB,IJB,IIE,IJE)
IKE = SIZE(PZMASS,3) - JPVEXT
IKB = 1 + JPVEXT  ! (is this correct?)

IF ((KII.EQ.0).AND.(KIJ.NE.0)) THEN


  II1 = IIB-1
  II2 = IIE+1
!JUAN   II1 = IIB-2
!JUAN   II2 = IIE+2
  IJ1 = IJB
  IJ2 = IJE


ELSE IF ((KIJ.EQ.0).AND.(KII.NE.0)) THEN
  II1 = IIB
  II2 = IIE
 IJ1 = IJB-1
 IJ2 = IJE+1
!JUAN    IJ1 = IJB-2
!JUAN    IJ2 = IJE+2



ELSE
 !callabortstop
CALL ABORT
  STOP 'Error in zdiffusetup'
ENDIF

DO JI = II1,II2
  DO JJ = IJ1,IJ2

    KKMIN(JI,JJ) = IKB

    DO JK = IKB,IKE
      JIR = JI + KII  ! location of the remote points for which the real index is to be computed
      JJR = JJ + KIJ
      ZHGT = PZMASS(JI,JJ,JK) ! height for which the index is needed

      IF (ZHGT.GT.PZMASS(JIR,JJR,JK)) THEN  ! local point higher than remote point; search upward
        IF (ZHGT.GT.PZMASS(JIR,JJR,IKE)) THEN
          PKIND(JI,JJ,JK) = IKE
          GOTO 20
        ENDIF
        DO JK2 = JK,IKE-1
          IF((ZHGT.GT.PZMASS(JIR,JJR,JK2)).AND.(ZHGT.LE.PZMASS(JIR,JJR,JK2+1))) THEN
            PKIND(JI,JJ,JK) = JK2+(ZHGT-PZMASS(JIR,JJR,JK2))/&
                              (PZMASS(JIR,JJR,JK2+1)-PZMASS(JIR,JJR,JK2))
          GOTO 20
          ENDIF
        ENDDO
  20	CONTINUE
      ELSE   ! local point lower than remote point; search downward
        IF (ZHGT.LT.PZMASS(JIR,JJR,IKB)) THEN
          PKIND(JI,JJ,JK) = IKB
          KKMIN(JI,JJ) = MAX(KKMIN(JI,JJ),JK+1)
          GOTO 25
        ELSE IF (ZHGT.EQ.PZMASS(JIR,JJR,IKB)) THEN
          PKIND(JI,JJ,JK) = IKB
          KKMIN(JI,JJ) = MAX(KKMIN(JI,JJ),JK)
          GOTO 25
        ENDIF
        DO JK2 = JK,IKB+1,-1
          IF((ZHGT.GT.PZMASS(JIR,JJR,JK2-1)).AND.(ZHGT.LE.PZMASS(JIR,JJR,JK2))) THEN
            PKIND(JI,JJ,JK) = JK2-1+(ZHGT-PZMASS(JIR,JJR,JK2-1))/&
                              (PZMASS(JIR,JJR,JK2)-PZMASS(JIR,JJR,JK2-1))
          GOTO 25
          ENDIF
        ENDDO
  25    CONTINUE
      ENDIF
    ENDDO
  ENDDO
ENDDO
IF (MINVAL(KKMIN) .EQ. 0 ) THEN
print *," zdiffusetup::PROBLEME MINVAL(KKMIN) .EQ. 0 "
call abort()
STOP
ELSE
!print *," zdiffusetup:: OK "
ENDIF
IF (MINVAL(INT(PKIND)) .EQ. 0 ) THEN
print *," zdiffusetup::PROBLEME MINVAL(INT(PKIND)) .EQ. 0 "
!PKIND = MAX (1.00001,PKIND)
call abort()
STOP
ELSE
!print *," zdiffusetup:: OK "
ENDIF
END SUBROUTINE INDINT_HALO2

END SUBROUTINE ZDIFFUSETUP

