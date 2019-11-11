!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!     ####################
      MODULE MODI_COMPUTE_SPECTRE
!     ####################
!
INTERFACE
!
      SUBROUTINE COMPUTE_SPECTRE(PDXHATM,PDYHATM, &
      PY,HFIELDSP,HOUTFILE)
!  
USE MODD_TYPE_DATE
!
!
REAL, INTENT(IN)                            :: PDXHATM  ! mean grid increment in the x direction
REAL, INTENT(IN)                            :: PDYHATM  ! mean grid increment in the y direction
!
REAL, DIMENSION(:,:,:), INTENT(IN)          :: PY       ! RHS of the equation
!                                                        
CHARACTER(LEN=*) , INTENT(IN)               :: HFIELDSP ! name of field to analyse
CHARACTER(LEN=*) , INTENT(IN)               :: HOUTFILE ! name of original file
!
END SUBROUTINE COMPUTE_SPECTRE
!
END INTERFACE
!
END MODULE MODI_COMPUTE_SPECTRE
!
!
!
!     ######spl
      SUBROUTINE COMPUTE_SPECTRE(PDXHATM,PDYHATM,   &
      PY,HFIELDSP,HOUTFILE)
!     ######################################################################
!
!****  *SPECTRE * - Compute variance spectra from grid field
!
!    PURPOSE
!    -------
!       This routine transpose a given grid field into spectral field, then
!     computes its variance spectra.
!
!   
!
!    AUTHOR
!    ------
!	A. Mary, R. Legrand          **ENM**
!       D. Ricard    **CNRM**            
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!  Philippe Wautelet: 10/01/2019: use NEWUNIT argument of OPEN
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SPECTRE, ONLY : CTYPEFILE,LSMOOTH,LSTAT
USE MODD_PARAMETERS
USE MODD_CONF
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_LUNIT_n,  ONLY: TLUOUT
USE MODE_FM
USE MODE_IO_ll
USE MODE_MSG
USE MODE_SPLITTINGZ_ll
!
USE MODI_FFT55
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
REAL, INTENT(IN)                            :: PDXHATM  ! mean grid increment in the x direction
REAL, INTENT(IN)                            :: PDYHATM  ! mean grid increment in the y direction
!
REAL, DIMENSION(:,:,:), INTENT(IN)          :: PY       ! RHS of the equation
!                                                        
CHARACTER(LEN=*) , INTENT(IN)               :: HFIELDSP ! name of field to analyse
CHARACTER(LEN=*) , INTENT(IN)               :: HOUTFILE ! name of original file
!
!*       0.2   declaration of local variables
!
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: ZY, ZY2 ! work array to store 
                                                    ! the RHS of the equation
!
INTEGER :: IIMAX        ! number of inner mass points along the x direction
INTEGER :: IJMAX        ! number of inner mass points along the y direction
INTEGER :: IKMAX        ! number of inner mass points along the k direction
INTEGER :: IKB          ! indice K for the first inner mass point along z
INTEGER :: IKU          ! size of the arrays along z
!
INTEGER :: JI,JJ,JK     ! loop indexes along x, y, z respectively
!
INTEGER :: ILOTX,ILOTY  ! number of data vectors along x, y resp. computed 
                        ! in parallel during the FFT process 
INTEGER :: INC1X,INC1Y  ! increment within each data vector for the FFT along
                        ! x, y resp.
INTEGER :: INC2X,INC2Y  ! increment between the start of one data vector and 
                        ! the next for the FFT along x,y resp.
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWORKX ! work array used by
! the FFT. It has been enlarged in order to be sufficient for 2D and 3D cases
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWORKY ! work array used by
! the FFT. It has been enlarged in order to be sufficient for 2D and 3D cases
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZBAND_X  ! intermediate array in X slices distribution
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZBAND_Y  ! intermediate array in Y slices distribution
!
INTEGER :: IINFO_ll           ! return code of parallel routine
!
INTEGER :: IIX,IJX,IIY,IJY ! dimensions of the extended x or y slices subdomain
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZBAND_YT  ! array in Y slices distribution transpose
!
! local variables for section 5
CHARACTER(len=9)                    :: YJNB       ! for error message
INTEGER                             :: JNB        ! loop index : adimensionned wavenumbers
INTEGER                             :: JERR       ! allocation or writing errors
INTEGER                             :: IMIN       ! minimum domain dimension
CHARACTER(LEN=25)                   :: YFMT       ! writing format in result file
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZVARP      ! variance parameter vector
REAL, DIMENSION(:),     ALLOCATABLE :: ZANB       ! wavenumbers vector
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZLGO       ! wavelengths and adimensionned wavenumbers vector
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZAP        ! normalized wavenumbers vector
INTEGER                             :: JLOGIMIN   ! size of wavenumber format
CHARACTER(LEN=80)                   :: YOUTFILE   ! outfile's name
REAL                                :: ZDEL       ! needed for the initialization of 
REAL                                :: ZANGLE     ! the arrays used by the FFT
REAL, DIMENSION(:), ALLOCATABLE     :: ZTRIGSX    ! arrays of sin or cos values for the FFT along x
REAL, DIMENSION(:), ALLOCATABLE     :: ZTRIGSY    ! arrays of sin or cos values for the FFT along y
INTEGER                             :: INN        ! temporary result for the computation of array TRIGS
INTEGER, DIMENSION(19)              :: IIFAXX     ! decomposition in prime numbers for the FFT along x
INTEGER, DIMENSION(19)              :: IIFAXY     ! decomposition in prime numbers for the FFT along y
INTEGER                             :: IIMAX_ll   ! Number of points of Global physical domain
INTEGER                             :: IJMAX_ll   ! Number of points of Global physical domain
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZSP        ! results table      
!              level 1 ! level 2 ! level 3  .....
!         ----!-------------------------------
!         k=1 !        !         !
!         ----!-------------------------------
!         k=2 !        !         !
!         ----!-------------------------------
!         ... !        !         !
!
INTEGER                         :: ILU
INTEGER                         :: ILUOUT    ! Logical unit number for
                                             ! output-listing   
INTEGER                         :: IRESP     ! return code in FM routines         
REAL                            :: ZMOY_C, ZMOY_S, ZVAR_C, ZVAR_S, ZVAR_S2 !computation of statistical moments
!-------------------------------------------------------------------------------
!
ILUOUT = TLUOUT%NLU
!
!*       1.    COMPUTE LOOP BOUNDS
!              -------------------
!
IIX=SIZE(PY,1)
IJX=SIZE(PY,2)
IIY=SIZE(PY,1)
IJY=SIZE(PY,2)
IIMAX = IIX-2*JPHEXT
IJMAX = IJY-2*JPHEXT
IIMAX_ll = IIX-2*JPHEXT
IJMAX_ll = IJY-2*JPHEXT
! test pour tableau 2D
IF (IIX==3 .OR. IJX==3) L2D=.TRUE.
!
IKU=SIZE(PY,3)
IKB=1+JPVEXT
IKMAX=IKU-2*JPVEXT
CALL SET_DIM_ll(IIMAX,IJMAX,IKMAX)
CALL INI_PARAZ_ll(IINFO_ll)
!
ALLOCATE(ZBAND_X(IIX,IJX,IKU))
ALLOCATE(ZBAND_Y(IIY,IJY,IKU))
ALLOCATE(ZWORKX(IIX,IJX,IKU))
ALLOCATE(ZWORKY(IIY,IJY,IKU))
IF (.NOT. L2D) THEN
  ALLOCATE(ZBAND_YT(IJY,IIY,IKU))
END IF
!
ALLOCATE(ZTRIGSX(3*IIMAX_ll))
ALLOCATE(ZTRIGSY(3*IJMAX_ll))

!-- Safety limitation
IF(PDXHATM/PDYHATM>1.2 .OR. PDYHATM/PDXHATM>1.2) THEN
  !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','COMPUTE_SPECTRE','delta X and delta Y must be close')
ENDIF
!
!-------------------------------------------------------------------------------
!*       2.    INITIALIZATION OF THE TRIGS AND IFAX ARRAYS FOR THE FFT 
!              -------------------------------------------------------
!
!        2.1 x lateral boundary conditions
!
CALL SET99(ZTRIGSX,IIFAXX,IIMAX_ll)
IF (IIFAXX(10) /=  IIMAX_ll) THEN
      WRITE(UNIT=ILUOUT,FMT="('   ERROR',/,                               &
      & '     THE FORM OF THE FFT REQUIRES'         ,/,&
      & '     THAT NIMAX MUST BE FACTORIZABLE'         ,/,&
      & '     AS A PRODUCT OF POWERS OF 2, 3 AND 5.')")
      !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','COMPUTE_SPECTRE','')
END IF 
!
!     extra trigs for shifted (co) sine transform (FFT55)
!
  INN=2*(IIMAX_ll)
  ZDEL=ASIN(1.0)/REAL(IIMAX_ll)
  DO JI=1,IIMAX_ll
    ZANGLE=REAL(JI)*ZDEL
    ZTRIGSX(INN+JI)=SIN(ZANGLE)
  END DO
!
!
!        2.2 y lateral boundary conditions
!
IF (.NOT. L2D) THEN 
  CALL SET99(ZTRIGSY,IIFAXY,IJMAX_ll)
  IF (IIFAXY(10) /=  IJMAX_ll) THEN
      WRITE(UNIT=ILUOUT,FMT="('   ERROR',/,                               &
      & '     THE FORM OF THE FFT REQUIRES'         ,/,&
      & '     THAT NJMAX MUST BE FACTORIZABLE'         ,/,&
      & '     AS A PRODUCT OF POWERS OF 2, 3 AND 5.')")
      !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','COMPUTE_SPECTRE','')
  END IF 
 !
 !     extra trigs for shifted (co) sine transform
 !
   INN=2*(IJMAX_ll)
   ZDEL=ASIN(1.0)/REAL(IJMAX_ll)
   DO JJ=1,IJMAX_ll
     ZANGLE=REAL(JJ)*ZDEL
     ZTRIGSY(INN+JJ)=SIN(ZANGLE)
   END DO
 !
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTE THE ARRAY INCREMENTS FOR THE FFT
!              ----------------------------------------
!
IF(.NOT. L2D) THEN
!
  ILOTX = IJX*IKU
  INC1X = 1
  INC2X = IIX
!
  ILOTY = IIY*IKU
  INC1Y = 1
  INC2Y = IJY
!
ELSE
!
  ILOTX = IKU
  INC1X = 1
  INC2X = IIX*IJX
ENDIF
!-------------------------------------------------------------------------------
!
!*      4.    FORM HOMOGENEOUS BOUNDARY CONDITIONS FOR A NONCYCLIC CASE
!             ---------------------------------------------------------
!
!
!*      4.1 copy the RHS in a local array REMAP functions will shift the indices for the FFT
!
ZY = PY
!
!	4.3 2way structure -> xslice structure, + data shift
!
ZBAND_X=0.
CALL REMAP_2WAY_X_ll(ZY,ZBAND_X,IINFO_ll)
!
!  
!-------------------------------------------------------------------------------
!*      5.    APPLY A REAL TO COMPLEX FFT
!             ---------------------------
!
!
CALL FFT55(ZBAND_X(1,1,IKB-1),ZWORKX,ZTRIGSX,IIFAXX,INC1X,INC2X,  &
             IIMAX,ILOTX,-1 )
!
!
ZBAND_Y=0.
CALL REMAP_X_Y_ll(ZBAND_X,ZBAND_Y,IINFO_ll)
!
IF (.NOT. L2D) THEN                                 
!
! array transposition I --> J
!
CALL FAST_TRANSPOSE(ZBAND_Y,ZBAND_YT,IIY,IJY,IKU)
!
CALL FFT55(ZBAND_YT(1,1,IKB-1),ZWORKY,ZTRIGSY,IIFAXY,INC1Y,INC2Y,  &
                  IJMAX,ILOTY,-1 )
!
END IF    
!
!-------------------------------------------------------------------------------
!
!*	6.   SPECTRA COMPUTATION
!            -------------------
!*	6.1 initialisation of values
!
!	computation of the minimum domain dimension
IF (IIX<IJY) THEN
  IMIN=IIX-2
ELSE
  IMIN=IJY-2
END IF
!
!	allocation errors
ALLOCATE(ZSP(IMIN-1,IKU))
ALLOCATE(ZANB(IMIN))
ALLOCATE(ZLGO(IMIN,2))
ALLOCATE(ZVARP(IJY,IIY,IKU))
ALLOCATE(ZAP(IJY,IIX))
! initialisation
ZSP=0.
ZANB=0.
ZLGO=0.
ZAP=0.
!
!*	6.2 algorithm
!	writing of the wavenumbers and wavelengths
DO JNB=1,IMIN
  ZANB(JNB)=real(JNB)/real(IMIN)
  ZLGO(JNB,1)=JNB
  ZLGO(JNB,2)=(PDXHATM+PDYHATM)/ZANB(JNB)
END DO
!
!
!	computation of variances, gathering of variances by wavenumbers
DO JI=1,IIX-2                                      
 DO JJ=1,IJY-2
   IF (JI/=1 .OR. JJ/=1) THEN
     ZAP(JJ,JI)=((real(JI-1)/real(IIX-2))**2+(real(JJ-1)/real(IJY-2))**2)**(1/2.)
     DO JNB=1,IMIN-1
       IF ( (JNB.EQ.1 .OR. ZANB(JNB) <= ZAP(JJ,JI)) .AND. ZAP(JJ,JI) < ZANB(JNB+1) ) THEN
         DO JK=1,IKU
            IF(JI==1 .OR. JJ==1) THEN
              ZVARP(JJ,JI,JK)=ZBAND_YT(JJ,JI,JK)**2/2.
            ELSE
              ZVARP(JJ,JI,JK)=ZBAND_YT(JJ,JI,JK)**2/4.
            ENDIF
            IF (LSMOOTH) THEN
              IF (JNB==IMIN-1) THEN 
                ZSP(JNB,JK)=ZSP(JNB,JK)+ZVARP(JJ,JI,JK)*(ZANB(JNB+1)-ZAP(JJ,JI))/(ZANB(JNB+1)-ZANB(JNB))
              ELSE ! introducing smoothing in distribution of spectral coefficients
                IF(ZANB(JNB) <= ZAP(JJ,JI)) THEN
                  !general case
                  ZSP(JNB,JK)=ZSP(JNB,JK)+ZVARP(JJ,JI,JK)*(ZANB(JNB+1)-ZAP(JJ,JI))/(ZANB(JNB+1)-ZANB(JNB))
                  ZSP(JNB+1,JK)=ZSP(JNB+1,JK)+ZVARP(JJ,JI,JK)*(ZAP(JJ,JI)-ZANB(JNB))/(ZANB(JNB+1)-ZANB(JNB))
                ELSE
                  !special case for some coefficients in the k=1 ring
                  ZSP(JNB,JK)=ZSP(JNB,JK)+ZVARP(JJ,JI,JK)
                END IF
              END IF
            ELSE
              ZSP(JNB,JK)=ZSP(JNB,JK)+ZVARP(JJ,JI,JK)
            ENDIF
          END DO
        END IF
     END DO
   END IF
 END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*	7.   STATISTICS
!            ----------
IF(LSTAT) THEN
  !*	7.1 Direct statistics
  DO JK=2,IKU-1
    ZMOY_C=0.
    ZVAR_S=0.
    ZVAR_S2=0.
    ZVAR_C=0.
    ZMOY_S=ZBAND_YT(1,1,JK)
    DO JI=2,SIZE(PY,1)-1
      DO JJ=2,SIZE(PY,2)-1
        ZMOY_C=ZMOY_C+PY(JI,JJ,JK)
        ZVAR_C=ZVAR_C+PY(JI,JJ,JK)**2
      ENDDO
    ENDDO
    DO JI=1,IIX-2                                      
      DO JJ=1,IJY-2
        IF(JI/=1 .OR. JJ/=1) THEN
          IF(JI==1 .OR. JJ==1) THEN
            ZVAR_S=ZVAR_S+ZBAND_YT(JJ,JI,JK)**2/2.
          ELSE
            ZVAR_S=ZVAR_S+ZBAND_YT(JJ,JI,JK)**2/4.
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    DO JNB=1,IMIN-1
      ZVAR_S2=ZVAR_S2+ZSP(JNB,JK)
    ENDDO
    ZMOY_C=ZMOY_C/(SIZE(PY,1)-2)/(SIZE(PY,2)-2)
    ZVAR_C=ZVAR_C/(SIZE(PY,1)-2)/(SIZE(PY,2)-2)-ZMOY_C**2
    print*,"Statistiques sur le champ, au niveau ",JK
    print*,"Moyenne du champ  (calculée à partir du champ en point de grille)=",ZMOY_C
    print*,"Moyenne du champ  (calculée à partir du champ en spectral)       =",ZMOY_S
    print*,"Variance du champ (calculée à partir du champ en point de grille)=",ZVAR_C
    print*,"Variance du champ (calculée à partir du champ en spectral)       =",ZVAR_S
    print*,"Variance du champ (calculée à partir du spectre)                 =",ZVAR_S2
    print*,"Résidus                                                          =",ZVAR_S-ZVAR_S2
  ENDDO
  !*	7.1 Second statistics
  !APPLY A COMPLEX TO REAL FFT
  IF (.NOT. L2D) THEN
    CALL FFT55( ZBAND_YT(1,1,IKB-1),ZWORKY,ZTRIGSY,IIFAXY,INC1Y,INC2Y,  &
                    IJMAX,ILOTY,+1 )
    ! array transposition J --> I
    CALL FAST_TRANSPOSE(ZBAND_YT,ZBAND_Y,IJY,IIY,IKU)
  ENDIF
  !
  ! Transposition Y-> X
  !
  ZBAND_X=0.
  CALL REMAP_Y_X_ll(ZBAND_Y,ZBAND_X,IINFO_ll)
  !
  !
  CALL FFT55( ZBAND_X(1,1,IKB-1),ZWORKX,ZTRIGSX,IIFAXX,INC1X,INC2X,  &
                  IIMAX,ILOTX,+1 )
  !
  !RETURN TO A NON HOMOGENEOUS NEUMAN CONDITION FOR NON-CYCLIC CASES
  ! Transposition + shift X -> 2way
  !
  CALL REMAP_X_2WAY_ll(ZBAND_X,ZY2,IINFO_ll)
  DO JK=2,IKU-1
    ZMOY_C=0.
    ZVAR_C=0.
    DO JI=2,SIZE(PY,1)-1
      DO JJ=2,SIZE(PY,2)-1
        ZMOY_C=ZMOY_C+PY(JI,JJ,JK)
        ZVAR_C=ZVAR_C+PY(JI,JJ,JK)**2
      ENDDO
    ENDDO
    ZMOY_C=ZMOY_C/(SIZE(PY,1)-2)/(SIZE(PY,2)-2)
    ZVAR_C=ZVAR_C/(SIZE(PY,1)-2)/(SIZE(PY,2)-2)-ZMOY_C**2
    print*,"Statistiques sur le champ, après aller-retour spectral, au niveau ",JK
    print*,"Moyenne du champ  (calculée à partir du champ en point de grille)=",ZMOY_C
    print*,"Variance du champ (calculée à partir du champ en point de grille)=",ZVAR_C
    print*,"Différence mini et maxi entre les champs (avant et après aller-retour):",&
      MINVAL(ZY(2:SIZE(PY,1)-1,2:SIZE(PY,2)-1,JK)-ZY2(2:SIZE(PY,1)-1,2:SIZE(PY,2)-1,JK)),&
      MAXVAL(ZY(2:SIZE(PY,1)-1,2:SIZE(PY,2)-1,JK)-ZY2(2:SIZE(PY,1)-1,2:SIZE(PY,2)-1,JK))
  ENDDO
ENDIF
!
!-------------------------------------------------------------------------------
!
!*	7.   WRITING IN RESULT FILE
!            -------------------
!    
!* 7.1 Name of result file
!
WRITE(YOUTFILE,FMT='(A,A,A)') TRIM(ADJUSTL(HOUTFILE)),"_",TRIM(ADJUSTL(HFIELDSP))
! 
!* 7.2 Write result in file
!
JLOGIMIN=int(log(real(IMIN-1)))
IF (IKU < 100) THEN
  WRITE(YFMT,FMT='(A,I1,A,I2,A)') "(I",JLOGIMIN,",F13.2,",IKU,"F25.16)"
ELSE
  WRITE(YFMT,FMT='(A,I1,A,I3,A)') "(I",JLOGIMIN,",F13.2,",IKU,"F25.16)"
ENDIF
YFMT=TRIM(ADJUSTL(YFMT))
!
OPEN(NEWUNIT=ILU, FILE=YOUTFILE, ACCESS='SEQUENTIAL', IOSTAT=JERR)
IF (JERR /= 0) THEN
  CALL PRINT_MSG(NVERB_FATAL,'IO','COMPUTE_SPECTRE','error when opening '//trim(YOUTFILE))
ENDIF
!
DO JNB=1,IMIN-1
  WRITE(UNIT=ILU,IOSTAT=JERR,FMT=YFMT) INT(ZLGO(JNB,1)),ZLGO(JNB,2),ZSP(JNB,:)
  IF (JERR /= 0) THEN
    WRITE(YJNB,'( I9 )') JNB
    CALL PRINT_MSG(NVERB_FATAL,'IO','COMPUTE_SPECTRE','error when writing JNB='//trim(YJNB)//' in file '//trim(YOUTFILE))
  ENDIF
END DO
!
CLOSE(UNIT=ILU, IOSTAT=JERR)
IF (JERR /= 0) THEN 
  CALL PRINT_MSG(NVERB_ERROR,'IO','COMPUTE_SPECTRE','error when closing '//trim(YOUTFILE))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*	7.   DEALLOCATION OF ARRAYS
!            -------------------
!
DEALLOCATE(ZBAND_X)
DEALLOCATE(ZBAND_Y)
IF (.NOT. L2D) THEN
  DEALLOCATE(ZBAND_YT)
END IF
DEALLOCATE(ZWORKX)
DEALLOCATE(ZWORKY)
DEALLOCATE(ZANB)
DEALLOCATE(ZLGO)
DEALLOCATE(ZSP)
DEALLOCATE(ZVARP)
DEALLOCATE(ZAP)
DEALLOCATE(ZTRIGSX)
DEALLOCATE(ZTRIGSY)

!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
CONTAINS
  SUBROUTINE FAST_TRANSPOSE(PX,PXT,KNI,KNJ,KNK)
    INTEGER                      :: KNI,KNJ,KNK ! 3D dimension of X and XT
    REAL, DIMENSION(KNI*KNJ,KNK) :: PX
    REAL, DIMENSION(KNJ*KNI,KNK) :: PXT
    !
    INTEGER                      :: IJI,II,IJ,IIJ ! index in array X and XT
    INTEGER                      :: JK
!
    DO JK=1,KNK
       ! PERMUTATION(PX,PXT)
       !CDIR NODEP
       !OCL NOVREC
       DO IJI = 1, KNJ*KNI
          ! I,J Indice in XT array from linearised index IJI
          II   = 1 +    (IJI-1)/KNJ
          IJ   = IJI - (II-1)*KNJ 
          ! linearised index in X
          IIJ = II + (IJ-1)*KNI 
          ! transposition
          PXT(IJI,JK) = PX(IIJ,JK)
!       
       END DO
    END DO
!    
END SUBROUTINE FAST_TRANSPOSE
!
!------------------------------------------------------------------------------  
END SUBROUTINE COMPUTE_SPECTRE
