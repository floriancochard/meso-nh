!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 solver 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ######spl
      SUBROUTINE FFT55(PA,PWORK,PTRIGS,KIFAX,KINC,KJUMP,KN,KLOT,KISIGN)
!     #################################################################
!
!!****  *FFT55 * - multiple fast real staggered (shifted) cosine transform
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to apply either a staggered cosine 
!     transform (ISIGN = -1) or a real periodic transform (ISIGN = +1).
!
!!**  METHOD
!!    ------
!!      The staggered cosine transform of length N is converted to real 
!!    periodic  transform by pre- and post- processing:
!!         ISIGN=+1: X(I)=SUM(J=0,...,N-1)(Z(J)*COS((2*I-1)*J*PI/(2*N))
!!
!!      The real periodic transform is performed by pruning redundant operations
!!    from complex transform:
!!       ISIGN=-1: Z(J)=(2/N)*SUM(I=1,...,N)(X(I)*COS((2*I-1)*J*PI/(2*N))
!!                           (scaling is 1/N for J=0)
!!      
!!      The input/output are ordered as follows: 
!!         data:            X(1) , X(2) , ... , X(N)
!!         coefficients:    Z(0) , Z(1) , ... , Z(N-1)  
!!
!!      Vectorization is achieved on CRAY by doing the transforms in parallel.
!!      
!!      N must be composed of factors 2,3 & 5 and must be even.
!!      -------------------------------------------------------
!!
!!      For this version enough space must be allowed for data vectors of 
!!    length (N+2).
!!     
!!    EXTERNAL
!!    --------
!!      Subroutine RFFTMLT : apply real-to-complex or complex-to-real Fast 
!!    Fourier Transform (FFT) on multiple input vectors.
!!      Subroutine FFT991  : equivalent to RFFTMLT
!!     
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (subroutine FFT55) 
!!      
!!    AUTHOR
!!    ------
!!	Clive Temperton       
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/02/87 
!!      Revision C. Temperton (december 91)
!!      Revision P. Jabouille (juin 96) replace the CRAY intrinsic function 
!!                   RFFTMLT by the arpege routine FFT991
!!      Revision J. Stein and P. Jabouille (juillet 96) extend the pre- 
!!                   and post-processing to the odd number
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
                                                 
REAL, DIMENSION(*), INTENT(IN) :: PTRIGS         ! previously prepared list of
                                                 ! trig function values
!
INTEGER, DIMENSION(19), INTENT(IN) :: KIFAX       ! previously prepared list of
                                                 ! factors of N 
!
INTEGER, INTENT(IN) :: KINC              ! increment within each data 'vector' 
                                         ! (e.g. INC = 1 for consectuvely stored
                                         ! vectors
!
INTEGER, INTENT(IN) :: KJUMP             ! increment between the start of each 
                                         ! data vector
!
INTEGER, INTENT(IN) :: KN                ! length of each data vector
!
INTEGER, INTENT(IN) :: KLOT              ! number of data vectors
!
INTEGER, INTENT(IN) :: KISIGN  ! +1 for transform from spectral to gridpoint
                               ! -1 for transform from gridpoint to spectral
!
REAL, DIMENSION(*), INTENT(OUT) :: PWORK ! area of size (2*N)*MIN(LOT,1020)
!
REAL, DIMENSION(*), INTENT(INOUT) :: PA ! input and output data
!
!*       0.2   declarations of local variables 
!
! 
INTEGER :: INH, INN, INBLOX, INVEX, ISTART, IJA, IJB, IIA,IIB, IIS
INTEGER :: IJA0, IJB0, IIA0, IIB0    ! cvj to promote vectorization
INTEGER :: IEVEN ! evenness flag
!
INTEGER :: JNB, JJ, JK ! loop control
!
REAL :: ZSCALE, ZSCALE1, ZSCALE2
!
REAL :: ZT1, ZT2, ZCO, ZRI
!  
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE LOOP BOUNDS
!              -------------------
IEVEN=MOD(KN,2)
INH=KN/2
INN=2*KN
INBLOX=1+(KLOT-1)/1020
INVEX=KLOT-(INBLOX-1)*1020
IF (KISIGN.EQ.1) THEN
!
!-------------------------------------------------------------------------------
!*       2.    KISIGN=+1, SPECTRAL TO GRIDPOINT TRANSFORM
!              ------------------------------------------
!
  ISTART=1
  DO JNB=1,INBLOX
!
!      2.1     preprocessing
!              -------------
!
    ZSCALE=0.5*FLOAT(KN)
! this loop works for odd and even case
    DO JK=1,(KN-1)/2
      IJA=JK+1
      IJB=KN+1-JK
!cvj  IIA=ISTART+(IJA-1)*KINC
!cvj  IIB=ISTART+(IJB-1)*KINC
      IIA0=ISTART+(IJA-1)*KINC
      IIB0=ISTART+(IJB-1)*KINC
      ZRI=ZSCALE*PTRIGS(INN+JK)
      ZCO=ZSCALE*PTRIGS(INN+KN-JK)
!
!ocl novrec
!CDIR NODEP
      DO JJ=1,INVEX
        IIA=IIA0+(JJ-1)*KJUMP
        IIB=IIB0+(JJ-1)*KJUMP
        ZT1 = ZRI*(PA(IIA)+PA(IIB)) - ZCO*(PA(IIB)-PA(IIA))
        ZT2 = ZRI*(PA(IIB)-PA(IIA)) + ZCO*(PA(IIA)+PA(IIB))
        PA(IIA) = ZT1
        PA(IIB) = ZT2
!cvj    IIA=IIA+KJUMP
!cvj    IIB=IIB+KJUMP
      ENDDO
    ENDDO
!
!cvjIIA=ISTART+INH*KINC
!cvjIIB=ISTART
    IIA0=ISTART+INH*KINC
    IIB0=ISTART
    ZSCALE1=2.0*ZSCALE
    ZSCALE2=2.0*ZSCALE*PTRIGS(INN+INH)
!
!ocl novrec
!CDIR NODEP
    DO JJ=1,INVEX
      IIA=IIA0+(JJ-1)*KJUMP
      IIB=IIB0+(JJ-1)*KJUMP
      IF(IEVEN == 0) PA(IIA)=ZSCALE2*PA(IIA)
      PA(IIB)=ZSCALE1*PA(IIB)
!cvj  IIA=IIA+KJUMP
!cvj  IIB=IIB+KJUMP
    ENDDO
!
!      2.2     periodic Fourier analysis
!              -------------------------
!
    IIS=-1
    IIA=ISTART
    CALL FFT991(PA(IIA),PWORK,PTRIGS,KIFAX,KINC,KJUMP,KN,INVEX,IIS)
!
!      2.3     postprocessing
!              --------------
    DO JK=2,KN,2
!cvj  IJA=ISTART+(JK-1)*KINC
      IJA0=ISTART+(JK-1)*KINC
!
!ocl novrec
!CDIR NODEP
      DO JJ=1,INVEX
        IJA=IJA0+(JJ-1)*KJUMP
        PA(IJA)=PA(IJA+KINC)-PA(IJA+2*KINC)
        PA(IJA+KINC)=PA(IJA+KINC)+PA(IJA+2*KINC)
!cvj    IJA=IJA+KJUMP
      ENDDO
    ENDDO
!
      ISTART=ISTART+INVEX*KJUMP
      INVEX=1020
  ENDDO
!
ELSE
!
!-------------------------------------------------------------------------------
!*      3.     KISIGN=-1, gridpoint to spectral transform
!              ------------------------------------------
!
  ISTART=1
  DO JNB=1,INBLOX
!
!      3.1     preprocessing
!              -------------
!
!cvjIIA=ISTART
!cvjIIB=IIA+(KN-1)*KINC
    IIA0=ISTART
    IIB0=IIA0+(KN-1)*KINC
!
!ocl novrec
!CDIR NODEP
    DO JJ=1,INVEX
      IIA=IIA0+(JJ-1)*KJUMP
      IIB=IIB0+(JJ-1)*KJUMP
      PA(IIA)=2.0*PA(IIA)
      IF(IEVEN == 0) PA(IIB+KINC)=2.0*PA(IIB)
!cvj  IIA=IIA+KJUMP
!cvj  IIB=IIB+KJUMP
    ENDDO
!
! this loop works for odd and even case
    DO JK=KN-2+IEVEN,2,-2
!cvj  IIA=ISTART+(JK-1)*KINC
      IIA0=ISTART+(JK-1)*KINC
!
!ocl novrec
!CDIR NODEP
      DO JJ=1,INVEX
        IIA=IIA0+(JJ-1)*KJUMP
        PA(IIA+2*KINC)=PA(IIA+KINC)-PA(IIA)
        PA(IIA+KINC)=PA(IIA+KINC)+PA(IIA)
!cvj    IIA=IIA+KJUMP
      ENDDO
    ENDDO
!
!      3.2     periodic Fourier synthesis
!              --------------------------
!
    IIA=ISTART
    IIS=1
    CALL FFT991(PA(IIA),PWORK,PTRIGS,KIFAX,KINC,KJUMP,KN,INVEX,IIS)
!
!      3.3     postprocessing
!              --------------
!
    ZSCALE=0.5/FLOAT(KN)
! this loop works for odd and even case
    DO JK=1,(KN-1)/2
      IIA=JK+1
      IIB=KN-JK+1
!cvj  IJA=ISTART+(IIA-1)*KINC
!cvj  IJB=ISTART+(IIB-1)*KINC
      IJA0=ISTART+(IIA-1)*KINC
      IJB0=ISTART+(IIB-1)*KINC
      ZRI=ZSCALE*PTRIGS(INN+JK)
      ZCO=ZSCALE*PTRIGS(INN+KN-JK)
!
!ocl novrec
!CDIR NODEP
      DO JJ=1,INVEX
        IJA=IJA0+(JJ-1)*KJUMP
        IJB=IJB0+(JJ-1)*KJUMP
        ZT1 = ZRI*(PA(IJA)-PA(IJB)) + ZCO*(PA(IJA)+PA(IJB))
        ZT2 = ZRI*(PA(IJA)+PA(IJB)) - ZCO*(PA(IJA)-PA(IJB))
        PA(IJA) = ZT1
        PA(IJB) = ZT2
!cvj    IJA=IJA+KJUMP
!cvj    IJB=IJB+KJUMP
      ENDDO
    ENDDO
!
!cvjIJA=ISTART+INH*KINC
!cvjIJB=ISTART
    IJA0=ISTART+INH*KINC
    IJB0=ISTART
    ZSCALE1=2.0*ZSCALE*PTRIGS(INN+INH)
!
!ocl novrec
!CDIR NODEP
    DO JJ=1,INVEX
      IJA=IJA0+(JJ-1)*KJUMP
      IJB=IJB0+(JJ-1)*KJUMP
      IF(IEVEN == 0) PA(IJA)=ZSCALE1*PA(IJA)
      PA(IJB)=ZSCALE*PA(IJB)
!cvj  IJA=IJA+KJUMP
!cvj  IJB=IJB+KJUMP
    ENDDO
!
    ISTART=ISTART+INVEX*KJUMP
    INVEX=1020
  ENDDO
!
ENDIF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE FFT55
