!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ########################
      MODULE MODE_TOOLSZ_ll
!     ########################
!!    Author
!!    ------
!     J. Escobar    * LA *
!
!!    Modifications
!!    -------------
!     Original 3/02/2009


  CONTAINS

!!$  SUBROUTINE SLIDE_COORD(KDIM_DATA,KDIM_PROC,THIS_PROC,KOR,KEND)
!!$
!!$    !!    Purpose
!!$    !
!!$    !   Compute for the processor=THIS_PROC the origine/end of slide in decomposing
!!$    !   an array of data of dimension=KDIM_DATA on KDIM_PROC
!!$    !
!!$    !!    Author
!!$    !!    ------
!!$    !     J. ESCOBAR    * LA *
!!$
!!$    IMPLICIT NONE
!!$    !
!!$    !*       0.1   declarations of arguments
!!$    !
!!$    INTEGER, INTENT(IN)  :: KDIM_DATA ! dimension of data to split
!!$    INTEGER, INTENT(IN)  :: KDIM_PROC ! numbers of processor to use in splitting
!!$    INTEGER, INTENT(IN)  :: THIS_PROC ! processor id from 1..NB_PROC
!!$    INTEGER, INTENT(OUT) :: KOR,KEND  ! Origine/End coordonate
!!$    !
!!$    !*       0.2   declarations of local variables
!!$    !
!!$    INTEGER             :: IDIM_SLIDE ! slide dimension ( without rest/delta )
!!$    INTEGER             :: IREST      ! number of point in surabondance to distribut 
!!$    INTEGER             :: IDELTAOR,IDELTAEND     ! offset in origine to apply 
!!$
!!$    IDIM_SLIDE   = KDIM_DATA/KDIM_PROC
!!$    IREST        = MOD(KDIM_DATA,KDIM_PROC)
!!$    IDELTAOR     = MIN(IREST,THIS_PROC-1)
!!$    IDELTAEND    = MIN(IREST,THIS_PROC)
!!$
!!$    KOR   = ( THIS_PROC - 1 ) * IDIM_SLIDE + 1 + IDELTAOR
!!$    KEND  =   THIS_PROC       * IDIM_SLIDE     + IDELTAEND
!!$
!!$  END SUBROUTINE SLIDE_COORD

  !###################################################################
  SUBROUTINE CARTESIANZ(TPROC,NB_PROC,X_DIM,Y_DIM,Z_DIM,X_DOMAINS,Y_DOMAINS,Z_DOMAINS,KORDER)
    !###################################################################
    !
    !!****  *CARTESIAN* - routine which splits a domain if NB_PROC 
    !  		      is not a prime number
    !
    !!    Purpose
    !!    -------
    !     this routine fills the elements of TPROC.
    !
    !!**  Method
    !!    ------
    !
    !!    External
    !!    --------
    !
    !!    Implicit Arguments
    !!    ------------------
    !
    !     Module MODD_STRUCTURE_ll
    !       type ZONE_LL
    !
    !!    Reference
    !!    ---------
    !
    !!    Author
    !!    ------
    !     J. ESCOBAR  * LA *
    !
    !!    Modifications
    !!    -------------
    !     Original 19/11/2007
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : ZONE_LL
    USE MODE_TOOLS_ll,    ONLY  : SLIDE_COORD
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    TYPE(ZONE_LL), INTENT(INOUT), DIMENSION(:),TARGET  :: TPROC                         ! split structur to file
    INTEGER      , INTENT(IN)                   :: NB_PROC                       ! number of processors
    INTEGER      , INTENT(IN)                   :: X_DIM,Y_DIM,Z_DIM             ! global data_grid dimension
    INTEGER      , INTENT(IN)                   :: X_DOMAINS,Y_DOMAINS,Z_DOMAINS ! processors_grid dimension
    INTEGER      , INTENT(IN)                   :: KORDER                        ! order of mapping of processor
    !
    !*       0.2   declarations of local variables
    !
    INTEGER :: IDOM,II,IJ,IK
    INTEGER :: XORDER,YORDER,ZORDER, IP_NUMBER
    !JUAN
    TYPE(ZONE_LL), POINTER  :: TP
    !JUAN
    !
    !-------------------------------------------------------------------------------
    !
    !*	 1.    COMPUTE THE AVERAGE DIMENSION
    !
    !
    !-------------------------------------------------------------------------------
    !
    !*	 2.    FILL THE FIELDS OF TPROC
    ! 
    IF      ( KORDER .EQ. 321 ) THEN ! ZYX MAPPING
       ZORDER = 1 ; YORDER = Z_DOMAINS ; XORDER = Z_DOMAINS * Y_DOMAINS
    ELSE IF ( KORDER .EQ. 312 ) THEN ! ZXY MAPPING
       ZORDER = 1 ; XORDER = Z_DOMAINS ; YORDER = Z_DOMAINS * X_DOMAINS
    ELSE IF ( KORDER .EQ. 213 ) THEN ! YXZ MAPPING / BSPLITTING
       YORDER = 1 ; XORDER = Y_DOMAINS ; ZORDER = Y_DOMAINS * X_DOMAINS
    ELSE IF ( KORDER .EQ. 123 ) THEN ! XYZ MAPPING / BSPLITTING
       XORDER = 1 ; YORDER = X_DOMAINS ; ZORDER = X_DOMAINS * Y_DOMAINS
    END IF
    IDOM = 0
    DO IK=1,Z_DOMAINS
       DO IJ=1,Y_DOMAINS
          DO II=1,X_DOMAINS
             !
             !       file processor number
             !
             IDOM = IDOM + 1
             !JUANZ TP => TPROC(IDOM)
             IP_NUMBER = 1 + ( II - 1 ) * XORDER + ( IJ - 1 ) * YORDER + ( IK - 1 ) * ZORDER
             TP => TPROC(IP_NUMBER)
             !JUANZ TP%NUMBER = IDOM
             TP%NUMBER = IP_NUMBER
             !
             !       compute/file x coordonate
             !
             CALL SLIDE_COORD(X_DIM,X_DOMAINS,II,TP%NXOR,TP%NXEND)
             !
             !       compute/file y coordonate
             !
             CALL SLIDE_COORD(Y_DIM,Y_DOMAINS,IJ,TP%NYOR,TP%NYEND)
             !
             !       compute/file Z coordonate
             !
             CALL SLIDE_COORD(Z_DIM,Z_DOMAINS,IK,TP%NZOR,TP%NZEND)
             !
          END DO
          !
       END DO

    END DO
    ! 
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE CARTESIANZ

  !     #################################################################
  SUBROUTINE SPLITZ(X_DIM,Y_DIM,Z_DIM,NB_PROC,TPROC,HSPLITTING,KZ_PROC, KX_DOMAINS,KY_DOMAINS)
    !     #################################################################
    !
    !!****  *SPLITZ* - routine which splits a domain in NB_PROC sub-domains
    !                     by using DEFINE_SPLITTING2.
    !
    !!    Purpose
    !!    -------
    !     this routine fills the fields of TPROC.
    !
    !!**  Method
    !!    ------
    !
    !!    External
    !!    --------
    !
    !!    Implicit Arguments
    !!    ------------------
    !
    !     Module MODD_STRUCTURE_ll
    !       type ZONE_LL
    !
    !     Module MODD_VAR_ll
    !       JPHEXT
    !
    !!    Reference
    !!    ---------
    !
    !!    Author
    !!    ------
    !     D. Lugato    * CERFACS *
    !
    !!    Modifications
    !!    -------------
    !     Original 01/05/98
    !     R.Guivarch 29/11/99 : x and y splitting : HSPLITTING
    !     J.Escobar 28/03/2019: check very small domain(0 size) 
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : ZONE_LL
    USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT
    !JUAN
    USE MODD_VAR_ll, ONLY : IP
    USE MODD_CONFZ , ONLY : NZ_VERB,NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two

    USE  MODE_SPLITTING_ll , ONLY : def_splitting2
    USE MODE_TOOLS_ll   ,    ONLY : SLIDE_COORD
    !JUAN
    !
    IMPLICIT NONE
    !
    !
    !*       0.1   declarations of arguments
    !
    INTEGER, INTENT(IN) :: NB_PROC,X_DIM,Y_DIM,Z_DIM
    CHARACTER*10, INTENT(IN) :: HSPLITTING ! kind of splitting
    TYPE(ZONE_LL), INTENT(OUT), DIMENSION(NB_PROC),TARGET  :: TPROC
    !
    !JUAN
    INTEGER, INTENT(IN) :: KZ_PROC ! number of proc for Z splitting
    INTEGER, INTENT(IN), OPTIONAL :: KX_DOMAINS,KY_DOMAINS

    !JUAN
    !
    !*       0.2   declarations of local variables
    !
    INTEGER                 :: X_DOMAINS,Y_DOMAINS,Z_DOMAINS,X_DOMAINS_NEW
    LOGICAL                 :: PREM
    INTEGER                 :: IK
    INTEGER                 :: IDOM 
    TYPE(ZONE_LL), POINTER  :: TP
    INTEGER                 :: NB_PROC_XY
    !
    !        0. CHECK NB_PROC/NZ_PROC
    PREM = .FALSE.
    IF ( MOD(NB_PROC,KZ_PROC) .NE. 0 ) THEN
       PRINT*
       WRITE(*,1000) NB_PROC, KZ_PROC
       PRINT*
1000   FORMAT("MODE_SPLITTINGZ::SPLITZ --> NB_PROC=", I4 ," NOT DIVISIBLE BY KZ_PROC=", I4)
       STOP
    ENDIF
    !
    !   Splitting in Z possible so
    !
    NB_PROC_XY = NB_PROC / KZ_PROC
    Z_DOMAINS  = KZ_PROC
    !
    !-------------------------------------------------------------------------------
    !
    !*	 1.    FIND THE SPLITTING in XY & Z 
    !
    !  CALL DEF_SPLITTINGZ(X_DOMAINS,Y_DOMAINS,Z_DOMAINS,X_DIM,Y_DIM,Z_DIM,NB_PROC,KZ_PROC,PREM)
    IF (PRESENT(KX_DOMAINS) ) THEN
       X_DOMAINS = KX_DOMAINS
       Y_DOMAINS = KY_DOMAINS
    ELSE
    CALL DEF_SPLITTING2(X_DOMAINS,Y_DOMAINS,X_DIM,Y_DIM,NB_PROC_XY,PREM)
    END IF
    !
    ! transpose the data X & Y
    X_DOMAINS_NEW = Y_DOMAINS
    Y_DOMAINS     = X_DOMAINS
    X_DOMAINS     = X_DOMAINS_NEW
    !
    !-------------------------------------------------------------------------------
    !
    !*	 2.    FILL THE FIELDS OF TPROC
    !
    IF(HSPLITTING.EQ."P2P1SPLITT") THEN
       IF ((PREM).AND.(NB_PROC_XY.GT.2)) THEN
          STOP "mode_toolsz_ll.f90::SPLITZ: NPROC PREMIER NON PREVUE !!! "
          !
          !   split x direction only on NB_PROC_XY - 1 processors
          !   and on reducted x-size = X_DIM - X_DIM/NB_PROC_XY -1
          !
          CALL DEF_SPLITTING2(X_DOMAINS,Y_DOMAINS,X_DIM - X_DIM/NB_PROC_XY -1,  &
               Y_DIM,NB_PROC_XY-1,PREM)
          !
          !     the last Z processor slide hold last all Y dim 
          !
          IDOM  = NB_PROC - KZ_PROC
          DO IK = 1,KZ_PROC
             IDOM = IDOM +  1
             TP => TPROC(IDOM)
             TP%NUMBER = IDOM
             ! X coordonate
             TP%NXOR  = X_DIM - X_DIM/NB_PROC_XY
             TP%NXEND = X_DIM
             ! Y coordonate
             TP%NYOR  = 1
             TP%NYEND = Y_DIM
             ! Z coordonate
             CALL SLIDE_COORD(Z_DIM,KZ_PROC,IK,TP%NZOR,TP%NZEND)
          ENDDO
          !
          ! cartesian splitting with NB_PROC-1
          !
          CALL CARTESIANZ(TPROC,NB_PROC_XY-1,X_DIM-X_DIM/NB_PROC_XY-1,   &
               Y_DIM,Z_DIM,X_DOMAINS,Y_DOMAINS,Z_DOMAINS,213)
          !
       ELSE ! X/P2 Y/P1 & Z splitting --> general case
          !
          CALL CARTESIANZ(TPROC,NB_PROC,X_DIM,Y_DIM,Z_DIM,X_DOMAINS,Y_DOMAINS,Z_DOMAINS,213)
       END IF
    ELSEIF(HSPLITTING.EQ."P1P2SPLITT") THEN
       ! X/P1 Y/P2 Z splitting 
    X_DOMAINS_NEW = Y_DOMAINS
    Y_DOMAINS     = X_DOMAINS
    X_DOMAINS     = X_DOMAINS_NEW
         CALL CARTESIANZ(TPROC,NB_PROC,X_DIM,Y_DIM,Z_DIM,X_DOMAINS,Y_DOMAINS,Z_DOMAINS,123)
    ELSEIF(HSPLITTING.EQ."XSPLITTING") THEN
       !
       X_DOMAINS=NB_PROC_XY
       Y_DOMAINS=1
       CALL CARTESIANZ(TPROC,NB_PROC,X_DIM,Y_DIM,Z_DIM,X_DOMAINS,Y_DOMAINS,Z_DOMAINS,312)
       !
    ELSE ! YSPLITTING
       !
       X_DOMAINS=1
       Y_DOMAINS=NB_PROC_XY
       CALL CARTESIANZ(TPROC,NB_PROC,X_DIM,Y_DIM,Z_DIM,X_DOMAINS,Y_DOMAINS,Z_DOMAINS,321)
       !
    END IF
    IF (      ( (1+TPROC(IP)%NXEND-TPROC(IP)%NXOR) == 0 )  &
         .OR. ( (1+TPROC(IP)%NYEND-TPROC(IP)%NYOR) == 0 )  ) THEN
       PRINT*, "/!\ SPLITZ: some proc have 0 size local domaine , to much processors used for domaine size ", &
            " IP="  ,IP , &
            " DIMX=",1+TPROC(IP)%NXEND-TPROC(IP)%NXOR, &
            " DIMY=",1+TPROC(IP)%NYEND-TPROC(IP)%NYOR
       CALL ABORT()
    END IF
    !
    !*       3.     shift from physical to extended domain
    !
    IF ( IAND(NZ_SPLITTING,2) > 0 ) THEN 
       IF (NZ_VERB .GE. 5 ) THEN
          IF ( IP .EQ. 1 )THEN
             PRINT*,"********************************************************************"
             PRINT*,"******************** HSPLITTING=",HSPLITTING," *************************"
             PRINT*,"********************************************************************"
             PRINT*,"NB_PROC=",NB_PROC," KZ_PROC=",KZ_PROC
             PRINT*,"X_DIM=", X_DIM ,", X_DOMAINS=",X_DOMAINS
             PRINT*,"Y_DIM=", Y_DIM ,", Y_DOMAINS=",Y_DOMAINS
             PRINT*,"Z_DIM=", Z_DIM ,", Z_DOMAINS=",Z_DOMAINS
             
             DO IK = 1,NB_PROC,MAX(1,NB_PROC-1)
                PRINT*," ============== NPROC=",IK,"========================"
                PRINT*,"NXOR=",TPROC(IK)%NXOR," NXEND=",TPROC(IK)%NXEND," TAILLE=",1+TPROC(IK)%NXEND-TPROC(IK)%NXOR
                PRINT*,"NYOR=",TPROC(IK)%NYOR," NYEND=",TPROC(IK)%NYEND," TAILLE=",1+TPROC(IK)%NYEND-TPROC(IK)%NYOR
                PRINT*,"NZOR=",TPROC(IK)%NZOR," NZEND=",TPROC(IK)%NZEND," TAILLE=",1+TPROC(IK)%NZEND-TPROC(IK)%NZOR
             END DO

          ENDIF
       END IF
    END IF
    !    STOP
    !
    ! Add 'Halo points' to global coordonne in X & Y direction
    !
    TPROC(:)%NXOR  = TPROC(:)%NXOR  + JPHEXT
    TPROC(:)%NXEND = TPROC(:)%NXEND + JPHEXT
    !
    TPROC(:)%NYOR  = TPROC(:)%NYOR  + JPHEXT
    TPROC(:)%NYEND = TPROC(:)%NYEND + JPHEXT
    !
    ! In Z direction Halo already intégrated in Z dimension 
    !
!!$    TPROC(:)%NZOR  = TPROC(:)%NZOR  + JPVEXT
!!$    TPROC(:)%NZEND = TPROC(:)%NZEND + JPVEXT
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE SPLITZ


  !     ########################################
  SUBROUTINE INI_PZZ( TP_TSPLITS, TPPZS)
    !     ########################################
    !
    !!****  *INI_PZ* - routine to initialize the physical 2way splitting
    !                  of the TPPROCONF variable
    ! 
    !!    Purpose
    !!    -------
    !     the purpose of this routine is to fill the arguments of the
    !     variable TPPROCONF$TSPLITS_B concerning physical subdomains
    !     with a given splitting TPPZS
    !
    !!**  Method
    !!    ------
    ! 
    !!    External
    !!    --------
    ! 
    !!    Implicit Arguments
    !!    ------------------
    !
    !     Module MODD_STRUCTURE_ll
    !        types PROCONF_ll, ZONE_ll
    !
    !     Module MODD_VAR_ll
    !        NPROC - Number of processors
    !
    !!    Reference
    !!    ---------
    ! 
    !!    Author
    !!    ------
    !     R. Guivarch               * CERFACS - ENSEEIHT *
    !!
    !!    Modifications
    !!    -------------
    !     Original 01/05/98
    !     R. Guivarch 01/08/98 arguments for grid-nesting
    !!
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : PROCONF_ll, ZONE_ll, MODELSPLITTING_ll
    USE MODD_VAR_ll, ONLY       : NPROC
    !
    IMPLICIT NONE
    !
    !
    !*       0.1   declarations of arguments
    !
    TYPE(MODELSPLITTING_ll), DIMENSION(:), INTENT(INOUT):: TP_TSPLITS     ! Physical Zone Splitting
    !
    TYPE(ZONE_ll), DIMENSION(:), INTENT(IN):: TPPZS     ! Physical Zone Splitting
    !
    !
    !*       0.2   declarations of local variables
    !
    INTEGER :: J ! loop control variable
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    FILL TPPROCONF%TSPLITS_B FOR EACH J :
    !              ---------------------------------------
    !
    DO J = 1, NPROC
       !
       TP_TSPLITS(J)%NUMBER = TPPZS(J)%NUMBER
       !
       TP_TSPLITS(J)%NXORP  = TPPZS(J)%NXOR
       TP_TSPLITS(J)%NXENDP = TPPZS(J)%NXEND
       TP_TSPLITS(J)%NDIMXP = TPPZS(J)%NXEND - TPPZS(J)%NXOR + 1
       !
       TP_TSPLITS(J)%NYORP  = TPPZS(J)%NYOR
       TP_TSPLITS(J)%NYENDP = TPPZS(J)%NYEND
       TP_TSPLITS(J)%NDIMYP = TPPZS(J)%NYEND - TPPZS(J)%NYOR + 1
       !
       TP_TSPLITS(J)%NZORP  = TPPZS(J)%NZOR
       TP_TSPLITS(J)%NZENDP = TPPZS(J)%NZEND
       TP_TSPLITS(J)%NDIMZP = TPPZS(J)%NZEND - TPPZS(J)%NZOR + 1
       !
    ENDDO
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE INI_PZZ
  !
  !##############################
  SUBROUTINE INI_EZZ( TPPROCONF )
    !############################
    !
    !!    Author
    !!    ------
    !     J. ESCOBAR              * LA - CNRS *
    ! 
    !!    Modifications
    !!    -------------
    !     Original 22/01/2008
    ! 
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY  : PROCONF_ll
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    TYPE(PROCONF_ll), POINTER :: TPPROCONF       ! splitting data structure
    !
    !-------------------------------------------------------------------------------
    !
    CALL INI_EZZZ( TPPROCONF%TBOUND_SXP1_YP2_Z, TPPROCONF%TSPLITS_SXP1_YP2_Z)
    CALL INI_EZZZ( TPPROCONF%TBOUND_SXP2_YP1_Z, TPPROCONF%TSPLITS_SXP2_YP1_Z)
    CALL INI_EZZZ( TPPROCONF%TBOUND_SX_YP2_ZP1, TPPROCONF%TSPLITS_SX_YP2_ZP1)
    CALL INI_EZZZ( TPPROCONF%TBOUND_SXP2_Y_ZP1, TPPROCONF%TSPLITS_SXP2_Y_ZP1)
    !-------------------------------------------------------------------------------
    !
  CONTAINS
    SUBROUTINE INI_EZZZ(TP_TBOUNDS,TP_TSPLITS)
      !
      !*       0.    DECLARATIONS
      !
      USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT ! Horizontal/Vertical External points number
      USE MODD_STRUCTURE_ll , ONLY : MODELSPLITTING_ll, LOCALISATION_ll
      USE MODD_VAR_ll       , ONLY : NPROC,  JPHALO
      !
      IMPLICIT NONE
      !
      !*       0.1   declarations of arguments
      !
      TYPE(LOCALISATION_ll)   , DIMENSION(:), POINTER :: TP_TBOUNDS
      TYPE(MODELSPLITTING_ll) , DIMENSION(:), POINTER :: TP_TSPLITS
      !
      !*       0.2   declarations of local variables
      !
      INTEGER                 :: J ! loop control variable
      INTEGER                 :: IFLAG
      IFLAG=0
      DO J = 1, NPROC
         !
         ! Origine indexe 
         !
         IF(TP_TBOUNDS(J)%WEST) THEN
            TP_TSPLITS(J)%NXORE = TP_TSPLITS(J)%NXORP - JPHEXT*IFLAG
         ELSE
            TP_TSPLITS(J)%NXORE = TP_TSPLITS(J)%NXORP - JPHALO*IFLAG
         ENDIF
         !
         IF(TP_TBOUNDS(J)%SOUTH) THEN
            TP_TSPLITS(J)%NYORE = TP_TSPLITS(J)%NYORP - JPHEXT*IFLAG
         ELSE
            TP_TSPLITS(J)%NYORE = TP_TSPLITS(J)%NYORP - JPHALO*IFLAG
         ENDIF
         !
         IF(TP_TBOUNDS(J)%BOTTOM) THEN
            TP_TSPLITS(J)%NZORE = TP_TSPLITS(J)%NZORP - JPVEXT*IFLAG
         ELSE
            TP_TSPLITS(J)%NZORE = TP_TSPLITS(J)%NZORP - JPHALO*IFLAG
         ENDIF
         !
         ! End indexe
         !
         IF(TP_TBOUNDS(J)%EAST) THEN
            TP_TSPLITS(J)%NXENDE = TP_TSPLITS(J)%NXENDP + JPHEXT*IFLAG
         ELSE
            TP_TSPLITS(J)%NXENDE = TP_TSPLITS(J)%NXENDP + JPHALO*IFLAG
         ENDIF
         !
         IF(TP_TBOUNDS(J)%NORTH) THEN
            TP_TSPLITS(J)%NYENDE = TP_TSPLITS(J)%NYENDP + JPHEXT*IFLAG
         ELSE
            TP_TSPLITS(J)%NYENDE = TP_TSPLITS(J)%NYENDP + JPHALO*IFLAG
         ENDIF
         !
         IF(TP_TBOUNDS(J)%TOP) THEN
            TP_TSPLITS(J)%NZENDE = TP_TSPLITS(J)%NZENDP + JPVEXT*IFLAG
         ELSE
            TP_TSPLITS(J)%NZENDE = TP_TSPLITS(J)%NZENDP + JPHALO*IFLAG
         ENDIF
         !
         !
         ! Size indexe
         !
         TP_TSPLITS(J)%NDIMXE = TP_TSPLITS(J)%NXENDE - TP_TSPLITS(J)%NXORE + 1
         TP_TSPLITS(J)%NDIMYE = TP_TSPLITS(J)%NYENDE - TP_TSPLITS(J)%NYORE + 1
         TP_TSPLITS(J)%NDIMZE = TP_TSPLITS(J)%NZENDE - TP_TSPLITS(J)%NZORE + 1
         !
      ENDDO
    END SUBROUTINE INI_EZZZ
  END SUBROUTINE INI_EZZ

  !     #######################################
  SUBROUTINE INI_BOUNDARIESZ( TPPROCONF )
    !####################################
    !
    !!****  *INI_BOUNDARIESZ* - routine to compute the localisation variable TBOUND
    !
    !!
    !!    Author
    !!    ------
    !!    J. ESCOBAR
    !!
    !!    Modifications
    !!    -------------
    !!    Original 22/01/2008
    !!
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY  : PROCONF_ll
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    TYPE(PROCONF_ll), POINTER :: TPPROCONF       ! splitting data structure
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    FILL TBOUND OF EACH PROCESSOR :
    !              ------------------------------
    !
    CALL INI_BOUNZ( TPPROCONF%TBOUND_SXP1_YP2_Z, TPPROCONF%TSPLITS_SXP1_YP2_Z)
    CALL INI_BOUNZ( TPPROCONF%TBOUND_SXP2_YP1_Z, TPPROCONF%TSPLITS_SXP2_YP1_Z)
    CALL INI_BOUNZ( TPPROCONF%TBOUND_SX_YP2_ZP1, TPPROCONF%TSPLITS_SX_YP2_ZP1)
    CALL INI_BOUNZ( TPPROCONF%TBOUND_SXP2_Y_ZP1, TPPROCONF%TSPLITS_SXP2_Y_ZP1)
    !
    !-------------------------------------------------------------------------------
    !
  CONTAINS
    SUBROUTINE INI_BOUNZ(TP_TBOUNDS,TP_TSPLITS)
      !
      !*       0.    DECLARATIONS
      !
      USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT ! Horizontal/Vertical External points number
      USE MODD_STRUCTURE_ll , ONLY : MODELSPLITTING_ll, LOCALISATION_ll
      USE MODD_VAR_ll       , ONLY : NPROC, DIMX, DIMY, DIMZ
      !
      IMPLICIT NONE
      !
      !*       0.1   declarations of arguments
      !
      TYPE(LOCALISATION_ll)   , DIMENSION(:), POINTER :: TP_TBOUNDS
      TYPE(MODELSPLITTING_ll) , DIMENSION(:), POINTER :: TP_TSPLITS
      !
      !*       0.2   declarations of local variables
      !
      INTEGER                 :: J ! loop control variable
      !JUANZ
      INTEGER                 :: JNUMBER
      !JUANZ
      !
      ALLOCATE(TP_TBOUNDS(NPROC))
      !
      DO J = 1, NPROC
         !
         JNUMBER =  TP_TSPLITS(J)%NUMBER
         !
         TP_TBOUNDS(JNUMBER) = &
              LOCALISATION_ll( .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE.)
         !
         IF( TP_TSPLITS(J)%NXORP .EQ. (1+JPHEXT)) &
              TP_TBOUNDS(JNUMBER)%WEST = .TRUE.
         !
         IF( TP_TSPLITS(J)%NXENDP  .EQ. (DIMX-JPHEXT)) &
              TP_TBOUNDS(JNUMBER)%EAST = .TRUE.
         !
         IF( TP_TSPLITS(J)%NYORP .EQ. (1+JPHEXT)) &
              TP_TBOUNDS(JNUMBER)%SOUTH = .TRUE.
         !
         IF( TP_TSPLITS(J)%NYENDP  .EQ. (DIMY-JPHEXT)) &
              TP_TBOUNDS(JNUMBER)%NORTH = .TRUE.
         !
         IF( TP_TSPLITS(J)%NZORP .EQ. (1+JPVEXT)) &
              TP_TBOUNDS(JNUMBER)%BOTTOM = .TRUE.
         !
         IF( TP_TSPLITS(J)%NZENDP  .EQ. (DIMZ-JPVEXT)) &
              TP_TBOUNDS(JNUMBER)%TOP = .TRUE.
      ENDDO
    END SUBROUTINE INI_BOUNZ
  END SUBROUTINE INI_BOUNDARIESZ

  !
  !     ##################################################
  SUBROUTINE CONSTRUCT_TRANSZ(TPCOMDATA, TPPROCONF )
    !     ##################################################
    !
    !!****  *CONSTRUCT_TRANS* - routine to construct
    !                           the transposition correspondants
    !
    !!    Purpose
    !!    -------
    !     the purpose of the routine is to fill the structured type variable
    !     TPCOMDATA with informations concerning the communications
    !     during transposition
    !
    !!**  Method
    !!    ------
    !     we compute for the local processor,
    !      - intersections between zones of the global domain in x-slice splitting
    !        and local zone in 2way splitting to find the send correspondant
    !        of the local processor for transposition 2way / x-slices
    ! 
    !      - intersections between physical zones of the global domain
    !        in 2way splitting and local zone in x-slices splitting
    !        to find the receive correspondant of the local processor
    !        for transposition 2way / x-slices
    ! 
    !      - intersections between zones of the global domain in y-slice splitting
    !        and local zone in x-slices splitting to find the send correspondant
    !        of the local processor for transposition x-slices / y-slices
    ! 
    !      - intersections between zones of the global domain
    !        in x-slices splitting and local zone in y-slices splitting
    !        to find the receive correspondant of the local processor
    !        for transposition x-slices / y-slices
    ! 
    !!    External
    !!    --------
    !     Module MODE_TOOLS_ll
    !        ADD_ZONE, INTERSECTION, GLOBAL2LOCAL, EXTRACT_ZONE, G2LX
    !
    !     Module MODE_CONSTRUCT_ll
    !        RESIZE_TRANS
    !
    !!    Implicit Arguments
    !!    ------------------
    !
    !     Module MODD_STRUCTURE_ll
    !        types ZONE_ll, PROC_COM_DATA_ll, PROCONF_ll
    !
    !     Module MODD_VAR_ll
    !        IP - Number of local processor=subdomain
    !        NPROC - Number of processors
    !
    !!    Reference
    !!    ---------
    ! 
    !!    Author
    !!    ------
    !     R. Guivarch               * CERFACS - ENSEEIHT *
    ! 
    !!    Modifications
    !!    -------------
    !     Original 01/05/98
    !     R. Guivarch 01/08/98 arguments for grid-nesting
    !!
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : ZONE_ll, PROC_COM_DATA_ll, PROCONF_ll
    USE MODD_VAR_ll, ONLY       : IP, NPROC , NP1 , NP2 , NMNH_ROW_WORLD , NMNH_COL_WORLD
    !
    USE MODE_TOOLS_ll, ONLY     : INTERSECTION, G2LX, GLOBAL2LOCAL, ADD_ZONE, &
         EXTRACT_ZONE
    !
    IMPLICIT NONE
    !
    !
    !*       0.1   declarations of arguments
    !
    TYPE(PROC_COM_DATA_ll), POINTER :: TPCOMDATA ! communications data
    TYPE(PROCONF_ll), POINTER       :: TPPROCONF ! splitting data
    !
    !*       0.2   declarations of local variables
    !
    INTEGER, PARAMETER              :: SEND_DIR = 0 , RECV_DIR = 1
    !
    !-------------------------------------------------------------------------------
    !
    !        1.1    CONSTRUCTION OF TRANSPOSITION B -> SX_YP2_ZP1-SLICES CORRESPONDANTS :
    !              -------------------------------------------------------------
    !
!!$    CALL CONSTRUCT_TRANSZZ(NPROC,1000,TPCOMDATA%TSEND_B_SX_YP2_ZP1,SEND_DIR,&
!!$         TPCOMDATA%TSEND_BOX_B_SX_YP2_ZP1,&
!!$         TPPROCONF%TSPLITS_B,TPPROCONF%TSPLITS_SX_YP2_ZP1,NMNH_ROW_WORLD,NP1)
!!$    !
!!$    CALL CONSTRUCT_TRANSZZ(NPROC,1000,TPCOMDATA%TRECV_B_SX_YP2_ZP1,RECV_DIR,&
!!$         TPCOMDATA%TRECV_BOX_B_SX_YP2_ZP1,&     
!!$         TPPROCONF%TSPLITS_SX_YP2_ZP1,TPPROCONF%TSPLITS_B,NMNH_ROW_WORLD,NP1)
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.    CONSTRUCTION OF TRANSPOSITION SXP1_YP2_Z -> SX_YP2_ZP1-SLICES CORRESPONDANTS :
    !              -------------------------------------------------------------
    !
    CALL CONSTRUCT_TRANSZZ(NPROC,3000,TPCOMDATA%TSEND_SXP1_YP2_Z_SX_YP2_ZP1,SEND_DIR,&
         TPCOMDATA%TSEND_BOX_SXP1_YP2_Z_SX_YP2_ZP1,&
         TPPROCONF%TSPLITS_SXP1_YP2_Z,TPPROCONF%TSPLITS_SX_YP2_ZP1,NMNH_ROW_WORLD,NP1)
    !
    CALL CONSTRUCT_TRANSZZ(NPROC,3000,TPCOMDATA%TRECV_SXP1_YP2_Z_SX_YP2_ZP1,RECV_DIR,&
         TPCOMDATA%TRECV_BOX_SXP1_YP2_Z_SX_YP2_ZP1,&
         TPPROCONF%TSPLITS_SX_YP2_ZP1,TPPROCONF%TSPLITS_SXP1_YP2_Z,NMNH_ROW_WORLD,NP1)
    !
    !-------------------------------------------------------------------------------
    !
    !*       3.    CONSTRUCTION OF
    !              TRANSPOSITION SX_YP2_ZP1-SLICES -> SXP2_Y_ZP1-SLICES CORRESPONDANTS :
    !              -------------------------------------------------
    !
    CALL CONSTRUCT_TRANSZZ(NPROC,4000,TPCOMDATA%TSEND_SX_YP2_ZP1_SXP2_Y_ZP1,SEND_DIR,&
         TPCOMDATA%TSEND_BOX_SX_YP2_ZP1_SXP2_Y_ZP1,&
         TPPROCONF%TSPLITS_SX_YP2_ZP1,TPPROCONF%TSPLITS_SXP2_Y_ZP1,NMNH_COL_WORLD,NP2)
    !
    CALL CONSTRUCT_TRANSZZ(NPROC,4000,TPCOMDATA%TRECV_SX_YP2_ZP1_SXP2_Y_ZP1,RECV_DIR,&
         TPCOMDATA%TRECV_BOX_SX_YP2_ZP1_SXP2_Y_ZP1,&
         TPPROCONF%TSPLITS_SXP2_Y_ZP1,TPPROCONF%TSPLITS_SX_YP2_ZP1,NMNH_COL_WORLD,NP2)
    !
    !-------------------------------------------------------------------------------
    !
    !*       4.    CONSTRUCTION OF
    !              TRANSPOSITION SXP2_Y_ZP1-SLICES -> SXP2_YP1_Z-SLICES CORRESPONDANTS :
    !              -------------------------------------------------
    !
    CALL CONSTRUCT_TRANSZZ(NPROC,5000,TPCOMDATA%TSEND_SXP2_Y_ZP1_SXP2_YP1_Z,SEND_DIR,&
         TPCOMDATA%TSEND_BOX_SXP2_Y_ZP1_SXP2_YP1_Z,&
         TPPROCONF%TSPLITS_SXP2_Y_ZP1,TPPROCONF%TSPLITS_SXP2_YP1_Z,NMNH_ROW_WORLD,NP1)
    !
    CALL CONSTRUCT_TRANSZZ(NPROC,5000,TPCOMDATA%TRECV_SXP2_Y_ZP1_SXP2_YP1_Z,RECV_DIR,&
         TPCOMDATA%TRECV_BOX_SXP2_Y_ZP1_SXP2_YP1_Z,&
         TPPROCONF%TSPLITS_SXP2_YP1_Z,TPPROCONF%TSPLITS_SXP2_Y_ZP1,NMNH_ROW_WORLD,NP1)
    !
    !*       3.1   extraction of x-slices splitting
    !
    !
    !
    !*       3.5   Switch to local coordinates
    !
    !CALL G2LX(TPPROCONF%TSPLITS_X(IP),TPCOMDATA%TSEND_TRANS_XY)
    !CALL G2LX(TPPROCONF%TSPLITS_Y(IP),TPCOMDATA%TRECV_TRANS_XY)
    !
    !
    !-------------------------------------------------------------------------------
    !
    !*       4.    DEALLOCATION OF LOCAL VARIABLES :
    !              -------------------------------
    !    
    !
    !-------------------------------------------------------------------------------
    !
  CONTAINS
    SUBROUTINE CONSTRUCT_TRANSZZ(KPROC,KMESG,TP_TRANS_FROM_TO,KDIR,&
         TBOX_FROM_TO,&
         TP_TSPLITS_FROM,TP_TSPLITS_TO,KDIM_COM,KDIM_SIZE)
      !
      USE MODD_STRUCTURE_ll , ONLY : CRSPD_ll, MODELSPLITTING_ll, BOX_ll
      USE MODD_VAR_ll       , ONLY : IP
      USE MODD_CONFZ        , ONLY : NZ_VERB
      IMPLICIT NONE
      !
      ! Argument
      !
      INTEGER                                         :: KPROC,KMESG
      TYPE(CRSPD_ll)                        , POINTER :: TP_TRANS_FROM_TO
      INTEGER                                         :: KDIR
      TYPE(BOX_ll)                          , POINTER :: TBOX_FROM_TO
      TYPE(MODELSPLITTING_ll),  DIMENSION(:), POINTER :: TP_TSPLITS_FROM,TP_TSPLITS_TO
      INTEGER                                         :: KDIM_COM,KDIM_SIZE
      !
      ! Local Variable
      !
      TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:)  :: TP_ZONE_FROM,TP_ZONE_TO,TZINTER
      INTEGER                                   :: JI,JII,NSEND,JB,JE,JINC,ITIC
      INTEGER                , DIMENSION(KPROC) :: ISEND
      !-------------------------------------------------------------------------------
      !
      !*       1.    ALLOCATE OF THE LOCAL VARIABLES :
      !              -------------------------------
      !
      ALLOCATE(TP_ZONE_FROM(KPROC),TP_ZONE_TO(KPROC),TZINTER(KPROC))
      !
      !-------------------------------------------------------------------------------
      !
      !*       2.    CONSTRUCTION OF TRANSPOSITION 2WAY -> X-SLICES CORRESPONDANTS :
      !              -------------------------------------------------------------
      !
      !*       2.1   extraction of physical 2way splitting
      !              from TPPROCONF structure
      !
      CALL EXTRACT_ZONEZ(TP_TSPLITS_FROM, TP_ZONE_FROM, TZINTER)
      !
      !*       2.2   extraction of x-slices splitting
      !              from TPPROCONF structure
      !
      CALL EXTRACT_ZONEZ(TP_TSPLITS_TO, TP_ZONE_TO, TZINTER)
      !
      !*       2.3   computation of intersection between local 2way zone
      !              and x-slices splitting -> send correspondants
      !
      ! find the from zone corresponding to this proc=IP
      DO  JI = 1, KPROC
         if ( TP_ZONE_FROM(JI)%NUMBER .EQ. IP ) EXIT
      END DO
      CALL INTERSECTIONZ(TP_ZONE_TO, KPROC, TP_ZONE_FROM(JI), TZINTER,&
                         TBOX_FROM_TO,KDIM_COM,KDIM_SIZE)
      !
      NULLIFY(TP_TRANS_FROM_TO)

      NSEND = 0
      ITIC = -1
      DO JI = 1, KPROC
         !JUAN add shift from IP to distribue send/recv
         ITIC = - ITIC
         JII= 1+ MOD(KPROC+ITIC*(JI/2)+IP-1,KPROC)
         IF(TZINTER(JII)%NUMBER.NE.0) THEN
	   NSEND = NSEND + 1
           ISEND(NSEND)=JII	   
	 ENDIF
      ENDDO
      IF (KDIR .EQ. SEND_DIR) THEN
      JB = NSEND ; JE = 1 ; JINC = -1
      IF (NZ_VERB .GE. 5 ) THEN
         IF (IP .EQ. (KPROC+1)/2 ) print*,"mode_toolsz_ll.f90:: SEND IP=",IP,":",(ISEND(JI),JI=JB,JE,JINC)
      END IF
      ELSE
      JB = 1 ; JE = NSEND ; JINC =  1
      IF (NZ_VERB .GE. 5 ) THEN
         IF (IP .EQ. (KPROC+1)/2 ) print*,"mode_toolsz_ll.f90:: RECV IP=",IP,":",(ISEND(JI),JI=JB,JE,JINC)
      END IF
      ENDIF
      DO JI = JB,JE,JINC
         !JUAN add shift from IP to distribue send/recv
         JII= ISEND(JI)
         IF(TZINTER(JII)%NUMBER.NE.0) THEN
            TZINTER(JII)%MSSGTAG = KMESG 
            CALL ADD_ZONE(TP_TRANS_FROM_TO, TZINTER(JII) )
!!$         IF(TZINTER(JI)%NUMBER.NE.0) THEN
!!$            !TZINTER(JI)%MSSGTAG = KMESG 
!!$            TZINTER(JI)%MSSGTAG = KMESG + JI*0000 + TZINTER(JI)%NUMBER*0
!!$            CALL ADD_ZONE(TP_TRANS_FROM_TO, TZINTER(JI) )
         ENDIF
      ENDDO
      !
      ! find the from zone corresponding to this proc=IP
      DO  JI = 1, KPROC
         if ( TP_TSPLITS_FROM(JI)%NUMBER .EQ. IP ) EXIT
      END DO
      CALL G2LXZ(TP_TSPLITS_FROM(JI),TP_TRANS_FROM_TO,TBOX_FROM_TO)
      !
      DEALLOCATE(TP_ZONE_FROM,TP_ZONE_TO,TZINTER)
      !
    END SUBROUTINE CONSTRUCT_TRANSZZ
  END SUBROUTINE CONSTRUCT_TRANSZ



  !     #################################################
  SUBROUTINE EXTRACT_ZONEZ( TPSPLITS, TPPZS, TPEZS )
    !     #################################################
    !
    !!****  *EXTRACT_ZONEZ* - routine to construct two splittings variables
    !!                       from a MODELSPLITTING_ll variable
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR   LA - CNRS
    ! 
    !!    Modifications
    !!    -------------
    !     Original 28/01/2008
    ! 
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : MODELSPLITTING_ll, ZONE_ll
    USE MODD_VAR_ll, ONLY : NPROC
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    TYPE(MODELSPLITTING_ll), DIMENSION(:), POINTER :: TPSPLITS
    !
    TYPE(ZONE_ll), DIMENSION(:), INTENT(OUT) :: TPPZS, TPEZS
    !
    !*       0.2   declarations of local variables
    !
    INTEGER :: J ! loop control variable
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    FILL TPPZS AND TPEZS FOR EACH J :
    !              -------------------------------
    !
    DO J = 1, NPROC
       !
       TPPZS(J) = ZONE_ll( 0, 0, 0, 0, 0, 0, 0, 0 )
       TPEZS(J) = ZONE_ll( 0, 0, 0, 0, 0, 0, 0, 0 )
       !
       TPPZS(J)%NUMBER = TPSPLITS(J)%NUMBER
       TPPZS(J)%NXOR   = TPSPLITS(J)%NXORP
       TPPZS(J)%NYOR   = TPSPLITS(J)%NYORP
       TPPZS(J)%NXEND  = TPSPLITS(J)%NXENDP
       TPPZS(J)%NYEND  = TPSPLITS(J)%NYENDP
       !JUAN
!!$       TPPZS(J)%NZOR   = TPSPLITS(J)%NZORP
!!$       TPPZS(J)%NZEND  = TPSPLITS(J)%NZENDP
       ! JUAN :: For Z splitting extended domain is needed in Z slide
       TPPZS(J)%NZOR   = TPSPLITS(J)%NZORE
       TPPZS(J)%NZEND  = TPSPLITS(J)%NZENDE
       !JUAN
       !
       TPEZS(J)%NUMBER = TPSPLITS(J)%NUMBER
       TPEZS(J)%NXOR   = TPSPLITS(J)%NXORE
       TPEZS(J)%NYOR   = TPSPLITS(J)%NYORE
       TPEZS(J)%NXEND  = TPSPLITS(J)%NXENDE
       TPEZS(J)%NYEND  = TPSPLITS(J)%NYENDE
       !JUAN
       TPEZS(J)%NZOR   = TPSPLITS(J)%NZORE
       TPEZS(J)%NZEND  = TPSPLITS(J)%NZENDE
       !JUAN
    ENDDO
    !
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE EXTRACT_ZONEZ

  !     ################################
  SUBROUTINE G2LXZ(TPSPLIT,TPCRSPD,TBOX)
    !     ################################
    !
    !!****  *G2LXZ* - routine to switch from global coordinates to local ones
    ! 
    !
    !!    Reference
    !!    ---------
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR
    !
    !!    Modifications
    !!    -------------
    !     Original 01/02/2008
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : MODELSPLITTING_ll, CRSPD_ll, ZONE_ll, BOX_ll
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    TYPE(MODELSPLITTING_ll), INTENT(IN) :: TPSPLIT ! x-slices or y-slices
    ! splitting
    !
    TYPE(CRSPD_ll), POINTER :: TPCRSPD ! CRSPD_ll to be switch
    TYPE(BOX_ll)  , POINTER :: TBOX
    !
    !*       0.2   declarations of local variables
    !
    ! intermediate variables to describe the list and the elements
    TYPE(ZONE_ll), POINTER :: TZZONE
    TYPE(CRSPD_ll), POINTER :: TZCRSPD
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    SWITCH :
    !              ------
    ! we list the variable TPCRSPD of type CRSPD_ll
    TZCRSPD => TPCRSPD
    DO WHILE (ASSOCIATED(TZCRSPD))
       TZZONE => TZCRSPD%TELT
       !
       !   we substract the origin of the local subdomain (physical subdomain)
       TZZONE%NXOR = TZZONE%NXOR - TPSPLIT%NXORP + 1
       TZZONE%NXEND = TZZONE%NXEND - TPSPLIT%NXORP + 1
       TZZONE%NYOR = TZZONE%NYOR - TPSPLIT%NYORP + 1
       TZZONE%NYEND = TZZONE%NYEND - TPSPLIT%NYORP + 1
       !JUAN
       TZZONE%NZOR = TZZONE%NZOR - TPSPLIT%NZORP + 1
       TZZONE%NZEND = TZZONE%NZEND - TPSPLIT%NZORP + 1
!!$       !JUAN
!!$       !JUAN
!!$       !   we substract the origin of the local subdomain (extended subdomain)
!!$       TZZONE%NXOR = TZZONE%NXOR - TPSPLIT%NXORE + 1
!!$       TZZONE%NXEND = TZZONE%NXEND - TPSPLIT%NXORE + 1
!!$       !
!!$       TZZONE%NYOR = TZZONE%NYOR - TPSPLIT%NYORE + 1
!!$       TZZONE%NYEND = TZZONE%NYEND - TPSPLIT%NYORE + 1
!!$       !JUAN
!!$       TZZONE%NZOR = TZZONE%NZOR - TPSPLIT%NZORE + 1
!!$       TZZONE%NZEND = TZZONE%NZEND - TPSPLIT%NZORE + 1
       !JUAN
       !
       TZCRSPD => TZCRSPD%TNEXT
    ENDDO
    TBOX%NXOR(:)  =  TBOX%NXOR(:)  -  TPSPLIT%NXORE + 1
    TBOX%NXEND(:) =  TBOX%NXEND(:) -  TPSPLIT%NXORE + 1
    !
    TBOX%NYOR(:)  =  TBOX%NYOR(:)  -  TPSPLIT%NYORE + 1
    TBOX%NYEND(:) =  TBOX%NYEND(:) -  TPSPLIT%NYORE + 1
    !
    TBOX%NZOR(:)  =  TBOX%NZOR(:)  -  TPSPLIT%NZORE + 1
    TBOX%NZEND(:) =  TBOX%NZEND(:) -  TPSPLIT%NZORE + 1
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE G2LXZ

  !     ####################################################
  SUBROUTINE INTERSECTIONZ( TPSPLIT, K, TPZONE, TPRES,&
                            TBOX_FROM_TO,KDIM_COM,KDIM_SIZE)
    !     ####################################################
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR
    !!
    !!    Modifications
    !!    -------------
    !     Original 23/01/2008
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : ZONE_ll, BOX_ll, ALLOC_BOX_LL
    USE MODD_VAR_ll, ONLY : DIMZ, NMNH_COMM_WORLD
    !
    IMPLICIT NONE
    !
    !
    !*       0.1   declarations of arguments
    !
    TYPE(ZONE_ll), DIMENSION(:), INTENT(IN) :: TPSPLIT ! Splitting of the domain
    !
    INTEGER, INTENT(IN) :: K ! Number of elements of TPSPLIT
    !
    TYPE(ZONE_ll), INTENT(IN) :: TPZONE ! Zone to be splitted
    !
    TYPE(ZONE_ll), DIMENSION(:), INTENT(OUT) :: TPRES ! Splitting of the zone
    !
     TYPE(BOX_ll)              , POINTER     :: TBOX_FROM_TO
    INTEGER                                 :: KDIM_COM , KDIM_SIZE
    !*       0.2   declarations of local variables
    !
    INTEGER :: J      ! loop control variable
    INTEGER :: JCNTS  ! Size of the current box
    INTEGER :: JSTRT  ! displacement in the pack box 
    INTEGER :: JBOX
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    LIST AND COMPUTE INTERSECTION BETWEEN TPSPLIT(J) AND TPZONE :
    !              -----------------------------------------------------------
    !
    
    CALL ALLOC_BOX_ll(TBOX_FROM_TO,KDIM_SIZE) 
    TBOX_FROM_TO%NCOM = KDIM_COM
    TBOX_FROM_TO%NBOX = KDIM_SIZE
!!$    CALL ALLOC_BOX_ll(TBOX_FROM_TO,K) 
!!$    TBOX_FROM_TO%NCOM = NMNH_COMM_WORLD
!!$    TBOX_FROM_TO%NBOX = K
    !
    JSTRT = 0
    JBOX  = 0
    !
    DO J = 1, K
       ! 
       !   Which subdomain is the owner of TPSPLIT(J)
       TPRES(J)%NUMBER = TPSPLIT(J)%NUMBER
       !
       !   Computation of the origin coordinate
       TPRES(J)%NXOR = MAX( TPZONE%NXOR, TPSPLIT(J)%NXOR )
       TPRES(J)%NYOR = MAX( TPZONE%NYOR, TPSPLIT(J)%NYOR )
       TPRES(J)%NZOR = MAX( TPZONE%NZOR, TPSPLIT(J)%NZOR )
       !
       !   Computation of the last coordinate
       TPRES(J)%NXEND = MIN( TPZONE%NXEND, TPSPLIT(J)%NXEND )
       TPRES(J)%NYEND = MIN( TPZONE%NYEND, TPSPLIT(J)%NYEND )
       TPRES(J)%NZEND = MIN( TPZONE%NZEND, TPSPLIT(J)%NZEND )
       !
       !   if the intersection is void, the result is nullified
       IF(  (TPRES(J)%NXEND - TPRES(J)%NXOR  < 0 ) .OR. &
            (TPRES(J)%NYEND - TPRES(J)%NYOR  < 0 ) .OR. &
            (TPRES(J)%NZEND - TPRES(J)%NZOR  < 0 ) ) THEN
          
          TPRES(J) = ZONE_ll ( 0, 0, 0, 0, 0, 0, 0, 0 )
          JCNTS    = 0
       ELSE                    
          JCNTS = ( TPRES(J)%NXEND - TPRES(J)%NXOR + 1 ) * &
                  ( TPRES(J)%NYEND - TPRES(J)%NYOR + 1 ) * &     
                  ( TPRES(J)%NZEND - TPRES(J)%NZOR + 1 )
!!$       ENDIF
          JBOX = JBOX + 1 
       ! 
       ! Fill the box_ll for mpialltoallv
       !
          TBOX_FROM_TO%NXOR(JBOX)  = TPRES(J)%NXOR
          TBOX_FROM_TO%NXEND(JBOX) = TPRES(J)%NXEND
       !
          TBOX_FROM_TO%NYOR(JBOX)  = TPRES(J)%NYOR
          TBOX_FROM_TO%NYEND(JBOX) = TPRES(J)%NYEND
       !
          TBOX_FROM_TO%NZOR(JBOX)  = TPRES(J)%NZOR
          TBOX_FROM_TO%NZEND(JBOX) = TPRES(J)%NZEND
       !
          TBOX_FROM_TO%NCNT(JBOX)  = JCNTS
          TBOX_FROM_TO%NSTRT(JBOX) = JSTRT
          TBOX_FROM_TO%NUMBER(JBOX) = TPRES(J)%NUMBER
       !
       JSTRT = JSTRT + JCNTS ! next box start pointer
       END IF
       !
    ENDDO
    TBOX_FROM_TO%NSIZE = JSTRT
    
  
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE INTERSECTIONZ

END MODULE MODE_TOOLSZ_ll
