!MNH_LIC Copyright 1996-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #################################################################
      SUBROUTINE WRITE_DIACHRO(TPDIAFILE,TPLUOUTDIA,HGROUP,HTYPE,     &
      KGRID,PDATIME,PVAR,PTRAJT,                                     &
      HTITRE,HUNITE,HCOMMENT,OICP,OJCP,OKCP,KIL,KIH,KJL,KJH,KKL,KKH, &
      PTRAJX,PTRAJY,PTRAJZ,PMASK)
!     #################################################################
!
!!****  *WRITE_DIACHRO* - Ecriture d'un enregistrement dans un fichier
!!                        diachronique (de nom de base HGROUP)
!!
!!    PURPOSE
!!    -------
!      
!
!!**  METHOD
!!    ------
!!      En fait pour un groupe donne HGROUP, on ecrit systematiquement
!       plusieurs enregistrements :
!       - 1: HGROUP.TYPE          (type d'informations a enregistrer)
!       - 2: HGROUP.DIM           (dimensions de toutes les matrices a 
!                                  enregistrer)
!       - 3: HGROUP.TITRE         (Nom des processus)
!       - 4: HGROUP.UNITE         (Unites pour chaque processus)
!       - 5: HGROUP.COMMENT       (Champ commentaire pour chaque processus)
!       - 6: HGROUP.TRAJT         (Temps)
!       - 7: HGROUP.PROCx         (Champ traite . 1 enr./ 1 processus)
!       - 8: HGROUP.DATIM         (Les differentes dates du modele)
!       et pour certains types d'informations on enregistre egalement
!       des coordonnees (HGROUP.TRAJX, HGROUP.TRAJY, HGROUP.TRAJZ)
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       08/01/96
!!      Updated   PM
!!      Modification (N. Asencio) 18/06/99  : the two first dimensions of PMASK
!!                   are linked to the horizontal grid, FMWRIT is called with 'XY' argument. 
!!                   In standard configuration of the budgets, the mask is written once 
!!                   outside this routine with FMWRIT call. Its record name is 'MASK_nnnn.MASK'
!!                   So optional PMASK is not used .
!!      Modification (J. Duron)   24/06/99  : add logical GPACK to disable the pack option,
!!                                            add the initialization of the dimensions of
!!                                          MASK array in MASK case with write outside the 
!!                                          routine.
!!      J.Escobar       02/10/2015 modif for JPHEXT(JPVEXT) variable
!!      D.Gazen+ G.Delautier 06/2016 modif for ncl files
!!      P. Wautelet     09/06/2017: name of the variable added to the name of the written field
!!                                  and better comment (true comment + units)
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BUDGET
USE MODD_CONF
USE MODE_FIELD
USE MODD_IO_ll,      ONLY : TFILEDATA
USE MODD_PARAMETERS, ONLY : JPHEXT
!
USE MODE_ll
USE MODE_FMWRIT
USE MODE_IO_ll
!
USE MODI_MENU_DIACHRO
!
IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------
TYPE(TFILEDATA),              INTENT(IN)          :: TPDIAFILE    ! file to write
TYPE(TFILEDATA),              INTENT(IN)          :: TPLUOUTDIA
CHARACTER(LEN=*),             INTENT(IN)          :: HGROUP, HTYPE
INTEGER,DIMENSION(:),         INTENT(IN)          :: KGRID
REAL,DIMENSION(:,:),          INTENT(IN)          :: PDATIME
REAL,DIMENSION(:,:,:,:,:,:),  INTENT(IN)          :: PVAR
REAL,DIMENSION(:,:),          INTENT(IN)          :: PTRAJT
CHARACTER(LEN=*),DIMENSION(:),INTENT(IN)          :: HTITRE, HUNITE, HCOMMENT
LOGICAL,                      INTENT(IN),OPTIONAL :: OICP, OJCP, OKCP
INTEGER,                      INTENT(IN),OPTIONAL :: KIL, KIH
INTEGER,                      INTENT(IN),OPTIONAL :: KJL, KJH
INTEGER,                      INTENT(IN),OPTIONAL :: KKL, KKH
REAL,DIMENSION(:,:,:),        INTENT(IN),OPTIONAL :: PTRAJX
REAL,DIMENSION(:,:,:),        INTENT(IN),OPTIONAL :: PTRAJY
REAL,DIMENSION(:,:,:),        INTENT(IN),OPTIONAL :: PTRAJZ
REAL,DIMENSION(:,:,:,:,:,:),  INTENT(IN),OPTIONAL :: PMASK
!
!*       0.1   Local variables
!              ---------------
CHARACTER(LEN=20) :: YCOMMENT
CHARACTER(LEN=3)  :: YJ
INTEGER   ::   ILENG, ILENTITRE, ILENUNITE, ILENCOMMENT
INTEGER   ::   ILUOUTDIA
INTEGER   ::   II, IJ, IK, IT, IN, IP, J, JJ
INTEGER   ::   INTRAJT, IKTRAJX, IKTRAJY, IKTRAJZ
INTEGER   ::   ITTRAJX, ITTRAJY, ITTRAJZ
INTEGER   ::   INTRAJX, INTRAJY, INTRAJZ
INTEGER   ::   IIMASK, IJMASK, IKMASK, ITMASK, INMASK, IPMASK
INTEGER   ::   ICOMPX, ICOMPY, ICOMPZ
INTEGER   ::   IIMAX_ll, IJMAX_ll ! size of the physical global domain
INTEGER,DIMENSION(:),ALLOCATABLE :: ITABCHAR
LOGICAL   ::   GPACK
TYPE(TFIELDDATA)  :: TZFIELD
!------------------------------------------------------------------------------
!
GPACK=LPACK
LPACK=.FALSE.
YCOMMENT='NOTHING'
!
ILUOUTDIA = TPLUOUTDIA%NLU
!
! BUG ...ca passe que si PRESENT(OICP) sinon OICP non defini 
! Question: doit-on mettre condition comme:
!  IF(HTYPE == 'CART' .AND. .NOT. PRESENT(OICP) .AND. .NOT. PRESENT(OJCP)) THEN

! en attendant correction on debranche avec un IF Present. ENDIF av
! RETURN
IF (PRESENT(OICP) .AND. PRESENT(OJCP)) THEN
  IF(HTYPE == 'CART' .AND. .NOT. OICP .AND. .NOT. OJCP) THEN
                              !for parallel execution, PVAR is distributed on several proc
    II=KIH-KIL+1
    IJ=KJH-KJL+1
  ELSE
    II = SIZE(PVAR,1)
    IJ = SIZE(PVAR,2)
  ENDIF
ELSE
    II = SIZE(PVAR,1)
    IJ = SIZE(PVAR,2)

ENDIF
IK = SIZE(PVAR,3)
IT = SIZE(PVAR,4)
IN = SIZE(PVAR,5)
IP = SIZE(PVAR,6)

INTRAJT=SIZE(PTRAJT,2)

IKTRAJX=0; IKTRAJY=0; IKTRAJZ=0
ITTRAJX=0; ITTRAJY=0; ITTRAJZ=0
INTRAJX=0; INTRAJY=0; INTRAJZ=0
IF(PRESENT(PTRAJX))THEN
  IKTRAJX=SIZE(PTRAJX,1)
  ITTRAJX=SIZE(PTRAJX,2)
  INTRAJX=SIZE(PTRAJX,3)
ENDIF
IF(PRESENT(PTRAJY))THEN
  IKTRAJY=SIZE(PTRAJY,1)
  ITTRAJY=SIZE(PTRAJY,2)
  INTRAJY=SIZE(PTRAJY,3)
ENDIF
IF(PRESENT(PTRAJZ))THEN
  IKTRAJZ=SIZE(PTRAJZ,1)
  ITTRAJZ=SIZE(PTRAJZ,2)
  INTRAJZ=SIZE(PTRAJZ,3)
ENDIF

IIMASK=0; IJMASK=0; IKMASK=0; ITMASK=0; INMASK=0; IPMASK=0
IF(HTYPE == 'MASK')THEN
!     MASK is written outside this routine but the dimensions must be initialized
!     the mask is defined on the extended domain
  CALL GET_GLOBALDIMS_ll (IIMAX_ll,IJMAX_ll)
  IIMASK=IIMAX_ll + 2 * JPHEXT
  IJMASK=IJMAX_ll + 2 * JPHEXT
  IF(PRESENT(PMASK))THEN
    IKMASK=SIZE(PMASK,3)
    ITMASK=SIZE(PMASK,4)
    INMASK=SIZE(PMASK,5)
    IPMASK=SIZE(PMASK,6)
  ELSE
    IKMASK=1
    ITMASK=NBUWRNB
    INMASK=NBUMASK
    IPMASK=1
  ENDIF
ENDIF

ILENTITRE = LEN(HTITRE)
ILENUNITE = LEN(HUNITE)
ILENCOMMENT = LEN(HCOMMENT)

ICOMPX=0; ICOMPY=0; ICOMPZ=0
IF(PRESENT(OICP))THEN
IF(OICP)THEN
  ICOMPX=1
ENDIF
IF(OJCP)THEN
  ICOMPY=1
ENDIF
IF(OKCP)THEN
  ICOMPZ=1
ENDIF
ENDIF
!
IF (NVERB>=5) THEN
  WRITE(ILUOUTDIA,*)' WRITE_DIACHRO: ',TRIM(TPDIAFILE%CNAME)//'.lfi'
ENDIF
!
! 1er enregistrement TYPE
!
TZFIELD%CMNHNAME   = TRIM(HGROUP)//'.TYPE'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = TRIM(HGROUP)//'.TYPE'
TZFIELD%CUNITS     = ''
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = TRIM(YCOMMENT)
TZFIELD%NGRID      = KGRID(1)
TZFIELD%NTYPE      = TYPECHAR
TZFIELD%NDIMS      = 0
TZFIELD%LTIMEDEP   = .FALSE.
CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,HTYPE)

IF (NVERB>=5) THEN
  WRITE(ILUOUTDIA,*)'  1st record (',TRIM(TZFIELD%CMNHNAME),'): OK'
ENDIF
!
! 2eme  enregistrement DIMENSIONS des MATRICES et LONGUEUR des TABLEAUX de CARACTERES et FLAGS de COMPRESSION sur les DIFFERENTS AXES
!
TZFIELD%CMNHNAME   = TRIM(HGROUP)//'.DIM'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = TRIM(HGROUP)//'.DIM'
TZFIELD%CUNITS     = ''
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = TRIM(YCOMMENT)
TZFIELD%NGRID      = KGRID(1)
TZFIELD%NTYPE      = TYPEINT
TZFIELD%NDIMS      = 1
TZFIELD%LTIMEDEP   = .FALSE.
SELECT CASE(HTYPE)
  CASE('CART','MASK','SPXY')
    ILENG = 34
    ALLOCATE(ITABCHAR(ILENG))
    ITABCHAR(1)=ILENTITRE; ITABCHAR(2)=ILENUNITE
    ITABCHAR(3)=ILENCOMMENT; ITABCHAR(4)=II
    ITABCHAR(5)=IJ; ITABCHAR(6)=IK
    ITABCHAR(7)=IT; ITABCHAR(8)=IN
    ITABCHAR(9)=IP; ITABCHAR(10)=KIL
    ITABCHAR(11)=KJL; ITABCHAR(12)=KKL
    ITABCHAR(13)=KIH; ITABCHAR(14)=KJH
    ITABCHAR(15)=KKH; ITABCHAR(16)=ICOMPX
    ITABCHAR(17)=ICOMPY; ITABCHAR(18)=ICOMPZ
    IF(HTYPE == 'MASK')THEN
!     ITABCHAR(10)=1; ITABCHAR(11)=1
!     ITABCHAR(13)=1; ITABCHAR(14)=1
      ITABCHAR(16)=1; ITABCHAR(17)=1
    ENDIF
    ITABCHAR(19)=INTRAJT; ITABCHAR(20)=IKTRAJX
    ITABCHAR(21)=IKTRAJY; ITABCHAR(22)=IKTRAJZ
    ITABCHAR(23)=ITTRAJX; ITABCHAR(24)=ITTRAJY
    ITABCHAR(25)=ITTRAJZ; ITABCHAR(26)=INTRAJX
    ITABCHAR(27)=INTRAJY; ITABCHAR(28)=INTRAJZ
    ITABCHAR(29)=IIMASK; ITABCHAR(30)=IJMASK
    ITABCHAR(31)=IKMASK; ITABCHAR(32)=ITMASK
    ITABCHAR(33)=INMASK; ITABCHAR(34)=IPMASK
    CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,ITABCHAR)
    DEALLOCATE(ITABCHAR)
    IF (NVERB>=5) THEN
      WRITE(ILUOUTDIA,*)' ILENTITRE,ILENUNITE,ILENCOMMENT ',ILENTITRE,ILENUNITE,ILENCOMMENT
    ENDIF
  CASE DEFAULT
    ILENG = 25 
    ALLOCATE(ITABCHAR(ILENG))
    ITABCHAR(1)=ILENTITRE; ITABCHAR(2)=ILENUNITE
    ITABCHAR(3)=ILENCOMMENT; ITABCHAR(4)=II
    ITABCHAR(5)=IJ; ITABCHAR(6)=IK
    ITABCHAR(7)=IT; ITABCHAR(8)=IN
    ITABCHAR(9)=IP
    ITABCHAR(10)=INTRAJT; ITABCHAR(11)=IKTRAJX
    ITABCHAR(12)=IKTRAJY; ITABCHAR(13)=IKTRAJZ
    ITABCHAR(14)=ITTRAJX; ITABCHAR(15)=ITTRAJY
    ITABCHAR(16)=ITTRAJZ; ITABCHAR(17)=INTRAJX
    ITABCHAR(18)=INTRAJY; ITABCHAR(19)=INTRAJZ
    ITABCHAR(20)=IIMASK; ITABCHAR(21)=IJMASK
    ITABCHAR(22)=IKMASK; ITABCHAR(23)=ITMASK
    ITABCHAR(24)=INMASK; ITABCHAR(25)=IPMASK
    CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,ITABCHAR)
    DEALLOCATE(ITABCHAR)
END SELECT
IF (NVERB>=5) THEN
  WRITE(ILUOUTDIA,*)'  2nd record (',TRIM(TZFIELD%CMNHNAME),'): OK'
ENDIF
!
! 3eme enregistrement TITRE
!
TZFIELD%CMNHNAME   = TRIM(HGROUP)//'.TITRE'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = TRIM(HGROUP)//'.TITRE'
TZFIELD%CUNITS     = ''
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = TRIM(YCOMMENT)
TZFIELD%NGRID      = KGRID(1)
TZFIELD%NTYPE      = TYPECHAR
TZFIELD%NDIMS      = 1
TZFIELD%LTIMEDEP   = .FALSE.
CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,HTITRE(1:IP))

IF (NVERB>=5) THEN
  WRITE(ILUOUTDIA,*)'  3rd record (',TRIM(TZFIELD%CMNHNAME),'): OK'
ENDIF
!
! 4eme enregistrement UNITE
!
TZFIELD%CMNHNAME   = TRIM(HGROUP)//'.UNITE'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = TRIM(HGROUP)//'.UNITE'
TZFIELD%CUNITS     = ''
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = TRIM(YCOMMENT)
TZFIELD%NGRID      = KGRID(1)
TZFIELD%NTYPE      = TYPECHAR
TZFIELD%NDIMS      = 1
TZFIELD%LTIMEDEP   = .FALSE.
CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,HUNITE(1:IP))

IF (NVERB>=5) THEN
  WRITE(ILUOUTDIA,*)'  4th record (',TRIM(TZFIELD%CMNHNAME),'): OK'
ENDIF
!
! 5eme enregistrement COMMENT
!
TZFIELD%CMNHNAME   = TRIM(HGROUP)//'.COMMENT'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = TRIM(HGROUP)//'.COMMENT'
TZFIELD%CUNITS     = ''
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = TRIM(YCOMMENT)
TZFIELD%NGRID      = KGRID(1)
TZFIELD%NTYPE      = TYPECHAR
TZFIELD%NDIMS      = 1
TZFIELD%LTIMEDEP   = .FALSE.
CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,HCOMMENT(1:IP))

IF (NVERB>=5) THEN
  WRITE(ILUOUTDIA,*)'  5th record (',TRIM(TZFIELD%CMNHNAME),'): OK'
ENDIF
!
! 6eme enregistrement PVAR
!
! Dans la mesure ou cette matrice risque d'etre tres volumineuse, on ecrira un 
! enregistrement par processus
!!!!!!!!!!!!!!!!  FUJI  compiler directive !!!!!!!!!!
!ocl scalar
!!!!!!!!!!!!!!!!  FUJI  compiler directive !!!!!!!!!!
DO J = 1,IP
  YJ = '   '
  IF(J < 10)WRITE(YJ,'(I1)')J ; YJ = ADJUSTL(YJ)
  IF(J >= 10 .AND. J < 100) THEN 
          WRITE(YJ,'(I2)')J ; YJ = ADJUSTL(YJ)
  ELSE IF(J >= 100 .AND. J < 1000) THEN 
          WRITE(YJ,'(I3)')J
  ENDIF
! BUG ...ca passe que si PRESENT(OICP) sinon OICP non defini 
IF (PRESENT(OICP) .AND. PRESENT(OJCP)) THEN
  IF(HTYPE == 'CART' .AND. .NOT. OICP .AND. .NOT. OJCP) THEN
    TZFIELD%CMNHNAME   = TRIM(HGROUP)//'.PROC'//YJ
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = TRIM(HUNITE(J))
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = TRIM(HTITRE(J))//' - '//TRIM(HCOMMENT(J))//' ('//TRIM(HUNITE(J))//')'
    TZFIELD%NGRID      = KGRID(J)
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 5
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD_BOX(TPDIAFILE,TZFIELD,'BUDGET',PVAR(:,:,:,:,:,J), &
                            KIL+JPHEXT,KIH+JPHEXT,KJL+JPHEXT,KJH+JPHEXT)
  ELSE
    TZFIELD%CMNHNAME   = TRIM(HGROUP)//'.PROC'//YJ
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = TRIM(HUNITE(J))
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = TRIM(HTITRE(J))//' - '//TRIM(HCOMMENT(J))//' ('//TRIM(HUNITE(J))//')'
    TZFIELD%NGRID      = KGRID(J)
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 5
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,PVAR(:,:,:,:,:,J))
  ENDIF
ELSE
    TZFIELD%CMNHNAME   = TRIM(HGROUP)//'.PROC'//YJ
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
  TZFIELD%CUNITS     = TRIM(HUNITE(J))
  TZFIELD%CDIR       = '--'
  TZFIELD%CCOMMENT   = TRIM(HTITRE(J))//' - '//TRIM(HCOMMENT(J))//' ('//TRIM(HUNITE(J))//')'
  TZFIELD%NGRID      = KGRID(J)
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 5
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,PVAR(:,:,:,:,:,J))
END IF
  IF (NVERB>=5) THEN
    WRITE(ILUOUTDIA,*)J,TRIM(TZFIELD%CMNHNAME)
  ENDIF
ENDDO
IF (NVERB>=5) THEN
  WRITE(ILUOUTDIA,*)'  6th record: OK'
ENDIF
!
! 7eme enregistrement TRAJT
!
TZFIELD%CMNHNAME   = TRIM(HGROUP)//'.TRAJT'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = TRIM(HGROUP)//'.TRAJT'
TZFIELD%CUNITS     = ''
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = TRIM(YCOMMENT)
TZFIELD%NGRID      = KGRID(1)
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .FALSE.
CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,PTRAJT)

IF (NVERB>=5) THEN
  WRITE(ILUOUTDIA,*)'  7th record (',TRIM(TZFIELD%CMNHNAME),'): OK'
ENDIF
!
! Dans certains cas
!
!
! 8eme enregistrement TRAJX
!
IF(PRESENT(PTRAJX))THEN
  TZFIELD%CMNHNAME   = TRIM(HGROUP)//'.TRAJX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = TRIM(HGROUP)//'.TRAJX'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = '--'
  TZFIELD%CCOMMENT   = TRIM(YCOMMENT)
  TZFIELD%NGRID      = KGRID(1)
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,PTRAJX)
ENDIF
!
!                        ou
!
IF(PRESENT(PMASK))THEN
  TZFIELD%CMNHNAME   = TRIM(HGROUP)//'.MASK'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = TRIM(HGROUP)//'.MASK'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = TRIM(YCOMMENT)
  TZFIELD%NGRID      = KGRID(1)
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 6
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,PMASK)
ENDIF
!
! 9eme enregistrement TRAJY
!
IF(PRESENT(PTRAJY))THEN
  TZFIELD%CMNHNAME   = TRIM(HGROUP)//'.TRAJY'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = TRIM(HGROUP)//'.TRAJY'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = '--'
  TZFIELD%CCOMMENT   = TRIM(YCOMMENT)
  TZFIELD%NGRID      = KGRID(1)
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,PTRAJY)
ENDIF
!
! 10eme enregistrement TRAJZ
!
IF(PRESENT(PTRAJZ))THEN
  TZFIELD%CMNHNAME   = TRIM(HGROUP)//'.TRAJZ'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = TRIM(HGROUP)//'.TRAJZ'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = '--'
  TZFIELD%CCOMMENT   = TRIM(YCOMMENT)
  TZFIELD%NGRID      = KGRID(1)
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,PTRAJZ)
ENDIF
!
! 11eme enregistrement PDATIME
!
TZFIELD%CMNHNAME   = TRIM(HGROUP)//'.DATIM'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = TRIM(HGROUP)//'.DATIM'
TZFIELD%CUNITS     = ''
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = TRIM(YCOMMENT)
TZFIELD%NGRID      = KGRID(1)
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .FALSE.
CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,PDATIME)
!
CALL MENU_DIACHRO(TPDIAFILE,HGROUP)
LPACK=GPACK
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
RETURN
END SUBROUTINE WRITE_DIACHRO
