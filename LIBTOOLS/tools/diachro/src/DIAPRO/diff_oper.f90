!     ######spl
      MODULE  MODI_DIFF_OPER
!     ##############################
!
INTERFACE
!
SUBROUTINE DIFF_OPER(K)
INTEGER :: K
END SUBROUTINE DIFF_OPER
!
END INTERFACE
!
END MODULE MODI_DIFF_OPER
!     #######################
      SUBROUTINE DIFF_OPER(K)
!     #######################
!
!!****  *DIFF_OPER* - 
!!
!!    PURPOSE
!!    -------
!      
!
!!**  METHOD
!!    ------
!!     
!!     N.A.
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ALLOC_FORDIACHRO
USE MODD_ALLOC2_FORDIACHRO
USE MODD_RESOLVCAR
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODD_FILES_DIACHRO
USE MODD_TIT
USE MODD_TYPE_AND_LH
USE MODD_MEMCV
USE MODN_NCAR

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

INTEGER  :: K
!
!*       0.1   Local variables
!              ---------------
INTEGER  :: JLOOPT, ILENT, ILENT2, ILENC, ILENFI, ILENTIM, ILENTIT
INTEGER  :: JA, JME, JME2, IP
INTEGER  :: INDK, INDKM1, INDNN
INTEGER  :: IFAC, INUMPM
INTEGER,SAVE  :: IDIFK, INIV1, INIV2
CHARACTER(LEN=800),SAVE :: YCAR
CHARACTER(LEN=100),SAVE :: YTITB3
CHARACTER(LEN=8) :: YTIM
REAL     :: ZSPVAL
!
!------------------------------------------------------------------------------
IDIFK=1
INIV1=0; INIV2=0
IF(LEV .OR. LPR .OR. LTK)THEN
  print *,' NON NON NON NON NON NON NON NON NON NON NON NON NON NON NON NON NON'
  print *,' **Diff_oper Operation IMPOSSIBLE . Seules sont autorisees les differences '
  print *,' entre 2 niveaux de modele (identiques ou non) ou entre 2 altitudes semblables'
  LPBREAD=.TRUE.
  RETURN
ENDIF
if(nverbia > 0)then
  print *,' ** Diff_oper NUMPM(MAX(K-1,1)),K-1,K ',NUMPM(MAX(K-1,1)),NUMPM(K-1),NUMPM(K),' K argument ',K
endif
IF(LFT .OR. LPVKT .OR. LFT1 .OR. LPVKT1 .AND. LSPVALT)THEN
  ZSPVAL=XSPVAL
  XSPVAL=XSPVALT
ENDIF
! Initialisation unquement si il n'y a pas eu de + ou - avant
INUMPM=0
DO JA=1,K-1
  IF(NUMPM(JA) == 1 .OR. NUMPM(JA) == 2)THEN
    INUMPM=INUMPM+1
  ENDIF
ENDDO
IF(INUMPM == 0)THEN
  YCAR(1:LEN(YCAR))=' '
  ILENC=0
  ILENC=ILENC+1
ELSE
  ILENC=LEN_TRIM(YCAR)
  ILENC=ILENC+1
  YCAR(ILENC:ILENC+2)=' , '
  ILENC=ILENC+3
ENDIF
if(nverbia > 0)then
  print *,' **Diff_oper ILENC entree ',ILENC
endif
YTIM(1:LEN(YTIM))=' '
IF(NUMFILECUR2 /= NUMFILECUR)THEN
  DO JA=1,NBFILES
    IF(NUMFILES(JA) == NUMFILECUR)THEN
      JME=JA
    ENDIF
    IF(NUMFILES(JA) == NUMFILECUR2)THEN
      JME2=JA
    ENDIF
  ENDDO
ENDIF
!
! Traitement d'un seul temps
!
IF(NBNDIA(K) /= 1 .OR. NBNDIA(K-1) /= 1 )THEN       
  LPBREAD=.TRUE.
  print *,' NB DE MASQUES (ou STATIONS) DEMANDES > 1 . PAS DE TRACE '
  IF(LFT .OR. LPVKT .OR. LFT1 .OR. LPVKT1 .AND. LSPVALT)THEN
    XSPVAL=ZSPVAL
  ENDIF
  RETURN
ENDIF     

INDK=NNDIA(1,K); INDKM1=NNDIA(1,K-1)

!Je mets directement 1 pour le 1er indice de xlvlzdia puisque= nblvlzdia(k,indk)
          if(nblvlzdia(k) == 1 .and. nblvlzdia(k-1) == 1 .AND. &
	     xlvlzdia(1,k) /= xlvlzdia(1,k-1) )then
            print *,' NON NON NON NON NON NON NON NON NON NON NON NON NON NON NON NON NON'
            print *,' **Diff_oper Operation IMPOSSIBLE . Seules sont autorisees les differences '
            print *,' entre 2 niveaux de modele (identiques ou non) ou entre 2 altitudes semblables'
	    print *,' **Altitudes demandees : xlvlzdia(1,k-1),xlvlzdia(1,k) ',xlvlzdia(1,k-1),xlvlzdia(1,k)
            LPBREAD=.TRUE.
            IF(LFT .OR. LPVKT .OR. LFT1 .OR. LPVKT1 .AND. LSPVALT)THEN
              XSPVAL=ZSPVAL
            ENDIF
            RETURN
          endif

IP=NPROCDIA(NBPROCDIA(K),K)

IF(NGRIDIA(IP) /= NGRIDIAM)THEN
  IF(CTYPE == 'MASK')THEN
    SELECT CASE(NGRIDIAM)
      CASE(1,2,3,5)
         SELECT CASE(NGRIDIA(IP))
           CASE(1,2,3,5)
             print *,' *** diff_oper Type MASK NGRIDIAM, NGRIDIA(IP), pas d interpolation ',NGRIDIAM,NGRIDIA(IP)
           CASE(4,6,7)
             print *,' *** diff_oper Type MASK NGRIDIAM, NGRIDIA(IP),&
& interpolation en K sur la grille du 1er processus traite (NGRIDIAM) ',NGRIDIAM,NGRIDIA(IP)
             CALL INTERP_GRIDS(K)
         END SELECT
      CASE(4,6,7)
         SELECT CASE(NGRIDIA(IP))
           CASE(1,2,3,5)
             print *,' *** diff_oper Type MASK NGRIDIAM, NGRIDIA(IP), interpolation en K&
& sur la grille du 1er processus traite (NGRIDIAM) ',NGRIDIAM,NGRIDIA(IP)
             CALL INTERP_GRIDS(K)
           CASE(4,6,7)
             print *,' *** diff_oper Type MASK NGRIDIAM, NGRIDIA(IP), pas d interpolation ',NGRIDIAM,NGRIDIA(IP)
         END SELECT
    END SELECT
  ELSE
    print *,' *** diff_oper NGRIDIAM, NGRIDIA(IP), interpolation ',NGRIDIAM,NGRIDIA(IP)
    CALL INTERP_GRIDS(K)
  ENDIF
ENDIF

if(nverbia >0)then
print *,' DIFF_OPER  INDK INDKM1 ',INDK,INDKM1
print *,' DIFF_OPER  NBTIMEDIA(K,INDK) NBTIMEDIA(K-1,INDKM1) ', &
NBTIMEDIA(K,INDK),NBTIMEDIA(K-1,INDKM1)
endif

!******************************** 1 seul temps

IF(NBTIMEDIA(K,INDK) == 1 .AND. NBTIMEDIA(K-1,INDKM1) == 1)THEN
! print *,' DIFF_OPER XVAR ',XVAR
! print *,' DIFF_OPER XVAR2 ',XVAR2

  IF(NBPROCDIA(K) /= 1 .OR. NBPROCDIA(K-1) /= 1 )THEN        !++++++++++
    LPBREAD=.TRUE.
    print *,' NB DE PROCESSUS DEMANDES > 1 . PAS DE TRACE '
    IF(LFT .OR. LPVKT .OR. LFT1 .OR. LPVKT1 .AND. LSPVALT)THEN
      XSPVAL=ZSPVAL
    ENDIF
    RETURN

  ELSE                                                       !++++++++++
!-----------------------------------------------------------------------
!   IF(K <= 2)THEN
    IF(NUMPM(K-1) == 0 .OR. NUMPM(MAX(K-1,1)) == 3)THEN
      if(nverbia > 0)then
        print *,' diff_oper1 K NUMPM(K-1) ',K,NUMPM(K-1)
      endif
      if(nblvlkdia(k,indk) == 1 .and. nblvlkdia(k-1,indk) == 1)then
        print *,' diff_oper1 Niveaux en K de part et d autre de MINUS '
        print *,' K1= ',NLVLKDIA(1,K-1,INDK),' K2= ',NLVLKDIA(1,K,INDK)
! INIV1=le 1er dans la directive <-> K-1 , INIV2=le 2e=courant <-> K
! Janv 2001
	IF(CTYPE == 'CART' .OR. CTYPE == 'MASK')THEN
            INIV1=NLVLKDIA(1,K-1,INDK)-NKL+1
            INIV2=NLVLKDIA(1,K,INDK)-NKL+1
	ELSE
! Janv 2001
          INIV1=NLVLKDIA(1,K-1,INDK)
          INIV2=NLVLKDIA(1,K,INDK)
! Janv 2001
	ENDIF
! Janv 2001
        IDIFK=2
      endif

      IF(LMUMVM .OR. LMUTVT .OR. LUMVM .OR. LUTVT .OR. &
         LSUMVM .OR. LSUTVT .OR. LMLSUMVM .OR. LMLSUTVT)THEN
        CGROUPS(K-1)=ADJUSTL(CGROUPS(K-1))
	ILENTIT=LEN_TRIM(CGROUPS(K-1))
        YCAR(ILENC:ILENC+ILENTIT-1)=CGROUPS(K-1)(1:ILENTIT)
      ELSE
        CTITRE2(NPROCDIA(1,K-1))=ADJUSTL(CTITRE2(NPROCDIA(1,K-1)))
	ILENTIT=LEN_TRIM(CTITRE2(NPROCDIA(1,K-1)))
        YCAR(ILENC:ILENC+ILENTIT-1)=CTITRE2(NPROCDIA(1,K-1))(1:ILENTIT)
      ENDIF
      YCAR=ADJUSTL(YCAR)
      ILENC=LEN_TRIM(YCAR)
      ILENC=ILENC+2
      YCAR(ILENC:ILENC)='('
      ILENC=ILENC+1
      IF(NUMFILECUR2 /= NUMFILECUR)THEN
        CFILEDIAS(JME2)=ADJUSTL(CFILEDIAS(JME2))
        ILENFI=LEN_TRIM(CFILEDIAS(JME2))
        YCAR(ILENC:ILENC+ILENFI-1)=CFILEDIAS(JME2)(1:ILENFI)
        ILENC=ILENC+ILENFI
        YCAR(ILENC:ILENC+1)=')('
        ILENC=ILENC+2
      ENDIF

      IF(IDIFK == 2 .AND. INIV1 /= INIV2)THEN
        YCAR(ILENC:ILENC+1)='K='
        ILENC=ILENC+2
        WRITE(YCAR(ILENC:ILENC+1),'(I2)')INIV1
        ILENC=ILENC+2
        YCAR(ILENC:ILENC+1)=')('
        ILENC=ILENC+2
      ENDIF

      IF(CTYPE == 'DRST' .OR. CTYPE == 'RSPL' .OR. CTYPE == 'RAPL')THEN
        WRITE(YTIM,'(F8.0)')XTRAJT2(NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),INDKM1)
      ELSE
        WRITE(YTIM,'(F8.0)')XTRAJT2(NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),1)
      ENDIF
      YTIM=ADJUSTL(YTIM)
!     print *,' YTIM ',YTIM
      ILENTIM=LEN_TRIM(YTIM)
      YCAR(ILENC:ILENC+ILENTIM-1)=YTIM(1:ILENTIM)
      YTIM(1:LEN(YTIM))=' '
      ILENC=ILENC+ILENTIM
      IF(LMINUS)THEN
        YCAR(ILENC:ILENC+3)=') - '
      ELSE IF(LPLUS)THEN
        YCAR(ILENC:ILENC+3)=') + '
      ENDIF
      ILENC=ILENC+4

    ELSE

      YCAR=ADJUSTL(CTITB3)
      ILENC=LEN_TRIM(YCAR)+1
      IF(LMINUS)THEN
        YCAR(ILENC:ILENC+2)=' - '
      ELSE IF(LPLUS)THEN
        YCAR(ILENC:ILENC+2)=' + '
      ENDIF
      ILENC=ILENC+3

    ENDIF
!-----------------------------------------------------------------------
    IF(LMUMVM .OR. LMUTVT .OR. LUMVM .OR. LUTVT .OR. &
       LSUMVM .OR. LSUTVT .OR. LMLSUMVM .OR. LMLSUTVT)THEN
      ILENTIT=LEN_TRIM(CGROUPS(K))
      YCAR(ILENC:ILENC+ILENTIT-1)=CGROUPS(K)(1:ILENTIT)
    ELSE
      CTITRE(NPROCDIA(1,K))=ADJUSTL(CTITRE(NPROCDIA(1,K)))
      ILENTIT=LEN_TRIM(CTITRE(NPROCDIA(1,K)))
      YCAR(ILENC:ILENC+ILENTIT-1)=CTITRE(NPROCDIA(1,K))(1:ILENTIT)
    ENDIF
    YCAR=ADJUSTL(YCAR)
    ILENC=LEN_TRIM(YCAR)
    ILENC=ILENC+2
    YCAR(ILENC:ILENC)='('
    ILENC=ILENC+1
!   print *,' AV 2eme IF'
    IF(NUMFILECUR2 /= NUMFILECUR)THEN
      CFILEDIAS(JME)=ADJUSTL(CFILEDIAS(JME))
      ILENFI=LEN_TRIM(CFILEDIAS(JME))
      YCAR(ILENC:ILENC+ILENFI-1)=CFILEDIAS(JME)(1:ILENFI)
      ILENC=ILENC+ILENFI
      YCAR(ILENC:ILENC+1)=')('
      ILENC=ILENC+2
    ENDIF

      IF(IDIFK == 2 .AND. INIV1 /= INIV2)THEN
        YCAR(ILENC:ILENC+1)='K='
        ILENC=ILENC+2
        WRITE(YCAR(ILENC:ILENC+1),'(I2)')INIV2
        ILENC=ILENC+2
        YCAR(ILENC:ILENC+1)=')('
        ILENC=ILENC+2
      ENDIF

      IF(CTYPE == 'DRST' .OR. CTYPE == 'RSPL' .OR. CTYPE == 'RAPL')THEN
        WRITE(YTIM,'(F8.0)')XTRAJT(NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),INDK)
      ELSE
        WRITE(YTIM,'(F8.0)')XTRAJT(NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1)
      ENDIF
    YTIM=ADJUSTL(YTIM)
!   print *,' YTIM ',YTIM,' ILENC ',ILENC
    ILENTIM=LEN_TRIM(YTIM)
    YCAR(ILENC:ILENC+ILENTIM-1)=YTIM(1:ILENTIM)
    print *,' YCAR ',YCAR(1:LEN_TRIM(YCAR))
!   YTIM(1:LEN(YTIM))=' '
    ILENC=ILENC+ILENTIM
!   print *,' ILENC ',ILENC
    YCAR(ILENC:ILENC)=')'
!   print *,' YCAR ',YCAR

!   IF(K <= 2)THEN
    IF(NUMPM(K-1) == 0 .OR. NUMPM(MAX(K-1,1)) == 3)THEN
      if(nverbia > 0)then
        print *,' diff_oper1 K NUMPM(K-1) ',K,NUMPM(K-1)
      endif
      if(nblvlkdia(k,indk) == 1 .and. nblvlkdia(k-1,indk) == 1)then
        print *,' diff_oper1-2 Niveaux en K de part et d autre de MINUS '
        print *,' K1= ',NLVLKDIA(1,K-1,INDK),' K2= ',NLVLKDIA(1,K,INDK)
! Janv 2001
	IF(CTYPE == 'CART' .OR. CTYPE == 'MASK')THEN
          INIV1=NLVLKDIA(1,K-1,INDK)-NKL+1
          INIV2=NLVLKDIA(1,K,INDK)-NKL+1
	ELSE
! Janv 2001
          INIV1=NLVLKDIA(1,K-1,INDK)
          INIV2=NLVLKDIA(1,K,INDK)
! Janv 2001
        ENDIF
! Janv 2001
        IDIFK=2
      endif
      LTITDEF=.FALSE.
!     YTITB3(1:LEN(YTITB3))=' '
!     YTITB3=ADJUSTL(CTITB3)
    ENDIF
!!! 1/3/04
    IF(CTITB3 /= 'DEFAULT')THEN
    ELSE
!!! 1/3/04
    CTITB3=ADJUSTL(YCAR(1:100))
    CTITB3=ADJUSTL(CTITB3)
!!! 1/3/04
    ENDIF
!!! 1/3/04
    print *,' CTITB3 ',CTITB3
    IF(LMINUS)THEN
      IFAC=-1
    ELSE IF(LPLUS)THEN
      IFAC=1
    ENDIF

    IF(IDIFK == 2)THEN
!!Mai 2003
      WHERE((XVAR(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),INDK, &
!     WHERE((XVAR(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
        NPROCDIA(NBPROCDIA(K),K)) == XSPVAL) .OR.  &
        (XVAR2(:,:,INIV1,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),INDKM1, &
!       (XVAR2(:,:,INIV1,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),1, &
        NPROCDIA(NBPROCDIA(K-1),K-1)) == XSPVAL))
        XVAR(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),INDK, &
!       XVAR(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
        NPROCDIA(NBPROCDIA(K),K))= XSPVAL
      ELSEWHERE
        XVAR(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),INDK,  &
!       XVAR(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1,  &
        NPROCDIA(NBPROCDIA(K),K))=  &
        XVAR2(:,:,INIV1,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),INDKM1, &
!       XVAR2(:,:,INIV1,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),1, &
        NPROCDIA(NBPROCDIA(K-1),K-1))+  &
        IFAC * XVAR(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),INDK, &
!       IFAC * XVAR(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
        NPROCDIA(NBPROCDIA(K),K))
      END WHERE
    ELSE
      WHERE((XVAR(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),INDK, &
!     WHERE((XVAR(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
        NPROCDIA(NBPROCDIA(K),K)) == XSPVAL) .OR.  &
        (XVAR2(:,:,:,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),INDKM1, &
!       (XVAR2(:,:,:,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),1, &
        NPROCDIA(NBPROCDIA(K-1),K-1)) == XSPVAL))
        XVAR(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),INDK, &
!       XVAR(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
        NPROCDIA(NBPROCDIA(K),K))= XSPVAL
      ELSEWHERE
        XVAR(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),INDK,  &
!       XVAR(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1,  &
        NPROCDIA(NBPROCDIA(K),K))=  &
        XVAR2(:,:,:,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),INDKM1, &
!       XVAR2(:,:,:,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),1, &
        NPROCDIA(NBPROCDIA(K-1),K-1))+  &
        IFAC * XVAR(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),INDK, &
!       IFAC * XVAR(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
!!Mai 2003
        NPROCDIA(NBPROCDIA(K),K))
      END WHERE
    ENDIF

    IF(ALLOCATED(XUMEM))THEN
      IF(IDIFK == 2)THEN
        WHERE((XU(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
          NPROCDIA(NBPROCDIA(K),K)) == XSPVAL) .OR.  &
          (XUMEM(:,:,INIV1,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),1, &
          NPROCDIA(NBPROCDIA(K-1),K-1)) == XSPVAL))
	  XU(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
	  NPROCDIA(NBPROCDIA(K),K))= XSPVAL
        ELSEWHERE
          XU(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1,  &
	  NPROCDIA(NBPROCDIA(K),K))=  &
          XUMEM(:,:,INIV1,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),1, &
	  NPROCDIA(NBPROCDIA(K-1),K-1))+  &
          IFAC * XU(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
	  NPROCDIA(NBPROCDIA(K),K))
        END WHERE
      ELSE
        WHERE((XU(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
          NPROCDIA(NBPROCDIA(K),K)) == XSPVAL) .OR.  &
          (XUMEM(:,:,:,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),1, &
          NPROCDIA(NBPROCDIA(K-1),K-1)) == XSPVAL))
	  XU(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
	  NPROCDIA(NBPROCDIA(K),K))= XSPVAL
        ELSEWHERE
          XU(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1,  &
	  NPROCDIA(NBPROCDIA(K),K))=  &
          XUMEM(:,:,:,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),1, &
	  NPROCDIA(NBPROCDIA(K-1),K-1))+  &
          IFAC * XU(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
	  NPROCDIA(NBPROCDIA(K),K))
        END WHERE
      ENDIF
    ENDIF
    IF(ALLOCATED(XVMEM))THEN
      IF(IDIFK == 2)THEN
        WHERE((XV(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
          NPROCDIA(NBPROCDIA(K),K)) == XSPVAL) .OR.  &
          (XVMEM(:,:,INIV1,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),1, &
          NPROCDIA(NBPROCDIA(K-1),K-1)) == XSPVAL))
	  XV(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
	  NPROCDIA(NBPROCDIA(K),K))= XSPVAL
        ELSEWHERE
          XV(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1,  &
	  NPROCDIA(NBPROCDIA(K),K))=  &
          XVMEM(:,:,INIV1,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),1, &
	  NPROCDIA(NBPROCDIA(K-1),K-1))+  &
          IFAC * XV(:,:,INIV2,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
	  NPROCDIA(NBPROCDIA(K),K))
        END WHERE
      ELSE
        WHERE((XV(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
          NPROCDIA(NBPROCDIA(K),K)) == XSPVAL) .OR.  &
          (XVMEM(:,:,:,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),1, &
          NPROCDIA(NBPROCDIA(K-1),K-1)) == XSPVAL))
	  XV(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
	  NPROCDIA(NBPROCDIA(K),K))= XSPVAL
        ELSEWHERE
          XV(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1,  &
	  NPROCDIA(NBPROCDIA(K),K))=  &
          XVMEM(:,:,:,NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),1, &
	  NPROCDIA(NBPROCDIA(K-1),K-1))+  &
          IFAC * XV(:,:,:,NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1, &
	  NPROCDIA(NBPROCDIA(K),K))
        END WHERE
      ENDIF
    ENDIF
  ENDIF                                                      !++++++++++

!********************************   plusieurs temps

ELSE

! Expression du temps en sequentiel   SSSSSSSSSSSS
  IF(.NOT.LTINCRDIA(K,INDK))THEN

    IF(NBTIMEDIA(K,INDK) == NBTIMEDIA(K-1,INDKM1))THEN
! Intervalle de temps de meme longueur pour les 2 champs OK

      IF(NBPROCDIA(K) /= 1 .OR. NBPROCDIA(K-1) /= 1 )THEN
	LPBREAD=.TRUE.
	print *,' NB DE PROCESSUS DEMANDES > 1 . PAS DE TRACE '
        IF(LFT .OR. LPVKT .OR. LFT1 .OR. LPVKT1 .AND. LSPVALT)THEN
          XSPVAL=ZSPVAL
        ENDIF
        RETURN
      ELSE
!       print *,' PASSAGE ICI'
        IF(LMUMVM .OR. LMUTVT .OR. LUMVM .OR. LUTVT .OR. &
           LSUMVM .OR. LSUTVT .OR. LMLSUMVM .OR. LMLSUTVT)THEN
          CGROUPS(K-1)=ADJUSTL(CGROUPS(K-1))
	  ILENTIT=LEN_TRIM(CGROUPS(K-1))
          YCAR(ILENC:ILENC+ILENTIT-1)=CGROUPS(K-1)(1:ILENTIT)
        ELSE
          CTITRE2(NPROCDIA(1,K-1))=ADJUSTL(CTITRE2(NPROCDIA(1,K-1)))
	  ILENTIT=LEN_TRIM(CTITRE2(NPROCDIA(1,K-1)))
          YCAR(ILENC:ILENC+ILENTIT-1)=CTITRE2(NPROCDIA(1,K-1))(1:ILENTIT)
        ENDIF
        YCAR=ADJUSTL(YCAR)
        ILENC=LEN_TRIM(YCAR)
        ILENC=ILENC+2
        YCAR(ILENC:ILENC)='('
        ILENC=ILENC+1
        IF(NUMFILECUR2 /= NUMFILECUR)THEN
	  CFILEDIAS(JME2)=ADJUSTL(CFILEDIAS(JME2))
	  ILENFI=LEN_TRIM(CFILEDIAS(JME2))
          YCAR(ILENC:ILENC+ILENFI-1)=CFILEDIAS(JME2)(1:ILENFI)
          ILENC=ILENC+ILENFI
          YCAR(ILENC:ILENC+1)=')('
          ILENC=ILENC+2
        ENDIF

      IF(IDIFK == 2 .AND. INIV1 /= INIV2)THEN
        YCAR(ILENC:ILENC+1)='K='
        ILENC=ILENC+2
        WRITE(YCAR(ILENC:ILENC+1),'(I2)')INIV1
        ILENC=ILENC+2
        YCAR(ILENC:ILENC+1)=')('
        ILENC=ILENC+2
      ENDIF

! Ecriture de la premiere serie de temps
!-----------------------------------------------------------------------
!       IF(K <=2)THEN
        IF(NUMPM(K-1) == 0 .OR. NUMPM(MAX(K-1,1)) == 3)THEN
          if(nverbia > 0)then
            print *,' diff_oper2 K NUMPM(K-1) ',K,NUMPM(K-1)
          endif
          if(nblvlkdia(k,indk) == 1 .and. nblvlkdia(k-1,indk) == 1)then
            print *,' diff_oper2 Niveaux en K de part et d autre de MINUS '
            print *,' K1= ',NLVLKDIA(1,K-1,INDK),' K2= ',NLVLKDIA(1,K,INDK)
! Janv 2001
	    IF(CTYPE == 'CART' .OR. CTYPE == 'MASK')THEN
              INIV1=NLVLKDIA(1,K-1,INDK)-NKL+1
              INIV2=NLVLKDIA(1,K,INDK)-NKL+1
	    ELSE
! Janv 2001
              INIV1=NLVLKDIA(1,K-1,INDK)
              INIV2=NLVLKDIA(1,K,INDK)
! Janv 2001
            ENDIF
! Janv 2001
            IDIFK=2
          endif

          IF(CTYPE == 'DRST' .OR. CTYPE == 'RSPL' .OR. CTYPE == 'RAPL')THEN
            WRITE(YTIM,'(F8.0)')XTRAJT2(NTIMEDIA(1,K-1,INDKM1),INDKM1)
          ELSE
            WRITE(YTIM,'(F8.0)')XTRAJT2(NTIMEDIA(1,K-1,INDKM1),1)
          ENDIF
          YTIM=ADJUSTL(YTIM)
          ILENTIM=LEN_TRIM(YTIM)
          YCAR(ILENC:ILENC+ILENTIM-1)=YTIM(1:ILENTIM)
          YTIM(1:LEN(YTIM))=' '
          ILENC=ILENC+ILENTIM
      IF(CTYPE == 'DRST' .OR. CTYPE == 'RSPL' .OR. CTYPE == 'RAPL')THEN
        INDNN=INDKM1
      ELSE
        INDNN=1
      ENDIF
!         IF(XTRAJT2(NTIMEDIA(1,K-1,INDKM1),1) /= XTRAJT2(NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),1))THEN
          IF(XTRAJT2(NTIMEDIA(1,K-1,INDKM1),INDNN) /= XTRAJT2(NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),INDNN))THEN
            YCAR(ILENC:ILENC+2)=' - '
	    ILENC=ILENC+3
              WRITE(YTIM,'(F8.0)')XTRAJT2(NTIMEDIA(NBTIMEDIA(K-1,INDKM1),K-1,INDKM1),INDNN)
            YTIM=ADJUSTL(YTIM)
            ILENTIM=LEN_TRIM(YTIM)
            YCAR(ILENC:ILENC+ILENTIM-1)=YTIM(1:ILENTIM)
            YTIM(1:LEN(YTIM))=' '
            ILENC=ILENC+ILENTIM

	  ELSE

	    ILENC=ILENC+1
	  ENDIF
          IF(LMINUS)THEN
	    YCAR(ILENC:ILENC+3)=') - '
          ELSE IF(LPLUS)THEN
	    YCAR(ILENC:ILENC+3)=') + '
          ENDIF
	  ILENC=ILENC+4

        ELSE
    
          YCAR=ADJUSTL(CTITB3)
          ILENC=LEN_TRIM(YCAR)+1
          IF(LMINUS)THEN
            YCAR(ILENC:ILENC+2)=' - '
          ELSE IF(LPLUS)THEN
            YCAR(ILENC:ILENC+2)=' + '
          ENDIF
          ILENC=ILENC+3
    
        ENDIF
!-----------------------------------------------------------------------
!  FIN Ecriture de la premiere serie de temps
        IF(LMUMVM .OR. LMUTVT .OR. LUMVM .OR. LUTVT .OR. &
           LSUMVM .OR. LSUTVT .OR. LMLSUMVM .OR. LMLSUTVT)THEN
          ILENTIT=LEN_TRIM(CGROUPS(K))
          YCAR(ILENC:ILENC+ILENTIT-1)=CGROUPS(K)(1:ILENTIT)
        ELSE
          CTITRE(NPROCDIA(1,K))=ADJUSTL(CTITRE(NPROCDIA(1,K)))
          ILENTIT=LEN_TRIM(CTITRE(NPROCDIA(1,K)))
          YCAR(ILENC:ILENC+ILENTIT-1)=CTITRE(NPROCDIA(1,K))(1:ILENTIT)
        ENDIF
        YCAR=ADJUSTL(YCAR)
        ILENC=LEN_TRIM(YCAR)
        ILENC=ILENC+2
        YCAR(ILENC:ILENC)='('
        ILENC=ILENC+1
        IF(NUMFILECUR2 /= NUMFILECUR)THEN
	  CFILEDIAS(JME)=ADJUSTL(CFILEDIAS(JME))
	  ILENFI=LEN_TRIM(CFILEDIAS(JME))
          YCAR(ILENC:ILENC+ILENFI-1)=CFILEDIAS(JME)(1:ILENFI)
          ILENC=ILENC+ILENFI
          YCAR(ILENC:ILENC+1)=')('
          ILENC=ILENC+2
        ENDIF

      IF(IDIFK == 2 .AND. INIV1 /= INIV2)THEN
        YCAR(ILENC:ILENC+1)='K='
        ILENC=ILENC+2
        WRITE(YCAR(ILENC:ILENC+1),'(I2)')INIV2
        ILENC=ILENC+2
        YCAR(ILENC:ILENC+1)=')('
        ILENC=ILENC+2
      ENDIF

! Ecriture de la deuxieme serie de temps
      IF(CTYPE == 'DRST' .OR. CTYPE == 'RSPL' .OR. CTYPE == 'RAPL')THEN
        WRITE(YTIM,'(F8.0)')XTRAJT(NTIMEDIA(1,K,INDK),INDK)
      ELSE
        WRITE(YTIM,'(F8.0)')XTRAJT(NTIMEDIA(1,K,INDK),1)
      ENDIF
        YTIM=ADJUSTL(YTIM)
        ILENTIM=LEN_TRIM(YTIM)
        YCAR(ILENC:ILENC+ILENTIM-1)=YTIM(1:ILENTIM)
        YTIM(1:LEN(YTIM))=' '
        ILENC=ILENC+ILENTIM

        IF(CTYPE == 'DRST' .OR. CTYPE == 'RSPL' .OR. CTYPE == 'RAPL')THEN
          INDNN=INDK
        ELSE
          INDNN=1
        ENDIF
!       IF(XTRAJT(NTIMEDIA(1,K,INDK),1) /= XTRAJT(NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1))THEN
        IF(XTRAJT(NTIMEDIA(1,K,INDK),INDNN) /= XTRAJT(NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),INDNN))THEN
	  YCAR(ILENC:ILENC+2)=' - '
	  ILENC=ILENC+3
          WRITE(YTIM,'(F8.0)')XTRAJT(NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),INDNN)
!         WRITE(YTIM,'(F8.0)')XTRAJT(NTIMEDIA(NBTIMEDIA(K,INDK),K,INDK),1)
          YTIM=ADJUSTL(YTIM)
          ILENTIM=LEN_TRIM(YTIM)
          YCAR(ILENC:ILENC+ILENTIM-1)=YTIM(1:ILENTIM)
          YTIM(1:LEN(YTIM))=' '
          ILENC=ILENC+ILENTIM

	ELSE

	  ILENC=ILENC+1
	ENDIF
	YCAR(ILENC:ILENC)=')'

!       IF(K <= 2)THEN
        IF(NUMPM(K-1) == 0 .OR. NUMPM(MAX(K-1,1)) == 3)THEN
          if(nverbia > 0)then
            print *,' diff_oper2 K NUMPM(K-1) ',K,NUMPM(K-1)
          endif
          if(nblvlkdia(k,indk) == 1 .and. nblvlkdia(k-1,indk) == 1)then
            print *,' diff_oper2-2 Niveaux en K de part et d autre de MINUS '
            print *,' K1= ',NLVLKDIA(1,K-1,INDK),' K2= ',NLVLKDIA(1,K,INDK)
! Janv 2001
	    IF(CTYPE == 'CART' .OR. CTYPE == 'MASK')THEN
              INIV1=NLVLKDIA(1,K-1,INDK)-NKL+1
              INIV2=NLVLKDIA(1,K,INDK)-NKL+1
	    ELSE
! Janv 2001
              INIV1=NLVLKDIA(1,K-1,INDK)
              INIV2=NLVLKDIA(1,K,INDK)
! Janv 2001
            ENDIF
! Janv 2001
            IDIFK=2
          endif
          LTITDEF=.FALSE.
!         YTITB3(1:LEN(YTITB3))=' '
!         YTITB3=ADJUSTL(CTITB3)
        ENDIF

!!! 1/3/04
    IF(CTITB3 /= 'DEFAULT')THEN
    ELSE
!!! 1/3/04
        CTITB3=ADJUSTL(YCAR(1:100))
        CTITB3=ADJUSTL(CTITB3)
!!! 1/3/04
    ENDIF
!!! 1/3/04
        print *,' CTITB3 ',CTITB3
        DO JLOOPT=1,NBTIMEDIA(K,INDK)
          IF(LMINUS)THEN
            IFAC=-1
          ELSE IF(LPLUS)THEN
            IFAC=1
          ENDIF
          IF(IDIFK == 2)THEN
! Mai 2003
            WHERE((XVAR(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),INDK,NPROCDIA(1,K))==XSPVAL)&
!           WHERE((XVAR(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))==XSPVAL)&
	    .OR. (XVAR2(:,:,INIV1,NTIMEDIA(JLOOPT,K-1,INDKM1),INDKM1, &
!           .OR. (XVAR2(:,:,INIV1,NTIMEDIA(JLOOPT,K-1,INDKM1),1, &
	    NPROCDIA(1,K-1))==XSPVAL))
	      XVAR(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),INDK,NPROCDIA(1,K))=XSPVAL
!             XVAR(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))=XSPVAL
            ELSEWHERE
  	      XVAR(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),INDK,NPROCDIA(1,K)) =  &
! 	      XVAR(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K)) =  &
  	      XVAR2(:,:,INIV1,NTIMEDIA(JLOOPT,K-1,INDKM1),INDKM1,NPROCDIA(1,K-1)) + &
! 	      XVAR2(:,:,INIV1,NTIMEDIA(JLOOPT,K-1,INDKM1),1,NPROCDIA(1,K-1)) + &
  	      IFAC * XVAR(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),INDK,NPROCDIA(1,K))
! 	      IFAC * XVAR(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))
            END WHERE
	  ELSE
            WHERE((XVAR(:,:,:,NTIMEDIA(JLOOPT,K,INDK),INDK,NPROCDIA(1,K))==XSPVAL)&
!           WHERE((XVAR(:,:,:,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))==XSPVAL)&
	    .OR. (XVAR2(:,:,:,NTIMEDIA(JLOOPT,K-1,INDKM1),INDKM1, &
!           .OR. (XVAR2(:,:,:,NTIMEDIA(JLOOPT,K-1,INDKM1),1, &
	    NPROCDIA(1,K-1))==XSPVAL))
	      XVAR(:,:,:,NTIMEDIA(JLOOPT,K,INDK),INDK,NPROCDIA(1,K))=XSPVAL
!             XVAR(:,:,:,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))=XSPVAL
            ELSEWHERE
  	      XVAR(:,:,:,NTIMEDIA(JLOOPT,K,INDK),INDK,NPROCDIA(1,K)) =  &
! 	      XVAR(:,:,:,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K)) =  &
  	      XVAR2(:,:,:,NTIMEDIA(JLOOPT,K-1,INDKM1),INDKM1,NPROCDIA(1,K-1)) + &
! 	      XVAR2(:,:,:,NTIMEDIA(JLOOPT,K-1,INDKM1),1,NPROCDIA(1,K-1)) + &
  	      IFAC * XVAR(:,:,:,NTIMEDIA(JLOOPT,K,INDK),INDK,NPROCDIA(1,K))
! 	      IFAC * XVAR(:,:,:,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))
            END WHERE
          ENDIF
          IF(ALLOCATED(XUMEM))THEN
            IF(IDIFK == 2)THEN
              WHERE((XU(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))==XSPVAL)&
	      .OR. (XUMEM(:,:,INIV1,NTIMEDIA(JLOOPT,K-1,INDKM1),1, &
	      NPROCDIA(1,K-1))==XSPVAL))
	        XU(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))=XSPVAL
              ELSEWHERE
  	        XU(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K)) =  &
  	        XUMEM(:,:,INIV1,NTIMEDIA(JLOOPT,K-1,INDKM1),1,NPROCDIA(1,K-1)) + &
  	        IFAC * XU(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))
              END WHERE
            ELSE
              WHERE((XU(:,:,:,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))==XSPVAL)&
	      .OR. (XUMEM(:,:,:,NTIMEDIA(JLOOPT,K-1,INDKM1),1, &
	      NPROCDIA(1,K-1))==XSPVAL))
	        XU(:,:,:,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))=XSPVAL
              ELSEWHERE
  	        XU(:,:,:,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K)) =  &
  	        XUMEM(:,:,:,NTIMEDIA(JLOOPT,K-1,INDKM1),1,NPROCDIA(1,K-1)) + &
  	        IFAC * XU(:,:,:,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))
              END WHERE
            ENDIF
          ENDIF
          IF(ALLOCATED(XVMEM))THEN
            IF(IDIFK == 2)THEN
              WHERE((XV(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))==XSPVAL)&
	        .OR. (XVMEM(:,:,INIV1,NTIMEDIA(JLOOPT,K-1,INDKM1),1, &
	        NPROCDIA(1,K-1))==XSPVAL))
	        XV(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))=XSPVAL
              ELSEWHERE
  	        XV(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K)) =  &
  	        XVMEM(:,:,INIV1,NTIMEDIA(JLOOPT,K-1,INDKM1),1,NPROCDIA(1,K-1)) + &
  	        IFAC * XV(:,:,INIV2,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))
              END WHERE
            ELSE
              WHERE((XV(:,:,:,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))==XSPVAL)&
	        .OR. (XVMEM(:,:,:,NTIMEDIA(JLOOPT,K-1,INDKM1),1, &
	        NPROCDIA(1,K-1))==XSPVAL))
	        XV(:,:,:,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))=XSPVAL
              ELSEWHERE
  	        XV(:,:,:,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K)) =  &
  	        XVMEM(:,:,:,NTIMEDIA(JLOOPT,K-1,INDKM1),1,NPROCDIA(1,K-1)) + &
  	        IFAC * XV(:,:,:,NTIMEDIA(JLOOPT,K,INDK),1,NPROCDIA(1,K))
              END WHERE
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ELSE
! Intervalle de temps de longueur differente  pour les 2 champs. On ne trace pas
      LPBREAD=.TRUE.
      print *,' INTERVALLE DE TEMPS DIFFERENT POUR LES 2 CHAMPS. PAS DE TRACE '
      IF(LFT .OR. LPVKT .OR. LFT1 .OR. LPVKT1 .AND. LSPVALT)THEN
        XSPVAL=ZSPVAL
      ENDIF
      RETURN
    ENDIF

!           A VERIFIER                SSSSSSSSSSSS
  ELSE
! Temps incremental
    IF(NBPROCDIA(K) /= 1 .OR. NBPROCDIA(K-1) /= 1 )THEN
      LPBREAD=.TRUE.
      print *,' NB DE PROCESSUS DEMANDES > 1 . PAS DE TRACE '
      IF(LFT .OR. LPVKT .OR. LFT1 .OR. LPVKT1 .AND. LSPVALT)THEN
        XSPVAL=ZSPVAL
      ENDIF
      RETURN
    ELSE
      ILENT=(NTIMEDIA(2,K,INDK)-NTIMEDIA(1,K,INDK))/NTIMEDIA(3,K,INDK)
      ILENT2=(NTIMEDIA(2,K-1,INDKM1)-NTIMEDIA(1,K-1,INDKM1))/NTIMEDIA(3,K-1,INDKM1)
      IF(ILENT2 /= ILENT)THEN
        LPBREAD=.TRUE.
        print *,' INTERVALLE DE TEMPS DIFFERENT POUR LES 2 CHAMPS. PAS DE TRACE '
	print *,' (',XTIMEDIA(1,K-1,INDKM1),' - (',XTIMEDIA(2,K-1,INDKM1),') ET (', &
	XTIMEDIA(1,K,INDK),' - (',XTIMEDIA(2,K,INDK),')'
        IF(LFT .OR. LPVKT .OR. LFT1 .OR. LPVKT1 .AND. LSPVALT)THEN
          XSPVAL=ZSPVAL
        ENDIF
        RETURN
      ELSE
!-----------------------------------------------------------------------
!       IF(K <= 2)THEN
        IF(NUMPM(K-1) == 0 .OR. NUMPM(MAX(K-1,1)) == 3)THEN
          if(nverbia > 0)then
            print *,' diff_oper3 K NUMPM(K-1) ',K,NUMPM(K-1)
          endif
      if(nblvlkdia(k,indk) == 1 .and. nblvlkdia(k-1,indk) == 1)then
        print *,' diff_oper3 Niveaux en K de part et d autre de MINUS '
        print *,' K1= ',NLVLKDIA(1,K-1,INDK),' K2= ',NLVLKDIA(1,K,INDK)
! INIV1=le 1er dans la directive <-> K-1 , INIV2=le 2e=courant <-> K
! Janv 2001
	IF(CTYPE == 'CART' .OR. CTYPE == 'MASK')THEN
          INIV1=NLVLKDIA(1,K-1,INDK)-NKL+1
          INIV2=NLVLKDIA(1,K,INDK)-NKL+1
	ELSE
! Janv 2001
          INIV1=NLVLKDIA(1,K-1,INDK)
          INIV2=NLVLKDIA(1,K,INDK)
! Janv 2001
        ENDIF
! Janv 2001
        IDIFK=2
      endif


          IF(LMUMVM .OR. LMUTVT .OR. LUMVM .OR. LUTVT .OR. &
            LSUMVM .OR. LSUTVT .OR. LMLSUMVM .OR. LMLSUTVT)THEN
            CGROUPS(K-1)=ADJUSTL(CGROUPS(K-1))
	    ILENTIT=LEN_TRIM(CGROUPS(K-1))
            YCAR(ILENC:ILENC+ILENTIT-1)=CGROUPS(K-1)(1:ILENTIT)
          ELSE
            CTITRE2(NPROCDIA(1,K-1))=ADJUSTL(CTITRE2(NPROCDIA(1,K-1)))
	    ILENTIT=LEN_TRIM(CTITRE2(NPROCDIA(1,K-1)))
            YCAR(ILENC:ILENC+ILENTIT-1)=CTITRE2(NPROCDIA(1,K-1))(1:ILENTIT)
          ENDIF
          YCAR=ADJUSTL(YCAR)
          ILENC=LEN_TRIM(YCAR)
          ILENC=ILENC+2
          YCAR(ILENC:ILENC)='('
          ILENC=ILENC+1
          IF(NUMFILECUR2 /= NUMFILECUR)THEN
	    CFILEDIAS(JME2)=ADJUSTL(CFILEDIAS(JME2))
	    ILENFI=LEN_TRIM(CFILEDIAS(JME2))
            YCAR(ILENC:ILENC+ILENFI-1)=CFILEDIAS(JME2)(1:ILENFI)
            ILENC=ILENC+ILENFI
            YCAR(ILENC:ILENC+1)=')('
            ILENC=ILENC+2
          ENDIF

      IF(IDIFK == 2 .AND. INIV1 /= INIV2) THEN
        YCAR(ILENC:ILENC+1)='K='
        ILENC=ILENC+2
        WRITE(YCAR(ILENC:ILENC+1),'(I2)')INIV1
        ILENC=ILENC+2
        YCAR(ILENC:ILENC+1)=')('
        ILENC=ILENC+2
      ENDIF

! Ecriture de la premiere serie de temps
          IF(CTYPE == 'DRST' .OR. CTYPE == 'RSPL' .OR. CTYPE == 'RAPL')THEN
            WRITE(YTIM,'(F8.0)')XTRAJT2(NTIMEDIA(1,K-1,INDKM1),INDKM1)
          ELSE
            WRITE(YTIM,'(F8.0)')XTRAJT2(NTIMEDIA(1,K-1,INDKM1),1)
          ENDIF
          YTIM=ADJUSTL(YTIM)
          ILENTIM=LEN_TRIM(YTIM)
          YCAR(ILENC:ILENC+ILENTIM-1)=YTIM(1:ILENTIM)
          YTIM(1:LEN(YTIM))=' '
          ILENC=ILENC+ILENTIM

          IF(CTYPE == 'DRST' .OR. CTYPE == 'RSPL' .OR. CTYPE == 'RAPL')THEN
            INDNN=INDKM1
          ELSE
            INDNN=1
          ENDIF
!         IF(XTRAJT2(NTIMEDIA(2,K-1,INDKM1),1) /= XTRAJT2(NTIMEDIA(1,K-1,INDKM1),1))THEN
          IF(XTRAJT2(NTIMEDIA(2,K-1,INDKM1),INDNN) /= XTRAJT2(NTIMEDIA(1,K-1,INDKM1),INDNN))THEN

	    YCAR(ILENC:ILENC+2)=' - '
	    ILENC=ILENC+3
            WRITE(YTIM,'(F8.0)')XTRAJT2(NTIMEDIA(2,K-1,INDKM1),INDNN)
!           WRITE(YTIM,'(F8.0)')XTRAJT2(NTIMEDIA(2,K-1,INDKM1),1)
            YTIM=ADJUSTL(YTIM)
            ILENTIM=LEN_TRIM(YTIM)
            YCAR(ILENC:ILENC+ILENTIM-1)=YTIM(1:ILENTIM)
            YTIM(1:LEN(YTIM))=' '
            ILENC=ILENC+ILENTIM

          ELSE
            ILENC=ILENC+1
          ENDIF

          IF(LMINUS)THEN
	    YCAR(ILENC:ILENC+3)=') - '
          ELSE IF(LPLUS)THEN
	    YCAR(ILENC:ILENC+3)=') + '
          ENDIF
	  ILENC=ILENC+4

        ELSE
    
          YCAR=ADJUSTL(CTITB3)
          ILENC=LEN_TRIM(YCAR)+1
          IF(LMINUS)THEN
            YCAR(ILENC:ILENC+2)=' - '
          ELSE IF(LPLUS)THEN
            YCAR(ILENC:ILENC+2)=' + '
          ENDIF
          ILENC=ILENC+3
    
        ENDIF
!-----------------------------------------------------------------------
!  FIN Ecriture de la premiere serie de temps
        IF(LMUMVM .OR. LMUTVT .OR. LUMVM .OR. LUTVT .OR. &
           LSUMVM .OR. LSUTVT .OR. LMLSUMVM .OR. LMLSUTVT)THEN
          ILENTIT=LEN_TRIM(CGROUPS(K))
          YCAR(ILENC:ILENC+ILENTIT-1)=CGROUPS(K)(1:ILENTIT)
        ELSE
          CTITRE(NPROCDIA(1,K))=ADJUSTL(CTITRE(NPROCDIA(1,K)))
          ILENTIT=LEN_TRIM(CTITRE(NPROCDIA(1,K)))
          YCAR(ILENC:ILENC+ILENTIT-1)=CTITRE(NPROCDIA(1,K))(1:ILENTIT)
        ENDIF
        YCAR=ADJUSTL(YCAR)
        ILENC=LEN_TRIM(YCAR)
        ILENC=ILENC+2
        YCAR(ILENC:ILENC)='('
        ILENC=ILENC+1
        IF(NUMFILECUR2 /= NUMFILECUR)THEN
	  CFILEDIAS(JME)=ADJUSTL(CFILEDIAS(JME))
	  ILENFI=LEN_TRIM(CFILEDIAS(JME))
          YCAR(ILENC:ILENC+ILENFI-1)=CFILEDIAS(JME)(1:ILENFI)
          ILENC=ILENC+ILENFI
          YCAR(ILENC:ILENC+1)=')('
          ILENC=ILENC+2
        ENDIF

      IF(IDIFK == 2 .AND. INIV1 /= INIV2)THEN
        YCAR(ILENC:ILENC+1)='K='
        ILENC=ILENC+2
        WRITE(YCAR(ILENC:ILENC+1),'(I2)')INIV2
        ILENC=ILENC+2
        YCAR(ILENC:ILENC+1)=')('
        ILENC=ILENC+2
      ENDIF


! Ecriture de la deuxieme serie de temps
        IF(CTYPE == 'DRST' .OR. CTYPE == 'RSPL' .OR. CTYPE == 'RAPL')THEN
          WRITE(YTIM,'(F8.0)')XTRAJT(NTIMEDIA(1,K,INDK),INDK)
        ELSE
          WRITE(YTIM,'(F8.0)')XTRAJT(NTIMEDIA(1,K,INDK),1)
        ENDIF
        YTIM=ADJUSTL(YTIM)
        ILENTIM=LEN_TRIM(YTIM)
        YCAR(ILENC:ILENC+ILENTIM-1)=YTIM(1:ILENTIM)
        YTIM(1:LEN(YTIM))=' '
        ILENC=ILENC+ILENTIM

        IF(CTYPE == 'DRST' .OR. CTYPE == 'RSPL' .OR. CTYPE == 'RAPL')THEN
          INDNN=INDK
        ELSE
          INDNN=1
        ENDIF
!       IF(XTRAJT(NTIMEDIA(2,K,INDK),1) /= XTRAJT(NTIMEDIA(1,K,INDK),1))THEN
        IF(XTRAJT(NTIMEDIA(2,K,INDK),INDNN) /= XTRAJT(NTIMEDIA(1,K,INDK),INDNN))THEN

	  YCAR(ILENC:ILENC+2)=' - '
	  ILENC=ILENC+3
          WRITE(YTIM,'(F8.0)')XTRAJT(NTIMEDIA(2,K,INDK),INDNN)
!         WRITE(YTIM,'(F8.0)')XTRAJT(NTIMEDIA(2,K,INDK),1)
          YTIM=ADJUSTL(YTIM)
          ILENTIM=LEN_TRIM(YTIM)
          YCAR(ILENC:ILENC+ILENTIM-1)=YTIM(1:ILENTIM)
          YTIM(1:LEN(YTIM))=' '
          ILENC=ILENC+ILENTIM

        ELSE
          ILENC=ILENC+1
        ENDIF

	YCAR(ILENC:ILENC)=')'
        
!       IF(K <= 2)THEN
        IF(NUMPM(K-1) == 0 .OR. NUMPM(MAX(K-1,1)) == 3)THEN
          if(nverbia > 0)then
            print *,' diff_oper3 K NUMPM(K-1) ',K,NUMPM(K-1)
          endif
          LTITDEF=.FALSE.
!         YTITB3(1:LEN(YTITB3))=' '
!         YTITB3=ADJUSTL(CTITB3)
        ENDIF
        if(nblvlkdia(k,indk) == 1 .and. nblvlkdia(k-1,indk) == 1)then
          print *,' diff_oper3-2 Niveaux en K de part et d autre de MINUS '
          print *,' K1= ',NLVLKDIA(1,K-1,INDK),' K2= ',NLVLKDIA(1,K,INDK)
! Janv 2001
	  IF(CTYPE == 'CART' .OR. CTYPE == 'MASK')THEN
            INIV1=NLVLKDIA(1,K-1,INDK)-NKL+1
            INIV2=NLVLKDIA(1,K,INDK)-NKL+1
	  ELSE
! Janv 2001
            INIV1=NLVLKDIA(1,K-1,INDK)
            INIV2=NLVLKDIA(1,K,INDK)
! Janv 2001
          ENDIF
! Janv 2001
            if(nverbia >0)then
            print *,' INIV1 INIV2 diff_oper',INIV1,INIV2
            endif
          IDIFK=2
        endif
!!! 1/3/04
    IF(CTITB3 /= 'DEFAULT')THEN
    ELSE
!!! 1/3/04

        CTITB3=ADJUSTL(YCAR(1:100))
        CTITB3=ADJUSTL(CTITB3)
!!! 1/3/04
    ENDIF
!!! 1/3/04
        print *,' CTITB3 ',CTITB3
        IF(LMINUS)THEN
          IFAC=-1
        ELSE IF(LPLUS)THEN
          IFAC=1
        ENDIF
! 220900
        IF(IDIFK == 2)THEN
          if(nverbia > 0)then
            print *,' **diff_oper IDIFK size(XVAR) ',IDIFK,size(xvar,1),&
            size(xvar,2),size(xvar,3),size(xvar,4),size(xvar,5),size(xvar,6)
            print *,' AV ',xvar(1,1,1,:,1,3)
            print *,' INDKM1 ',INDKM1,' K ',K,' NPROCDIA(1,K) ',NPROCDIA(1,K)
          endif
! Mai 2003
          WHERE((XVAR(:,:,INIV2,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	  NTIMEDIA(3,K,INDK),INDK,NPROCDIA(1,K)) == XSPVAL) &
!         NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K)) == XSPVAL) &
	  .OR. (XVAR2(:,:,INIV1,NTIMEDIA(1,K-1,INDKM1): &
	  NTIMEDIA(2,K-1,INDKM1):&
	  NTIMEDIA(3,K-1,INDKM1),INDKM1,NPROCDIA(1,K-1)) == XSPVAL))
!         NTIMEDIA(3,K-1,INDKM1),1,NPROCDIA(1,K-1)) == XSPVAL))
	    XVAR(:,:,INIV2,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
            NTIMEDIA(3,K,INDK),INDK,NPROCDIA(1,K))=XSPVAL
!           NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))=XSPVAL
          ELSEWHERE
	    XVAR(:,:,INIV2,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	    NTIMEDIA(3,K,INDK),INDK,NPROCDIA(1,K))= &
!           NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))= &
	    XVAR2(:,:,INIV1,NTIMEDIA(1,K-1,INDKM1):NTIMEDIA(2,K-1,INDKM1):&
	    NTIMEDIA(3,K-1,INDKM1),INDKM1,NPROCDIA(1,K-1)) +  &
!           NTIMEDIA(3,K-1,INDKM1),1,NPROCDIA(1,K-1)) +  &
	    IFAC * XVAR(:,:,INIV2,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	    NTIMEDIA(3,K,INDK),INDK,NPROCDIA(1,K))
!           NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))
          END WHERE
          if(nverbia > 0)then
            print *,' AP ',xvar(1,1,1,:,1,3)
          endif
	ELSE
          WHERE((XVAR(:,:,:,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	  NTIMEDIA(3,K,INDK),INDK,NPROCDIA(1,K)) == XSPVAL) &
!         NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K)) == XSPVAL) &
	  .OR. (XVAR2(:,:,:,NTIMEDIA(1,K-1,INDKM1): &
	  NTIMEDIA(2,K-1,INDKM1):&
	  NTIMEDIA(3,K-1,INDKM1),INDKM1,NPROCDIA(1,K-1)) == XSPVAL))
!         NTIMEDIA(3,K-1,INDKM1),1,NPROCDIA(1,K-1)) == XSPVAL))
	    XVAR(:,:,:,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
            NTIMEDIA(3,K,INDK),INDK,NPROCDIA(1,K))=XSPVAL
!           NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))=XSPVAL
          ELSEWHERE
	    XVAR(:,:,:,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	    NTIMEDIA(3,K,INDK),INDK,NPROCDIA(1,K))= &
!           NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))= &
	    XVAR2(:,:,:,NTIMEDIA(1,K-1,INDKM1):NTIMEDIA(2,K-1,INDKM1):&
	    NTIMEDIA(3,K-1,INDKM1),INDKM1,NPROCDIA(1,K-1)) +  &
!           NTIMEDIA(3,K-1,INDKM1),1,NPROCDIA(1,K-1)) +  &
	    IFAC * XVAR(:,:,:,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	    NTIMEDIA(3,K,INDK),INDK,NPROCDIA(1,K))
!           NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))
          END WHERE
        ENDIF
        IF(ALLOCATED(XUMEM))THEN
          IF(IDIFK == 2)THEN
            WHERE((XU(:,:,INIV2,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	    NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K)) == XSPVAL) &
	    .OR. (XUMEM(:,:,INIV1,NTIMEDIA(1,K-1,INDKM1): &
	    NTIMEDIA(2,K-1,INDKM1):&
	    NTIMEDIA(3,K-1,INDKM1),1,NPROCDIA(1,K-1)) == XSPVAL))
	      XU(:,:,INIV2,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
              NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))=XSPVAL
            ELSEWHERE
	      XU(:,:,INIV2,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	      NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))= &
	      XUMEM(:,:,INIV1,NTIMEDIA(1,K-1,INDKM1):NTIMEDIA(2,K-1,INDKM1):&
	      NTIMEDIA(3,K-1,INDKM1),1,NPROCDIA(1,K-1)) +  &
	      IFAC * XU(:,:,INIV2,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	      NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))
            END WHERE
          ELSE
            WHERE((XU(:,:,:,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	    NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K)) == XSPVAL) &
	    .OR. (XUMEM(:,:,:,NTIMEDIA(1,K-1,INDKM1): &
	    NTIMEDIA(2,K-1,INDKM1):&
	    NTIMEDIA(3,K-1,INDKM1),1,NPROCDIA(1,K-1)) == XSPVAL))
	      XU(:,:,:,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
              NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))=XSPVAL
            ELSEWHERE
	      XU(:,:,:,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	      NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))= &
	      XUMEM(:,:,:,NTIMEDIA(1,K-1,INDKM1):NTIMEDIA(2,K-1,INDKM1):&
	      NTIMEDIA(3,K-1,INDKM1),1,NPROCDIA(1,K-1)) +  &
	      IFAC * XU(:,:,:,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	      NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))
            END WHERE
	  ENDIF
	ENDIF
        IF(ALLOCATED(XVMEM))THEN
          IF(IDIFK == 2)THEN
            WHERE((XV(:,:,INIV2,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	    NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K)) == XSPVAL) &
	    .OR. (XVMEM(:,:,INIV1,NTIMEDIA(1,K-1,INDKM1): &
	    NTIMEDIA(2,K-1,INDKM1):&
	    NTIMEDIA(3,K-1,INDKM1),1,NPROCDIA(1,K-1)) == XSPVAL))
	      XV(:,:,INIV2,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
              NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))=XSPVAL
            ELSEWHERE
	      XV(:,:,INIV2,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	      NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))= &
	      XVMEM(:,:,INIV1,NTIMEDIA(1,K-1,INDKM1):NTIMEDIA(2,K-1,INDKM1):&
	      NTIMEDIA(3,K-1,INDKM1),1,NPROCDIA(1,K-1)) +  &
	      IFAC * XV(:,:,INIV2,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	      NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))
            END WHERE
          ELSE
            WHERE((XV(:,:,:,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	    NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K)) == XSPVAL) &
	    .OR. (XVMEM(:,:,:,NTIMEDIA(1,K-1,INDKM1): &
	    NTIMEDIA(2,K-1,INDKM1):&
	    NTIMEDIA(3,K-1,INDKM1),1,NPROCDIA(1,K-1)) == XSPVAL))
	      XV(:,:,:,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
              NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))=XSPVAL
            ELSEWHERE
	      XV(:,:,:,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	      NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))= &
	      XVMEM(:,:,:,NTIMEDIA(1,K-1,INDKM1):NTIMEDIA(2,K-1,INDKM1):&
	      NTIMEDIA(3,K-1,INDKM1),1,NPROCDIA(1,K-1)) +  &
	      IFAC * XV(:,:,:,NTIMEDIA(1,K,INDK):NTIMEDIA(2,K,INDK): &
	      NTIMEDIA(3,K,INDK),1,NPROCDIA(1,K))
            END WHERE
	  ENDIF
	ENDIF

      ENDIF

    ENDIF
!           A VERIFIER                SSSSSSSSSSSS
  ENDIF

!******************************** A VERIFIER
ENDIF

!
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
!IF(K == NSUPERDIA .AND. YTITB3 /= ' ' .AND. YTITB3 /= 'DEFAULT')THEN
! CTITB3=YTITB3
!ENDIF

IF(LFT .OR. LPVKT .OR. LFT1 .OR. LPVKT1 .AND. LSPVALT)THEN
  XSPVAL=ZSPVAL
ENDIF
RETURN
END SUBROUTINE DIFF_OPER
