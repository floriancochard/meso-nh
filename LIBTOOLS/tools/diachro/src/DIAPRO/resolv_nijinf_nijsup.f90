!     ######spl
      SUBROUTINE RESOLV_NIJINF_NIJSUP
!     ###############################
!
!!****  *RESOLV_NIJINF_NIJSUP* -  Affectation des valeurs de NIINF, NISUP,
!!                                NJINF et NJSUP dans les 2 cas CH et CV
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      None
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       16/01/95
!!      Updated   PM 
!-------------------------------------------------------------------------------
USE MODD_DIM1
USE MODD_TYPE_AND_LH
USE MODD_PARAMETERS
USE MODD_RESOLVCAR
USE MODN_PARA
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CH Positionnement NIINF, NJINF, NISUP, NJSUP
! Defaut : NIINF=NIL, NJINF=NJL, NISUP=NIH, NJSUP=NJH
! Sinon valeurs fournies par l'utilisateur dans les limites (NIL,NJL NIH,
! NJH)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		if(nverbia > 0)then
                  print *,' **resolv_ni... LCH LCV LCHXY LFT LPVKT ',LCH,LCV,LCHXY,LFT,LPVKT
                endif

		IF((LCH.AND..NOT.LCV) .OR. (LCHXY.AND..NOT.LCV))THEN

		  IF(NIINF == 0)THEN
		    NIINF=NIL
		    IF(NIINF == 1)NIINF=NIINF+JPHEXT
                  ELSE IF(NIINF /=0)THEN
		    print *,' NIINF DEMANDE NIL NIH ', &
				NIINF,NIL,NIH
		    IF(NIINF < NIL .OR. NIINF > NIH)THEN
		      NIINF=NIL
		      IF(NIINF == 1)THEN
		        NIINF=NIINF+JPHEXT
		      ENDIF
                      print *,' NIINF MODIFIE ', NIINF
		    ENDIF
		  ENDIF

		  IF(NJINF == 0)THEN
		    NJINF=NJL
		    IF(NJINF == 1)NJINF=NJINF+JPHEXT
                  ELSE IF(NJINF /=0)THEN
		    print *,' NJINF DEMANDE NJL NJH ', &
				NJINF,NJL,NJH
		    IF(NJINF < NJL .OR. NJINF > NJH)THEN
		      NJINF=NJL
		      IF(NJINF == 1)THEN
		        NJINF=NJINF+JPHEXT
		      ENDIF
                      print *,' NJINF MODIFIE ', NJINF
		    ENDIF
		  ENDIF

		  IF(NISUP == 0)THEN
		    NISUP=NIH
		    IF(NISUP > NIMAX+JPHEXT)NISUP=NIMAX+JPHEXT
                  ELSE IF(NISUP /=0)THEN
		    print *,' NISUP DEMANDE NIL NIH ', &
			NISUP,NIL,NIH
		    IF(NISUP < NIL .OR. NISUP > NIH)THEN
		      NISUP=NIH
		      IF(NISUP > NIMAX+JPHEXT)THEN
			NISUP=NIMAX+JPHEXT
		      ENDIF
                      print *,' NISUP MODIFIE ', NISUP
		    ENDIF
		  ENDIF

		  IF(NJSUP == 0)THEN
		    NJSUP=NJH
		    IF(NJSUP > NJMAX+JPHEXT)NJSUP=NJMAX+JPHEXT
                  ELSE IF(NJSUP /=0)THEN
		    print *,' NJSUP DEMANDE NJL NJH ', &
		               NJSUP,NJL,NJH
		    IF(NJSUP < NJL .OR. NJSUP > NJH)THEN
		      NJSUP=NJH
		      IF(NJSUP > NJMAX+JPHEXT)THEN
			NJSUP=NJMAX+JPHEXT
		      ENDIF
                      print *,' NJSUP MODIFIE ', NJSUP
		    ENDIF
		  ENDIF

!           	  print *,' NIINF,NISUP,NJINF,NJSUP,NKL,NKH ', &
!   		            NIINF,NISUP,NJINF,NJSUP,NKL,NKH
    
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CV Positionnement NIINF, NJINF, NISUP, NJSUP
! CV Positionnement LHORIZ et LVERTI
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                ELSE IF(LCV)THEN
                  LHORIZ=.FALSE.; LVERTI=.TRUE.
		  NIINF=NIL
		  NJINF=NJL
		  NISUP=NIH
		  NJSUP=NJH

!           	  print *,' NIINF,NISUP,NJINF,NJSUP,NKL,NKH ', &
!   		            NIINF,NISUP,NJINF,NJSUP,NKL,NKH
		ENDIF
    
		if(nverbia > 0)then
                  print *,' **resolv_nii.. NIINF,NISUP,NJINF,NJSUP,NIL,NIH,NJL,NJH,NKL,NKH ', &
    		            NIINF,NISUP,NJINF,NJSUP,NIL,NIH,NJL,NJH,NKL,NKH
                endif
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
!
RETURN
END SUBROUTINE  RESOLV_NIJINF_NIJSUP
