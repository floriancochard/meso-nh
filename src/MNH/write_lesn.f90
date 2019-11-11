!MNH_LIC Copyright 2000-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!######################
MODULE MODI_WRITE_LES_n
!######################
!
INTERFACE
!
      SUBROUTINE  WRITE_LES_n(TPDIAFILE,HLES_AVG)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA),  INTENT(IN) :: TPDIAFILE! file to write
CHARACTER(LEN=1), INTENT(IN) :: HLES_AVG ! flag to perform the averages
!                                        ! or normalizations
END SUBROUTINE WRITE_LES_n
!
END INTERFACE
!
END MODULE MODI_WRITE_LES_n

!     ######################
      SUBROUTINE  WRITE_LES_n(TPDIAFILE,HLES_AVG)
!     ######################
!
!
!!****  *WRITE_LES_n* writes the LES final diagnostics for model _n 
!!                         
!!
!!    PURPOSE
!!    -------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         07/02/00
!!                       01/02/01 (D. Gazen) add module MODD_NSV for NSV variable
!!                       06/11/02 (V. Masson) some minor bugs
!!                       01/04/03 (V. Masson) idem
!!                       10/10/09 (P. Aumond) Add user multimaskS
!!                          11/15 (C.Lac) Add production terms of TKE
!!                    10/2016 (C.Lac) Add droplet deposition
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!!!                     02/2019 (C. Lac) Add rain fraction as a LES diagnostic

!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LES
USE MODD_LES_n
USE MODD_FIELD_n
USE MODD_CONF_n
USE MODD_PARAM_n
USE MODD_TURB_n
USE MODD_GRID_n
USE MODD_NSV, ONLY : NSV
USE MODD_PARAM_ICE, ONLY : LDEPOSC
USE MODD_PARAM_C2R2, ONLY : LDEPOC
!
USE MODE_ll
!
USE MODE_LES_DIACHRO
USE MODE_LES_SPEC_N
USE MODE_MODELN_HANDLER
!
USE MODI_WRITE_LES_BUDGET_n
USE MODI_WRITE_LES_RT_BUDGET_n
USE MODI_WRITE_LES_SV_BUDGET_n
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(TFILEDATA),  INTENT(IN) :: TPDIAFILE! file to write
CHARACTER(LEN=1), INTENT(IN) :: HLES_AVG ! flag to perform the averages
!                                        ! or normalizations
!
!
!*      0.2  declaration of local variables
!
INTEGER :: IMASK
!
INTEGER :: JSV       ! scalar loop counter
INTEGER :: JI        ! loop counter
INTEGER :: JPDF      ! pdf loop counter
!
CHARACTER(len=9), DIMENSION(NLES_MASKS) :: YSUBTITLE
CHARACTER(len=5)                        :: YGROUP
!
REAL, DIMENSION(:,:,:), ALLOCATABLE     :: ZAVG_PTS_ll
REAL, DIMENSION(:,:,:), ALLOCATABLE     :: ZUND_PTS_ll
REAL                                    :: ZCART_PTS_ll
INTEGER                                 :: IMI ! Current model inde
!
!
!-------------------------------------------------------------------------------
!
IF (.NOT. LLES) RETURN
!
IF (HLES_AVG=='A'                                                       &
     .AND. (XLES_TEMP_MEAN_START==XUNDEF .OR. XLES_TEMP_MEAN_END==XUNDEF)) RETURN
IF (HLES_AVG=='E' .AND. CLES_NORM_TYPE=='NONE'                          ) RETURN
IF (HLES_AVG=='H' .AND. (CLES_NORM_TYPE=='NONE'                          &
     .OR. XLES_TEMP_MEAN_START/=XUNDEF .OR. XLES_TEMP_MEAN_END/=XUNDEF)) RETURN
!
!*      1.   Initializations
!            ---------------
!
IMI = GET_CURRENT_MODEL_INDEX()
!
!
!*      1.1  Normalization variables
!            -----------------------
!
IF (CLES_NORM_TYPE/='NONE' ) THEN
  ALLOCATE(XLES_NORM_M  (NLES_TIMES))
  ALLOCATE(XLES_NORM_S  (NLES_TIMES))
  ALLOCATE(XLES_NORM_K  (NLES_TIMES))
  ALLOCATE(XLES_NORM_RHO(NLES_TIMES))
  ALLOCATE(XLES_NORM_RV (NLES_TIMES))
  ALLOCATE(XLES_NORM_SV (NLES_TIMES,NSV))
  ALLOCATE(XLES_NORM_P  (NLES_TIMES))
  !
  IF (CLES_NORM_TYPE=='CONV') THEN
    WHERE (XLES_WSTAR(:)>0.)
      XLES_NORM_M(:)   = XLES_BL_HEIGHT(:)
      XLES_NORM_S(:)   = XLES_NORM_M(:) / XLES_WSTAR(:)
      XLES_NORM_K(:)   = XLES_Q0(:) / XLES_WSTAR(:)
      XLES_NORM_RHO(:) = XLES_MEAN_RHO(1,:,1)
      XLES_NORM_RV(:)  = XLES_E0(:) / XLES_WSTAR(:)
      XLES_NORM_P(:)   = XLES_MEAN_RHO(1,:,1) * XLES_WSTAR(:)**2
    ELSEWHERE
      XLES_NORM_M(:)   = 0.
      XLES_NORM_S(:)   = 0.
      XLES_NORM_K(:)   = 0.
      XLES_NORM_RHO(:) = 0.
      XLES_NORM_RV(:)  = 0.
      XLES_NORM_P(:)   = 0.
    END WHERE
    DO JSV=1,NSV
      WHERE (XLES_WSTAR(:)>0.)
        XLES_NORM_SV(:,JSV)= XLES_SV0(:,JSV) / XLES_WSTAR(:)
      ELSEWHERE
        XLES_NORM_SV(:,JSV)= 0.
      END WHERE
    END DO
  ELSE IF (CLES_NORM_TYPE=='EKMA') THEN
    WHERE (XLES_USTAR(:)>0.)
      XLES_NORM_M(:)   = XLES_BL_HEIGHT(:)
      XLES_NORM_S(:)   = XLES_NORM_M(:) / XLES_USTAR(:)
      XLES_NORM_K(:)   = XLES_Q0(:) / XLES_USTAR(:)
      XLES_NORM_RHO(:) = XLES_MEAN_RHO(1,:,1)
      XLES_NORM_RV(:)  = XLES_E0(:) / XLES_USTAR(:)
      XLES_NORM_P(:)   = XLES_MEAN_RHO(1,:,1) * XLES_USTAR(:)**2
    ELSEWHERE
      XLES_NORM_M(:)   = 0.
      XLES_NORM_S(:)   = 0.
      XLES_NORM_K(:)   = 0.
      XLES_NORM_RHO(:) = 0.
      XLES_NORM_RV(:)  = 0.
      XLES_NORM_P(:)   = 0.
    END WHERE
    DO JSV=1,NSV
      WHERE (XLES_USTAR(:)>0.)
        XLES_NORM_SV(:,JSV)= XLES_SV0(:,JSV) / XLES_USTAR(:)
      ELSEWHERE
        XLES_NORM_SV(:,JSV)= 0.
      END WHERE
    END DO
  ELSE IF (CLES_NORM_TYPE=='MOBU') THEN
    XLES_NORM_M(:) = XLES_MO_LENGTH(:)
    WHERE (XLES_USTAR(:)>0.)
      XLES_NORM_S(:)   = XLES_NORM_M(:) / XLES_USTAR(:)
      XLES_NORM_K(:)   = XLES_Q0(:) / XLES_USTAR(:)
      XLES_NORM_RHO(:) = XLES_MEAN_RHO(1,:,1)
      XLES_NORM_RV(:)  = XLES_E0(:) / XLES_USTAR(:)
      XLES_NORM_P(:)   = XLES_MEAN_RHO(1,:,1) * XLES_USTAR(:)**2
    ELSEWHERE
      XLES_NORM_S(:)   = 0.
      XLES_NORM_K(:)   = 0.
      XLES_NORM_RHO(:) = 0.
      XLES_NORM_RV(:)  = 0.
      XLES_NORM_P(:)   = 0.
    END WHERE
    DO JSV=1,NSV
      WHERE (XLES_USTAR(:)>0.)
        XLES_NORM_SV(:,JSV)= XLES_SV0(:,JSV) / XLES_USTAR(:)
      ELSEWHERE
        XLES_NORM_SV(:,JSV)= 0.
      END WHERE
    END DO
  END IF
END IF
!
!*      1.2  Initializations for WRITE_DIACHRO
!            ---------------------------------
!
NLES_CURRENT_TIMES=NLES_TIMES
!
ALLOCATE(XLES_CURRENT_TRAJT(NLES_TIMES,1))
XLES_CURRENT_TRAJT(:,:) = XLES_TRAJT(:,:)
ALLOCATE(XLES_CURRENT_Z(NLES_K))
XLES_CURRENT_Z(:) = XLES_Z(:)
ALLOCATE(XLES_CURRENT_DATIME(16,NLES_TIMES))
XLES_CURRENT_DATIME(:,:) = XLES_DATIME(:,:)
!
XLES_CURRENT_ZS = XLES_ZS
!
NLES_CURRENT_IINF=NLESn_IINF(IMI)
NLES_CURRENT_ISUP=NLESn_ISUP(IMI)
NLES_CURRENT_JINF=NLESn_JINF(IMI)
NLES_CURRENT_JSUP=NLESn_JSUP(IMI)
!
XLES_CURRENT_DOMEGAX=XDXHAT(1)
XLES_CURRENT_DOMEGAY=XDYHAT(1)
!
!
!
!*      2.   (z,t) profiles (all masks)
!            --------------
IMASK = 1
YSUBTITLE(IMASK) = " (cart)"
IF (LLES_NEB_MASK) THEN
  IMASK=IMASK+1
  YSUBTITLE(IMASK) = " (neb)"
  IMASK=IMASK+1
  YSUBTITLE(IMASK) = " (clear)"
END IF
IF (LLES_CORE_MASK) THEN
  IMASK=IMASK+1
  YSUBTITLE(IMASK) = " (core)"
  IMASK=IMASK+1
  YSUBTITLE(IMASK) = " (env)"
END IF
IF (LLES_MY_MASK) THEN
   DO JI=1,NLES_MASKS_USER
        IMASK=IMASK+1
        YSUBTITLE(IMASK) = " (user)"
   END DO
END IF
IF (LLES_CS_MASK) THEN
  IMASK=IMASK+1
  YSUBTITLE(IMASK) = " (cs1)"
  IMASK=IMASK+1
  YSUBTITLE(IMASK) = " (cs2)"
  IMASK=IMASK+1
  YSUBTITLE(IMASK) = " (cs3)"
END IF
!
!*      2.0  averaging diagnostics
!            ---------------------
!
IF (HLES_AVG==' ' .OR. HLES_AVG=='A') THEN
  ALLOCATE(ZAVG_PTS_ll (NLES_K,NLES_TIMES,NLES_MASKS))
  ALLOCATE(ZUND_PTS_ll (NLES_K,NLES_TIMES,NLES_MASKS))
  !
  ZAVG_PTS_ll(:,:,:) = NLES_AVG_PTS_ll(:,:,:)
  ZUND_PTS_ll(:,:,:) = NLES_UND_PTS_ll(:,:,:)
  ZCART_PTS_ll       = (NLESn_ISUP(IMI)-NLESn_IINF(IMI)+1) * (NLESn_JSUP(IMI)-NLESn_JINF(IMI)+1)
  !
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"AVG_PTS  ",YSUBTITLE(:), &
  "number of points used for averaging"//YSUBTITLE(:),"1",ZAVG_PTS_ll,HLES_AVG)
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"AVG_PTSF",YSUBTITLE(:), &
  "fraction of points used for averaging"//YSUBTITLE(:),"1",ZAVG_PTS_ll/ZCART_PTS_ll,HLES_AVG)
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"UND_PTS  ",YSUBTITLE(:), &
  "number of points below orography"//YSUBTITLE(:),"1",ZUND_PTS_ll,HLES_AVG)
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"UND_PTSF",YSUBTITLE(:), &
  "fraction of points below orography"//YSUBTITLE(:),"1",ZUND_PTS_ll/ZCART_PTS_ll,HLES_AVG)
  !
  DEALLOCATE(ZAVG_PTS_ll)
  DEALLOCATE(ZUND_PTS_ll)
END IF
!
!
!*      2.1  mean quantities
!            ---------------
!
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_U  ",YSUBTITLE(:), &
  "Mean U Profile"//YSUBTITLE(:),"m s-1",XLES_MEAN_U,HLES_AVG)

CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_V  ",YSUBTITLE(:), &
  "Mean V Profile"//YSUBTITLE(:),"m s-1",XLES_MEAN_V,HLES_AVG)

CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_W  ",YSUBTITLE(:), &
  "Mean W Profile"//YSUBTITLE(:),"m s-1",XLES_MEAN_W,HLES_AVG)

CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_PRE",YSUBTITLE(:), &
  "Mean pressure Profile"//YSUBTITLE(:),"Pa",XLES_MEAN_P,HLES_AVG)

CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_DP",YSUBTITLE(:), &
  "Mean Dyn production TKE Profile"//YSUBTITLE(:),"m2 s-3",XLES_MEAN_DP,HLES_AVG)

CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_TP",YSUBTITLE(:), &
  "Mean Thermal  production TKE Profile "//YSUBTITLE(:),"m2 s-3",XLES_MEAN_TP,HLES_AVG)

CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_TR",YSUBTITLE(:), &
  "Mean transport production TKE Profile"//YSUBTITLE(:),"m2 s-3",XLES_MEAN_TR,HLES_AVG)

CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_DISS",YSUBTITLE(:), &
  "Mean Dissipation TKE Profile"//YSUBTITLE(:),"m2 s-3",XLES_MEAN_DISS,HLES_AVG)

CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_LM",YSUBTITLE(:), &
  "Mean mixing length Profile"//YSUBTITLE(:),"m",XLES_MEAN_LM,HLES_AVG)

CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_RHO",YSUBTITLE(:), &
  "Mean density Profile"//YSUBTITLE(:),"kg m-3",XLES_MEAN_RHO,HLES_AVG)

CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_TH ",YSUBTITLE(:),&
  "Mean potential temperature Profile"//YSUBTITLE(:),"K",XLES_MEAN_Th,HLES_AVG)

CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_MF ",YSUBTITLE(:),&
  "Mass-flux Profile"//YSUBTITLE(:),"m s-1",XLES_MEAN_Mf,HLES_AVG)

IF (LUSERC) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_THL",YSUBTITLE(:), &
   "Mean liquid potential temperature Profile"//YSUBTITLE(:),"K",XLES_MEAN_Thl,HLES_AVG)

IF (LUSERV) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_THV",YSUBTITLE(:), &
   "Mean virtual potential temperature Profile"//YSUBTITLE(:),"K",XLES_MEAN_Thv,HLES_AVG)

IF (LUSERC) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_RT ",YSUBTITLE(:), &
  "Mean Rt Profile"//YSUBTITLE(:),"kg kg-1",XLES_MEAN_Rt,HLES_AVG)

IF (LUSERV) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_RV ",YSUBTITLE(:), &
  "Mean Rv Profile"//YSUBTITLE(:),"kg kg-1",XLES_MEAN_Rv,HLES_AVG)

IF (LUSERV) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_REHU ",YSUBTITLE(:), &
  "Mean Rh Profile"//YSUBTITLE(:),"percent",XLES_MEAN_Rehu,HLES_AVG)

IF (LUSERV) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_QS ",YSUBTITLE(:), &
  "Mean Qs Profile"//YSUBTITLE(:),"kg kg-1",XLES_MEAN_Qs,HLES_AVG)

IF (LUSERC) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_KHT ",YSUBTITLE(:),&
  "Eddy-diffusivity (temperature) Profile"//YSUBTITLE(:),"m2 s-1",XLES_MEAN_KHt,HLES_AVG)

IF (LUSERC) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_KHR ",YSUBTITLE(:),&
  "Eddy-diffusivity (wvapor) Profile"//YSUBTITLE(:),"m2 s-1",XLES_MEAN_KHr,HLES_AVG)

IF (LUSERC) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_RC ",YSUBTITLE(:), &
  "Mean Rc Profile"//YSUBTITLE(:),"kg kg-1",XLES_MEAN_Rc,HLES_AVG)

IF (LUSERC) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_CF ",YSUBTITLE(:), &
  "Mean Cf Profile"//YSUBTITLE(:),"1",XLES_MEAN_Cf,HLES_AVG)

IF (LUSERC) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_INDCF ",YSUBTITLE(:), &
  "Mean Cf>1-6 Profile (0 ou 1)"//YSUBTITLE(:),"1",XLES_MEAN_INDCf,HLES_AVG)

IF (LUSERC) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_INDCF2 ",YSUBTITLE(:), &
  "Mean Cf>1-5 Profile (0 ou 1)"//YSUBTITLE(:),"1",XLES_MEAN_INDCf2,HLES_AVG)

IF (LUSERR) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_RR ",YSUBTITLE(:), &
  "Mean Rr Profile"//YSUBTITLE(:),"kg kg-1",XLES_MEAN_Rr,HLES_AVG)

IF (LUSERR) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_RF ",YSUBTITLE(:), &
  "Mean RF Profile"//YSUBTITLE(:),"1",XLES_MEAN_RF,HLES_AVG)

IF (LUSERI) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_RI ",YSUBTITLE(:), &
  "Mean Ri Profile"//YSUBTITLE(:),"kg kg-1",XLES_MEAN_Ri,HLES_AVG)

IF (LUSERS) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_RS ",YSUBTITLE(:), &
  "Mean Rs Profile"//YSUBTITLE(:),"kg kg-1",XLES_MEAN_Rs,HLES_AVG)

IF (LUSERG) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_RG ",YSUBTITLE(:), &
  "Mean Rg Profile"//YSUBTITLE(:),"kg kg-1",XLES_MEAN_Rg,HLES_AVG)

IF (LUSERH) &
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEAN_RH ",YSUBTITLE(:), &
  "Mean Rh Profile"//YSUBTITLE(:),"kg kg-1",XLES_MEAN_Rh,HLES_AVG)

IF (NSV>0) &
CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"MEAN_SV ",YSUBTITLE(:), &
  "Mean Sv Profiles"//YSUBTITLE(:),"kg kg-1",XLES_MEAN_Sv,HLES_AVG)

CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEANWIND",YSUBTITLE(:), &
  "Profile of Mean Modulus of Wind"//YSUBTITLE(:),"m s-1",XLES_MEAN_WIND,HLES_AVG)
!
CALL LES_DIACHRO_MASKS(TPDIAFILE,"MEANMSFX",YSUBTITLE(:),  &
     "Total updraft mass flux"//YSUBTITLE(:),"kg m-2 s-1",XLES_RESOLVED_MASSFX   ,HLES_AVG)
!
IF (LLES_PDF) THEN
  CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"PDF_TH ",YSUBTITLE(:), &
  "Pdf potential temperature Profiles"//YSUBTITLE(:),"1",XLES_PDF_TH,HLES_AVG)
  CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"PDF_W ",YSUBTITLE(:), &
  "Pdf vertical velocity Profiles"//YSUBTITLE(:),"1",XLES_PDF_W,HLES_AVG)
  CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"PDF_THV ",YSUBTITLE(:), &
  "Pdf virtual pot. temp. Profiles"//YSUBTITLE(:),"1",XLES_PDF_THV,HLES_AVG)
    
  IF (LUSERV) THEN
   CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"PDF_RV ",YSUBTITLE(:), &
     "Pdf Rv Profiles"//YSUBTITLE(:),"1",XLES_PDF_RV,HLES_AVG)
  END IF

  IF (LUSERC) THEN
   CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"PDF_RC ",YSUBTITLE(:), &
   "Pdf Rc Profiles"//YSUBTITLE(:),"1",XLES_PDF_RC,HLES_AVG)

   CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"PDF_RT ",YSUBTITLE(:), &
   "Pdf Rt Profiles"//YSUBTITLE(:),"1",XLES_PDF_RT,HLES_AVG)

   CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"PDF_THL ",YSUBTITLE(:), &
   "Pdf Thl Profiles"//YSUBTITLE(:),"1",XLES_PDF_THL,HLES_AVG)
  END IF
  IF (LUSERR) &
  CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"PDF_RR ",YSUBTITLE(:), &
  "Pdf Rr Profiles"//YSUBTITLE(:),"1",XLES_PDF_RR,HLES_AVG)

  IF (LUSERI) &
  CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"PDF_RI ",YSUBTITLE(:), &
  "Pdf Ri Profiles"//YSUBTITLE(:),"1",XLES_PDF_RI,HLES_AVG)

  IF (LUSERS) &
  CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"PDF_RS ",YSUBTITLE(:), &
  "Pdf Rs Profiles"//YSUBTITLE(:),"1",XLES_PDF_RS,HLES_AVG)

  IF (LUSERG) &
  CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"PDF_RG ",YSUBTITLE(:), &
  "Pdf Rg Profiles"//YSUBTITLE(:),"1",XLES_PDF_RG,HLES_AVG)

END IF
!
!*      2.2  resolved quantities
!            -------------------
!
IF (LLES_RESOLVED) THEN
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_U2  ",YSUBTITLE(:), &
     "Resolved <u2> variance "//YSUBTITLE(:),"m2 s-2",XLES_RESOLVED_U2,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_V2  ",YSUBTITLE(:), &
     "Resolved <v2> variance"//YSUBTITLE(:),"m2 s-2",XLES_RESOLVED_V2,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_W2  ",YSUBTITLE(:), &
     "Resolved <w2> variance"//YSUBTITLE(:),"m2 s-2",XLES_RESOLVED_W2,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_UV  ",YSUBTITLE(:), &
     "Resolved <uv> Flux"//YSUBTITLE(:),"m2 s-2",XLES_RESOLVED_UV,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WU  ",YSUBTITLE(:), &
   "Resolved <wu> Flux"//YSUBTITLE(:),"m2 s-2",XLES_RESOLVED_WU,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WV  ",YSUBTITLE(:), &
     "Resolved <wv> Flux"//YSUBTITLE(:),"m2 s-2",XLES_RESOLVED_WV,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_KE  ",YSUBTITLE(:), &
     "Resolved TKE Profile"//YSUBTITLE(:),"m2 s-2",XLES_RESOLVED_Ke,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_P2  ",YSUBTITLE(:), &
     "Resolved pressure variance"//YSUBTITLE(:),"Pa2",XLES_RESOLVED_P2,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_UPZ ",YSUBTITLE(:), &
     "Resolved <up> horizontal Flux"//YSUBTITLE(:),"Pa s-1",XLES_RESOLVED_UP,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_VPZ ",YSUBTITLE(:), &
     "Resolved <vp> horizontal Flux"//YSUBTITLE(:),"Pa s-1",XLES_RESOLVED_VP,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WPZ ",YSUBTITLE(:), &
     "Resolved <wp> vertical Flux"//YSUBTITLE(:),"Pa s-1",XLES_RESOLVED_WP,HLES_AVG)

  IF (LUSERV) &
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_THTV ",YSUBTITLE(:), &
     "Resolved potential temperature - virtual potential temperature covariance"//YSUBTITLE(:), &
     "K2",XLES_RESOLVED_ThThv,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_TLTV ",YSUBTITLE(:), &
     "Resolved liquid potential temperature - virtual potential temperature covariance"//YSUBTITLE(:), &
     "K2",XLES_RESOLVED_ThlThv,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_TH2 ",YSUBTITLE(:), &
     "Resolved potential temperature variance"//YSUBTITLE(:),"K2",XLES_RESOLVED_Th2,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_THL2",YSUBTITLE(:), &
     "Resolved liquid potential temperature variance"//YSUBTITLE(:),"K2",XLES_RESOLVED_Thl2,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_UTH ",YSUBTITLE(:), &
     "Resolved <uth> horizontal Flux"//YSUBTITLE(:),"m K s-1",XLES_RESOLVED_UTh,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_VTH ",YSUBTITLE(:), &
     "Resolved <vth> horizontal Flux"//YSUBTITLE(:),"m K s-1",XLES_RESOLVED_VTh,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WTH ",YSUBTITLE(:), &
     "Resolved <wth> vertical Flux"//YSUBTITLE(:),"m K s-1",XLES_RESOLVED_WTh,HLES_AVG)

  IF (LUSERC) THEN
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_UTHL",YSUBTITLE(:), &
       "Resolved <uthl> horizontal Flux"//YSUBTITLE(:),"m K s-1",XLES_RESOLVED_UThl,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_VTHL",YSUBTITLE(:), &
       "Resolved <vthl> horizontal Flux"//YSUBTITLE(:),"m K s-1",XLES_RESOLVED_VThl,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WTHL",YSUBTITLE(:), &
       "Resolved <wthl> vertical Flux "//YSUBTITLE(:),"m K s-1",XLES_RESOLVED_WThl,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_RT2 ",YSUBTITLE(:), &
     "Resolved total water variance"//YSUBTITLE(:),"kg2 kg-2",XLES_RESOLVED_Rt2,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WRT ",YSUBTITLE(:), &
       "Resolved <wrt> vertical Flux "//YSUBTITLE(:),"m kg kg-1 s-1",XLES_RESOLVED_WRt,HLES_AVG)
  END IF

  IF (LUSERV) THEN

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_UTHV",YSUBTITLE(:), &
       "Resolved <uthv> horizontal Flux"//YSUBTITLE(:),"m K s-1",XLES_RESOLVED_UThv,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_VTHV",YSUBTITLE(:), &
       "Resolved <vthl> horizontal Flux"//YSUBTITLE(:),"m K s-1",XLES_RESOLVED_VThv,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WTHV",YSUBTITLE(:), &
       "Resolved <wthv> vertical Flux "//YSUBTITLE(:),"m K s-1",XLES_RESOLVED_WThv,HLES_AVG)
  END IF
!
  IF (LUSERV) THEN
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_RV2 ",YSUBTITLE(:), &
       "Resolved water vapor variance"//YSUBTITLE(:),"kg2 kg-2",XLES_RESOLVED_Rv2,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_THRV",YSUBTITLE(:), &
       "Resolved <thrv> covariance"//YSUBTITLE(:),"K kg kg-1",XLES_RESOLVED_ThRv,HLES_AVG)

    IF (LUSERC) &
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_TLRV",YSUBTITLE(:), &
       "Resolved <thlrv> covariance"//YSUBTITLE(:),"K kg kg-1",XLES_RESOLVED_ThlRv,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_TVRV",YSUBTITLE(:), &
       "Resolved <thvrv> covariance"//YSUBTITLE(:),"K kg kg-1",XLES_RESOLVED_ThvRv,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_URV ", YSUBTITLE(:), &
       "Resolved <urv> horizontal flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_RESOLVED_URv,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_VRV ", YSUBTITLE(:), &
       "Resolved <vrv> horizontal flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_RESOLVED_VRv,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WRV ", YSUBTITLE(:), &
       "Resolved <wrv> vertical flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_RESOLVED_WRv,HLES_AVG)
  END IF

  IF (LUSERC) THEN
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_RC2 ", YSUBTITLE(:), &
       "Resolved cloud water variance"//YSUBTITLE(:),"kg2 kg-2",XLES_RESOLVED_Rc2,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_THRC", YSUBTITLE(:), &
       "Resolved <thrc> covariance"//YSUBTITLE(:),"K kg kg-1",XLES_RESOLVED_ThRc,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_TLRC", YSUBTITLE(:), &
       "Resolved <thlrc> covariance"//YSUBTITLE(:),"K kg kg-1",XLES_RESOLVED_ThlRc,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_TVRC", YSUBTITLE(:), &
       "Resolved <thvrc> covariance"//YSUBTITLE(:),"K kg kg-1",XLES_RESOLVED_ThvRc,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_URC ", YSUBTITLE(:), &
       "Resolved <urc> horizontal flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_RESOLVED_URc,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_VRC ", YSUBTITLE(:), &
       "Resolved <vrc> horizontal flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_RESOLVED_VRc,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WRC ", YSUBTITLE(:), &
       "Resolved <wrc> vertical flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_RESOLVED_WRc,HLES_AVG)
  END IF

  IF (LUSERI) THEN
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_RI2 ", YSUBTITLE(:), &
       "Resolved cloud ice variance"//YSUBTITLE(:),"kg2 kg-2",XLES_RESOLVED_Ri2,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_THRI", YSUBTITLE(:), &
       "Resolved <thri> covariance"//YSUBTITLE(:),"K kg kg-1",XLES_RESOLVED_ThRi,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_TLRI", YSUBTITLE(:), &
       "Resolved <thlri> covariance"//YSUBTITLE(:),"K kg kg-1",XLES_RESOLVED_ThlRi,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_TVRI", YSUBTITLE(:), &
       "Resolved <thvri> covariance"//YSUBTITLE(:),"K kg kg-1",XLES_RESOLVED_ThvRi,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_URI ", YSUBTITLE(:), &
       "Resolved <uri> horizontal flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_RESOLVED_URi,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_VRI ", YSUBTITLE(:), &
       "Resolved <vri> horizontal flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_RESOLVED_VRi,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WRI ", YSUBTITLE(:), &
       "Resolved <wri> vertical flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_RESOLVED_WRi,HLES_AVG)
  END IF

  IF (LUSERR) THEN
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WRR ", YSUBTITLE(:), &
       "Resolved <wrr> vertical flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_RESOLVED_WRr,HLES_AVG)
    
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"INPRR3D ", YSUBTITLE(:), &
       "Precipitation flux"//YSUBTITLE(:),"m s-1",XLES_INPRR3D,HLES_AVG)
        
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"MAXINPR3D ", YSUBTITLE(:), &
       "Max Precip flux"//YSUBTITLE(:),"m s-1",XLES_MAX_INPRR3D,HLES_AVG)
        
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"EVAP3D ", YSUBTITLE(:), &
       "Evaporation profile"//YSUBTITLE(:),"kg kg-1 s-1",XLES_EVAP3D,HLES_AVG)
  ENDIF
  IF (NSV>0) THEN
    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"RES_SV2 ", YSUBTITLE(:), &
       "Resolved scalar variables variances"//YSUBTITLE(:),"kg2 kg-2",XLES_RESOLVED_Sv2,HLES_AVG)

    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"RES_THSV", YSUBTITLE(:), &
       "Resolved <ThSv> variance"//YSUBTITLE(:),"K kg kg-1",XLES_RESOLVED_ThSv,HLES_AVG)

    IF (LUSERC) &
    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"RES_TLSV", YSUBTITLE(:), &
       "Resolved <ThlSv> variance"//YSUBTITLE(:),"K kg kg-1",XLES_RESOLVED_ThlSv,HLES_AVG)

    IF (LUSERV) &
    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"RES_TVSV", YSUBTITLE(:), &
       "Resolved <ThvSv> variance"//YSUBTITLE(:),"K kg kg-1",XLES_RESOLVED_ThvSv,HLES_AVG)

    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"RES_USV ", YSUBTITLE(:), &
       "Resolved <uSv> horizontal flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_RESOLVED_USv,HLES_AVG)

    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"RES_VSV ", YSUBTITLE(:), &
       "Resolved <vSv> horizontal flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_RESOLVED_VSv,HLES_AVG)

    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"RES_WSV ", YSUBTITLE(:), &
       "Resolved <wSv> vertical flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_RESOLVED_WSv,HLES_AVG)
  END IF

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_U3  ",YSUBTITLE(:),  &
       "Resolved <w3>"//YSUBTITLE(:),"m3 s-3",XLES_RESOLVED_U3,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_V3  ",YSUBTITLE(:),  &
       "Resolved <w3>"//YSUBTITLE(:),"m3 s-3",XLES_RESOLVED_V3,HLES_AVG)
    
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_W3  ",YSUBTITLE(:),  &
       "Resolved <w3>"//YSUBTITLE(:),"m3 s-3",XLES_RESOLVED_W3,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_U4  ",YSUBTITLE(:),  &
       "Resolved <w3>"//YSUBTITLE(:),"m4 s-4",XLES_RESOLVED_U4,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_V4  ",YSUBTITLE(:),  &
       "Resolved <w3>"//YSUBTITLE(:),"m4 s-4",XLES_RESOLVED_V4,HLES_AVG)
    
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_W4  ",YSUBTITLE(:),  &
       "Resolved <w3>"//YSUBTITLE(:),"m4 s-4",XLES_RESOLVED_W4,HLES_AVG)


  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WTL2",YSUBTITLE(:),  &
       "Resolved <wThl2>"//YSUBTITLE(:),"m K2 s-1",XLES_RESOLVED_WThl2,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_W2TL",YSUBTITLE(:),  &
       "Resolved <w2Thl>"//YSUBTITLE(:),"m2 K s-2",XLES_RESOLVED_W2Thl,HLES_AVG)

  IF (LUSERV) THEN
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WRV2",YSUBTITLE(:),  &
         "Resolved <wRv2>"//YSUBTITLE(:),"m kg2 kg-2 s-1",XLES_RESOLVED_WRv2,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_W2RV",YSUBTITLE(:),  &
         "Resolved <w2Rv>"//YSUBTITLE(:),"m2 kg kg-1 s-2",XLES_RESOLVED_W2Rv,HLES_AVG)
     
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WRT2",YSUBTITLE(:),  &
         "Resolved <wRt2>"//YSUBTITLE(:),"m kg2 kg-2 s-1",XLES_RESOLVED_WRt2,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_W2RT",YSUBTITLE(:),  &
         "Resolved <w2Rt>"//YSUBTITLE(:),"m2 kg kg-1 s-2",XLES_RESOLVED_W2Rt,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RE_WTLRV",YSUBTITLE(:),  &
         "Resolved <wThlRv>"//YSUBTITLE(:),"m K kg kg-1 s-1",XLES_RESOLVED_WThlRv,HLES_AVG)
   
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RE_WTLRT",YSUBTITLE(:),  &
         "Resolved <wThlRt>"//YSUBTITLE(:),"m K kg kg-1 s-1",XLES_RESOLVED_WThlRt,HLES_AVG)
  END IF

  IF (LUSERC) THEN
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WRC2",YSUBTITLE(:),  &
         "Resolved <wRc2>"//YSUBTITLE(:),"m kg2 kg-2 s-1",XLES_RESOLVED_WRc2,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_W2RC",YSUBTITLE(:),  &
         "Resolved <w2Rc>"//YSUBTITLE(:),"m2 kg kg-1 s-2",XLES_RESOLVED_W2Rc,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RE_WTLRC",YSUBTITLE(:),  &
         "Resolved <wThlRc>"//YSUBTITLE(:),"m K kg kg-1 s-1",XLES_RESOLVED_WThlRc,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RE_WRVRC",YSUBTITLE(:),  &
         "Resolved <wRvRc>"//YSUBTITLE(:),"m kg2 kg-2 s-1",XLES_RESOLVED_WRvRc,HLES_AVG)
  END IF

  IF (LUSERI) THEN
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WRI2",YSUBTITLE(:),  &
         "Resolved <wRi2>"//YSUBTITLE(:),"m kg2 kg-2 s-1",XLES_RESOLVED_WRi2,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_W2RI",YSUBTITLE(:),  &
         "Resolved <w2Ri>"//YSUBTITLE(:),"m2 kg kg-1 s-2",XLES_RESOLVED_W2Ri,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RE_WTLRI",YSUBTITLE(:),  &
         "Resolved <wThlRi>"//YSUBTITLE(:),"m K kg kg-1 s-1",XLES_RESOLVED_WThlRi,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"RE_WRVRI",YSUBTITLE(:),  &
         "Resolved <wRvRi>"//YSUBTITLE(:),"m kg2 kg-2 s-1",XLES_RESOLVED_WRvRi,HLES_AVG)
  END IF

  IF (NSV>0) THEN
    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"RES_WSV2",YSUBTITLE(:),  &
         "Resolved <wSv2>"//YSUBTITLE(:),"m kg2 kg-2 s-1",XLES_RESOLVED_WSv2,HLES_AVG)

    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"RES_W2SV",YSUBTITLE(:),  &
         "Resolved <w2Sv>"//YSUBTITLE(:),"m2 kg kg-1 s-2",XLES_RESOLVED_W2Sv,HLES_AVG)

    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"RE_WTLSV",YSUBTITLE(:),  &
         "Resolved <wThlSv>"//YSUBTITLE(:),"m K kg kg-1 s-1",XLES_RESOLVED_WThlSv,HLES_AVG)

    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"RE_WRVSV",YSUBTITLE(:),  &
         "Resolved <wRvSv>"//YSUBTITLE(:),"m kg2 kg-2 s-1",XLES_RESOLVED_WRvSv,HLES_AVG)
  END IF

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_TLPZ",YSUBTITLE(:),  &
       "Resolved <Thldp/dz>"//YSUBTITLE(:),"K Pa m-1",XLES_RESOLVED_ThlPz,HLES_AVG)

  IF (LUSERV) &
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_RVPZ",YSUBTITLE(:),  &
       "Resolved <Rvdp/dz>"//YSUBTITLE(:),"kg2 kg-2 Pa m-1",XLES_RESOLVED_RvPz,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_RCPZ",YSUBTITLE(:),  &
       "Resolved <Rcdp/dz>"//YSUBTITLE(:),"kg2 kg-2 Pa m-1",XLES_RESOLVED_RcPz,HLES_AVG)

  IF (LUSERI) &
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_RIPZ",YSUBTITLE(:),  &
       "Resolved <Ridp/dz>"//YSUBTITLE(:),"kg2 kg-2 Pa m-1",XLES_RESOLVED_RiPz,HLES_AVG)

  IF (NSV>0) THEN
    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"RES_SVPZ",YSUBTITLE(:),  &
         "Resolved <Svdp/dz>"//YSUBTITLE(:),"kg2 kg-2 Pa m-1",XLES_RESOLVED_SvPz,HLES_AVG)
  END IF

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_UKE ", YSUBTITLE(:), &
       "Resolved flux of resolved kinetic energy"//YSUBTITLE(:),"m3 s-3",XLES_RESOLVED_UKe,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_VKE ", YSUBTITLE(:), &
       "Resolved flux of resolved kinetic energy"//YSUBTITLE(:),"m3 s-3",XLES_RESOLVED_VKe,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RES_WKE ", YSUBTITLE(:), &
       "Resolved flux of resolved kinetic energy"//YSUBTITLE(:),"m3 s-3",XLES_RESOLVED_WKe,HLES_AVG)

END IF
!
!
!*      2.3  subgrid quantities
!            ------------------
!
IF (LLES_SUBGRID) THEN

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_TKE ",YSUBTITLE(:), &
       "Subgrid TKE"//YSUBTITLE(:),"m2 s-2",XLES_SUBGRID_Tke,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_U2  ", YSUBTITLE(:), &
       "Subgrid <u2> variance"//YSUBTITLE(:),"m2 s-2",XLES_SUBGRID_U2,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_V2  ", YSUBTITLE(:), &
       "Subgrid <v2> variance"//YSUBTITLE(:),"m2 s-2",XLES_SUBGRID_V2,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_W2  ", YSUBTITLE(:), &
       "Subgrid <w2> variance"//YSUBTITLE(:),"m2 s-2",XLES_SUBGRID_W2,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_UV  ", YSUBTITLE(:), &
       "Subgrid <uv> flux"//YSUBTITLE(:),"m2 s-2",XLES_SUBGRID_UV,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_WU  ", YSUBTITLE(:), &
       "Subgrid <wu> flux"//YSUBTITLE(:),"m2 s-2",XLES_SUBGRID_WU,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_WV  ", YSUBTITLE(:), &
       "Subgrid <wv> flux"//YSUBTITLE(:),"m2 s-2",XLES_SUBGRID_WV,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_THL2", YSUBTITLE(:), &
       "Subgrid liquid potential temperature variance"//YSUBTITLE(:),"K2",XLES_SUBGRID_Thl2,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_UTHL", YSUBTITLE(:), &
       "Subgrid hor. flux of liquid potential temperature"//YSUBTITLE(:),"m K s-1",XLES_SUBGRID_UThl,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_VTHL", YSUBTITLE(:), &
       "Subgrid hor. flux of liquid potential temperature"//YSUBTITLE(:),"m K s-1",XLES_SUBGRID_VThl,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_WTHL", YSUBTITLE(:), &
       "Subgrid vert. flux of liquid potential temperature"//YSUBTITLE(:),"m K s-1",XLES_SUBGRID_WThl,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_WP  ",YSUBTITLE(:), &
     "Subgrid <wp> vertical Flux"//YSUBTITLE(:),"mPa/s",XLES_SUBGRID_WP,HLES_AVG)
!!
!!
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"THLUP_MF",YSUBTITLE(:), &
     "Subgrid <thl> of updraft"//YSUBTITLE(:),"K",XLES_SUBGRID_THLUP_MF,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RTUP_MF ",YSUBTITLE(:), &
     "Subgrid <rt> of updraft"//YSUBTITLE(:),"kg kg-1",XLES_SUBGRID_RTUP_MF,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RVUP_MF ",YSUBTITLE(:), &
     "Subgrid <rv> of updraft"//YSUBTITLE(:),"kg kg-1",XLES_SUBGRID_RVUP_MF,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RCUP_MF ",YSUBTITLE(:), &
     "Subgrid <rc> of updraft"//YSUBTITLE(:),"kg kg-1",XLES_SUBGRID_RCUP_MF,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"RIUP_MF ",YSUBTITLE(:), &
     "Subgrid <ri> of updraft"//YSUBTITLE(:),"kg kg-1",XLES_SUBGRID_RIUP_MF,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"WUP_MF  ",YSUBTITLE(:), &
     "Subgrid <w> of updraft"//YSUBTITLE(:),"m s-1",XLES_SUBGRID_WUP_MF,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"MAFLX_MF",YSUBTITLE(:), &
     "Subgrid <MF> of updraft"//YSUBTITLE(:),"kg m-2 s-1",XLES_SUBGRID_MASSFLUX,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"DETR_MF ",YSUBTITLE(:), &
     "Subgrid <detr> of updraft"//YSUBTITLE(:),"kg m-3 s-1",XLES_SUBGRID_DETR,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"ENTR_MF ",YSUBTITLE(:), &
     "Subgrid <entr> of updraft"//YSUBTITLE(:),"kg m-3 s-1",XLES_SUBGRID_ENTR,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"FRCUP_MF",YSUBTITLE(:), &
     "Subgrid <FracUp> of updraft"//YSUBTITLE(:),"1",XLES_SUBGRID_FRACUP,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"THVUP_MF",YSUBTITLE(:), &
     "Subgrid <thv> of updraft"//YSUBTITLE(:),"K",&
                                XLES_SUBGRID_THVUP_MF,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"WTHL_MF ",YSUBTITLE(:), &
     "Subgrid <wthl> of mass flux convection scheme"//YSUBTITLE(:),"m K s-1",&
                                XLES_SUBGRID_WTHLMF,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"WRT_MF  ",YSUBTITLE(:), &
     "Subgrid <wrt> of mass flux convection scheme"//YSUBTITLE(:),"m kg kg-1 s-1",&
                                XLES_SUBGRID_WRTMF,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"WTHV_MF ",YSUBTITLE(:), &
     "Subgrid <wthv> of mass flux convection scheme"//YSUBTITLE(:),"m K s-1",&
                                XLES_SUBGRID_WTHVMF,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"WU_MF   ",YSUBTITLE(:), &
     "Subgrid <wu> of mass flux convection scheme"//YSUBTITLE(:),"m2 s-2",&
                                XLES_SUBGRID_WUMF,HLES_AVG)
!     
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"WV_MF   ",YSUBTITLE(:), &
     "Subgrid <wv> of mass flux convection scheme"//YSUBTITLE(:),"m2 s-2",&
                                XLES_SUBGRID_WVMF,HLES_AVG)
!!     

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_PHI3",YSUBTITLE(:), &
     "Subgrid Phi3 function"//YSUBTITLE(:),"1",XLES_SUBGRID_PHI3,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_LMIX",YSUBTITLE(:), &
     "Subgrid Mixing Length"//YSUBTITLE(:),"1",XLES_SUBGRID_LMix,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_LDIS",YSUBTITLE(:), &
     "Subgrid Dissipation Length"//YSUBTITLE(:),"1",XLES_SUBGRID_LDiss,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_KM  ",YSUBTITLE(:), &
     "Eddy diffusivity for momentum"//YSUBTITLE(:),"m2 s-1",XLES_SUBGRID_Km,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_KH  ",YSUBTITLE(:), &
     "Eddy diffusivity for heat"//YSUBTITLE(:),"m2 s-1",XLES_SUBGRID_Kh,HLES_AVG)
!
  IF (LUSERV) THEN
     CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_WTHV", YSUBTITLE(:), &
       "Subgrid vert. flux of liquid potential temperature"//YSUBTITLE(:),"m K s-1",XLES_SUBGRID_WThv,HLES_AVG)
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_RT2 ", YSUBTITLE(:), &
       "Subgrid total water variance"//YSUBTITLE(:),"kg2 kg-2",XLES_SUBGRID_Rt2,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_TLRT", YSUBTITLE(:), &
       "Subgrid <thlrt> covariance"//YSUBTITLE(:),"K kg kg-1",XLES_SUBGRID_ThlRt,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_URT ", YSUBTITLE(:), &
       "Subgrid total water horizontal flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_SUBGRID_URt,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_VRT ", YSUBTITLE(:), &
       "Subgrid total water horizontal flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_SUBGRID_VRt,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_WRT ", YSUBTITLE(:), &
       "Subgrid total water vertical flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_SUBGRID_WRt,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_PSI3",YSUBTITLE(:), &
       "Subgrid Psi3 function"//YSUBTITLE(:),"1",XLES_SUBGRID_PSI3,HLES_AVG)
  END IF

  IF (LUSERC) THEN
    CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_RC2 ", YSUBTITLE(:), &
       "Subgrid cloud water variance"//YSUBTITLE(:),"kg2 kg-2",XLES_SUBGRID_Rc2,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_URC ", YSUBTITLE(:), &
       "Subgrid cloud water horizontal flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_SUBGRID_URc,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_VRC ", YSUBTITLE(:), &
       "Subgrid cloud water horizontal flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_SUBGRID_VRc,HLES_AVG)

    CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_WRC ", YSUBTITLE(:), &
       "Subgrid cloud water vertical flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_SUBGRID_WRc,HLES_AVG)
  END IF

  IF (NSV>0) THEN
    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"SBG_USV ", YSUBTITLE(:), &
       "Subgrid <uSv> horizontal flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_SUBGRID_USv,HLES_AVG)

    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"SBG_VSV ", YSUBTITLE(:), &
       "Subgrid <vSv> horizontal flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_SUBGRID_VSv,HLES_AVG)

    CALL LES_DIACHRO_SV_MASKS(TPDIAFILE,"SBG_WSV ", YSUBTITLE(:), &
       "Subgrid <wSv> vertical flux"//YSUBTITLE(:),"m kg kg-1 s-1",XLES_SUBGRID_WSv,HLES_AVG)
  END IF

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_UTKE", YSUBTITLE(:), &
     "Subgrid flux of subgrid kinetic energy"//YSUBTITLE(:),"m3 s-3",XLES_SUBGRID_UTke,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_VTKE", YSUBTITLE(:), &
     "Subgrid flux of subgrid kinetic energy"//YSUBTITLE(:),"m3 s-3",XLES_SUBGRID_VTke,HLES_AVG)

  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_WTKE",YSUBTITLE(:),  &
     "Subgrid flux of subgrid kinetic energy"//YSUBTITLE(:),"m3 s-3",XLES_SUBGRID_WTke,HLES_AVG)
!
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_W2TL",YSUBTITLE(:),  &
     "Subgrid flux of subgrid kinetic energy"//YSUBTITLE(:),"m2 K s-2",XLES_SUBGRID_W2Thl,HLES_AVG)
!
  CALL LES_DIACHRO_MASKS(TPDIAFILE,"SBG_WTL2",YSUBTITLE(:),  &
     "Subgrid flux of subgrid kinetic energy"//YSUBTITLE(:),"m K2 s-1",XLES_SUBGRID_WThl2,HLES_AVG)
!
END IF
!
!*      2.4  Updraft quantities
!            ------------------
!
IF (LLES_UPDRAFT) THEN
  CALL LES_DIACHRO(TPDIAFILE,"UP_FRAC ",  &
       "Updraft fraction","1",XLES_UPDRAFT,HLES_AVG)

  CALL LES_DIACHRO(TPDIAFILE,"UP_W    ",  &
       "Updraft W mean value","m s-1",XLES_UPDRAFT_W,HLES_AVG)

  CALL LES_DIACHRO(TPDIAFILE,"UP_TH   ",  &
       "Updraft potential temperature mean value","K",XLES_UPDRAFT_Th,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO(TPDIAFILE,"UP_THL  ",  &
       "Updraft liquid potential temperature mean value","K",XLES_UPDRAFT_Thl,HLES_AVG)

  IF (LUSERV) &
  CALL LES_DIACHRO(TPDIAFILE,"UP_THV  ",  &
       "Updraft virutal potential temperature mean value","K",XLES_UPDRAFT_Thv,HLES_AVG)

  CALL LES_DIACHRO(TPDIAFILE,"UP_KE   ",  &
       "Updraft resolved TKE mean value","m2 s-2",XLES_UPDRAFT_Ke,HLES_AVG)

  CALL LES_DIACHRO(TPDIAFILE,"UP_TKE  ",  &
       "Updraft subgrid TKE mean value","m2 s-2",XLES_UPDRAFT_Tke,HLES_AVG)

  IF (LUSERV) &
  CALL LES_DIACHRO(TPDIAFILE,"UP_RV   ",  &
       "Updraft water vapor mean value","kg kg-1",XLES_UPDRAFT_Rv,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO(TPDIAFILE,"UP_RC   ",  &
       "Updraft cloud water mean value","kg kg-1",XLES_UPDRAFT_Rc,HLES_AVG)

  IF (LUSERR) &
  CALL LES_DIACHRO(TPDIAFILE,"UP_RR   ",  &
       "Updraft rain mean value","kg kg-1",XLES_UPDRAFT_Rr,HLES_AVG)

  IF (LUSERI) &
  CALL LES_DIACHRO(TPDIAFILE,"UP_RI   ",  &
       "Updraft ice mean value","kg kg-1",XLES_UPDRAFT_Ri,HLES_AVG)

  IF (LUSERS) &
  CALL LES_DIACHRO(TPDIAFILE,"UP_RS   ",  &
       "Updraft snow mean value","kg kg-1",XLES_UPDRAFT_Rs,HLES_AVG)

  IF (LUSERG) &
  CALL LES_DIACHRO(TPDIAFILE,"UP_RG   ",  &
       "Updraft graupel mean value","kg kg-1",XLES_UPDRAFT_Rg,HLES_AVG)

  IF (LUSERH) &
  CALL LES_DIACHRO(TPDIAFILE,"UP_RH   ",  &
       "Updraft hail mean value","kg kg-1",XLES_UPDRAFT_Rh,HLES_AVG)

  IF (NSV>0) &
  CALL LES_DIACHRO_SV(TPDIAFILE,"UP_SV   ",  &
       "Updraft scalar variables mean values","kg kg-1",XLES_UPDRAFT_Sv,HLES_AVG)
  !
  CALL LES_DIACHRO(TPDIAFILE,"UP_TH2 ",  &
       "Updraft resolved Theta variance ","K2",XLES_UPDRAFT_Th2,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO(TPDIAFILE,"UP_THL2",  &
       "Updraft resolved Theta_l variance ","K2",XLES_UPDRAFT_Thl2,HLES_AVG)

  IF (LUSERV) &
  CALL LES_DIACHRO(TPDIAFILE,"UP_THTV",  &
       "Updraft resolved Theta Theta_v covariance ","K2",XLES_UPDRAFT_ThThv,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO(TPDIAFILE,"UP_TLTV",  &
       "Updraft resolved Theta_l Theta_v covariance ","K2",XLES_UPDRAFT_ThlThv,HLES_AVG)

  CALL LES_DIACHRO(TPDIAFILE,"UP_WTH  ",  &
       "Updraft resolved WTh flux","m K s-1",XLES_UPDRAFT_WTh,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO(TPDIAFILE,"UP_WTHL ",  &
       "Updraft resolved WThl flux","m K s-1",XLES_UPDRAFT_WThl,HLES_AVG)

  IF (LUSERV) &
  CALL LES_DIACHRO(TPDIAFILE,"UP_WTHV ",  &
       "Updraft resolved WThv flux","m K s-1",XLES_UPDRAFT_WThv,HLES_AVG)
  !
  IF (LUSERV) THEN
    CALL LES_DIACHRO(TPDIAFILE,"UP_RV2  ",  &
       "Updraft resolved water vapor variance","kg2 kg-2",XLES_UPDRAFT_Rv2,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"UP_THRV ",  &
       "Updraft resolved <thrv> covariance","K kg kg-1",XLES_UPDRAFT_ThRv,HLES_AVG)

    IF (LUSERC) &
    CALL LES_DIACHRO(TPDIAFILE,"UP_THLRV",  &
     "Updraft resolved <thlrv> covariance","K kg kg-1",XLES_UPDRAFT_ThlRv,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"UP_THVRV",  &
       "Updraft resolved <thvrv> covariance","K kg kg-1",XLES_UPDRAFT_ThvRv,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"UP_WRV  ",  &
       "Updraft resolved <wrv> vertical flux","m kg kg-1 s-1",XLES_UPDRAFT_WRv,HLES_AVG)
  END IF

  IF (LUSERC) THEN
    CALL LES_DIACHRO(TPDIAFILE,"UP_RC2  ",  &
       "Updraft resolved cloud water variance","kg2 kg-2",XLES_UPDRAFT_Rc2,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"UP_THRC ",  &
       "Updraft resolved <thrc> covariance","K kg kg-1",XLES_UPDRAFT_ThRc,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"UP_THLRC", &
       "Updraft resolved <thlrc> covariance","K kg kg-1",XLES_UPDRAFT_ThlRc,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"UP_THVRC", &
       "Updraft resolved <thvrc> covariance","K kg kg-1",XLES_UPDRAFT_ThvRc,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"UP_WRC  ",  &
       "Updraft resolved <wrc> vertical flux","m kg kg-1 s-1",XLES_UPDRAFT_WRc,HLES_AVG)
  END IF

  IF (LUSERI) THEN
    CALL LES_DIACHRO(TPDIAFILE,"UP_RI2  ", &
       "Updraft resolved cloud ice variance","kg2 kg-2",XLES_UPDRAFT_Ri2,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"UP_THRI ",  &
       "Updraft resolved <thri> covariance","K kg kg-1",XLES_UPDRAFT_ThRi,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"UP_THLRI",  &
       "Updraft resolved <thlri> covariance","K kg kg-1",XLES_UPDRAFT_ThlRi,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"UP_THVRI",  &
       "Updraft resolved <thvri> covariance","K kg kg-1",XLES_UPDRAFT_ThvRi,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"UP_WRI  ",  &
       "Updraft resolved <wri> vertical flux","m kg kg-1 s-1",XLES_UPDRAFT_WRi,HLES_AVG)
  END IF

  IF (NSV>0) THEN
    CALL LES_DIACHRO_SV(TPDIAFILE,"UP_SV2  ",  &
       "Updraft resolved scalar variables variances","kg2 kg-2",XLES_UPDRAFT_Sv2,HLES_AVG)

    CALL LES_DIACHRO_SV(TPDIAFILE,"UP_THSV ", &
       "Updraft resolved <ThSv> variance","K kg kg-1",XLES_UPDRAFT_ThSv,HLES_AVG)

    IF (LUSERC) &
    CALL LES_DIACHRO_SV(TPDIAFILE,"UP_THLSV",  &
       "Updraft resolved <ThlSv> variance","K kg kg-1",XLES_UPDRAFT_ThlSv,HLES_AVG)

    IF (LUSERV) &
    CALL LES_DIACHRO_SV(TPDIAFILE,"UP_THVSV",  &
       "Updraft resolved <ThvSv> variance","K kg kg-1",XLES_UPDRAFT_ThvSv,HLES_AVG)

    CALL LES_DIACHRO_SV(TPDIAFILE,"UP_WSV  ",  &
       "Updraft resolved <wSv> vertical flux","m kg kg-1 s-1",XLES_UPDRAFT_WSv,HLES_AVG)
  END IF
END IF
!                
!
!*      2.5  Downdraft quantities
!            --------------------
!
IF (LLES_DOWNDRAFT) THEN
   CALL LES_DIACHRO(TPDIAFILE,"DW_FRAC ",  &
       "Downdraft fraction","1",XLES_DOWNDRAFT,HLES_AVG)

  CALL LES_DIACHRO(TPDIAFILE,"DW_W    ", &
       "Downdraft W mean value","m s-1",XLES_DOWNDRAFT_W,HLES_AVG)

  CALL LES_DIACHRO(TPDIAFILE,"DW_TH   ",  &
       "Downdraft potential temperature mean value","K",XLES_DOWNDRAFT_Th,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO(TPDIAFILE,"DW_THL  ", &
       "Downdraft liquid potential temperature mean value","K",XLES_DOWNDRAFT_Thl,HLES_AVG)

  IF (LUSERV) &
  CALL LES_DIACHRO(TPDIAFILE,"DW_THV  ",  &
       "Downdraft virtual potential temperature mean value","K",XLES_DOWNDRAFT_Thv,HLES_AVG)


  CALL LES_DIACHRO(TPDIAFILE,"DW_KE   ", &
       "Downdraft resolved TKE mean value","m2 s-2",XLES_DOWNDRAFT_Ke,HLES_AVG)

  CALL LES_DIACHRO(TPDIAFILE,"DW_TKE  ",  &
       "Downdraft subgrid TKE mean value","m2 s-2",XLES_DOWNDRAFT_Tke,HLES_AVG)

  IF (LUSERV) &
  CALL LES_DIACHRO(TPDIAFILE,"DW_RV   ",  &
       "Downdraft water vapor mean value","kg kg-1",XLES_DOWNDRAFT_Rv,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO(TPDIAFILE,"DW_RC   ",  &
       "Downdraft cloud water mean value","kg kg-1",XLES_DOWNDRAFT_Rc,HLES_AVG)

  IF (LUSERR) &
  CALL LES_DIACHRO(TPDIAFILE,"DW_RR   ",  &
       "Downdraft rain mean value","kg kg-1",XLES_DOWNDRAFT_Rr,HLES_AVG)

  IF (LUSERI) &
  CALL LES_DIACHRO(TPDIAFILE,"DW_RI   ",  &
       "Downdraft ice mean value","kg kg-1",XLES_DOWNDRAFT_Ri,HLES_AVG)

  IF (LUSERS) &
  CALL LES_DIACHRO(TPDIAFILE,"DW_RS   ",  &
       "Downdraft snow mean value","kg kg-1",XLES_DOWNDRAFT_Rs,HLES_AVG)

  IF (LUSERG) &
  CALL LES_DIACHRO(TPDIAFILE,"DW_RG   ",  &
       "Downdraft graupel mean value","kg kg-1",XLES_DOWNDRAFT_Rg,HLES_AVG)

  IF (LUSERH) &
  CALL LES_DIACHRO(TPDIAFILE,"DW_RH   ",  &
       "Downdraft hail mean value","kg kg-1",XLES_DOWNDRAFT_Rh,HLES_AVG)

  IF (NSV>0) &
  CALL LES_DIACHRO_SV(TPDIAFILE,"DW_SV   ", &
       "Downdraft scalar variables mean values","kg kg-1",XLES_DOWNDRAFT_Sv,HLES_AVG)
  !
  CALL LES_DIACHRO(TPDIAFILE,"DW_TH2 ",  &
       "Downdraft resolved Theta variance ","K2",XLES_DOWNDRAFT_Th2,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO(TPDIAFILE,"DW_THL2",  &
       "Downdraft resolved Theta_l variance ","K2",XLES_DOWNDRAFT_Thl2,HLES_AVG)

  IF (LUSERV) &
  CALL LES_DIACHRO(TPDIAFILE,"DW_THTV ",  &
       "Downdraft resolved Theta Theta_v covariance ","K2",XLES_DOWNDRAFT_ThThv,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO(TPDIAFILE,"DW_TLTV ",  &
       "Downdraft resolved Theta_l Theta_v covariance ","K2",XLES_DOWNDRAFT_ThlThv,HLES_AVG)

  CALL LES_DIACHRO(TPDIAFILE,"DW_WTH  ",  &
       "Downdraft resolved WTh flux","m K s-1",XLES_DOWNDRAFT_WTh,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO(TPDIAFILE,"DW_WTHL ",  &
       "Downdraft resolved WThl flux","m K s-1",XLES_DOWNDRAFT_WThl,HLES_AVG)

  IF (LUSERV) &
  CALL LES_DIACHRO(TPDIAFILE,"DW_WTHV ",  &
       "Downdraft resolved WThv flux","m K s-1",XLES_DOWNDRAFT_WThv,HLES_AVG)
  !
  IF (LUSERV) THEN
    CALL LES_DIACHRO(TPDIAFILE,"DW_RV2  ",  &
       "Downdraft resolved water vapor variance","kg2 kg-2",XLES_DOWNDRAFT_Rv2,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"DW_THRV ",  &
       "Downdraft resolved <thrv> covariance","K kg kg-1",XLES_DOWNDRAFT_ThRv,HLES_AVG)

    IF (LUSERC) &
    CALL LES_DIACHRO(TPDIAFILE,"DW_THLRV",  &
       "Downdraft resolved <thlrv> covariance","K kg kg-1",XLES_DOWNDRAFT_ThlRv,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"DW_THVRV",  &
       "Downdraft resolved <thvrv> covariance","K kg kg-1",XLES_DOWNDRAFT_ThvRv,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"DW_WRV  ",  &
       "Downdraft resolved <wrv> vertical flux","m kg kg-1 s-1",XLES_DOWNDRAFT_WRv,HLES_AVG)
  END IF

  IF (LUSERC) THEN
    CALL LES_DIACHRO(TPDIAFILE,"DW_RC2  ",  &
       "Downdraft resolved cloud water variance","kg2 kg-2",XLES_DOWNDRAFT_Rc2,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"DW_THRC ",  &
       "Downdraft resolved <thrc> covariance","K kg kg-1",XLES_DOWNDRAFT_ThRc,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"DW_THLRC",  &
       "Downdraft resolved <thlrc> covariance","K kg kg-1",XLES_DOWNDRAFT_ThlRc,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"DW_THVRC",  &
       "Downdraft resolved <thvrc> covariance","K kg kg-1",XLES_DOWNDRAFT_ThvRc,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"DW_WRC  ",  &
       "Downdraft resolved <wrc> vertical flux","m kg kg-1 s-1",XLES_DOWNDRAFT_WRc,HLES_AVG)
  END IF

  IF (LUSERI) THEN
    CALL LES_DIACHRO(TPDIAFILE,"DW_RI2  ",  &
       "Downdraft resolved cloud ice variance","kg2 kg-2",XLES_DOWNDRAFT_Ri2,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"DW_THRI ",  &
       "Downdraft resolved <thri> covariance","K kg kg-1",XLES_DOWNDRAFT_ThRi,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"DW_THLRI", &
       "Downdraft resolved <thlri> covariance","K kg kg-1",XLES_DOWNDRAFT_ThlRi,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"DW_THVRI",  &
       "Downdraft resolved <thvri> covariance","K kg kg-1",XLES_DOWNDRAFT_ThvRi,HLES_AVG)

    CALL LES_DIACHRO(TPDIAFILE,"DW_WRI  ", &
       "Downdraft resolved <wri> vertical flux","m kg kg-1 s-1",XLES_DOWNDRAFT_WRi,HLES_AVG)
  END IF

  IF (NSV>0) THEN
    CALL LES_DIACHRO_SV(TPDIAFILE,"DW_SV2  ", &
       "Downdraft resolved scalar variables variances","kg2 kg-2",XLES_DOWNDRAFT_Sv2,HLES_AVG)

   CALL LES_DIACHRO_SV(TPDIAFILE,"DW_THSV ",  &
       "Downdraft resolved <ThSv> variance","K kg kg-1",XLES_DOWNDRAFT_ThSv,HLES_AVG)

    IF (LUSERC) &
    CALL LES_DIACHRO_SV(TPDIAFILE,"DW_THLSV",  &
       "Downdraft resolved <ThlSv> variance","K kg kg-1",XLES_DOWNDRAFT_ThlSv,HLES_AVG)

    IF (LUSERV) &
    CALL LES_DIACHRO_SV(TPDIAFILE,"DW_THVSV",  &
       "Downdraft resolved <ThvSv> variance","K kg kg-1",XLES_DOWNDRAFT_ThvSv,HLES_AVG)

    CALL LES_DIACHRO_SV(TPDIAFILE,"DW_WSV  ", &
       "Downdraft resolved <wSv> vertical flux","m kg kg-1 s-1",XLES_DOWNDRAFT_WSv,HLES_AVG)
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
!*      3.   surface normalization parameters
!            --------------------------------
!
IF (HLES_AVG==' ' .OR. HLES_AVG=='A') THEN

  CALL LES_DIACHRO_SURF(TPDIAFILE,"Q0      ",  &
     "Sensible heat flux at the surface","m K s-1",XLES_Q0,HLES_AVG)

  IF (LUSERV) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"E0      ",  &
     "Latent heat flux at the surface","kg kg-1 m s-1",XLES_E0,HLES_AVG)
     !writes  sw and lw flux and dthrad at all levels
  CALL LES_DIACHRO(TPDIAFILE,"SWU      ",  &
     "sw_up ","W m-2 ",XLES_SWU,HLES_AVG)

  CALL LES_DIACHRO(TPDIAFILE,"SWD      ",  &
     "sw_down ","W m-2 ",XLES_SWD,HLES_AVG)
     
  CALL LES_DIACHRO(TPDIAFILE,"LWU      ",  &
     "lw_up ","W m-2 ",XLES_LWU,HLES_AVG)

  CALL LES_DIACHRO(TPDIAFILE,"LWD      ",  &
     "lw_down ","W m-2 ",XLES_LWD,HLES_AVG)

  CALL LES_DIACHRO(TPDIAFILE,"DTHRADSW      ",  &
     "dthrad_sw ","K s-1 ",XLES_DTHRADSW,HLES_AVG)

  CALL LES_DIACHRO(TPDIAFILE,"DTHRADLW      ",  &
     "dthrad_lw ","K s-1 ",XLES_DTHRADLW,HLES_AVG)
!writes mean_effective radius at all levels
  CALL LES_DIACHRO(TPDIAFILE,"RADEFF      ",  &
     "mean effective radius ","microm ",XLES_RADEFF,HLES_AVG)

  IF (NSV>0) &
  CALL LES_DIACHRO_SURF_SV(TPDIAFILE,"SV0     ",  &
     "Scalar variable fluxes at the surface","kg kg-1 m s-1",XLES_SV0,HLES_AVG)

  CALL LES_DIACHRO_SURF(TPDIAFILE,"U*      ",  &
     "Friction velocity","m s-1",XLES_USTAR,HLES_AVG)

  CALL LES_DIACHRO_SURF(TPDIAFILE,"W*      ",  &
     "Convective velocity","m s-1",XLES_WSTAR,HLES_AVG)

  CALL LES_DIACHRO_SURF(TPDIAFILE,"BL_H    ",  &
     "Boundary Layer Height","m",XLES_BL_HEIGHT,HLES_AVG)

  CALL LES_DIACHRO_SURF(TPDIAFILE,"L_MO    ",  &
     "Monin-Obukhov length","m",XLES_MO_LENGTH,HLES_AVG)

  CALL LES_DIACHRO_SURF(TPDIAFILE,"INT_TKE    ",  &
     "Vertical integrated tke","m2.s-2",XLES_INT_TKE,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"ZCB    ",  &
     "Cloud base Height","m",XLES_ZCB,HLES_AVG)   

  IF (LUSERC) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"ZCFTOT    ",  &
     "Total Cloud cover","1",XLES_CFtot,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"ZCF2TOT    ",  &
     "Total Cloud cove 2r","1",XLES_CF2tot,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"LWP    ",  &
     "Liquid Water path","kg m-2",XLES_LWP,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"LWPVAR ",  &
     "Liquid Water path variance","kg m-4",XLES_LWPVAR,HLES_AVG)

  IF (LUSERR) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"RWP    ",  &
     "Rain Water path","kg m-2",XLES_RWP,HLES_AVG)

  IF (LUSERI) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"IWP    ",  &
     "Ice Water path","kg m-2",XLES_IWP,HLES_AVG)

  IF (LUSERS) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"SWP    ",  &
     "Snow Water path","kg m-2",XLES_SWP,HLES_AVG)

  IF (LUSERG) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"GWP    ",  &
     "Graupel Water path","kg m-2",XLES_GWP,HLES_AVG)

  IF (LUSERH) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"HWP    ",  &
     "Hail Water path","kg m-2",XLES_HWP,HLES_AVG)

  IF (LUSERR) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"PREC_FRAC    ",  &
  "Fract of col where rain at surface","",XLES_PRECFR,HLES_AVG)

  IF (LUSERR) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"INST_PREC    ",  &
     "Inst precip rate","mm day-1",XLES_INPRR,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"INST_SEDIM   ",  &
     "Inst cloud precip rate","mm day-1",XLES_INPRC,HLES_AVG)

  IF (LUSERC .AND. (LDEPOSC .OR. LDEPOC)) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"INST_DEPOS   ",  &
     "Inst cloud deposi rate","mm day-1",XLES_INDEP,HLES_AVG)

  IF (LUSERR) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"RAIN_PREC    ",  &
     "inst pr. rate over rainy grid cells","mm day-1",XLES_RAIN_INPRR,HLES_AVG)
     
  IF (LUSERR) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"ACCU_PREC    ",  &
     "Accu precip rate","mm",XLES_ACPRR,HLES_AVG)   

  IF (LUSERC) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"ZMAXCF    ",  &
     "Height of Cloud fraction max","m",XLES_ZMAXCF,HLES_AVG)

  IF (LUSERC) &
  CALL LES_DIACHRO_SURF(TPDIAFILE,"ZMAXCF2   ",  &
     "Height of Cloud fraction2max","m",XLES_ZMAXCF2,HLES_AVG)

END IF
!
!-------------------------------------------------------------------------------
!
!*      4.   LES budgets
!            -----------
!
CALL WRITE_LES_BUDGET_n(TPDIAFILE,HLES_AVG)
IF (LUSERV) CALL WRITE_LES_RT_BUDGET_n(TPDIAFILE,HLES_AVG)
IF (NSV>0)  CALL WRITE_LES_SV_BUDGET_n(TPDIAFILE,HLES_AVG)
!
!-------------------------------------------------------------------------------
!
!*      5.   (ni,z,t) and (nj,z,t) 2points correlations
!            ------------------------------------------
!
IF (HLES_AVG==' ' .OR. HLES_AVG=='A') THEN
  IF (NSPECTRA_K>0) THEN
    CALL LES_DIACHRO_2PT(TPDIAFILE,"UU   ","U*U     2 points correlations", &
  "m2 s-2",XCORRi_UU,    XCORRj_UU,HLES_AVG)
    CALL LES_DIACHRO_2PT(TPDIAFILE,"VV   ","V*V     2 points correlations", &
  "m2 s-2",XCORRi_VV,    XCORRj_VV,HLES_AVG)
    CALL LES_DIACHRO_2PT(TPDIAFILE,"WW   ","W*W     2 points correlations", &
  "m2 s-2",XCORRi_WW,    XCORRj_WW,HLES_AVG)
    CALL LES_DIACHRO_2PT(TPDIAFILE,"UV   ","U*V     2 points correlations", &
  "m2 s-2",XCORRi_UV,    XCORRj_UV,HLES_AVG)
    CALL LES_DIACHRO_2PT(TPDIAFILE,"WU   ","W*U     2 points correlations", &
  "m2 s-2",XCORRi_WU,    XCORRj_WU,HLES_AVG)
    CALL LES_DIACHRO_2PT(TPDIAFILE,"WV   ","W*V     2 points correlations", &
  "m2 s-2",XCORRi_WV,    XCORRj_WV,HLES_AVG)
    CALL LES_DIACHRO_2PT(TPDIAFILE,"THTH ","Th*Th   2 points correlations", &
  "K2   ",XCORRi_ThTh,  XCORRj_ThTh,HLES_AVG)
    IF (LUSERC) &
    CALL LES_DIACHRO_2PT(TPDIAFILE,"TLTL ","Thl*Thl 2 points correlations", &
  "K2   ",XCORRi_ThlThl,XCORRj_ThlThl,HLES_AVG)
    CALL LES_DIACHRO_2PT(TPDIAFILE,"WTH  ","W*Th    2 points correlations", &
  "m K s-1 ",XCORRi_WTh,   XCORRj_WTh,HLES_AVG)
    IF (LUSERC) &
    CALL LES_DIACHRO_2PT(TPDIAFILE,"WTHL ","W*Thl   2 points correlations", &
  "m K s-1 ",XCORRi_WThl,  XCORRj_WThl,HLES_AVG)
    !
    IF (LUSERV) THEN
      CALL LES_DIACHRO_2PT(TPDIAFILE,"RVRV ","rv*rv   2 points correlations", &
  "kg2 kg-2 ",XCORRi_RvRv,  XCORRj_RvRv,HLES_AVG)
      CALL LES_DIACHRO_2PT(TPDIAFILE,"THRV ","th*rv   2 points correlations", &
  "K kg kg-1  ",XCORRi_ThRv,  XCORRj_ThRv,HLES_AVG)
      IF (LUSERC) &
      CALL LES_DIACHRO_2PT(TPDIAFILE,"TLRV ","thl*rv  2 points correlations", &
  "K kg kg-1  ",XCORRi_ThlRv, XCORRj_ThlRv,HLES_AVG)
      CALL LES_DIACHRO_2PT(TPDIAFILE,"WRV  ","W*rv    2 points correlations", &
  "m kg s-1 kg-1",XCORRi_WRv,   XCORRj_WRv,HLES_AVG)
    END IF
    IF (LUSERC) THEN
      CALL LES_DIACHRO_2PT(TPDIAFILE,"RCRC ","rc*rc   2 points correlations", &
  "kg2 kg-2 ",XCORRi_RcRc,  XCORRj_RcRc,HLES_AVG)
      CALL LES_DIACHRO_2PT(TPDIAFILE,"THRC ","th*rc   2 points correlations", &
  "K kg kg-1  ",XCORRi_ThRc,  XCORRj_ThRc,HLES_AVG)
      CALL LES_DIACHRO_2PT(TPDIAFILE,"TLRC ","thl*rc  2 points correlations", &
  "K kg kg-1  ",XCORRi_ThlRc, XCORRj_ThlRc,HLES_AVG)
      CALL LES_DIACHRO_2PT(TPDIAFILE,"WRC  ","W*rc    2 points correlations", &
  "m kg s-1 kg-1",XCORRi_WRc,   XCORRj_WRc,HLES_AVG)
    END IF
    IF (LUSERI) THEN
      CALL LES_DIACHRO_2PT(TPDIAFILE,"RCRC ","ri*ri   2 points correlations", &
  "kg2 kg-2 ",XCORRi_RiRi,  XCORRj_RiRi,HLES_AVG)
      CALL LES_DIACHRO_2PT(TPDIAFILE,"THRC ","th*ri   2 points correlations", &
  "K kg kg-1  ",XCORRi_ThRi,  XCORRj_ThRi,HLES_AVG)
      CALL LES_DIACHRO_2PT(TPDIAFILE,"TLRC ","thl*ri  2 points correlations", &
  "K kg kg-1  ",XCORRi_ThlRi, XCORRj_ThlRi,HLES_AVG)
      CALL LES_DIACHRO_2PT(TPDIAFILE,"WRC  ","W*ri    2 points correlations", &
  "m kg s-1 kg-1",XCORRi_WRi,   XCORRj_WRi,HLES_AVG)
    END IF
    DO JSV=1,NSV
      WRITE (YGROUP,FMT="(A2,I3.3)") "SS",JSV
      CALL LES_DIACHRO_2PT(TPDIAFILE,YGROUP,"Sv*Sv   2 points correlations", &
  "kg2 kg-2 ",XCORRi_SvSv(:,:,:,JSV),  XCORRj_SvSv(:,:,:,JSV),HLES_AVG)
    END DO
    DO JSV=1,NSV
      WRITE (YGROUP,FMT="(A2,I3.3)") "WS",JSV
      CALL LES_DIACHRO_2PT(TPDIAFILE,YGROUP,"W*Sv    2 points correlations", &
 "m kg s-1 kg-1",XCORRi_WSv(:,:,:,JSV),   XCORRj_WSv(:,:,:,JSV),HLES_AVG)
    END DO
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
!*      6.   spectra and time-averaged profiles (if first call to WRITE_LES_n)
!            ----------------------------------
!
IF (HLES_AVG==' ') CALL LES_SPEC_n(TPDIAFILE)
!
!-------------------------------------------------------------------------------
!
!*      7.   deallocations
!            -------------
!
DEALLOCATE(XLES_CURRENT_TRAJT )
DEALLOCATE(XLES_CURRENT_Z     )
DEALLOCATE(XLES_CURRENT_DATIME)

IF (CLES_NORM_TYPE/='NONE' ) THEN
  DEALLOCATE(XLES_NORM_M  )
  DEALLOCATE(XLES_NORM_S  )
  DEALLOCATE(XLES_NORM_K  )
  DEALLOCATE(XLES_NORM_RHO)
  DEALLOCATE(XLES_NORM_RV )
  DEALLOCATE(XLES_NORM_SV )
  DEALLOCATE(XLES_NORM_P  )
END IF
!
END SUBROUTINE WRITE_LES_n 
