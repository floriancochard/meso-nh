!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /srv/cvsroot/MNH-VX-Y-Z/src/MNH/modn_budget.f90,v $ $Revision: 1.2.2.1.2.1.2.1.10.1.2.3 $ $Date: 2014/01/09 15:01:56 $
!-----------------------------------------------------------------
!     ##################
      MODULE MODN_BUDGET
!     ##################
!
!!****  *MODN_BUDGET* - declaration of namelist NAM_BUDGET
!!
!!    PURPOSE
!!    -------
!       The purpose of this  module is to specify  the namelists NAM_GENBUDGET,
!     NAM_BURU, NAM_BURV, NAM_BURW, NAM_BURTH, NAM_BURTKE, NAM_BURSV, 
!     NAM_BURRV, NAM_BURRC, NAM_BURR, NAM_BURI, NAM_BURS, NAM_BURG, NAM_BURH 
!     which concern the budgets activation and the choice of processes.    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_BUDGET : contains declaration of budget variables
!!
!!         CBUTYPE      : type of desired budget
!!                          'CART' for cartesian box configuration
!!                          'MASK' for budget zone defined by a mask 
!!                          'NONE'  ' for no budget
!!         NBUMOD       : model number in which budget is calculated
!!         XBULEN       : length of the budget temporal average in seconds
!!         NBUKL, NBUKH : lowest and highest K indice values
!!         LBU_KCP      : switch for compression in the K direction
!!                       .TRUE. = compression in the K direction
!!                       .FALSE. = no compression in the K direction
!!         XBUWRI       : Period in seconds when the budget is written 
!!                       on FM-files
!!         
!!         Variables used by the cartesian box case ('CART') only 
!!         
!!         NBUIL, NBUIH : lowest and highest I indice values of cartesian box 
!!         NBUJL, NBUJH : lowest and highest J indice values of cartesian box
!!         LBU_ICP      : switch for compression in I direction
!!                       .TRUE. = compression in the I direction
!!                       .FALSE. = no compression in the I direction
!!         LBU_JCP      : switch for compression in the J direction 
!!                       .TRUE. = compression in the J direction
!!                       .FALSE. = no compression in the J direction
!!         
!!         Variable used by the  mask case ('MASK') only
!!         
!!         NBUMASK      : number of MASK zones for which budgets are performed
!!
!!         Logicals for budgets activations:
!!         
!!         LBU_RU       : logical for budget of RU (wind component along x)
!!                        .TRUE. = budget of RU         
!!                        .FALSE. = no budget of RU 
!!         LBU_RV       : logical for budget of RV (wind component along y)
!!                        .TRUE. = budget of RV         
!!                        .FALSE. = no budget of RV 
!!         LBU_RW        : logical for budget of RW (wind component along z)
!!                        .TRUE. = budget of RW         
!!                        .FALSE. = no budget of RW 
!!         LBU_RTH      : logical for budget of RTH (potential temperature)
!!                        .TRUE. = budget of RTH        
!!                        .FALSE. = no budget of RTH
!!         LBU_RTKE     : logical for budget of RTKE (turbulent kinetic energy)
!!                        .TRUE. = budget of RTKE       
!!                        .FALSE. = no budget of RTKE
!!         LBU_RRV      : logical for budget of RRV (water vapor)
!!                        .TRUE. = budget of RRV 
!!                        .FALSE. = no budget of RRV 
!!         LBU_RRC      : logical for budget of RRC (cloud water)
!!                        .TRUE. = budget of RRC 
!!                        .FALSE. = no budget of RRC 
!!         LBU_RRR      : logical for budget of RRR (rain water)
!!                        .TRUE. = budget of RRR 
!!                        .FALSE. = no budget of RRR 
!!         LBU_RRI      : logical for budget of RRI (ice)
!!                        .TRUE. = budget of RRI 
!!                        .FALSE. = no budget of RRI 
!!         LBU_RRS      : logical for budget of RRS (snow)
!!                        .TRUE. = budget of RRS 
!!                        .FALSE. = no budget of RRS 
!!         LBU_RRG      : logical for budget of RRG (graupel)
!!                        .TRUE. = budget of RRG 
!!                        .FALSE. = no budget of RRG 
!!         LBU_RRH      : logical for budget of RRH (hail)
!!                        .TRUE. = budget of RRH 
!!                        .FALSE. = no budget of RRH 
!!         LBU_RSV      : logical for budget of RSV1 (scalar variable)
!!                        .TRUE. = budget of RSV   
!!         
!!         
!!                         Variables specific to a budget
!!     Their values are integer from 0 to JPBUPROMAX (max. of possible 
!!     processes). The mechanism for activating a given process is detailled 
!!     in Book 3. Each variable name is related to the location of the 
!!     process in the run and may be decomposed in 2 parts: the first part
!!     gives the name of the process implicated and the second part gives 
!!     the name of the budget of interest. For example, NADVXU is linked to 
!!     the budget of the variable RU in the process ADVection  along X. The 
!!     syntaxic code is detailled above:
!!
!!              Left part of the variable name = implicated process 
!!       
!!     ADVX  : advection along X (subroutine ADVECTION)        
!!     ADVY  : advection along Y (subroutine ADVECTION)        
!!     ADVZ  : advection along Z (subroutine ADVECTION)        
!!     CURV  : curvature terms (subroutine DYN_SOURCES)    
!!     COR   : coriolis terms (subroutine DYN_SOURCES)   
!!     GRAV  : gravity term (subroutine DYN_SOURCES) 
!!     DIF   : numerical diffusion (subroutine NUM_DIFF)    
!!     REL   : relaxation (subroutine RELAXATION)  
!!     HTURB : horizontal turbulence (subroutine TURB)        
!!     VTURB : vertical turbulence (subroutine TURB)         
!!     PRES  : pressure term (subroutine PRESSURE)
!!     COND  : evaporation/condensation (subroutine CLOUD_EXPLI)
!!     REVA  : rain evaporation (subroutine CLOUD_EXPLI)
!!     ACCR  : accretion (subroutine CLOUD_EXPLI)
!!     AUTO  : autoconversion (subroutine CLOUD_EXPLI)
!!     SEDI  : sedimentation (subroutine CLOUD_EXPLI)       
!!         
!!              Right part of the variable name = budget of interest
!!
!!                                                          Budget number  
!!     U     : budget of RU                                      1
!!     V     : budget of RV                                      2
!!     W     : budget of RW                                      3
!!     TH    : budget of RTH                                     4    
!!     TKE   : budget of RTKE (turbulent kinetic energy)         5
!!     RV    : budget of RRV variable (water)                    6
!!     RC    : budget of RRC variable (water)                    7   
!!     RR    : budget of RRR variable (water)                    8  
!!     RI    : budget of RRI variable (water)                    9    
!!     RS    : budget of RRS variable (water)                    10    
!!     RG    : budget of RRG variable (water)                    11    
!!     RH    : budget of RRH variable (water)                    12    
!!     SV1   : budget of RTRACER1 Scalar Variable                13
!!              Right part of the variable name = budget of interest
!!         
!!         
!!     It is now easy to deduce the process and the budget from every variable 
!!     name:
!!                    Budget of RU      
!!         
!!         NADVXU, NADVYU, NADVZU, NCURVU, NCORU, NDIFU, NRELU, NHTURBU,  
!!         NVTURBU, NPRESU,NMAFLU
!!         
!!                    Budget of RV          
!!         
!!         NADVXV, NADVYV, NADVZV, NCURVV, NCORV, NDIFV, NRELV, NHTURBV,  
!!         NVTURBV, NPRESV, NMAFLV
!!         
!!                    Budget of RW       
!!         
!!         NADVXW, NADVYW, NADVZW, NCURVW, NCORW, NGRAVW, NDIFW, NRELW,  
!!         NHTURBW, NVTURBW, NPRESW
!!         
!!                    Budget of RTH         
!!         
!!         NADVXTH, NADVYTH, NADVZTH, NPREFTH, NDIFTH, NRELTH, NHTURBTH, 
!!         NVTURBTH, NMAFLTH, NREVATH, NCONDTH        
!!         
!!                    Budget of RTKE       
!!         
!!         NADVXTKE, NADVYTKE, NADVZTKE, NDIFTKE, NDPTKE, NTPTKE, NDISSTKE,
!!         NTRTKE
!!         
!!                    Budget of RSV1        
!!         
!!         NADVXSV1, NADVYSV1,  NADVZSV1, NDIFSV1, NHTURBSV1, NVTURBSV1 
!!         
!!                    Budget of RRV
!!         
!!         NADVXRV, NADVYRV, NADVZRV, NDIFRV, NRELRV, NHTURBRV, NVTURBRV,
!!         NREVARV, NCONDRV, NMAFLRV
!!
!!                    Budget of RRC
!!         
!!         NADVXRC, NADVYRC, NADVZRC, NDIFRC, NHTURBRC, NVTURBRC, NACCRRC,
!!         NAUTORC, NCONDRC
!!
!!                    Budget of RRR
!!         
!!         NADVXRR, NADVYRR, NADVZRR, NDIFRR, NACCRRR, NAUTORR, NREVARR, 
!!         NSEDIRR
!!
!!                    Budget of RRI
!!         
!!         NADVXRI, NADVYRI, NADVZRI, NDIFRI
!!
!!                    Budget of RRS
!!         
!!         NADVXRS, NADVYRS, NADVZRS, NDIFRS
!!
!!                    Budget of RRG 
!!         
!!         NADVXRG, NADVYRG, NADVZRG, NDIFRG
!!
!!                    Budget of RRH
!!         
!!         NADVXRH, NADVYRH, NADVZRH, NDIFRH
!!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_BUDGET)
!!      Asencio N. et al., 1994, "Le projet de modele non-hydrostatique 
!!    commun CNRM-LA, specifications techniques", Note CNRM/GMME, 26, 139p,
!!    (chapters 2 and 3)
!!       
!!    AUTHOR
!!    ------
!!	P. Hereil   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/03/95                      
!!      J. Stein    29/06/95  new processes' list
!!      J.-P. Pinty 11/01/97  add several SVx
!!      J.-P. Pinty 18/02/97  add forcing and ice
!!      J.-P. Pinty 25/09/00  add budget terms for C2R2
!!      D. Gazen    22/01/01  add NCHEMSV
!!      C.Lac           04/2016  negative contribution to the budget splitted between advection, turbulence and microphysics for KHKO/C2R2
!!      C. Barthe        /16  add budget terms for LIMA
!!      C.Lac        10/2016  Add droplet deposition
!!      S. Riette   11/2016 New budgets for ICE3/ICE4
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_BUDGET
!
IMPLICIT NONE
! 
NAMELIST/NAM_BUDGET/CBUTYPE, NBUMOD, XBULEN, NBUKL, NBUKH, LBU_KCP, XBUWRI, &
                    NBUIL, NBUIH, NBUJL, NBUJH, LBU_ICP, LBU_JCP, NBUMASK 
!
NAMELIST/NAM_BU_RU/LBU_RU, NASSEU, NNESTU, NADVU, NFRCU, NNUDU, &
                   NCURVU, NCORU, NDIFU, NRELU, NDRAGU, NHTURBU, NVTURBU, NMAFLU, NPRESU  
!
NAMELIST/NAM_BU_RV/LBU_RV, NASSEV, NNESTV, NADVV, NFRCV, NNUDV, &
                   NCURVV, NCORV, NDIFV, NRELV, NDRAGV, NHTURBV, NVTURBV, NMAFLV, NPRESV  

NAMELIST/NAM_BU_RW/LBU_RW, NASSEW, NNESTW, NADVW, NFRCW, NNUDW, &
                   NCURVW, NCORW, NGRAVW, NDIFW, NRELW, NHTURBW, NVTURBW, NPRESW
!
NAMELIST/NAM_BU_RTH/LBU_RTH, NASSETH, NNESTTH, NADVTH, NFRCTH, &
                   NNUDTH, NPREFTH, NDIFTH, NRELTH, NRADTH, NDCONVTH, NHTURBTH, &
                   NVTURBTH, NDISSHTH, NNEGATH, NREVATH, NCONDTH, NHENUTH, NHONTH, &
                   NSFRTH, NDEPSTH, NDEPGTH,NRIMTH, NACCTH, NCFRZTH, NWETGTH, &
                   NDRYGTH, NGMLTTH, NIMLTTH, NBERFITH, NCDEPITH, NWETHTH, NHMLTTH, &
                   NMAFLTH, NNETURTH, NNEADVTH,NNECONTH, NDRYHTH, NADJUTH, NCORRTH, &
                   NHINDTH, NHINCTH, NHONHTH, NHONCTH, NHONRTH, NCEDSTH, NSEDITH
!
NAMELIST/NAM_BU_RTKE/LBU_RTKE, NASSETKE, NADVTKE,    &
                     NFRCTKE, NDIFTKE, NRELTKE, NDRAGTKE,                           &
                     NDPTKE, NTPTKE, NDISSTKE, NTRTKE
!
NAMELIST/NAM_BU_RRV/LBU_RRV, NASSERV, NNESTRV, NADVRV, NFRCRV, &
                    NNUDRV, NDIFRV, NRELRV, NDCONVRV, NHTURBRV, NVTURBRV, NNEGARV, &
                    NREVARV, NCONDRV, NHENURV, NDEPSRV, NDEPGRV, NCDEPIRV, NMAFLRV, &
                    NNETURRV, NNEADVRV,NNECONRV, NADJURV, NCORRRV, NHINDRV, NHONHRV, NCEDSRV
! 
NAMELIST/NAM_BU_RRC/LBU_RRC, NASSERC, NNESTRC, NADVRC, NFRCRC, &
                    NDIFRC, NRELRC, NDCONVRC, NHTURBRC, NVTURBRC, NNEGARC, NACCRRC, &
                    NAUTORC, NCONDRC, NHONRC, NRIMRC, NWETGRC, NDRYGRC, NIMLTRC,   &
                    NBERFIRC, NCDEPIRC, NHENURC, NSEDIRC, NWETHRC, NNETURRC, &
                    NNEADVRC,NNECONRC, NDRYHRC, NADJURC, NCORRRC, NCMELRC, &
                    NHINCRC, NHONCRC, NCEDSRC, NREVARC, NDEPORC,NDEPOTRRC, &
                    NCORRRC, NR2C1RC, NCVRCRC
! 
NAMELIST/NAM_BU_RRR/LBU_RRR, NASSERR, NNESTRR, NADVRR, NFRCRR, &
                    NDIFRR, NRELRR, NNEGARR, NACCRRR, NAUTORR, NREVARR, NSEDIRR,    &
                    NSFRRR, NACCRR, NCFRZRR, NWETGRR, NDRYGRR, NGMLTRR, NWETHRR,    &
                    NHMLTRR, NDRYHRR, NCORRRR, NCMELRR,NHONRRR, NCORRRR, NR2C1RR, NCVRCRR
! 
NAMELIST/NAM_BU_RRI/LBU_RRI, NASSERI, NNESTRI, NADVRI, NFRCRI, &
                    NDIFRI, NRELRI, NDCONVRI, NHTURBRI, NVTURBRI, NNEGARI, NSEDIRI, &
                    NHENURI, NHONRI, NAGGSRI, NAUTSRI, NCFRZRI, NWETGRI, NDRYGRI,   &
                    NIMLTRI, NBERFIRI, NCDEPIRI, NWETHRI, NDRYHRI, NADJURI, NCORRRI, &
                    NHINDRI, NHINCRI, NHONHRI, NHONCRI, NCNVIRI, NCNVSRI, &
                    NHMSRI, NHMGRI, NCEDSRI, NCORRRI
! 
NAMELIST/NAM_BU_RRS/LBU_RRS, NASSERS, NNESTRS, NADVRS, NFRCRS, &
                    NDIFRS, NRELRS, NNEGARS, NSEDIRS, NDEPSRS, NAGGSRS, NAUTSRS,    &
                    NRIMRS, NACCRS, NCMELRS, NWETGRS, NDRYGRS, NWETHRS, NDRYHRS,    &
                    NCORRRS, NCNVIRS, NCNVSRS, NHMSRS, NCORRRS
! 
NAMELIST/NAM_BU_RRG/LBU_RRG, NASSERG, NNESTRG, NADVRG, NFRCRG, &
                    NDIFRG, NRELRG, NNEGARG, NSEDIRG, NSFRRG, NDEPGRG, NRIMRG, NACCRG,    &
                    NCMELRG, NCFRZRG, NWETGRG, NDRYGRG, NGMLTRG, NWETHRG, &
                    NDRYHRG, NCORRRG, NHGCVRG, NGHCVRG,NHONRRG, NHMGRG, NCOHGRG 
! 
NAMELIST/NAM_BU_RRH/LBU_RRH, NASSERH, NNESTRH, NADVRH, NFRCRH, &
                    NDIFRH, NRELRH, NNEGARH, NSEDIRH, NWETGRH, NWETHRH, NDRYHRH, NHMLTRH, &
                    NCORRRH, NHGCVRH, NGHCVRH, NCOHGRH, NHMLTRH
! 
NAMELIST/NAM_BU_RSV/ LBU_RSV, NASSESV, NNESTSV, NADVSV, NFRCSV, &
                     NDIFSV, NRELSV, NDCONVSV, NVTURBSV, NHTURBSV, NCHEMSV, NMAFLSV,       &
                     NNEGASV,                                                              & 
                     NAUTOQC, NACCRQC, NRIMQC, NWETGQC, NDRYGQC, NIMLTQC, NBERFIQC,        &
                     NDEPIQC, NINDQC, NSEDIQC, NNEUTQC,                                    &
                     NAUTOQR, NACCRQR, NREVAQR, NACCQR, NCFRZQR, NWETGQR, NDRYGQR,         &
                     NGMLTQR, NSEDIQR, NNEUTQR,                                            &
                     NAGGSQI, NAUTSQI, NCFRZQI, NWETGQI, NDRYGQI, NIMLTQI, NBERFIQI,       &
                     NDEPIQI, NNIISQI, NSEDIQI, NNEUTQI,                                   &
                     NDEPSQS, NAGGSQS, NAUTSQS, NRIMQS, NACCQS, NCMELQS, NWETGQS,          &
                     NDRYGQS, NNIISQS, NSEDIQS, NNEUTQS,                                   &
                     NDEPGQG, NRIMQG, NACCQG, NCMELQG, NCFRZQG, NWETGQG, NDRYGQG,          &
                     NGMLTQG, NINDQG, NSEDIQG, NNEUTQG,NDEPOSV,NDEPOTRSV
! must add budget for hail
!
END MODULE MODN_BUDGET
