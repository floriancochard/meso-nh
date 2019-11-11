!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!!    ################################ 
      MODULE MODI_INIT_FOREFIRE_n
!!    ################################ 

      INTERFACE
      
      	SUBROUTINE INIT_FOREFIRE_n(KMODEL, KLUOUT, IP, KYEAR, KMONTH, KDAY, PTIME, DT)
      	
      		IMPLICIT NONE
   			INTEGER, INTENT(IN)    	:: KMODEL    	! number of the model in nest hierarchy
      		INTEGER, INTENT(IN)    	:: KLUOUT    	! output listing
				INTEGER, INTENT(IN)		:: IP				! mpi rank of the process
				INTEGER, INTENT(IN)  	:: KYEAR     	! current year (UTC)
				INTEGER, INTENT(IN)  	:: KMONTH   	! current month (UTC)
				INTEGER, INTENT(IN)  	:: KDAY      	! current day (UTC)
				REAL, INTENT(IN)  		:: PTIME     	! current time since midnight (UTC, s)
				REAL, INTENT(IN)  		:: DT     		! time-step of MNH (s)
				
      	END SUBROUTINE INIT_FOREFIRE_n
      	
      END INTERFACE 
      
      END MODULE MODI_INIT_FOREFIRE_n

!     ##################################################
      SUBROUTINE INIT_FOREFIRE_n(KMODEL, KLUOUT, IP, KYEAR, KMONTH, KDAY, PTIME, DT)
!     ##################################################

!!    AUTHORS
!!    -------
!!      P. Tulet  CNRM
!!      X. Pialat SPE

USE MODD_FOREFIRE
USE MODD_FOREFIRE_n
USE MODD_FOREFIRE_FORTRAN_API

USE MODD_DYN_n, 		ONLY: NSTOP
USE MODD_GRID,      	ONLY: XLATORI, XLONORI
USE MODD_GRID_n,    	ONLY: XXHAT, XYHAT, XZHAT, XZS, XZZ
USE MODD_NSV,       	ONLY: NSV_CHEMBEG, NSV_CHEMEND, NSV_FF, NSV_FFBEG, NSV_FFEND, CSV
USE MODD_CH_MNHC_n, 	ONLY: LUSECHEM
USE MODD_CH_M9_n,   	ONLY: CNAMES

   INTEGER, INTENT(IN)    	:: KMODEL    	! number of the model in nest hierarchy
   INTEGER, INTENT(IN)    	:: KLUOUT    	! output listing
	INTEGER, INTENT(IN)		:: IP				! mpi rank of the process
	INTEGER, INTENT(IN)  	:: KYEAR     	! current year (UTC)
	INTEGER, INTENT(IN)  	:: KMONTH    	! current month (UTC)
	INTEGER, INTENT(IN)  	:: KDAY      	! current day (UTC)
	REAL, INTENT(IN)  		:: PTIME     	! current time since midnight (UTC, s)
	REAL, INTENT(IN)  		:: DT     		! time-step of MNH (s)
	
	INTEGER                 :: IINFO_ll       ! return code of parallel routine
	INTEGER						:: JSV, JCHEMV
	
	LOGICAL 						:: DUMPEDGRID 

	WRITE (KLUOUT,*)  'COUPLING WITH FOREFIRE IN MODEL ', KMODEL
	
	!******************************!
	!* COMMON PART FOR ALL MODELS *!
	!******************************!
	
	! INITIALISATION OF THE FOREFIRE SIMULATION
	!------------------------------------------
	CALL MNH_TO_FF_INIT(PTIME)
   
   ! SIZE OF THE MATRICES TO RETRIEVE THE SURFACE PROPERTIES
   !--------------------------------------------------------
   FF_NX = SIZE(XXHAT)
   FF_NY = SIZE(XYHAT)
   FF_NZ = SIZE(XZHAT)
   FF_MATRIXSIZE = FF_NX*FF_NY
   
   ! REFERENCE FOR THE OUTPUTS
   !--------------------------
	FFNMODEL = KMODEL
	PROCID = IP
	FFOUTPUTSUPDATE = FFOUTUPS(KMODEL)
	IF ( FFOUTPUTSUPDATE .LE. 0 ) THEN
		FFOUTPUTSUPDATE = 1000000.
	ELSE
		CALL MNH_DUMP_DOUBLEARRAY(FFNMODEL, PROCID, sZgrid, FF_TIME, XZZ, FF_NX*FF_NY*FF_NZ, FF_NX, FF_NY, FF_NZ, 1)
	ENDIF
   FFREFERENCETIME = PTIME
   FF_TIME = 0.
   FFNUMOUT = 0;

   ! ARE FOREFIRE 3D OUTPUTS REQUESTED BY THE USER ?
   !------------------------------------------------
   
   IF ( FLOWOUT(FFNMODEL) .GT. 0 ) THEN
      FF3DOUTPUTSFLOW = .TRUE.
   ELSE
      FF3DOUTPUTSFLOW = .FALSE.
   END IF
   IF ( PHYSOUT(FFNMODEL) .GT. 0 ) THEN
      FF3DOUTPUTSPHYS = .TRUE.
   ELSE
      FF3DOUTPUTSPHYS = .FALSE.
   END IF
   IF ( CHEMOUT(FFNMODEL) .GT. 0 ) THEN
      FF3DOUTPUTSCHEM = .TRUE.
   ELSE
      FF3DOUTPUTSCHEM = .FALSE.
   END IF
   	
	!************************************************!
	!* SPECIFIC PART IS THE MODEL IS COUPLED TO MNH *!
	!************************************************!
  
   IF ( FFCOUPLING ) THEN 
   	! CREATING THE FOREFIRE DOMAIN
   	!-----------------------------
   	CALL MNH_TO_FF_CREATEDOMAIN(IP, KYEAR, KMONTH, KDAY, PTIME, XLATORI, XLONORI &
   			, FF_NX, XXHAT, FF_NY, XYHAT, FF_NZ, FF_NX*FF_NY*FF_NZ, XZZ, DT)

   	! SETTING THE FIRE SIMULATION TO PARALLEL
   	!----------------------------------------
   	CALL MNH_PUT_INT(sParallel, 1)

   		! COMMUNICATING THE NUMBER OF MNH ITERATIONS
   	!-------------------------------------------
   	CALL MNH_PUT_INT(sNumIte, NSTOP)
	
	   ! COMMUNICATION OF THE OROGRAPHY TO FOREFIRE
	   !-------------------------------------------
	   CALL MNH_PUT_DOUBLEARRAY(sAltitude, FF_TIME, XZS, FF_MATRIXSIZE, 1)
	
	   ! ALLOCATION OF THE FLUXES STEMMING FROM FOREFIRE
	   !------------------------------------------------
	   ALLOCATE(FF_HEATFLUX(FF_NX, FF_NY))
	 	CALL MNH_TO_FF_CHECKLAYER(sHeatFlux)
	   ALLOCATE(FF_VAPORFLUX(FF_NX, FF_NY))
	 	CALL MNH_TO_FF_CHECKLAYER(sVaporFlux)
	   ALLOCATE(FF_SVFLUXES(FF_NX, FF_NY, NSV_FF))
	   ALLOCATE(sScalarVariables(NSV_FF))
	   DO JSV = 1, NSV_FF
	 		CSV(NSV_FFBEG+JSV-1) = TRIM(FFSVNAMES(JSV))
	 		sScalarVariables(JSV) = cast_char_to_c(FFSVNAMES(JSV))
	 		CALL MNH_TO_FF_CHECKLAYER(sScalarVariables(JSV))
	 	ENDDO
	
		! RETRIEVING INFORMATION IN CASE OF CHEMICAL COUPLING
		!----------------------------------------------------
		IF ( LFFCHEM ) THEN
			IF ( .NOT.LUSECHEM ) THEN
				WRITE(*,*) 'PROBLEM: YOU SWITCHED CHEMISTRY ON IN FOREFIRE WHILE OFF IN MNH'
			ENDIF
			ALLOCATE(FF_CHEMINDICES(NFFCHEMVAR))
	   	ALLOCATE(FF_CVFLUXES(FF_NX, FF_NY, NFFCHEMVAR))
	   	ALLOCATE(sChemicalVariables(NFFCHEMVAR))
			DO JSV = 1, NFFCHEMVAR
	 			sChemicalVariables(JSV) = cast_char_to_c(FFCVNAMES(JSV))
	 			CALL MNH_TO_FF_CHECKLAYER(sChemicalVariables(JSV))
		 		DO JCHEMV = 1, NSV_CHEMEND-NSV_CHEMBEG
		   		IF ( TRIM(CNAMES(JCHEMV)) == TRIM(FFCVNAMES(JSV)) ) THEN
		   			FF_CHEMINDICES(JSV) = JCHEMV
		   		END IF
		 		ENDDO
		 	ENDDO
			ALLOCATE(FF_CHEMINDOUT(NFFCHEMVAROUT))
			DO JSV = 1, NFFCHEMVAROUT
		 		DO JCHEMV = 1, NSV_CHEMEND-NSV_CHEMBEG
		   		IF ( TRIM(CNAMES(JCHEMV)) == TRIM(FFCONAMES(JSV)) ) THEN
		   			FF_CHEMINDOUT(JSV) = JCHEMV
		   		END IF
		 		ENDDO
		 	ENDDO
		END IF
	   
	   ! RETRIEVING THE MAXIMUM NUMBER OF FIRENODES PER MNH CELL
	   !--------------------------------------------------------
	   CALL MNH_GET_INT(sNumFNMax, NBFNMAX)
	   FF_PARALMATRIXSIZE = FF_MATRIXSIZE*NBFNMAX
	
	   ! ALLOCATION OF THE MATRICES FOR FOREFIRE PARALLELISATION
	   !--------------------------------------------------------
	   ALLOCATE(FFNODES_POSX(FF_NX,FF_NY,NBFNMAX))
	   ALLOCATE(FFNODES_POSY(FF_NX,FF_NY,NBFNMAX))
	   ALLOCATE(FFNODES_VELX(FF_NX,FF_NY,NBFNMAX))
	   ALLOCATE(FFNODES_VELY(FF_NX,FF_NY,NBFNMAX))
	   ALLOCATE(FFNODES_TIME(FF_NX,FF_NY,NBFNMAX))
	   ALLOCATE(FFNODES_ID(FF_NX,FF_NY,NBFNMAX))
	   ALLOCATE(FFOUTERWINDU(FF_NX,FF_NY))
	   ALLOCATE(FFOUTERWINDV(FF_NX,FF_NY))
	   FFOUTERWINDU(:,:) = 0.
	   FFOUTERWINDV(:,:) = 0.
	   CALL MNH_GET_INT(sMNHmult, FFMULT)
		
		! PARALLELIZATION PROCESS
		!------------------------
		CALL FOREFIRE_SEND_PARAL_n(IINFO_ll)
   ENDIF


END SUBROUTINE INIT_FOREFIRE_n
