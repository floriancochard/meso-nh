!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_nesting.f90, Version:1.3, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     ###################
      MODULE MODD_NESTING
!     ###################
!
!!****  *MODD_NESTING* - declaration of gridnesting configuration variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the variables
!     which concern the gridnesting configuration of all models.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_PARAMETERS  :
!!         JPMODELMAX : Maximum allowed  number of nested models
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_NESTING)
!!       
!!    AUTHOR
!!    ------
!!	J.P. Lafore   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/08/95
!!      updated     29/07/96  (J.P. Lafore) MY_NAME(m) introduction          

!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
                            ! resolution RATIO between models m and its father NDAD(m)
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NDXRATIO_ALL        ! in x-direction 
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NDYRATIO_ALL        ! in y-direction 
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NDTRATIO            ! in Time 
!
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NXOR_ALL, NYOR_ALL  ! horizontal position (i,j) of the
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NXEND_ALL,NYEND_ALL ! ORigin and END of model m 
                                                     ! relative to its father NDAD(m)    
!
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NDAD ! model number of the father of each model "m"
REAL,SAVE,     DIMENSION(JPMODELMAX) :: XWAY ! model m interactive nesting level with its father NDAD(m)
!
                                                            !   MeSsaGes concerning 
INTEGER,SAVE,  DIMENSION(JPMODELMAX,JPMODELMAX) :: NMSG_IF  ! var. Interpolation at Flux
INTEGER,SAVE,  DIMENSION(JPMODELMAX,JPMODELMAX) :: NMSG_IS  ! and Scalar location
INTEGER,SAVE,  DIMENSION(JPMODELMAX,JPMODELMAX) :: NMSG_AVR ! AVeRage
INTEGER,SAVE,  DIMENSION(JPMODELMAX,JPMODELMAX) :: NMSG_END ! timestep END
                                                            !   MeSsaGes concerning
INTEGER,SAVE,  DIMENSION(JPMODELMAX,JPMODELMAX) :: NMSG_AVR_END ! AVeRage END
!
CHARACTER(LEN=28),SAVE,   DIMENSION(JPMODELMAX) :: CMY_NAME,CDAD_NAME
                                                  ! names of the initial FM-Files
                                                  ! then generic names of output FM-Files
                                                  ! of each model "m"
                                                  ! and of its DAD model
                                                  ! (read and written on the LFI parts)
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NDT_2_WAY ! number of times the time step
              ! of model n used for the relaxation time of the 2_WAY grid-nesting
              ! interaction  i.e. Tau = NDT_2_WAY * XTSTEP
END MODULE MODD_NESTING
