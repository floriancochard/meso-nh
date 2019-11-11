!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_lunitn.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     ###################
      MODULE MODD_LUNIT1
!     ###################
!
!!****  *MODD_LUNIT1* - declaration of names and logical unit numbers of files 
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the names 
!     for the initial Meso-NH files  
!     and also the  generic names  for the output files for model n.    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_LUNITn)
!!          
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      05/05/94  
!!      Modification  20/10/94 (J.Stein) add the output files                    
!!      Modification  10/03/95 (I.Mallet)   add the coupling files names 
!!      Modification  25/09/95 (J.Stein) add the output diachronic file                    
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
! 
CHARACTER(LEN=28),SAVE :: CINIFILE      ! Name of the input FM-file
CHARACTER(LEN=24),SAVE :: COUTFILE      ! Generic name of the output FM-files
CHARACTER(LEN=28),SAVE :: CFMDIAC       ! diachronic output FM-file 
!
CHARACTER(LEN=16),SAVE :: CLUOUT        ! Name of output_listing file
CHARACTER(LEN=28),SAVE,DIMENSION(JPCPLFILEMAX) :: CCPLFILE ! Names of the 
                                                           ! coupling FM-files
!
END MODULE MODD_LUNIT1
