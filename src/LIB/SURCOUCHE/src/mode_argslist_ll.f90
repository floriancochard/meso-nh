!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!!    #######################
      MODULE MODE_ARGSLIST_ll
!!    #######################
!
!!****  *MODE_ARGSLIST_ll *-
!
!!    Routines Of The User Interface
!!    ------------------------------
!
!    SUBROUTINES : ADD1DFIELD_ll, ADD2DFIELD_ll, ADD3DFIELD_ll
!                  DEL1DFIELD_ll, DEL2DFIELD_ll, DEL3DFIELD_ll
!                  CLEANLIST_ll
!
!!    Purpose
!!    -------
!     This module manages a list of fields. Fields may be added
!     (routines ADD1DFIELD_ll, ADD2DFIELD_ll, ADD3DFIELD_ll) or deleted 
!     (routines DEL1DFIELD_ll, DEL2DFIELD_ll, DEL3DFIELD_ll) from
!     the list. The list can then be used in routines where
!     all fields in the list are handled identically, e.g.
!     the distribution of the data or the update of the halos.
!     The list may contain only one type of fields (1D, 2D or 3D).
! 
!!    Reference
!!     ---------
!
!     User interface for the Meso-NH parallel package
!
!!    Authors
!!    -------
!
!     Ph. Kloos                 * CERFACS - CNRM *
!
!!    Implicit Arguments
!!     ------------------
!
!     Module MODD_ARGSLIST_ll
!         LIST_ll : list of 1D/2D/3D fields
!         LIST1D_ll : list of 1D fields
!
!!    Modifications
!!    -------------
!     Original    May 19, 1998
!
!-------------------------------------------------------------------------------
!
! Implicit arguments
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll, LIST1D_ll
!
  CONTAINS
!
!!    ##############################################
      SUBROUTINE ADD1DFIELD_ll(HDIR, TPLIST, PFIELD)
!!    ##############################################
!
!!****  *ADD1DFIELD_ll* -
!
!!    Purpose
!!    -------
!     This routine is used to add a 1D field (PFIELD) to a list of 
!     1D fields (TPLIST). 
!
!!    Reference
!!    ---------
!
!     User interface for the Meso-NH parallel package
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_ARGSLIST_ll
!         LIST1D_ll : list of 1D fields
!!    Author
!!    ------
!
!     Ph. Kloos                 * CERFACS - CNRM *
!
!!    Modifications
!!    -------------
!     Original    May 19, 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!*       0.1   declarations of arguments
!
  IMPLICIT NONE
!
  TYPE(LIST1D_ll), POINTER     :: TPLIST ! list of fields
  REAL, DIMENSION(:), TARGET   :: PFIELD ! field to be added
                                         ! to the list of fields
  CHARACTER(LEN=1), INTENT(IN) :: HDIR ! direction of the field ("X" or "Y")
!
!*       0.2   declarations of local variables
!
  TYPE(LIST1D_ll), POINTER :: TZLIST
!
!-------------------------------------------------------------------------------
!
!*       1.    Test value of HDIR
!
  IF (HDIR /= "X" .AND. HDIR /= "Y") THEN
    WRITE(*,*) 'Error ADD1DFIELD : Bad HDIR argument'
    STOP
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.    In case TPLIST has not been already used
!              ----------------------------------------
  IF (.NOT. ASSOCIATED(TPLIST)) THEN ! TPLIST is empty
!
    ALLOCATE(TPLIST)
    TPLIST%ARRAY1D => PFIELD
    NULLIFY(TPLIST%NEXT)
    TPLIST%NCARD = 1
    TPLIST%CDIR = HDIR
!
  ELSE 
!
!-------------------------------------------------------------------------------
!
!*       3.    TPLIST already contains fields;
!*             add the PFIELD field at the end of the list
!              -------------------------------------------------
!
    TZLIST => TPLIST
    DO WHILE (ASSOCIATED(TZLIST%NEXT))
      TZLIST => TZLIST%NEXT
    ENDDO
!
    ALLOCATE(TZLIST%NEXT)
    TZLIST => TZLIST%NEXT
    TZLIST%ARRAY1D => PFIELD
    TZLIST%NCARD = 0
    TZLIST%CDIR = HDIR
    NULLIFY(TZLIST%NEXT)
!
    TPLIST%NCARD = TPLIST%NCARD + 1
!
  ENDIF
! 
!-------------------------------------------------------------------------------
!
      END SUBROUTINE ADD1DFIELD_ll
!
!!    ########################################
      SUBROUTINE ADD2DFIELD_ll(TPLIST, PFIELD)
!!    ########################################
!
!!****  *ADD2DFIELD_ll* -
!
!!    Purpose
!!    -------
!     This routine is used to add a 2D field (PFIELD) to a list of 
!     2D fields (TPLIST). 
! 
!!    Reference
!!    ---------
! 
!     User interface for the Meso-NH parallel package
! 
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_ARGSLIST :
!        LIST_ll : list of fields
!
!!    Author
!!    ------
!
!     Ph. Kloos                 * CERFACS - CNRM *
!
!!    Modifications
!!    -------------
!     Original    May 19, 1998
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!*       0.1   declarations of arguments
!
  IMPLICIT NONE
!
  TYPE(LIST_ll), POINTER       :: TPLIST ! list of fields
  REAL, DIMENSION(:,:), TARGET :: PFIELD ! field to be added
                                         ! to the list of fields
!
!*       0.2   declarations of local variables
!
  TYPE(LIST_ll), POINTER :: TZLIST
!
!-------------------------------------------------------------------------------
!
!*       1.    In case TPLIST has not been already used
!              ----------------------------------------
  IF (.NOT. ASSOCIATED(TPLIST)) THEN ! TPLIST is empty
!
    ALLOCATE(TPLIST)
    NULLIFY(TPLIST%ARRAY1D)
    NULLIFY(TPLIST%ARRAY3D)
    TPLIST%ARRAY2D => PFIELD
    NULLIFY(TPLIST%NEXT)
    TPLIST%NCARD = 1
    TPLIST%L1D = .FALSE.
    TPLIST%L2D = .TRUE.
    TPLIST%L3D = .FALSE.
!
  ELSE 
!
!-------------------------------------------------------------------------------
!
!*       2.    TPLIST already contains fields;
!*             add the PFIELD field at the end of the list
!
    TZLIST => TPLIST
    DO WHILE (ASSOCIATED(TZLIST%NEXT))
      TZLIST => TZLIST%NEXT
    ENDDO
!
    ALLOCATE(TZLIST%NEXT)
    TZLIST => TZLIST%NEXT
    NULLIFY(TZLIST%ARRAY1D)
    NULLIFY(TZLIST%ARRAY3D)
    TZLIST%ARRAY2D => PFIELD
    TZLIST%NCARD = 0
    TZLIST%L1D = .FALSE.
    TZLIST%L2D = .TRUE.
    TZLIST%L3D = .FALSE.
    NULLIFY(TZLIST%NEXT)
!
    TPLIST%NCARD = TPLIST%NCARD + 1
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE ADD2DFIELD_ll
!
!!    ########################################
      SUBROUTINE ADD3DFIELD_ll(TPLIST, PFIELD)
!!    ########################################
!
!!****  *ADD3DFIELD_ll* -
!
!!    Purpose
!!    -------
!     This routine is used to add a 3D field (PFIELD) to a list of 
!     3D fields (TPLIST). 
!
!!    Reference
!!    ---------
!
!     User interface for the Meso-NH parallel package
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_ARGSLIST :
!         LIST_ll : list of fields
!
!!    Author
!!    ------
!
!     Ph. Kloos                 * CERFACS - CNRM *
!
!!    Modifications
!!    -------------
!     Original    May 19, 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!*       0.1   declarations of arguments
!
  IMPLICIT NONE
!
  TYPE(LIST_ll), POINTER         :: TPLIST   ! list of fields
  REAL, DIMENSION(:,:,:), TARGET :: PFIELD   ! field to be added to the list
!                                              of fields
!
!*       0.2   declarations of local variables
!
  TYPE(LIST_ll), POINTER :: TZLIST
!
!-------------------------------------------------------------------------------
!
!*       1.    In case TPLIST has not been already used
!              ----------------------------------------
!
  IF (.NOT. ASSOCIATED(TPLIST)) THEN ! TPLIST is empty
!
    ALLOCATE(TPLIST)
    NULLIFY(TPLIST%ARRAY1D)
    NULLIFY(TPLIST%ARRAY2D)
    TPLIST%ARRAY3D => PFIELD
    NULLIFY(TPLIST%NEXT)
    TPLIST%NCARD = 1
    TPLIST%L1D = .FALSE.
    TPLIST%L2D = .FALSE.
    TPLIST%L3D = .TRUE.
!
  ELSE 
!
!-------------------------------------------------------------------------------
!
!*       2.    TPLIST already contains fields;
!*             add the PFIELD field at the end of the list
!              -------------------------------------------
!
    TZLIST => TPLIST
    DO WHILE (ASSOCIATED(TZLIST%NEXT))
      TZLIST => TZLIST%NEXT
    ENDDO  
!
    ALLOCATE(TZLIST%NEXT)
    TZLIST => TZLIST%NEXT
    NULLIFY(TZLIST%ARRAY1D)
    NULLIFY(TZLIST%ARRAY2D)
    TZLIST%ARRAY3D => PFIELD
    TZLIST%NCARD = 0
    TZLIST%L1D = .FALSE.
    TZLIST%L2D = .FALSE.
    TZLIST%L3D = .TRUE.
    NULLIFY(TZLIST%NEXT)
!
    TPLIST%NCARD = TPLIST%NCARD + 1
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE ADD3DFIELD_ll
!
!!    ###############################################
      SUBROUTINE DEL1DFIELD_ll(TPLIST, PFIELD, KINFO)
!!    ###############################################
!
!!****  *DEL1DFIELD_ll* -
!
!!    Purpose
!!    -------
!     This routine deletes the PFIELD 1D array from the TPLIST list of fields
!
!!    Reference
!!    ---------
!
!     User interface for the Meso-NH parallel package
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_ARGSLIST
!        LIST1D_ll : type for a list of 1D field
!
!!    Author
!!    ------
!
!     Ph. Kloos                 * CERFACS - CNRM *
!     Didier Gazen              * LA *
!     Ronan Guivarch            * CERFACS - ENSEEIHT *
!
!!    Modifications
!!    -------------
!     Original    May 19, 1998
!     Didier Gazen (ajout DEALLOCATE)
!     Ronan Guivarch LIST_ll -> LIST1D_ll 23/02/00
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!*       0.1   declarations of arguments
!
  IMPLICIT NONE
!
  TYPE(LIST1D_ll), POINTER     :: TPLIST ! list of fields
  REAL, DIMENSION(:), TARGET   :: PFIELD ! field to be deleted from the list
                                         ! of fields
  INTEGER, INTENT(OUT)         :: KINFO  ! return status : 
                                         ! 0 if PFIELD has been found 
                                         ! in TPLIST, 1 otherwise.
! 
!*       0.2   declarations of local variables
!
  INTEGER                  :: NCARD ! number of elements in TPLIST
  TYPE(LIST1D_ll), POINTER :: TZCURRENT, TZPREV, TZTEMP 
!
!-------------------------------------------------------------------------------
!
  KINFO=1
!
!*       1.    Remove PFIELD from TPLIST 
!              -------------------------
!
  IF (ASSOCIATED(TPLIST)) THEN
    ! TPLIST is not empty : good
    NULLIFY(TZPREV)
    TZCURRENT => TPLIST
    NCARD = TPLIST%NCARD
    DO WHILE(ASSOCIATED(TZCURRENT))
      ! look for all occurrences of PFIELD in TPLIST
      !
      IF (.NOT. ASSOCIATED(TZCURRENT%ARRAY1D,PFIELD)) THEN
        ! PFIELD not found, try next element
        TZPREV => TZCURRENT
        TZCURRENT => TZCURRENT%NEXT
      ELSE
        ! PFIELD found in TZCURRENT
        TZTEMP=>TZCURRENT     
        ! Delete TZCURRENT from TPLIST
        IF (ASSOCIATED(TZPREV)) THEN
          ! not the first element of TPLIST
          TZPREV%NEXT => TZCURRENT%NEXT
        ELSE
          ! first element of TPLIST
          TPLIST => TZCURRENT%NEXT
        END IF
!
        TZCURRENT   => TZCURRENT%NEXT
        NCARD = NCARD-1
        IF (NCARD>0) TPLIST%NCARD = NCARD
!
        DEALLOCATE(TZTEMP)
        KINFO=0
      END IF
    END DO
  END IF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE DEL1DFIELD_ll
!
!!    ###############################################
      SUBROUTINE DEL2DFIELD_ll(TPLIST, PFIELD, KINFO)
!!    ###############################################
!
!!****  *DEL2DFIELD_ll* -
!
!!    Purpose
!!    -------
!     This routine deletes the PFIELD 2D array from the TPLIST list of fields
!     
!!    Reference
!!    ---------
!
!     User interface for the Meso-NH parallel package
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_ARGSLIST
!        LIST_ll : type for a list of fields
!
!!    Author
!!    ------
!
!     Ph. Kloos                 * CERFACS - CNRM *
!     Didier Gazen              * LA *
!
!!    Modifications
!!    -------------
!     Original    May 19, 1998
!     Didier Gazen (ajout DEALLOCATE)
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!*       0.1   declarations of arguments
!
  IMPLICIT NONE
!
  TYPE(LIST_ll), POINTER       :: TPLIST ! list of fields
  REAL, DIMENSION(:,:), TARGET :: PFIELD ! field to be deleted from the list
                                         ! of fields
!
  INTEGER, INTENT(OUT)       :: KINFO  ! return status : 0 if PFIELD
                                       ! has been found  in TPLIST, 1 otherwise.
!
!
!*       0.2   declarations of local variables
!
  INTEGER :: NCARD ! number of elements in TPLIST
  TYPE(LIST_ll), POINTER :: TZCURRENT, TZPREV, TZTEMP 
!
!-------------------------------------------------------------------------------
!
  KINFO=1
!
!*       1.    Remove PFIELD from TPLIST 
!              -------------------------
!
  IF (ASSOCIATED(TPLIST)) THEN
    ! TPLIST is not empty : good
    NULLIFY(TZPREV)
    TZCURRENT => TPLIST
    NCARD = TPLIST%NCARD
    DO WHILE(ASSOCIATED(TZCURRENT))
      ! look for all occurrences of PFIELD in TPLIST
      !
      IF (.NOT. ASSOCIATED(TZCURRENT%ARRAY2D,PFIELD)) THEN
        ! PFIELD not found, try next element
        TZPREV => TZCURRENT
        TZCURRENT => TZCURRENT%NEXT
      ELSE
        ! PFIELD found in TZCURRENT
        TZTEMP=>TZCURRENT
        ! Delete TZCURRENT from TPLIST
        IF (ASSOCIATED(TZPREV)) THEN
          ! not the first element of TPLIST
          TZPREV%NEXT => TZCURRENT%NEXT
        ELSE
          ! first element of TPLIST
          TPLIST => TZCURRENT%NEXT
        END IF
!      
        TZCURRENT   => TZCURRENT%NEXT
        NCARD = NCARD-1
        IF (NCARD>0) TPLIST%NCARD = NCARD
!
        DEALLOCATE(TZTEMP)
        KINFO=0
      END IF
    END DO
  END IF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE DEL2DFIELD_ll
!
!!    ###############################################
      SUBROUTINE DEL3DFIELD_ll(TPLIST, PFIELD, KINFO)
!!    ###############################################
!!
!!    Purpose
!!    -------
!     This routine deletes the PFIELD 3D array from the TPLIST list of fields 
!     
!!    Reference
!!    ---------
!
!     User interface for the Meso-NH parallel package
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_ARGSLIST
!        LIST_ll : type for a list of fields
!
!!    Author
!!    ------
!
!     Ph. Kloos                 * CERFACS - CNRM *
!     Didier Gazen              * LA *
!
!!    Modifications
!!    -------------
!     Original    May 19, 1998
!     Didier Gazen (ajout DEALLOCATE)
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!*       0.1   declarations of arguments
!
  IMPLICIT NONE
!
  TYPE(LIST_ll), POINTER         :: TPLIST ! list of fields
  REAL, DIMENSION(:,:,:), TARGET :: PFIELD ! field to be deleted
                                           ! from the list of fields
  INTEGER, INTENT(OUT)           :: KINFO  ! return status : 
                                           ! 0 if PFIELD has been found 
                                           ! in TPLIST, 1 otherwise
!
!*       0.2   declarations of local variables
!
  INTEGER :: NCARD ! number of elements in TPLIST
  TYPE(LIST_ll), POINTER :: TZCURRENT, TZPREV, TZTEMP 
!
!-------------------------------------------------------------------------------
!
  KINFO=1
!
!*       1.    Remove PFIELD from TPLIST 
!              -------------------------
!
  IF (ASSOCIATED(TPLIST)) THEN
    ! TPLIST is not empty : good
    NULLIFY(TZPREV)
    TZCURRENT => TPLIST
    NCARD = TPLIST%NCARD
    DO WHILE(ASSOCIATED(TZCURRENT))
      ! look for all occurrences of PFIELD in TPLIST
      !
      IF (.NOT. ASSOCIATED(TZCURRENT%ARRAY3D,PFIELD)) THEN
        ! PFIELD not found, try next element
        TZPREV => TZCURRENT
        TZCURRENT => TZCURRENT%NEXT
      ELSE
        ! PFIELD found in TZCURRENT
        TZTEMP=>TZCURRENT
        ! Delete TZCURRENT from TPLIST
        IF (ASSOCIATED(TZPREV)) THEN
          ! not the first element of TPLIST
          TZPREV%NEXT => TZCURRENT%NEXT
        ELSE
          ! first element of TPLIST
          TPLIST => TZCURRENT%NEXT
        END IF
!      
        TZCURRENT   => TZCURRENT%NEXT
        NCARD = NCARD-1
        IF (NCARD>0) TPLIST%NCARD = NCARD
!
        DEALLOCATE(TZTEMP)
        KINFO=0
      END IF
    END DO
  END IF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE DEL3DFIELD_ll
!
!!    ###############################
      SUBROUTINE CLEANLIST_ll(TPLIST)
!!    ###############################
!
!!****  *CLEANLIST_ll* -
!
!!    Purpose
!!    -------
!     This routine frees a TPLIST by deallocating all elements
!     created by the ADDxDFIELD routines. Arrays associated to each
!     element are untouched.
!     
!!    Reference
!!    ---------
!
!     User interface for the Meso-NH parallel package
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_ARGSLIST
!        LIST_ll : type for a list of fields
!
!!    Author
!!    ------
!
!     Didier Gazen              * LA *
!
!!    Modifications
!!    -------------
!     Original    May 19, 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
!*       0.1   declarations of arguments
!
  IMPLICIT NONE 
!
  TYPE(LIST_ll),  POINTER :: TPLIST ! List of fields
!
!*       0.2   declarations of local variables
!
  TYPE(LIST_ll), POINTER :: TZTEMP
!
!------------------------------------------------------------------------------
!
!*       1.    Dealloacte one by one the elements of TPLIST
!              --------------------------------------------
!
  IF (ASSOCIATED(TPLIST)) THEN
    DO WHILE(ASSOCIATED(TPLIST))
      TZTEMP => TPLIST
      TPLIST => TPLIST%NEXT
      DEALLOCATE(TZTEMP)
    END DO
    NULLIFY(TPLIST)
  END IF
!
!------------------------------------------------------------------------------
!
      END SUBROUTINE CLEANLIST_ll
!
END MODULE MODE_ARGSLIST_ll
