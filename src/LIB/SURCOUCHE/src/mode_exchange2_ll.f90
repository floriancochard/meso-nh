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
!Correction :
!  J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-----------------------------------------------------------------

!     ########################
      MODULE MODE_EXCHANGE2_ll
!     ########################
!!
!!    Purpose
!!    -------
! 
!     The purpose of this module is the implementation of communication routine
!     UPDATE_HALO2 that updates the second layer of the halo
!!
!!    Glossary
!!    --------
!       For short, we will refer as "halo2" the zone corresponding
!     to the second layer of the halo to be updated.
!!
!!    Routines Of The User Interface
!!    ------------------------------
! 
!     SUBROUTINES : UPDATE_HALO2_ll, INIT_HALO2_ll
! 
!!    Reference
!!    ---------
!
!     User Interface for Meso-NH parallel package
!     Ph. Kloos, L. Giraud, R. Guivarch, D. Lugato
!
!!    Authors
!!    -------
!
!     R. Guivarch, D. Lugato    * CERFACS *
!     Ph. Kloos                 * CERFACS - CNRM *
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_ARGSLIST_ll
!      types HALO2LIST_ll, LIST_ll
!
!     Module MODD_STRUCTURE_ll
!       types CRSPD_ll, ZONE_ll
!
!     Module MODD_VAR_ll
!       IP - Number of local processor=subdomain
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       NHALO2_COM - MPI communicator for halo 2
!       NCOMBUFFSIZE2 - buffer size
!       MPI_PRECISION - mpi precision
!       NNEXTTAG, NMAXTAG - variable to define message tag
!
!!    Modifications
!!    -------------
!       Original    May 19, 1998
!       R. Guivarch June 24, 1998 _ll
!       R. Guivarch June 29, 1998 MPI_PRECISION
!       N. Gicquel  October 30, 1998 COPY_CRSPD2
!       J.Escobar 10/02/2012 : Bug , in MPI_RECV replace 
!            MPI_STATUSES_IGNORE with MPI_STATUS_IGNORE
! 
!-------------------------------------------------------------------------------
!
!
  CONTAINS
!
!-----------------------------------------------------------------------
!
!     ######################################################################
      SUBROUTINE INIT_HALO2_ll(TPHALO2LIST, KNBVAR, KDIMX, KDIMY, KDIMZ)
!     ######################################################################
!
!!****  *INIT_HALO2_ll* initialise the second layer of the halo
!!
!!
!!    Purpose
!!    -------
!       The purpose of this routine is to allocate and initialise the 
!     TPHALO2LIST variable which contains the second layer of the
!     halo for each variable.
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_ARGSLIST_ll
!       type HALO2LIST_ll
!!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     P. Kloos                 * CERFACS - CNRM *
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
  TYPE(HALO2LIST_ll), POINTER :: TPHALO2LIST ! list of HALO2_lls
  INTEGER                     :: KNBVAR      ! number of HALO2_lls to allocate
  INTEGER        :: KDIMX, KDIMY, KDIMZ      ! dimensions of the HALO2_lls
!
!
!*       0.2   Declarations of local variables :
!
  TYPE(HALO2LIST_ll), POINTER :: TZHALO2LIST
  INTEGER :: JJ ! loop counter
!
!-------------------------------------------------------------------------------
!
!*       1.    Allocate the list of HALO2_lls
!
  ALLOCATE(TPHALO2LIST)
  TZHALO2LIST => TPHALO2LIST
!
  DO JJ=1, KNBVAR
!
!*       1.1   Allocate the current HALO2_ll
!
    ALLOCATE(TZHALO2LIST%HALO2) 
    ALLOCATE(TZHALO2LIST%HALO2%WEST(KDIMY, KDIMZ))
    ALLOCATE(TZHALO2LIST%HALO2%EAST(KDIMY, KDIMZ))
    ALLOCATE(TZHALO2LIST%HALO2%SOUTH(KDIMX, KDIMZ))
    ALLOCATE(TZHALO2LIST%HALO2%NORTH(KDIMX, KDIMZ))
    ALLOCATE(TZHALO2LIST%NEXT)
    TZHALO2LIST%HALO2%WEST=0.
    TZHALO2LIST%HALO2%EAST=0.
    TZHALO2LIST%HALO2%SOUTH=0.
    TZHALO2LIST%HALO2%NORTH=0.
!
!*       1.2   Go to the next HALO2_ll, or terminate the list
!
    IF (JJ < KNBVAR) THEN
      TZHALO2LIST => TZHALO2LIST%NEXT
    ELSE
      DEALLOCATE(TZHALO2LIST%NEXT)
      NULLIFY(TZHALO2LIST%NEXT)
    ENDIF
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE INIT_HALO2_ll
!
!     ######################################################
      SUBROUTINE UPDATE_HALO2_ll(TPLIST, TPLISTHALO2, KINFO)
!     ######################################################
!
!!****  *UPDATE_HALO2_ll* - routine to update the second layer halo
!!
!!    Purpose
!!    -------
!       This routine updates the halo with the values computed by the 
!     neighbor subdomains. The fields to be updated are in the
!     TPLIST list of fields. Before UPDATE_HALO is called, TPLIST
!     has been filled with the fields to be communicated
!
!!**  Method
!!    ------
!       First the processors send their internal halos to their
!     neighboring processors/subdomains, and then they receive
!     their external halos from their neighboring processors.
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_ARGSLIST_ll
!       type HALO2LIST_ll, LIST_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       NHALO2_COM - MPI communicator for halo 2
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     P. Kloos                 * CERFACS - CNRM *
!!    J.Escobar 21/03/014: add mppd_check for all updated field
! 
!-------------------------------------------------------------------------------!
!*       0.    DECLARATIONS
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll, HALO2LIST_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, NHALO2_COM
!
  USE MODE_MPPDB
!
!*       0.1   declarations of arguments
!
  TYPE(LIST_ll), POINTER      :: TPLIST      ! pointer to the list of 
                                             ! fields to be sent
  TYPE(HALO2LIST_ll), POINTER :: TPLISTHALO2 ! pointer to the list of
                                             ! halo2 to be received
  INTEGER                     :: KINFO       ! return status
!
  TYPE(LIST_ll), POINTER :: TZFIELD
!
  INTEGER                :: ICOUNT
  CHARACTER*2            :: YCOUNT
!
!-------------------------------------------------------------------------------
!
!*       1.   SEND / RECV THE INTERNAL HALO TO/FROM THE NEIGHBORING PROCESSORS
!             ----------------------------------------------------
!
  CALL SEND_RECV_CRSPD2(TCRRT_COMDATA%TSEND_HALO2, TCRRT_COMDATA%TRECV_HALO2, &
                        TPLIST, TPLISTHALO2, NHALO2_COM, KINFO)

!------------------------------------------------------------------------------
!
!*       2.  ZONES TO SEND TO THE PROC ITSELF
!            --------------------------------
!
   CALL COPY_CRSPD2(TCRRT_COMDATA%TSEND_HALO2, TCRRT_COMDATA%TRECV_HALO2, &
                    TPLIST, TPLISTHALO2, KINFO)
!
!JUAN MPP_CHECK2D/3D
!
   IF (MPPDB_INITIALIZED) THEN
      TZFIELD => TPLIST
      ICOUNT=0
      DO WHILE (ASSOCIATED(TZFIELD))
         ICOUNT=ICOUNT+1
         WRITE(YCOUNT,'(I2)') ICOUNT
         IF (TZFIELD%L2D) THEN
            CALL MPPDB_CHECK2D(TZFIELD%ARRAY2D,"UPDATE_HALO2_ll::TAB2D("//YCOUNT//")",PRECISION)
         ELSEIF(TZFIELD%L3D) THEN
            CALL MPPDB_CHECK3D(TZFIELD%ARRAY3D,"UPDATE_HALO2_ll::TAB2D("//YCOUNT//")",PRECISION)
         END IF
         TZFIELD => TZFIELD%NEXT
      END DO
   END IF
!
!----------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_HALO2_ll
!
!     ##########################################################################
      SUBROUTINE COPY_CRSPD2(TPSENDCRSPD, TPRECVCRSPD, TPSENDLIST, TPRECVLIST, &
                             KINFO)
!     ##########################################################################
!
!!****  *COPY_CRSPD2* - routine to copy zones that a proc sends to itself
!!
!!    Purpose
!!    -------
!!  
!      Copy the field sendtplist of the zone sendcrspd to the field recvtplist
!      of the zones recvcrspd.    
!
!!**  Method
!!    ------
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       type CRSPD_ll
!
!     Module MODD_ARGSLIST_ll
!       type HALO2LIST_ll, LIST_ll
!
!     Module MODD_VAR_ll
!       IP - Number of local processor=subdomain
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     N. Gicquel               * CERFACS - CNRM *
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : CRSPD_ll
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll, HALO2LIST_ll
  USE MODD_VAR_ll, ONLY : IP
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(CRSPD_ll), POINTER      :: TPSENDCRSPD, TPRECVCRSPD
  TYPE(LIST_ll), POINTER       :: TPSENDLIST
  TYPE(HALO2LIST_ll), POINTER  :: TPRECVLIST
  INTEGER                      :: KINFO
!
!*       0.2   declarations of local variables
!
  TYPE(CRSPD_ll), POINTER      :: TZSEND, TZRECV
!
!------------------------------------------------------------------------------
!
  TZSEND => TPSENDCRSPD
  DO WHILE (ASSOCIATED(TZSEND))
    IF (TZSEND%TELT%NUMBER == IP) THEN
      TZRECV => TPRECVCRSPD 
      DO WHILE (ASSOCIATED(TZRECV))
        IF (TZRECV%TELT%NUMBER == IP&
             & .AND.TZSEND%TELT%MSSGTAG == TZRECV%TELT%MSSGTAG) THEN
          CALL COPY_ZONE2(TZSEND%TELT, TZRECV%TELT, TPSENDLIST, &
                          TPRECVLIST, KINFO)
        ENDIF
        TZRECV => TZRECV%TNEXT
      ENDDO
    ENDIF
    TZSEND => TZSEND%TNEXT
  ENDDO
!
!------------------------------------------------------------------------------
!
      END SUBROUTINE COPY_CRSPD2
!
!     ####################################################################
      SUBROUTINE COPY_ZONE2(TPSEND, TPRECV, TPSENDLIST, TPRECVLIST, KINFO)
!     ####################################################################
!
!!****  *COPY_ZONE2* - 
!!
!!    Purpose
!!    -------
!       This routine copies the values of the fields in the TPSENDLIST to the
!       halo 2 fields of TPRECVLIST according the TPSEND and TPRECV zones
!
!!**  Method
!!    ------
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       type ZONE_ll
!
!     Module MODD_ARGSLIST_ll
!       type HALO2LIST_ll, LIST_ll
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     P. Kloos                 * CERFACS - CNRM *
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll, HALO2LIST_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(ZONE_ll)                :: TPSEND, TPRECV 
  TYPE(LIST_ll), POINTER       :: TPSENDLIST
  TYPE(HALO2LIST_ll), POINTER  :: TPRECVLIST
  INTEGER                      :: KINFO
!
!*       0.2   declarations of local variables
!
  TYPE(LIST_ll), POINTER       :: TZLIST
  TYPE(HALO2LIST_ll), POINTER  :: TZHALO
  REAL, DIMENSION(:,:), POINTER :: TZTAB2D 
  INTEGER :: IIBS, IIES, IJBS, IJES, IIBR, IIER, IJBR, IJER, IKES, IKBS, &
             IKBR, IKER
!
  INTEGER, PARAMETER :: ISENDNORTH=2, &
                        ISENDWEST=4, &
                        ISENDSOUTH=6, &
                        ISENDEAST=8, &
                        ISENDCYCNORTH=12, &
                        ISENDCYCWEST=14, &
                        ISENDCYCSOUTH=16, &
                        ISENDCYCEAST=18
!
!-------------------------------------------------------------------------------
!
  IIBS = TPSEND%NXOR
  IIES = TPSEND%NXEND
  IJBS = TPSEND%NYOR
  IJES = TPSEND%NYEND
  IKBS = TPSEND%NZOR
  IKES = TPSEND%NZEND
!
  IIBR = TPRECV%NXOR
  IIER = TPRECV%NXEND
  IJBR = TPRECV%NYOR
  IJER = TPRECV%NYEND
  IKBR = TPRECV%NZOR
  IKER = TPRECV%NZEND
!
  TZLIST => TPSENDLIST
  TZHALO => TPRECVLIST 
  DO WHILE (ASSOCIATED(TZLIST)) 
    SELECT CASE(TPSEND%MSSGTAG)
      CASE(ISENDNORTH, ISENDCYCNORTH)
        TZTAB2D => TZHALO%HALO2%SOUTH
        TZTAB2D(IIBR:IIER, IKBR:IKER) = &
                                     TZLIST%ARRAY3D(IIBS:IIES, IJBS, IKBS:IKES)
      CASE(ISENDSOUTH, ISENDCYCSOUTH)
        TZTAB2D => TZHALO%HALO2%NORTH
        TZTAB2D(IIBR:IIER, IKBR:IKER) = &
                                     TZLIST%ARRAY3D(IIBS:IIES, IJBS, IKBS:IKES)
      CASE(ISENDWEST, ISENDCYCWEST)
        TZTAB2D => TZHALO%HALO2%EAST
        TZTAB2D(IJBR:IJER, IKBR:IKER) = &
                                     TZLIST%ARRAY3D(IIBS, IJBS:IJES, IKBS:IKES)
      CASE(ISENDEAST, ISENDCYCEAST)
        TZTAB2D => TZHALO%HALO2%WEST
        TZTAB2D(IJBR:IJER, IKBR:IKER) = &
                                     TZLIST%ARRAY3D(IIBS, IJBS:IJES, IKBS:IKES)
    END SELECT    
    TZLIST => TZLIST%NEXT
    TZHALO => TZHALO%NEXT
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE COPY_ZONE2
!
!     ######################################################
      SUBROUTINE FILLOUT_ZONE2(TPHALO2LIST, TPZONE, PBUFFER)
!     ######################################################
!
!!****  *FILLOUT_ZONE2* -
!!
!!    Purpose
!!    -------
!       This routine receives the data of the fields of the TPHALO2LIST
!     list of fields situated in the TPZONE ZONE_ll.
!
!!**  Method
!!    ------
!       First the data are received in a buffer. Then each field
!     of the TPHALO2LIST list of fields is filled in at the
!     location pointed by the zone.
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       type ZONE_ll
!
!     Module MODD_ARGSLIST_ll
!       type HALO2LIST_ll
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     P. Kloos                 * CERFACS - CNRM *
!
!-------------------------------------------------------------------------------
!
!
!*       0.    DECLARATIONS
!              ------------
!
  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
  USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
  USE MODD_MPIF
!
  IMPLICIT NONE

!  INCLUDE 'mpif.h'
!
!*       0.1   declarations of arguments
!
  TYPE(ZONE_ll)          :: TPZONE      ! ZONE_ll to be received
  TYPE(HALO2LIST_ll), POINTER :: TPHALO2LIST ! list of halo2 to be received
  REAL, DIMENSION(:) :: PBUFFER     ! reception buffer for unpacking data
!
!  INTEGER, DIMENSION(MPI_STATUS_SIZE) ::  KSTATUS ! status of received message
!
!*       0.2   declarations of local variables
!
  INTEGER :: JI,JJ,JK,JINC           ! loop counters
  TYPE(HALO2LIST_ll), POINTER :: TZHALO2  ! temporary list of halo2
!
  INTEGER, PARAMETER :: ISENDNORTH=2, &
             ISENDWEST=4, &
             ISENDSOUTH=6, &
             ISENDEAST=8, &
             ISENDCYCNORTH=12, &
             ISENDCYCWEST=14, & 
             ISENDCYCSOUTH=16, &
             ISENDCYCEAST=18
!
!-------------------------------------------------------------------------------
!
!*       1.    Set MPI message tags
!              --------------------
!  ISENDNORTH = 2
!  ISENDWEST  = 4
!  ISENDSOUTH = 6
!  ISENDEAST  = 8 
!
!  ISENDCYCNORTH = 12 ! !  ISENDCYCWEST  = 14 !  | In case of cyclic 

!  ISENDCYCSOUTH = 16 !  | boundary conditions
!  ISENDCYCEAST  = 18 ! /
 
!
!*       2.    Store the received message
!              --------------------------
!
!*       2.1   See at which side of the halo the message has to be stored
!
  SELECT CASE(TPZONE%MSSGTAG)
!
!*       2.1.1 Message is to be put in the south halo
!
  CASE(ISENDNORTH, ISENDCYCNORTH) 
!
!*       2.1.1.1   Go over the TPHALO2LIST list of halo2
!
    TZHALO2 => TPHALO2LIST
    JINC=0
    DO WHILE (ASSOCIATED(TZHALO2))
!
!*       2.1.1.2   Fill out the buffer in the south part of the halo2
!
      DO JK=TPZONE%NZOR, TPZONE%NZEND
          DO JI=TPZONE%NXOR, TPZONE%NXEND
!
            JINC = JINC + 1
            TZHALO2%HALO2%SOUTH(JI,JK) = PBUFFER(JINC)
!
          ENDDO
      ENDDO
!
!*       2.1.1.3   Go to the next halo2 in the list
!
      TZHALO2 => TZHALO2%NEXT
!
    ENDDO
!
!*       2.1.2 Message is coming from the east
!
  CASE(ISENDWEST, ISENDCYCWEST)
!
!*       2.1.2.1   Go over the TPHALO2LIST list of halo2
!
    TZHALO2 => TPHALO2LIST
    JINC=0
    DO WHILE (ASSOCIATED(TZHALO2))
!
!*       2.1.2.2   Fill out the buffer in the east part of the halo2
!
      DO JK=TPZONE%NZOR, TPZONE%NZEND
        DO JJ=TPZONE%NYOR, TPZONE%NYEND
!
            JINC = JINC + 1
            TZHALO2%HALO2%EAST(JJ,JK) = PBUFFER(JINC)
!
        ENDDO
      ENDDO
!
!*       2.1.2.3   Go to the next halo2 in the list
!
      TZHALO2 => TZHALO2%NEXT
!
    ENDDO
!
!*       2.1.3 Message is coming from the north
!
  CASE(ISENDSOUTH, ISENDCYCSOUTH)
!
!*       2.1.3.1   Go over the TPHALO2LIST list of halo2
!
    TZHALO2 => TPHALO2LIST
    JINC=0
    DO WHILE (ASSOCIATED(TZHALO2))
!
!*       2.1.3.2   Fill out the buffer in the north part of the halo2
!
      DO JK=TPZONE%NZOR, TPZONE%NZEND
          DO JI=TPZONE%NXOR, TPZONE%NXEND
!
            JINC = JINC + 1
            TZHALO2%HALO2%NORTH(JI,JK) = PBUFFER(JINC)
!
          ENDDO
      ENDDO
!
!*       2.1.3.3   Go to the next halo2 in the list
!
      TZHALO2 => TZHALO2%NEXT
!
    ENDDO
!
!*       2.1.4 Message is coming from the west
!
  CASE(ISENDEAST, ISENDCYCEAST)
!
!*       2.1.4.1   Go over the TPHALO2LIST list of halo2
!
    TZHALO2 => TPHALO2LIST
    JINC=0
    DO WHILE (ASSOCIATED(TZHALO2))
!
!*       2.1.4.2   Fill out the buffer in the west part of the halo2
!
      DO JK=TPZONE%NZOR, TPZONE%NZEND
        DO JJ=TPZONE%NYOR, TPZONE%NYEND
!
            JINC = JINC + 1
            TZHALO2%HALO2%WEST(JJ,JK) = PBUFFER(JINC)
!
        ENDDO
      ENDDO
!
!*       2.1.4.3   Go to the next halo2 in the list
!
      TZHALO2 => TZHALO2%NEXT
!
    ENDDO
!
  END SELECT
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE FILLOUT_ZONE2
!
!     ########################################################################
      SUBROUTINE SEND_RECV_CRSPD2(TPCRSPDSEND, TPCRSPDRECV, TPFIELDLISTSEND, &
                                  TPFIELDLISTRECV, KMPI_COMM, KINFO)
!     ########################################################################
!
!!****  *SEND_RECV_CRSPD2* -
!
!!    Purpose
!!    -------
!       This routine sends the data of the TPFIELDLISTSEND list 
!       to the correspondants of the TPCRSPDSEND list 
!       and receives the data of the TPFIELDLISTRECV list
!       from the correspondants  of the TPCRSPDRECV list
!
!!**  Method
!!    ------
!
!     This routine is based on the following rule : 
!        one sent for one received (if it is possible);
!     The algorithm is the following :
!
!     while (there is some messages to send)
!       OR  (there is some messages to receive)
!     do
!       if (there is some messages to send)
!         send the first of the list
!       end if
!       if (there is some messages to receive)
!         try to receive
!       end if
!     done
!
!     The receptions are non-blocking and don't follow necessarly
!     the order of the TPCRSPDRECV list
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       type CRSPD_ll, ZONE_ll
!
!     Module MODD_ARGSLIST_ll
!       type LIST_ll, HALO2LIST_ll
!
!     Module MODD_VAR_ll
!       IP - Number of local processor=subdomain
!       NCOMBUFFSIZE2 - buffer size
!       MPI_PRECISION - mpi precision
!       NNEXTTAG, NMAXTAG - variable to define message tag
!               
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     P. Kloos                 * CERFACS - CNRM *
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : CRSPD_ll, ZONE_ll
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll, HALO2LIST_ll
!
  USE MODD_VAR_ll, ONLY : NCOMBUFFSIZE2, IP, MPI_PRECISION, NNEXTTAG, NMAXTAG
  USE MODE_EXCHANGE_ll, ONLY : FILLIN_BUFFERS
  USE MODD_MPIF
!JUANZ
  USE MODD_CONFZ, ONLY : LMNH_MPI_BSEND
!JUANZ
!
  USE MODD_VAR_ll, ONLY : MNH_STATUSES_IGNORE
!
  IMPLICIT NONE
!
!  INCLUDE 'mpif.h'
!
!*       0.1   declarations of arguments
!
  TYPE(CRSPD_ll), POINTER :: TPCRSPDSEND, TPCRSPDRECV
  TYPE(LIST_ll), POINTER :: TPFIELDLISTSEND
  TYPE(HALO2LIST_ll), POINTER :: TPFIELDLISTRECV 
  INTEGER :: KMPI_COMM
  INTEGER :: KINFO
!
!*       0.2   declarations of local variables
!
  INTEGER :: JINC, JI, JJ, JK  ! Loop and counter variables
  INTEGER :: FOUND, KERROR
!
!JUAN
!if defined (MNH_MPI_ISEND)
  REAL, DIMENSION (:,:), ALLOCATABLE,TARGET :: TZBUFFER ! Buffers for info
                                                      ! received
!!$#else
!!$  REAL, DIMENSION (NCOMBUFFSIZE2), TARGET :: TZBUFFER ! Buffers for info
!!$                                                      ! received
!!$#endif
!JUAN
!
 ! INTEGER IRECVSTATUS(MPI_STATUS_SIZE) ! Status of completed receive request
  LOGICAL :: IRECVFLAG, ISENDFLAG
  INTEGER :: IMSGTAG, ISENDERPROC
!
  INTEGER :: IRECVNB, ISENDNB ! Total numbers of receive and send to do
  INTEGER :: IRECVDONE ! RECEIVE COMPLETED (receive and treated)
!
  TYPE(CRSPD_ll), POINTER :: TPMAILSEND, TPMAILRECV
  TYPE(ZONE_ll), POINTER :: TZZONESEND
!
  TYPE(LIST_ll), POINTER :: TZFIELDLISTSEND
  TYPE(HALO2LIST_ll), POINTER :: TZFIELDLISTRECV
  INTEGER, SAVE :: ITAGOFFSET = 0
! JUAN
!if defined (MNH_MPI_ISEND)
INTEGER,PARAMETER                                     :: MPI_MAX_REQ = 1024
!INTEGER,SAVE,DIMENSION(MPI_MAX_REQ)                  :: REQ_TAB
!INTEGER,SAVE,DIMENSION(MPI_STATUS_SIZE,MPI_MAX_REQ)  :: STATUS_TAB
INTEGER,ALLOCATABLE,DIMENSION(:)                      :: REQ_TAB
INTEGER                                               :: NB_REQ,NFIRST_REQ_RECV
!endif
! JUAN
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISATIONS
!              ---------------
!
  TZFIELDLISTSEND => TPFIELDLISTSEND
  TZFIELDLISTRECV => TPFIELDLISTRECV
!
  IRECVDONE = 0
!
  IF (ASSOCIATED(TPCRSPDRECV)) THEN
    IRECVNB = TPCRSPDRECV%NCARDDIF
  ELSE
    IRECVNB = 0
  ENDIF
  IF (ASSOCIATED(TPCRSPDSEND)) THEN
    ISENDNB = TPCRSPDSEND%NCARDDIF
  ELSE
    ISENDNB = 0
  ENDIF

!JUAN
!if defined (MNH_MPI_ISEND)
  IF (LMNH_MPI_BSEND) THEN
     ALLOCATE(TZBUFFER(NCOMBUFFSIZE2,1))
  ELSE
     ALLOCATE(TZBUFFER(NCOMBUFFSIZE2,ISENDNB+IRECVNB))
     NB_REQ = 0
     ALLOCATE(REQ_TAB(ISENDNB+IRECVNB))
  END IF
!endif
!JUAN
! 
  TPMAILRECV => TPCRSPDRECV
  TPMAILSEND => TPCRSPDSEND
!
!NZJUAN  CALL MPI_BARRIER(KMPI_COMM, KERROR)
!
!-------------------------------------------------------------------------------
!
!*       2.    MAIN LOOP
!              ---------
!
  DO WHILE (ASSOCIATED(TPMAILSEND))
!
!*       2.1  if there is still something to send
!
     IF (ASSOCIATED(TPMAILSEND)) THEN
        TZZONESEND => TPMAILSEND%TELT
        IF (TZZONESEND%NUMBER /= IP) THEN 
           JINC = 0
!JUAN
!if defined(MNH_MPI_ISEND)
 IF ( .NOT. LMNH_MPI_BSEND) THEN
           NB_REQ = NB_REQ + 1
           CALL FILLIN_BUFFERS(TZFIELDLISTSEND, TZZONESEND, TZBUFFER(:,NB_REQ), JINC)
 else
           CALL FILLIN_BUFFERS(TZFIELDLISTSEND, TZZONESEND, TZBUFFER(:,1), JINC)
 endif
!JUAN
!if defined(MNH_MPI_BSEND)
 IF (LMNH_MPI_BSEND) THEN
           CALL MPI_BSEND(TZBUFFER, JINC, MPI_PRECISION, &
                TZZONESEND%NUMBER - 1, TZZONESEND%MSSGTAG + ITAGOFFSET, &
                KMPI_COMM, KERROR)
else
!JUAN
!if defined(MNH_MPI_ISEND)
           CALL MPI_ISEND(TZBUFFER(1,NB_REQ), JINC, MPI_PRECISION, &
                TZZONESEND%NUMBER - 1, TZZONESEND%MSSGTAG + ITAGOFFSET, &
                KMPI_COMM, REQ_TAB(NB_REQ), KERROR)

 endif
!JUAN 

        ENDIF
        TPMAILSEND => TPMAILSEND%TNEXT
     ENDIF
  ENDDO

!NZJUAN  CALL MPI_BARRIER(KMPI_COMM, KERROR)

! JUAN
!if defined (MNH_MPI_ISEND)
  IF ( .NOT. LMNH_MPI_BSEND) THEN
     NFIRST_REQ_RECV = NB_REQ
  endif
  ! JUAN

  DO WHILE (ASSOCIATED(TPMAILRECV)) 
     IF (TPMAILRECV%TELT%NUMBER  == IP) THEN
        TPMAILRECV => TPMAILRECV%TNEXT
     ELSE
! JUAN
!if defined (MNH_MPI_ISEND)
 IF ( .NOT. LMNH_MPI_BSEND) THEN
        NB_REQ = NB_REQ + 1
        CALL MPI_IRECV(TZBUFFER(1,NB_REQ), NCOMBUFFSIZE2, MPI_PRECISION, &
             TPMAILRECV%TELT%NUMBER-1, &
             TPMAILRECV%TELT%MSSGTAG + ITAGOFFSET, &
             KMPI_COMM, REQ_TAB(NB_REQ), KERROR)
else
        CALL MPI_RECV(TZBUFFER, NCOMBUFFSIZE2, MPI_PRECISION, &
             TPMAILRECV%TELT%NUMBER-1, &
             TPMAILRECV%TELT%MSSGTAG + ITAGOFFSET, &
             KMPI_COMM, MPI_STATUS_IGNORE, KERROR)
        CALL FILLOUT_ZONE2(TZFIELDLISTRECV, TPMAILRECV%TELT, TZBUFFER(:,1))
endif
! JUAN
        TPMAILRECV => TPMAILRECV%TNEXT
     ENDIF
     
  ENDDO

! JUAN
!if defined (MNH_MPI_ISEND)
 IF ( .NOT. LMNH_MPI_BSEND) THEN
    CALL MPI_WAITALL(NB_REQ,REQ_TAB,MNH_STATUSES_IGNORE,KINFO) 
    
    TPMAILRECV => TPCRSPDRECV
    NB_REQ = NFIRST_REQ_RECV
    
    DO WHILE (ASSOCIATED(TPMAILRECV)) 
       IF (TPMAILRECV%TELT%NUMBER  == IP) THEN
          TPMAILRECV => TPMAILRECV%TNEXT
       ELSE
          NB_REQ = NB_REQ + 1
          CALL FILLOUT_ZONE2(TZFIELDLISTRECV, TPMAILRECV%TELT, TZBUFFER(:,NB_REQ))
          TPMAILRECV => TPMAILRECV%TNEXT
       ENDIF
       
    ENDDO
    
    DEALLOCATE(REQ_TAB)
 endif
!JUAN 
! 
 DEALLOCATE(TZBUFFER)
  ITAGOFFSET = MOD((ITAGOFFSET + NNEXTTAG), NMAXTAG)
!
!NZJUAN  CALL MPI_BARRIER(KMPI_COMM, KERROR)
!
!-------------------------------------------------------------------------------!
      END SUBROUTINE SEND_RECV_CRSPD2
!
END MODULE MODE_EXCHANGE2_ll 
