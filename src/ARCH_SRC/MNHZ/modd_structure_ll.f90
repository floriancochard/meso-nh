!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------

!      ########################
       MODULE MODD_STRUCTURE_ll
!      ########################
!
!!****  *MODD_STRUCTURE_lll* - declaration of parallel structure
!
!!     Purpose
!!     -------
!
!      The purpose of this module is the definition of the parallel data
!      structured types.
!
!!     Reference
!!     ---------
!
!      Surcouche de parallelisation : Structure de donnees
!      R. Guivarch, Ph. Kloos, D. Lugato
!
!!     Authors
!!     -------
!
!      R. Guivarch, D. Lugato    * CERFACS - ENSEEIHT*
!      Ph. Kloos                 * CERFACS - CNRM *
!
!!    Implicit Arguments
!!    ------------------
!
!     MODD_ARGSLIST_ll
!          type LIST_ll
!     MODD_PARAMETERS_ll
!          NMAXRIM, maximum number of different RIM sizes
!
!!    Modifications
!!    -------------
!
!     Original 04/05/98
!     Juan     19/08/2005: distinction Halo NORD/SUD & EST/WEST
!
!-------------------------------------------------------------------------------
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODD_PARAMETERS_ll, ONLY : NMAXRIM
!
!-------------------------------------------------------------------------------
!
! Routines for allocation of structured types
!
INTERFACE ALLOC
  MODULE PROCEDURE LCRSPD_ALLOC, COMDATA_ALLOC, P2CDATA_ALLOC
END INTERFACE
!
!-------------------------------------------------------------------------------
!
!     ############
      TYPE ZONE_ll
!     ############
!
!!****  *Type ZONE_ll* - basic type to described a geometric domain ;
!                        it is used in the communications between the
!                        processors
!!
!!    Purpose
!!    -------
!     This type contains natural data such as a zone's number
!     to identify it, the coordinate of the origine point and
!     coordinate of the last point.
!
!     The zone's number most often corresponds to a processor number ;
!     this processor is the destination or the source of the
!     communication.
!
!     A fourth data, less natural is the Message Tag ; more than one
!     zone can correspond to the same processor (see the reference)
!     and we must separate each communication.
!
!-------------------------------------------------------------------------------
!
  SEQUENCE
!
  INTEGER :: NUMBER       ! zone's number
  INTEGER :: MSSGTAG      ! Message Tag
  INTEGER :: NXOR, NXEND  ! coordinate of the origine point
  INTEGER :: NYOR, NYEND  ! coordinate of the last point
  INTEGER :: NZOR, NZEND  !
!!$  INTEGER :: NXOR, NYOR, NZOR    ! coordinate of the origine point
!!$  INTEGER :: NXEND, NYEND, NZEND ! coordinate of the last point

!
!-------------------------------------------------------------------------------
!
      END TYPE ZONE_ll
!
!     ############
      TYPE ONEDX_ll
!     ############
!
!!****  *Type ONEDX_ll* - type used for the 1DX updates WEST/ESAT
!
!!    Purpose
!!    -------
!     This type contains two arrays for listing the neighbors
!     that the process has to communicate with : one for the
!     west neighbors and one for the east ones 
!
!-------------------------------------------------------------------------------
!
  INTEGER, DIMENSION(:), POINTER :: NSEND_WEST, NSEND_EAST
  INTEGER                        :: NBSEND_WEST, NBSEND_EAST
  INTEGER                        :: NRECV_WEST, NRECV_EAST

!
!-------------------------------------------------------------------------------
!
      END TYPE ONEDX_ll
!
!     ############
      TYPE ONEDY_ll
!     ############
!
!!****  *Type ONEDY_ll* - type used for the 1DY updates NORTH/SOUTH
!
!!    Purpose
!!    -------
!     This type contains two arrays for listing the neighbors
!     that the process has to communicate with :  one for the
!     south neighbors and one for the north ones 
!
!-------------------------------------------------------------------------------
!
  INTEGER, DIMENSION(:), POINTER :: NSEND_SOUTH, NSEND_NORTH
  INTEGER                        :: NBSEND_SOUTH, NBSEND_NORTH
  INTEGER                        :: NRECV_SOUTH, NRECV_NORTH

!
!-------------------------------------------------------------------------------
!
      END TYPE ONEDY_ll

!
!     #############
      TYPE CRSPD_ll
!     #############
!
!!****  *Type CRSPD_ll* - this type describes a list of zones
!!
!!    Purpose
!!    -------
!     This type is useful to list dynamically a set of zones ;
!     A processor will use it to list all the zones to be sent or to be
!     received
!
!-------------------------------------------------------------------------------
!
  SEQUENCE
!
  INTEGER                 :: NCARD    ! number of zones in the list
  INTEGER                 :: NCARDDIF ! number of zones in the list
                                      ! whose NUMBER field is different
                                      ! from the number of the current
                                      ! proc (IP)
!
  TYPE(ZONE_ll)           :: TELT     ! element
  TYPE(CRSPD_ll), POINTER :: TNEXT    ! pointer to the next element
!
!-------------------------------------------------------------------------------
!
      END TYPE CRSPD_ll
!
!
!     ##############
      TYPE LCRSPD_ll
!     ##############
!
!!**** *Type LCRSPD_ll* - list of CRSPD_ll
!!
!!    PURPOSE
!!    -------
!     This type implement a list which elements are of type CRSPD_ll
!     It is used to list the communication structures for the LB
!     in the grid-nesting.
!
  INTEGER                  :: NCARD  ! number of elements in the list
                                     ! (=number of different NWIDTHs)
                                     !
  INTEGER                  :: NWIDTH ! width of the zone of the coarse
                                     ! mesh used for the interpolation
                                     !
  TYPE(CRSPD_ll), POINTER  :: TCRSPD ! list of zones to be communicated
  TYPE(LIST_ll),  POINTER  :: TLIST  ! list of fields to be communicated
  TYPE(LCRSPD_ll), POINTER :: TNEXT  ! next element in the list
!
!-------------------------------------------------------------------------------
!
      END TYPE LCRSPD_ll


!     ####################
      TYPE LOCALISATION_ll
!     ####################
!
!!****  *Type LOCALISATION_ll* - localisation of a subdomain
!                                or a processor according
!                                to the boundaries of the model
!!
!!    PURPOSE
!!    -------
!     To perform communications and many specific calculation
!     we must know where the subdomain is located according to
!     the boundaries ;
!     The four logicals allow to give this information
!
!-------------------------------------------------------------------------------
!
!JUAN
  LOGICAL :: NORTH, SOUTH, EAST, WEST, TOP, BOTTOM
!JUAN
!
!-------------------------------------------------------------------------------
!
      END TYPE LOCALISATION_ll
!
!     ######################
      TYPE MODELSPLITTING_ll
!     ######################
!
!!****  *MODELSPLITTING_ll* - structure which contains for one processor
!                             its subdomain of the model
!!
!!    PURPOSE
!!    -------
!     This type contains data about the subdomain attributed to
!     a processor ; this subdomain is composed of two parts :
!     the physical domain and the extended domain which
!     is the physical domain + the halos
!
!     (P stands for physical and E for extended subdomains)
!
!-------------------------------------------------------------------------------
!
  INTEGER :: NUMBER    ! number of the subdomain
!
  INTEGER :: NXORP     ! coordinate of the origin
  INTEGER :: NXENDP    ! coordinate of the last
  INTEGER :: NDIMXP    ! x-dimension of the physical subdomain
!
  INTEGER :: NYORP     !  point of the physical subdomain
  INTEGER :: NYENDP    !  point of the physical subdomain
  INTEGER :: NDIMYP    ! y-dimension of the physical subdomain
!
  INTEGER :: NZORP     !  point of the physical subdomain
  INTEGER :: NZENDP     !  point of the physical subdomain
  INTEGER :: NDIMZP    ! z-dimension of the physical subdomain
!
  INTEGER :: NXORE     ! coordinate of the origin
  INTEGER :: NXENDE    ! coordinate of the last
  INTEGER :: NDIMXE    ! x-dimension of the extended subdomain
!
  INTEGER :: NYORE     !  point of the extended subdomain
  INTEGER :: NYENDE    !  point of the extended subdomain
  INTEGER :: NDIMYE    ! y-dimension of the extended subdomain
!
  INTEGER :: NZORE     !  point of the extended subdomain
  INTEGER :: NZENDE    !  point of the extended subdomain
  INTEGER :: NDIMZE    ! z-dimension of the extended subdomain
!
!-------------------------------------------------------------------------------
!
      END TYPE MODELSPLITTING_ll
!
!JUAN MPIALLDATA
!     ############
      TYPE BOX_ll
!     ############

  SEQUENCE
!
  INTEGER                        :: NUMBER       ! zone's number
  INTEGER                        :: MSSGTAG      ! Message Tag
  INTEGER                        :: NBOX         ! number of box/proc
  INTEGER                        :: NSIZE        ! size of buffer to allocate
  INTEGER, DIMENSION(:), POINTER :: NXOR=>NULL(), NXEND=>NULL()  ! coordinate of the origine point
  INTEGER, DIMENSION(:), POINTER :: NYOR=>NULL(), NYEND=>NULL()  ! coordinate of the last point
  INTEGER, DIMENSION(:), POINTER :: NZOR=>NULL(), NZEND=>NULL()  !
  INTEGER, DIMENSION(:), POINTER :: NCNT=>NULL()         ! size of the box
  INTEGER, DIMENSION(:), POINTER :: NSTRT=>NULL()        ! displacement in the pack buffer for mpiallallv
!
!-------------------------------------------------------------------------------
!
      END TYPE BOX_ll
!JUAN MPIALLDATA
!
!     ###############
      TYPE PROCONF_ll
!     ###############
!
!!****  *PROCONF_ll* - structure which contains all the shared
!                      informations about the knowledge of the
!                      distribution of a model upon the processors
!!
!!    PURPOSE
!!    -------
!     This type contains data to know all the decompositions of a model
!     according to different splittings ; these arrays will be allocated
!     with NPROC the number of processors.
!
!     Another information is the array of localisation
!
!     PROCONF_ll is a recursive type ; two pointers allow to go up to
!     the data of the parent model or to go down to the child models.
!
!-------------------------------------------------------------------------------
!
! number of the model
!
  INTEGER              :: NUMBER
!
! array of the 2-way splitting of the domain on the processors
!
  TYPE(MODELSPLITTING_ll), DIMENSION(:), POINTER :: TSPLITS_B
!
! array of the x-slices splitting of the domain on the processors
!
  TYPE(MODELSPLITTING_ll), DIMENSION(:), POINTER :: TSPLITS_X
!
! array of the y-slices splitting of the domain on the processors
!
  TYPE(MODELSPLITTING_ll), DIMENSION(:), POINTER :: TSPLITS_Y
!
!JUAN Z SPLITTING
  TYPE(MODELSPLITTING_ll),  DIMENSION(:),POINTER :: TSPLITS_SXP2_YP1_Z=>NULL()
  TYPE(MODELSPLITTING_ll),  DIMENSION(:),POINTER :: TSPLITS_SX_YP2_ZP1=>NULL()
  TYPE(MODELSPLITTING_ll),  DIMENSION(:),POINTER :: TSPLITS_SXP2_Y_ZP1=>NULL()
!JUAN
! array for the position of the processors according to
! the model boundaries

!
  TYPE(LOCALISATION_ll), DIMENSION(:), POINTER :: TBOUND
!
!JUAN
  TYPE(LOCALISATION_ll), DIMENSION(:), POINTER :: TBOUND_SXP2_YP1_Z=>NULL()
  TYPE(LOCALISATION_ll), DIMENSION(:), POINTER :: TBOUND_SX_YP2_ZP1=>NULL()
  TYPE(LOCALISATION_ll), DIMENSION(:), POINTER :: TBOUND_SXP2_Y_ZP1=>NULL()
!JUAN
! Pointer to the parent model (recursivity)
!
  TYPE(PROCONF_ll), POINTER :: TPARENT
!
! list of the child models (recursivity)
!
  TYPE(LPROCONF_ll), POINTER :: TCHILDREN
!
!-------------------------------------------------------------------------------
!
      END TYPE PROCONF_ll
!
!     ################
      TYPE LPROCONF_ll
!     ################
!
!!****  *LPROCONF_ll* - this type describes a list of PROCONF_ll
!!
!!    PURPOSE
!!    -------
!     This type is useful to list dynamically a set of PROCONF_ll ;
!     A model can have many child models and this set can be described
!     using this type
!
!-------------------------------------------------------------------------------
!
  INTEGER                    :: NCARD ! number of elements in the list
  TYPE(PROCONF_ll)           :: TELT
  TYPE(LPROCONF_ll), POINTER :: TNEXT
!
!-------------------------------------------------------------------------------
!
      END TYPE LPROCONF_ll
!
!     #########################
      TYPE PARENT2CHILD_DATA_ll
!     #########################
!
!!****  *PARENT2CHILD_DATA_ll*- a variable of this type contains
!                               the informations about data
!                               that a processor on the parent model
!                               sends or received from another processor
!                               on one of the child model
!
!-------------------------------------------------------------------------------
!
  INTEGER                  :: NUMBER
  TYPE(CRSPD_ll), POINTER  :: TFEEDBACK_COORD
  TYPE(LCRSPD_ll), POINTER :: TSEND_1WAY_LS, TRECV_2WAY_LS
  TYPE(LCRSPD_ll), POINTER :: TSEND_1WAY_LBXW, TSEND_1WAY_LBXE, &
                              TSEND_1WAY_LBYS, TSEND_1WAY_LBYN
!
!-------------------------------------------------------------------------------
!
      END TYPE PARENT2CHILD_DATA_ll
!
!     ##########################
      TYPE LPARENT2CHILD_DATA_ll
!     ##########################
!
!!****  *LPARENT2CHILD_DATA_ll* - this type describes
!                                 a list of PARENT2CHILD_DATA_ll
!!
!!    PURPOSE
!!    -------
!     This type is useful to list dynamically a set of PARENT2CHILD_DATA_ll ;
!     A model can have many child models and this set can be described
!     using this type
!
!-------------------------------------------------------------------------------
!
  INTEGER                              :: NCARD ! number of elements in the list
  TYPE(PARENT2CHILD_DATA_ll)           :: TELT
  TYPE(LPARENT2CHILD_DATA_ll), POINTER :: TNEXT
!
!-------------------------------------------------------------------------------
!
      END TYPE LPARENT2CHILD_DATA_ll
!
!     #####################
      TYPE PROC_COM_DATA_ll
!     #####################
!
!!****  * PROC_COM_DATA_ll* - this type contains for one processor all the data
!                             for all kind of the communications
!!
!!    PURPOSE
!!    -------
!     All processor has a variable of this type ; this variable contains
!     all the data usefull for halos and transpositions communications
!     grid nesting in future.
!
!     It's a recursive type ; two pointeurs allow to go up to the equivalent
!     data of parent model or to go down to the child models.
!
!-------------------------------------------------------------------------------
!
! Grid Nesting
!
  TYPE(PROC_COM_DATA_ll), POINTER    :: TPARENT ! pointer to the parent model
!
! number of the model
!
  INTEGER              :: NUMBER
!
! information on LB sizes
!
  INTEGER                     :: NDIMRIMLBX ! dimension of array NRIMLBX
  INTEGER, DIMENSION(NMAXRIM) :: NRIMLBX    ! array with the different
                                            ! dimensions of LBs in x
  INTEGER                     :: NDIMRIMLBY ! dimension of NRIMLBY
  INTEGER, DIMENSION(NMAXRIM) :: NRIMLBY    ! array with the different
                                            ! dimensions of LBs in y
!
! Pointers to the processor splittings in PROCONF_ll variable
! for the convenience
!
  TYPE(MODELSPLITTING_ll), POINTER :: TSPLIT_B=>NULL()
  TYPE(MODELSPLITTING_ll), POINTER :: TSPLIT_X=>NULL()
  TYPE(MODELSPLITTING_ll), POINTER :: TSPLIT_Y=>NULL()

!JUAN Z SPLITTING
  TYPE(MODELSPLITTING_ll), POINTER :: TSPLIT_SXP2_YP1_Z=>NULL()
  TYPE(MODELSPLITTING_ll), POINTER :: TSPLIT_SX_YP2_ZP1=>NULL()
  TYPE(MODELSPLITTING_ll), POINTER :: TSPLIT_SXP2_Y_ZP1=>NULL()
!JUAN
!
! subsets of correspondants for the boundaries communications
!
  TYPE(CRSPD_ll), POINTER  :: TSEND_BOUNDX, TRECV_BOUNDX, &
                              TSEND_BOUNDY, TRECV_BOUNDY, &
                              TSEND_BOUNDXY, TRECV_BOUNDXY
!
! subsets of correspondants for the halos communications
!
  TYPE(CRSPD_ll), POINTER  :: TSEND_HALO1, TRECV_HALO1, &
                              TSEND_HALO2, TRECV_HALO2
!
! subsets of correspondants for the transpositions communications
!
  TYPE(CRSPD_ll), POINTER  :: TSEND_TRANS_BX, TRECV_TRANS_BX, &
                              TSEND_TRANS_XY, TRECV_TRANS_XY
!JUAN
  TYPE(CRSPD_ll), POINTER  :: TSEND_SXP2_YP1_Z_SX_YP2_ZP1=>NULL(),&
                              TRECV_SXP2_YP1_Z_SX_YP2_ZP1=>NULL(),&
                              TSEND_SX_YP2_ZP1_SXP2_Y_ZP1=>NULL(),&
                              TRECV_SX_YP2_ZP1_SXP2_Y_ZP1=>NULL(),&
                              TSEND_SXP2_Y_ZP1_SXP2_YP1_Z=>NULL(),&
                              TRECV_SXP2_Y_ZP1_SXP2_YP1_Z=>NULL()
  TYPE(CRSPD_ll), POINTER  :: TSEND_B_SX_YP2_ZP1=>NULL(),&
                              TRECV_B_SX_YP2_ZP1=>NULL(),&
                              TSEND_SXP2_Y_ZP1_B=>NULL(),&
                              TRECV_SXP2_Y_ZP1_B=>NULL()
!JUAN MPIALLDATA
  TYPE(BOX_ll), POINTER    :: TSEND_BOX_B_SX_YP2_ZP1=>NULL(),&
                              TRECV_BOX_B_SX_YP2_ZP1=>NULL(),&
                              TSEND_BOX_SX_YP2_ZP1_SXP2_Y_ZP1=>NULL(),&
                              TRECV_BOX_SX_YP2_ZP1_SXP2_Y_ZP1=>NULL(),&
                              TSEND_BOX_SXP2_Y_ZP1_B=>NULL(),&
                              TRECV_BOX_SXP2_Y_ZP1_B=>NULL(),&
                              TSEND_BOX_SXP2_Y_ZP1_SXP2_YP1_Z=>NULL(),&
                              TRECV_BOX_SXP2_Y_ZP1_SXP2_YP1_Z=>NULL()
!JUAN MPIALLDATA
!JUAN
!
! subsets of variables for the 1D updates
!
  TYPE(ONEDX_ll), POINTER :: HALO1DX
  TYPE(ONEDY_ll), POINTER :: HALO1DY
!
! subsets of correspondants for the gridnesting communications
! (as a child)
!
  TYPE(LCRSPD_ll), POINTER :: TRECV_1WAY_LS, TSEND_2WAY_LS
  TYPE(LCRSPD_ll), POINTER :: TRECV_LBXW, TRECV_LBXE, &
                              TRECV_LBYS, TRECV_LBYN
!
  INTEGER :: NLSDIMX, NLSDIMY
!
! list of subsets for the gridnesting communications
! (as a parent)
!
  TYPE(LPARENT2CHILD_DATA_ll), POINTER       :: TP2C_DATA
  TYPE(LPROC_COM_DATA_ll), POINTER           :: TCHILDREN
!
!-------------------------------------------------------------------------------
!
      END TYPE PROC_COM_DATA_ll
!
!     ######################
      TYPE LPROC_COM_DATA_ll
!     ######################
!
!!****  *LPROC_COM_DATA_ll* - this type describes a list of PROC_COM_DATA_ll
!!
!!    PURPOSE
!!    -------
!     This type is useful to list dynamically a set of PROC_COM_DATA_ll ;
!     A model can have many child models and this set can be described
!     using this type
!
!-------------------------------------------------------------------------------
!
  INTEGER                          :: NCARD
  TYPE(PROC_COM_DATA_ll)           :: TELT
  TYPE(LPROC_COM_DATA_ll), POINTER :: TNEXT
!
!-------------------------------------------------------------------------------
!
      END TYPE LPROC_COM_DATA_ll
!
  CONTAINS
!
!     #################################
      SUBROUTINE LCRSPD_ALLOC(TPLCRSPD)
!     #################################
!
!!****  *LCRSPD_ALLOC * - allocate the memory zone for a variable
!!                        of type LCRSPD_ALLOC
!!
!!    Purpose
!!    -------
!       The purpose of this routine is allocate the memory zone for a variable
!     of type LCRSPD_ALLOC and ensure that the pointers are set to null
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(LCRSPD_ll), POINTER :: TPLCRSPD
!
!-------------------------------------------------------------------------------
!
!*       1.    ALLOCATE AND NULLIFY
!             ---------------------
!
  ALLOCATE(TPLCRSPD)
  NULLIFY(TPLCRSPD%TCRSPD)
  NULLIFY(TPLCRSPD%TLIST)
  NULLIFY(TPLCRSPD%TNEXT)
!  
!-------------------------------------------------------------------------------
!
      END SUBROUTINE LCRSPD_ALLOC
!
!     ############################
      SUBROUTINE COMDATA_ALLOC(TP)
!     ############################
!
!!****  *COMDATA_ALLOC * - allocate the memory zone for a variable
!!                         of type PROC_COM_DATA_ll
!!
!!    Purpose
!!    -------
!       The purpose of this routine is allocate the memory zone for a variable
!     of type PROC_COM_DATA_ll and ensure that the pointers are set to null and
!     that the data are equal to zero
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TP
!
!-------------------------------------------------------------------------------
!
!*       1.    ALLOCATE AND NULLIFY
!             ---------------------
!
  ALLOCATE(TP)
  NULLIFY(TP%TPARENT)
  TP%NUMBER = 0
  TP%NDIMRIMLBX = 0
  TP%NRIMLBX = 0
  TP%NDIMRIMLBY = 0
  TP%NRIMLBY = 0
  NULLIFY(TP%TSPLIT_B)
  NULLIFY(TP%TSPLIT_X)
  NULLIFY(TP%TSPLIT_Y)
!JUAN
  NULLIFY(TP%TSPLIT_SXP2_YP1_Z)
  NULLIFY(TP%TSPLIT_SX_YP2_ZP1)
  NULLIFY(TP%TSPLIT_SXP2_Y_ZP1)
!JUAN
  NULLIFY(TP%TSEND_BOUNDX)
  NULLIFY(TP%TSEND_BOUNDY)
  NULLIFY(TP%TSEND_BOUNDXY)
  NULLIFY(TP%TRECV_BOUNDX)
  NULLIFY(TP%TRECV_BOUNDY)
  NULLIFY(TP%TRECV_BOUNDXY)
  NULLIFY(TP%TSEND_HALO1)
  NULLIFY(TP%TRECV_HALO1)
  NULLIFY(TP%TSEND_HALO2)
  NULLIFY(TP%TRECV_HALO2)
  NULLIFY(TP%TSEND_TRANS_BX)
  NULLIFY(TP%TRECV_TRANS_BX)
  NULLIFY(TP%TSEND_TRANS_XY)
  NULLIFY(TP%TRECV_TRANS_XY)
  NULLIFY(TP%HALO1DX)
  NULLIFY(TP%HALO1DY)
  NULLIFY(TP%TRECV_1WAY_LS)
  NULLIFY(TP%TSEND_2WAY_LS)
  NULLIFY(TP%TRECV_LBXW)
  NULLIFY(TP%TRECV_LBXE)
  NULLIFY(TP%TRECV_LBYS)
  NULLIFY(TP%TRECV_LBYN)
  TP%NLSDIMX=0
  TP%NLSDIMY=0
  NULLIFY(TP%TP2C_DATA)
  NULLIFY(TP%TCHILDREN)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE COMDATA_ALLOC
!
!     ############################
      SUBROUTINE P2CDATA_ALLOC(TP)
!     ############################
!
!!****  *P2CDATA_ALLOC * - allocate the memory zone for a variable
!!                         of type PARENT2CHILD_DATA_ll
!!
!!    Purpose
!!    -------
!       The purpose of this routine is allocate the memory zone for a variable
!     of type PARENT2CHILD_DATA_ll and ensure that the pointers are set to null
!     and that the data are equal to zero
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(PARENT2CHILD_DATA_ll), POINTER :: TP
!
!-------------------------------------------------------------------------------
!
!*       1.    ALLOCATE AND NULLIFY
!             ---------------------
!
  ALLOCATE(TP)
  TP%NUMBER=0
  NULLIFY(TP%TSEND_1WAY_LS)
  NULLIFY(TP%TRECV_2WAY_LS)
  NULLIFY(TP%TFEEDBACK_COORD)
  NULLIFY(TP%TSEND_1WAY_LBXW)
  NULLIFY(TP%TSEND_1WAY_LBXE)
  NULLIFY(TP%TSEND_1WAY_LBYS)
  NULLIFY(TP%TSEND_1WAY_LBYN)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE P2CDATA_ALLOC

!     ############################
      SUBROUTINE ALLOC_BOX_ll(TP,K)
!     ############################
!
!!****  *P2CDATA_ALLOC * - allocate the memory zone for a variable
!!                         of type PARENT2CHILD_DATA_ll
!!
!!    Purpose
!!    -------
!       The purpose of this routine is allocate the memory zone for a variable
!     of type PARENT2CHILD_DATA_ll and ensure that the pointers are set to null
!     and that the data are equal to zero
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(BOX_ll), POINTER :: TP
  INTEGER               :: K
!
!-------------------------------------------------------------------------------
!
!*       1.    ALLOCATE AND NULLIFY
!             ---------------------
!
  ALLOCATE(TP)
  TP%NBOX = K
  ALLOCATE(TP%NXOR(K),TP%NXEND(K))
  TP%NXOR(:)  = 0 
  TP%NXEND(:) = 0 
  ALLOCATE(TP%NYOR(K),TP%NYEND(K))
  TP%NYOR(:)  = 0 
  TP%NYEND(:) = 0 
  ALLOCATE(TP%NZOR(K),TP%NZEND(K))
  TP%NZOR(:)  = 0 
  TP%NZEND(:) = 0 
  ALLOCATE(TP%NCNT(K),TP%NSTRT(K))
  TP%NCNT(:)  = 0 
  TP%NSTRT(:) = 0 

!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ALLOC_BOX_LL
!
END MODULE MODD_STRUCTURE_ll
