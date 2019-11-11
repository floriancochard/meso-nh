!SURFEX_LIC Copyright 1994-2014 Meteo-France 
!SURFEX_LIC This is part of the SURFEX software governed by the CeCILL-C  licence
!SURFEX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!SURFEX_LIC for details. version 1.
!     ################################################################
      SUBROUTINE UPDATE_NHALO1D( NHALO, PFIELD1D, KISIZE_ll, KJSIZE_ll, KXOR, KXEND, KYOR, KYEND, HREC )
!     ################################################################
!
!!****  *UPDATE_NHALO1D* - routine to update a halo of size NHALO for 1D field PFIELD1D 
!!                         of physical size KSIZE 
!!                         WARNING : if NHALO > KISIZE_ll or NHALO > KJSIZE_ll, this routine will crash
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	M.Moge   *LA - CNRS*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/2015 
!!        M.Moge    07/2015 bug fix : TAG computation and ICARDDIF computation 
!!                          + use a temp field PFIELDTMP for SEND_RECV_FIELD
!!        M.Moge    08/2015 calling ABORT if local subdomain is of size < NHALO
!!                          (this causes problems on the boundary of the domain)
!!        M.Moge    08/2015 bug fix : changing the computation of IISIZE
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR, ONLY : NUNDEF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODE_ll
USE MODE_EXCHANGE_ll, ONLY : SEND_RECV_FIELD
USE MODE_SPLITTING_ll, ONLY : SPLIT2
USE MODD_VAR_ll, ONLY : NPROC, IP, YSPLITTING, NMNH_COMM_WORLD
USE MODD_MPIF
USE MODD_STRUCTURE_ll, ONLY : ZONE_ll, CRSPD_ll
USE MODE_TOOLS_ll, ONLY : INTERSECTION
USE MODD_PARAMETERS, ONLY : JPHEXT, JPVEXT
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER,                      INTENT(IN)    :: NHALO      ! size of HALO to update
REAL, DIMENSION(:),       INTENT(INOUT)    :: PFIELD1D   ! field to be updated
INTEGER,                      INTENT(IN)    :: KISIZE_ll   ! global physical size of the domain in X direction
INTEGER,                      INTENT(IN)    :: KJSIZE_ll   ! global physical size of the domain in Y direction
INTEGER,                      INTENT(IN)    :: KXOR   ! origin  and 
INTEGER,                      INTENT(IN)    :: KXEND  ! end of local subdomain
INTEGER,                      INTENT(IN)    :: KYOR   ! in global domain
INTEGER,                      INTENT(IN)    :: KYEND  ! defined by KISIZE_ll and KJSIZE_ll
CHARACTER(LEN=6),       INTENT(IN) :: HREC        ! name of the parameter, 'XX' or 'YY'
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: PFIELDTMP   ! PFIELDTMP will get the communicated halos
!
!* other variables
!
INTEGER     :: JI,JJ         ! loop controls
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: IINFO_ll
!
! structures for the partitionning
!
TYPE(ZONE_ll), DIMENSION(:),ALLOCATABLE :: TZSPLITTING_PHYS  !physical splitting of the field
TYPE(ZONE_ll), DIMENSION(:),ALLOCATABLE :: TZSPLITTING_EXT   !extended splitting of the field
TYPE(ZONE_ll) :: TZFIELD_ll ! global field
!
! structures for the communications
!
TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:)  :: TZSEND, TZRECV
TYPE(CRSPD_ll), POINTER      :: TZCRSPDSEND, TZCRSPDRECV
TYPE(CRSPD_ll), ALLOCATABLE, DIMENSION(:), TARGET :: TZCRSPDSENDTAB, TZCRSPDRECVTAB
!
INTEGER :: J
INTEGER :: INBMSG
INTEGER :: ICARD
INTEGER :: ICARDDIF
INTEGER :: IISIZE, IJSIZE, IKSIZE	! local sizes of the field
INTEGER :: IXOR_ll, IYOR_ll, IZOR_ll, IXEND_ll, IYEND_ll, IZEND_ll  ! origin and end of local physical subdomain
REAL, DIMENSION(:,:,:), ALLOCATABLE    :: PFIELD3D   ! field to be updated
INTEGER , DIMENSION(NPROC) :: IXORARRAY_ALL, IYORARRAY_ALL, IXENDARRAY_ALL, IYENDARRAY_ALL  ! array containing origin and end of each subdomain in local coordinates
INTEGER , DIMENSION(NPROC) :: IDISPLS !displacement array for MPI_ALLGATHERV
INTEGER , DIMENSION(NPROC) :: IRECVCOUNTS !nteger array containing the number of elements that are to be received from each process 
!------------------------------------------------------------------------------
!
!*       1.    Coherence check
!              --------------------------------------------------------------
!
! verify that the local sizes of PFIELD3D, the halo size (NHALO) 
! and the global physical sizes (KISIZE_ll,KJSIZE_ll) are OK
!IF (  ) THEN
!
!  RETURN
!ENDIF
!!
!!	we assume that the sizes are correct
!!
ALLOCATE(TZSPLITTING_PHYS(NPROC),TZSPLITTING_EXT(NPROC))
!------------------------------------------------------------------------------
!
!*       1.    Partitionning of the field
!              --------------------------------------------------------------
!  
! Si le sous-domaine local est plus petit que NHALO, alors on aura des problemes sur les sous-domaines
! voisins des sous-domaines au bord :
!   si le sous-domaine local est a moins de NHALO points du bord mais n'est pas au bord,
!   alors il y a un probleme car on ne recevra pas de valeurs pour les points du HALO situes sur le halo externe
! Donc on fait un WARNING et un ABORT
!
IF ( NHALO > KXEND - KXOR + 1 .OR. NHALO > KYEND - KYOR + 1 ) THEN
  WRITE(*,*) "ERROR in UPDATE_NHALO1D : size of local subdomain is (", KXEND - KXOR + 1,",",KYEND - KYOR + 1, &
       ") which is less than NHALO=",NHALO
  WRITE(*,*) "Try with less MPI processes or a larger domain"
  CALL ABORT
ENDIF
!
! physical splitting of the field
!
IXOR_ll = KXOR + NHALO - JPHEXT
IXEND_ll = KXEND + NHALO - JPHEXT
IYOR_ll = KYOR + NHALO - JPHEXT
IYEND_ll = KYEND + NHALO - JPHEXT
IZOR_ll = 1
IZEND_ll = 1
!
IISIZE = IXEND_ll - IXOR_ll + 1 + 2 * NHALO
IJSIZE = IYEND_ll - IYOR_ll + 1 + 2 * NHALO
IKSIZE = 1
ALLOCATE( PFIELDTMP(IISIZE,IJSIZE,IKSIZE) )
ALLOCATE( PFIELD3D(IISIZE,IJSIZE,IKSIZE) )
IF ( HREC == 'XX' ) THEN
  DO JI = 1, IJSIZE
    PFIELD3D(:,JI,1) = PFIELD1D(:)
  ENDDO
ELSE IF ( HREC == 'YY' ) THEN
  DO JI = 1, IISIZE
    PFIELD3D(JI,:,1) = PFIELD1D(:)
  ENDDO
ENDIF
!
! we get the splitting of the physical domain with a MPI_AllGatherv
!
DO JI = 1, NPROC
  IDISPLS( JI ) = JI-1
ENDDO
IRECVCOUNTS(:) = 1
CALL MPI_ALLGATHERV( KXOR, 1, MPI_INTEGER, IXORARRAY_ALL, IRECVCOUNTS, IDISPLS, MPI_INTEGER, NMNH_COMM_WORLD, IINFO_ll)
CALL MPI_ALLGATHERV( KXEND, 1, MPI_INTEGER, IXENDARRAY_ALL, IRECVCOUNTS, IDISPLS, MPI_INTEGER, NMNH_COMM_WORLD, IINFO_ll)
CALL MPI_ALLGATHERV( KYOR, 1, MPI_INTEGER, IYORARRAY_ALL, IRECVCOUNTS, IDISPLS, MPI_INTEGER, NMNH_COMM_WORLD, IINFO_ll)
CALL MPI_ALLGATHERV( KYEND, 1, MPI_INTEGER, IYENDARRAY_ALL, IRECVCOUNTS, IDISPLS, MPI_INTEGER, NMNH_COMM_WORLD, IINFO_ll)
!
DO JI = 1, NPROC
  TZSPLITTING_PHYS(JI)%NUMBER = JI
  TZSPLITTING_PHYS(JI)%NXOR = IXORARRAY_ALL(JI) + NHALO - JPHEXT
  TZSPLITTING_PHYS(JI)%NXEND = IXENDARRAY_ALL(JI) + NHALO - JPHEXT
  TZSPLITTING_PHYS(JI)%NYOR = IYORARRAY_ALL(JI) + NHALO - JPHEXT
  TZSPLITTING_PHYS(JI)%NYEND = IYENDARRAY_ALL(JI) + NHALO - JPHEXT
  TZSPLITTING_PHYS(JI)%NZOR = 1
  TZSPLITTING_PHYS(JI)%NZEND = 1
ENDDO
!
! extended splitting of the field
!
DO JI = 1, NPROC
  TZSPLITTING_EXT(JI)%NUMBER = JI
  TZSPLITTING_EXT(JI)%NXOR = TZSPLITTING_PHYS(JI)%NXOR - NHALO
  TZSPLITTING_EXT(JI)%NXEND = TZSPLITTING_PHYS(JI)%NXEND + NHALO
  TZSPLITTING_EXT(JI)%NYOR = TZSPLITTING_PHYS(JI)%NYOR - NHALO
  TZSPLITTING_EXT(JI)%NYEND = TZSPLITTING_PHYS(JI)%NYEND + NHALO
  TZSPLITTING_EXT(JI)%NZOR = 1
  TZSPLITTING_EXT(JI)%NZEND = 1
ENDDO
!
!------------------------------------------------------------------------------
!
!*       2.    Creation of global ZONE_ll for the whole field
!              --------------------------------------------------------------
!
  TZFIELD_ll%NXOR = 1
  TZFIELD_ll%NYOR = 1
  TZFIELD_ll%NZOR = 1
    TZFIELD_ll%NXEND = KISIZE_ll + 2 * NHALO
  IF( KJSIZE_ll > 1 ) THEN
    TZFIELD_ll%NYEND = KJSIZE_ll + 2 * NHALO
  ELSE
    TZFIELD_ll%NYEND = KJSIZE_ll
  ENDIF
  TZFIELD_ll%NZEND = 1
!------------------------------------------------------------------------------
!
!*       3.    Preparing the structures for the communications
!              --------------------------------------------------------------
!
  !
  ! ######## initializing the structures for the SEND ########
  ! we take the intersection of the local physical subdomain with all extended subdomains
  !
  ALLOCATE(TZSEND(NPROC))
  CALL INTERSECTION( TZSPLITTING_EXT, NPROC, TZSPLITTING_PHYS(IP), TZSEND)
  ! il faut initialiser le TAG de manière a avoir un meme tag unique pour le send et le recv :
  !   on concatene le num du proc qui envoie et le num du proc qui recoit
  DO J = 1, NPROC
    IF ( TZSEND(J)%NUMBER > 0 ) THEN
      TZSEND(J)%MSSGTAG = IP * 10**(FLOOR(LOG10(real(TZSEND(J)%NUMBER)))+1) + TZSEND(J)%NUMBER
    ENDIF
  ENDDO
  ! switching to local coordinates
  DO J = 1, NPROC
    IF ( TZSEND(J)%NUMBER > 0 ) THEN
       TZSEND(J)%NXOR = TZSEND(J)%NXOR - IXOR_ll + NHALO + 1
       TZSEND(J)%NXEND = TZSEND(J)%NXEND - IXOR_ll + NHALO + 1
       TZSEND(J)%NYOR = TZSEND(J)%NYOR - IYOR_ll + 1
       TZSEND(J)%NYEND = TZSEND(J)%NYEND - IYOR_ll + 1
       IF( KJSIZE_ll > 1 ) THEN
         TZSEND(J)%NYOR = TZSEND(J)%NYOR + NHALO
         TZSEND(J)%NYEND = TZSEND(J)%NYEND + NHALO
       ENDIF
       TZSEND(J)%NZOR = 1
       TZSEND(J)%NZEND = 1
    ENDIF
  ENDDO
  ! we do not need the Z dimension
  DO J = 1, NPROC
    IF ( TZSEND(J)%NUMBER > 0 ) THEN
      TZSEND(J)%NZOR = 1
      TZSEND(J)%NZEND = 1
    ENDIF
  ENDDO
  ! switching from an array of CRSPD_ll to a CRSPD_ll pointer
  INBMSG = 0
  DO J = 1, NPROC
    IF ( TZSEND(J)%NUMBER > 0 ) THEN
      INBMSG = INBMSG+1
    ENDIF
  ENDDO
  ALLOCATE( TZCRSPDSENDTAB(INBMSG) )
  ICARD = 0
  ICARDDIF = 0
  DO J = 1, NPROC
    IF ( TZSEND(J)%NUMBER > 0 ) THEN
      ICARD = ICARD+1
      IF ( TZSEND(J)%NUMBER /= IP ) THEN
        ICARDDIF = ICARDDIF+1
      ENDIF
      TZCRSPDSENDTAB(ICARD)%TELT = TZSEND(J)
      IF ( ICARD == INBMSG ) THEN
        TZCRSPDSENDTAB(ICARD)%TNEXT => NULL()
      ELSE
        TZCRSPDSENDTAB(ICARD)%TNEXT => TZCRSPDSENDTAB(ICARD+1)
      ENDIF
    ENDIF
  ENDDO
  DO J = 1, ICARD
    TZCRSPDSENDTAB(J)%NCARD = ICARD-J+1
    IF ( TZCRSPDSENDTAB(J)%TELT%NUMBER > IP ) THEN
      TZCRSPDSENDTAB(J)%NCARDDIF = ICARDDIF-J+2
    ELSE
      TZCRSPDSENDTAB(J)%NCARDDIF = ICARDDIF-J+1
    ENDIF
  ENDDO
  IF (ICARD > 0) THEN
    TZCRSPDSEND => TZCRSPDSENDTAB(1)
  ELSE
    TZCRSPDSEND => NULL()
  ENDIF
  !
  ! ######## initializing the structures for the RECV ########
  ! we take the intersection of the local extend subdomain with all physical subdomains
  !
  ALLOCATE(TZRECV(NPROC))
  CALL INTERSECTION( TZSPLITTING_PHYS, NPROC, TZSPLITTING_EXT(IP), TZRECV)
  ! il faut initialiser le TAG de manière a avoir un meme tag unique pour le send et le recv :
  !   on concatene le num du proc qui envoie et le num du proc qui recoit
  DO J = 1, NPROC
    IF ( TZRECV(J)%NUMBER > 0 ) THEN
        TZRECV(J)%MSSGTAG = TZRECV(J)%NUMBER * 10**(FLOOR(LOG10(real(IP)))+1) + IP
    ENDIF
  ENDDO
  ! switching to local coordinates
  DO J = 1, NPROC
    IF ( TZRECV(J)%NUMBER > 0 ) THEN
       TZRECV(J)%NXOR = TZRECV(J)%NXOR - IXOR_ll + NHALO + 1
       TZRECV(J)%NXEND = TZRECV(J)%NXEND - IXOR_ll + NHALO + 1
       TZRECV(J)%NYOR = TZRECV(J)%NYOR - IYOR_ll + 1
       TZRECV(J)%NYEND = TZRECV(J)%NYEND - IYOR_ll + 1
      IF( KJSIZE_ll > 1 ) THEN
       TZRECV(J)%NYOR = TZRECV(J)%NYOR + NHALO
       TZRECV(J)%NYEND = TZRECV(J)%NYEND + NHALO
      ENDIF
       TZRECV(J)%NZOR = 1
       TZRECV(J)%NZEND = 1
    ENDIF
  ENDDO
  ! we do not need the Z dimension
  DO J = 1, NPROC
    IF ( TZRECV(J)%NUMBER > 0 ) THEN
      TZRECV(J)%NZOR = 1
      TZRECV(J)%NZEND = 1
    ENDIF
  ENDDO
  ! switching from an array of CRSPD_ll to a CRSPD_ll pointer
  INBMSG = 0
  DO J = 1, NPROC
    IF ( TZRECV(J)%NUMBER > 0 ) THEN
      INBMSG = INBMSG+1
    ENDIF
  ENDDO
  ALLOCATE( TZCRSPDRECVTAB(INBMSG) )
  ICARD = 0
  ICARDDIF = 0
  DO J = 1, NPROC
    IF ( TZRECV(J)%NUMBER > 0 ) THEN
      ICARD = ICARD+1
      IF ( TZRECV(J)%NUMBER /= IP ) THEN
        ICARDDIF = ICARDDIF+1
      ENDIF
      TZCRSPDRECVTAB(ICARD)%TELT = TZRECV(J)
      IF ( ICARD == INBMSG ) THEN
        TZCRSPDRECVTAB(ICARD)%TNEXT => NULL()
      ELSE
        TZCRSPDRECVTAB(ICARD)%TNEXT => TZCRSPDRECVTAB(ICARD+1)
      ENDIF
    ENDIF
  ENDDO
  DO J = 1, ICARD
    TZCRSPDRECVTAB(J)%NCARD = ICARD-J+1
    IF ( TZCRSPDRECVTAB(J)%TELT%NUMBER > IP ) THEN
      TZCRSPDRECVTAB(J)%NCARDDIF = ICARDDIF-J+2
    ELSE
      TZCRSPDRECVTAB(J)%NCARDDIF = ICARDDIF-J+1
    ENDIF
  ENDDO
  IF (ICARD > 0) THEN
    TZCRSPDRECV => TZCRSPDRECVTAB(1)
  ELSE
    TZCRSPDRECV => NULL()
  ENDIF
!
!------------------------------------------------------------------------------
!
!*       5.    Communications
!              ------------------------------------------------------------
!
PFIELDTMP(:,:,:) = PFIELD3D(:,:,:)
CALL SEND_RECV_FIELD( TZCRSPDSEND, TZCRSPDRECV, PFIELD3D, PFIELDTMP, IINFO_ll)
!
IF ( HREC == 'XX' ) THEN
  PFIELD1D(:) = PFIELDTMP(:,NHALO+1,1)
ELSE IF ( HREC == 'YY' ) THEN
  PFIELD1D(:) = PFIELDTMP(NHALO+1,:,1)
ENDIF
DEALLOCATE(PFIELD3D)
DEALLOCATE(PFIELDTMP)
DEALLOCATE(TZSEND)
DEALLOCATE(TZRECV)
DEALLOCATE(TZSPLITTING_PHYS,TZSPLITTING_EXT)
IF (LHOOK) CALL DR_HOOK('UPDATE_NHALO1D',1,ZHOOK_HANDLE)
!---------------------------------------------------------------------------
!
END SUBROUTINE UPDATE_NHALO1D
