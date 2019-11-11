/*
 *  $Id$
 *
 *  (C) 1993 by Argonne National Laboratory and Mississipi State University.
 *      All rights reserved.  See COPYRIGHT in top-level directory.
 */

/* Atomic types */

#ifndef _DMPIATOM_INCLUDE
#define _DMPIATOM_INCLUDE

#include "mpi.h"

/*****************************************************************************
*                           MPI ATOMIC DATA STRUCTURES                       *
*****************************************************************************/

/*****************************************************************************
*  We place "cookies" into the data structures to improve error detection and*
*  reporting of invalid objects.  In order to make this flexible, the        *
*  cookies are defined as macros.                                            *
*  If MPIR_HAS_COOKIES is not defined, then the "cookie" fields are not      *
*  set or tested                                                             *
*****************************************************************************/

#define MPIR_HAS_COOKIES
#ifdef MPIR_HAS_COOKIES
#define MPIR_COOKIE unsigned long cookie;
#define MPIR_SET_COOKIE(obj,value) (obj)->cookie = (value);
#else
#define MPIR_COOKIE
#define MPIR_SET_COOKIE(obj,value)
#endif
/****************************************************************************/

/* Ensure that this is at least 32 bits */
typedef unsigned long MPIR_uint32;

typedef struct _MPIR_HBT_node {
  MPIR_COOKIE                    /* Cookie to help detect valid items */
  void *value;
  int   keyval;
  short balance;
  struct _MPIR_HBT_node *left;
  struct _MPIR_HBT_node *right;
} MPIR_HBT_node;
#define MPIR_HBT_NODE_COOKIE 0x03b740de

typedef struct _MPIR_HBT {
  MPIR_COOKIE                    /* Cookie to help detect valid items */
  unsigned int   height;
  int            ref_count;
  MPIR_HBT_node *root;
} *MPIR_HBT;
#define MPIR_HBT_COOKIE 0x03b7c007

/*
 * In order to provide better compile time error checking for the 
 * implementation, we use a union to store the actual copy/delete functions
 * for the different languages
 */
#ifndef ANSI_ARGS
#if defined(__STDC__) || defined(__cplusplus) || defined(HAVE_PROTOTYPES)
#define ANSI_ARGS(a) a
#else
#define ANSI_ARGS(a) ()
#endif
#endif

typedef union {
    int (*c_copy_fn) ANSI_ARGS(( MPI_Comm, int, void *, void *, void *, 
				 int * ));
    void (*f77_copy_fn) ANSI_ARGS(( MPI_Comm, int *, int *, int *, int *, 
				   int *, int * ));
} MPIR_Copy_fn;
typedef union {
    int (*c_delete_fn) ANSI_ARGS(( MPI_Comm, int, void *, void * ));
    void (*f77_delete_fn) ANSI_ARGS(( MPI_Comm, int *, int *, void *, int * ));
} MPIR_Delete_fn;

typedef struct  {
    MPIR_COOKIE                    /* Cookie to help detect valid items */
    MPIR_Copy_fn copy_fn;
    MPIR_Delete_fn delete_fn;
    void  *extra_state;
    int    FortranCalling;        /* Used to indicate whether Fortran or
				     C calling conventions are used for
				     copy_fn (attribute_in is passed by 
				     value in C, must be passed by reference
				     in Fortran); the underlying code
				     must also understand what a 
				     Fortran logical looks like */
    int    ref_count;
    int    permanent;             /* Used to mark the permanent attributes of
				     MPI_COMM_WORLD */
} MPIR_Attr_key;
#define MPIR_ATTR_COOKIE 0xa774c003

#define MPIR_GROUP_COOKIE 0xea01beaf
struct MPIR_GROUP {
    MPIR_COOKIE             /* Cookie to help detect valid item */
    int np;			        /* Number of processes in group */
    int local_rank;         /* My rank in the group (if I belong) */
    int ref_count;          /* Number of references to this group */
    int N2_next;            /* Next power of 2 from np */
    int N2_prev;            /* Previous power of 2 from np */
    int permanent;          /* Permanent group */
    int *lrank_to_grank;    /* Mapping from local to "global" ranks */
    int *set_mark;          /* Used for set marking/manipulation on groups */
};

/* 
   Error handlers must survive being deleted and set to MPI_ERRHANDLER_NULL,
   the reference count is for knowing how many communicators still have this
   error handler active 
 */
struct MPIR_Errhandler {
    MPIR_COOKIE                    /* Cookie to help detect valid items */
    MPI_Handler_function *routine;
    int                  ref_count;
    };
#define MPIR_ERRHANDLER_COOKIE 0xe443a2dd


typedef unsigned long MPIR_CONTEXT;
#define  MPIR_CONTEXT_TYPE MPI_UNSIGNED_LONG

#define  MPIR_WORLD_PT2PT_CONTEXT 0
#define  MPIR_WORLD_COLL_CONTEXT  1
#define  MPIR_SELF_PT2PT_CONTEXT  2
#define  MPIR_SELF_COLL_CONTEXT   3
#define  MPIR_FIRST_FREE_CONTEXT  4

typedef enum { MPIR_INTRA=1, MPIR_INTER } MPIR_COMM_TYPE;

#ifdef ANSI_ARGS
#undef ANSI_ARGS
#endif

#if defined(__STDC__) || defined(__cplusplus) || defined(HAVE_PROTOTYPES)
#define ANSI_ARGS(a) a
#else
#define ANSI_ARGS(a) ()
#endif

typedef struct _MPIR_COLLOPS *MPIR_COLLOPS;

/*
   The local_rank field is used to reduce unnecessary memory references
   when doing send/receives.  It must equal local_group->local_rank.

   lrank_to_grank is group->lrank_to_grank; this is also used to 
   reduce memory refs.  (it is IDENTICAL, not just a copy; the "group"
   owns the array.)

   These have been ordered so that the most common elements are 
   near the top, in hopes of improving cache utilization.

   For a normal intra-communicator the group and local_group are identical
   The group differs from the local_group only in an inter-communicator
 */
#define MPIR_COMM_COOKIE 0xea02beaf
struct MPIR_COMMUNICATOR {
    MPIR_COOKIE                   /* Cookie to help detect valid item */
    /* Most common data from group is cached here */
    int           np;             /* size of (remote) group */
    int           local_rank;     /* rank in local_group of this process */
    int           *lrank_to_grank;/* mapping for group */
    MPIR_CONTEXT   send_context;  /* context to send messages */
    MPIR_CONTEXT   recv_context;  /* context to recv messages */
    void          *ADIctx;        /* Context (if any) for abstract device */

    /* This stuff is needed for the communicator implemenation, but less
       often than the above items */
    MPIR_COMM_TYPE comm_type;	  /* inter or intra */
    MPI_Group     group;	  /* group associated with communicator */
    MPI_Group     local_group;    /* local group */
    MPI_Comm      comm_coll;      /* communicator for collective ops */

    int            ref_count;     /* number of references to communicator */
    void          *comm_cache;	  /* Hook for communicator cache */
    MPIR_HBT      attr_cache;     /* Hook for attribute cache */
    int           use_return_handler;   /* Allows us to override error_handler
					   when the MPI implementation
					   calls MPI routines */
#ifdef NEW_POINTERS
    MPI_Errhandler error_handler; /* Error handler structure */
#else
    struct MPIR_Errhandler *error_handler;  /* Error handler structure */
#endif
    int            permanent;      /* Is this a permanent object? */
    void          *mutex;          /* Local for threaded versions */

    /*** BEGIN HETEROGENEOUS ONLY ***/
    int           is_homogeneous; /* True if communicator's members
				     are homogeneous in data representation
				     (not used yet) */
    int           msgrep;         /* Message representation form for 
				     this communicator.  This is either
				     0 == MPIR_MSGREP_RECEIVER (all the same)
				     or != 0 (sender or XDR) */
    /* Note that point-to-point information on message representations
       is managed directly by the device and is not duplicated in the
       communicator */
    /*** END HETEROGENEOUS ONLY ***/

    /* These are used to support collective operations in this context */
    void          *adiCollCtx;
    MPIR_COLLOPS  collops;
#ifndef MPI_ADI2
    void          *ADIBarrier;
    void          *ADIReduce;
    void          *ADIScan;
    void          *ADIBcast;
    void          *ADICollect;

    /* These are only required to allow debuggers a way to locate
     * all of the communicators in the code, and provide a print name
     * for each. (The user may be able to set this name, at some point).
     */
    struct MPIR_COMMUNICATOR *comm_next; /* A chain through all 
					    communicators */
    char 		     *comm_name; /* A print name for this 
					    communicator */
#endif
};

/*
 * The list of all communicators in the program.
 */
typedef struct _MPIR_Comm_list {
  int	 	   	sequence_number;
  struct MPIR_COMMUNICATOR * comm_first;
} MPIR_Comm_list ;

extern MPIR_Comm_list MPIR_All_communicators;

#ifdef FOO
/* not a bad idea, but not used consistently */
typedef enum {
  MPIR_NO = 0, MPIR_YES
} MPIR_BOOL;
#endif

#ifdef FOO
/* 
   MPIR_FORT_INT is a special integer type for systems where sizeof(int) !=
   sizeof(Fortran REAL).  See the Fortran standard for why this is required.

   MPIR_LOGICAL is required for Fortran, since while the length of a 
   Fortran LOGICAL == Fortran INTEGER, the true/false values are 
   unspecified 
 */
typedef enum {
    MPIR_INT, MPIR_FLOAT, MPIR_DOUBLE, MPIR_COMPLEX, MPIR_LONG, MPIR_SHORT,
    MPIR_CHAR, MPIR_BYTE, MPIR_UCHAR, MPIR_USHORT, MPIR_ULONG, MPIR_UINT,
    MPIR_CONTIG, MPIR_VECTOR, MPIR_HVECTOR, MPIR_INDEXED,
    MPIR_HINDEXED, MPIR_STRUCT, MPIR_DOUBLE_COMPLEX, MPIR_PACKED, 
	MPIR_UB, MPIR_LB, MPIR_LONGDOUBLE, MPIR_LONGLONGINT, 
    MPIR_LOGICAL, MPIR_FORT_INT 
} MPIR_NODETYPE;
#endif

/* These are used by some of the datatype routines */
#ifdef FOO
#ifndef MPIR_TRUE
#define MPIR_TRUE  1
#define MPIR_FALSE 0
#endif
#endif

/* These are used in the Group manipulation routines */
#define MPIR_UNMARKED 0
#define MPIR_MARKED   1

/* #include "mpi_bc.h" */
#endif
