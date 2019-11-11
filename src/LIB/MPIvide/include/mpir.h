/*
 *  $Id$
 *
 *  (C) 1993 by Argonne National Laboratory and Mississipi State University.
 *      All rights reserved.  See COPYRIGHT in top-level directory.
 */

/* include file for the MPIR implementation of MPI */


#ifndef _MPIR_INCLUDE
#define _MPIR_INCLUDE

#include <stdio.h>
#include "mpi.h"
#include "mpi_bc.h"
#include "dmpi.h"
#include "dm.h"

#ifndef ANSI_ARGS
#if defined(__STDC__) || defined(__cplusplus) || defined(HAVE_PROTOTYPES)
#define ANSI_ARGS(a) a
#else
#define ANSI_ARGS(a) ()
#endif
#endif

#ifndef VOLATILE
#if defined(__STDC__) || defined(__cplusplus)
#define VOLATILE volatile
#else
#define VOLATILE
#endif
#endif
/*****************************************************************************
*                           MPIR DATA STRUCTURES                             *
*****************************************************************************/

/*****************************************************************************/
/* Datatypes.  The contiguous, predefined datatypes are handled separately   */
/* to demonstrate that the added functionality has low cost                  */
/*****************************************************************************/

#include "../mpid/ch2/datatype.h"
#ifdef FOO
#define MPIR_DATATYPE_COOKIE 0xea31beaf
struct MPIR_DATATYPE {
    MPIR_NODETYPE dte_type; /* type of datatype element this is */
    MPIR_COOKIE             /* Cookie to help detect valid item */
    int          committed; /* whether committed or not */
    int          is_contig; /* whether entirely contiguous */
    int              basic; /* Is this a basic type */
    int          permanent; /* Is this a permanent type */
    MPI_Aint        ub, lb; /* upper/lower bound of type */
    MPI_Aint real_ub, real_lb; /* values WITHOUT TYPE_UB/TYPE_LB */
    int             has_ub; /* Indicates that the datatype has a TYPE_UB */
    int             has_lb; /* Indicates that the datatype has a TYPE_LB */
    int             extent; /* extent of this datatype */
    int               size; /* size of type */
    int           elements; /* number of basic elements */
    int          ref_count; /* nodes depending on this node */
    int              align; /* alignment needed for start of datatype */
    int              count; /* replication count */
    MPI_Aint        stride; /* stride, for VECTOR and HVECTOR types */
    MPI_Aint      *indices; /* array of indices, for (H)INDEXED, STRUCT */
    int           blocklen; /* blocklen, for VECTOR and HVECTOR types */
    int         *blocklens; /* array of blocklens for (H)INDEXED, STRUCT */
    MPI_Datatype old_type;  /* type this type is built of, if 1 */
    MPI_Datatype *old_types;/* array of types, for STRUCT */
    MPI_Datatype flattened; /* Flattened version, if available */
};

/* Holds translation from integer to MPI_Datatype */
#define MPIR_MAX_DATATYPE_ARRAY 256
/* defined in src/env/initutil.c */
/* MPIR_datatypes[0] is always 0, so null->null requires no special test */
extern MPI_Datatype MPIR_datatypes[MPIR_MAX_DATATYPE_ARRAY];

/* Translate between index and datatype pointer */
#define MPIR_GET_REAL_DATATYPE(a) \
  {if(MPIR_TEST_PREDEF_DATATYPE(a)) a = MPIR_datatypes[(MPI_Aint)(a)];}
#define MPIR_TEST_PREDEF_DATATYPE(a) \
    ((MPI_Aint)(a)<MPIR_MAX_DATATYPE_ARRAY && (MPI_Aint)(a) >0)
#define MPIR_DATATYPE_CONTIG(a) \
    (MPIR_TEST_PREDEF_DATATYPE(a) || (a)->is_contig)
/* For ONLY the predefined datatypes, the size MAY be encoded in the 
   value of the datatype */
#define MPIR_DATATYPE_SIZE(a) (1 + ( (MPI_Aint)(a)&0xf ) )
#endif

/*****************************************************************************/
/* Requests come in many flavors                                             */
/*****************************************************************************/

typedef enum {
    MPIR_SEND,
    MPIR_RECV,
    MPIR_PERSISTENT_SEND,
    MPIR_PERSISTENT_RECV,
    MPIR_USER
} MPIR_OPTYPE;

/* MPIR_COMMON is a subset of the handle structures that contain common
   elements.  MAKE SURE THAT ANY CHANGES IN THE HANDLES PRESERVES
   THIS ORDER! */
#define MPIR_REQUEST_COOKIE 0xe0a1beaf
typedef struct {
    MPIR_OPTYPE handle_type;
    MPIR_COOKIE                 /* Cookie to help detect valid item */
    int         contextid;
    int         partner;        /* source or destination, depending on type */
    int         tag;
    int         completer;      /* index of routine to complete, 
				   or 0 if done */
    int         errval;         /* Holds any error code; 0 for none */

    MPI_Datatype datatype;
    MPI_Comm    comm;           /* The communicator this request is 
				   on. Needed for MPI_Start in
				   multiprotocol and heterogeneous
				   communication */
    int         persistent;     /* is handle persistent (created with 
				   MPI_..._init)? */
    int         active;         /* relavent only for persistent handles,
				   indicates that the handle is active */
    int         msgrep;         /* Message representation; used to indicate
				   XDR used */

    void *bufadd;       /* address of buffer */
    int buflen;         /* length of buffer at bufadd */
    int count;          /* requested count */

    char *bufpos;       /* position of next byte to be transferred */
    int  totallen;      /* length in local bytes of message */

    } MPIR_COMMON;

typedef struct {
    MPIR_OPTYPE handle_type;    
    MPIR_COOKIE             /* Cookie to help detect valid item */
    int         contextid;  /* context id */
    int         dest;       /* destination process for message */
    int         tag;        /* tag */
    int         completer;      /* index of routine to complete, 
				   or 0 if done */
    int         errval;         /* Holds any error code; 0 for none */

    MPI_Datatype datatype;  /* basic or derived datatype */
    MPI_Comm    comm;
    int         persistent;
    int         active;
    int         msgrep;         /* Message representation; used to indicate
				   XDR used */
    void *bufadd;       /* address of buffer */
    int buflen;         /* length of buffer at bufadd */
    int count;          /* requested count */
    char *bufpos;       /* position of next byte to be transferred */
                        /* This is only used to point to a temporary buffer
			   to send out of. The buffer pointed to by this
			   is freed at the end of blocking sends and 
			   waits on non-blocking sends */
    int  totallen;      /* length in local bytes of message */

    MPIR_Mode mode;     /* mode: STANDARD, READY, SYNC, BUFFERED */
    int  lrank;                 /* Rank in sending group */
    void *bsend;        /* Pointer to structure managed for buffered sends */

    MPID_SHANDLE dev_shandle;   /* device's version of send handle */
} MPIR_SHANDLE;

typedef struct {
    MPIR_OPTYPE handle_type;    
    MPIR_COOKIE             /* Cookie to help detect valid item */
    int         contextid;  /* context id */
    int         source;     /* source process message */
    int         tag;        /* tag */
    int         completer;      /* index of routine to complete, 
				   or 0 if done */
    int         errval;         /* Holds any error code; 0 for none */

    MPI_Datatype datatype;  /* basic or derived datatype */
    MPI_Comm    comm;
    int         persistent;
    int         active;
    int         msgrep;         /* Message representation; used to indicate
				   XDR used */
    void *bufadd;       /* address of buffer */
    int buflen;         /* length of buffer at bufadd */
    int count;          /* requested count */
    char *bufpos;       /* position of next byte to be transferred */
    int  totallen;          /* length in local bytes of message */

    void *p;            /* Pointer to unexpected data */
    int  len;           /* Length of this data ? */
    int  actcount;      /* number of items actually read */

    MPID_RHANDLE dev_rhandle;   /* device's version of recv handle */

    /* For persistant receives, we may need to restore some of the fields
       after the operation completes */
    int perm_source, perm_tag, perm_len;

} MPIR_RHANDLE;

/* This is an "extension" handle and is NOT part of the MPI standard.
   Defining it, however, introduces no problems with the standard, and
   it allows us to easily extend the request types.

   Note that this is not yet compatible with the essential fields of
   MPIR_COMMON.
 */
typedef struct {
    MPIR_OPTYPE handle_type;    
    MPIR_COOKIE                 /* Cookie to help detect valid item */
    int         completer;      /* index of routine to complete, 
				   or 0 if done */
    int         active;
    int         (*create_ureq) ANSI_ARGS((MPI_Request));
    int         (*free_ureq)   ANSI_ARGS((MPI_Request));
    int         (*wait_ureq)   ANSI_ARGS((MPI_Request));
    int         (*test_ureq)   ANSI_ARGS((MPI_Request));
    int         (*start_ureq)  ANSI_ARGS((MPI_Request));
    int         (*cancel_ureq) ANSI_ARGS((MPI_Request));
    void        *private_data;
} MPIR_UHANDLE;

#define MPIR_HANDLES_DEFINED

union MPIR_HANDLE {
    MPIR_OPTYPE type;           /* send or receive */
    MPIR_COMMON  chandle;       /* common fields */
    MPIR_SHANDLE shandle;
    MPIR_RHANDLE rhandle;
    MPIR_SHANDLE persistent_shandle;
    MPIR_RHANDLE persistent_rhandle;
    MPIR_UHANDLE uhandle;
};


/*
   MPIR manages queues of unexpected messages and posted receives.
   Originally, these were arranged as general queue elements, but 
   they have been modified to meet the needs of the MPIR implementation.

   The current implementation does not make use of the delivery-order
   fields.

   Because of the importance of this queue, we've added explicit support
   for MPI.  This is indicated by the fields:
   context_id
   tag, tagmask
   lsrc, srcmask

   (note that the last two may depart if we split the queues by source).
   The mask fields allow us to replace

   itag == tag || tag == MPI_ANY_TAG

   with 

   (itag & tagmask) == tag

   saving us a compare and branch at the cost of a load and bitwise-and.
   Generalizations of the tagmask could provide users with additional
   functionality (for MPICH EXTENSIONS).
 */
typedef enum {
    MPIR_QSHANDLE,
    MPIR_QRHANDLE
} MPIR_QEL_TYPE;

typedef struct _MPIR_QEL {  /* queue elements */
    void *ptr;              /* queues can contain anything  */
    int  context_id, 
         tag, tagmask,
         lsrc, srcmask;
    MPIR_QEL_TYPE qel_type; /* but we may need to know what */
    struct _MPIR_QEL *prev; /* previous queue element */
    struct _MPIR_QEL *next; /* next queue element */
    struct _MPIR_QEL *deliv_next;   /* next in delivery order */
} MPIR_QEL;

typedef struct {        /* header for queues of things like handles */
    MPIR_QEL *first;        /* first queue element */
    MPIR_QEL *last;         /* last queue element */
    MPIR_QEL *deliv_first;  /* first in delivery order */
    MPIR_QEL *deliv_last;   /* last in delivery order */
    int      currlen;       /* current length */
    int      maxlen;        /* maximum since beginning of run */
    MPID_THREAD_DS_LOCK_DECLARE   /* Used for controlling access by threads */
} MPIR_QHDR;

/* >>>> is this still in use <<<< */
typedef struct _MPIR_FDTEL {    /* flat datatype elements */
    MPIR_NODETYPE      type;    /* one of the elementary data types */
    int                disp;    /* displacement */
    struct _MPIR_FDTEL *next;   /* pointer to next element */
} MPIR_FDTEL;


/*****************************************************************************
*                                GLOBAL DATA                                 *
*****************************************************************************/

/* memory management for fixed-size blocks */
extern void *MPIR_shandles;     /* sbcnst MPIR_SHANDLES */
extern void *MPIR_rhandles;     /* sbcnst MPIR_RHANDLES */
extern void *MPIR_errhandlers;  /* sbcnst Error handlers */
extern void *MPIR_dtes;   /* sbcnst datatype elements */
extern void *MPIR_qels;   /* sbcnst queue elements */
extern void *MPIR_fdtels; /* sbcnst flat datatype elements */
extern void *MPIR_hbts;   /* sbcnst height balanced tree roots for cacheing */
extern void *MPIR_hbt_els;/* sbcnst height balanced tree nodes for cacheing */
extern void *MPIR_topo_els;/* sbcnst topology elements */

/* queues */
extern MPIR_QHDR MPIR_posted_recvs;
extern MPIR_QHDR MPIR_unexpected_recvs;

/* Predefined function tables for collective routines, the device
 * can also use its own, but these are the defaults.
 */
#ifdef FOO
defined in mpicoll.h
extern MPIR_COLLOPS * MPIR_inter_collops;   /* Simply raises appropriate error */
extern MPIR_COLLOPS * MPIR_intra_collops;   /* Do the business using pt2pt     */
#endif
/* MPIR routines are defined in mpiimpl.h */

/* MPIR process id (from device) */
extern int MPIR_tid;

/* Fortran logical values */
#ifndef _CRAY
extern int MPIR_F_TRUE, MPIR_F_FALSE;
#define MPIR_TO_FLOG(a) ((a) ? MPIR_F_TRUE : MPIR_F_FALSE)
/* 
   Note on true and false.  This code is only an approximation.
   Some systems define either true or false, and allow some or ALL other
   patterns for the other.  This is just like C, where 0 is false and 
   anything not zero is true.  Modify this test as necessary for your
   system.
 */
#define MPIR_FROM_FLOG(a) ( (a) == MPIR_F_TRUE ? 1 : 0 )

#else
/* CRAY Vector processors only; these are defined in /usr/include/fortran.h 
   Thanks to lmc@cray.com */
#define MPIR_TO_FLOG(a) (_btol(a))
#define MPIR_FROM_FLOG(a) ( _ltob(&(a)) )    /*(a) must be a pointer */
#endif

/* MPIR_F_MPI_BOTTOM is the address of the Fortran MPI_BOTTOM value */
extern void *MPIR_F_MPI_BOTTOM;

/* MPIR_F_PTR checks for the Fortran MPI_BOTTOM and provides the value 
   MPI_BOTTOM if found 
   See src/pt2pt/addressf.c for why MPIR_F_PTR(a) is just (a)
*/
/*  #define MPIR_F_PTR(a) (((a)==(MPIR_F_MPI_BOTTOM))?MPI_BOTTOM:a) */
#define MPIR_F_PTR(a) (a)

/*  
 * These are hooks for Fortran characters.
 * MPID_FCHAR_T is the type of a Fortran character argument
 * MPID_FCHAR_LARG is the "other" argument that some Fortran compilers use
 * MPID_FCHAR_STR gives the pointer to the characters
 */
#ifdef MPID_CHARACTERS_ARE_CRAYPVP
typedef <whatever> MPID_FCHAR_T;
#define MPID_FCHAR_STR(a) (a)->characters   <or whatever>
#define MPID_FCHAR_LARG(d) 
#else
typedef char *MPID_FCHAR_T;
#define MPID_FCHAR_STR(a) a
#define MPID_FCHAR_LARG(d) ,d
#endif

/*
 * This allows us to define the C type that corresponds to a Fortran
 * integer with e.g., -DMPIR_FORT_INT_T=long
 */
#ifndef MPIR_FORT_INT_T
#define MPIR_FORT_INT_T int
#endif

/* 
   Message tag ranges.  This may be overridden in mpid.h .
 */
#ifndef MPID_TAG_UB
#define MPID_TAG_UB (1<<((sizeof(int)*8) -2))
#endif

/* Message encodings - How messages are enecoded once they
   reach the device. MPI tries to put messages into the receiver's 
   format if the conversion is easy, i.e. switching bytes.
   If the conversion could be to difficult, then the messages use
   another form of encoding. (By default XDR) */

#define MPIR_MSGREP_UNKNOWN	-1
/* Encoded in the receiver's native format */
#define MPIR_MSGREP_RECEIVER	0
/* Encoded with XDR */
#define MPIR_MSGREP_XDR		1
/* Encoded in the sender's native format */
#define MPIR_MSGREP_SENDER	2

/* 
   Formats to use.  This is a little different from the above; these
   specify the actual conversion to perform 

   These are not the final forms, since we'll eventually want to allow
   forms that do swap-and-extend etc., but for now, these are ok.

   To aid in debugging, we make the non-XDR form different in value.
   Currently, these are used in src/dmpi/dmpi.c dmpipk.c and pkutil.c for
   the "destination type"
 */
#define MPIR_MSGFORM_XDR 1
#define MPIR_MSGFORM_SWAP 10
#define MPIR_MSGFORM_OK 0

/*
   These macros allow inlining of some common cases
 */     
#ifdef MPID_HAS_HETERO
#define  MPIR_SEND_SETUP_BUFFER( request, shandle ) \
   {mpi_errno = MPIR_Send_setup( request );} 
#else        
#define  MPIR_SEND_SETUP_BUFFER( request, shandle ) \
   if (shandle.datatype->is_contig) {\
   shandle.active       = 1;\
   shandle.dev_shandle.bytes_as_contig =\
   shandle.count * shandle.datatype->size;\
   if (shandle.dev_shandle.bytes_as_contig > 0 && shandle.bufadd == 0)\
       mpi_errno = MPI_ERR_BUFFER;\
   shandle.dev_shandle.start = shandle.bufadd;\
   shandle.bufpos		 = 0;}\
   else{mpi_errno = MPIR_Send_setup( request );} 
#endif
#ifdef MPID_HAS_HETERO
#define MPIR_RECV_SETUP_BUFFER( request, rhandle ) \
    {mpi_errno = MPIR_Receive_setup( request );}
#else
#define MPIR_RECV_SETUP_BUFFER( request, rhandle ) \
    if (rhandle.datatype->is_contig) {\
    rhandle.active       = 1;\
    rhandle.dev_rhandle.start = rhandle.bufadd;\
    rhandle.dev_rhandle.bytes_as_contig =\
      rhandle.count * rhandle.datatype->extent;\
    if (rhandle.dev_rhandle.bytes_as_contig > 0 && \
	rhandle.bufadd == 0) \
	mpi_errno = MPI_ERR_BUFFER;\
    rhandle.bufpos                      = 0;}\
    else {mpi_errno = MPIR_Receive_setup( request );}
#endif

/* 
   These restore any possible wild-cards in persistant receives.  These
   are used in the various Test/Wait calls
 */
#define MPIR_RESET_RECVINIT(rhandle) {\
     (rhandle)->source = (rhandle)->perm_source;\
     (rhandle)->tag    = (rhandle)->perm_tag;}
#define MPIR_RESET_PERSISTENT(request) {\
    (request)->chandle.active    = 0;\
    MPID_Clr_completed( MPID_Ctx( (request) ), (request) );\
    if ((request)->type == MPIR_RECV) {\
	MPID_Reuse_recv_handle( (request)->rhandle.comm->ADIctx,\
			      &(request)->rhandle.dev_rhandle );\
	MPIR_RESET_RECVINIT(&(request)->rhandle);\
	}\
    else {\
	MPID_Reuse_send_handle( (request)->shandle.comm->ADIctx,\
			      &(request)->shandle.dev_shandle );\
	}}
#define MPIR_RESET_PERSISTENT_RECV(request) {\
		(request)->chandle.active = 0;\
		MPID_Clr_completed( MPID_Ctx( request ), request );\
		MPID_Reuse_recv_handle( (request)->rhandle.comm->ADIctx, \
				        &(request)->rhandle.dev_rhandle );\
		MPIR_RESET_RECVINIT(&(request)->rhandle);}
#define MPIR_RESET_PERSISTENT_SEND(request) {\
	      (request)->chandle.active = 0;\
	      MPID_Clr_completed( MPID_Ctx(request), request );\
	      MPID_Reuse_send_handle( (request)->shandle.comm->ADIctx, \
				      &(request)->shandle.dev_shandle );}

#ifdef FOO

/* 
   Standardized error testing

   Many of the MPI routines take arguments of the same type.  These
   macros provide tests for these objects.

   It is intended that the tests for a valid opaque object such as 
   a communicator can check to insure that the object is both a communicator
   and that it is valid (hasn't been freed).  They can also test for
   null pointers.

   These are not used yet; we are still looking for the best ways to 
   define them.

   The intent is to use them in this manner:

   if (MPIR_TEST_...() || MPIR_TEST_... || ... ) 
        return MPIR_ERROR( comm, mpi_errno, "Error in MPI_routine" );

   The hope is, that in the NO_ERROR_CHECKING case, the optimizer will
   be smart enough to remove the code.
 */
#ifdef MPIR_NO_ERROR_CHECKING
#define MPIR_TEST_SEND_TAG(comm,tag)      0
#define MPIR_TEST_RECV_TAG(comm,tag)      0
#define MPIR_TEST_SEND_RANK(comm,rank)    0
#define MPIR_TEST_RECV_RANK(comm,rank)    0
#define MPIR_TEST_COUNT(comm,count)       0
#define MPIR_TEST_OP(comm,op)             0
#define MPIR_TEST_GROUP(comm,group)       0
#define MPIR_TEST_COMM(comm,comm1)        0
#define MPIR_TEST_REQUEST(comm,request)   0
#define MPIR_TEST_IS_DATATYPE(comm,datatype) 0
#define MPIR_TEST_DATATYPE(comm,datatype) 0
#define MPIR_TEST_ERRHANDLER(comm,errhandler) 0
#define MPIR_TEST_ALIAS(b1,b2)            0
#define MPIR_TEST_ARG(arg)                0

#else
#ifdef MPIR_HAS_COOKIES
#define MPIR_TEST_COOKIE(val,value) || ( ((val)->cookie != (value)) )
#else 
#define MPIR_TEST_COOKIE(val,value) 
#endif


#define MPIR_TEST_SEND_TAG(comm,tag) \
    ( ((tag) < 0 ) && (mpi_errno = MPI_ERR_TAG ))
    /* This requires MPI_ANY_TAG == -1 */
#define MPIR_TEST_RECV_TAG(comm,tag) \
    ( ((tag) < MPI_ANY_TAG) &&  (mpi_errno = MPI_ERR_TAG ))
    /* This exploits MPI_ANY_SOURCE==-2, MPI_PROC_NULL==-1 */
#define MPIR_TEST_SEND_RANK(comm,rank) \
    ( ((rank) < MPI_PROC_NULL || (rank) >= (comm)->np)\
           && (mpi_errno = MPI_ERR_RANK))
    /* This requires min(MPI_PROC_NULL,MPI_ANY_SOURCE)=-2 */
#define MPIR_TEST_RECV_RANK(comm,rank) \
    (((rank) < -2 || (rank) >= (comm)->np) && \
     (mpi_errno = MPI_ERR_RANK))
#define MPIR_TEST_COUNT(comm,count) ( ((count) < 0) && \
				     (mpi_errno = MPI_ERR_COUNT))
#define MPIR_TEST_OP(comm,op)       \
    ( (!(op) MPIR_TEST_COOKIE(op,MPIR_OP_COOKIE)) && (mpi_errno = MPI_ERR_OP ))
#define MPIR_TEST_GROUP(comm,group) \
    ( (!(group) MPIR_TEST_COOKIE(group,MPIR_GROUP_COOKIE)) && \
       (mpi_errno = MPI_ERR_GROUP ))
#define MPIR_TEST_COMM(comm,comm1)  \
    ( (!(comm1) MPIR_TEST_COOKIE(comm1,MPIR_COMM_COOKIE)) \
     && (mpi_errno = MPI_ERR_COMM ))
#define MPIR_TEST_REQUEST(comm,request) \
 ( (!(request) MPIR_TEST_COOKIE(&((request)->chandle),MPIR_REQUEST_COOKIE)) \
     && (mpi_errno = MPI_ERR_REQUEST))

#ifdef MPIR_HAS_COOKIES
#define MPIR_TEST_IS_DATATYPE(comm,datatype) \
    ( (!(datatype) || \
       (!MPIR_TEST_PREDEF_DATATYPE(datatype) && \
	((datatype)->cookie!=MPIR_DATATYPE_COOKIE))) \
     && (mpi_errno = MPI_ERR_TYPE ))
#else
#define MPIR_TEST_IS_DATATYPE(comm,datatype) \
    ( (!(datatype) ) && (mpi_errno = MPI_ERR_TYPE ))
#endif
#define MPIR_TEST_DATATYPE(comm,datatype) \
    (MPIR_TEST_IS_DATATYPE(comm,datatype) || \
  (!MPIR_TEST_PREDEF_DATATYPE(datatype) && !(datatype)->committed && \
   (mpi_errno = (MPI_ERR_TYPE | MPIR_ERR_UNCOMMITTED))))

#define MPIR_TEST_ERRHANDLER(comm,errhandler) \
    ( ( (!(errhandler) MPIR_TEST_COOKIE(errhandler,MPIR_ERRHANDLER_COOKIE)) \
       && (mpi_errno = MPI_ERR_ARG )))
#define MPIR_TEST_HBT_NODE(comm,node) \
    ( ( !(node) MPIR_TEST_COOKIE(node,MPIR_HBT_NODE_COOKIE)) \
      && (mpi_errno = MPI_ERR_INTERN))
#define MPIR_TEST_HBT(comm,hbt) \
    ( ( !(hbt) MPIR_TEST_COOKIE(hbt,MPIR_HBT_COOKIE)) \
      && (mpi_errno = MPI_ERR_INTERN))

#define MPIR_TEST_ALIAS(b1,b2)      \
    ( ((b1)==(b2)) && (mpi_errno = MPI_ERR_BUFFER_ALIAS ))
#define MPIR_TEST_ARG(arg)  (!(arg) && (mpi_errno = MPI_ERR_ARG) )
#endif 

#endif

#endif
