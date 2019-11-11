#ifndef MPIR_ERROR

/* 
 * Macros to simplify error handling
 * MPIR_ERROR - used primarily in "return MPIR_ERROR( comm, MPI_ERR_ccc, "")
 * MPIR_RETURN(comm,mpi_errno,"") use instead of return mpi_errno (calls
 *  MPIR_ERROR if mpi_errno != MPI_SUCCESS).
 * MPIR_CALL(fcn,comm,msg) calls an MPI routine and returns if error
 * 
 * The following are to allow MPI routines to call other MPI routines and
 * get the "correct" error behavior (i.e., return up to the outermost caller).
 * MPIR_ERROR_DECL - declaration (holds previous state)
 * MPIR_ERROR_PUSH(comm) - Change error handler on comm
 * MPIR_ERROR_POP(comm) - Change error handler on comm
 * MPIR_CALL_POP(fcn,comm,msg) - like MPIR_CALL, but also does MPIR_ERROR_POP.
 */

#ifndef ANSI_ARGS
#if defined(__STDC__) || defined(__cplusplus) || defined(HAVE_PROTOTYPES)
#define ANSI_ARGS(a) a
#else
#define ANSI_ARGS(a) ()
#endif
#endif

/* Generic error handling code.  This handles inserting the file and line
   number (in MPI) where the error occured.  In addition, it
   checks the error handler and calls the appropriate one.  Finally, 
   it returns the errorcode as its value.
 */
int MPIR_Error ANSI_ARGS(( struct MPIR_COMMUNICATOR *, int, char *, char *, int ));

#define MPIR_ERROR(comm,code,string) \
    MPIR_Error( comm, code, string, __FILE__, __LINE__ )
#define MPIR_RETURN(comm,code,string) \
    return (code) ? MPIR_ERROR(comm,code,string) : code
#define MPIR_ERROR_DECL int mpi_comm_err_ret
#define MPIR_ERROR_PUSH(comm) {mpi_comm_err_ret = (comm)->use_return_handler;\
        (comm)->use_return_handler = 1;}
#define MPIR_ERROR_POP(comm) (comm)->use_return_handler = mpi_comm_err_ret
#define MPIR_RETURN_POP(comm,code,string) \
    {MPIR_ERROR_POP(comm);MPIR_RETURN(comm,code,string);}

/* 
 * This routine can be called to call an MPI function and call the
 * appropriate error handler.
 */
#define MPIR_CALL(fcn,comm,msg) {if ((mpi_errno = fcn) != 0) \
				 return MPIR_ERROR(comm,mpi_errno,msg);}
#define MPIR_CALL_POP(fcn,comm,msg) {if ((mpi_errno = fcn) != 0) {\
	MPIR_ERROR_POP(comm); return MPIR_ERROR(comm,mpi_errno,msg);}}

/*
 * These can be called to handle allocating storage that might fail, returning
 * a NULL value in that case.  Note that "internal" routines like trmalloc
 * and MPID_SBalloc should just return a null, so that the appropriate 
 * error handler can be invoked
 */
#define MPIR_ALLOC(ptr,fcn,comm,code,msg) \
   {if (!((ptr) = fcn)) {return MPIR_ERROR(comm,code,msg);}}
/* MPIR_FALLOC is for Fortran */
#define MPIR_FALLOC(ptr,fcn,comm,code,msg) \
   {if (!((ptr) = fcn)) {*__ierr = MPIR_ERROR(comm,code,msg);return;}}
#define MPIR_ALLOC_POP(ptr,fcn,comm,code,msg) \
   {if (!((ptr) = fcn)) {MPIR_RETURN_POP(comm,code,msg);}}

#define MPIR_MAX_ARGS 10
extern void *(MPIR_errargs[MPIR_MAX_ARGS]);
extern int    MPIR_errargcnt;

#define MPIR_ERROR_PUSH_ARG(ptr) MPIR_errargs[MPIR_errargcnt++] = (void*)(ptr)

/* Here we define some additional error information values.  These need to be
   or'ed into the appropriate MPI error class (from mpi_errno.h) 
 */
#define MPIR_ERR_CLASS_BITS 8
#define MPIR_ERR_CLASS_MASK 0xff

/* 
   The various error codes are in the second 8 bits.  We reserve the 
   remaining 16 bits to indicate special error handling, for example,
   to indicate that runtime data for the message is available
 */
#define MPIR_ERR_CODE_BITS 16
#define MPIR_ERR_CODE_MASK 0xff00

/* Here are error CODE bits, to be or'ed with the error CLASS.
   In addition, some MPI_ERR types are defined that are CODES, not CLASSES 
 */

/* These are all error CODES mapped onto some of the error CLASSES.
   The error CLASS is the low byte; the code is in the second byte 
 */

/* MPI_ERR_COMM */
#define MPIR_ERR_COMM_NULL    (1 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_COMM_NULL    (MPIR_ERR_COMM_NULL | MPI_ERR_COMM)
                                    /* NULL communicator argument 
				       passed to function */
#define MPIR_ERR_COMM_INTER   (2 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_COMM_INTER   (MPIR_ERR_COMM_INTER | MPI_ERR_COMM)
			            /* Intercommunicator is not allowed 
				       in function */
#define MPIR_ERR_COMM_INTRA   (3 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_COMM_INTRA   (MPIR_ERR_COMM_INTRA | MPI_ERR_COMM)      
                                    /* Intracommunicator is not allowed 
				       in function */
#define MPIR_ERR_COMM_CORRUPT (4 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_COMM_CORRUPT  (MPIR_ERR_COMM_CORRUPT | MPI_ERR_COMM)

/* MPI_ERR_GROUP */
#define MPIR_ERR_GROUP_NULL   (1 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_GROUP_NULL    (MPIR_ERR_GROUP_NULL | MPI_ERR_GROUP)

#define MPIR_ERR_GROUP_CORRUPT (2 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_GROUP_CORRUPT (MPIR_ERR_GROUP_CORRUPT | MPI_ERR_GROUP)

/* MPI_ERR_TYPE */
#define MPIR_ERR_UNCOMMITTED  (1 << MPIR_ERR_CLASS_BITS) 
                                    /* Uncommitted datatype */  
#define MPI_ERR_UNCOMMITTED   (MPIR_ERR_UNCOMMITTED | MPI_ERR_TYPE)

#define MPIR_ERR_TYPE_NULL    (2 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_TYPE_NULL     (MPIR_ERR_TYPE_NULL | MPI_ERR_TYPE)

#define MPIR_ERR_TYPE_CORRUPT (3 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_TYPE_CORRUPT  (MPIR_ERR_TYPE_CORRUPT | MPI_ERR_TYPE)

#define MPIR_ERR_PERM_TYPE    (4 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_PERM_TYPE    (MPIR_ERR_PERM_TYPE | MPI_ERR_TYPE)      
                                    /* Can't free a perm type */
/* MPI_ERR_OP */
#define MPIR_ERR_NOT_DEFINED  (1 << MPIR_ERR_CLASS_BITS)
                                    /* Operation not defined for this 
				      datatype */
#define MPI_ERR_OP_NOT_DEFINED (MPIR_ERR_NOT_DEFINED | MPI_ERR_OP)

#define MPIR_ERR_OP_NULL      (2 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_OP_NULL       (MPIR_ERR_OP_NULL | MPI_ERR_OP)

/* MPI_ERR_ARG */
#define MPIR_ERR_ERRORCODE    (1 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_ERRORCODE    (MPIR_ERR_ERRORCODE | MPI_ERR_ARG)
			            /* Invalid error code */
#define MPIR_ERR_NULL         (2 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_NULL         (MPIR_ERR_NULL | MPI_ERR_ARG)
                                    /* Null parameter */
#define MPIR_ERR_PERM_KEY     (4 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_PERM_KEY     (MPIR_ERR_PERM_KEY | MPI_ERR_ARG)
                                    /* Can't free a perm key */
#define MPIR_ERR_PERM_OP      (6 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_PERM_OP      (MPIR_ERR_PERM_OP | MPI_ERR_ARG)      
                                    /* Can't free a permanent operator */
#define MPIR_ERR_FORTRAN_ADDRESS_RANGE (7 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_FORTRAN_ADDRESS_RANGE \
                             (MPIR_ERR_FORTRAN_ADDRESS_RANGE | MPI_ERR_ARG)
           /* Address of location given to MPI_ADDRESS does not fit in 
	      Fortran int */
#define MPIR_ERR_PERM_GROUP      (8 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_PERM_GROUP      (MPIR_ERR_PERM_GROUP | MPI_ERR_ARG)      

#define MPIR_ERR_KEYVAL          (9 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_KEYVAL           (MPIR_ERR_KEYVAL | MPI_ERR_ARG)      
				 /* Can't free a permanent group */   

#define MPIR_ERR_ERRHANDLER_NULL (10 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_ERRHANDLER_NULL  (MPIR_ERR_ERRHANDLER_NULL | MPI_ERR_ARG)

#define MPIR_ERR_ERRHANDLER_CORRUPT  (11 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_ERRHANDLER_CORRUPT (MPIR_ERR_ERRHANDLER_CORRUPT | MPI_ERR_ARG)

/* MPI_ERR_BUFFER */
#define MPIR_ERR_BUFFER_EXISTS (1 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_BUFFER_EXISTS (MPIR_ERR_BUFFER_EXISTS | MPI_ERR_BUFFER)

#define MPIR_ERR_USER_BUFFER_EXHAUSTED (2 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_USER_BUFFER_EXHAUSTED \
                            (MPIR_ERR_USER_BUFFER_EXHAUSTED | MPI_ERR_BUFFER)
                                    /* BSend with insufficent buffer space */
#define MPIR_ERR_BUFFER_ALIAS (3 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_BUFFER_ALIAS (MPIR_ERR_BUFFER_ALIAS | MPI_ERR_BUFFER)
                                    /* User has aliased an argument */
#define MPIR_ERR_BUFFER_SIZE (4 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_BUFFER_SIZE (MPIR_ERR_BUFFER_SIZE | MPI_ERR_BUFFER)

/* MPI_ERR_COUNT */
#define MPIR_ERR_COUNT_ARRAY_NEG (1 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_COUNT_ARRAY_NEG (MPIR_ERR_COUNT_ARRAY_NEG | MPI_ERR_COUNT)

/* MPI_ERR_OTHER */
#define MPIR_ERR_LIMIT        (1 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_LIMIT        (MPIR_ERR_LIMIT | MPI_ERR_OTHER)
                                    /* limit reached */
#define MPIR_ERR_NOMATCH      (2 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_NOMATCH      (MPIR_ERR_NOMATCH | MPI_ERR_OTHER)
                                    /* no recv posted for ready send */

#define MPIR_ERR_INIT         (3 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_INIT         (MPIR_ERR_INIT | MPI_ERR_OTHER)
                                    /* MPI_INIT already called */
#define MPIR_ERR_PRE_INIT     (4 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_PRE_INIT     (MPIR_ERR_PRE_INIT | MPI_ERR_OTHER)
                                    /* MPI_INIT has not been called */

#define MPIR_ERR_MPIRUN       (5 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_MPIRUN        (MPIR_ERR_MPIRUN | MPI_ERR_OTHER)

#define MPIR_ERR_BAD_INDEX    (6 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_BAD_INDEX     (MPIR_ERR_BAD_INDEX | MPI_ERR_OTHER)

#define MPIR_ERR_INDEX_EXHAUSTED (7 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_INDEX_EXHAUSTED (MPIR_ERR_INDEX_EXHAUSTED | MPI_ERR_OTHER)

#define MPIR_ERR_INDEX_FREED  (8 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_INDEX_FREED   (MPIR_ERR_INDEX_FREED | MPI_ERR_OTHER)

#define MPIR_ERR_BUFFER_TOO_SMALL (9 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_BUFFER_TOO_SMALL (MPIR_ERR_BUFFER_TOO_SMALL | MPI_ERR_OTHER)

/* MPI_ERR_INTERN */
#define MPIR_ERR_EXHAUSTED    (1 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_EXHAUSTED    (MPIR_ERR_EXHAUSTED | MPI_ERR_INTERN)
                                    /* Memory exhausted */
#define MPIR_ERR_ONE_CHAR     (2 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_ONE_CHAR     (MPIR_ERR_ONE_CHAR | MPI_ERR_INTERN)

#define MPIR_ERR_MSGREP_SENDER (3 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_MSGREP_SENDER  (MPIR_ERR_MSGREP_SENDER | MPI_ERR_INTERN)

#define MPIR_ERR_MSGREP_UNKNOWN (4 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_MSGREP_UNKNOWN  (MPIR_ERR_MSGREP_UNKNOWN | MPI_ERR_INTERN)

/* MPI_ERR_REQUEST */
#define MPIR_ERR_REQUEST_NULL (1 << MPIR_ERR_CLASS_BITS)
#define MPI_ERR_REQUEST_NULL  (MPIR_ERR_REQUEST_NULL | MPI_ERR_REQUEST)

/* 
   Standardized argument testing

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
        return MPIR_ERROR( comm, mpi_errno, "MPI_routine" );

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

/*
 * Some compilers may complain about (a=b) tests.  They should be upgraded
 * to do what the GNU gcc compiler does: complain about a=b but not about
 * (a=b).  If you ABSOLUTELY must shut up a hostile compiler, change
 * (a=b) to ((a=b),1). 
 */

/*
 * Tag tests.  In order to detect "tag too large", we need to check the
 * tag value against the maximum tag value.  This is MPID_TAG_UB from
 * the device (which should normally be a compile time constant, but
 * could be a global variable if the device wants that option)
 */

#define MPIR_TEST_SEND_TAG(comm,tag) \
    ((((tag) < 0 ) && (MPIR_ERROR_PUSH_ARG(&tag),mpi_errno = MPI_ERR_TAG )) ||\
     (((tag) > MPID_TAG_UB) && (MPIR_ERROR_PUSH_ARG(&tag),mpi_errno = MPI_ERR_TAG)))
    /* This requires MPI_ANY_TAG == -1 */
#define MPIR_TEST_RECV_TAG(comm,tag) \
    (( ((tag) < MPI_ANY_TAG) && \
         (MPIR_ERROR_PUSH_ARG(&tag),mpi_errno = MPI_ERR_TAG )) || \
    (((tag) > MPID_TAG_UB) && (MPIR_ERROR_PUSH_ARG(&tag),mpi_errno = MPI_ERR_TAG)))
    /* This exploits MPI_ANY_SOURCE==-2, MPI_PROC_NULL==-1 */
#define MPIR_TEST_SEND_RANK(comm,rank) \
    ( ((rank) < MPI_PROC_NULL || (rank) >= (comm)->np)\
           && (MPIR_ERROR_PUSH_ARG(&rank),mpi_errno = MPI_ERR_RANK))
    /* This requires min(MPI_PROC_NULL,MPI_ANY_SOURCE)=-2 */
#define MPIR_TEST_RECV_RANK(comm,rank) \
    (((rank) < -2 || (rank) >= (comm)->np) && \
     (MPIR_ERROR_PUSH_ARG(&rank),mpi_errno = MPI_ERR_RANK))
#define MPIR_TEST_COUNT(comm,count) ( ((count) < 0) && \
				     (mpi_errno = MPI_ERR_COUNT))
#ifdef NEW_POINTERS
#define MPIR_TEST_OP(comm,op) 'fixme'
#else
#define MPIR_TEST_OP(comm,op)       \
    ( (!(op) MPIR_TEST_COOKIE(op,MPIR_OP_COOKIE)) && (mpi_errno = MPI_ERR_OP ))
#endif
#ifdef NEW_POINTERS
#define MPIR_TEST_GROUP(comm,group) 'fixme'
#else
#define MPIR_TEST_GROUP(comm,group) \
    ( (!(group) MPIR_TEST_COOKIE(group,MPIR_GROUP_COOKIE)) && \
       (mpi_errno = MPI_ERR_GROUP ))
#endif
#ifdef NEW_POINTERS
#define MPIR_TEST_COMM(comm,comm1)  'fixme'
#else
#define MPIR_TEST_COMM(comm,comm1)  \
    ( (!(comm1) MPIR_TEST_COOKIE(comm1,MPIR_COMM_COOKIE)) \
     && (mpi_errno = MPI_ERR_COMM ))
#endif
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
#define MPIR_TEST_DATATYPE(comm,datatype) 'fixme'
/*     (!(datatype) && (mpi_errno = MPI_ERR_TYPE)) */
#ifdef FOO
#define MPIR_TEST_DATATYPE(comm,datatype) \
    (MPIR_TEST_IS_DATATYPE(comm,datatype) || \
  (!MPIR_TEST_PREDEF_DATATYPE(datatype) && !(datatype)->committed && \
   (mpi_errno = (MPI_ERR_TYPE | MPIR_ERR_UNCOMMITTED))))
#endif

#ifdef NEW_POINTERS
#define MPIR_TEST_ERRHANDLER(comm,errhandler) 'fixme'
#else
#define MPIR_TEST_ERRHANDLER(comm,errhandler) \
    ( ( (!(errhandler) MPIR_TEST_COOKIE(errhandler,MPIR_ERRHANDLER_COOKIE)) \
       && (mpi_errno = MPI_ERR_ARG )))
#endif
#define MPIR_TEST_HBT_NODE(comm,node) \
    ( ( !(node) MPIR_TEST_COOKIE(node,MPIR_HBT_NODE_COOKIE)) \
      && (mpi_errno = MPI_ERR_INTERN))
#define MPIR_TEST_HBT(comm,hbt) \
    ( ( !(hbt) MPIR_TEST_COOKIE(hbt,MPIR_HBT_COOKIE)) \
      && (mpi_errno = MPI_ERR_INTERN))

#define MPIR_TEST_ALIAS(b1,b2)      \
    ( ((b1)==(b2)) && (mpi_errno = (MPI_ERR_BUFFER | MPIR_ERR_BUFFER_ALIAS) ))
#define MPIR_TEST_ARG(arg)  (!(arg) && (mpi_errno = MPI_ERR_ARG) )
#endif 

/* 
   Here are the definitions of the actual error messages; this is also needed
   by end-users (MPI error names are visible to all)
 */
#include "mpi_errno.h"

#endif
