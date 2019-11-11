/*
 *  $Id$
 *
 *  (C) 1993 by Argonne National Laboratory and Mississipi State University.
 *      All rights reserved.  See COPYRIGHT in top-level directory.
 */

/* user include file for MPI programs, with no dependencies */

#ifndef _MPI_INCLUDE_BC
#define _MPI_INCLUDE_BC

/* assorted defined constants */
/* Error codes */
#include "mpi_error.h"            

#ifdef FOO
/* Generic error handling code.  This handles inserting the file and line
   number (in MPI) where the error occured.  In addition, it
   checks the error handler and calls the appropriate one.  Finally, 
   it returns the errorcode as its value.
 */
#define MPIR_ERROR(comm,code,string) \
    MPIR_Error( comm, code, string, __FILE__, __LINE__ )

/* Here we define some additional error information values.  These need to be
   or'ed into the appropriate MPI error class (from mpi_errno.h) */
#define MPIR_ERR_CLASS_BITS 8
#define MPIR_ERR_CLASS_MASK 0xff
           /* Uncommitted datatype - MPI_ERR_TYPE */  
#define MPIR_ERR_UNCOMMITTED  (1 << MPIR_ERR_CLASS_BITS) 

           /* Operation not defined for this datatype - MPI_ERR_OP */
#define MPIR_ERR_NOT_DEFINED  (1 << MPIR_ERR_CLASS_BITS)

           /* BSend with insufficent buffer space */
#define MPIR_ERR_USER_BUFFER_EXHAUSTED (1 << MPIR_ERR_CLASS_BITS)

           /* User has aliased an argument */
#define MPIR_ERR_BUFFER_ALIAS (2 << MPIR_ERR_CLASS_BITS)

           /* Address of location given to MPI_ADDRESS does fit in 
	      Fortran int */
#define MPIR_ERR_FORTRAN_ADDRESS_RANGE (4 << MPIR_ERR_CLASS_BITS)
#endif

/* communication modes */
typedef enum { 
    MPIR_MODE_STANDARD = 0, 
    MPIR_MODE_READY = 0x1, 
    MPIR_MODE_SYNCHRONOUS = 0x2, 
    MPIR_MODE_BUFFERED = 0x3,
    MPIR_MODE_SYNC_ACK = 0x4,
    MPIR_MODE_RECV = 0x5         /* This mode is not used by the ADI */
} MPIR_Mode;

/* Value of tag in status for a cancelled message */
#define MPIR_MSG_CANCELLED (-3)

/* This is the only global state in MPI */
extern int MPIR_Has_been_initialized;

#endif

