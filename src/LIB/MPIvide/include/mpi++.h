/*
 *  $Id$
 *
 *  (C) 1993 by Argonne National Laboratory and Mississipi State University.
 *      All rights reserved.  See COPYRIGHT in top-level directory.
 */

/* user include file for MPI++ programs */

#ifndef _MPIPP_INCLUDE
#define _MPIPP_INCLUDE

#include <mpi++P.h>

// External declarations 
extern MPIX_Comm_world   MPIX_COMM_WORLD;
extern MPIX_Comm         MPIX_COMM_NULL;
extern MPIX_Group        MPIX_GROUP_NULL, MPIX_GROUP_EMPTY;
extern MPIX_Datatype     MPIX_CHAR, MPIX_SHORT, MPIX_INT, MPIX_LONG;
extern MPIX_Datatype     MPIX_UNSIGNED_CHAR, MPIX_UNSIGNED_SHORT;
extern MPIX_Datatype     MPIX_UNSIGNED, MPIX_UNSIGNED_LONG, MPIX_FLOAT;
extern MPIX_Datatype     MPIX_DOUBLE, MPIX_LONG_DOUBLE, MPIX_BYTE, MPIX_PACKED;
extern MPIX_Op           MPIX_SUM, MPIX_MAX, MPIX_MIN, MPIX_PROD, MPIX_BXOR;
extern MPIX_Op           MPIX_LAND, MPIX_BAND, MPIX_LOR, MPIX_BOR, MPIX_LXOR;

// Various defines to allow use of "MPI" instead of "MPIX" naming
#ifdef  MPI_COMM_WORLD
#undef  MPI_COMM_WORLD
#endif
#define MPI_COMM_WORLD      MPIX_COMM_WORLD

#ifdef  MPI_COMM_NULL
#undef  MPI_COMM_NULL
#endif
#define MPI_COMM_NULL       MPIX_COMM_NULL

#ifdef  MPI_GROUP_NULL
#undef  MPI_GROUP_NULL
#endif
#define MPI_GROUP_NULL      MPIX_GROUP_NULL

#ifdef  MPI_GROUP_EMPTY
#undef  MPI_GROUP_EMPTY
#endif
#define MPI_GROUP_EMPTY     MPIX_GROUP_EMPTY

#ifdef  MPI_Comm
#undef  MPI_Comm
#endif
#define MPI_Comm            MPIX_Comm

#ifdef  MPI_Comm_world
#undef  MPI_Comm_world
#endif
#define MPI_Comm_world      MPIX_Comm_world

#ifdef  MPI_Group
#undef  MPI_Group
#endif
#define MPI_Group           MPIX_Group

#ifdef  MPI_Datatype
#undef  MPI_Datatype
#endif
#define MPI_Datatype        MPIX_Datatype

#ifdef  MPI_Op
#undef  MPI_Op
#endif
#define MPI_Op              MPIX_Op

#ifdef  MPI_CHAR
#undef  MPI_CHAR
#endif
#define MPI_CHAR            MPIX_CHAR

#ifdef  MPI_SHORT
#undef  MPI_SHORT
#endif
#define MPI_SHORT           MPIX_SHORT

#ifdef  MPI_INT
#undef  MPI_INT
#endif
#define MPI_INT             MPIX_INT

#ifdef  MPI_LONG
#undef  MPI_LONG
#endif
#define MPI_LONG            MPIX_LONG

#ifdef  MPI_UNSIGNED_CHAR
#undef  MPI_UNSIGNED_CHAR
#endif
#define MPI_UNSIGNED_CHAR   MPIX_UNSIGNED_CHAR

#ifdef  MPI_UNSIGNED_SHORT
#undef  MPI_UNSIGNED_SHORT
#endif
#define MPI_UNSIGNED_SHORT  MPIX_UNSIGNED_SHORT

#ifdef  MPI_UNSIGNED
#undef  MPI_UNSIGNED
#endif
#define MPI_UNSIGNED        MPIX_UNSIGNED

#ifdef  MPI_UNSIGNED_LONG
#undef  MPI_UNSIGNED_LONG
#endif
#define MPI_UNSIGNED_LONG   MPIX_UNSIGNED_LONG

#ifdef  MPI_FLOAT
#undef  MPI_FLOAT
#endif
#define MPI_FLOAT           MPIX_FLOAT

#ifdef  MPI_DOUBLE
#undef  MPI_DOUBLE
#endif
#define MPI_DOUBLE          MPIX_DOUBLE

#ifdef  MPI_LONG_DOUBLE
#undef  MPI_LONG_DOUBLE
#endif
#define MPI_LONG_DOUBLE     MPIX_LONG_DOUBLE

#ifdef  MPI_BYTE
#undef  MPI_BYTE
#endif
#define MPI_BYTE            MPIX_BYTE

#ifdef  MPI_PACKED
#undef  MPI_PACKED
#endif
#define MPI_PACKED          MPIX_PACKED

#ifdef  MPI_SUM
#undef  MPI_SUM
#endif
#define MPI_SUM             MPIX_SUM

#ifdef  MPI_MAX
#undef  MPI_MAX
#endif
#define MPI_MAX             MPIX_MAX

#ifdef  MPI_MIN
#undef  MPI_MIN
#endif
#define MPI_MIN             MPIX_MIN

#ifdef  MPI_ROD
#undef  MPI_ROD
#endif
#define MPI_ROD             MPIX_PROD

#ifdef  MPI_XOR
#undef  MPI_XOR
#endif
#define MPI_XOR             MPIX_BXOR

#ifdef  MPI_AND
#undef  MPI_AND
#endif
#define MPI_AND             MPIX_LAND

#ifdef  MPI_AND
#undef  MPI_AND
#endif
#define MPI_AND             MPIX_BAND

#ifdef  MPI_LOR
#undef  MPI_LOR
#endif
#define MPI_LOR             MPIX_LOR

#ifdef  MPI_BOR
#undef  MPI_BOR
#endif
#define MPI_BOR             MPIX_BOR

#ifdef  MPI_XOR
#undef  MPI_XOR
#endif
#define MPI_XOR             MPIX_LXOR

#endif
