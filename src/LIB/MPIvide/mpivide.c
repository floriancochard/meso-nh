/* 
MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
MNH_LIC for details. version 1.
*/ 
#include <string.h>
#include "mpi.h"

/* Variables defined in meso-nh code */
#ifdef FUJI
#if MNH_REAL == 4
  #define MPI_PRECISION MPI_REAL
  #define MPI_2PRECISION MPI_2REAL
#else
  #define MPI_PRECISION MPI_DOUBLE_PRECISION
  #define MPI_2PRECISION MPI_2DOUBLE_PRECISION
#endif
#else
  #define MPI_PRECISION MPI_REAL
  #define MPI_2PRECISION MPI_2REAL
#endif


/* MPI_INTEGER is defined in mpi.h */

#define MPI_INTEGER8 MPI_LONG_LONG_INT

#ifdef FUJI
#if MNH_INT == 8 
#define SIZEINTEGER 8
#define SIZEINTEGER8 8
#define SIZELOGICAL 8
#else
#define SIZEINTEGER 4
#define SIZEINTEGER8 8
#define SIZELOGICAL 4
#endif
#if MNH_REAL == 4
#define SIZEPRECISION 4
#define SIZE2PRECISION 8
#else
#define SIZEPRECISION 8 
#define SIZE2PRECISION 16 
#endif
#else
#define SIZEINTEGER 8 
#define SIZEINTEGER8 8
#define SIZEPRECISION 8
#define SIZE2PRECISION 16
#endif

#define MPI_DOUBLEDOUBLE 999
#define SIZE_DOUBLEDOUBLE SIZE2PRECISION

#if MNH_INT == 8 
#define int long long
#endif
 
static int initflag = 0; /* ADDON Didier */

void disppass(fct)
char *fct;
{
  /* printf("MPIVIDE::Passage dans %s \n", fct); */
}

#pragma weak mpi_cart_sub__ = mpi_cart_sub
#pragma weak mpi_cart_sub_  = mpi_cart_sub
void mpi_cart_sub
( comm, remains_dims, comm_new, __ierr )
int              *comm;
int              *remains_dims;
int              *comm_new;
int              *__ierr;
{
    disppass("cart_sub");
    *comm_new = MPI_COMM_NULL + 1; 
    *__ierr = 0;
}

#pragma weak mpi_wtime__ =  mpi_wtime
#pragma weak mpi_wtime_  =  mpi_wtime
double mpi_wtime()
{
   return 1.0 ; 
}

#pragma weak mpi_cart_coords__ = mpi_cart_coords
#pragma weak mpi_cart_coords_  = mpi_cart_coords
void mpi_cart_coords
( comm, rank, maxdims, coords, __ierr )
int              *comm;
int              *rank;
int              *maxdims;
int              *coords;
int              *__ierr;
{
    disppass("cart_coords");
    coords[0] = 0 ;
    coords[1] = 0 ;
    *__ierr = 0;
}

#pragma weak mpi_cart_create__ = mpi_cart_create
#pragma weak mpi_cart_create_  = mpi_cart_create
void mpi_cart_create
( comm_old, ndims, dims, l_periods, l_reorder, comm_cart, __ierr )
int              *comm_old;
int              *ndims;
int              *dims;
int              *l_periods;
int              *l_reorder;
int              *comm_cart;
int              *__ierr;
{
    disppass("cart_create");
    *comm_cart = MPI_COMM_NULL + 1; 
    *__ierr = 0;
}


#pragma weak mpi_alltoallv__ = mpi_alltoallv
#pragma weak mpi_alltoallv_  = mpi_alltoallv
void mpi_alltoallv(void *sendbuf, int *sendcounts,
            int *sdispls, int *sendtype,
            void *recvbuf, int *recvcounts,
            int *rdispls, int *recvtype, int *comm, int *__ierr)
{
    int size = SIZE2PRECISION;
    disppass("alltoallv");
    switch(*sendtype)
    {
      case MPI_INTEGER:
        size = SIZEINTEGER;
        break;
      case MPI_PRECISION:
        size = SIZEPRECISION;
        break;
      case MPI_2PRECISION:
        size = SIZE2PRECISION;
        break;
      case MPI_DOUBLEDOUBLE:
        size = SIZE_DOUBLEDOUBLE;
        break;
      case MPI_LOGICAL:
        size = SIZELOGICAL ;
        break;
    }
    memcpy(recvbuf, sendbuf, (*recvcounts)*size);

    *__ierr = 0;
}

#pragma weak mpi_comm_split__ = mpi_comm_split
#pragma weak mpi_comm_split_  = mpi_comm_split
void mpi_comm_split(int *comm, int color, int key,
		   int *newcomm, int *__ierr )
{
    disppass("comm_split");
    *newcomm = MPI_COMM_NULL + 2; 
    *__ierr = 0;
}

#pragma weak mpi_allgatherv__ =  mpi_allgatherv
#pragma weak mpi_allgatherv_  =  mpi_allgatherv
void mpi_allgatherv
( sendbuf, sendcount,  sendtype, recvbuf, recvcounts, displs,   recvtype, comm, __ierr )
void             *sendbuf;
int              *sendcount;
int              *sendtype;
void             *recvbuf;
int              *recvcounts;
int              *displs;
int              *recvtype;
int              *comm;
int              *__ierr;
{
   int size = SIZE2PRECISION;
   disppass("allgatherv");

   switch(*sendtype)
    {
      case MPI_INTEGER:
        size = SIZEINTEGER;
        break;
      case MPI_INTEGER8:
        size = SIZEINTEGER8;
        break;
      case MPI_PRECISION:
        size = SIZEPRECISION;
        break;
      case MPI_2PRECISION:
        size = SIZE2PRECISION;
        break;
      case MPI_DOUBLEDOUBLE:
        size = SIZE_DOUBLEDOUBLE;
        break;
      case MPI_LOGICAL:
        size = SIZELOGICAL ;
        break;
    }
    memcpy(recvbuf, sendbuf, (*recvcounts)*size);
    *__ierr = 0;
}

#pragma weak mpi_gather__  = mpi_gather
#pragma weak mpi_gather_   = mpi_gather
void mpi_gather
( sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, __ierr )

void             *sendbuf;
int              *sendcount;
int              *sendtype;
void             *recvbuf;
int              *recvcount;
int              *recvtype;
int              *root;
int              *comm;
int              *__ierr;
{
    int size = SIZE2PRECISION;
    disppass("gather");
    switch(*sendtype)
    {
      case MPI_INTEGER:
        size = SIZEINTEGER;
        break;
      case MPI_PRECISION:
        size = SIZEPRECISION;
        break;
      case MPI_2PRECISION:
        size = SIZE2PRECISION;
        break;
      case MPI_DOUBLEDOUBLE:
        size = SIZE_DOUBLEDOUBLE;
        break;
      case MPI_DOUBLE:
        size = 8 ;
        break;
      case MPI_LOGICAL:
        size = SIZELOGICAL ;
        break;
    }
    memcpy(recvbuf, sendbuf, (*recvcount)*size);

    *__ierr = 0;
}

#pragma weak mpi_gatherv__  = mpi_gatherv
#pragma weak mpi_gatherv_   = mpi_gatherv
void mpi_gatherv
( sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, __ierr )

void             *sendbuf;
int              *sendcount;
int              *sendtype;
void             *recvbuf;
int              *recvcounts;
int              *displs;
int              *recvtype;
int              *root;
int              *comm;
int              *__ierr;
{
    int size = SIZE2PRECISION;
    disppass("gatherv");
    switch(*sendtype)
    {
      case MPI_INTEGER:
        size = SIZEINTEGER;
        break;
      case MPI_PRECISION:
        size = SIZEPRECISION;
        break;
      case MPI_2PRECISION:
        size = SIZE2PRECISION;
        break;
      case MPI_DOUBLEDOUBLE:
        size = SIZE_DOUBLEDOUBLE;
        break;
      case MPI_LOGICAL:
        size = SIZELOGICAL ;
        break;
    }
    memcpy(recvbuf, sendbuf, (*recvcounts)*size);

    *__ierr = 0;
}

#pragma weak mpi_comm_get_parent__ =  mpi_comm_get_parent
#pragma weak mpi_comm_get_parent_  =  mpi_comm_get_parent
void mpi_comm_get_parent(int *parent , int *__ierr)
{
  *parent = 0 ;
  *__ierr = 0;
}

#pragma weak mpi_info_create__ =  mpi_info_create
#pragma weak mpi_info_create_  =  mpi_info_create
void mpi_info_create(int *info , int *__ierr)
{
  *info = 0 ;
  *__ierr = 0;
}

#pragma weak mpi_info_set__ =  mpi_info_set
#pragma weak mpi_info_set_  =  mpi_info_set
void mpi_info_set(int *info ,  char *key, char *value, int *__ierr)
{
  *__ierr = 0;
}

#pragma weak mpi_info_get__ =  mpi_info_get
#pragma weak mpi_info_get_  =  mpi_info_get
void mpi_info_get(int *info ,  char *key, int *valuelen , char *value, int *__ierr)
{
  *__ierr = 0;
}

#pragma weak mpi_comm_spawn__ = mpi_comm_spawn
#pragma weak mpi_comm_spawn_  = mpi_comm_spawn
void mpi_comm_spawn(char *command, char *argv[], int *maxprocs,
            int* info, int *root, int *comm,
		    int *intercomm, int array_of_errcodes[] ,  int *__ierr )
{
  *__ierr = 0;
}

#pragma weak mpi_intercomm_merge__ = mpi_intercomm_merge
#pragma weak mpi_intercomm_merge_  = mpi_intercomm_merge
void mpi_intercomm_merge(int *intercomm, int *high, int *newintracomm ,  int *__ierr )
{
  *__ierr = 0;
}

/* ############################################################################################
############################################################################################
############################################################################################ */


#pragma weak mpi_allgather__  = mpi_allgather
#pragma weak mpi_allgather_   = mpi_allgather
void mpi_allgather
( sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, __ierr )
void             *sendbuf;
int              *sendcount;
int              *sendtype;
void             *recvbuf;
int              *recvcount;
int              *recvtype;
int              *comm;
int              *__ierr;
{
    int size = SIZE2PRECISION;
    disppass("allgather");
    switch(*sendtype)
    {
      case MPI_INTEGER:
        size = SIZEINTEGER;
        break;
      case MPI_PRECISION:
        size = SIZEPRECISION;
        break;
      case MPI_2PRECISION:
        size = SIZE2PRECISION;
        break;
      case MPI_DOUBLEDOUBLE:
        size = SIZE_DOUBLEDOUBLE;
        break;
      case MPI_LOGICAL:
        size = SIZELOGICAL ;
        break;
    }
    memcpy(recvbuf, sendbuf, (*recvcount)*size);

    *__ierr = 0;
}

#pragma weak mpi_op_create__ =  mpi_op_create
#pragma weak mpi_op_create_  =  mpi_op_create
void mpi_op_create
( function, commute, op, __ierr )
void             *function;
int              *commute;
int              *op;
int              *__ierr;
{
    *__ierr = 0;
}

#pragma weak mpi_type_contiguous__ = mpi_type_contiguous
#pragma weak mpi_type_contiguous_  = mpi_type_contiguous
void mpi_type_contiguous
( count, oldtype, newtype, __ierr )
int              *count;
int              *oldtype;
int              *newtype;
int              *__ierr;
{
  *newtype =  MPI_DOUBLEDOUBLE;
   *__ierr = 0;
}

#pragma weak mpi_type_commit__ = mpi_type_commit
#pragma weak mpi_type_commit_  = mpi_type_commit
void mpi_type_commit
( newtype, __ierr )
int              *newtype;
int              *__ierr;
{
  *newtype =  MPI_DOUBLEDOUBLE;
   *__ierr = 0;
}

#pragma weak mpi_allreduce__ =  mpi_allreduce
#pragma weak mpi_allreduce_  =  mpi_allreduce
void mpi_allreduce
( sendbuf, recvbuf, count, datatype, op, comm, __ierr )
void             *sendbuf;
void             *recvbuf;
int              *count;
int              *datatype;
int              *comm;
int              *op;
int              *__ierr;
{
    int size = SIZE2PRECISION;
 
    disppass("allreduce"); 
    switch(*datatype)  
    {
      case MPI_INTEGER:
        size = SIZEINTEGER;
        break;
      case MPI_PRECISION:
        size = SIZEPRECISION;
        break;
      case MPI_2PRECISION:
        size = SIZE2PRECISION;
        break;
      case MPI_DOUBLE:
        size = 8 ;
        break;
      case MPI_LOGICAL:
        size = SIZELOGICAL ;
        break;
    }
    memcpy(recvbuf, sendbuf, (*count)*size); 
    *__ierr = 0;
}

#pragma weak mpi_bcast__ =  mpi_bcast
#pragma weak mpi_bcast_  =  mpi_bcast
void mpi_bcast
( buffer, count, datatype, root, comm, __ierr )
void             *buffer;
int              *count;
int              *datatype;
int              *root;
int              *comm;
int              *__ierr;
{
    disppass("bcast");
    *__ierr = 0;
}

#pragma weak mpi_comm_group__ =  mpi_comm_group
#pragma weak mpi_comm_group_  =  mpi_comm_group
void mpi_comm_group
( comm, group, __ierr )
int              *comm;
int              *group;
int              *__ierr;
{
    disppass("comm_group");
    *__ierr = 0;
}


#pragma weak mpi_group_incl__ = mpi_group_incl
#pragma weak mpi_group_incl_  = mpi_group_incl
void mpi_group_incl
( group, n, ranks, group_out, __ierr )
int              *group, *group_out;
int              *n, *ranks;
int              *__ierr;
{
    disppass("group_incl"); 
    *__ierr = 0;
}


#pragma weak mpi_comm_create__ = mpi_comm_create
#pragma weak mpi_comm_create_  = mpi_comm_create
void mpi_comm_create
( comm, group, comm_out, __ierr )
int              *comm;
int              *group;
int              *comm_out;
int              *__ierr;
{
    disppass("comm_create");
    *comm_out = MPI_COMM_NULL + 1; 
    *__ierr = 0;
}

#pragma weak mpi_comm_free__ = mpi_comm_free
#pragma weak mpi_comm_free_  = mpi_comm_free
void mpi_comm_free
( comm, __ierr )
int              *comm;
int              *__ierr;
{
    disppass("comm_free");
    *__ierr = 0;
}

#pragma weak mpi_group_free__ = mpi_group_free
#pragma weak mpi_group_free_  = mpi_group_free
void mpi_group_free
( group, __ierr )
int              *group;
int              *__ierr;
{
    disppass("group_free");
    *__ierr = 0;
}


#pragma weak mpi_group_excl__ = mpi_group_excl
#pragma weak mpi_group_excl_  = mpi_group_excl
void mpi_group_excl
( group, n, ranks, newgroup, __ierr )
int              *group, *newgroup;
int              *n, *ranks;
int              *__ierr;
{
    disppass("group_excl");
    *__ierr = 0;
}

#pragma weak mpi_comm_rank__ = mpi_comm_rank
#pragma weak mpi_comm_rank_  = mpi_comm_rank
void mpi_comm_rank
( comm, rank, __ierr )
int              *comm;
int              *rank;
int              *__ierr;
{
    disppass("comm_rank");
    *rank = 0;
    *__ierr = 0;
}

#pragma weak mpi_comm_compare__ = mpi_comm_compare
#pragma weak mpi_comm_compare_  = mpi_comm_compare
void mpi_comm_compare
( comm1, comm2, result, __ierr )
int  *comm1;
int  *comm2;
int       *result;
int *__ierr;
{
    disppass("comm_compare");
    *result = MPI_CONGRUENT; /* ADDON Didier */
    *__ierr = 0;
}

#pragma weak mpi_comm_dup__ = mpi_comm_dup
#pragma weak mpi_comm_dup_  = mpi_comm_dup
void mpi_comm_dup
( comm, comm_out, __ierr )
int *comm, *comm_out;
int *__ierr;
{
    disppass("comm_dup");
    *comm_out = *comm;   /* ADDON Didier */
    *__ierr = 0;
}

#pragma weak mpi_comm_size__ = mpi_comm_size
#pragma weak mpi_comm_size_  = mpi_comm_size
void mpi_comm_size
( comm, size, __ierr )
int *comm;
int *size;
int *__ierr;
{
    disppass("comm_size");
    *size = 1;
    *__ierr = 0;
}


#pragma weak mpi_send__ = mpi_send
#pragma weak mpi_send_  = mpi_send
void mpi_send
( buf, count, datatype, dest, tag, comm, __ierr )
void             *buf;
int *count,*dest,*tag;
int     *datatype;
int         *comm;
int *__ierr;
{
    disppass("send");
    *__ierr = 0;
}
#pragma weak mpi_bsend__ = mpi_bsend
#pragma weak mpi_bsend_  = mpi_bsend
void mpi_bsend
( buf, count, datatype, dest, tag, comm, __ierr )
void             *buf;
int *count,*dest,*tag;
int     *datatype;
int         *comm;
int *__ierr;
{
    disppass("bsend");
    *__ierr = 0;
}

#pragma weak mpi_buffer_attach__ = mpi_buffer_attach
#pragma weak mpi_buffer_attach_  = mpi_buffer_attach
void mpi_buffer_attach
( )
{
disppass("buffer_attach");
}

/* JUAN void fclose( iunit )
int *iunit;
{

printf("fclose juanito unit=%d\n",*iunit);
} */


#pragma weak mpi_probe__ = mpi_probe
#pragma weak mpi_probe_  = mpi_probe
void mpi_probe
( source, tag, comm, status, __ierr )
int *source;
int *tag;
int    *comm;
int  *status;
int *__ierr;
{
    disppass("probe");
    *__ierr = 0;
}

#pragma weak mpi_iprobe__ = mpi_iprobe
#pragma weak mpi_iprobe_  = mpi_iprobe
void mpi_iprobe
( source, tag, comm, flag, status, __ierr )
int*source;
int*tag;
int         *flag;
MPI_Comm    *comm;
MPI_Status  *status;
int *__ierr;
{
    disppass("iprobe");
    *__ierr = 0;
}

#pragma weak mpi_get_count__ = mpi_get_count
#pragma weak mpi_get_count_  = mpi_get_count
void mpi_get_count
( status, datatype, count, __ierr )
int   *status;
int *datatype;
int          *count;
int *__ierr;
{
    disppass("get_count");
    *__ierr = 0;
}

#pragma weak mpi_recv__ = mpi_recv
#pragma weak mpi_recv_  = mpi_recv
void mpi_recv
( buf, count, datatype, source, tag, comm, status, __ierr )
void             *buf;
int              *count,*source,*tag;
int     *datatype;
int         *comm;
int       *status;
int *__ierr;
{
    disppass("recv");
    *__ierr = 0;
}

#pragma weak mpi_initialized__ = mpi_initialized
#pragma weak mpi_initialized_  = mpi_initialized
void mpi_initialized
( flag, __ierr )
int  *flag;
int *__ierr;
{
    int lflag;
    disppass("initialized");
    *__ierr = 0;
    *flag = initflag; /* ADDON Didier */
    /*  *flag = MPIR_TO_FLOG(lflag); */
}

#pragma weak mpi_init__ = mpi_init
#pragma weak mpi_init_  = mpi_init
void mpi_init
( ierr )
int *ierr;
{
  disppass("init");
  initflag = 1; /* ADDON Didier */
  *ierr = 0;
}

#pragma weak mpi_abort__ = mpi_abort
#pragma weak mpi_abort_  = mpi_abort
void mpi_abort
( comm, errorcode, __ierr )
int         *comm;
int *errorcode;
int *__ierr;
{
    disppass("abort");
    *__ierr = 0;
}

#pragma weak mpi_finalize__ = mpi_finalize
#pragma weak mpi_finalize_  = mpi_finalize
void mpi_finalize
(__ierr )
int *__ierr;
{
    disppass("finalize");
    *__ierr = 0;
}

#pragma weak mpi_isend__ = mpi_isend
#pragma weak mpi_isend_  = mpi_isend
void mpi_isend
( buf, count, datatype, dest, tag, comm, request, __ierr )
void             *buf;
int*count;
MPI_Datatype     *datatype;
int*dest;
int*tag;
MPI_Comm         *comm;
MPI_Request      *request;
int *__ierr;
{
    disppass("isend");
    *__ierr = 0;
}

#pragma weak mpi_irecv__ = mpi_irecv
#pragma weak mpi_irecv_  = mpi_irecv
void mpi_irecv
( buf, count, datatype, source, tag, comm, request, __ierr )
void             *buf;
int              *count;
MPI_Datatype     * datatype;
int              *source;
int              *tag;
MPI_Comm         *comm;
MPI_Request      *request;
int              *__ierr;
{
    disppass("irecv");
    *__ierr = 0;
}

#pragma weak mpi_testany__ = mpi_testany
#pragma weak mpi_testany_  = mpi_testany
void mpi_testany
( count, array_of_requests, index, flag, status, __ierr )
int*count;
MPI_Request array_of_requests[];
int         *index, *flag;
MPI_Status  *status;
int *__ierr;
{
    disppass("testany");
    *__ierr = 0;
}

#pragma weak mpi_waitall__ = mpi_waitall
#pragma weak mpi_waitall_  = mpi_waitall
void mpi_waitall
(count, array_of_requests, array_of_statuses, __ierr )
int*count;
MPI_Request array_of_requests[];
MPI_Status  array_of_statuses[];
int *__ierr;
{
    disppass("waitall");
    *__ierr = 0;
}

#pragma weak mpi_barrier__ = mpi_barrier
#pragma weak mpi_barrier_  = mpi_barrier
void mpi_barrier
( comm, __ierr )
MPI_Comm *comm;
int *__ierr;
{
    disppass("barrier");
    *__ierr = 0;
}

