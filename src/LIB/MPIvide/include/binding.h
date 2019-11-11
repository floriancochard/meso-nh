#ifndef __MPI_BINDINGS
#define __MPI_BINDINGS

#include "mpi.h"

#if defined(__STDC__) || defined(__cplusplus) || defined(HAVE_PROTOTYPES)
int MPI_Send(void*, int, MPI_Datatype, int, int, MPI_Comm);
int MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
int MPI_Get_count(MPI_Status *, MPI_Datatype, int *);
int MPI_Bsend(void*, int, MPI_Datatype, int, int, MPI_Comm);
int MPI_Ssend(void*, int, MPI_Datatype, int, int, MPI_Comm);
int MPI_Rsend(void*, int, MPI_Datatype, int, int, MPI_Comm);
int MPI_Buffer_attach( void*, int);
int MPI_Buffer_detach( void*, int*);
int MPI_Isend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Ibsend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Issend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Irsend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Irecv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Wait(MPI_Request *, MPI_Status *);
int MPI_Test(MPI_Request *, int *, MPI_Status *);
int MPI_Request_free(MPI_Request *);
int MPI_Waitany(int, MPI_Request *, int *, MPI_Status *);
int MPI_Testany(int, MPI_Request *, int *, int *, MPI_Status *);
int MPI_Waitall(int, MPI_Request *, MPI_Status *);
int MPI_Testall(int, MPI_Request *, int *, MPI_Status *);
int MPI_Waitsome(int, MPI_Request *, int *, int *, MPI_Status *);
int MPI_Testsome(int, MPI_Request *, int *, int *, MPI_Status *);
int MPI_Iprobe(int, int, MPI_Comm, int *flag, MPI_Status *);
int MPI_Probe(int, int, MPI_Comm, MPI_Status *);
int MPI_Cancel(MPI_Request *);
int MPI_Test_cancelled(MPI_Status *, int *flag);
int MPI_Send_init(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Bsend_init(void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *);
int MPI_Ssend_init(void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *);
int MPI_Rsend_init(void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *);
int MPI_Recv_init(void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *);
int MPI_Start(MPI_Request *);
int MPI_Startall(int, MPI_Request *);
int MPI_Sendrecv(void *, int, MPI_Datatype,int, int, void *, int, 
		 MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
int MPI_Sendrecv_replace(void*, int, MPI_Datatype, 
			 int, int, int, int, MPI_Comm, MPI_Status *);
int MPI_Type_contiguous(int, MPI_Datatype, MPI_Datatype *);
int MPI_Type_vector(int, int, int, MPI_Datatype, MPI_Datatype *);
int MPI_Type_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *);
int MPI_Type_indexed(int, int *, int *, MPI_Datatype, MPI_Datatype *);
int MPI_Type_hindexed(int, int *, MPI_Aint *, MPI_Datatype, MPI_Datatype *);
int MPI_Type_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *);
int MPI_Address(void*, MPI_Aint *);
int MPI_Type_extent(MPI_Datatype, MPI_Aint *);

/* See the 1.1 version of the Standard; I think that the standard is in 
   error; however, it is the standard */
/* int MPI_Type_size(MPI_Datatype, MPI_Aint *size); */
int MPI_Type_size(MPI_Datatype, int *);
int MPI_Type_count(MPI_Datatype, int *);
int MPI_Type_lb(MPI_Datatype, MPI_Aint*);
int MPI_Type_ub(MPI_Datatype, MPI_Aint*);
int MPI_Type_commit(MPI_Datatype *);
int MPI_Type_free(MPI_Datatype *);
int MPI_Get_elements(MPI_Status *, MPI_Datatype, int *);
int MPI_Pack(void* inbuf, int, MPI_Datatype, void *outbuf, 
	     int outsize, int *position,  MPI_Comm);
int MPI_Unpack(void* inbuf, int insize, int *position, void *outbuf, 
	       int , MPI_Datatype, MPI_Comm);
int MPI_Pack_size(int, MPI_Datatype, MPI_Comm, 
		  int *);
int MPI_Barrier(MPI_Comm );
int MPI_Bcast(void*fer, int, MPI_Datatype, int root, 
	      MPI_Comm );
int MPI_Gather(void* , int, MPI_Datatype, 
	       void*, int, MPI_Datatype, 
	       int root, MPI_Comm); 
int MPI_Gatherv(void* , int, MPI_Datatype, 
		void*, int *recvcounts, int *displs, 
		MPI_Datatype, int root, MPI_Comm); 
int MPI_Scatter(void* , int, MPI_Datatype, 
		void*, int, MPI_Datatype, 
		int root, MPI_Comm);
int MPI_Scatterv(void* , int *sendcounts, int *displs, 
		 MPI_Datatype, void*, int, 
		 MPI_Datatype, int root, MPI_Comm);
int MPI_Allgather(void* , int, MPI_Datatype, 
		  void*, int, MPI_Datatype, 
		  MPI_Comm);
int MPI_Allgatherv(void* , int, MPI_Datatype, 
		   void*, int *recvcounts, int *displs, 
		   MPI_Datatype, MPI_Comm);
int MPI_Alltoall(void* , int, MPI_Datatype, 
		 void*, int, MPI_Datatype, 
		 MPI_Comm);
int MPI_Alltoallv(void* , int *sendcounts, int *sdispls, 
		  MPI_Datatype, void*, int *recvcounts, 
		  int *rdispls, MPI_Datatype, MPI_Comm);
int MPI_Reduce(void* , void*, int, 
	       MPI_Datatype, MPI_Op op, int root, MPI_Comm);
int MPI_Op_create(MPI_User_function *, int, MPI_Op *);
int MPI_Op_free( MPI_Op *);
int MPI_Allreduce(void* , void*, int, 
		  MPI_Datatype, MPI_Op op, MPI_Comm);
int MPI_Reduce_scatter(void* , void*, int *recvcounts, 
		       MPI_Datatype, MPI_Op op, MPI_Comm);
int MPI_Scan(void* , void*, int, MPI_Datatype, 
	     MPI_Op op, MPI_Comm );
int MPI_Group_size(MPI_Group group, int *);
int MPI_Group_rank(MPI_Group group, int *rank);
int MPI_Group_translate_ranks (MPI_Group group1, int n, int *ranks1, 
			       MPI_Group group2, int *ranks2);
int MPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result);
int MPI_Comm_group(MPI_Comm, MPI_Group *);
int MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
int MPI_Group_intersection(MPI_Group group1, MPI_Group group2, 
			   MPI_Group *newgroup);
int MPI_Group_difference(MPI_Group group1, MPI_Group group2, 
			 MPI_Group *newgroup);
int MPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
int MPI_Group_excl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
int MPI_Group_range_incl(MPI_Group group, int n, int ranges[][3], 
			 MPI_Group *newgroup);
int MPI_Group_range_excl(MPI_Group group, int n, int ranges[][3], 
			 MPI_Group *newgroup);
int MPI_Group_free(MPI_Group *);
int MPI_Comm_size(MPI_Comm, int *);
int MPI_Comm_rank(MPI_Comm, int *rank);
int MPI_Comm_compare(MPI_Comm, MPI_Comm, int *result);
int MPI_Comm_dup(MPI_Comm, MPI_Comm *newcomm);
int MPI_Comm_create(MPI_Comm, MPI_Group group, MPI_Comm *newcomm);
int MPI_Comm_split(MPI_Comm, int color, int key, MPI_Comm *newcomm);
int MPI_Comm_free(MPI_Comm *comm);
int MPI_Comm_test_inter(MPI_Comm, int *flag);
int MPI_Comm_remote_size(MPI_Comm, int *);
int MPI_Comm_remote_group(MPI_Comm, MPI_Group *);
int MPI_Intercomm_create(MPI_Comm local_comm, int local_leader, 
			 MPI_Comm peer_comm, int remote_leader, 
			 int, MPI_Comm *newintercomm);
int MPI_Intercomm_merge(MPI_Comm intercomm, int high, MPI_Comm *newintracomm);
int MPI_Keyval_create(MPI_Copy_function *copy_fn, 
		      MPI_Delete_function *delete_fn, 
		      int *keyval, void* extra_state);
int MPI_Keyval_free(int *keyval);
int MPI_Attr_put(MPI_Comm, int keyval, void* attribute_val);
int MPI_Attr_get(MPI_Comm, int keyval, void *attribute_val, int *flag);
int MPI_Attr_delete(MPI_Comm, int keyval);
int MPI_Topo_test(MPI_Comm, int *);
int MPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods,
		    int reorder, MPI_Comm *comm_cart);
int MPI_Dims_create(int nnodes, int ndims, int *dims);
int MPI_Graph_create(MPI_Comm, int, int *, int *, int, MPI_Comm *);
int MPI_Graphdims_get(MPI_Comm, int *nnodes, int *nedges);
int MPI_Graph_get(MPI_Comm, int, int, int *, int *);
int MPI_Cartdim_get(MPI_Comm, int *ndims);
int MPI_Cart_get(MPI_Comm, int maxdims, int *dims, int *periods,
		 int *coords);
int MPI_Cart_rank(MPI_Comm, int *coords, int *rank);
int MPI_Cart_coords(MPI_Comm, int rank, int maxdims, int *coords);
int MPI_Graph_neighbors_count(MPI_Comm, int rank, int *nneighbors);
int MPI_Graph_neighbors(MPI_Comm, int rank, int maxneighbors,
			int *neighbors);
int MPI_Cart_shift(MPI_Comm, int direction, int disp, 
		   int *rank_source, int *rank_dest);
int MPI_Cart_sub(MPI_Comm, int *remain_dims, MPI_Comm *newcomm);
int MPI_Cart_map(MPI_Comm, int ndims, int *dims, int *periods, 
		 int *newrank);
int MPI_Graph_map(MPI_Comm, int, int *, int *, int *);
int MPI_Get_processor_name(char *name, int *result_len);
int MPI_Errhandler_create(MPI_Handler_function *function, 
			  MPI_Errhandler *errhandler);
int MPI_Errhandler_set(MPI_Comm, MPI_Errhandler errhandler);
int MPI_Errhandler_get(MPI_Comm, MPI_Errhandler *errhandler);
int MPI_Errhandler_free(MPI_Errhandler *errhandler);
int MPI_Error_string(int errorcode, char *string, int *result_len);
int MPI_Error_class(int errorcode, int *errorclass);
double MPI_Wtime(void);
double MPI_Wtick(void);
#ifndef MPI_Wtime
double PMPI_Wtime(void);
double PMPI_Wtick(void);
#endif
int MPI_Init(int *argc, char ***argv);
int MPI_Finalize(void);
int MPI_Initialized(int *flag);
int MPI_Abort(MPI_Comm, int errorcode);
/* MPI-2 communicator naming functions */
int MPI_Comm_set_name(MPI_Comm, char *);
int MPI_Comm_get_name(MPI_Comm, char **);
#ifdef HAVE_NO_C_CONST
/* Default Solaris compiler does not accept const but does accept prototypes */
int MPI_Pcontrol(int level, ...);
#else
int MPI_Pcontrol(const int level, ...);
#endif

int MPI_NULL_COPY_FN ( MPI_Comm oldcomm, int keyval, void *extra_state, 
		       void *attr_in, void *attr_out, int *flag );
int MPI_NULL_DELETE_FN ( MPI_Comm, int keyval, void *attr, 
			 void *extra_state );
int MPI_DUP_FN ( MPI_Comm, int keyval, void *extra_state, void *attr_in,
		 void *attr_out, int *flag );

/* Here are the bindings of the profiling routines */
#if !defined(MPI_BUILD_PROFILING)
int PMPI_Send(void*, int, MPI_Datatype, int, int, 
	     MPI_Comm);
int PMPI_Recv(void*, int, MPI_Datatype, int, 
	     int, MPI_Comm, MPI_Status *);
int PMPI_Get_count(MPI_Status *, MPI_Datatype, int *);
int PMPI_Bsend(void*, int, MPI_Datatype, int, int, 
	      MPI_Comm);
int PMPI_Ssend(void*, int, MPI_Datatype, int, int, 
	      MPI_Comm);
int PMPI_Rsend(void*, int, MPI_Datatype, int, int, 
	      MPI_Comm);
int PMPI_Buffer_attach( void* buffer, int size);
int PMPI_Buffer_detach( void* buffer, int* size);
int PMPI_Isend(void*, int, MPI_Datatype, int, int, 
	      MPI_Comm, MPI_Request *);
int PMPI_Ibsend(void*, int, MPI_Datatype, int, 
	       int, MPI_Comm, MPI_Request *);
int PMPI_Issend(void*, int, MPI_Datatype, int, 
	       int, MPI_Comm, MPI_Request *);
int PMPI_Irsend(void*, int, MPI_Datatype, int, 
	       int, MPI_Comm, MPI_Request *);
int PMPI_Irecv(void*, int, MPI_Datatype, int, 
	      int, MPI_Comm, MPI_Request *);
int PMPI_Wait(MPI_Request *, MPI_Status *);
int PMPI_Test(MPI_Request *, int *flag, MPI_Status *);
int PMPI_Request_free(MPI_Request *);
int PMPI_Waitany(int, MPI_Request *, int *, MPI_Status *);
int PMPI_Testany(int, MPI_Request *, int *, int *, MPI_Status *);
int PMPI_Waitall(int, MPI_Request *, 
		MPI_Status *);
int PMPI_Testall(int, MPI_Request *, int *flag, 
		MPI_Status *);
int PMPI_Waitsome(int, MPI_Request *, int *, 
		 int *, MPI_Status *);
int PMPI_Testsome(int, MPI_Request *, int *, 
		 int *, MPI_Status *);
int PMPI_Iprobe(int, int, MPI_Comm, int *flag, 
	       MPI_Status *);
int PMPI_Probe(int, int, MPI_Comm, MPI_Status *);
int PMPI_Cancel(MPI_Request *);
int PMPI_Test_cancelled(MPI_Status *, int *flag);
int PMPI_Send_init(void*, int, MPI_Datatype, int, 
		  int, MPI_Comm, MPI_Request *);
int PMPI_Bsend_init(void*, int, MPI_Datatype, int, 
		   int, MPI_Comm, MPI_Request *);
int PMPI_Ssend_init(void*, int, MPI_Datatype, int, 
		   int, MPI_Comm, MPI_Request *);
int PMPI_Rsend_init(void*, int, MPI_Datatype, int, 
		   int, MPI_Comm, MPI_Request *);
int PMPI_Recv_init(void*, int, MPI_Datatype, int, 
		  int, MPI_Comm, MPI_Request *);
int PMPI_Start(MPI_Request *);
int PMPI_Startall(int, MPI_Request *);
int PMPI_Sendrecv(void *, int, MPI_Datatype, 
		 int, int, void *, int, 
		 MPI_Datatype, int, int, 
		 MPI_Comm, MPI_Status *);
int PMPI_Sendrecv_replace(void*, int, MPI_Datatype, 
			 int, int, int, int, 
			 MPI_Comm, MPI_Status *);
int PMPI_Type_contiguous(int, MPI_Datatype, 
			MPI_Datatype *);
int PMPI_Type_vector(int, int, int, 
		    MPI_Datatype, MPI_Datatype *);
int PMPI_Type_hvector(int, int, MPI_Aint, 
		     MPI_Datatype, MPI_Datatype *);
int PMPI_Type_indexed(int, int *, 
		     int *, MPI_Datatype, 
		     MPI_Datatype *);
int PMPI_Type_hindexed(int, int *, 
		      MPI_Aint *, MPI_Datatype, 
		      MPI_Datatype *);
int PMPI_Type_struct(int, int *, 
		    MPI_Aint *, 
		    MPI_Datatype *, MPI_Datatype *);
int PMPI_Address(void*, MPI_Aint *);
int PMPI_Type_extent(MPI_Datatype, MPI_Aint *);

/* See the 1.1 version of the Standard; I think that the standard is in 
   error; however, it is the standard */
/* int PMPI_Type_size(MPI_Datatype, MPI_Aint *); */
int PMPI_Type_size(MPI_Datatype, int *);
int PMPI_Type_count(MPI_Datatype, int *);
int PMPI_Type_lb(MPI_Datatype, MPI_Aint*);
int PMPI_Type_ub(MPI_Datatype, MPI_Aint*);
int PMPI_Type_commit(MPI_Datatype *);
int PMPI_Type_free(MPI_Datatype *);
int PMPI_Get_elements(MPI_Status *, MPI_Datatype, int *);
int PMPI_Pack(void* inbuf, int, MPI_Datatype, void *outbuf, 
	     int outsize, int *position,  MPI_Comm);
int PMPI_Unpack(void*, int, int *, void *, 
	       int, MPI_Datatype, MPI_Comm);
int PMPI_Pack_size(int, MPI_Datatype, MPI_Comm, 
		  int *);
int PMPI_Barrier(MPI_Comm );
int PMPI_Bcast(void* buffer, int, MPI_Datatype, int root, 
	      MPI_Comm );
int PMPI_Gather(void* , int, MPI_Datatype, 
	       void*, int, MPI_Datatype, 
	       int root, MPI_Comm); 
int PMPI_Gatherv(void* , int, MPI_Datatype, 
		void*, int *recvcounts, int *displs, 
		MPI_Datatype, int root, MPI_Comm); 
int PMPI_Scatter(void* , int, MPI_Datatype, 
		void*, int, MPI_Datatype, 
		int root, MPI_Comm);
int PMPI_Scatterv(void* , int *sendcounts, int *displs, 
		 MPI_Datatype, void*, int, 
		 MPI_Datatype, int root, MPI_Comm);
int PMPI_Allgather(void* , int, MPI_Datatype, 
		  void*, int, MPI_Datatype, 
		  MPI_Comm);
int PMPI_Allgatherv(void* , int, MPI_Datatype, 
		   void*, int *recvcounts, int *displs, 
		   MPI_Datatype, MPI_Comm);
int PMPI_Alltoall(void* , int, MPI_Datatype, 
		 void*, int, MPI_Datatype, 
		 MPI_Comm);
int PMPI_Alltoallv(void* , int *sendcounts, int *sdispls, 
		  MPI_Datatype, void*, int *recvcounts, 
		  int *rdispls, MPI_Datatype, MPI_Comm);
int PMPI_Reduce(void* , void*, int, 
	       MPI_Datatype, MPI_Op op, int root, MPI_Comm);
int PMPI_Op_create(MPI_User_function *, int, MPI_Op *);
int PMPI_Op_free( MPI_Op *);
int PMPI_Allreduce(void* , void*, int, 
		  MPI_Datatype, MPI_Op op, MPI_Comm);
int PMPI_Reduce_scatter(void* , void*, int *recvcounts, 
		       MPI_Datatype, MPI_Op op, MPI_Comm);
int PMPI_Scan(void* , void*, int, MPI_Datatype, 
	     MPI_Op op, MPI_Comm );
int PMPI_Group_size(MPI_Group group, int *);
int PMPI_Group_rank(MPI_Group group, int *rank);
int PMPI_Group_translate_ranks (MPI_Group group1, int n, int *ranks1, 
			       MPI_Group group2, int *ranks2);
int PMPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result);
int PMPI_Comm_group(MPI_Comm, MPI_Group *);
int PMPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
int PMPI_Group_intersection(MPI_Group group1, MPI_Group group2, 
			   MPI_Group *newgroup);
int PMPI_Group_difference(MPI_Group group1, MPI_Group group2, 
			 MPI_Group *newgroup);
int PMPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
int PMPI_Group_excl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
int PMPI_Group_range_incl(MPI_Group group, int n, int ranges[][3], 
			 MPI_Group *newgroup);
int PMPI_Group_range_excl(MPI_Group group, int n, int ranges[][3], 
			 MPI_Group *newgroup);
int PMPI_Group_free(MPI_Group *);
int PMPI_Comm_size(MPI_Comm, int *);
int PMPI_Comm_rank(MPI_Comm, int *rank);
int PMPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result);
int PMPI_Comm_dup(MPI_Comm, MPI_Comm *newcomm);
int PMPI_Comm_create(MPI_Comm, MPI_Group group, MPI_Comm *newcomm);
int PMPI_Comm_split(MPI_Comm, int color, int key, MPI_Comm *newcomm);
int PMPI_Comm_free(MPI_Comm *comm);
int PMPI_Comm_test_inter(MPI_Comm, int *flag);
int PMPI_Comm_remote_size(MPI_Comm, int *);
int PMPI_Comm_remote_group(MPI_Comm, MPI_Group *);
int PMPI_Intercomm_create(MPI_Comm local_comm, int local_leader, 
			 MPI_Comm peer_comm, int remote_leader, 
			 int, MPI_Comm *newintercomm);
int PMPI_Intercomm_merge(MPI_Comm intercomm, int high, MPI_Comm *newintracomm);
int PMPI_Keyval_create(MPI_Copy_function *copy_fn, 
		      MPI_Delete_function *delete_fn, 
		      int *keyval, void* extra_state);
int PMPI_Keyval_free(int *keyval);
int PMPI_Attr_put(MPI_Comm, int keyval, void* attribute_val);
int PMPI_Attr_get(MPI_Comm, int keyval, void *attribute_val, int *flag);
int PMPI_Attr_delete(MPI_Comm, int keyval);
int PMPI_Topo_test(MPI_Comm, int *);
int PMPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods,
		    int reorder, MPI_Comm *comm_cart);
int PMPI_Dims_create(int nnodes, int ndims, int *dims);
int PMPI_Graph_create(MPI_Comm, int, int *, int *, int, MPI_Comm *);
int PMPI_Graphdims_get(MPI_Comm, int *nnodes, int *nedges);
int PMPI_Graph_get(MPI_Comm, int, int, int *, int *);
int PMPI_Cartdim_get(MPI_Comm, int *ndims);
int PMPI_Cart_get(MPI_Comm, int maxdims, int *dims, int *periods,
		 int *coords);
int PMPI_Cart_rank(MPI_Comm, int *coords, int *rank);
int PMPI_Cart_coords(MPI_Comm, int rank, int maxdims, int *coords);
int PMPI_Graph_neighbors_count(MPI_Comm, int rank, int *nneighbors);
int PMPI_Graph_neighbors(MPI_Comm, int rank, int maxneighbors,
			int *neighbors);
int PMPI_Cart_shift(MPI_Comm, int direction, int disp, 
		   int *rank_source, int *rank_dest);
int PMPI_Cart_sub(MPI_Comm, int *remain_dims, MPI_Comm *newcomm);
int PMPI_Cart_map(MPI_Comm, int ndims, int *dims, int *periods, 
		 int *newrank);
int PMPI_Graph_map(MPI_Comm, int, int *, int *, int *);
int PMPI_Get_processor_name(char *name, int *result_len);
int PMPI_Errhandler_create(MPI_Handler_function *function, 
			  MPI_Errhandler *errhandler);
int PMPI_Errhandler_set(MPI_Comm, MPI_Errhandler errhandler);
int PMPI_Errhandler_get(MPI_Comm, MPI_Errhandler *errhandler);
int PMPI_Errhandler_free(MPI_Errhandler *errhandler);
int PMPI_Error_string(int errorcode, char *string, int *result_len);
int PMPI_Error_class(int errorcode, int *errorclass);
/* Wtime done above */
int PMPI_Init(int *argc, char ***argv);
int PMPI_Finalize(void);
int PMPI_Initialized(int *flag);
int PMPI_Abort(MPI_Comm, int);
/* MPI-2 communicator naming functions */
int PMPI_Comm_set_name(MPI_Comm, char *);
int PMPI_Comm_get_name(MPI_Comm, char **);
#ifdef HAVE_NO_C_CONST
/* Default Solaris compiler does not accept const but does accept prototypes */
int PMPI_Pcontrol(int level, ...);
#else
int PMPI_Pcontrol(const int level, ...);
#endif
#endif  /* MPI_BUILD_PROFILING */

#else 
extern double MPI_Wtime();
extern double MPI_Wtick();
#ifndef MPI_Wtime
extern double PMPI_Wtime();
extern double PMPI_Wtick();
#endif

extern int MPI_NULL_COPY_FN(), MPI_NULL_DELETE_FN();
#endif

#endif
