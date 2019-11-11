/*
 *  $Id$
 *
 *  (C) 1993 by Argonne National Laboratory and Mississipi State University.
 *      All rights reserved.  See COPYRIGHT in top-level directory.
 */

/* user include file for MPI++ programs */

#ifndef _MPIPPP_INCLUDE
#define _MPIPPP_INCLUDE

#ifndef __STDC__
#define __STDC__
#endif

extern "C" {
#include <mpi.h>

// Functions that do not belong to a specific class 
int MPI_Address(void* location, MPI_Aint *address);
}


// Forward class declarations
class MPIX_Group;
class MPIX_Comm;
class MPIX_Comm_world;  // Only instance is MPI_COMM_WORLD 
class MPIX_I_Comm;
class MPIX_Intercomm;
class MPIX_Datatype;
class MPIX_Op;


//----------------------------------------------------------//
// MPIX_Op                                                  //
//----------------------------------------------------------//
class MPIX_Op {
friend class MPIX_Comm;
friend class MPIX_Comm_world;
  
public:
  // Constructor
  MPIX_Op (void) { }
  
  // Constructor
  MPIX_Op (MPI_Op oldop) { op = oldop; }
  
  // Initialization of an op
  Create ( MPI_User_function *function, int commute ) {
	return (MPI_Op_create( function, commute, &op ));
  }
  
  // Destructor
  ~MPIX_Op (void) { }
  
  // Type destructor
  int Free (void) {
      return (MPI_Op_free (&op));
  }
  
private:
  MPI_Op op;

  // Standard operations not provided by this class

  // New op objects cannot created from existing op objects
  // MPIX_Op(const MPIX_Op&) { }

  // Op objects cannot be assigned from other op objects
  // MPIX_Op& operator=(const MPIX_Op&) { }
};


//----------------------------------------------------------//
// MPIX_Datatype                                            //
//----------------------------------------------------------//
class MPIX_Datatype {
friend class MPIX_I_Comm;
friend class MPIX_Comm;
friend class MPIX_Comm_world;

public:
  // Constructor
  MPIX_Datatype (void) { }

  // Constructor
  MPIX_Datatype (MPI_Datatype datatype) { type = datatype; }
  
  // Destructor
  ~MPIX_Datatype (void) { }
  
  // Type initializers
  Contiguous ( int count, MPIX_Datatype& oldtype ) {
	return (MPI_Type_contiguous(count, oldtype.type, &type));
  }
  Vector ( int count, int blocklen, int stride, MPIX_Datatype& oldtype ) {
	return (MPI_Type_vector(count, blocklen, stride, oldtype.type, &type));
  }
  Hvector ( int count, int blocklen, MPI_Aint stride, MPIX_Datatype& oldtype ) {
	return (MPI_Type_hvector(count, blocklen, stride, oldtype.type, &type));
  }
  Indexed ( int count, int blocklens[], int displs[], MPIX_Datatype& oldtype ) {
	return (MPI_Type_indexed ( count, blocklens, displs, oldtype.type, &type));
  }
  Hindexed ( int count, int blocklens[], 
			 MPI_Aint displs[], MPIX_Datatype& oldtype ) {
	return (MPI_Type_hindexed ( count, blocklens, displs, oldtype.type, &type));
  }
  Struct(int count,int blocklens[],MPI_Aint displs[],MPIX_Datatype oldtype[]);
  
  // Type accessors
  Extent ( MPI_Aint& extent ) {
	return (MPI_Type_extent ( type, &extent ));
  }
  Size ( int & size ) {
	return (MPI_Type_size ( type, &size ));
  }
  Lb ( MPI_Aint& lb ) {
	return (MPI_Type_lb ( type, &lb ));
  }
  Ub ( MPI_Aint& ub ) {
	return (MPI_Type_ub ( type, &ub ));
  }
  Commit (void) {
	return (MPI_Type_commit (&type));
  }
  
  // Type destructor
  int Free (void) {
	return (MPI_Type_free (&type));
  }
  
private:
  MPI_Datatype type;

  // Standard operations not provided by this class

  // New datatype objects cannot created from existing datatype objects
  // MPIX_Datatype(const MPIX_Datatype&) { }

  // Datatype objects cannot be assigned from other datatype objects
  // MPIX_Datatype& operator=(const MPIX_Datatype&) { }
};


//----------------------------------------------------------//
// MPIX_Group                                               //
//----------------------------------------------------------//
class MPIX_Group {
friend class MPIX_Comm;
friend class MPIX_I_Comm;
friend class MPIX_Comm_world;

public:
  // Constructor
  MPIX_Group (void) { group = (MPI_Group)0; }

  // Destructor
  ~MPIX_Group (void) { }
  
  // Group accessors
  int Size(int& size_out) { return(MPI_Group_size ( group, &size_out )); }
  int Rank(int& rank_out) { return(MPI_Group_rank ( group, &rank_out)); }
  int Translate_ranks(int n, int r1[], MPIX_Group& group2, int r2[]) { 
	return (MPI_Group_translate_ranks(group, n, r1, group2.group, r2)); 
  }
  
  // Group constructors 
  int Union(MPIX_Group group2, MPIX_Group& group_out) {
	return (MPI_Group_union(group,group2.group,&(group_out.group)));
  }
  int Intersection(MPIX_Group group2, MPIX_Group& group_out) {
	return (MPI_Group_intersection(group,group2.group,&(group_out.group)));
  }
  int Difference(MPIX_Group group2, MPIX_Group& group_out) {
	return (MPI_Group_difference(group,group2.group,&(group_out.group)));
  }
  int Incl(int n, int ranks[], MPIX_Group& group_out) {
	return (MPI_Group_incl(group,n,ranks,&(group_out.group)));
  }
  int Excl(int n, int ranks[], MPIX_Group& group_out) {
	return (MPI_Group_excl(group,n,ranks,&(group_out.group)));
  }
  int Range_incl(int n, int ranges[][3], MPIX_Group& group_out) {
	return (MPI_Group_range_incl(group,n,ranges,&(group_out.group)));
  }
  int Range_excl(int n, int ranges[][3], MPIX_Group& group_out) {
	return (MPI_Group_range_excl(group,n,ranges,&(group_out.group)));
  }
  int Compare(MPIX_Group& group2, int& result) {
	return (MPI_Group_compare(group,group2.group,&result));
  }
  
  // Group destructors
  int Free(void) { return (MPI_Group_free(&group)); }

private:
  MPI_Group group;

  // Operations not provided by this class

  // Initialization of a group from an existing group is not allowed
  // MPIX_Group(const MPIX_Group&) { }

  // Cannot assign one group to another
  // MPIX_Group& operator=(const MPIX_Group&) { }
};


//----------------------------------------------------------//
// MPIX_I_Comm                                              //
//----------------------------------------------------------//
class MPIX_I_Comm {
friend class MPIX_Comm;
friend class MPIX_Comm_world;
  
public:
  // Constructor
  MPIX_I_Comm (void) { comm = (MPI_Comm)0; }
  
  // Destructor
  ~MPIX_I_Comm (void) { }
  
  // Accessors
  int Size(int& size_out) { 
	return(MPI_Comm_size ( comm, &size_out )); 
  }
  int Rank(int& rank_out) { 
	return(MPI_Comm_rank ( comm, &rank_out )); 
  }
  int Group(MPIX_Group& group_out) { 
	return (MPI_Comm_group(comm, &(group_out.group)));
  }
  
  // Initializers
  int Dup(MPIX_I_Comm& comm_out) {
	return (MPI_Comm_dup(comm, &(comm_out.comm)));
  }

  // Destructors
  int Free(void) { return (MPI_Comm_free(&comm)); }


  // Sends
  int Send(void *buf,int cnt, MPIX_Datatype type, int dest, int tag) { 
	return (MPI_Send(buf,cnt,type.type,dest,tag,comm)); 
  }
  int Bsend(void* buf, int cnt, MPIX_Datatype type, int dest, int tag) {
	return (MPI_Bsend(buf,cnt,type.type,dest,tag,comm));
  }
  int Ssend(void* buf, int cnt, MPIX_Datatype type, int dest, int tag) {
	return (MPI_Ssend(buf,cnt,type.type,dest,tag,comm));
  }
  int Rsend(void* buf, int cnt, MPIX_Datatype type, int dest, int tag) {
	return (MPI_Rsend(buf,cnt,type.type,dest,tag,comm));
  }
  int Isend(void *buf, int cnt, MPIX_Datatype type, int dest, int tag,
			MPI_Request& request) { 
	return (MPI_Isend(buf,cnt,type.type,dest,tag,comm,&request)); 
  }
  int Ibsend(void *buf, int cnt, MPIX_Datatype type, int dest, int tag,
			 MPI_Request& request) { 
	return (MPI_Ibsend(buf,cnt,type.type,dest,tag,comm,&request)); 
  }
  int Issend(void *buf, int cnt, MPIX_Datatype type, int dest, int tag,
			 MPI_Request& request) { 
	return (MPI_Issend(buf,cnt,type.type,dest,tag,comm,&request)); 
  }
  int Irsend(void *buf, int cnt, MPIX_Datatype type, int dest, int tag,
			 MPI_Request& request) { 
	return (MPI_Irsend(buf,cnt,type.type,dest,tag,comm,&request)); 
  }
  int Send_init(void* buf, int cnt, MPIX_Datatype type, int dest, 
				int tag, MPI_Request& request) {
	return (MPI_Send_init(buf,cnt,type.type,dest,tag,comm,&request));
  }
  int Bsend_init(void* buf, int cnt, MPIX_Datatype type, int dest, 
				 int tag, MPI_Request& request) {
	return (MPI_Bsend_init(buf,cnt,type.type,dest,tag,comm,&request));
  }
  int Ssend_init(void* buf, int cnt, MPIX_Datatype type, int dest, 
				 int tag, MPI_Request& request) {
	return (MPI_Ssend_init(buf,cnt,type.type,dest,tag,comm,&request));
  }
  int Rsend_init(void* buf, int cnt, MPIX_Datatype type, int dest, 
				 int tag, MPI_Request& request) {
	return (MPI_Rsend_init(buf,cnt,type.type,dest,tag,comm,&request));
  }


  // Receives
  int Recv(void *buf, int cnt, MPIX_Datatype type, int source, 
		   int tag, MPI_Status& status) { 
	return (MPI_Recv(buf,cnt,type.type,source,tag,comm, &status)); 
  }
  int Irecv(void *buf, int cnt, MPIX_Datatype type, int source, 
			int tag, MPI_Request& request) { 
	return (MPI_Irecv(buf,cnt,type.type,source,tag,comm, &request)); 
  }
  int Recv_init(void* buf, int cnt, MPIX_Datatype type, int source, 
				int tag, MPI_Request &request) {
	return (MPI_Recv_init(buf,cnt,type.type,source,tag,comm, &request));
  }

  // Pack and Unpack operartions
  int Pack(void* inbuf, int incnt, MPIX_Datatype type, void *outbuf, 
	       int outsize, int& position) {
	return (MPI_Pack(inbuf,incnt,type.type,outbuf,outsize,&position,comm));
  }
  int Unpack(void* inbuf, int insize, int& position, void *outbuf, 
			 int outcnt, MPIX_Datatype type) {
	return (MPI_Unpack(inbuf,insize,&position,outbuf,outcnt,type.type,comm));
  }
  int Pack_size(int incount, MPIX_Datatype type, int& size) {
	return (MPI_Pack_size(incount, type.type, comm, &size));
  }

  // Completion operations
  int Wait(MPI_Request& request, MPI_Status& status) { 
	return (MPI_Wait(&request, &status)); 
  }
  int Waitall(int count, MPI_Request reqs[], MPI_Status status[]) {
	return (MPI_Waitall(count,reqs,status));
  }
  int Waitany(int count, MPI_Request reqs[], int& index, MPI_Status& status) {
	return (MPI_Waitany(count, reqs, &index, &status));
  }
  int Waitsome(int incnt, MPI_Request reqs[], 
			   int& outcnt, int& index, MPI_Status& status) {
	return (MPI_Waitsome(incnt, reqs, &outcnt, &index, &status));
  }      
  // Test needs to be fixed to be draft compliant
  int Test(MPI_Request& request,int& flag, MPI_Status& status) { 
	return (MPI_Test(&request, &flag, &status)); 
  }
  int Testany(int count, MPI_Request reqs[], int& index, 
			  int& flag, MPI_Status& status) {
	return (MPI_Testany(count, reqs, &index, &flag, &status));
  }
  int Testall(int count, MPI_Request reqs[], int& flag, MPI_Status& status) {
	return (MPI_Testall(count, reqs, &flag, &status));
  }
  int Testsome(int incnt, MPI_Request reqs[], int& outcnt, 
			   int& index, MPI_Status& status) {
	return (MPI_Testsome(incnt, reqs, &outcnt, &index, &status));
  }
  int Test_cancelled(MPI_Status& status, int& flag) {
	return (MPI_Test_cancelled(&status, &flag));
  }
  int Probe(int source, int tag, MPI_Status& status) {
	return (MPI_Probe(source,tag,comm,&status));
  }
  int Iprobe(int source, int tag, int& flag, MPI_Status& status) {
	return (MPI_Iprobe(source,tag,comm,&flag,&status));
  }

  // Need all of the persistent communication request functions
  // Also need the functions from 3.11

private:
  MPI_Comm comm;

  // Standard operations not provided by this class

  // Creating a new communicator from an existing communicator
  // MPIX_Comm(const MPIX_Comm&) { }
  
  // You cannot assign a value to MPIX_Comm's
  // MPIX_Comm& operator=(const MPIX_Comm&) { }
};


//----------------------------------------------------------//
// MPIX_Comm                                                //
//----------------------------------------------------------//
class MPIX_Comm: public MPIX_I_Comm {
public:
  // Constructor
  MPIX_Comm (void) { }

  // Initializers
  int Create(MPIX_Group& group, MPIX_Comm& comm_out) {
	return (MPI_Comm_create(comm, group.group, &(comm_out.comm)));
  }
  int Split(int color, int key, MPIX_Comm& comm_out) {
	return (MPI_Comm_split(comm, color, key, &(comm_out.comm)));
  }


  // Inter-communicator initializer
  int Intercomm_create(int ll, MPIX_Comm& pc, int rl, int tag, 
					 MPIX_Comm& comm_out) {
	return (MPI_Intercomm_create(comm,ll,pc.comm,rl,tag,&(comm_out.comm)));
  }


  // Collective operations
  int Barrier(void) { return (MPI_Barrier(comm)); }
  int Bcast(void *buf, int cnt, MPIX_Datatype type, int root) { 
	return (MPI_Bcast(buf, cnt, type.type, root, comm)); 
  }
  int Gather(void *sbuf,int scnt, MPIX_Datatype stype, void *rbuf,
			 int rcnt, MPIX_Datatype rtype, int root) {
	return (MPI_Gather(sbuf,scnt,stype.type,rbuf,rcnt,rtype.type,root,comm)); 
  }
  int Gatherv(void *sbuf,int scnt, MPIX_Datatype stype, void *rbuf,
			  int *rcnts, int *displs, MPIX_Datatype rtype, int root) {
	return (MPI_Gatherv(sbuf,scnt,stype.type,rbuf,rcnts,displs,rtype.type,root,comm));
  }
  int Scatter(void *sbuf, int scnt, MPIX_Datatype stype, void *rbuf,
			  int rcnt, MPIX_Datatype rtype, int root) { 
	return (MPI_Scatter(sbuf,scnt,stype.type,rbuf,rcnt,rtype.type,root,comm)); 
  }
  int Scatterv(void *sbuf, int *scnts, int *displs, 
			   MPIX_Datatype stype, void *rbuf, int rcnt, 
			   MPIX_Datatype rtype, int root) {
	return (MPI_Scatterv(sbuf,scnts,displs,stype.type,rbuf,rcnt,rtype.type,root,
						 comm)); 
  }
  int Allgather(void *sbuf,int scnt, MPIX_Datatype stype, void *rbuf,
				int rcnt, MPIX_Datatype rtype) { 
	return (MPI_Allgather(sbuf,scnt,stype.type,rbuf,rcnt,rtype.type,comm)); 
  }
  int Allgatherv(void *sbuf,int scnt, MPIX_Datatype stype, void *rbuf,
				 int *rcnts, int *displs, MPIX_Datatype rtype) {
	return (MPI_Allgatherv(sbuf,scnt,stype.type,rbuf,rcnts,displs,rtype.type,comm)); 
  }
  int Alltoall(void *sbuf, int scnt, MPIX_Datatype stype, void *rbuf,
			   int rcnt, MPIX_Datatype rtype) {
	return (MPI_Alltoall(sbuf,scnt,stype.type,rbuf,rcnt,rtype.type,comm)); 
  }
  int Alltoallv(void *sbuf, int *scnts, int *sdispls, 
				MPIX_Datatype stype, void *rbuf, int *rcnts, 
				int *rdispls, MPIX_Datatype rtype) {
	return (MPI_Alltoallv(sbuf,scnts,sdispls,stype.type,rbuf,rcnts,
						  rdispls,rtype.type,comm)); 
  }
  int Reduce(void *in,void *out,int cnt,
			 MPIX_Datatype type, MPIX_Op op, int root) {
	return (MPI_Reduce ( in, out, cnt, type.type, op.op, root, comm )); 
  }
  int Allreduce(void *in, void *out, int cnt, 
				MPIX_Datatype type, MPIX_Op op) {
	return (MPI_Allreduce (in, out, cnt, type.type, op.op, comm)); 
  }
  int Reduce_scatter(void *sbuf, void *rbuf, int *rcnts, 
					 MPIX_Datatype type, MPIX_Op op)	{ 
	return (MPI_Reduce_scatter(sbuf,rbuf,rcnts,type.type,op.op,comm)); 
  }
  int Scan(void *sbuf, void *rbuf, int cnt, 
		   MPIX_Datatype type, MPIX_Op op) {
	return (MPI_Scan(sbuf,rbuf,cnt,type.type,op.op,comm)); 
  }
  
  
  // Topology functions
  int Dims_create(int nnodes, int ndims, int dims[]) {
	return (MPI_Dims_create(nnodes,ndims,dims));
  }
  int Cart_create(int ndims, int dims[], int periods[], int reorder,
				  MPIX_Comm& comm_out) {
	return (MPI_Cart_create(comm,ndims,dims,periods,reorder,&(comm_out.comm)));
  }
  int Cart_coords(int rank, int ndims, int coords[]) {
	return (MPI_Cart_coords(comm,rank,ndims,coords));
  }
  int Cart_rank(int coords[],int &rank) {
	return (MPI_Cart_rank(comm,coords,&rank));
  }
  int Topo_test(int& status) {
	return (MPI_Topo_test(comm,&status));
  }
  int Graph_create(int nnodes, int& index, int& edges,
				   int reorder, MPIX_Comm& comm_graph) {
	return (MPI_Graph_create(comm,nnodes,&index,&edges,reorder,&(comm_graph.comm)));
  }
};


// MPI_COMM_WORLD - specific class
class MPIX_Comm_world : public MPIX_Comm {
private:
  static int init;

public:
  int Init(int& argc, char **argv);
  int Finalize (void) { return (MPI_Finalize()); }
};

#endif


