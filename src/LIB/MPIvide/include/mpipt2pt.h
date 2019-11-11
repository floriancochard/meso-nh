#ifndef _MPI_PT2PT
#define _MPI_PT2PT

#ifndef ANSI_ARGS
#if defined(__STDC__) || defined(__cplusplus) || defined(HAVE_PROTOTYPES)
#define ANSI_ARGS(a) a
#else
#define ANSI_ARGS(a) ()
#endif
#endif

#if defined(USE_STDARG) && \
   (defined(__STDC__) || defined(__cplusplus) || defined(HAVE_PROTOTYPES))
/* This SAME test must be used in pt2pt/mperror.c */
#define MPIR_USE_STDARG
#endif

int MPIR_Pack  ANSI_ARGS(( struct MPIR_COMMUNICATOR *, int, void *, int, 
			   struct MPIR_DATATYPE *, void *, int, int *));
int MPIR_Pack_size  ANSI_ARGS(( int, struct MPIR_DATATYPE *, 
				struct MPIR_COMMUNICATOR *, int, int *));
#ifdef MPI_ADI2
int MPIR_Unpack ANSI_ARGS(( struct MPIR_COMMUNICATOR *, void *, int, int, 
			    struct MPIR_DATATYPE *, 
			    MPID_Msgrep_t, void *, int *, int * ));
#else
int MPIR_Unpack ANSI_ARGS(( struct MPIR_COMMUNICATOR *, void *, int, int, 
			    struct MPIR_DATATYPE *,
			    int, void *, int *, int * ));
#endif
int MPIR_UnPackMessage ANSI_ARGS(( char *, int, MPI_Datatype, int, 
				   MPI_Request, int * ));

/* These are used in the debugging-enabled version */
#ifdef MPI_ADI2
void MPIR_Sendq_init ANSI_ARGS(( void ));
void MPIR_Sendq_finalize ANSI_ARGS(( void ));
void MPIR_Remember_send ANSI_ARGS(( MPIR_SHANDLE *, void *, int, MPI_Datatype,
				    int, int, struct MPIR_COMMUNICATOR * ));
void MPIR_Forget_send ANSI_ARGS(( MPIR_SHANDLE * ));
#endif

/* Datatype service routines */
int MPIR_Type_free ANSI_ARGS(( struct MPIR_DATATYPE ** ));
#ifdef FOO
void MPIR_Type_free_struct ANSI_ARGS(( struct MPIR_DATATYPE * ));
#endif
struct MPIR_DATATYPE *MPIR_Type_dup ANSI_ARGS(( struct MPIR_DATATYPE * ));
int MPIR_Type_permanent ANSI_ARGS(( struct MPIR_DATATYPE * ));
void MPIR_Free_perm_type ANSI_ARGS(( MPI_Datatype ));
void MPIR_Free_struct_internals ANSI_ARGS(( struct MPIR_DATATYPE * ));
void MPIR_Type_get_limits ANSI_ARGS(( struct MPIR_DATATYPE *, MPI_Aint *, MPI_Aint *));
#ifndef MPI_ADI2
int MPIR_Send_init ANSI_ARGS(( void *, int, struct MPIR_DATATYPE *, int, int, 
			       MPI_Comm, MPI_Request, MPIR_Mode, int ));
int MPIR_Recv_init ANSI_ARGS(( void *, int, struct MPIR_DATATYPE *, int, int, 
			       MPI_Comm, MPI_Request, int ));
#endif
extern MPI_Handler_function MPIR_Errors_are_fatal;
extern MPI_Handler_function MPIR_Errors_return;
extern MPI_Handler_function MPIR_Errors_warn;
/* Since these are declared as handler functions, we do not
   redeclare them */
#ifdef MPIR_USE_STDARG
/* gcc requires an explicit declaration when checking for strict prototypes */
void MPIR_Errors_are_fatal ANSI_ARGS(( MPI_Comm *, int *, ... ));
void MPIR_Errors_return ANSI_ARGS(( MPI_Comm *, int *, ... ));
void MPIR_Errors_warn ANSI_ARGS(( MPI_Comm *, int *, ... ));
#else
#ifdef FOO
/* Otherwise, just accept the typedef declaration */
void MPIR_Errors_are_fatal ANSI_ARGS(( MPI_Comm *, int *, char *, 
				      char *, int *));
void MPIR_Errors_return ANSI_ARGS(( MPI_Comm *, int *, char *, char *, int *));
void MPIR_Errors_warn ANSI_ARGS(( MPI_Comm *, int *, char *, char *, int *));
#endif
#endif /* MPIR_USE_STDARG */

#endif
