/*
 *  $Id$
 *
 *  (C) 1993 by Argonne National Laboratory and Mississipi State University.
 *      All rights reserved.  See COPYRIGHT in top-level directory.
 */

/* mpir version of device interface */

/* #include "mpir.h" */

#ifndef _DMPI_INCLUDE
#define _DMPI_INCLUDE

#define  DMPI_mark_send_completed(DMPI_send_handle) \
          (DMPI_send_handle)->completer = 0;
#define  DMPI_mark_recv_completed(DMPI_recv_handle) \
          (DMPI_recv_handle)->completer = 0;
#define  DMPI_Clr_send_completed(DMPI_send_handle) \
          (DMPI_send_handle)->completer = 1;
#define  DMPI_Clr_recv_completed(DMPI_recv_handle) \
          (DMPI_recv_handle)->completer = 1;
#define DMPI_mpid_recv_handle_from_rhandle( a, b ) b = &((a)->dev_rhandle)

#define DMPI_search_unexpected_queue(src,tag,ctxt_id,found,flag,unex) \
MPID_THREAD_LOCK(0,0);\
if (MPIR_unexpected_recvs.first) {\
MPIR_search_unexpected_queue( src, tag, ctxt_id, found, flag, unex );}\
else *(found) = 0;\
MPID_THREAD_UNLOCK(0,0);
 
/* Defined but do nothing */
#define DMPI_check_mpi(blocking)
#endif


