/* Copyright (c) 2010, RWTH Aachen University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or 
 * without modification, are permitted provided that the following
 * conditions are met:
 *   * Redistributions of source code must retain the above 
 *     copyright notice, this list of conditions and the following
 *     disclaimer.
 *   * Redistributions in binary form must reproduce the above 
 *     copyright notice, this list of conditions and the following 
 *     disclaimer in the documentation and/or other materials 
 *     provided with the distribution.
 *   * Neither the name of the RWTH Aachen University nor the
 *     names of its contributors may be used to endorse or promote 
 *     products derived from this software without specific prior 
 *     written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL RWTH 
 * AACHEN UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF 
 * SUCH DAMAGE.
 *
 * Coded by Matthias Petschow (petschow@aices.rwth-aachen.de),
 * August 2010, Version 0.6
 *
 * This code was the result of a collaboration between 
 * Matthias Petschow and Paolo Bientinesi. When you use this 
 * code, kindly reference a paper related to this work.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <semaphore.h>
#include "mpi.h"
#include "pmrrr.h"
#include "plarrv.h"
#include "global.h"
#include "queue.h"
#include "counter.h"
#include "rrr.h"
#include "structs.h"
#include "process_task.h"


#define THREE            3.0
#define FOUR             4.0


static inline 
rrr_t* compute_new_rrr(cluster_t *cl, int tid, proc_t *procinfo,
		       val_t *Wstruct, vec_t *Zstruct,
		       tol_t *tolstruct, long double *work, int *iwork);

static inline 
int refine_eigvals(cluster_t *cl, int rf_begin, int rf_end,
		   int tid, proc_t *procinfo,
		   rrr_t *RRR, val_t *Wstruct, vec_t *Zstruct,
		   tol_t *tolstruct, counter_t *num_left, 
		   workQ_t *workQ, long double *work, 
		   int *iwork);

static inline 
int communicate_refined_eigvals(cluster_t *cl, proc_t *procinfo,
				int tid, val_t *Wstruct, rrr_t *RRR);

static inline 
int test_comm_status(cluster_t *cl, val_t *Wstruct);

static inline 
int create_subtasks(cluster_t *cl, int tid, proc_t *procinfo,
		    rrr_t *RRR, val_t *Wstruct, vec_t *Zstruct,
		    workQ_t *workQ,
		    counter_t *num_left);




int PMR_process_c_task(cluster_t *cl, int tid, proc_t *procinfo,
		       val_t *Wstruct, vec_t *Zstruct, 
		       tol_t *tolstruct, workQ_t *workQ, 
		       counter_t *num_left, long double *work, int *iwork)
{
  /* From inputs */
  int   depth      = cl->depth;
  int   left_pid   = cl->left_pid;
  int   right_pid  = cl->right_pid;
  int   pid        = procinfo->pid;
  int   n          = Wstruct->n;

  /* Others */
  rrr_t *RRR;
  int   rf_begin, rf_end;
  int   status;

  /* Protection against infinitely deep trees */
  assert(depth < n);

  /* Check if task only need to be split into subtasks */
  if (cl->wait_until_refined == true) {
    status = test_comm_status(cl, Wstruct);
    if (status == COMM_COMPLETE) {
      create_subtasks(cl, tid, procinfo, cl->RRR, Wstruct, Zstruct,
		      workQ, num_left);
      return(C_TASK_PROCESSED);
    } else {
      return(C_TASK_NOT_PROCESSED);
    }
  }

  /* Otherwise: compute new rrr, refine part own cluster,
   * communicate the refined eigenvalues if necessary,
   * and create subtasks if possible */

  RRR = compute_new_rrr(cl, tid, procinfo, Wstruct, Zstruct,
			tolstruct, work, iwork);

  /* Refine eigenvalues 'rf_begin' to 'rf_end' */
  if (left_pid != right_pid) {
    rf_begin = imax(cl->begin, cl->proc_W_begin);
    rf_end   = imin(cl->end,   cl->proc_W_end);
  } 
  if (pid == left_pid ) rf_begin = cl->begin;
  if (pid == right_pid) rf_end   = cl->end;

  refine_eigvals(cl, rf_begin, rf_end, tid, procinfo, RRR,
		 Wstruct, Zstruct, tolstruct, num_left,
		 workQ, work, iwork);

  /* Communicate results: non-blocking */
  status = COMM_COMPLETE;
  if (left_pid != right_pid) {

    status = communicate_refined_eigvals(cl, procinfo, tid,
					 Wstruct, RRR);
    /* status = COMM_INCOMPLETE if communication not finished */
  }

  if (status == COMM_COMPLETE) {
    
    create_subtasks(cl, tid, procinfo, RRR, Wstruct, Zstruct,
		    workQ, num_left);

    return(C_TASK_PROCESSED);
  } else {
    return(C_TASK_NOT_PROCESSED);
  }

} /* end process_c_task */





static inline 
rrr_t* compute_new_rrr(cluster_t *cl, int tid, proc_t *procinfo,
		       val_t *Wstruct, vec_t *Zstruct,
		       tol_t *tolstruct, long double *work, int *iwork)
{
  /* From inputs */
  int              cl_begin    = cl->begin;
  int              cl_end      = cl->end;
  int              cl_size     = cl_end - cl_begin + 1;
  int              depth       = cl->depth;
  int              bl_begin    = cl->bl_begin;
  int              bl_end      = cl->bl_end;
  int              bl_size     = bl_end - bl_begin + 1;
  long double           bl_spdiam   = cl->bl_spdiam;
  rrr_t            *RRR_parent = cl->RRR;

  long double *restrict Werr        = Wstruct->Werr;
  long double *restrict Wgap        = Wstruct->Wgap;
  int    *restrict Windex      = Wstruct->Windex;
  long double *restrict Wshifted    = Wstruct->Wshifted;

  long double           pivmin      = tolstruct->pivmin;

  /* New RRR */
  long double           *restrict D,         *restrict L;
  long double           *restrict DL,        *restrict DLL;
  long double           *restrict D_parent,  *restrict L_parent;
  long double           *DL_parent,          *DLL_parent;
  long double           left_gap, right_gap, tau, fudge;
  rrr_t            *RRR;
  long double           RQtol = 2*LDBL_EPSILON;
  long double           savegap;

  /* Others */
  int              i, k, p, info;
  long double           tmp;
  int              offset, IONE=1;

  /* Allocate memory for new representation for cluster */
  D   = (long double *) malloc(bl_size * sizeof(long double));
  assert(D != NULL);
  
  L   = (long double *) malloc(bl_size * sizeof(long double));
  assert(L != NULL);
  
  DL  = (long double *) malloc(bl_size * sizeof(long double));
  assert(DL != NULL);
  
  DLL = (long double *) malloc(bl_size * sizeof(long double));
  assert(DLL != NULL);

  /* Recompute DL and DLL */
  D_parent = RRR_parent->D;
  L_parent = RRR_parent->L;
  for (i=0; i<bl_size-1; i++) {
    tmp    = D_parent[i]*L_parent[i];
    DL[i]  = tmp;
    DLL[i] = tmp*L_parent[i];
  }
  DL_parent  = DL;
  DLL_parent = DLL;

  /* to shift as close as possible refine extremal eigenvalues */
  for (k=0; k<2; k++) {
    if (k == 0) {
      p              = Windex[cl_begin];
      savegap        = Wgap[cl_begin];
      Wgap[cl_begin] = 0.0;
    } else {
      p              = Windex[cl_end  ];
      savegap        = Wgap[cl_end];
      Wgap[cl_end]   = 0.0;
    }
    
    offset  = Windex[cl_begin] - 1;

    xdrrb_(&bl_size, D_parent, DLL_parent, &p, &p, &RQtol,
	    &RQtol, &offset, &Wshifted[cl_begin], &Wgap[cl_begin],
	    &Werr[cl_begin], work, iwork, &pivmin, &bl_spdiam,
	    &bl_size, &info);
    assert( info == 0 );

    if (k == 0) {
      Wgap[cl_begin] = fmaxl(0, (Wshifted[cl_begin+1]-Werr[cl_begin+1])
			    - (Wshifted[cl_begin]+Werr[cl_begin]) );
    } else {
      Wgap[cl_end]   = savegap;
    }
  } /* end k */

  left_gap  = cl->lgap;
  right_gap = Wgap[cl_end];

  /* Compute new RRR and store it in D and L */
  xdrrf_(&bl_size, D_parent, L_parent, DL_parent,
	  &IONE, &cl_size, &Wshifted[cl_begin], &Wgap[cl_begin],
	  &Werr[cl_begin], &bl_spdiam, &left_gap, &right_gap,
	  &pivmin, &tau, D, L, work, &info);
  assert(info == 0);

  /* Update shift and store it */
  tmp = L_parent[bl_size-1] + tau;
  L[bl_size-1] = tmp;

  /* Compute D*L and D*L*L */
  for (i=0; i<bl_size-1; i++) {
    tmp    = D[i]*L[i];
    DL[i]  = tmp;
    DLL[i] = tmp*L[i];
  }

  /* New RRR of cluster is usually created at the parent level and
   * initialized to parent RRR, now reset to contain new RRR */
  if (RRR_parent->copied_parent_rrr == true) {
    free(RRR_parent->D);
    free(RRR_parent->L);
  }
  RRR = PMR_reset_rrr(RRR_parent, D, L, DL, DLL, bl_size, depth+1);
  
  /* Update shifted eigenvalues */
  for (k=cl_begin; k<=cl_end; k++) {
    fudge  = THREE * LDBL_EPSILON * fabsl( Wshifted[k] );
    Wshifted[k] -= tau;
    fudge += FOUR * LDBL_EPSILON * fabsl( Wshifted[k] );
    Werr[k] += fudge;
  }

  /* Assure that structure is not freed while it is processed */
  PMR_increment_rrr_dependencies(RRR);

  return(RRR);
} /* end compute_new_rrr */





/* 
 * Refine eigenvalues with respect to new rrr 
 */
static inline 
int refine_eigvals(cluster_t *cl, int rf_begin, int rf_end,
		   int tid, proc_t *procinfo, rrr_t *RRR, 
		   val_t *Wstruct, vec_t *Zstruct,
		   tol_t *tolstruct, counter_t *num_left,
		   workQ_t *workQ, long double *work,
		   int *iwork)
{
  /* From inputs */
  int              rf_size   = rf_end-rf_begin+1;
  int              bl_begin  = cl->bl_begin;
  int              bl_end    = cl->bl_end;
  int              bl_size   = bl_end - bl_begin + 1;
  long double           bl_spdiam = cl->bl_spdiam;

  int              nthreads  = procinfo->nthreads;

  long double *restrict D         = RRR->D;
  long double *restrict L         = RRR->L;
  long double *restrict DLL       = RRR->DLL;

  long double *restrict W         = Wstruct->W;
  long double *restrict Werr      = Wstruct->Werr;
  long double *restrict Wgap      = Wstruct->Wgap;
  int    *restrict Windex    = Wstruct->Windex;
  long double *restrict Wshifted  = Wstruct->Wshifted;

  int              nz        = Zstruct->nz;

  long double           pivmin    = tolstruct->pivmin;
  long double           rtol1     = tolstruct->rtol1;
  long double           rtol2     = tolstruct->rtol2;

  /* Others */
  int              info, i, p, q, offset;
  long double           sigma, savegap;
  int              MIN_REFINE_CHUNK = fmax(2,nz/(4*nthreads));
  int              left, own_part, others_part, num_tasks;
  int              ts_begin, ts_end, chunk, count;
  task_t           *task;
  sem_t            sem;
  int              num_iter;

  /* Determine if refinement should be split into tasks */
  left = PMR_get_counter_value(num_left);
  own_part = (int) fmax( ceil( (double) left / nthreads ),
			 MIN_REFINE_CHUNK);

  if (own_part < rf_size) {

    others_part = rf_size - own_part;
    num_tasks   = iceil(rf_size, own_part) - 1; /* >1 */
    chunk       = others_part/num_tasks;        /* floor */

    sem_init(&sem, 0, 0);
    ts_begin = rf_begin;
    p        = Windex[rf_begin];
    for (i=0; i<num_tasks; i++) {
      ts_end = ts_begin + chunk - 1;
      q      = p        + chunk - 1;

      task = PMR_create_r_task(ts_begin, ts_end, D, DLL, p, q, 
			       bl_size, bl_spdiam, tid, &sem);
     
      if (ts_begin <= ts_end)
	PMR_insert_task_at_back(workQ->r_queue, task);
      else
	sem_post(&sem); /* case chunk=0 */

      ts_begin = ts_end + 1;
      p        = q      + 1;
    }
    ts_end = rf_end;
    q      = Windex[rf_end];
    offset = Windex[ts_begin] - 1;

    /* Call bisection routine to refine the values */
    if (ts_begin <= ts_end) {
      xdrrb_(&bl_size, D, DLL, &p, &q, &rtol1, &rtol2, &offset, 
	      &Wshifted[ts_begin], &Wgap[ts_begin], &Werr[ts_begin],
	      work, iwork, &pivmin, &bl_spdiam, &bl_size, &info);
      assert( info == 0 );
    }

    /* Empty "all" r-queue refine tasks before waiting */
    num_iter = PMR_get_num_tasks(workQ->r_queue);
    for (i=0; i<num_iter; i++) {
      task = PMR_remove_task_at_front(workQ->r_queue);
      if (task != NULL) {
	if (task->flag == REFINE_TASK_FLAG) {
	  PMR_process_r_task((refine_t *) task->data, procinfo, 
			     Wstruct, tolstruct, work, iwork);
	  free(task);
	} else {
	  PMR_insert_task_at_back(workQ->r_queue, task);
	}
      } /* if task */
    } /* end for i */
    
    /* Barrier: wait until all created tasks finished */
    count = num_tasks;
    while (count > 0) {
      while (sem_wait(&sem) != 0) { };
      count--;
    }
    sem_destroy(&sem);

    /* Edit right gap at splitting point */
    ts_begin = rf_begin;
    for (i=0; i<num_tasks; i++) {
      ts_end = ts_begin + chunk - 1;
      
      Wgap[ts_end] = fmaxl(0.0, Wshifted[ts_end + 1] - Werr[ts_end + 1]
			  - Wshifted[ts_end] - Werr[ts_end]);
      
      ts_begin = ts_end + 1;
    }

  } else {
    /* Refinement of cluster without creating tasks */
    
    /* 'p' and 'q' are local (within block) indices of
     * the first/last eigenvalue of the cluster */
    p = Windex[rf_begin];
    q = Windex[rf_end];
    
    offset = Windex[rf_begin] - 1;    /* = p - 1 */
    
    if (p == q) {
      savegap = Wgap[rf_begin];
      Wgap[rf_begin] = 0.0;
    }  
    
    /* Bisection routine to refine the values */
    xdrrb_(&bl_size, D, DLL, &p, &q, &rtol1, &rtol2, &offset, 
	    &Wshifted[rf_begin], &Wgap[rf_begin], &Werr[rf_begin],
	    work, iwork, &pivmin, &bl_spdiam, &bl_size, &info);
    assert( info == 0 );
    
    if (p == q) {
      Wgap[rf_begin] = savegap;
    }  
  
  } /* end refine with or without creating tasks */

  sigma     = L[bl_size-1];
  
  /* refined eigenvalues with all shifts applied in W */
  for (i=rf_begin; i<=rf_end; i++) {
    W[i] = Wshifted[i] + sigma;
  }

  return(0);
} /* end refine_eigvals */








static inline 
int communicate_refined_eigvals(cluster_t *cl, proc_t *procinfo,
				int tid, val_t *Wstruct, rrr_t *RRR)
{
  /* From inputs */
  int              cl_begin     = cl->begin;
  int              cl_end       = cl->end;
  int              bl_begin     = cl->bl_begin;
  int              bl_end       = cl->bl_end;
  int              proc_W_begin = cl->proc_W_begin;
  int              proc_W_end   = cl->proc_W_end;
  int              left_pid     = cl->left_pid;
  int              right_pid    = cl->right_pid;
  int              num_messages;

  int              pid          = procinfo->pid;

  long double *restrict W            = Wstruct->W;
  long double *restrict Werr         = Wstruct->Werr;
  long double *restrict Wgap         = Wstruct->Wgap;
  long double *restrict Wshifted     = Wstruct->Wshifted;
  int    *restrict iproc        = Wstruct->iproc;

  /* Others */
  int              p, i_msg, u, k, i;
  int              my_begin, my_end, my_size;
  int              other_begin, other_end, other_size;
  long double           sigma;
  int              status, communication_done;
  MPI_Request      *requests;
  MPI_Status       *stats;
  comm_t           *comm;
  bool             proc_involved;

  my_begin = imax(cl_begin, proc_W_begin);
  my_end   = imin(cl_end,   proc_W_end);
  if (pid == left_pid ) my_begin = cl_begin;
  if (pid == right_pid) my_end   = cl_end;
  my_size  = my_end - my_begin + 1;

  num_messages = 0;
  for (i=left_pid; i<=right_pid; i++) {
    for (k=cl_begin; k<=cl_end; k++) {
      if (iproc[k] == i) {
	num_messages += 4;
	break;
      }
    }    
  }

  requests = (MPI_Request *) malloc( num_messages *
					  sizeof(MPI_Request) );
  stats    = (MPI_Status  *) malloc( num_messages * 
					  sizeof(MPI_Status) );

  i_msg = 0;
  for (p=left_pid; p<=right_pid; p++) {

    proc_involved = false;
    for (k=cl_begin; k<=cl_end; k++) {
      if (iproc[k] == p) {
	proc_involved = true;
	break;
      }
    }

    if (p != pid && proc_involved == true) {

      /* send message to process p (non-blocking) */
      MPI_Isend(&Wshifted[my_begin], my_size, MPI_LONG_DOUBLE, p,
		my_begin, procinfo->comm, &requests[4*i_msg]);

      MPI_Isend(&Werr[my_begin], my_size, MPI_LONG_DOUBLE, p,
		my_begin, procinfo->comm, &requests[4*i_msg+1]);

      /* Find eigenvalues in of process p */
      other_size = 0;
      for (k=cl_begin; k<=cl_end; k++) {
	if (other_size == 0 && iproc[k] == p) {
	  other_begin = k;
	  other_end   = k;
	  other_size++;
	  u = k+1;
	  while (u <=cl_end && iproc[u] == p) {
	    other_end++;
	    other_size++;
	    u++;
	  }
	}
      }
      if (p == left_pid) {
	other_begin = cl_begin;
	u = cl_begin;
	while (iproc[u] == -1) {
	  other_size++;
	  u++;
	}
      }
      if (p == right_pid) {
	other_end = cl_end;
	u = cl_end;
	while (iproc[u] == -1) {
	  other_size++;
	  u--;
	}
      }

      /* receive message from process p (non-blocking) */
      MPI_Irecv(&Wshifted[other_begin], other_size, MPI_LONG_DOUBLE,	p,
		other_begin, procinfo->comm, &requests[4*i_msg+2]);

      MPI_Irecv(&Werr[other_begin], other_size, MPI_LONG_DOUBLE, p,
		other_begin, procinfo->comm, &requests[4*i_msg+3]);
     
      i_msg++;
    }

  } /* end for p */
  num_messages = 4*i_msg; /* messages actually send */

  status = MPI_Testall(num_messages, requests, 
		       &communication_done, stats);
  assert(status == MPI_SUCCESS);

  if (communication_done == true) {

    sigma     = RRR->L[bl_end-bl_begin];
    for (k=cl_begin; k<cl_end; k++) {
      W[k]    = Wshifted[k] + sigma;
      Wgap[k] = fmaxl(0, Wshifted[k+1]-Werr[k+1] 
                        - (Wshifted[k]-Werr[k]));
    }
    W[cl_end] = Wshifted[cl_end] + sigma;

    free(requests);
    free(stats);
    status = COMM_COMPLETE;
  } else {
   
    comm = (comm_t *) malloc( sizeof(comm_t) );
    assert(comm != NULL);

    comm->num_messages      = num_messages;
    comm->requests          = requests;
    comm->stats             = stats;
    cl->wait_until_refined  = true;
    cl->messages            = comm;
    
    status = COMM_INCOMPLETE;
  }
  
  return(status);
} /* end communicate_refined_eigvals */






static inline 
int test_comm_status(cluster_t *cl, val_t *Wstruct)
{
  int         cl_begin            = cl->begin;
  int         cl_end              = cl->end;
  int         bl_begin            = cl->bl_begin;
  int         bl_end              = cl->bl_end;
  rrr_t       *RRR                = cl->RRR;
  comm_t      *comm               = cl->messages;
  int         num_messages        = comm->num_messages;
  MPI_Request *requests           = comm->requests;
  MPI_Status  *stats              = comm->stats;
  long double      *restrict W         = Wstruct->W;
  long double      *restrict Werr      = Wstruct->Werr;
  long double      *restrict Wgap      = Wstruct->Wgap;
  long double      *restrict Wshifted  = Wstruct->Wshifted;

  int         status, k, communication_done;
  long double      sigma;

  /* Test if communication complete */
  status = MPI_Testall(num_messages, requests, 
		       &communication_done, stats);
  assert(status == MPI_SUCCESS);

  if (communication_done == true) {

    cl->wait_until_refined = false;

    sigma     = RRR->L[bl_end-bl_begin];
    for (k=cl_begin; k<cl_end; k++) {
      W[k]    = Wshifted[k] + sigma;
      Wgap[k] = fmaxl(0, Wshifted[k+1]-Werr[k+1] 
		        - (Wshifted[k]-Werr[k]));
    }
    W[cl_end] = Wshifted[cl_end] + sigma;
    
    free(comm);
    free(requests);
    free(stats);
    status = COMM_COMPLETE;
  } else {
    status = COMM_INCOMPLETE;
  }

  return(status);
} /* test_comm_status */





static inline 
int create_subtasks(cluster_t *cl, int tid, proc_t *procinfo, 
		    rrr_t *RRR, val_t *Wstruct, vec_t *Zstruct,
		    workQ_t *workQ, counter_t *num_left)
{
  /* From inputs */
  int              cl_begin  = cl->begin;
  int              cl_end    = cl->end;
  int              depth     = cl->depth;
  int              bl_begin  = cl->bl_begin;
  int              bl_end    = cl->bl_end;
  int              bl_size   = bl_end - bl_begin + 1;
  long double           bl_spdiam = cl->bl_spdiam;
  long double           lgap;

  int              pid       = procinfo->pid;
  int              nproc     = procinfo->nproc;
  int              nthreads  = procinfo->nthreads;
  bool           proc_involved=true;

  long double *restrict Wgap      = Wstruct->Wgap;
  long double *restrict Wshifted  = Wstruct->Wshifted;
  int    *restrict iproc     = Wstruct->iproc;

  int              ldz       = Zstruct->ldz;
  double *restrict Z         = Zstruct->Z;
  int    *restrict Zindex    = Zstruct->Zindex;

  /* others */
  int    i, l, k, j, jj;
  int    max_size;
  task_t *task;
  bool   task_inserted;
  int    new_first, new_last, new_size, new_ftt1, new_ftt2;
  int    sn_first, sn_last, sn_size;
  rrr_t  *RRR_parent;
  int    new_lpid, new_rpid;
  long double *restrict D_parent;
  long double *restrict L_parent;
  int    my_first, my_last;
  bool   copy_parent_rrr;


  max_size = fmax(1, PMR_get_counter_value(num_left) /
		     (fmin(depth+1,4)*nthreads) );
  task_inserted = true;
  new_first = cl_begin;
  for (i=cl_begin; i<=cl_end; i++) {    

    if ( i == cl_end )
      new_last = i;
    else if ( Wgap[i] >= MIN_RELGAP*fabsl(Wshifted[i]) )
      new_last = i;
    else
      continue;

    new_size = new_last - new_first + 1;

    if (new_size == 1) {
      /* singleton was found */
      
      if (new_first==cl_begin || task_inserted==true) {
	/* initialize new singleton task */
	sn_first = new_first;
	sn_last  = new_first;
	sn_size  = 1;
      } else {
	/* extend singleton task by one */
	sn_last++;
	sn_size++;
      }
      
      /* insert task if ... */
      if (i==cl_end || sn_size>=max_size ||
	    Wgap[i+1] < MIN_RELGAP*fabsl(Wshifted[i+1])) {

	/* Check if process involved in s-task */
	proc_involved = false;
	for (k=sn_first; k<=sn_last; k++) {
	  if (iproc[k] == pid) {
	    proc_involved = true;
	    break;
	  }
	}
	if (proc_involved == false) {
	  task_inserted = true;
	  new_first = i + 1;
	  continue;
	}

	/* Insert task as process is involved */
	if (sn_first == cl_begin) {
	  lgap = cl->lgap;
	} else {
	  lgap = Wgap[sn_first-1];
	}
	
	PMR_increment_rrr_dependencies(RRR);
	
	task = PMR_create_s_task(sn_first, sn_last, depth+1, bl_begin,
				 bl_end, bl_spdiam, lgap, RRR);
	
	PMR_insert_task_at_back(workQ->s_queue, task);
	  
	task_inserted = true;
      } else {
	task_inserted = false;
      }
      
    } else {
      /* cluster was found */

      /* check if process involved in processing the new cluster */
      new_lpid = nproc-1;
      new_rpid = -1;
      for (l=new_first; l<=new_last; l++) {
	if (iproc[l] != -1) {
	  new_lpid = imin(new_lpid, iproc[l]);
	  new_rpid = imax(new_rpid, iproc[l]);
	  }
      }
      if (new_lpid > pid || new_rpid < pid) {
	task_inserted = true;
	new_first = i + 1;
	continue;
      }

      /* find gap to the left */
      if (new_first == cl_begin) {
	lgap = cl->lgap;
      } else {
	lgap = Wgap[new_first - 1];
      }
    
      /* determine where to store the parent rrr needed by the
       * cluster to find its new rrr */
      my_first = imax(new_first, cl->proc_W_begin);
      my_last  = imin(new_last,  cl->proc_W_end);
      //      if ( my_first == my_last ) {      
      if ( my_last - my_first + 1 < 4 ) {
	/* less than 4 eigenvalues of cluster belongs to process */
	copy_parent_rrr = true;
      } else {
	/* store parent rrr in Z at column new_ftt */
	copy_parent_rrr = true;
	for (j = my_first; j < my_last; j++) {
	  if (Zindex[j+1] == Zindex[j] + 1) {
	    new_ftt1 = Zindex[j];
	    break;
	  }
	}
	jj = j + 2;
	for (j = jj; j < my_last; j++) {
	  if (Zindex[j+1] == Zindex[j] + 1) {
	    new_ftt2 = Zindex[j];
	    /* Found enough space to store RRR */
	    copy_parent_rrr = false;
	    break;
	  }
	}
      }

      if (copy_parent_rrr == true) {
	/* Copy parent RRR into alloceted arrays and mark them
	 * for freeing later */
	D_parent = (long double *) malloc(bl_size * sizeof(long double));
	assert(D_parent != NULL);
	
	L_parent = (long double *) malloc(bl_size * sizeof(long double));
	assert(L_parent != NULL);

	memcpy(D_parent, RRR->D, bl_size*sizeof(long double));
	memcpy(L_parent, RRR->L, bl_size*sizeof(long double));

	RRR_parent = PMR_create_rrr(D_parent, L_parent, NULL, 
				    NULL, bl_size, depth);
	PMR_set_copied_parent_rrr_flag(RRR_parent, true);

      } else {
	memcpy((long double*) &Z[new_ftt1*ldz], RRR->D, 
	       bl_size*sizeof(long double));
	memcpy((long double*) &Z[new_ftt2*ldz], RRR->L, 
	       bl_size*sizeof(long double));

	RRR_parent = PMR_create_rrr( (long double*) &Z[new_ftt1*ldz],
				    (long double*) &Z[new_ftt2*ldz],
				    NULL, NULL, bl_size, depth);
      }
      
      /* Create the task for the cluster and put it in the queue */ 
      task = PMR_create_c_task(new_first, new_last, depth+1, 
			       bl_begin, bl_end, bl_spdiam, lgap, 
			       cl->proc_W_begin, cl->proc_W_end, 
			       new_lpid, new_rpid, RRR_parent);

      if (new_lpid != new_rpid)
	PMR_insert_task_at_back(workQ->r_queue, task);
      else
	PMR_insert_task_at_back(workQ->c_queue, task);

      task_inserted = true;
      
    } /* if singleton or cluster found */

    new_first = i + 1;
  } /* end i */
  
  /* set flag in RRR that last singleton is created */
  PMR_set_parent_processed_flag(RRR);
  
  /* clean up */
  PMR_try_destroy_rrr(RRR);
  free(cl);

  return(0);
} /* end create_subtasks */
