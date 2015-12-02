/* Parallel computation of eigenvalues and symmetric tridiagonal 
 * matrix T, given by its diagonal elements D and its super-/sub-
 * diagonal elements E.
 *
 * Copyright (c) 2010, RWTH Aachen University
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
 * August 2010, Version 0.7
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
#include <pthread.h>
#include "mpi.h"
#include "pmrrr.h" 
#include "plarre.h"
#include "global.h"
#include "structs.h" 


#define ONE                1.0
#define HUNDRED          100.0
#define HALF               0.5
#define FOURTH             0.25


static void *eigval_subset_thread_a(void *argin);
static void *eigval_subset_thread_r(void *argin);
static void clean_up_plarre(long double*, long double*, int*, int*, int*);


static
int eigval_approx_proc(proc_t *procinfo, int ifirst, int ilast, 
			   int n, long double *D, long double *E, long double *E2,  
			   int *Windex, int *iblock, long double *gersch, tol_t *tolstruct, 
			   long double *W, long double *Werr, long double *Wgap, long double *work,
			   int *iwork);

static
int eigval_root_proc(proc_t *procinfo, int ifirst, int ilast, 
			   int n, long double *D, long double *E, long double *E2,  
			   int *Windex, int *iblock, long double *gersch, tol_t *tolstruct, 
			   long double *W, long double *Werr, long double *Wgap, long double *work,
			 int *iwork);

static
int eigval_refine_proc(proc_t *procinfo, int ifirst, int ilast, 
			   int n, long double *D, long double *E, long double *E2,  
			   int *Windex, int *iblock, long double *gersch, tol_t *tolstruct, 
			   long double *W, long double *Werr, long double *Wgap, long double *work,
			   int *iwork);

static 
auxarg1_t *create_auxarg1(int, long double*, long double*, long double*, int, int, 
			  int, int, int, int*, long double,  long double,
			  long double*, long double*, long double*, int*, int*);
static 
void retrieve_auxarg1(auxarg1_t*, int*, long double**, long double**, long double**,
		      int*, int*, int*, int*, int*, int**, long double*, 
		      long double*, long double**, long double**, long double**, int**, 
		      int**);
static
auxarg2_t *create_auxarg2(int, long double*, long double*, int, int, long double*,
			      long double*,long double*,int*,long double, long double, long double, long double);
static
void retrieve_auxarg2(auxarg2_t*, int*, long double**, long double**, int*,
			  int*, long double**, long double**, long double**, int**, long double*, long double*, long double*,
		      long double*);

static int cmp(const void*, const void*);




/* Routine to compute eigenvalues */
int plarre(proc_t *procinfo, char *jobz, char *range, in_t *Dstruct, 
	       val_t *Wstruct, tol_t *tolstruct, int *nzp, int *offsetp)
{
  /* input variables */
  int              pid    = procinfo->pid;
  int              nproc  = procinfo->nproc;
  bool             wantZ  = (jobz[0]  == 'V' || jobz[0]  == 'v');
  bool             cntval = (jobz[0]  == 'C' || jobz[0]  == 'c');
  int              n      = Dstruct->n;
  long double *restrict D      = Dstruct->D;
  long double *restrict E      = Dstruct->E;
  int    *restrict isplit = Dstruct->isplit;
  long double           *vl    = Wstruct->vl;
  long double           *vu    = Wstruct->vu;
  int              *il    = Wstruct->il;
  int              *iu    = Wstruct->iu;
  long double *restrict W      = Wstruct->W;
  long double *restrict Werr   = Wstruct->Werr;
  long double *restrict Wgap   = Wstruct->Wgap;
  int    *restrict Windex = Wstruct->Windex;
  int    *restrict iblock = Wstruct->iblock;
  long double *restrict gersch = Wstruct->gersch;

  /* constants */
  int             IZERO = 0,   IONE = 1;
  long double          DZERO = 0.0;

  /* work space */
  long double          *E2;
  long double         *work;
  int             *iwork;

  /* compute geschgorin disks and spectral diameter */
  long double          gl, gu, eold, emax, eabs;

  /* compute splitting points */
  int             bl_begin, bl_end, bl_size;

  /* distribute work among processes */
  int             ifirst, ilast, ifirst_tmp, ilast_tmp;
  int             chunk, isize, iil, iiu;

  /* gather results */
  int             *rcount, *rdispl;

  /* others */
  int             info, i, j, jbl, idummy;
  long double          tmp1, dummy;
  bool             sorted;
  enum range_enum {allrng=1, valrng=2, indrng=3} irange;
  long double          intervals[2];
  int             negcounts[2];
  long double          sigma;

  if (range[0] == 'A' || range[0] == 'a') {
    irange = allrng;
  } else if (range[0] == 'V' || range[0] == 'v') {
    irange = valrng;
  } else if (range[0] == 'I' || range[0] == 'i') {
    irange = indrng;
  } else {
    return(1);
  }

  /* allocate work space */
  E2     = (long double *) malloc(     n * sizeof(long double) );
  assert(E2 != NULL);
  work   = (long double *) malloc(   4*n * sizeof(long double) );
  assert(work != NULL);
  iwork  = (int *)    malloc(   3*n * sizeof(int) );
  assert(iwork != NULL);
  rcount = (int *)    malloc( nproc * sizeof(int) );
  assert(rcount != NULL);
  rdispl = (int *)    malloc( nproc * sizeof(int) );
  assert(rdispl != NULL);

  /* Compute square of off-diagonal elements */
  for (i=0; i<n-1; i++) {
    E2[i] = E[i]*E[i];
  }

  /* compute geschgorin disks and spectral diameter */
  gl     = D[0];
  gu     = D[0];
  eold   =  0.0;
  emax   =  0.0;
  E[n-1] =  0.0;

  for (i=0; i<n; i++) {
    eabs = fabsl(E[i]);
    if (eabs >= emax) emax = eabs;
    tmp1 = eabs + eold;
    gersch[2*i] = D[i] - tmp1;
    gl = fminl(gl, gersch[2*i]);
    gersch[2*i+1] = D[i] + tmp1;
    gu = fmaxl(gu, gersch[2*i+1]);
    eold = eabs;
  }
  /* min. pivot allowed in the Sturm sequence of T */
  tolstruct->pivmin = DBL_MIN * fmaxl((long double) 1.0, emax*emax);
  /* estimate of spectral diameter */
  Dstruct->spdiam = gu - gl;

  /* compute splitting points with threshold "split" */
  xdrra_(&n, D, E, E2, &tolstruct->split, &Dstruct->spdiam,
  	  &Dstruct->nsplit, isplit, &info);
  assert(info == 0);

  if (irange == allrng || irange == indrng) {
    *vl = gl;
    *vu = gu;
  }

  /* set eigenvalue indices in case of all or subset by value has
   * to be computed; thereby convert all problem to subset by index
   * computation */
  if (irange == allrng) {
    *il = 1;
    *iu = n;
  } else if (irange == valrng) {
    intervals[0] = *vl; intervals[1] = *vu;
    
    /* find negcount at boundaries 'vl' and 'vu';
     * needs work of dim(n) and iwork of dim(n) */
    xdebz_(&IONE, &IZERO, &n, &IONE, &IONE, &IZERO,
  	    &DZERO, &DZERO, &tolstruct->pivmin, D, E, E2, &idummy,
  	    intervals, &dummy, &idummy, negcounts, work,
  	    iwork, &info);
    assert(info == 0);
    
    /* update negcounts of whole matrix with negcounts found for block */
    *il = negcounts[0] + 1;
    *iu = negcounts[1];
  }

  if (cntval && irange == valrng) {
    /* clean up and return */
    *nzp = iceil(*iu-*il+1, nproc);
    clean_up_plarre(E2, work, iwork, rcount, rdispl);
    return(0);
  }


  /* loop over unreduced blocks */  
  bl_begin  = 0;
  
  for (jbl=0; jbl<Dstruct->nsplit; jbl++) {
    
    bl_end  = isplit[jbl] - 1;
    bl_size = bl_end - bl_begin + 1;
    
    /* deal with 1x1 block immediately */
    if (bl_size == 1) {
      E[bl_end] = 0.0;
      W[bl_begin]      = D[bl_begin];
      Werr[bl_begin]   = 0.0;
      Werr[bl_begin]   = 0.0;
      iblock[bl_begin] = jbl + 1;
      Windex[bl_begin] = 1;
      bl_begin  = bl_end + 1;
      continue;
    }

    /* Indix range of block */
    iil = 1;
    iiu = bl_size;

    /* each process computes a subset of the eigenvalues of the block */
    ifirst_tmp = iil;
    for (i=0; i<nproc; i++) {
      chunk  = (iiu-iil+1)/nproc + (i < (iiu-iil+1)%nproc);
      if (i == nproc-1) {
	ilast_tmp = iiu;
      } else {
	ilast_tmp = ifirst_tmp + chunk - 1;
	ilast_tmp = imin(ilast_tmp, iiu);
      }
      if (i == pid) {
	ifirst    = ifirst_tmp;
	ilast     = ilast_tmp;
	isize     = ilast - ifirst + 1;
	*offsetp = ifirst - iil;
	*nzp      = isize;
      }
      rcount[i]  = ilast_tmp - ifirst_tmp + 1;
      rdispl[i]  = ifirst_tmp - iil;
      ifirst_tmp = ilast_tmp + 1;
      ifirst_tmp = imin(ifirst_tmp, iiu + 1);
    }
    
    /* approximate eigenvalues of input assigned to process */
    if (isize != 0) {      
      info = eigval_approx_proc(procinfo, ifirst, ilast,
				    bl_size, &D[bl_begin], &E[bl_begin], &E2[bl_begin], 
				    &Windex[bl_begin], &iblock[bl_begin], &gersch[2*bl_begin],
				    tolstruct, &W[bl_begin], &Werr[bl_begin], &Wgap[bl_begin], 
				    work, iwork);
      assert(info == 0);    
    }

    /* compute root representation of block */
    info = eigval_root_proc(procinfo, ifirst, ilast,
				  bl_size, &D[bl_begin], &E[bl_begin], &E2[bl_begin], 
				  &Windex[bl_begin], &iblock[bl_begin], &gersch[2*bl_begin],
				  tolstruct, &W[bl_begin], &Werr[bl_begin], &Wgap[bl_begin], 
				  work, iwork);
    assert(info == 0);    

    /* refine eigenvalues assigned to process w.r.t root */
    if (isize != 0) {
      info = eigval_refine_proc(procinfo, ifirst, ilast,
				    bl_size, &D[bl_begin], &E[bl_begin], &E2[bl_begin], 
				    &Windex[bl_begin], &iblock[bl_begin], &gersch[2*bl_begin],
				    tolstruct, &W[bl_begin], &Werr[bl_begin], &Wgap[bl_begin], 
				    work, iwork);
      assert(info == 0);    
    }
    
    memcpy(work, &W[bl_begin], isize * sizeof(long double) );
    MPI_Allgatherv(work, isize, MPI_LONG_DOUBLE, &W[bl_begin], rcount, rdispl,
		   MPI_LONG_DOUBLE, procinfo->comm);
    
    memcpy(work, &Werr[bl_begin], isize * sizeof(long double) );
    MPI_Allgatherv(work, isize, MPI_LONG_DOUBLE, &Werr[bl_begin], rcount, rdispl,
		   MPI_LONG_DOUBLE, procinfo->comm);
    
    memcpy(iwork, &Windex[bl_begin], isize * sizeof(int) );
    MPI_Allgatherv(iwork, isize, MPI_INT, &Windex[bl_begin], rcount, rdispl,
		   MPI_INT, procinfo->comm);
    
    /* Ensure that within block eigenvalues sorted */
    sorted = false;
    while (sorted == false) {
    	sorted = true;
    	for (j=bl_begin; j < bl_end; j++) {
    	  if (W[j+1] < W[j]) {
    	    sorted = false;
    	    tmp1 = W[j];
    	    W[j] = W[j+1];
    	    W[j+1] = tmp1;
    	    tmp1 = Werr[j];
    	    Werr[j] = Werr[j+1];
    	    Werr[j+1] = tmp1;
    	  }
    	}
    }
    
    /* Set indices index correctly */
    for (j=bl_begin; j <= bl_end; j++)
      iblock[j] = jbl + 1;
    
    /* Recompute gaps within the blocks */
    for (j = bl_begin; j < bl_end; j++) {
      Wgap[j] = fmaxl(0.0, (W[j+1] - Werr[j+1]) - (W[j] + Werr[j]) );
    }
    sigma = E[bl_end];
    Wgap[bl_end] = fmaxl(0.0, (gu - sigma) - (W[bl_end] + Werr[bl_end]) );

    /* Compute UNSHIFTED eigenvalues */
    if (!wantZ) {
      sigma = E[bl_end];
      for (i = bl_begin; i <= bl_end; i++) {
	W[i]   += sigma;
      }
    }
    
    /* Proceed with next block */
    bl_begin  = bl_end  + 1;
  }
  /* end of loop over unreduced blocks */    
  
  /* free memory */
  clean_up_plarre(E2, work, iwork, rcount, rdispl);
  
  return(0);
}
  



/*
 * Free's on allocated memory of plarre routine
 */
static  
void clean_up_plarre(long double *E2, long double *work, int *iwork, 
		     int *rcount, int *rdispl)
{
  free(E2);
  free(work);
  free(iwork);
  free(rcount);
  free(rdispl);
}




static 
int eigval_approx_proc(proc_t *procinfo, int ifirst, int ilast, 
			   int n, long double *D, long double *E, long double *E2,  
			   int *Windex, int *iblock, long double *gersch, tol_t *tolstruct, 
			   long double *W, long double *Werr, long double *Wgap, long double *work,
			   int *iwork)
{
  /* Input parameter */
  int              pid = procinfo->pid;
  int              isize        = ilast-ifirst+1;
  long double       pivmin       = tolstruct->pivmin;

  /* double gl, gu, wl, wu; */
  long double wl, wu;

  /* Tolerances */
  long double bsrtol;

  /* /\* Multithreading *\/ */
  int            nthreads;
  int              max_nthreads = procinfo->nthreads;
  int            iifirst, iilast, chunk;
  pthread_t      *threads;
  pthread_attr_t attr;
  auxarg1_t      *auxarg1;
  void           *status;

  /* Others */
  int    nsplit, *isplit;
  int    info, m, i, j;
  long double dummy;
  
  /* Allocate workspace */
  isplit = (int *) malloc( n * sizeof(int) );
  assert(isplit != NULL);
  threads = (pthread_t *) malloc( max_nthreads * sizeof(pthread_t) );
  assert(threads != NULL);

  /* This is an unreduced block */
  nsplit = 1;
  isplit[0] = n;
  
  if (max_nthreads > 1) {
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
  }

  /* Set tolerance parameters */
  bsrtol = sqrt(LDBL_EPSILON);    


  /* APPROXIMATE EIGENVALUES */

  /* compute approximations of the eigenvalues with muliple threads */
  /* equivalent to: */
  /* dlarrd_("I", "B", &n, &dummy, &dummy, &ifirst, &ilast, gersch, */
  /*         &bsrtol, D, E, E2, &pivmin, &nsplit, isplit, &m, W, Werr, */
  /*         &wl, &wu, iblock, Windex, work, iwork, &info); */
  /* assert(info == 0); */
  /* assert(m == ilast-ifirst+1); */
  
  nthreads = max_nthreads;
  while (nthreads > 1 && isize / nthreads < 2)
    nthreads--;

  if (nthreads > 1) {
    
    /* each threads computes W[iifirst:iilast] and places them in
     * work[0:n-1]; the corresponding errors in work[n:2*n-1];
     * the blocks they belong in iwork[0:n-1]; and their indices in
     * iwork[n:2*n-1]; */
    
    iifirst = ifirst;
    chunk = isize / nthreads;
    for (i=1; i<nthreads; i++) {
      
      iilast = iifirst + chunk - 1;

      auxarg1 = create_auxarg1(n, D, E, E2, ifirst, ilast, iifirst, iilast,
			       nsplit, isplit, bsrtol, pivmin, gersch,
			       &work[0], &work[n], &iwork[n], &iwork[0]);
      
      info = pthread_create(&threads[i], &attr,
			    eigval_subset_thread_a,
			    (void *) auxarg1);
      assert(info == 0);
      
      iifirst = iilast + 1;
    }
    iilast = ilast;

    auxarg1 = create_auxarg1(n, D, E, E2, ifirst, ilast, iifirst, iilast,
			     nsplit, isplit, bsrtol, pivmin, gersch,
			     &work[0], &work[n], &iwork[n], &iwork[0]);
    
    status = eigval_subset_thread_a( (void *) auxarg1 );
    assert(status == NULL);
    
    /* join threads */
    for (i=1; i<nthreads; i++) {
      info = pthread_join(threads[i], &status);
      assert(info == 0 && status == NULL);
    }
    
    /* m counts the numbers of eigenvalues computed by process */
    m = isize;
    for (j=0; j<isize; j++) {
      W[j]      = work[j];
      Werr[j]   = work[j+n];
      iblock[j] = iwork[j];
      Windex[j] = iwork[j+n];
    }
    
  } else {
    /* no multithreaded computation */
    
    xdrrd_("I", "B", &n, &dummy, &dummy, &ifirst, &ilast, gersch,
  	    &bsrtol, D, E, E2, &pivmin, &nsplit, isplit, &m, W, Werr,
  	    &wl, &wu, iblock, Windex, work, iwork, &info);
    assert(info == 0);
    assert(m == ilast-ifirst+1);
  }

  /* clean up */
  free(threads);
  free(isplit);

  if (max_nthreads > 1) {
    pthread_attr_destroy(&attr);
  }

  return(0);
}



static 
int eigval_root_proc(proc_t *procinfo, int ifirst, int ilast, 
			   int n, long double *D, long double *E, long double *E2,  
			   int *Windex, int *iblock, long double *gersch, tol_t *tolstruct, 
			   long double *W, long double *Werr, long double *Wgap, long double *work,
			   int *iwork)
{
  /* Input parameter */
  int              pid = procinfo->pid;
  /* int              isize        = ilast-ifirst+1; */
  long double       pivmin       = tolstruct->pivmin;

  /* Tolerances */
  long double rtl;

  /* Create random vector to perturb rrr, same seed */
  int    two_n = 2*n;
  int    iseed[4] = {1,1,1,1};
  long double *randvec;

  long double isleft, isright, spdiam;
  long double sigma, s1, s2;
  int    sgndef, cnt, negcnt_lft, negcnt_rgt;
  long double tau;

  int    jtry, off_L, off_invD;
  long double Dpivot, Dmax;
  bool   noREP;

  int   info, i, j;
  int   IONE = 1, ITWO = 2;
  long double tmp, tmp1, tmp2;
  long double gl, gu;

  /* Set tolerance parameters (need to be same as in refine function) */
  rtl    = sqrt(LDBL_EPSILON);
  
  /* Allocate workspace */
  randvec = (long double *) malloc( 2*n * sizeof(long double) );
  assert(randvec != NULL);

  /* create random vector to perturb rrr and broadcast it */
  xdrnv_(&ITWO, iseed, &two_n, randvec);
  
  /* store shift of initial RRR, here set to zero */
  E[n-1] = 0.0;

  /* find outer bounds GL, GU for block and spectral diameter */
  gl = D[0];
  gu = D[0];
  for (i = 0; i < n; i++) {
    gl = fminl(gl, gersch[2*i]  );
    gu = fmaxl(gu, gersch[2*i+1]);
  }
  spdiam = gu - gl;
  
  /* find approximation of extremal eigenvalues of the block
   * xdrrk computes one eigenvalue of tridiagonal matrix T
   * tmp1 and tmp2 one hold the eigenvalue and error, respectively */
  xdrrk_(&n, &IONE, &gl, &gu, D, E2,
	  &pivmin, &rtl, &tmp1, &tmp2, &info);
  assert(info == 0);  /* if info=-1 => eigenvalue did not converge */
    
  isleft = fmaxl(gl, tmp1-tmp2 - HUNDRED*LDBL_EPSILON*fabsl(tmp1-tmp2) );
    
  xdrrk_(&n, &n, &gl, &gu, D, E2,
  	    &pivmin, &rtl, &tmp1, &tmp2, &info);
  assert(info == 0);  /* if info=-1 => eigenvalue did not converge */
    
  isright = fminl(gu, tmp1+tmp2 + HUNDRED*LDBL_EPSILON*fabsl(tmp1+tmp2) );
  
  spdiam = isright - isleft;
  
  /* compute negcount at points s1 and s2 */
  s1 = isleft  + HALF   * spdiam;
  s2 = isright - FOURTH * spdiam;  /* not needed currently */

  /* compute negcount at points s1 and s2 */
  /* cnt = number of eigenvalues in (s1,s2] = count_right - count_left
   * negcnt_lft = number of eigenvalues smaller equals than s1
   * negcnt_rgt = number of eigenvalues smaller equals than s2 */
  xdrrc_("T", &n, &s1, &s2, D, E, &pivmin,
	  &cnt, &negcnt_lft, &negcnt_rgt, &info);
  assert(info == 0);
  
  /* if more of the desired eigenvectors are in the left part shift left
   * and the other way around */
  if ( negcnt_lft >= n - negcnt_lft ) {
    /* shift left */
    sigma = isleft;
    sgndef = ONE;
  } else {
    /* shift right */
    sigma = isright;
    sgndef = -ONE;
  }

  /* define increment to perturb initial shift to find RRR
   * with not too much element growth */
  tau = spdiam*LDBL_EPSILON*n + 2.0*pivmin;


  /* try to find initial RRR of block:
   * need work space of 3*n here to store D, L, D^-1 of possible
   * representation:
   * D_try      = work[0  :  n-1]
   * L_try      = work[n  :2*n-1]
   * inv(D_try) = work[2*n:3*n-1] */

  off_L    = n;
  off_invD = 2*n;
    
  for (jtry = 0; jtry < MAX_TRY_RRR; jtry++) {

    Dpivot  = D[0] - sigma;
    work[0] = Dpivot;
    Dmax    = fabsl( work[0] );
    j = 0;

    for (i = 0; i < n-1; i++) {
      work[i+off_invD] = 1.0 / work[i];
      tmp = E[j] * work[i+off_invD];
      work[i+off_L] = tmp;
      Dpivot = (D[j+1] - sigma) - tmp*E[j];
      work[i+1] = Dpivot;
      Dmax = fmaxl(Dmax, fabsl(Dpivot) );
      j++;
    }
      
    /* except representation only if not too much element growth */
    if (Dmax > MAX_GROWTH*spdiam) {
      noREP = true;
    } else {
      noREP = false;
    }
      
    if (noREP == true) {
      /* if all eigenvalues are desired shift is made definite to use DQDS
       * so we should not end here */
      if (jtry == MAX_TRY_RRR-2) {
	if (sgndef == ONE) { /* floating point comparison okay here */
	  sigma = gl - FUDGE_FACTOR*spdiam*LDBL_EPSILON*n
	    - FUDGE_FACTOR*2.0*pivmin;
	} else {
	  sigma = gu + FUDGE_FACTOR*spdiam*LDBL_EPSILON*n
	    + FUDGE_FACTOR*2.0*pivmin;
	}
      } else if (jtry == MAX_TRY_RRR-1) {
	fprintf(stderr,"No initial representation could be found.\n");
	exit(3);
      } else {
	sigma -= sgndef*tau;
	tau   *= 2.0;
	continue;
      }
    } else {   /* found representation */
      break;
    }
  }
  /* end trying to find initial RRR of block */

  /* save initial RRR and corresponding shift */
  memcpy(D, &work[0],  n    * sizeof(long double) );
  memcpy(E, &work[n], (n-1) * sizeof(long double) );
  E[n-1] = sigma;
  /* work[0:4*n-1] can now be used again for anything */

  /* perturb root rrr by small relative amount, first make sure
   * that at least two values are actually disturbed enough,
   * which might not be necessary */
  while( fabsl(randvec[0])*RAND_FACTOR < 1.0 )
    randvec[0] *= 2.0;
  while( fabsl(randvec[n-1])  *RAND_FACTOR < 1.0 )
    randvec[n-1]   *= 2.0;

  for (i=0; i<n-1; i++) {
    D[i] *= 1.0 + DBL_EPSILON*RAND_FACTOR*randvec[i];
    E[i] *= 1.0 + DBL_EPSILON*RAND_FACTOR*randvec[i+n];
  }
  D[n-1] *= 1.0 + DBL_EPSILON*RAND_FACTOR*randvec[n-1];

  /* clean up */
  free(randvec);

  return(0);
}




static 
int eigval_refine_proc(proc_t *procinfo, int ifirst, int ilast, 
			   int n, long double *D, long double *E, long double *E2,  
			   int *Windex, int *iblock, long double *gersch, tol_t *tolstruct, 
			   long double *W, long double *Werr, long double *Wgap, long double *work,
			   int *iwork)
{
  /* Input parameter */
  int              pid = procinfo->pid;
  int              isize        = ilast-ifirst+1;
  long double       pivmin       = tolstruct->pivmin;

  /* double gl, gu, wl, wu; */
  long double gl, gu;

  /* Multithreading */
  int            nthreads;
  int              max_nthreads = procinfo->nthreads;
  int            iifirst, iilast, chunk;
  pthread_t      *threads;
  pthread_attr_t attr;
  auxarg2_t      *auxarg2;
  void           *status;

  /* Others */
  int    nsplit, *isplit;
  long double spdiam;
  int    i_low, i_upp;
  long double sigma;

  int    off_DE2, offset;
  int    rf_begin, rf_end;

  int    info, i;

  /* Allocate space */
  threads = (pthread_t *) malloc( max_nthreads * sizeof(pthread_t) );
  assert(threads != NULL);
  isplit = (int *) malloc( n * sizeof(int) );
  assert(isplit != NULL);

  /* This is an unreduced block */
  nsplit = 1;
  isplit[0] = n;
  
  /* Prepare multi-threading */
  if (max_nthreads > 1) {
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
  }

  /* find outer bounds GL, GU for block and spectral diameter */
  gl = D[0];
  gu = D[0];
  for (i = 0; i < n; i++) {
    gl = fminl(gl, gersch[2*i]  );
    gu = fmaxl(gu, gersch[2*i+1]);
  }
  spdiam = gu - gl;

  /* REFINE EIGENVALUES i_low:i_upp WITH REPECT TO RRR */
  
  i_low = Windex[0];
  i_upp = Windex[isize-1];
  sigma = E[n-1];

  /* calculate gaps */
  for (i=0; i<isize-1; i++) {
    Wgap[i] = fmaxl(0.0, (W[i+1] - Werr[i+1]) - (W[i] + Werr[i]) );
  }
  Wgap[isize-1] = fmaxl(0.0, gu - (W[isize-1] + Werr[isize-1]) );
    
  /* shift eigenvalues to be consistent with dqds
   * and compute eigenvalues of SHIFTED matrix */
  for (i=0; i<isize; i++) {
    W[i]    -= sigma;
    Werr[i] += fabsl(W[i])*LDBL_EPSILON;
  }

  /* work  for sequential xdrrb = work[0:2*n-1]
   * iwork for sequential xdrrb = iwork[0:2*n-1]
   * DE2 = work[2*n:3*n-1] strting at bl_begin */
  off_DE2 = 2*n;
    
  /* compute DE2 at store it in work[bl_begin+2*n:bl_end-1+2*n] */
  for (i=0; i<n; i++) {
    work[i+off_DE2] = D[i]*E[i]*E[i];
  }
    
  nthreads = max_nthreads;
  while (nthreads > 1 && isize/nthreads < 2) {
    nthreads--;
  }

  if (nthreads > 1) {

    rf_begin = 0;
    chunk    = isize / nthreads;
    for (i=1; i<nthreads; i++) {
      
      rf_end = rf_begin + chunk - 1;
            
      auxarg2 = create_auxarg2(n, D,
			       &work[off_DE2],
			       rf_begin, rf_end, W, Werr, Wgap, Windex,
			       tolstruct->rtol1, tolstruct->rtol2,
			       pivmin, spdiam);
      
      info = pthread_create(&threads[i], &attr,
  			      eigval_subset_thread_r,
			    (void *) auxarg2);
      assert(info == 0);
      
      rf_begin = rf_end + 1;
    }
    rf_end = isize-1;

    auxarg2 = create_auxarg2(n, D,
			     &work[off_DE2],
			     rf_begin, rf_end, W, Werr, Wgap, Windex,
			     tolstruct->rtol1, tolstruct->rtol2,
			     pivmin, spdiam);
      
    status = eigval_subset_thread_r( (void *) auxarg2 );
    assert(status == NULL);
    
    /* join threads */
    for (i=1; i<nthreads; i++) {
      info = pthread_join(threads[i], &status);
      assert(info == 0 && status == NULL);
    }
    /* should update gaps at splitting points here, but the gaps
     * will be recomputed anyway */
      
  } else {
    
    offset = i_low-1;
    
    /* refine eigenvalues found by xdrrb for i_low:i_upp */
    xdrrb_(&n, D, &work[off_DE2], &i_low,
	    &i_upp, &tolstruct->rtol1, &tolstruct->rtol2, &offset, W, Wgap, 
	    Werr, work, iwork, &pivmin, &spdiam, &n, &info);
    assert(info == 0);
    /* needs work of dim(2*n) and iwork of dim(2*n) */
  }
  /* xdrrb computes gaps correctly, but not last one;
   * this is ignored since the gaps are recomputed anyway */
  
  /* clean up */
  free(threads);
  free(isplit);
  
  if (max_nthreads > 1) {
    pthread_attr_destroy(&attr);
  }
  
  return(0);
}




static 
void *eigval_subset_thread_a(void *argin)
{
  /* from input argument */
  int    n, il, iu, my_il, my_iu;
  long double *D, *E, *E2, *gersch;
  long double bsrtol, pivmin;
  int    nsplit, *isplit;

  /* others */
  int    info;
  long double dummy1, dummy2;
  int    num_vals;
  long double *W_tmp, *Werr_tmp, *W, *Werr;
  int    *iblock_tmp, *Windex_tmp, *iblock, *Windex;
  long double *work;
  int    *iwork;
  
  retrieve_auxarg1((auxarg1_t *) argin, &n, &D, &E, &E2,
		   &il, &iu, &my_il, &my_iu, &nsplit,
		   &isplit, &bsrtol, &pivmin, &gersch,
		   &W, &Werr, &Windex, &iblock);

  /* allocate memory */
  W_tmp = (long double *) malloc( n * sizeof(long double) );
  assert(W_tmp != NULL);
  
  Werr_tmp = (long double *) malloc( n * sizeof(long double) );
  assert(Werr_tmp != NULL);
  
  Windex_tmp = (int *) malloc( n * sizeof(int) );
  assert(Windex_tmp != NULL);

  iblock_tmp = (int *) malloc( n * sizeof(int) );
  assert(iblock_tmp != NULL);

  work  = (long double *) malloc( 4*n * sizeof(long double) );
  assert (work != NULL);

  iwork = (int *) malloc( 3*n * sizeof(int) );
  assert (iwork != NULL);

  /* compute eigenvalues 'my_il' to 'my_iu', put into temporary arrays */
  xdrrd_("I", "B", &n, &dummy1, &dummy2, &my_il, &my_iu, gersch,
  	  &bsrtol, D, E, E2, &pivmin, &nsplit, isplit, &num_vals,
  	  W_tmp, Werr_tmp, &dummy1, &dummy2, iblock_tmp, Windex_tmp,
  	  work, iwork, &info);

  assert(info == 0);

  /* copy computed values in W, Werr, Windex, iblock (which are work space) */
  memcpy(&W[my_il-il],      W_tmp,      num_vals * sizeof(long double) );
  memcpy(&Werr[my_il-il],   Werr_tmp,   num_vals * sizeof(long double) );
  memcpy(&Windex[my_il-il], Windex_tmp, num_vals * sizeof(int)    );
  memcpy(&iblock[my_il-il], iblock_tmp, num_vals * sizeof(int)    );
  
  free(W_tmp);
  free(Werr_tmp);
  free(Windex_tmp);
  free(iblock_tmp);
  free(work);
  free(iwork);

  return(NULL);
}




static 
auxarg1_t *create_auxarg1(int n, long double *D, long double *E, long double *E2,
			  int il, int iu, int my_il, int my_iu, 
			  int nsplit, int *isplit, long double bsrtol, 
			  long double pivmin, long double *gersch, long double *W, 
			  long double *Werr, int *Windex, int *iblock)
{
  auxarg1_t *arg;

  arg = (auxarg1_t *) malloc( sizeof(auxarg1_t) );
  assert(arg != NULL);

  arg->n       = n;
  arg->D       = D;
  arg->E       = E;
  arg->E2      = E2;
  arg->il      = il;
  arg->iu      = iu;
  arg->my_il   = my_il;
  arg->my_iu   = my_iu;
  arg->nsplit  = nsplit;
  arg->isplit  = isplit;
  arg->bsrtol  = bsrtol;
  arg->pivmin  = pivmin;
  arg->gersch  = gersch;
  arg->W       = W;
  arg->Werr    = Werr;
  arg->Windex  = Windex;
  arg->iblock  = iblock;

  return(arg);
}




static 
void retrieve_auxarg1(auxarg1_t *arg, int *n, long double **D, long double **E,
		      long double **E2, int *il, int *iu, int *my_il, 
		      int *my_iu, int *nsplit, int **isplit, 
		      long double *bsrtol, long double *pivmin, long double **gersch, 
		      long double **W, long double **Werr, int **Windex, 
		      int **iblock)
{
  *n      = arg->n;
  *D      = arg->D;
  *E      = arg->E;
  *E2     = arg->E2;
  *il     = arg->il;
  *iu     = arg->iu;
  *my_il  = arg->my_il;
  *my_iu  = arg->my_iu;
  *nsplit = arg->nsplit;
  *isplit = arg->isplit;
  *bsrtol = arg->bsrtol;
  *pivmin = arg->pivmin;
  *gersch = arg->gersch;
  *W      = arg->W;
  *Werr   = arg->Werr;
  *Windex = arg->Windex;
  *iblock = arg->iblock;

  free(arg);
}




static 
void *eigval_subset_thread_r(void *argin)
{
  /* from input argument */
  int          bl_size, rf_begin, rf_end;
  long double       *D, *DE2;
  long double       rtol1, rtol2, pivmin;
  long double       bl_spdiam;
  val_t        *Wstruct;

  /* others */
  int          info, offset;
  long double       *W, *Werr, *Wgap;
  int          *Windex;
  long double       *work;
  int          *iwork;

  retrieve_auxarg2((auxarg2_t *) argin, &bl_size, &D, &DE2,
		   &rf_begin, &rf_end, &W, &Werr, &Wgap, &Windex, &rtol1, &rtol2,
		   &pivmin, &bl_spdiam);

  /* malloc work space */
  work = (long double *) malloc( 2*bl_size * sizeof(long double) );
  assert(work != NULL);
  
  iwork = (int *)   malloc( 2*bl_size * sizeof(int) );
  assert(iwork != NULL);

  /* special case of only one eigenvalue */
  if (rf_begin == rf_end)
    Wgap[rf_begin] = 0.0;
 
  offset = Windex[rf_begin] - 1;
  
  /* call bisection routine to refine the eigenvalues */
  xdrrb_(&bl_size, D, DE2, &Windex[rf_begin], &Windex[rf_end],
	  &rtol1, &rtol2, &offset, &W[rf_begin], &Wgap[rf_begin],
	  &Werr[rf_begin], work, iwork, &pivmin, &bl_spdiam,
	  &bl_size, &info);
  assert(info == 0);

  /* clean up */
  free(work);
  free(iwork);

  return(NULL);
}




static 
auxarg2_t *create_auxarg2(int bl_size, long double *D, long double *DE2,
			  int rf_begin, int rf_end, long double *W, long double *Werr,
			  long double *Wgap, int *Windex,    
			  long double rtol1, long double rtol2, long double pivmin, 
			  long double bl_spdiam)
{
  auxarg2_t *arg;

  arg = (auxarg2_t *) malloc( sizeof(auxarg2_t) );
  assert(arg != NULL);

  arg->bl_size   = bl_size;
  arg->D         = D;
  arg->DE2       = DE2;
  arg->rf_begin  = rf_begin;
  arg->rf_end    = rf_end;
  arg->W   = W;
  arg->Werr   = Werr;
  arg->Wgap   = Wgap;
  arg->Windex   = Windex;
  arg->rtol1     = rtol1;
  arg->rtol2     = rtol2;
  arg->pivmin    = pivmin;
  arg->bl_spdiam = bl_spdiam;

  return(arg);
}




static 
void retrieve_auxarg2(auxarg2_t *arg, int *bl_size, long double **D,
		      long double **DE2, int *rf_begin, int *rf_end,
		      long double **W, long double **Werr, long double **Wgap, int **Windex, 
		      long double *rtol1, long double *rtol2,
		      long double *pivmin, long double *bl_spdiam)
{
  *bl_size   = arg->bl_size;
  *D         = arg->D;
  *DE2       = arg->DE2;
  *rf_begin  = arg->rf_begin;
  *rf_end    = arg->rf_end;
  *W   = arg->W;
  *Werr   = arg->Werr;
  *Wgap   = arg->Wgap;
  *Windex   = arg->Windex;
  *rtol1     = arg->rtol1;
  *rtol2     = arg->rtol2;
  *pivmin    = arg->pivmin;
  *bl_spdiam = arg->bl_spdiam;

  free(arg);
}



/*
 * Compare function for using qsort() on an array
 * of doubles
 */
static 
int cmp(const void *a1, const void *a2)
{
  long double arg1 = *(long double *)a1;
  long double arg2 = *(long double *)a2;

  if (arg1 < arg2)
    return(-1);
  else
    return(1);
}
