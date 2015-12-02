/* dlarrr.f -- translated by f2c (version 20061008) */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#define TRUE_ (1)
#define FALSE_ (0)

/* Subroutine */ 
int xdrrr_(int *n, long double *d__, long double *e, 
	int *info)
{
    /* System generated locals */
    int i__1;
    long double d__1;

    /* Builtin functions */
    // long double sqrt(long double);

    /* Local variables */
    int i__;
    long double eps, tmp, tmp2, rmin;
    // extern long double odmch_(char *);
    long double offdig, safmin;
    int yesrel;
    long double smlnum, offdig2;


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */


/*  Purpose */
/*  ======= */

/*  Perform tests to decide whether the symmetric tridiagonal matrix T */
/*  warrants expensive computations which guarantee high relative accuracy */
/*  in the eigenvalues. */

/*  Arguments */
/*  ========= */

/*  N       (input) INT */
/*          The order of the matrix. N > 0. */

/*  D       (input) LONG DOUBLE PRECISION array, dimension (N) */
/*          The N diagonal elements of the tridiagonal matrix T. */

/*  E       (input/output) LONG DOUBLE PRECISION array, dimension (N) */
/*          On entry, the first (N-1) entries contain the subdiagonal */
/*          elements of the tridiagonal matrix T; E(N) is set to ZERO. */

/*  INFO    (output) INT */
/*          INFO = 0(default) : the matrix warrants computations preserving */
/*                              relative accuracy. */
/*          INFO = 1          : the matrix warrants computations guaranteeing */
/*                              only absolute accuracy. */

/*  Further Details */
/*  =============== */

/*  Based on contributions by */
/*     Beresford Parlett, University of California, Berkeley, USA */
/*     Jim Demmel, University of California, Berkeley, USA */
/*     Inderjit Dhillon, University of Texas, Austin, USA */
/*     Osni Marques, LBNL/NERSC, USA */
/*     Christof Voemel, University of California, Berkeley, USA */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     As a default, do NOT go for relative-accuracy preserving computations. */
    /* Parameter adjustments */
    --e;
    --d__;

    /* Function Body */
    *info = 1;
    safmin = LDBL_MIN; // safmin = odmch_("Safe minimum");
    eps = LDBL_EPSILON; // eps = odmch_("Precision");
    smlnum = safmin / eps;
    rmin = sqrt(smlnum);
/*     Tests for relative accuracy */

/*     Test for scaled diagonal dominance */
/*     Scale the diagonal entries to one and check whether the sum of the */
/*     off-diagonals is less than one */

/*     The sdd relative error bounds have a 1/(1- 2*x) factor in them, */
/*     x = max(OFFDIG + OFFDIG2), so when x is close to 1/2, no relative */
/*     accuracy is promised.  In the notation of the code fragment below, */
/*     1/(1 - (OFFDIG + OFFDIG2)) is the condition number. */
/*     We don't think it is worth going into "sdd mode" unless the relative */
/*     condition number is reasonable, not 1/macheps. */
/*     The threshold should be compatible with other thresholds used in the */
/*     code. We set  OFFDIG + OFFDIG2 <= .999 =: RELCOND, it corresponds */
/*     to losing at most 3 decimal digits: 1 / (1 - (OFFDIG + OFFDIG2)) <= 1000 */
/*     instead of the current OFFDIG + OFFDIG2 < 1 */

    yesrel = TRUE_;
    offdig = 0.;
    tmp = sqrt((fabsl(d__[1])));
    if (tmp < rmin) {
	yesrel = FALSE_;
    }
    if (! yesrel) {
	goto L11;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	tmp2 = sqrt((d__1 = d__[i__], fabsl(d__1)));
	if (tmp2 < rmin) {
	    yesrel = FALSE_;
	}
	if (! yesrel) {
	    goto L11;
	}
	offdig2 = (d__1 = e[i__ - 1], fabsl(d__1)) / (tmp * tmp2);
	if (offdig + offdig2 >= .999) {
	    yesrel = FALSE_;
	}
	if (! yesrel) {
	    goto L11;
	}
	tmp = tmp2;
	offdig = offdig2;
/* L10: */
    }
L11:
    if (yesrel) {
	*info = 0;
	return 0;
    } else {
    }


/*     *** MORE TO BE IMPLEMENTED *** */


/*     Test if the lower bidiagonal matrix L from T = L D L^T */
/*     (zero shift facto) is well conditioned */


/*     Test if the upper bidiagonal matrix U from T = U D U^T */
/*     (zero shift facto) is well conditioned. */
/*     In this case, the matrix needs to be flipped and, at the end */
/*     of the eigenvector computation, the flip needs to be applied */
/*     to the computed eigenvectors (and the support) */


    return 0;

/*     END OF XDRRR */

} /* xdrrr_ */
