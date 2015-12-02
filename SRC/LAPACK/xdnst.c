/* dlanst.f -- translated by f2c (version 20061008) */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

/* Table of constant values */
static int c__1 = 1;

long double xdnst_(char *norm, int *n, long double *d__, long double *e)
{
    /* System generated locals */
    int i__1;
    long double ret_val, d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    // long double sqrt(long double);

    /* Local variables */
    int i__;
    long double sum, scale;
    extern int xlsame_(char *, char *);
    long double anorm;
    extern /* Subroutine */ int xdssq_(int *, long double *, int *, 
	    long double *, long double *);


/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  XDNST  returns the value of the one norm,  or the Frobenius norm, or */
/*  the  infinity norm,  or the  element of  largest absolute value  of a */
/*  real symmetric tridiagonal matrix A. */

/*  Description */
/*  =========== */

/*  XDNST returns the value */

/*     XDNST = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
/*              ( */
/*              ( norm1(A),         NORM = '1', 'O' or 'o' */
/*              ( */
/*              ( normI(A),         NORM = 'I' or 'i' */
/*              ( */
/*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e' */

/*  where  norm1  denotes the  one norm of a matrix (maximum column sum), */
/*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and */
/*  normF  denotes the  Frobenius norm of a matrix (square root of sum of */
/*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm. */

/*  Arguments */
/*  ========= */

/*  NORM    (input) CHARACTER*1 */
/*          Specifies the value to be returned in XDNST as described */
/*          above. */

/*  N       (input) INT */
/*          The order of the matrix A.  N >= 0.  When N = 0, XDNST is */
/*          set to zero. */

/*  D       (input) LONG DOUBLE PRECISION array, dimension (N) */
/*          The diagonal elements of A. */

/*  E       (input) LONG DOUBLE PRECISION array, dimension (N-1) */
/*          The (n-1) sub-diagonal or super-diagonal elements of A. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --e;
    --d__;

    /* Function Body */
    if (*n <= 0) {
	anorm = 0.;
    } else if (xlsame_(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	anorm = (d__1 = d__[*n], fabsl(d__1));
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__2 = anorm, d__3 = (d__1 = d__[i__], fabsl(d__1));
	    anorm = fmaxl(d__2,d__3);
/* Computing MAX */
	    d__2 = anorm, d__3 = (d__1 = e[i__], fabsl(d__1));
	    anorm = fmaxl(d__2,d__3);
/* L10: */
	}
    } else if (xlsame_(norm, "O") || *(unsigned char *)
	    norm == '1' || xlsame_(norm, "I")) {

/*        Find norm1(A). */

	if (*n == 1) {
	    anorm = fabsl(d__[1]);
	} else {
/* Computing MAX */
	    d__3 = fabsl(d__[1]) + fabsl(e[1]), d__4 = (d__1 = e[*n - 1], fabsl(
		    d__1)) + (d__2 = d__[*n], fabsl(d__2));
	    anorm = fmaxl(d__3,d__4);
	    i__1 = *n - 1;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing MAX */
		d__4 = anorm, d__5 = (d__1 = d__[i__], fabsl(d__1)) + (d__2 = e[
			i__], fabsl(d__2)) + (d__3 = e[i__ - 1], fabsl(d__3));
		anorm = fmaxl(d__4,d__5);
/* L20: */
	    }
	}
    } else if (xlsame_(norm, "F") || xlsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	if (*n > 1) {
	    i__1 = *n - 1;
	    xdssq_(&i__1, &e[1], &c__1, &scale, &sum);
	    sum *= 2;
	}
	xdssq_(n, &d__[1], &c__1, &scale, &sum);
	anorm = scale * sqrt(sum);
    }

    ret_val = anorm;
    return ret_val;

/*     End of XDNST */

} /* xdnst_ */
