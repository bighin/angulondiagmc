/* njgraf.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    char namsub[6];
} nam_;

#define nam_1 nam_

union {
    struct {
	integer j23[72]	/* was [24][3] */, arrow[72]	/* was [24][3] */, 
		line[120]	/* was [60][2] */, lcol[120]	/* was [60][2]
		 */;
	logical tabs[24];
	integer nbtr;
    } _1;
    struct {
	logical sum6j[20], t6j[20], jt[12], js[12];
	integer inver[60], jnsum[12], jinv[12], n6jn[20], in6j[20], jsumt[120]
			/* was [20][6] */;
    } _2;
} tree_;

#define tree_1 (tree_._1)
#define tree_2 (tree_._2)

struct {
    logical cut;
} cutdig_;

#define cutdig_1 cutdig_

union {
    struct {
	integer jdiag[144]	/* was [48][3] */, arr[144]	/* was [48][3]
		 */, tab1[120]	/* was [60][2] */, il[48], ih[48], npoint[24],
		 nbnode, ifirst, ilast, iparts, ipartl, npart, icross, nfree, 
		itfree[20], nfin, nc;
    } _1;
    struct {
	integer j12[576]	/* was [4][12][12] */;
    } _2;
    struct {
	integer jtem4[240]	/* was [12][20] */, jtem5[240]	/* was [12][
		20] */, jtem6[12], nsum6j[20], j6sum[20];
    } _3;
} graph_;

#define graph_1 (graph_._1)
#define graph_2 (graph_._2)
#define graph_3 (graph_._3)

struct {
    integer nograf;
    real tgraf;
    integer nogen;
    real tgen;
    integer norac;
    real trac;
} timeg_;

#define timeg_1 timeg_

struct {
    integer m, n, j1[60], j2[36]	/* was [12][3] */, j3[36]	/* 
	    was [12][3] */;
    logical free[60];
} couple_;

#define couple_1 couple_

union {
    struct {
	integer j6c, j7c, j8c, j9c, jwc, j6[180], j7[180], j8[180], j9[40], 
		kw[120]	/* was [6][20] */, jdel, ldel[40]	/* was [20][2]
		 */;
	logical sumvar[60];
	integer mp;
    } _1;
    struct {
	integer j6c, j7c, j8c, j9c, jwc, j6[180], j7[180], j8[180], j9[40], 
		jw[120]	/* was [6][20] */, jdel, ldel[40]	/* was [20][2]
		 */;
	logical sumvar[60];
	integer mp;
    } _2;
} argu_;

#define argu_1 (argu_._1)
#define argu_2 (argu_._2)

union {
    struct {
	integer i6c, i7c, i8c, i9c, idel, iwc;
    } _1;
    struct {
	integer eqv_pad[7];
    } _2;
} const_;

#define const_1 (const_._1)
#define const_2 (const_._2)

union {
    struct {
	integer ibug1, ibug2, ibug3, ibug4, ibug5, ibug6;
    } _1;
    struct {
	integer ibug1, ibug2, ibug3;
    } _2;
} debug_;

#define debug_1 (debug_._1)
#define debug_2 (debug_._2)

struct {
    real gam[100];
} facts_;

#define facts_1 facts_

struct {
    integer jkp[6]	/* was [2][3] */, jarr[6]	/* was [2][3] */, it2,
	     it3, it5;
} keep_;

#define keep_1 keep_

struct {
    integer j6cc, j7cc, j8cc, j9cc, jwcc, jdelc;
} dim_;

#define dim_1 dim_

struct {
    integer ial[48], if1, if2, node;
} build_;

#define build_1 build_

union {
    struct {
	integer i__, j, k, l, m, n;
    } _1;
    struct {
	integer ist[6];
    } _2;
} jrac_;

#define jrac_1 (jrac_._1)
#define jrac_2 (jrac_._2)

struct {
    integer ird, ipd, incp, iccp, incb, iccb, inct, icct;
} inform_;

#define inform_1 inform_

struct {
    integer j6p[40], j7p[40], j8p[40], j9p[40], jword[120]	/* was [6][20]
	     */, nlsum, nbj[10], nb6j[10], k6cp[10], k7cp[10], k8cp[10], k9cp[
	    10], jsum6[12], jsum4[240]	/* was [12][20] */, jsum5[240]	/* 
	    was [12][20] */, inv6j[20];
} sumarg_;

#define sumarg_1 sumarg_

struct {
    integer nzero, jzero[20];
} zer_;

#define zer_1 zer_

/* Table of constant values */

static integer c__8 = 8;
static integer c__0 = 0;
static integer c__48 = 48;
static integer c__12 = 12;
static integer c__1 = 1;
static integer c__24 = 24;
static integer c__2 = 2;
static integer c__10 = 10;
static integer c__3 = 3;
static integer c__4 = 4;

#define namsub (nam_1.namsub)
#define nbtr (tree_1.nbtr)
#define cut (cutdig_1.cut)
#define ih (graph_1.ih)
#define nbnode (graph_1.nbnode)
#define ifirst (graph_1.ifirst)
#define ilast (graph_1.ilast)
#define ipartl (graph_1.ipartl)
#define nfree (graph_1.nfree)
#define nfin (graph_1.nfin)
#define nc (graph_1.nc)
#define nograf (timeg_1.nograf)
#define tgraf (timeg_1.tgraf)
#define m (couple_1.m)
#define j6c (argu_1.j6c)
#define j7c (argu_1.j7c)
#define j8c (argu_1.j8c)
#define j9c (argu_1.j9c)
#define jwc (argu_1.jwc)
#define jdel (argu_1.jdel)
#define mp (argu_1.mp)
#define i6c (const_1.i6c)
#define i7c (const_1.i7c)
#define i8c (const_1.i8c)
#define i9c (const_1.i9c)
#define idel (const_1.idel)
#define iwc (const_1.iwc)
#define ibug3 (debug_1.ibug3)
/*     ******************************                                    ABBY0162 */
/* Subroutine */ int njgraf_(real *recup, logical *fail)
{
    /* Initialized data */

    static char name__[6] = "NJGRAF";

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer jf, jf1;
    real tt0, tt1;
    integer ncp;
    logical find;
    integer jpol, jump;
    extern /* Subroutine */ int zero_(integer *, integer *, logical *), 
	    cut1l_(logical *), cut2l_(logical *), order_(void), cutnl_(
	    logical *), bubble_(integer *, logical *), diagrm_(integer *), 
	    search_(logical *), second_(real *), settab_(logical *), triang_(
	    logical *), setdim_(void), sprate_(integer *), gensum_(real *), 
	    square_(void), lolpop_(logical *), printj_(char *, integer *, 
	    ftnlen), polygn_(integer *);

/*     ******************************                                    ABBY0164 */
/*                                                                       ABBY0165 */
/*                                                                       ABBY0166 */
/*  ***THIS IS THE MAIN PROGRAM.IT HANDLES  ALL THE ANALYSIS OF THE      ABBY0167 */
/*  ***RECOUPLING COEFFICIENT WITHOUT REFERRING EXPLICITLY TO THE VALUES ABBY0168 */
/*  ***OF ANGULAR MOMENTA WHICH ARE IN J1(J),EXCEPT FOR ZERO IN CASE FREEABBY0169 */
/*  ***=.FALSE. .LIKE NJSYM IT PREPARES ARRAYS OF ARGUMENTS FOR PHASE    ABBY0170 */
/*  ***FACTORS,(2*J+1) FACTORS AND 6J COEFFICIENTS TO BE COMPUTED IN     ABBY0171 */
/*  ***GENSUM,WHICH CAN ALSO BE CALLED SEPARATELY WHEN ONLY THE NUMERICALABBY0172 */
/*  ***VALUES OF ANGULAR MOMENTA CHANGE.THESE VARIABLE ANGULAR MOMENTA   ABBY0173 */
/*  ***SHOULD BE DECLARED FREE(J)=.TRUE.,SO THAT THE FORMULA PREPARED FORABBY0174 */
/*  ***GENSUM SHOULD BE CORRECT WHEN  J1 IS NOT ZERO.                    ABBY0175 */
/*  ***FAIL WILL BE TRUE WHEN THE RECOUPLING COEFFICIENT IS ZERO BECAUSE ABBY0176 */
/*  ***OF UNSATISFIED DELTA OR OTHER SIMILAR CAUSES.                     ABBY0177 */
/*                                                                       ABBY0178 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY0179 */
/*                                                                       ABBY0180 */
/*                                                                       ABBY0184 */
/*                                                                       ABBY0186 */
/*                                                                       ABBY0188 */
/*                                                                       ABBY0190 */
/*                                                                       ABBY0199 */
/*                                                                       ABBY0202 */
/*                                                                       ABBY0205 */
/*                                                                       ABBY0206 */
/*                                                                       ABBY0208 */
/*                                                                       ABBY0210 */
/*     COMMON BLOCK /FACTS / USED TO STORE LN(I!) IN THE NEW RACAH       ABBY0211 */
/*     ROUTINE WRITTEN BY STAN SCOTT.                                    ABBY0212 */
/*                                                                       ABBY0213 */
/*                                                                       ABBY0215 */
/*                                                                       ABBY0217 */
/*                                                                       ABBY0218 */
/* ***BUILDING UP OF THE UNSTRUCTURED GRAPH                              ABBY0219 */
/*                                                                       ABBY0220 */
    if (ibug3 == 1) {
	second_(&tt0);
    }
    ++nograf;
    *fail = FALSE_;
    j6c = 0;
    j7c = 0;
    j8c = 0;
    j9c = 0;
    jwc = 0;
    jdel = 0;
    i6c = 1;
    i7c = 1;
    i8c = 1;
    i9c = 1;
    iwc = 1;
    idel = 1;
    setdim_();
    nfin = 0;
    cut = FALSE_;
    settab_(fail);
    ++m;
    if (*fail) {
	goto L7;
    }
    --m;
    jf = 0;
    jf1 = 0;
/*                                                                       ABBY0246 */
/*   ***LOCATING AND HANDLING OF ZEROS                                   ABBY0247 */
/*                                                                       ABBY0248 */
    zero_(&jf1, &jf, fail);
    if (*fail) {
	goto L7;
    }
    mp = m;
    if (nbtr == 0) {
	goto L6;
    }
    jump = 1;
/*                                                                       ABBY0254 */
/*   ***BUILDING OF A FLAT DIAGRAM OUT OF THE UNSTRUCTURED GRAPH.        ABBY0255 */
/*   ***THERE MAY BE SEVERAL FLAT DIAGRAMS OUT OF THE ORIGINAL           ABBY0256 */
/*   ***GRAPH,IN CASE OF POSSIBLE CUTS.THEN THE FLAT DIAGRAMS            ABBY0257 */
/*   ***WILL HAVE FREE ENDS.                                             ABBY0258 */
/*                                                                       ABBY0259 */
L1:
    diagrm_(&jump);
/* Computing MAX */
    i__1 = 0, i__2 = nfree - 2;
    nfin = max(i__1,i__2);
/*                                                                       ABBY0262 */
    if (nfin != 0) {
	jump = 3;
/*                                                                       ABBY0265 */
/*  ***HANDLING OF FREE ENDS IF A CUT WAS FOUND                          ABBY0266 */
/*                                                                       ABBY0267 */
	cutnl_(fail);
	if (*fail) {
	    goto L7;
	}
    } else {
	jump = 2;
	if (nfree == 1) {
	    cut1l_(fail);
	    if (*fail) {
		goto L7;
	    }
	} else {
	    if (nfree > 1) {
		cut2l_(fail);
		if (*fail) {
		    goto L7;
		}
	    }
	}
    }
/*                                                                       ABBY0282 */
    nbtr += nfin;
    if (nbtr != 0) {
	cut = TRUE_;
    }
/*                                                                       ABBY0285 */
/*  ***ANALYSIS OF THE FLAT DIAGRAM.                                     ABBY0286 */
/*  ***CLOSED CIRCUITS OF INCREASING ORDER NC ARE SEARCHED,ANALYSED,AND  ABBY0287 */
/*  ***TAKEN OUT OF THE FLAT DIAGRAM,THUS REDUCING THE NUMBER OF NODES,  ABBY0288 */
/*  ***NBNODE.                                                           ABBY0289 */
/*                                                                       ABBY0291 */
    nc = 0;
L10:
    ++nc;
    search_(&find);
    if (! find) {
	goto L10;
    }
    ncp = nc - 2;
    jpol = 0;
    if (m == mp && nc > 3) {
	setdim_();
    }
    if (ipartl > 2) {
	polygn_(&jpol);
    }
    switch (nc) {
	case 1:  goto L11;
	case 2:  goto L12;
	case 3:  goto L13;
	case 4:  goto L14;
    }
L11:
    lolpop_(fail);
    if (*fail) {
	goto L7;
    }
    goto L15;
L12:
    bubble_(&jpol, fail);
    if (*fail) {
	goto L7;
    }
    goto L15;
L13:
    triang_(fail);
    if (*fail) {
	goto L7;
    }
    goto L15;
L14:
    square_();
L15:
    nbnode += -2;
    if (nbnode == 0) {
	goto L9;
    }
    ifirst = ih[0];
    ilast = ih[nbnode - 1];
/*                                                                       ABBY0315 */
/* ***PRINTJ IS AN ALL PURPOSE PRINTING SUBROUTINE CALLED FROM MANY PLACESABBY0316 */
/*                                                                       ABBY0317 */
    printj_(namsub, &c__8, (ftnlen)6);
    if (nbnode == nfin) {
	goto L9;
    }
    nc = ncp;
/*                                                                       ABBY0321 */
/*  ***PROCEED TO OTHER CIRCUITS OF ORDER NC-1                           ABBY0322 */
/*                                                                       ABBY0323 */
    goto L10;
L9:
    if (nbtr == 0) {
	goto L6;
    }
    if (jump == 3) {
	order_();
    }
/*                                                                       ABBY0327 */
/* ***AT THIS STAGE,THE FLAT DIAGRAM HAS BEEN REDUCED TO NODES           ABBY0328 */
/* ***INVOLVING FREE ENDS.PROCEED TO BUILD OTHER FLAT DIAGRAMS           ABBY0329 */
/* ***IF NECESSARY.                                                      ABBY0330 */
/*                                                                       ABBY0331 */
    goto L1;
/*                                                                       ABBY0333 */
/*  ***ALL PARTS OF THE ORIGINAL GRAPH HAVE BEEN REDUCED.                ABBY0334 */
/*                                                                       ABBY0335 */
L7:
    *recup = 0.f;
    --m;
    return 0;
L6:
    printj_(name__, &c__0, (ftnlen)6);
/*                                                                       ABBY0340 */
/*                                                                       ABBY0341 */
/*  ***PREPARATION OF THE RESULTS,AND SEPARATION IN SEVERAL SUMS         ABBY0342 */
/*  *** IF CUTS HAVE BEEN DETECTED,ALSO IN THE FLAT DIAGRAM ITSELF       ABBY0343 */
/*                                                                       ABBY0344 */
    second_(&tt1);
    tgraf = tgraf + tt1 - tt0;
    sprate_(&m);
    --m;
/*                                                                       ABBY0349 */
/*                                                                       ABBY0350 */
/*  ***GENSUM COMPUTES THE NUMERICAL VALUE OF THE RECOUPLING             ABBY0351 */
/*  ***COEFFICIENT.                                                      ABBY0352 */
/*                                                                       ABBY0353 */
    gensum_(recup);
/*                                                                       ABBY0355 */
    return 0;
} /* njgraf_ */

#undef namsub
#undef nbtr
#undef cut
#undef ih
#undef nbnode
#undef ifirst
#undef ilast
#undef ipartl
#undef nfree
#undef nfin
#undef nc
#undef nograf
#undef tgraf
#undef m
#undef j6c
#undef j7c
#undef j8c
#undef j9c
#undef jwc
#undef jdel
#undef mp
#undef i6c
#undef i7c
#undef i8c
#undef i9c
#undef idel
#undef iwc
#undef ibug3

#define namsub (nam_1.namsub)
#define jdiag (graph_1.jdiag)
#define arr (graph_1.arr)
#define tab1 (graph_1.tab1)
#define il (graph_1.il)
#define ih (graph_1.ih)
#define npoint (graph_1.npoint)
#define nbnode (graph_1.nbnode)
#define ifirst (graph_1.ifirst)
#define ilast (graph_1.ilast)
#define j6c (argu_1.j6c)
#define j8c (argu_1.j8c)
#define j9c (argu_1.j9c)
#define jwc (argu_1.jwc)
#define j6 (argu_1.j6)
#define j8 (argu_1.j8)
#define j9 (argu_1.j9)
#define kw (argu_1.kw)
#define mp (argu_1.mp)

/*                                                                       ABBY0358 */
/*                                                                       ABBY0359 */
/*     ****************************                                      ABBY0360 */
/* Subroutine */ int bubble_(integer *jpol, logical *fail)
{
    /* Initialized data */

    static char name__[6] = "BUBBLE";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    integer i__, k, l, i1, i2, k2, l1, k23, it, il1, it1, it2;
    extern /* Subroutine */ int delta_(integer *, integer *, logical *), 
	    phase_(integer *, integer *, integer *), phase2_(integer *);


#define kw_ref(a_1,a_2) kw[(a_2)*6 + a_1 - 7]
#define arr_ref(a_1,a_2) arr[(a_2)*48 + a_1 - 49]
#define tab1_ref(a_1,a_2) tab1[(a_2)*60 + a_1 - 61]
#define jdiag_ref(a_1,a_2) jdiag[(a_2)*48 + a_1 - 49]

/*     ****************************                                      ABBY0362 */
/*                                                                       ABBY0363 */
/*  ***REDUCES A CIRCUIT OF ORDER 2,GIVING DELTA FUNCTION AND PHASE      ABBY0364 */
/*  ***FACTORS.                                                          ABBY0365 */
/*                                                                       ABBY0366 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY0367 */
/*                                                                       ABBY0368 */
/*                                                                       ABBY0371 */
/*                                                                       ABBY0373 */
/*                                                                       ABBY0375 */
/*                                                                       ABBY0377 */
/*                                                                       ABBY0384 */
/*                                                                       ABBY0386 */
    s_copy(namsub, name__, (ftnlen)6, (ftnlen)6);
    k2 = 2;
    k23 = 3;
    i1 = 1;
    i2 = 1;
    it1 = npoint[0];
    it2 = npoint[1];
/*                                                                       ABBY0394 */
    if (it2 == ilast) {
	if (it1 != ifirst) {
	    it2 = it1;
	    it1 = ilast;
	}
	i1 = -1;
	k23 = 2;
	i2 = 2;
    }
/*                                                                       ABBY0404 */
    phase_(&it1, jdiag, &c__48);
    k = (i__1 = (arr_ref(it2, 1) * 3 + (arr_ref(it2, 2) << 1) + arr_ref(it2, 
	    3)) / 2, abs(i__1)) + 1;
    if (k != 4) {
	phase2_(&jdiag_ref(it2, k));
    }
    if (nbnode == 2) {
	return 0;
    }
    il1 = il[it2 - 1] + i1;
    it = ih[il1 - 1];
    arr_ref(it, k23) = arr_ref(it1, k23);
    l = jdiag_ref(it1, k23);
    l1 = jdiag_ref(it, k23);
    jdiag_ref(it, k23) = l;
/*                                                                       ABBY0415 */
/*                                                                       ABBY0416 */
    if (*jpol != 1) {
	delta_(&l, &l1, fail);
	if (*fail) {
	    return 0;
	}
    } else {
	--mp;
	kw_ref(2, jwc) = l;
	j6[j6c - 2] = l;
	j6[j6c - 1] = l;
	if (k == 2) {
	    j8[j8c - 1] = l;
	}
    }
/*                                                                       ABBY0427 */
    tab1_ref(l, i2) = it;
/*                                                                       ABBY0429 */
    if (it1 != ilast) {
	if (it2 == ilast) {
	    tab1_ref(l, 1) = ih[1];
	    il1 = 2;
	    k2 = 1;
	}
/*                                                                       ABBY0436 */
	i__1 = nbnode;
	for (i__ = il1; i__ <= i__1; ++i__) {
	    it = ih[i__ - 1];
	    il[it - 1] = i__ - k2;
	    ih[i__ - k2 - 1] = it;
/* L5: */
	}
/*                                                                       ABBY0442 */
    }
/*                                                                       ABBY0444 */
/* L6: */
    j9[j9c] = l;
    j9c += 2;
    j9[j9c - 1] = l;
/*                                                                       ABBY0448 */
    return 0;
} /* bubble_ */

#undef namsub
#undef jdiag
#undef arr
#undef tab1
#undef il
#undef ih
#undef npoint
#undef nbnode
#undef ifirst
#undef ilast
#undef j6c
#undef j8c
#undef j9c
#undef jwc
#undef j6
#undef j8
#undef j9
#undef kw
#undef mp


#undef jdiag_ref
#undef tab1_ref
#undef arr_ref
#undef kw_ref

#define jdiag (graph_1.jdiag)
#define arr (graph_1.arr)

/*                                                                       ABBY0451 */
/*                                                                       ABBY0452 */
/*     **********************                                            ABBY0453 */
/* Subroutine */ int change_(integer *l, integer *k)
{
    integer jp, jar;
    extern /* Subroutine */ int phase_(integer *, integer *, integer *);


#define arr_ref(a_1,a_2) arr[(a_2)*48 + a_1 - 49]
#define jdiag_ref(a_1,a_2) jdiag[(a_2)*48 + a_1 - 49]

/*     **********************                                            ABBY0455 */
/*                                                                       ABBY0456 */
/*     EXCHANGES THE FREE ENDS IN EITHER FIRST OR LAST TRIAD OF JDIAG.   ABBY0457 */
/*                                                                       ABBY0458 */
/*     IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY0459 */
/*                                                                       ABBY0460 */
/*                                                                       ABBY0463 */
/*                                                                       ABBY0465 */
/*                                                                       ABBY0469 */
    phase_(l, jdiag, &c__48);
    jp = jdiag_ref(*l, *k);
    jdiag_ref(*l, *k) = jdiag_ref(*l, 1);
    jdiag_ref(*l, 1) = jp;
    jar = arr_ref(*l, *k);
    arr_ref(*l, *k) = arr_ref(*l, 1);
    arr_ref(*l, 1) = jar;
/*                                                                       ABBY0477 */
    return 0;
} /* change_ */

#undef jdiag
#undef arr


#undef jdiag_ref
#undef arr_ref


/*                                                                       ABBY0480 */
/*                                                                       ABBY0481 */
/*     *****************************************                         ABBY0482 */
/* Subroutine */ int chvar_(integer *jp, integer *nbc, integer *kbc, logical *
	jt, integer *jinv, integer *nsum)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, kb, jk;

/*     *****************************************                         ABBY0484 */
/*                                                                       ABBY0485 */
/*   ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                               ABBY0486 */
/*                                                                       ABBY0487 */
/*   ***CHANGE THE ORDER OF SUMMATION VARIABLE TO BE ABLE TO PERFORM     ABBY0488 */
/*   ***SEPARATELY THE SUMMATIONS IN GENSUM.                             ABBY0489 */
/*                                                                       ABBY0490 */
/*                                                                       ABBY0493 */
/*                                                                       ABBY0494 */
    /* Parameter adjustments */
    --jp;
    --jinv;
    --jt;

    /* Function Body */
    kb = *kbc + 1;
/*                                                                       ABBY0496 */
    if (kb <= *nbc) {
	i__1 = *nbc;
	for (i__ = kb; i__ <= i__1; ++i__) {
	    jk = jp[i__];
	    if (jt[jk]) {
		++(*kbc);
		jp[i__] = jp[*kbc];
		jp[*kbc] = jinv[jk];
	    }
/* L1: */
	}
    }
/*                                                                       ABBY0507 */
    return 0;
} /* chvar_ */
#define jdiag (graph_1.jdiag)
#define arr (graph_1.arr)
#define il (graph_1.il)
#define ih (graph_1.ih)
#define nbnode (graph_1.nbnode)
#define itfree (graph_1.itfree)
#define m (couple_1.m)
#define j9c (argu_1.j9c)
#define j9 (argu_1.j9)

/*                                                                       ABBY0510 */
/*                                                                       ABBY0511 */
/*     **********************                                            ABBY0512 */
/* Subroutine */ int cut1l_(logical *fail)
{
    /* Initialized data */

    static char name__[6] = "CUT1L ";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, j, j0, it, il1, ilp;
    extern /* Subroutine */ int zero_(integer *, integer *, logical *), 
	    delta_(integer *, integer *, logical *), printj_(char *, integer *
	    , ftnlen);


#define arr_ref(a_1,a_2) arr[(a_2)*48 + a_1 - 49]
#define jdiag_ref(a_1,a_2) jdiag[(a_2)*48 + a_1 - 49]

/*     **********************                                            ABBY0514 */
/*                                                                       ABBY0515 */
/*  ***CUT ON ONE LINE,THAT WAS LEFT AS A FREE END IN JDIAG.PUTS         ABBY0516 */
/*  ***CORRESPONDING DELTA IN J23.                                       ABBY0517 */
/*                                                                       ABBY0518 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY0519 */
/*                                                                       ABBY0520 */
/*                                                                       ABBY0523 */
/*                                                                       ABBY0525 */
/*                                                                       ABBY0527 */
/*                                                                       ABBY0529 */
/*                                                                       ABBY0533 */
/*                                                                       ABBY0536 */
/*                                                                       ABBY0538 */
/*                                                                       ABBY0540 */
/*                                                                       ABBY0541 */
    it = itfree[0];
    j0 = jdiag_ref(it, 1);
    delta_(&j0, &m, fail);
    if (*fail) {
	goto L2;
    }
    delta_(&jdiag_ref(it, 3), &jdiag_ref(it, 2), fail);
    if (*fail) {
	goto L2;
    }
    jdiag_ref(it + 1, 3) = jdiag_ref(it, 3);
/*                                                                       ABBY0549 */
    if (arr_ref(it, 2) == arr_ref(it, 3)) {
	arr_ref(it + 1, 3) = 1;
	arr_ref(it - 1, 2) = -1;
    } else {
	if (arr_ref(it, 2) < arr_ref(it, 3)) {
	    arr_ref(it + 1, 3) = -1;
	    arr_ref(it - 1, 2) = 1;
	}
    }
/*                                                                       ABBY0559 */
    ++j9c;
    j9[j9c - 1] = jdiag_ref(it, 3);
    j = 2;
    zero_(&j, &j0, fail);
    if (*fail) {
	goto L2;
    }
    il1 = il[it];
/*                                                                       ABBY0566 */
    i__1 = nbnode;
    for (i__ = il1; i__ <= i__1; ++i__) {
	it = ih[i__ - 1];
	ilp = i__ - 1;
	il[it - 1] = ilp;
	ih[ilp - 1] = it;
/* L1: */
    }
/*                                                                       ABBY0573 */
    --nbnode;
/*                                                                       ABBY0575 */
L2:
    printj_(name__, &c__12, (ftnlen)6);
    return 0;
} /* cut1l_ */

#undef jdiag
#undef arr
#undef il
#undef ih
#undef nbnode
#undef itfree
#undef m
#undef j9c
#undef j9


#undef jdiag_ref
#undef arr_ref

#define j23 (tree_1.j23)
#define arrow (tree_1.arrow)
#define line (tree_1.line)
#define lcol (tree_1.lcol)
#define jdiag (graph_1.jdiag)
#define arr (graph_1.arr)
#define tab1 (graph_1.tab1)
#define itfree (graph_1.itfree)
#define j9c (argu_1.j9c)
#define j9 (argu_1.j9)

/*                                                                       ABBY0579 */
/*                                                                       ABBY0580 */
/*     **********************                                            ABBY0581 */
/* Subroutine */ int cut2l_(logical *fail)
{
    /* Initialized data */

    static char name__[6] = "CUT2L ";

    integer k1, l1, l2, k2, lc1, lc2, it1, it2, jt1, jt2;
    extern /* Subroutine */ int delta_(integer *, integer *, logical *), 
	    phase2_(integer *), otherj_(integer *, integer *, integer *, 
	    integer *, integer *), printj_(char *, integer *, ftnlen);


#define j23_ref(a_1,a_2) j23[(a_2)*24 + a_1 - 25]
#define arr_ref(a_1,a_2) arr[(a_2)*48 + a_1 - 49]
#define tab1_ref(a_1,a_2) tab1[(a_2)*60 + a_1 - 61]
#define line_ref(a_1,a_2) line[(a_2)*60 + a_1 - 61]
#define lcol_ref(a_1,a_2) lcol[(a_2)*60 + a_1 - 61]
#define jdiag_ref(a_1,a_2) jdiag[(a_2)*48 + a_1 - 49]
#define arrow_ref(a_1,a_2) arrow[(a_2)*24 + a_1 - 25]

/*     **********************                                            ABBY0583 */
/*                                                                       ABBY0584 */
/*  ***CUT ON TWO LINES THAT WERE LEFT AS FREE ENDS IN JDIAG.PUTS        ABBY0585 */
/*  ***CORRESPONDING DELTA IN J23.                                       ABBY0586 */
/*                                                                       ABBY0587 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY0588 */
/*                                                                       ABBY0589 */
/*                                                                       ABBY0592 */
/*                                                                       ABBY0594 */
/*                                                                       ABBY0596 */
/*                                                                       ABBY0598 */
/*                                                                       ABBY0599 */
/*                                                                       ABBY0603 */
/*                                                                       ABBY0606 */
/*                                                                       ABBY0609 */
/*                                                                       ABBY0611 */
/*                                                                       ABBY0612 */
    it1 = itfree[0];
    it2 = itfree[1];
    jt1 = jdiag_ref(it1, 1);
    jt2 = jdiag_ref(it2, 1);
    delta_(&jt1, &jt2, fail);
    if (*fail) {
	goto L1;
    }
    if (arr_ref(it1, 1) == arr_ref(it2, 1)) {
	phase2_(&jt1);
    }
    arr_ref(it2, 1) = -arr_ref(it1, 1);
    jdiag_ref(it2, 1) = jt1;
    tab1_ref(jt1, 2) = it2;
    j9[j9c] = jt1;
    j9c += 2;
    j9[j9c - 1] = jt1;
    otherj_(&c__0, &jt1, &l1, &lc1, &k1);
    otherj_(&c__0, &jt2, &l2, &lc2, &k2);
    j23_ref(l2, lc2) = jt1;
    line_ref(jt1, k1) = l2;
    lcol_ref(jt1, k1) = lc2;
    arrow_ref(l2, lc2) = -arrow_ref(l1, lc1);
/*                                                                       ABBY0632 */
L1:
    printj_(name__, &c__12, (ftnlen)6);
/*                                                                       ABBY0634 */
    return 0;
} /* cut2l_ */

#undef j23
#undef arrow
#undef line
#undef lcol
#undef jdiag
#undef arr
#undef tab1
#undef itfree
#undef j9c
#undef j9


#undef arrow_ref
#undef jdiag_ref
#undef lcol_ref
#undef line_ref
#undef tab1_ref
#undef arr_ref
#undef j23_ref

#define nbtr (tree_1.nbtr)
#define jdiag (graph_1.jdiag)
#define arr (graph_1.arr)
#define il (graph_1.il)
#define ih (graph_1.ih)
#define nbnode (graph_1.nbnode)
#define nfree (graph_1.nfree)
#define itfree (graph_1.itfree)
#define nfin (graph_1.nfin)
#define j9c (argu_1.j9c)
#define j9 (argu_1.j9)
#define jkp (keep_1.jkp)
#define jarr (keep_1.jarr)
#define it2 (keep_1.it2)
#define it3 (keep_1.it3)
#define it5 (keep_1.it5)

/*                                                                       ABBY0637 */
/*                                                                       ABBY0638 */
/*     **********************                                            ABBY0639 */
/* Subroutine */ int cutnl_(logical *fail)
{
    /* Initialized data */

    static char name__[6] = "CUTNL ";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, j, k, jt, it, il1, it1, it4, nfr, ntf, ilp;
    extern /* Subroutine */ int delta_(integer *, integer *, logical *), 
	    phase_(integer *, integer *, integer *), phase2_(integer *), 
	    printj_(char *, integer *, ftnlen);


#define arr_ref(a_1,a_2) arr[(a_2)*48 + a_1 - 49]
#define jkp_ref(a_1,a_2) jkp[(a_2)*2 + a_1 - 3]
#define jarr_ref(a_1,a_2) jarr[(a_2)*2 + a_1 - 3]
#define jdiag_ref(a_1,a_2) jdiag[(a_2)*48 + a_1 - 49]

/*     **********************                                            ABBY0641 */
/*                                                                       ABBY0642 */
/*  ***THIS SUBROUTINE  EXAMINES THE CASE WHERE THERE ARE MORE THAN      ABBY0643 */
/*  ***TWO FREE ENDS,BUT THEY ARE CONTIGUOUS,SO THAT THE GRAPH CAN       ABBY0644 */
/*  ***BE CUT WITHOUT DESTROYING THE FLAT STRUCTURE.                     ABBY0645 */
/*                                                                       ABBY0646 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY0647 */
/*                                                                       ABBY0648 */
/*                                                                       ABBY0651 */
/*                                                                       ABBY0653 */
/*                                                                       ABBY0655 */
/*                                                                       ABBY0657 */
/*                                                                       ABBY0660 */
/*                                                                       ABBY0664 */
/*                                                                       ABBY0666 */
/*                                                                       ABBY0669 */
/*                                                                       ABBY0671 */
/*                                                                       ABBY0672 */
    ntf = itfree[(2144 + (0 + ((nfree - 1) << 2)) - 2144) / 4] - itfree[0];
    if (ntf > nfree) {
	goto L8;
    }
    it2 = itfree[0];
    it3 = itfree[nfree - 1];
    it1 = it2 - 1;
    it4 = it3 + 1;
/*                                                                       ABBY0679 */
    if (ntf != nfree) {
	jt = jdiag_ref(it2, 3);
	delta_(&jt, &jdiag_ref(it3, 2), fail);
/*                                                                       ABBY0683 */
	if (*fail) {
	    goto L8;
	}
/*                                                                       ABBY0685 */
	if (arr_ref(it2, 3) == arr_ref(it3, 2)) {
	    phase2_(&jt);
	    arr_ref(it2, 3) = -arr_ref(it2, 3);
	    arr_ref(it1, 2) = -arr_ref(it1, 2);
	}
/*                                                                       ABBY0691 */
	jdiag_ref(it3, 2) = jt;
	jdiag_ref(it4, 3) = jt;
	j9[j9c] = jt;
	j9c += 2;
	j9[j9c - 1] = jt;
	nbtr += nfree;
	it5 = 0;
	goto L6;
    }
    nfr = 0;
/*                                                                       ABBY0702 */
    i__1 = it3;
    for (it5 = it2; it5 <= i__1; ++it5) {
	++nfr;
	if (itfree[nfr - 1] > it5) {
	    goto L4;
	}
/* L3: */
    }
/*                                                                       ABBY0707 */
L4:
    jkp_ref(1, 1) = jdiag_ref(it5, 1);
    jarr_ref(1, 1) = -arr_ref(it5, 1);
    jkp_ref(1, 2) = jdiag_ref(it2, 3);
    jarr_ref(1, 2) = -arr_ref(it2, 3);
    jkp_ref(1, 3) = jdiag_ref(it3, 2);
    jarr_ref(1, 3) = -arr_ref(it3, 2);
/*                                                                       ABBY0714 */
    for (j = 1; j <= 3; ++j) {
	jkp_ref(2, j) = jdiag_ref(it5, j);
	jarr_ref(2, j) = arr_ref(it5, j);
/* L5: */
    }
/*                                                                       ABBY0719 */
    jdiag_ref(it5, 2) = jdiag_ref(it3, 2);
    arr_ref(it5, 2) = arr_ref(it3, 2);
    jdiag_ref(it5, 3) = jdiag_ref(it2, 3);
    arr_ref(it5, 3) = arr_ref(it2, 3);
    ilp = il[it2 - 1];
    il[it5 - 1] = ilp;
    ih[ilp - 1] = it5;
    nbtr = nbtr + nfree + 2;
    phase_(&it5, jdiag, &c__48);
    k = (i__1 = (arr_ref(it5, 1) * 3 + (arr_ref(it5, 2) << 1) + arr_ref(it5, 
	    3)) / 2 + 1, abs(i__1));
    if (k != 4) {
	phase2_(&jdiag_ref(it5, k));
    }
L6:
    il1 = il[it4 - 1];
/*                                                                       ABBY0732 */
    i__1 = nbnode;
    for (i__ = il1; i__ <= i__1; ++i__) {
	it = ih[i__ - 1];
	ilp = i__ - nfree;
	il[it - 1] = ilp;
	ih[ilp - 1] = it;
/* L7: */
    }
/*                                                                       ABBY0739 */
    nbnode -= nfree;
    nfin = 0;
/*                                                                       ABBY0742 */
L8:
    printj_(name__, &c__8, (ftnlen)6);
/*                                                                       ABBY0744 */
    return 0;
} /* cutnl_ */

#undef nbtr
#undef jdiag
#undef arr
#undef il
#undef ih
#undef nbnode
#undef nfree
#undef itfree
#undef nfin
#undef j9c
#undef j9
#undef jkp
#undef jarr
#undef it2
#undef it3
#undef it5


#undef jdiag_ref
#undef jarr_ref
#undef jkp_ref
#undef arr_ref

#define cut (cutdig_1.cut)
#define j1 (couple_1.j1)
#define free (couple_1.free)
#define j6c (argu_1.j6c)
#define j7c (argu_1.j7c)
#define j8c (argu_1.j8c)
#define j9c (argu_1.j9c)
#define jwc (argu_1.jwc)
#define j6 (argu_1.j6)
#define j7 (argu_1.j7)
#define j8 (argu_1.j8)
#define j9 (argu_1.j9)
#define kw (argu_1.kw)
#define jdel (argu_1.jdel)
#define ldel (argu_1.ldel)
#define sumvar (argu_1.sumvar)
#define ibug3 (debug_1.ibug3)
#define j6cc (dim_1.j6cc)
#define j7cc (dim_1.j7cc)
#define j8cc (dim_1.j8cc)
#define j9cc (dim_1.j9cc)
#define jwcc (dim_1.jwcc)
#define jdelc (dim_1.jdelc)

/*                                                                       ABBY0747 */
/*                                                                       ABBY0748 */
/*     ****************************                                      ABBY0749 */
/* Subroutine */ int delta_(integer *ja, integer *jb, logical *fail)
{
    /* Format strings */
    static char fmt_1000[] = "(/2x,\002FROM DELTA\002,2x,\002JA=\002,i2,l2,5"
	    "x,\002JB=\002,i2,l2)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    integer i__, j, j61, j71, j81, j91, jw1, jdel1;

    /* Fortran I/O blocks */
    static cilist io___58 = { 0, 6, 0, fmt_1000, 0 };



#define kw_ref(a_1,a_2) kw[(a_2)*6 + a_1 - 7]
#define ldel_ref(a_1,a_2) ldel[(a_2)*20 + a_1 - 21]

/*     ****************************                                      ABBY0751 */
/*                                                                       ABBY0752 */
/*  ***TEST FOR DELTA(JA,JB).IF THEY ARE SUMMATION VARIABLES,THE SECOND  ABBY0753 */
/*  ***IS CHANGED INTO THE FIRST EVERYWHERE.IF THEY ARE FIXED,THEIR      ABBY0754 */
/*  ***VALUE IS CHECKED,AND FAIL PUT TO .TRUE. IF THEY DIFFER.           ABBY0755 */
/*                                                                       ABBY0756 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY0757 */
/*                                                                       ABBY0758 */
/*                                                                       ABBY0761 */
/*                                                                       ABBY0763 */
/*                                                                       ABBY0766 */
/*                                                                       ABBY0771 */
/*                                                                       ABBY0773 */
/*                                                                       ABBY0774 */
/*                                                                       ABBY0775 */
/* L1000: */
    if (ibug3 == 1) {
	s_wsfe(&io___58);
	do_fio(&c__1, (char *)&(*ja), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&sumvar[*ja - 1], (ftnlen)sizeof(logical));
	do_fio(&c__1, (char *)&(*jb), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&sumvar[*jb - 1], (ftnlen)sizeof(logical));
	e_wsfe();
    }
    if (sumvar[*ja - 1] && sumvar[*jb - 1]) {
	goto L2;
    }
    if (free[*ja - 1] || free[*jb - 1]) {
	++jdel;
	ldel_ref(jdel, 1) = *ja;
	ldel_ref(jdel, 2) = *jb;
	sumvar[*ja - 1] = FALSE_;
	sumvar[*jb - 1] = FALSE_;
	return 0;
    }
/*                                                                       ABBY0786 */
    if (j1[*ja - 1] != j1[*jb - 1]) {
	*fail = TRUE_;
    }
    cut = TRUE_;
    return 0;
/*                                                                       ABBY0790 */
L2:
    if (j6c != j6cc) {
	j61 = j6cc + 1;
/*                                                                       ABBY0793 */
	i__1 = j6c;
	for (i__ = j61; i__ <= i__1; ++i__) {
	    if (j6[i__ - 1] == *jb) {
		j6[i__ - 1] = *ja;
	    }
/* L3: */
	}
/*                                                                       ABBY0797 */
    }
/*                                                                       ABBY0799 */
    if (j7c != j7cc) {
	j71 = j7cc + 1;
/*                                                                       ABBY0802 */
	i__1 = j7c;
	for (i__ = j71; i__ <= i__1; ++i__) {
	    if (j7[i__ - 1] == *jb) {
		j7[i__ - 1] = *ja;
	    }
/* L5: */
	}
    }
/*                                                                       ABBY0807 */
    if (j8c != j8cc) {
	j81 = j8cc + 1;
/*                                                                       ABBY0810 */
	i__1 = j8c;
	for (i__ = j81; i__ <= i__1; ++i__) {
	    if (j8[i__ - 1] == *jb) {
		j8[i__ - 1] = *ja;
	    }
/* L7: */
	}
    }
/*                                                                       ABBY0815 */
    if (j9c != j9cc) {
	j91 = j9cc + 1;
/*                                                                       ABBY0818 */
	i__1 = j9c;
	for (i__ = j91; i__ <= i__1; ++i__) {
	    if (j9[i__ - 1] == *jb) {
		j9[i__ - 1] = *ja;
	    }
/* L9: */
	}
    }
/*                                                                       ABBY0823 */
    if (jwc != jwcc) {
	jw1 = jwcc + 1;
/*                                                                       ABBY0826 */
	i__1 = jwc;
	for (i__ = jw1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 6; ++j) {
		if (kw_ref(j, i__) == *jb) {
		    kw_ref(j, i__) = *ja;
		}
/* L13: */
	    }
/* L14: */
	}
    }
/*                                                                       ABBY0833 */
    if (jdel != jdelc) {
	jdel1 = jdelc + 1;
/*                                                                       ABBY0836 */
	i__1 = jdel;
	for (i__ = jdel1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 2; ++j) {
		if (ldel_ref(i__, j) == *jb) {
		    ldel_ref(i__, j) = *ja;
		}
/* L16: */
	    }
/* L17: */
	}
/*                                                                       ABBY0842 */
	sumvar[*jb - 1] = FALSE_;
    }
/*                                                                       ABBY0845 */
    return 0;
} /* delta_ */

#undef cut
#undef j1
#undef free
#undef j6c
#undef j7c
#undef j8c
#undef j9c
#undef jwc
#undef j6
#undef j7
#undef j8
#undef j9
#undef kw
#undef jdel
#undef ldel
#undef sumvar
#undef ibug3
#undef j6cc
#undef j7cc
#undef j8cc
#undef j9cc
#undef jwcc
#undef jdelc


#undef ldel_ref
#undef kw_ref

#define j23 (tree_1.j23)
#define arrow (tree_1.arrow)
#define tabs (tree_1.tabs)
#define nbtr (tree_1.nbtr)
#define jdiag (graph_1.jdiag)
#define arr (graph_1.arr)
#define nbnode (graph_1.nbnode)
#define ifirst (graph_1.ifirst)
#define ilast (graph_1.ilast)
#define mp (argu_1.mp)
#define ial (build_1.ial)
#define if1 (build_1.if1)
#define if2 (build_1.if2)
#define node (build_1.node)

/*                                                                       ABBY0848 */
/*                                                                       ABBY0849 */
/*     ***********************                                           ABBY0850 */
/* Subroutine */ int diagrm_(integer *jump)
{
    /* Initialized data */

    static char name__[6] = "DIAGRM";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, j, l, i1, k1, k2, k3, l1, l2, jb, nb, lc, nd, kp, lp, jt, 
	    ich, jar, ict, nbp;
    extern /* Subroutine */ int way_(integer *, integer *, integer *, integer 
	    *, integer *), intab_(void), phase_(integer *, integer *, integer 
	    *);
    integer ntime;
    extern /* Subroutine */ int neibor_(integer *, integer *, integer *), 
	    otherj_(integer *, integer *, integer *, integer *, integer *), 
	    printj_(char *, integer *, ftnlen);


#define j23_ref(a_1,a_2) j23[(a_2)*24 + a_1 - 25]
#define arr_ref(a_1,a_2) arr[(a_2)*48 + a_1 - 49]
#define jdiag_ref(a_1,a_2) jdiag[(a_2)*48 + a_1 - 49]
#define arrow_ref(a_1,a_2) arrow[(a_2)*24 + a_1 - 25]

/*     ***********************                                           ABBY0852 */
/*                                                                       ABBY0853 */
/*     THIS SUBROUTINE BUILDS UP A FLAT DIAGRAM FROM THE TRIADS J23 AND  ABBY0854 */
/*  ***PLACES THEM IN JDIAG.ARROWS ARE IN ARR (INTEGER).THE DIAGRAM IS   ABBY0855 */
/*  ***BUILT SO AS TO MAXIMIZE THE NUMBER OF TRIADS INVOLVED,WITHIN  A   ABBY0856 */
/*  ***ONE-STEP-FORWARD-CHECK PROCESS.IF THE DIAGRAM DOES NOT            ABBY0857 */
/*  ***INCLUDE ALL THE NBTR TRIADS,IT WILL HAVE 'FREE ENDS'.JDIAG HAS    ABBY0858 */
/*  ***DIMENSION DOUBLE THAT OF J23,BECAUSE THE PATH MAY PROCEED EITHER  ABBY0859 */
/*  ***WAY.                                                              ABBY0860 */
/*                                                                       ABBY0861 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY0862 */
/*                                                                       ABBY0863 */
/*                                                                       ABBY0866 */
/*                                                                       ABBY0868 */
/*                                                                       ABBY0870 */
/*                                                                       ABBY0872 */
/*                                                                       ABBY0875 */
/*                                                                       ABBY0878 */
/*                                                                       ABBY0882 */
/*                                                                       ABBY0885 */
/*                                                                       ABBY0887 */
/*                                                                       ABBY0888 */
/*                                                                       ABBY0889 */
/*  ***INITIALIZATION                                                    ABBY0890 */
/*                                                                       ABBY0891 */
    if (*jump > 2) {
	goto L17;
    }
    if (*jump < 2) {
	nb = 0;
    }
L1:
    ++nb;
    if (tabs[nb - 1]) {
	goto L1;
    }
    node = nbtr;
    ilast = nbtr;
/*                                                                       ABBY0898 */
    for (j = 1; j <= 3; ++j) {
	jdiag_ref(node, j) = j23_ref(nb, j);
	arr_ref(node, j) = arrow_ref(nb, j);
/* L2: */
    }
/*                                                                       ABBY0903 */
    tabs[nb - 1] = TRUE_;
/*                                                                       ABBY0905 */
    i__1 = mp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ial[i__ - 1] = 0;
/* L15: */
    }
/*                                                                       ABBY0909 */
    if1 = jdiag_ref(node, 1);
    if2 = jdiag_ref(node, 3);
    ial[if1 - 1] = 1;
    ial[if2 - 1] = 1;
L17:
    ntime = 0;
    i1 = 1;
    k1 = 1;
    k2 = 2;
    k3 = 3;
L3:
    jb = jdiag_ref(node, k2);
    otherj_(&c__0, &jb, &l, &lc, &kp);
    neibor_(&lc, &l1, &l2);
    way_(&l, &l1, &l2, &ich, &nd);
    node += i1;
    tabs[l - 1] = TRUE_;
    jdiag_ref(node, k3) = j23_ref(l, lc);
    arr_ref(node, k3) = arrow_ref(l, lc);
    ict = ich * i1;
/*                                                                       ABBY0928 */
    if (ich <= 0) {
	lp = l1;
	l1 = l2;
	l2 = lp;
    }
/*                                                                       ABBY0934 */
    if (ict <= 0) {
	phase_(&l, j23, &c__24);
    }
    jdiag_ref(node, k1) = j23_ref(l, l1);
    arr_ref(node, k1) = arrow_ref(l, l1);
    jdiag_ref(node, k2) = j23_ref(l, l2);
    arr_ref(node, k2) = arrow_ref(l, l2);
    j = j23_ref(l, l1);
    ++ial[j - 1];
    j = j23_ref(l, l2);
    ++ial[j - 1];
    if (nd < 1) {
	goto L3;
    }
    ++ntime;
    ilast = max(node,ilast);
    ifirst = min(node,nbtr);
    nbp = ial[if1 - 1] + ial[if2 - 1];
    if (nbp > 3 || ntime > 1) {
	nbnode = ilast - ifirst + 1;
	nbtr -= nbnode;
/*                                                                       ABBY0952 */
/*  ***DEFINITION OF FREE ENDS AND OTHER QUANTITIES.                     ABBY0953 */
/*                                                                       ABBY0954 */
	intab_();
	printj_(name__, &c__12, (ftnlen)6);
	goto L50;
    }
/*                                                                       ABBY0959 */
    if (nbp > 2) {
	if (ial[if1 - 1] <= ial[if2 - 1]) {
	    jt = jdiag_ref(nbtr, 1);
	    jar = arr_ref(nbtr, 1);
	    jdiag_ref(nbtr, 1) = jdiag_ref(nbtr, 3);
	    arr_ref(nbtr, 1) = arr_ref(nbtr, 3);
	    jdiag_ref(nbtr, 3) = jt;
	    arr_ref(nbtr, 3) = jar;
	    phase_(&nbtr, jdiag, &c__48);
	}
    }
/*                                                                       ABBY0971 */
    node = nbtr;
    i1 = -1;
    k2 = 3;
    k3 = 2;
    goto L3;
/*                                                                       ABBY0977 */
L50:
    return 0;
} /* diagrm_ */

#undef j23
#undef arrow
#undef tabs
#undef nbtr
#undef jdiag
#undef arr
#undef nbnode
#undef ifirst
#undef ilast
#undef mp
#undef ial
#undef if1
#undef if2
#undef node


#undef arrow_ref
#undef jdiag_ref
#undef arr_ref
#undef j23_ref

#define norac (timeg_1.norac)
#define trac (timeg_1.trac)
#define gam (facts_1.gam)
#define i__ (jrac_1.i__)
#define j (jrac_1.j)
#define k (jrac_1.k)
#define l (jrac_1.l)
#define m (jrac_1.m)
#define n (jrac_1.n)

/*                                                                       ABBY0980 */
/*                                                                       ABBY0981 */
/*     **********************                                            ABBY0982 */
/* Subroutine */ int dracah_(real *rac)
{
    /* Initialized data */

    static real zero = 0.f;
    static real one = 1.f;
    static real two = 2.f;

    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    double exp(doublereal), pow_ri(real *, integer *);

    /* Local variables */
    integer j1, j2, j3, j4, j5, j6, j7, ki, kk;
    real tt0, tt1;
    integer ntr, numin, numax;
    extern /* Subroutine */ int second_(real *);
    real deltat;
    integer icount;

/*     **********************                                            ABBY0984 */
/*                                                                       ABBY0985 */
/*     IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY0986 */
/*                                                                       ABBY0987 */
/*  ***SUBROUTINE TO CALCULATE RACAH COEFFICIENTS                        ABBY0988 */
/*  ***THE ARGUMENTS I,J,K,L,M,N SHOULD BE TWICE THEIR ACTUAL VALUE      ABBY0989 */
/*  ***WORKS FOR INTEGER AND HALF-INTEGER VALUES OF ANGULAR MOMENTA.     ABBY0990 */
/*  ***THE ROUTINE MAKES USE OF THE GAM ARRAY, THUS SUBROUTINE FACTT     ABBY0991 */
/*  ***MUST BE CALLED BEFORE THIS ROUTINE IS USED.                       ABBY0992 */
/*  ***WRITTEN BY N S SCOTT.                                             ABBY0993 */
/*                                                                       ABBY0994 */
/*                                                                       ABBY0996 */
/*                                                                       ABBY1001 */
/*                                                                       ABBY1003 */
/*                                                                       ABBY1004 */
/* IN THE REAL VERSION PUT HERE IF(IBUG3.EQ.1)                           ABBY1005 */
    second_(&tt0);
    ++norac;
/*                                                                       ABBY1008 */
/*                                                                       ABBY1009 */
    j1 = i__ + j + m;
    j2 = k + l + m;
    j3 = i__ + k + n;
    j4 = j + l + n;
/* Computing MAX */
    i__1 = max(i__,j);
    if ((max(i__1,m) << 1) - j1 > 0 || j1 % 2 != 0) {
	goto L2;
    }
/* Computing MAX */
    i__1 = max(k,l);
    if ((max(i__1,m) << 1) - j2 > 0 || j2 % 2 != 0) {
	goto L2;
    }
/* Computing MAX */
    i__1 = max(i__,k);
    if ((max(i__1,n) << 1) - j3 > 0 || j3 % 2 != 0) {
	goto L2;
    }
/* Computing MAX */
    i__1 = max(j,l);
    if ((max(i__1,n) << 1) - j4 > 0 || j4 % 2 != 0) {
	goto L2;
    }
    goto L1;
L2:
    *rac = zero;
    return 0;
/*                                                                       ABBY1021 */
L1:
    j1 /= 2;
    j2 /= 2;
    j3 /= 2;
    j4 /= 2;
    j5 = (i__ + j + k + l) / 2;
    j6 = (i__ + l + m + n) / 2;
    j7 = (j + k + m + n) / 2;
/* Computing MAX */
    i__1 = max(j1,j2), i__1 = max(i__1,j3);
    numin = max(i__1,j4) + 1;
/* Computing MIN */
    i__1 = min(j5,j6);
    numax = min(i__1,j7) + 1;
    *rac = one;
    icount = 0;
/*                                                                       ABBY1034 */
    if (numin == numax) {
	goto L4;
    }
    ++numin;
/*                                                                       ABBY1037 */
    i__1 = numax;
    for (kk = numin; kk <= i__1; ++kk) {
	ki = numax - icount;
	*rac = one - *rac * (real) (ki * (j5 - ki + 2) * (j6 - ki + 2) * (j7 
		- ki + 2)) / (real) ((ki - 1 - j1) * (ki - 1 - j2) * (ki - 1 
		- j3) * (ki - 1 - j4));
	++icount;
/* L3: */
    }
/*                                                                       ABBY1044 */
    --numin;
L4:
    r__1 = -one;
    i__1 = j5 + numin + 1;
    *rac = *rac * pow_ri(&r__1, &i__1) * exp(gam[numin] - gam[numin - j1 - 1] 
	    - gam[numin - j2 - 1] - gam[numin - j3 - 1] - gam[numin - j4 - 1] 
	    - gam[j5 + 2 - numin - 1] - gam[j6 + 2 - numin - 1] - gam[j7 + 2 
	    - numin - 1] + (gam[j1 + 1 - i__ - 1] + gam[j1 + 1 - j - 1] + gam[
	    j1 + 1 - m - 1] - gam[j1 + 1] + gam[j2 + 1 - k - 1] + gam[j2 + 1 
	    - l - 1] + gam[j2 + 1 - m - 1] - gam[j2 + 1] + gam[j3 + 1 - i__ - 
	    1] + gam[j3 + 1 - k - 1] + gam[j3 + 1 - n - 1] - gam[j3 + 1] + 
	    gam[j4 + 1 - j - 1] + gam[j4 + 1 - l - 1] + gam[j4 + 1 - n - 1] - 
	    gam[j4 + 1]) / two);
/* ----IN THE REAL VERSION , SUPPRESS THE COMMENT SIGN                    ABBY1052 */
/*     IF(IBUG3.NE.1) RETURN                                             ABBY1053 */
    for (ntr = 1; ntr <= 10; ++ntr) {
	second_(&tt1);
	deltat = tt1 - tt0;
	if (deltat < .001f) {
	    goto L5;
	}
	trac += deltat;
	goto L6;
L5:
	;
    }
/*                                                                       ABBY1061 */
L6:
    return 0;
} /* dracah_ */

#undef norac
#undef trac
#undef gam
#undef i__
#undef j
#undef k
#undef l
#undef m
#undef n

#define gam (facts_1.gam)

/*                                                                       ABBY1064 */
/*                                                                       ABBY1065 */
/*     ****************                                                  ABBY1066 */
/* Subroutine */ int factt_(void)
{
    /* Initialized data */

    static real thirty = 30.f;
    static real one = 1.f;
    static real two = 2.f;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    integer i__;
    real x;

/*     ****************                                                  ABBY1068 */
/*                                                                       ABBY1069 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY1070 */
/*                                                                       ABBY1071 */
/*  ***CALCULATES THE LOGS OF FACTORIALS REQUIRED BY THE RACAH           ABBY1072 */
/*  ***COEFFICIENT ROUTINE DRACAH                                        ABBY1073 */
/*  *** WRITTEN BY N.S. SCOTT                                            ABBY1074 */
/*                                                                       ABBY1075 */
/*                                                                       ABBY1079 */
/*                                                                       ABBY1080 */
    gam[0] = one;
    gam[1] = one;
    x = two;
/*                                                                       ABBY1084 */
    for (i__ = 3; i__ <= 30; ++i__) {
	gam[i__ - 1] = gam[i__ - 2] * x;
	x += one;
/* L10: */
    }
/*                                                                       ABBY1089 */
    for (i__ = 1; i__ <= 30; ++i__) {
	gam[i__ - 1] = log(gam[i__ - 1]);
/* L20: */
    }
/*                                                                       ABBY1093 */
    x = thirty;
/*                                                                       ABBY1095 */
    for (i__ = 31; i__ <= 100; ++i__) {
	gam[i__ - 1] = gam[i__ - 2] + log(x);
	x += one;
/* L30: */
    }
/*                                                                       ABBY1100 */
    return 0;
} /* factt_ */

#undef gam

#define j12 (graph_2.j12)
#define nogen (timeg_1.nogen)
#define tgen (timeg_1.tgen)
#define m (couple_1.m)
#define j1 (couple_1.j1)
#define j6c (argu_2.j6c)
#define j7c (argu_2.j7c)
#define j8c (argu_2.j8c)
#define j9c (argu_2.j9c)
#define jwc (argu_2.jwc)
#define j6 (argu_2.j6)
#define j7 (argu_2.j7)
#define j8 (argu_2.j8)
#define j9 (argu_2.j9)
#define jw (argu_2.jw)
#define jdel (argu_2.jdel)
#define ldel (argu_2.ldel)
#define mp (argu_2.mp)
#define ibug3 (debug_1.ibug3)
#define ist (jrac_2.ist)
#define ipd (inform_1.ipd)
#define j6p (sumarg_1.j6p)
#define j7p (sumarg_1.j7p)
#define j8p (sumarg_1.j8p)
#define j9p (sumarg_1.j9p)
#define jword (sumarg_1.jword)
#define nlsum (sumarg_1.nlsum)
#define nbj (sumarg_1.nbj)
#define nb6j (sumarg_1.nb6j)
#define k6cp (sumarg_1.k6cp)
#define k7cp (sumarg_1.k7cp)
#define k8cp (sumarg_1.k8cp)
#define k9cp (sumarg_1.k9cp)
#define jsum6 (sumarg_1.jsum6)
#define jsum4 (sumarg_1.jsum4)
#define jsum5 (sumarg_1.jsum5)
#define inv6j (sumarg_1.inv6j)

/*                                                                       ABBY1103 */
/*                                                                       ABBY1104 */
/*     ************************                                          ABBY1105 */
/* Subroutine */ int gensum_(real *recup)
{
    /* Initialized data */

    static real zero = 0.f;
    static real one = 1.f;
    static real epsil = 1e-10f;
    static integer mxcsvr = 4;

    /* Format strings */
    static char fmt_302[] = "(\002   SUM NR.\002,i3)";
    static char fmt_303[] = "(\002 NO SUMMATION. RECOUPLING COEFFICIENT=\002"
	    ",g15.8)";
    static char fmt_304[] = "(\002 RECOUPLING COEFFICIENT=\002,g15.8)";
    static char fmt_305[] = "(6f5.1,10x,g15.8)";
    static char fmt_306[] = "(\002 NUMBER OF INDEPENDENT SUMS:\002,i3)";
    static char fmt_307[] = "(\002 SUM NR.\002,i2,\002 SUM VALUE=\002,g15.8"
	    ",\002 RECUP=\002,g15.8)";
    static char fmt_309[] = "(\002   NOT INVOLVING SUMMATION VARIABLE\002)";
    static char fmt_400[] = "(//\002 PRINT OUT FROM SUBROUTINE GENSUM\002/"
	    "/\002 VALUES OF ANGULAR MOMENTA IN *REAL* FORMAT\002/,(14f5.1))";
    static char fmt_401[] = "(/\002 RACAH W FUNCTIONS(6J)\002/\002 ARGUMENTS"
	    " IN *REAL* FORMAT\002,18x,\002VALUE\002)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sqrt(doublereal);

    /* Local variables */
    integer i__, j, k, i1, i2, i3, i4, i5, l1, l3, l2;
    real x1;
    integer k1, jb, jj, ik, km, ll, mm, ip, nj, j6f, j7f, j8f, j9f, jj1, jj2, 
	    i1t, jt1;
    real xj1[60];
    integer jt2, jt4, ik1;
    real tt0, tt1;
    integer ik2, no1, ix2, jnd, ias, mat[144]	/* was [12][12] */, jkm, jmn, 
	    nfs, jwj, jmx, nps, jwr;
    real spr, sqr;
    integer j6cp, j7cp, j8cp, j9cp;
    logical noel[12];
    integer jmin, jmax, jmnp[5], jwrd, nolp, jmxp[5], jsum[40]	/* was [2][20]
	     */, nsum;
    real stor;
    integer ij6cp, jmin1, jmax1, jsum2[12], jsum3[12], nsum1;
    real stor1;
    logical ldiag[12];
    integer ichan, ipair[4]	/* was [2][2] */, maxlp[12];
    real wstor[20];
    extern /* Subroutine */ int dracah_(real *), second_(real *);
    integer iastor, jwtest[20];

    /* Fortran I/O blocks */
    static cilist io___121 = { 0, 0, 0, fmt_400, 0 };
    static cilist io___122 = { 0, 0, 0, fmt_306, 0 };
    static cilist io___123 = { 0, 0, 0, fmt_401, 0 };
    static cilist io___128 = { 0, 0, 0, fmt_309, 0 };
    static cilist io___131 = { 0, 0, 0, fmt_305, 0 };
    static cilist io___138 = { 0, 0, 0, fmt_303, 0 };
    static cilist io___146 = { 0, 0, 0, fmt_302, 0 };
    static cilist io___198 = { 0, 0, 0, fmt_305, 0 };
    static cilist io___204 = { 0, 0, 0, fmt_307, 0 };
    static cilist io___205 = { 0, 0, 0, fmt_304, 0 };



#define j12_ref(a_1,a_2,a_3) j12[((a_3)*12 + (a_2))*4 + a_1 - 53]
#define jw_ref(a_1,a_2) jw[(a_2)*6 + a_1 - 7]
#define mat_ref(a_1,a_2) mat[(a_2)*12 + a_1 - 13]
#define ldel_ref(a_1,a_2) ldel[(a_2)*20 + a_1 - 21]
#define jsum_ref(a_1,a_2) jsum[(a_2)*2 + a_1 - 3]
#define jsum4_ref(a_1,a_2) jsum4[(a_2)*12 + a_1 - 13]
#define jsum5_ref(a_1,a_2) jsum5[(a_2)*12 + a_1 - 13]
#define ipair_ref(a_1,a_2) ipair[(a_2)*2 + a_1 - 3]
#define jword_ref(a_1,a_2) jword[(a_2)*6 + a_1 - 7]

/*     ************************                                          ABBY1107 */
/*                                                                       ABBY1108 */
/*                                                                       ABBY1109 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY1110 */
/*                                                                       ABBY1111 */
/*  ***CARRIES OUT THE SUMMATION OVER COEFFICIENTS DEFINED BY THE ARRAYS ABBY1112 */
/*  ***THE ARRAYS J6,J7,J8,LDEL AND JW TO GIVE RECUP                     ABBY1113 */
/*  ***THE ENTRY IS EITHER MADE FROM NJGRAF OR DIRECTLY ASSUMING THAT THEABBY1114 */
/*  ***ARRAYS J6,...,JW HAVE ALREADY BEEN DETERMINED BY A PREVIOUS       ABBY1115 */
/*  ***ENTRY TO NJGRAF AND THAT THE SUMMATION IS REQUIRED FOR ANOTHER SETABBY1116 */
/*  ***OF J VALUES DEFINED BY THE ARRAY J1                               ABBY1117 */
/*                                                                       ABBY1118 */
/*  ***RECUP IS THE RECOUPLING COEFFICIENT                               ABBY1119 */
/*                                                                       ABBY1120 */
/*  ***SUBROUTINE CALLED: DRACAH                                         ABBY1121 */
/*                                                                       ABBY1122 */
/*                                                                       ABBY1126 */
/*                                                                       ABBY1128 */
/*                                                                       ABBY1136 */
/*                                                                       ABBY1141 */
/*                                                                       ABBY1144 */
/*                                                                       ABBY1148 */
/*                                                                       ABBY1150 */
/*                                                                       ABBY1154 */
/*  ***FORMAT STATEMENTS USED IN GENSUM                                  ABBY1155 */
/*                                                                       ABBY1156 */
/* L302: */
/* L303: */
/* L304: */
/* L305: */
/* L306: */
/* L307: */
/* L308: */
/* L309: */
/* L400: */
/* L401: */
/*                                                                       ABBY1169 */
/*                                                                       ABBY1170 */
/*                                                                       ABBY1171 */
/* ***EVALUATES ALL TERMS IN J6,J7,J8,J9,LDEL,AND JW WHICH DO NOT INVOLVEABBY1172 */
/*   ***A SUMMATION.THE RESULT IS STORED IN RECUP AND IASTOR             ABBY1173 */
/*                                                                       ABBY1174 */
    ++nogen;
    second_(&tt0);
    if (ibug3 == 1) {
/*                                                                       ABBY1178 */
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xj1[i__ - 1] = (j1[i__ - 1] - 1.f) / 2.f;
/* L139: */
	}
/*                                                                       ABBY1182 */
	io___121.ciunit = ipd;
	s_wsfe(&io___121);
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&xj1[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	io___122.ciunit = ipd;
	s_wsfe(&io___122);
	do_fio(&c__1, (char *)&nlsum, (ftnlen)sizeof(integer));
	e_wsfe();
	io___123.ciunit = ipd;
	s_wsfe(&io___123);
	e_wsfe();
    }
/*                                                                       ABBY1187 */
    mm = m + 1;
    j1[mm - 1] = 1;
/*                                                                       ABBY1190 */
/*  ***TEST DELTA FUNCTIONS                                              ABBY1191 */
/*                                                                       ABBY1192 */
    j1[mm - 1] = 1;
    if (jdel <= 0) {
	goto L180;
    }
/*                                                                       ABBY1195 */
    i__1 = jdel;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = ldel_ref(i__, 1);
	i2 = ldel_ref(i__, 2);
	if (i1 > mm || i2 > mm) {
	    if (i1 > mm) {
		j1[i1 - 1] = j1[i2 - 1];
	    }
	    if (i2 > mm) {
		j1[i2 - 1] = j1[i1 - 1];
	    }
	} else {
	    if (j1[i1 - 1] != j1[i2 - 1]) {
		second_(&tt1);
		*recup = zero;
		tgen = tgen + tt1 - tt0;
		return 0;
	    }
	}
/* L141: */
    }
/*                                                                       ABBY1211 */
L180:
    *recup = one;
    if (jwc != 0) {
/*                                                                       ABBY1214 */
/*  ***MULTIPLY RECUP BY ALL RACAH COEFFICIENTS WHICH DO NOT INVOLVE A   ABBY1215 */
/*  ***SUMMATION                                                         ABBY1216 */
/*                                                                       ABBY1217 */
	if (ibug3 == 1) {
	    io___128.ciunit = ipd;
	    s_wsfe(&io___128);
	    e_wsfe();
	}
/*                                                                       ABBY1219 */
	i__1 = jwc;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (inv6j[i__ - 1] > 0) {
		goto L7;
	    }
	    for (j = 1; j <= 6; ++j) {
		i1 = jw_ref(j, i__);
		ist[j - 1] = j1[i1 - 1] - 1;
/* L3: */
	    }
/*                                                                       ABBY1226 */
	    dracah_(&x1);
	    if (ibug3 == 1) {
		io___131.ciunit = ipd;
		s_wsfe(&io___131);
		for (k = 1; k <= 6; ++k) {
		    do_fio(&c__1, (char *)&xj1[jw_ref(k, i__) - 1], (ftnlen)
			    sizeof(real));
		}
		do_fio(&c__1, (char *)&x1, (ftnlen)sizeof(real));
		e_wsfe();
	    }
	    *recup *= x1;
/*                                                                       ABBY1230 */
L7:
	    ;
	}
/*                                                                       ABBY1232 */
    }
/*                                                                       ABBY1234 */
    sqr = 1.f;
/*                                                                       ABBY1236 */
    if (j6c != 0) {
	i__1 = j6c;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i1 = j6[i__ - 1];
	    sqr *= j1[i1 - 1];
/* L12: */
	}
    }
/*                                                                       ABBY1243 */
    spr = 1.f;
/*                                                                       ABBY1245 */
    if (j9c != 0) {
	i__1 = j9c;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i1 = j9[i__ - 1];
	    spr *= j1[i1 - 1];
/* L144: */
	}
    }
/*                                                                       ABBY1252 */
    *recup *= sqrt(sqr / spr);
    if (dabs(*recup) < epsil) {
	goto L145;
    }
    iastor = 0;
/*                                                                       ABBY1256 */
    if (j7c != 0) {
	i__1 = j7c;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i1 = j7[i__ - 1];
	    iastor = iastor + j1[i1 - 1] - 1;
/* L17: */
	}
    }
/*                                                                       ABBY1263 */
    if (j8c != 0) {
	i__1 = j8c;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i1 = j8[i__ - 1];
	    iastor += (j1[i1 - 1] - 1) << 1;
/* L22: */
	}
    }
/*                                                                       ABBY1270 */
    if (nlsum <= 0) {
	iastor /= 2;
/*                                                                       ABBY1273 */
/*  ***NO SUMMATION INVOLVED.END OF COMPUTATION                          ABBY1274 */
/*                                                                       ABBY1275 */
	stor1 = one;
	stor = one;
	if (iastor % 2 == 1) {
	    *recup = -(*recup);
	}
	if (ibug3 == 1) {
	    io___138.ciunit = ipd;
	    s_wsfe(&io___138);
	    do_fio(&c__1, (char *)&(*recup), (ftnlen)sizeof(real));
	    e_wsfe();
	}
	second_(&tt1);
	tgen = tgen + tt1 - tt0;
	return 0;
/*                                                                       ABBY1283 */
    }
/*                                                                       ABBY1285 */
/*                                                                       ABBY1286 */
/*  ***EVALUATION OF THE PART INVOLVING SUMMATIONS.                      ABBY1287 */
/*                                                                       ABBY1288 */
/*                                                                       ABBY1289 */
    nfs = 0;
    jwr = 0;
    j6f = 0;
    j7f = 0;
    j8f = 0;
    j9f = 0;
    nps = 0;
L25:
    ++nps;
    if (ibug3 == 1) {
	io___146.ciunit = ipd;
	s_wsfe(&io___146);
	do_fio(&c__1, (char *)&nps, (ftnlen)sizeof(integer));
	e_wsfe();
    }
/*                                                                       ABBY1299 */
/*                                                                       ABBY1300 */
/*   *** LOOP ON THE DISCONNECTED SUMMATIONS                             ABBY1301 */
/*                                                                       ABBY1302 */
/*                                                                       ABBY1303 */
    ias = 0;
    nsum = nbj[nps - 1] - nfs;
    jwrd = nb6j[nps - 1] - jwr;
    j6cp = k6cp[nps - 1];
    j7cp = k7cp[nps - 1];
    j8cp = k8cp[nps - 1];
    j9cp = k9cp[nps - 1];
/*                                                                       ABBY1311 */
/*                                                                       ABBY1312 */
/*     ***THE RANGE OF VALUES OF EACH SUMMATION VARIABLE IS              ABBY1313 */
/*     ***DEFINED BY ESTABLISHING A MATRIX OF THE LINKS BETWEEN          ABBY1314 */
/*     ***VARIABLES.MAT(I,J) CONTAINS:                                   ABBY1315 */
/*        I=J    NUMBER OF POSSIBLE VALUES OF I DUE TO TRIANGULAR        ABBY1316 */
/*               RELATIONS WITH NON-VARIABLES,I.E. CONSTANTS.            ABBY1317 */
/*       I.GT.J  NUMBER OF LINKS BETWEEN I AND J THROUGH CONSTANTS       ABBY1318 */
/*       I.LT.J  VALUE OF THE CONSTANT,IF THE ABOVE IS 1.IF NOT,         ABBY1319 */
/*               THESE VALUES ARE SRORED IN J12(L,I,J) WHERE THERE       ABBY1320 */
/*               IS ROOM FOR MXCSVR SUCH VALUES (L.LE.4)                 ABBY1321 */
/*                                                                       ABBY1322 */
/*                                                                       ABBY1323 */
    i__1 = nsum;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nsum;
	for (j = 1; j <= i__2; ++j) {
	    mat_ref(i__, j) = 0;
/* L152: */
	}
/* L52: */
    }
/*                                                                       ABBY1329 */
    i__1 = nsum;
    for (i1 = 1; i1 <= i__1; ++i1) {
	i1t = i1 + nfs;
	i2 = jsum6[i1t - 1];
	i__2 = i2;
	for (i3 = 1; i3 <= i__2; ++i3) {
	    i__ = jsum5_ref(i1t, i3);
	    j = jsum4_ref(i1t, i3);
	    switch (j) {
		case 1:  goto L54;
		case 2:  goto L55;
		case 3:  goto L56;
		case 4:  goto L57;
		case 5:  goto L58;
		case 6:  goto L59;
	    }
/*                                                                       ABBY1337 */
/*  ***THE ROWS OF THE IPAIR ARRAYS GIVE LIMITS OF SUMMATION IMPOSED     ABBY1338 */
/*                                                                       ABBY1339 */
/*                                                                       ABBY1340 */
L54:
	    ipair_ref(1, 1) = jword_ref(2, i__);
	    ipair_ref(1, 2) = jword_ref(5, i__);
	    ipair_ref(2, 1) = jword_ref(3, i__);
	    ipair_ref(2, 2) = jword_ref(6, i__);
	    goto L60;
/*                                                                       ABBY1346 */
L55:
	    ipair_ref(1, 1) = jword_ref(1, i__);
	    ipair_ref(1, 2) = jword_ref(5, i__);
	    ipair_ref(2, 1) = jword_ref(4, i__);
	    ipair_ref(2, 2) = jword_ref(6, i__);
	    goto L60;
/*                                                                       ABBY1352 */
L56:
	    ipair_ref(1, 1) = jword_ref(1, i__);
	    ipair_ref(1, 2) = jword_ref(6, i__);
	    ipair_ref(2, 1) = jword_ref(4, i__);
	    ipair_ref(2, 2) = jword_ref(5, i__);
	    goto L60;
/*                                                                       ABBY1358 */
L57:
	    ipair_ref(1, 1) = jword_ref(2, i__);
	    ipair_ref(1, 2) = jword_ref(6, i__);
	    ipair_ref(2, 1) = jword_ref(3, i__);
	    ipair_ref(2, 2) = jword_ref(5, i__);
	    goto L60;
/*                                                                       ABBY1364 */
L58:
	    ipair_ref(1, 1) = jword_ref(1, i__);
	    ipair_ref(1, 2) = jword_ref(2, i__);
	    ipair_ref(2, 1) = jword_ref(3, i__);
	    ipair_ref(2, 2) = jword_ref(4, i__);
	    goto L60;
/*                                                                       ABBY1370 */
L59:
	    ipair_ref(1, 1) = jword_ref(1, i__);
	    ipair_ref(1, 2) = jword_ref(3, i__);
	    ipair_ref(2, 1) = jword_ref(2, i__);
	    ipair_ref(2, 2) = jword_ref(4, i__);
/*                                                                       ABBY1375 */
L60:
	    for (i4 = 1; i4 <= 2; ++i4) {
		km = 0;
		for (i5 = 1; i5 <= 2; ++i5) {
		    if (ipair_ref(i4, i5) > mp) {
			++km;
		    }
/* L62: */
		}
/*                                                                       ABBY1381 */
		jj1 = ipair_ref(i4, 1);
		jj2 = ipair_ref(i4, 2);
		if (km == 1) {
		    goto L67;
		}
		if (km > 1) {
		    goto L63;
		}
/*                                                                       ABBY1386 */
/*  ***ONE VARIABLE LINKED TO TWO CONSTANTS.FIX THE DIAGONAL MAT(I,I)    ABBY1387 */
/*                                                                       ABBY1388 */
		jt1 = j1[jj1 - 1] - 1;
		jt2 = j1[jj2 - 1] - 1;
		jmin = (i__3 = jt1 - jt2, abs(i__3));
		jmax = jt1 + jt2;
/*                                                                       ABBY1393 */
		if (mat_ref(i1, i1) > 1) {
/*                                                                       ABBY1395 */
/*  ***IF THERE ARE SEVERAL COUPLES OF CONSTANTS ,TAKE THE MORE          ABBY1396 */
/*  ***STRINGENT COMBINATION                                             ABBY1397 */
/*                                                                       ABBY1398 */
/* Computing MAX */
		    i__3 = jmin, i__4 = jsum_ref(1, i1);
		    jmin = max(i__3,i__4);
/* Computing MIN */
		    i__3 = jmax, i__4 = jsum_ref(2, i1);
		    jmax = min(i__3,i__4);
		    if (jmax >= jmin) {
			jsum_ref(1, i1) = jmin;
			jsum_ref(2, i1) = jmax;
			mat_ref(i1, i1) = (jmax - jmin) / 2 + 1;
			goto L63;
		    } else {
			*recup = zero;
			goto L110;
		    }
		} else {
		    if (mat_ref(i1, i1) < 1) {
/*                                                                       ABBY1412 */
/*  ***FIRST TIME                                                        ABBY1413 */
/*                                                                       ABBY1414 */
			mat_ref(i1, i1) = (jmax - jmin) / 2 + 1;
			jsum_ref(1, i1) = jmin;
			jsum_ref(2, i1) = jmax;
/*                                                                       ABBY1418 */
		    }
		}
/*                                                                       ABBY1421 */
		goto L63;
/*                                                                       ABBY1423 */
/*                                                                       ABBY1424 */
/*                                                                       ABBY1425 */
/*                                                                       ABBY1426 */
/*  ***ONE VARIABLE LINKED TO ONE CONSTANT AND ONE VARIABLE  NON DIAGONALABBY1427 */
/*  ***ELEMENT                                                           ABBY1428 */
/*                                                                       ABBY1429 */
L67:
		jt1 = min(jj1,jj2);
		jt2 = max(jj1,jj2) - mp;
		if (jt2 > i1) {
		    goto L63;
		}
		jt4 = j1[jt1 - 1] - 1;
		k = mat_ref(i1, jt2);
		if (k == 0) {
		    goto L107;
		}
/*                                                                       ABBY1436 */
		i__3 = k;
		for (ll = 1; ll <= i__3; ++ll) {
		    if (jt4 == j12_ref(ll, jt2, i1)) {
			goto L63;
		    }
/* L71: */
		}
/*                                                                       ABBY1440 */
L107:
		++k;
		if (k > mxcsvr) {
		    goto L63;
		}
		mat_ref(i1, jt2) = k;
		j12_ref(k, jt2, i1) = jt4;
/*                                                                       ABBY1445 */
L63:
		;
	    }
/* L65: */
	}
/* L66: */
    }
/*                                                                       ABBY1449 */
/*  ***REDUCE THE DIAGONAL ELEMENTS BY TAKING INTO ACCOUNT THE NON       ABBY1450 */
/*  ***DIAGONAL ELEMENTS,AND KEEP THE LATTER ONLY IF NEEDED              ABBY1451 */
/*                                                                       ABBY1452 */
L150:
    ichan = 0;
/*                                                                       ABBY1454 */
    i__1 = nsum;
    for (i__ = 1; i__ <= i__1; ++i__) {
	noel[i__ - 1] = TRUE_;
	i1 = i__ - 1;
	if (i1 == 0) {
	    goto L170;
	}
	i__2 = i1;
	for (j = 1; j <= i__2; ++j) {
	    if (mat_ref(i__, j) == 0 || mat_ref(j, j) == 0) {
		goto L72;
	    }
	    ik1 = i__;
	    ik2 = j;
	    jkm = 1;
	    goto L201;
L171:
	    noel[i__ - 1] = FALSE_;
L72:
	    ;
	}
L170:
	if (i__ == nsum) {
	    goto L74;
	}
	i2 = i__ + 1;
/*                                                                       ABBY1470 */
	i__2 = nsum;
	for (j = i2; j <= i__2; ++j) {
	    if (mat_ref(j, i__) == 0 || mat_ref(j, j) == 0) {
		goto L73;
	    }
	    ik1 = j;
	    ik2 = i__;
	    jkm = 2;
	    goto L201;
L172:
L73:
	    ;
	}
L74:
	;
    }
/*                                                                       ABBY1480 */
    if (ichan != 0) {
	goto L150;
    }
    goto L220;
/*                                                                       ABBY1483 */
/*                                                                       ABBY1484 */
L201:
    jmin1 = 0;
    jmax1 = 1000;
    k = mat_ref(ik1, ik2);
/*                                                                       ABBY1488 */
    i__1 = k;
    for (l1 = 1; l1 <= i__1; ++l1) {
/*                                                                       ABBY1490 */
	l3 = mat_ref(j, j);
	jj1 = jsum_ref(1, j);
	jnd = j12_ref(l1, ik2, ik1);
	jmin = 1000;
	jmax = 0;
	jmnp[l1 - 1] = 0;
	jmxp[l1 - 1] = 1000;
/*                                                                       ABBY1498 */
	i__2 = l3;
	for (l2 = 1; l2 <= i__2; ++l2) {
/*                                                                       ABBY1500 */
	    jmn = (i__3 = jnd - jj1, abs(i__3));
	    jmx = jnd + jj1;
	    jmin = min(jmn,jmin);
	    jmax = max(jmx,jmax);
/* Computing MAX */
	    i__3 = jmn, i__4 = jmnp[l1 - 1];
	    jmnp[l1 - 1] = max(i__3,i__4);
/* Computing MIN */
	    i__3 = jmx, i__4 = jmxp[l1 - 1];
	    jmxp[l1 - 1] = min(i__3,i__4);
	    jj1 += 2;
/*                                                                       ABBY1508 */
/* L204: */
	}
/*                                                                       ABBY1510 */
	jmin1 = max(jmin1,jmin);
	jmax1 = min(jmax1,jmax);
/*                                                                       ABBY1513 */
/* L203: */
    }
/*                                                                       ABBY1515 */
    if (mat_ref(i__, i__) == 0) {
	jsum_ref(1, i__) = jmin1;
	jsum_ref(2, i__) = jmax1;
	mat_ref(i__, i__) = (jmax1 - jmin1) / 2 + 1;
	++ichan;
	goto L206;
    }
/*                                                                       ABBY1523 */
    if (jsum_ref(1, i__) < jmin1) {
	jsum_ref(1, i__) = jmin1;
	++ichan;
    }
/*                                                                       ABBY1528 */
    if (jsum_ref(2, i__) > jmax1) {
	jsum_ref(2, i__) = jmax1;
	++ichan;
    }
/*                                                                       ABBY1533 */
L206:
    k1 = 0;
/*                                                                       ABBY1535 */
    i__1 = k;
    for (l1 = 1; l1 <= i__1; ++l1) {
	if (jmnp[l1 - 1] <= jsum_ref(1, i__) && jmxp[l1 - 1] >= jsum_ref(2, 
		i__)) {
	    goto L207;
	}
	++k1;
	j12_ref(k1, ik2, ik1) = j12_ref(l1, ik2, ik1);
L207:
	;
    }
/*                                                                       ABBY1541 */
    if (k1 != k) {
	mat_ref(ik1, ik2) = k1;
	++ichan;
    }
/*                                                                       ABBY1546 */
    mat_ref(ik2, ik1) = j12_ref(1, ik2, ik1);
    switch (jkm) {
	case 1:  goto L171;
	case 2:  goto L172;
    }
/*                                                                       ABBY1549 */
/*                                                                       ABBY1550 */
/*                                                                       ABBY1551 */
/*  ***CARRY OUT THE SUMMATIONS.                                         ABBY1552 */
/*                                                                       ABBY1553 */
L220:
    i__1 = nsum;
    for (i__ = 1; i__ <= i__1; ++i__) {
	jsum3[i__ - 1] = 1;
	ldiag[i__ - 1] = FALSE_;
	if (mat_ref(i__, i__) == 1) {
	    ldiag[i__ - 1] = TRUE_;
	}
/* L230: */
    }
/*                                                                       ABBY1559 */
    i__1 = jwrd;
    for (i__ = 1; i__ <= i__1; ++i__) {
	jwtest[i__ - 1] = 1;
/* L231: */
    }
/*                                                                       ABBY1563 */
    stor = zero;
    stor1 = one;
    nolp = 0;
    ip = 1;
L240:
    ++nolp;
/*                                                                       ABBY1569 */
/*                                                                       ABBY1570 */
/*  ***FIND THE RANGE OF JSUM2(NOLP)                                     ABBY1571 */
/*  ***NOLP IS THE INDEX  OF THE SUMMATION VARIABLE                      ABBY1572 */
/*                                                                       ABBY1573 */
    jmin = jsum_ref(1, nolp);
    jmax = jsum_ref(2, nolp);
    if (noel[nolp - 1]) {
	goto L241;
    }
    no1 = nolp - 1;
/*                                                                       ABBY1578 */
    i__1 = no1;
    for (nj = 1; nj <= i__1; ++nj) {
	if (mat_ref(nolp, nj) == 1) {
	    jj1 = mat_ref(nj, nolp);
	    jj2 = jsum2[nj - 1];
/* Computing MAX */
	    i__3 = jmin, i__4 = (i__2 = jj2 - jj1, abs(i__2));
	    jmin = max(i__3,i__4);
/* Computing MIN */
	    i__2 = jmax, i__3 = jj1 + jj2;
	    jmax = min(i__2,i__3);
	} else {
	    if (mat_ref(nolp, nj) > 1) {
		k = mat_ref(nolp, nj);
		jj2 = jsum2[nj - 1];
/*                                                                       ABBY1589 */
		i__2 = k;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    jj1 = j12_ref(i__, nj, nolp);
/* Computing MAX */
		    i__4 = jmin, i__5 = (i__3 = jj2 - jj1, abs(i__3));
		    jmin = max(i__4,i__5);
/* Computing MIN */
		    i__3 = jmax, i__4 = jj1 + jj2;
		    jmax = min(i__3,i__4);
/* L245: */
		}
/*                                                                       ABBY1595 */
	    }
	}
/*                                                                       ABBY1598 */
/* L242: */
    }
/*                                                                       ABBY1600 */
L241:
    jsum2[nolp - 1] = jmin;
    maxlp[nolp - 1] = jmax;
    if (ldiag[nolp - 1]) {
	jsum3[nolp - 1] = 0;
    }
    if (nolp < nsum) {
	goto L240;
    }
/*                                                                       ABBY1605 */
    i__1 = jmax;
    for (jj = jmin; jj <= i__1; jj += 2) {
	jsum2[nsum - 1] = jj;
/*                                                                       ABBY1608 */
/*                                                                       ABBY1609 */
/*  ***DETERMINE WHICH RACAH COEFFICIENTS NEED RE-EVALUATING AND         ABBY1610 */
/*  ***SET JWTEST APPROPRIATELY                                          ABBY1611 */
/*                                                                       ABBY1612 */
	i__2 = nsum;
	for (j = ip; j <= i__2; ++j) {
	    if (jsum3[j - 1] <= 0) {
		goto L114;
	    }
	    i2 = jsum6[j - 1];
/*                                                                       ABBY1616 */
	    i__3 = i2;
	    for (i1 = 1; i1 <= i__3; ++i1) {
		i3 = jsum5_ref(j, i1);
		jwtest[i3 - 1] = 1;
/*                                                                       ABBY1620 */
/* L113: */
	    }
/*                                                                       ABBY1622 */
L114:
	    ;
	}
/*                                                                       ABBY1624 */
	i__2 = jwrd;
	for (j = 1; j <= i__2; ++j) {
	    if (jwtest[j - 1] == 0) {
		goto L98;
	    }
	    jwj = j + jwr;
/*                                                                       ABBY1628 */
	    for (i__ = 1; i__ <= 6; ++i__) {
		if (jword_ref(i__, jwj) <= mp) {
		    i1 = jword_ref(i__, jwj);
		    ist[i__ - 1] = j1[i1 - 1] - 1;
		} else {
		    i1 = jword_ref(i__, jwj) - mp - nfs;
		    ist[i__ - 1] = jsum2[i1 - 1];
		}
/* L90: */
	    }
/*                                                                       ABBY1638 */
	    dracah_(&x1);
	    wstor[j - 1] = x1;
	    if (ibug3 == 1) {
		for (i__ = 1; i__ <= 6; ++i__) {
		    xj1[i__ - 1] = ist[i__ - 1] / 2.f;
/* L99: */
		}
/*                                                                       ABBY1645 */
		io___198.ciunit = ipd;
		s_wsfe(&io___198);
		for (i__ = 1; i__ <= 6; ++i__) {
		    do_fio(&c__1, (char *)&xj1[i__ - 1], (ftnlen)sizeof(real))
			    ;
		}
		do_fio(&c__1, (char *)&x1, (ftnlen)sizeof(real));
		e_wsfe();
	    }
L98:
	    ;
	}
/*                                                                       ABBY1649 */
/*                                                                       ABBY1650 */
/*  ***FORM PRODUCT OF RACAH COEFFICIENTS,(2J+1) FACTORS AND (-1)        ABBY1651 */
/*  ***FACTORS IN STOR1                                                  ABBY1652 */
/*                                                                       ABBY1653 */
	i__2 = jwrd;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    stor1 *= wstor[i__ - 1];
/* L126: */
	}
/*                                                                       ABBY1657 */
/*  ***IASTOR CONTAINS THE POWER OF (-1)WHICH IS COMMON TO ALL TERMS     ABBY1658 */
/*                                                                       ABBY1659 */
	ix2 = 0;
	ij6cp = 1;
	if (j6cp != j6f) {
	    jb = j6f + 1;
/*                                                                       ABBY1664 */
	    i__2 = j6cp;
	    for (i__ = jb; i__ <= i__2; ++i__) {
		i1 = j6p[i__ - 1] - nfs;
		ij6cp *= jsum2[i1 - 1] + 1;
/* L128: */
	    }
	}
/*                                                                       ABBY1670 */
	if (j9cp != j9f) {
	    jb = j9f + 1;
/*                                                                       ABBY1673 */
	    i__2 = j9cp;
	    for (i__ = jb; i__ <= i__2; ++i__) {
		i1 = j9p[i__ - 1] - nfs;
		ij6cp /= jsum2[i1 - 1] + 1;
/* L147: */
	    }
	}
/*                                                                       ABBY1679 */
	stor1 *= sqrt((real) ij6cp);
/*                                                                       ABBY1681 */
	if (j7cp != j7f) {
	    jb = j7f + 1;
/*                                                                       ABBY1684 */
	    i__2 = j7cp;
	    for (i__ = jb; i__ <= i__2; ++i__) {
		i1 = j7p[i__ - 1] - nfs;
		ix2 += jsum2[i1 - 1];
/* L131: */
	    }
	}
/*                                                                       ABBY1690 */
	if (j8cp != j8f) {
	    jb = j8f + 1;
/*                                                                       ABBY1693 */
	    i__2 = j8cp;
	    for (i__ = jb; i__ <= i__2; ++i__) {
		i1 = j8p[i__ - 1] - nfs;
		ix2 += jsum2[i1 - 1] << 1;
/* L134: */
	    }
	}
/*                                                                       ABBY1699 */
	if (ix2 % 2 == 1) {
	    ias = -1;
	    ++ix2;
	}
/*                                                                       ABBY1704 */
	ix2 /= 2;
/*                                                                       ABBY1706 */
/*                                                                       ABBY1707 */
/*  ***ADD TERM INTO STOR AND RESET STOR1 TO 1 READY FOR NEXT TERM       ABBY1708 */
/*                                                                       ABBY1709 */
	if (ix2 % 2 == 1) {
	    stor1 = -stor1;
	}
	stor += stor1;
	stor1 = one;
	nsum1 = nsum - 1;
	if (nsum1 == 0) {
	    goto L260;
	}
/*                                                                       ABBY1715 */
	i__2 = nsum1;
	for (ik = 1; ik <= i__2; ++ik) {
	    jsum3[ik - 1] = 0;
/* L261: */
	}
/*                                                                       ABBY1719 */
	i__2 = jwrd;
	for (ik = 1; ik <= i__2; ++ik) {
	    jwtest[ik - 1] = 0;
/* L262: */
	}
/*                                                                       ABBY1723 */
L260:
	;
    }
/*                                                                       ABBY1725 */
L250:
    --nolp;
/*                                                                       ABBY1727 */
    if (nolp != 0) {
	if (ldiag[nolp - 1]) {
	    goto L250;
	}
	jsum3[nolp - 1] = 1;
	jsum2[nolp - 1] += 2;
	if (jsum2[nolp - 1] > maxlp[nolp - 1]) {
	    goto L250;
	}
	ip = nolp;
/*                                                                       ABBY1734 */
/*                                                                       ABBY1735 */
/*    ***PROCEED TO NEXT VARIABLE                                        ABBY1736 */
/*                                                                       ABBY1737 */
	goto L240;
/*                                                                       ABBY1739 */
    }
/*                                                                       ABBY1741 */
    *recup *= stor;
    if (ibug3 == 1) {
	io___204.ciunit = ipd;
	s_wsfe(&io___204);
	do_fio(&c__1, (char *)&nps, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&stor, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*recup), (ftnlen)sizeof(real));
	e_wsfe();
    }
    if (dabs(*recup) < epsil) {
	goto L145;
    }
    jwr = jwrd + jwr;
    nfs = nsum + nfs;
    j6f = j6cp;
    j7f = j7cp;
    j8f = j8cp;
    j9f = j9cp;
    iastor += ias;
/*                                                                       ABBY1752 */
/*                                                                       ABBY1753 */
/*  ***PROCEED TO NEXT SUM                                               ABBY1754 */
/*                                                                       ABBY1755 */
    if (nps < nlsum) {
	goto L25;
    }
    iastor /= 2;
    if (iastor % 2 != 0) {
	*recup = -(*recup);
    }
    if (ibug3 == 1) {
	io___205.ciunit = ipd;
	s_wsfe(&io___205);
	do_fio(&c__1, (char *)&(*recup), (ftnlen)sizeof(real));
	e_wsfe();
    }
    second_(&tt1);
    tgen = tgen + tt1 - tt0;
L110:
    return 0;
/*                                                                       ABBY1763 */
/*  ***NO SUMMATIONS. CHECK THAT THERE ARE NO INCONSISTENCIES. THEN      ABBY1764 */
/*  ***MULTIPLY BY (-1) FACTOR AND EXIT                                  ABBY1765 */
/*                                                                       ABBY1766 */
L145:
    *recup = zero;
    return 0;
} /* gensum_ */

#undef j12
#undef nogen
#undef tgen
#undef m
#undef j1
#undef j6c
#undef j7c
#undef j8c
#undef j9c
#undef jwc
#undef j6
#undef j7
#undef j8
#undef j9
#undef jw
#undef jdel
#undef ldel
#undef mp
#undef ibug3
#undef ist
#undef ipd
#undef j6p
#undef j7p
#undef j8p
#undef j9p
#undef jword
#undef nlsum
#undef nbj
#undef nb6j
#undef k6cp
#undef k7cp
#undef k8cp
#undef k9cp
#undef jsum6
#undef jsum4
#undef jsum5
#undef inv6j


#undef jword_ref
#undef ipair_ref
#undef jsum5_ref
#undef jsum4_ref
#undef jsum_ref
#undef ldel_ref
#undef mat_ref
#undef jw_ref
#undef j12_ref

#define jdiag (graph_1.jdiag)
#define tab1 (graph_1.tab1)
#define il (graph_1.il)
#define ih (graph_1.ih)
#define ifirst (graph_1.ifirst)
#define ilast (graph_1.ilast)
#define nfree (graph_1.nfree)
#define itfree (graph_1.itfree)
#define m (couple_1.m)
#define ial (build_1.ial)

/*                                                                       ABBY1770 */
/*                                                                       ABBY1771 */
/*     ****************                                                  ABBY1772 */
/* Subroutine */ int intab_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, j, k, it, ifr, itt;


#define tab1_ref(a_1,a_2) tab1[(a_2)*60 + a_1 - 61]
#define jdiag_ref(a_1,a_2) jdiag[(a_2)*48 + a_1 - 49]

/*     ****************                                                  ABBY1774 */
/*                                                                       ABBY1775 */
/*  ***THIS SUBROUTINE CALLED AT THE END OF DIAGRM,FIXES THE ARRAYS IH   ABBY1776 */
/*  ***AND IL-SO TO SPEAK HARDWARE AND LOGICAL ADDRESSES OF TRIADS IN    ABBY1777 */
/*  ***JDIAG.ALSO DETERMINES THE NUMBER OF FREE ENDS NFREE AND THEIR     ABBY1778 */
/*  ***LOCATION ITFREE.                                                  ABBY1779 */
/*                                                                       ABBY1780 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY1781 */
/*                                                                       ABBY1782 */
/*                                                                       ABBY1785 */
/*                                                                       ABBY1787 */
/*                                                                       ABBY1789 */
/*                                                                       ABBY1793 */
/*                                                                       ABBY1796 */
/*                                                                       ABBY1797 */
/*                                                                       ABBY1798 */
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L1: */
	ial[i__ - 1] = 1;
    }
/*                                                                       ABBY1801 */
    i__1 = ilast;
    for (i__ = ifirst; i__ <= i__1; ++i__) {
	j = jdiag_ref(i__, 1);
	k = ial[j - 1];
	tab1_ref(j, k) = i__;
	ial[j - 1] = k + 1;
/* L3: */
    }
/*                                                                       ABBY1808 */
    ifr = ifirst - 1;
/*                                                                       ABBY1810 */
    i__1 = ilast;
    for (i__ = ifirst; i__ <= i__1; ++i__) {
	it = i__ - ifr;
	il[i__ - 1] = it;
	ih[it - 1] = i__;
/* L4: */
    }
/*                                                                       ABBY1816 */
    j = jdiag_ref(ifirst, 3);
    k = ial[j - 1];
    if (k > 1) {
	tab1_ref(j, 2) = tab1_ref(j, 1);
    }
    tab1_ref(j, 1) = ifirst;
    ial[j - 1] = 3;
    j = jdiag_ref(ilast, 2);
    tab1_ref(j, 2) = ilast;
    ial[j - 1] = 3;
    nfree = 0;
/*                                                                       ABBY1826 */
    i__1 = ilast;
    for (i__ = ifirst; i__ <= i__1; ++i__) {
	j = jdiag_ref(i__, 1);
	if (ial[j - 1] != 3) {
	    ++nfree;
	    itt = ilast + nfree;
	    tab1_ref(j, 2) = itt;
	    il[itt - 1] = nfree * 1000;
	    itfree[nfree - 1] = i__;
	}
/* L7: */
    }
/*                                                                       ABBY1837 */
    return 0;
} /* intab_ */

#undef jdiag
#undef tab1
#undef il
#undef ih
#undef ifirst
#undef ilast
#undef nfree
#undef itfree
#undef m
#undef ial


#undef jdiag_ref
#undef tab1_ref

#define namsub (nam_1.namsub)
#define jdiag (graph_1.jdiag)
#define arr (graph_1.arr)
#define il (graph_1.il)
#define ih (graph_1.ih)
#define npoint (graph_1.npoint)
#define nbnode (graph_1.nbnode)
#define ilast (graph_1.ilast)
#define m (couple_1.m)
#define j6c (argu_1.j6c)
#define j9c (argu_1.j9c)
#define j6 (argu_1.j6)
#define j9 (argu_1.j9)

/*                                                                       ABBY1840 */
/*                                                                       ABBY1841 */
/*     ***********************                                           ABBY1842 */
/* Subroutine */ int lolpop_(logical *fail)
{
    /* Initialized data */

    static char name__[6] = "LOLPOP";
    static integer kp[3] = { 2,3,1 };
    static integer ks[3] = { 0,1,-1 };

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    integer i__, k, l, i1, k3, k1, i2, l1, i3, k2, it, il1, il2, ilp;
    extern /* Subroutine */ int delta_(integer *, integer *, logical *), 
	    phase2_(integer *);


#define arr_ref(a_1,a_2) arr[(a_2)*48 + a_1 - 49]
#define jdiag_ref(a_1,a_2) jdiag[(a_2)*48 + a_1 - 49]

/*     ***********************                                           ABBY1844 */
/*                                                                       ABBY1845 */
/*  ***REDUCES A LOOP WITH ONE LINE AND ONE NODE IN THE FLAT GRAPH.      ABBY1846 */
/*                                                                       ABBY1847 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY1848 */
/*                                                                       ABBY1849 */
/*                                                                       ABBY1852 */
/*                                                                       ABBY1854 */
/*                                                                       ABBY1856 */
/*                                                                       ABBY1858 */
/*                                                                       ABBY1860 */
/*                                                                       ABBY1862 */
/*                                                                       ABBY1866 */
/*                                                                       ABBY1870 */
/*                                                                       ABBY1872 */
/*                                                                       ABBY1876 */
/*                                                                       ABBY1877 */
/*                                                                       ABBY1878 */
    s_copy(namsub, name__, (ftnlen)6, (ftnlen)6);
    i1 = npoint[0];
    k3 = 2;
    if (i1 == ilast) {
	k3 = 3;
    }
    l = jdiag_ref(i1, k3);
    delta_(&l, &m, fail);
    if (*fail) {
	return 0;
    }
    k = kp[k3 - 1];
    if (arr_ref(i1, k) < 0) {
	phase2_(&jdiag_ref(i1, k));
    }
    k1 = ks[k3 - 1];
    il1 = il[i1 - 1] + k1;
    i2 = ih[il1 - 1];
    l1 = jdiag_ref(i2, 1);
    delta_(&l1, &jdiag_ref(i2, k3), fail);
    if (*fail) {
	return 0;
    }
    if (arr_ref(i2, k3) == k1) {
	phase2_(&l1);
    }
    il2 = il[i2 - 1] + k1;
    i3 = ih[il2 - 1];
    k2 = k3 + k1;
    jdiag_ref(i3, k2) = l1;
    arr_ref(i3, k2) = arr_ref(i2, 1);
    ++j9c;
    j9[j9c - 1] = l1;
    ++j6c;
    j6[j6c - 1] = jdiag_ref(i1, 1);
    if (k3 == 3) {
	return 0;
    }
/*                                                                       ABBY1905 */
    i__1 = nbnode;
    for (i__ = 3; i__ <= i__1; ++i__) {
	it = ih[i__ - 1];
	ilp = i__ - 2;
	il[it - 1] = ilp;
	ih[ilp - 1] = it;
/* L1: */
    }
/*                                                                       ABBY1912 */
    return 0;
} /* lolpop_ */

#undef namsub
#undef jdiag
#undef arr
#undef il
#undef ih
#undef npoint
#undef nbnode
#undef ilast
#undef m
#undef j6c
#undef j9c
#undef j6
#undef j9


#undef jdiag_ref
#undef arr_ref


/*                                                                       ABBY1915 */
/*                                                                       ABBY1916 */
/*     ***************************                                       ABBY1917 */
/* Subroutine */ int neibor_(integer *lc, integer *l1, integer *l2)
{
/*     ***************************                                       ABBY1919 */
/*                                                                       ABBY1920 */
/*  ***GIVES THE POSITIONS OF THE OTHER TWO ARGUMENTS IN THE TRIAD.      ABBY1921 */
/*                                                                       ABBY1922 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY1923 */
/*                                                                       ABBY1924 */
    if (*lc < 2) {
	*l1 = 2;
	*l2 = 3;
    } else {
	if (*lc == 2) {
	    *l1 = 3;
	    *l2 = 1;
	} else {
	    *l1 = 1;
	    *l2 = 2;
	}
    }
    return 0;
} /* neibor_ */
#define nbtr (tree_1.nbtr)
#define jdiag (graph_1.jdiag)
#define arr (graph_1.arr)
#define il (graph_1.il)
#define ih (graph_1.ih)
#define nbnode (graph_1.nbnode)
#define ilast (graph_1.ilast)
#define nfree (graph_1.nfree)
#define itfree (graph_1.itfree)
#define nfin (graph_1.nfin)
#define mp (argu_1.mp)
#define jkp (keep_1.jkp)
#define jarr (keep_1.jarr)
#define it2 (keep_1.it2)
#define it3 (keep_1.it3)
#define it5 (keep_1.it5)
#define ial (build_1.ial)
#define if1 (build_1.if1)
#define if2 (build_1.if2)
#define node (build_1.node)

/*                                                                       ABBY1939 */
/*                                                                       ABBY1940 */
/*     ****************                                                  ABBY1941 */
/* Subroutine */ int order_(void)
{
    /* Initialized data */

    static char name__[6] = "ORDER ";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, j, k, i1, i2, n0, nb, nf, ik, nm, jt, it, nbt, nft, nfr, isw,
	     nbt1, nbtt;
    extern /* Subroutine */ int change_(integer *, integer *), printj_(char *,
	     integer *, ftnlen);


#define arr_ref(a_1,a_2) arr[(a_2)*48 + a_1 - 49]
#define jkp_ref(a_1,a_2) jkp[(a_2)*2 + a_1 - 3]
#define jarr_ref(a_1,a_2) jarr[(a_2)*2 + a_1 - 3]
#define jdiag_ref(a_1,a_2) jdiag[(a_2)*48 + a_1 - 49]

/*     ****************                                                  ABBY1943 */
/*                                                                       ABBY1944 */
/*  ***THIS SUBROUTINE ORDERS THE TRIADS WHICH WERE LEFT WITH FREE ENDS  ABBY1945 */
/*  ***AS CONSEQUENCE OF CUTTING,SO THAT THE NEW GRAPH WILL START THERE. ABBY1946 */
/*                                                                       ABBY1947 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY1948 */
/*                                                                       ABBY1949 */
/*                                                                       ABBY1952 */
/*                                                                       ABBY1954 */
/*                                                                       ABBY1956 */
/*                                                                       ABBY1958 */
/*                                                                       ABBY1961 */
/*                                                                       ABBY1964 */
/*                                                                       ABBY1968 */
/*                                                                       ABBY1971 */
/*                                                                       ABBY1973 */
/*                                                                       ABBY1974 */
    i__1 = mp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ial[i__ - 1] = 0;
/* L10: */
    }
/*                                                                       ABBY1978 */
    if (nfin != 0) {
	nbt1 = nbtr - 1;
	nbt = nbt1 + nfin;
	nbtt = nbt + 1;
	nb = 0;
	goto L31;
    }
/*                                                                       ABBY1986 */
    nf = nbtr - itfree[0];
/*                                                                       ABBY1988 */
    if (it5 == 0) {
	nbt1 = nbtr - 1;
	n0 = 0;
	nft = nfree;
	isw = 2;
	goto L100;
    }
/*                                                                       ABBY1996 */
    nft = it5 - it2;
    nm = nft + nbtr + 1;
    nbt1 = nbtr;
/*                                                                       ABBY2000 */
    for (j = 1; j <= 3; ++j) {
	jdiag_ref(nbtr, j) = jkp_ref(1, j);
	arr_ref(nbtr, j) = jarr_ref(1, j);
	jdiag_ref(nm, j) = jkp_ref(2, j);
	arr_ref(nm, j) = jarr_ref(2, j);
/* L21: */
    }
/*                                                                       ABBY2007 */
    jt = jdiag_ref(nm, 1);
    n0 = 0;
    isw = 1;
    goto L100;
/*                                                                       ABBY2012 */
L22:
    n0 = nft;
    nbt1 += n0;
    nft = it3 - it5;
    isw = 3;
    goto L100;
/*                                                                       ABBY2018 */
L24:
    ++nft;
/*                                                                       ABBY2020 */
L23:
    node = nbt1 + nft;
    change_(&node, &c__2);
    goto L40;
/*                                                                       ABBY2024 */
/*                                                                       ABBY2025 */
L31:
    i__1 = nbnode;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = ih[i__ - 1];
	if (il[i1 - 1] > ilast) {
	    goto L35;
	}
	i2 = nbt1 + i__;
	if (i1 > nbtt) {
	    goto L33;
	}
	if (i1 == i2) {
	    goto L32;
	}
	if (il[i2 - 1] <= nbnode) {
	    goto L35;
	}
/*                                                                       ABBY2033 */
L33:
	for (j = 1; j <= 3; ++j) {
	    jdiag_ref(i2, j) = jdiag_ref(i1, j);
	    arr_ref(i2, j) = arr_ref(i1, j);
/* L34: */
	}
/*                                                                       ABBY2038 */
	il[i1 - 1] = ilast + i__;
L32:
	++nb;
	il[i2 - 1] = 0;
/*                                                                       ABBY2042 */
L35:
	;
    }
/*                                                                       ABBY2044 */
    if (nb != nfin) {
	goto L31;
    }
    node = nbt;
L40:
    if1 = jdiag_ref(nbtr, 1);
    if2 = jdiag_ref(nbtr, 3);
/*                                                                       ABBY2049 */
    i__1 = node;
    for (i__ = nbtr; i__ <= i__1; ++i__) {
	for (k = 1; k <= 3; ++k) {
	    j = jdiag_ref(i__, k);
	    ++ial[j - 1];
/* L50: */
	}
/* L51: */
    }
/*                                                                       ABBY2056 */
    ilast = node;
    printj_(name__, &c__8, (ftnlen)6);
/*                                                                       ABBY2059 */
    return 0;
/*                                                                       ABBY2061 */
L100:
    if (nf <= 0) {
	nfr = n0;
	i1 = 1;
    } else {
	nfr = nft + 1;
	i1 = -1;
    }
/*                                                                       ABBY2069 */
    i__1 = nft;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ik = nfr + i1 * i__;
	it = itfree[ik - 1];
	k = nbt1 + ik;
/*                                                                       ABBY2074 */
	for (j = 1; j <= 3; ++j) {
	    jdiag_ref(k, j) = jdiag_ref(it, j);
	    arr_ref(k, j) = arr_ref(it, j);
/* L3: */
	}
/*                                                                       ABBY2079 */
/* L4: */
    }
/*                                                                       ABBY2081 */
    switch (isw) {
	case 1:  goto L22;
	case 2:  goto L23;
	case 3:  goto L24;
    }
    return 0;
} /* order_ */

#undef nbtr
#undef jdiag
#undef arr
#undef il
#undef ih
#undef nbnode
#undef ilast
#undef nfree
#undef itfree
#undef nfin
#undef mp
#undef jkp
#undef jarr
#undef it2
#undef it3
#undef it5
#undef ial
#undef if1
#undef if2
#undef node


#undef jdiag_ref
#undef jarr_ref
#undef jkp_ref
#undef arr_ref

#define line (tree_1.line)
#define lcol (tree_1.lcol)
#define tabs (tree_1.tabs)

/*                                                                       ABBY2084 */
/*                                                                       ABBY2085 */
/*     *********************************                                 ABBY2086 */
/* Subroutine */ int otherj_(integer *lin, integer *j, integer *lo, integer *
	lco, integer *k)
{

#define line_ref(a_1,a_2) line[(a_2)*60 + a_1 - 61]
#define lcol_ref(a_1,a_2) lcol[(a_2)*60 + a_1 - 61]

/*     *********************************                                 ABBY2088 */
/*                                                                       ABBY2089 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY2090 */
/*                                                                       ABBY2091 */
/*  ***GIVES THE OTHER TRIAD WHERE A GIVEN J OCCURS AND ITS POSITION.    ABBY2092 */
/*                                                                       ABBY2093 */
/*                                                                       ABBY2094 */
/*                                                                       ABBY2097 */
/*                                                                       ABBY2102 */
/*                                                                       ABBY2103 */
    *lo = line_ref(*j, 1);
    if (*lo == *lin || tabs[*lo - 1]) {
	*k = 1;
	*lo = line_ref(*j, 2);
	*lco = lcol_ref(*j, 2);
    } else {
	*k = 2;
	*lco = lcol_ref(*j, 1);
    }
/*                                                                       ABBY2113 */
    return 0;
} /* otherj_ */

#undef line
#undef lcol
#undef tabs


#undef lcol_ref
#undef line_ref

#define j7c (argu_1.j7c)
#define j7 (argu_1.j7)

/*                                                                       ABBY2116 */
/*                                                                       ABBY2117 */
/*     ***************************                                       ABBY2118 */
/* Subroutine */ int phase_(integer *l, integer *jm, integer *ndim)
{
    /* System generated locals */
    integer jm_dim1, jm_offset;


#define jm_ref(a_1,a_2) jm[(a_2)*jm_dim1 + a_1]

/*     ***************************                                       ABBY2120 */
/*                                                                       ABBY2121 */
/*  ***PHASE FACTOR ARISING FROM NON-CYCLIC PERMUTATION OF ARGUMENTS IN  ABBY2122 */
/*  ***TRIAD L.JM MAY BE EITHER J23 OR JDIAG.                            ABBY2123 */
/*                                                                       ABBY2124 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY2125 */
/*                                                                       ABBY2126 */
/*                                                                       ABBY2129 */
/*                                                                       ABBY2134 */
/*                                                                       ABBY2135 */
    /* Parameter adjustments */
    jm_dim1 = *ndim;
    jm_offset = 1 + jm_dim1;
    jm -= jm_offset;

    /* Function Body */
    j7[j7c] = jm_ref(*l, 1);
    j7[j7c + 1] = jm_ref(*l, 2);
    j7c += 3;
    j7[j7c - 1] = jm_ref(*l, 3);
/*                                                                       ABBY2140 */
    return 0;
} /* phase_ */

#undef j7c
#undef j7


#undef jm_ref

#define j8c (argu_1.j8c)
#define j8 (argu_1.j8)

/*                                                                       ABBY2143 */
/*                                                                       ABBY2144 */
/*     ********************                                              ABBY2145 */
/* Subroutine */ int phase2_(integer *j)
{
/*     ********************                                              ABBY2147 */
/*                                                                       ABBY2148 */
/*  ***ADDS A PHASE FACTOR (-1)**2J                                      ABBY2149 */
/*                                                                       ABBY2150 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY2151 */
/*                                                                       ABBY2152 */
/*                                                                       ABBY2155 */
/*                                                                       ABBY2159 */
/*                                                                       ABBY2160 */
    ++j8c;
    j8[j8c - 1] = *j;
/*                                                                       ABBY2163 */
    return 0;
} /* phase2_ */

#undef j8c
#undef j8

#define jdiag (graph_1.jdiag)
#define arr (graph_1.arr)
#define tab1 (graph_1.tab1)
#define npoint (graph_1.npoint)
#define iparts (graph_1.iparts)
#define ipartl (graph_1.ipartl)
#define npart (graph_1.npart)
#define nc (graph_1.nc)
#define j6c (argu_1.j6c)
#define jwc (argu_1.jwc)
#define j6 (argu_1.j6)
#define kw (argu_1.kw)
#define sumvar (argu_1.sumvar)
#define mp (argu_1.mp)

/*                                                                       ABBY2166 */
/*                                                                       ABBY2167 */
/*     ***********************                                           ABBY2168 */
/* Subroutine */ int polygn_(integer *jpol)
{
    /* Initialized data */

    static char name__[6] = "POLYGN";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, jb, jc, je, nc1, nc2, it1, it2, nbc, jar;
    extern /* Subroutine */ int phase2_(integer *), printj_(char *, integer *,
	     ftnlen);


#define kw_ref(a_1,a_2) kw[(a_2)*6 + a_1 - 7]
#define arr_ref(a_1,a_2) arr[(a_2)*48 + a_1 - 49]
#define tab1_ref(a_1,a_2) tab1[(a_2)*60 + a_1 - 61]
#define jdiag_ref(a_1,a_2) jdiag[(a_2)*48 + a_1 - 49]

/*     ***********************                                           ABBY2170 */
/*                                                                       ABBY2171 */
/*  ***THIS ROUTINE REDUCES A CIRCUIT OF ARBITRARY ORDER NC.IT EXCHANGES ABBY2172 */
/*  ***NODES ON THE FLAT DIAGRAM UNTIL THE DISTANCE ON THE AXIS BETWEEN  ABBY2173 */
/*  ***NODES EQUEALS ONE.EACH EXCHANGE INTRODUCES A SUMMATION VARIABLE   ABBY2174 */
/*  ***AND A 6J SYMBOL.THE CIRCUIT HAS A MAXIMUM OF NPART=2 DISCONNECTED ABBY2175 */
/*  ***PARTS ON THE AXIS.                                                ABBY2176 */
/*                                                                       ABBY2177 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY2178 */
/*                                                                       ABBY2179 */
/*                                                                       ABBY2182 */
/*                                                                       ABBY2186 */
/*                                                                       ABBY2190 */
/*                                                                       ABBY2193 */
/*                                                                       ABBY2195 */
/*                                                                       ABBY2196 */
    nc1 = nc + 1;
    nc2 = nc;
    nbc = ipartl - 2;
/*                                                                       ABBY2200 */
L10:
    i__1 = nbc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	it2 = npoint[nc1 - i__ - 1];
	it1 = npoint[nc2 - i__ - 1];
	jb = jdiag_ref(it1, 1);
	jc = jdiag_ref(it2, 1);
	jdiag_ref(it1, 1) = jc;
	jdiag_ref(it2, 1) = jb;
	jar = arr_ref(it1, 1);
	arr_ref(it1, 1) = arr_ref(it2, 1);
	arr_ref(it2, 1) = jar;
	je = jdiag_ref(it1, 2);
	++mp;
	sumvar[mp - 1] = TRUE_;
	jdiag_ref(it1, 2) = mp;
	jdiag_ref(it2, 3) = mp;
/*                                                                       ABBY2216 */
	if (tab1_ref(jb, 1) == it1) {
	    tab1_ref(jb, 1) = it2;
	} else {
	    tab1_ref(jb, 2) = it2;
	}
/*                                                                       ABBY2222 */
	if (tab1_ref(jc, 1) == it2) {
	    tab1_ref(jc, 1) = it1;
	} else {
	    tab1_ref(jc, 2) = it1;
	}
/*                                                                       ABBY2228 */
	if (arr_ref(it1, 2) <= 0) {
	    phase2_(&je);
	    arr_ref(it1, 2) = 1;
	    arr_ref(it2, 3) = -1;
	}
/*                                                                       ABBY2234 */
	++jwc;
	kw_ref(1, jwc) = jb;
	kw_ref(2, jwc) = mp;
	kw_ref(3, jwc) = je;
	kw_ref(4, jwc) = jc;
	kw_ref(5, jwc) = jdiag_ref(it2, 2);
	kw_ref(6, jwc) = jdiag_ref(it1, 3);
	j6[j6c] = mp;
	j6c += 2;
	j6[j6c - 1] = mp;
/* L8: */
    }
/*                                                                       ABBY2246 */
    nc -= nbc;
/*                                                                       ABBY2248 */
    if (nc > 4) {
	nbc = iparts - 2;
	nc1 = iparts + 1;
	nc2 = iparts;
	goto L10;
    }
/*                                                                       ABBY2255 */
    if (npart != 1) {
	npoint[2] = npoint[nc1 - 1];
	npoint[3] = npoint[nc1];
    }
/*                                                                       ABBY2260 */
    if (nc == 2) {
	*jpol = 1;
    }
    printj_(name__, &c__10, (ftnlen)6);
/*                                                                       ABBY2263 */
    return 0;
} /* polygn_ */

#undef jdiag
#undef arr
#undef tab1
#undef npoint
#undef iparts
#undef ipartl
#undef npart
#undef nc
#undef j6c
#undef jwc
#undef j6
#undef kw
#undef sumvar
#undef mp


#undef jdiag_ref
#undef tab1_ref
#undef arr_ref
#undef kw_ref

#define j23 (tree_1.j23)
#define arrow (tree_1.arrow)
#define line (tree_1.line)
#define lcol (tree_1.lcol)
#define tabs (tree_1.tabs)
#define nbtr (tree_1.nbtr)
#define jdiag (graph_1.jdiag)
#define arr (graph_1.arr)
#define tab1 (graph_1.tab1)
#define il (graph_1.il)
#define ih (graph_1.ih)
#define npoint (graph_1.npoint)
#define nbnode (graph_1.nbnode)
#define ifirst (graph_1.ifirst)
#define ilast (graph_1.ilast)
#define iparts (graph_1.iparts)
#define ipartl (graph_1.ipartl)
#define npart (graph_1.npart)
#define icross (graph_1.icross)
#define nfree (graph_1.nfree)
#define nfin (graph_1.nfin)
#define nc (graph_1.nc)
#define m (couple_1.m)
#define n (couple_1.n)
#define j1 (couple_1.j1)
#define j6c (argu_1.j6c)
#define j7c (argu_1.j7c)
#define j8c (argu_1.j8c)
#define j9c (argu_1.j9c)
#define jwc (argu_1.jwc)
#define j6 (argu_1.j6)
#define j7 (argu_1.j7)
#define j8 (argu_1.j8)
#define j9 (argu_1.j9)
#define kw (argu_1.kw)
#define jdel (argu_1.jdel)
#define ldel (argu_1.ldel)
#define sumvar (argu_1.sumvar)
#define mp (argu_1.mp)
#define i6c (const_1.i6c)
#define i7c (const_1.i7c)
#define i8c (const_1.i8c)
#define i9c (const_1.i9c)
#define idel (const_1.idel)
#define iwc (const_1.iwc)
#define ibug3 (debug_1.ibug3)
#define nzero (zer_1.nzero)
#define jzero (zer_1.jzero)

/*                                                                       ABBY2266 */
/*                                                                       ABBY2267 */
/*     ***************************                                       ABBY2268 */
/* Subroutine */ int printj_(char *names, integer *jp, ftnlen names_len)
{
    /* Initialized data */

    static char iblank[8] = "        ";
    static char ifree[8] = "FREE END";
    static char ip[1] = "+";
    static char im[1] = "-";
    static char nsettb[6] = "SETTAB";
    static char i6[4] = "I6= ";
    static char i7[4] = "I7= ";
    static char i8[4] = "I8= ";
    static char i9[4] = "I9= ";
    static char ij1[4] = "J1= ";

    /* Format strings */
    static char fmt_1000[] = "(/10x,\002NBNODE=\002,i3,10x,\002NBTR=\002,i3,"
	    "10x,\002NFIN=\002,i3,/10x,\002IFIRST=\002,i3,10x,\002ILAST=\002,"
	    "i3,9x,\002NFREE=\002,i3)";
    static char fmt_1001[] = "(//7x,\002IL\002,3x,\002IH\002,14x,\002JDIA"
	    "G\002//)";
    static char fmt_1002[] = "(28x,3(a1,2x))";
    static char fmt_1003[] = "(7x,i2,3x,i2,2x,a8,2x,3i3/)";
    static char fmt_1004[] = "(/5x,\002TAB1\002/)";
    static char fmt_1005[] = "(4(i3,\002)\002,2x,i3,i5,5x))";
    static char fmt_1006[] = "(/2x,\002SUMVAR=\002,15(i3,l1))";
    static char fmt_1010[] = "(//10x,\002J23\002,10x,\002NBTR1=\002,i3//)";
    static char fmt_1012[] = "(18x,3(a1,2x))";
    static char fmt_1013[] = "(i9,i5,2x,3i3/)";
    static char fmt_1014[] = "(/3x,\002J  L1 K1  L2 K2\002)";
    static char fmt_1015[] = "(4(i4,\002)\002,i3,i3,i4,i3))";
    static char fmt_1020[] = "(/3x,a4,3x,3(20i3/))";
    static char fmt_1021[] = "(/3x,\002DELTA=\002,7(i5,i3))";
    static char fmt_1022[] = "(/3x,\002KW(ARG. OF 6J)\002,6i3)";
    static char fmt_1030[] = "(//2x,\002NC=\002,i2,4x,\002NPART=\002,i2,4x"
	    ",\002IPARTL=\002,i2,4x,\002IPARTS=\002,i2,4x,\002ICROSS=\002,i2,"
	    "4x,/2x,\002NPOINT=\002,20i3)";
    static char fmt_1040[] = "(//2x,\002NZERO=\002,i2,5x,12(i4,\002)\002,i3))"
	    ;
    static char fmt_1050[] = "(///3x,\002PRINT OUT AFTER CALLING SUBROUTINE"
	    " \002,a7)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    integer i__, j, k, mm;
    char is[1*3];
    integer it, jt;
#define ix ((integer *)&const_1)
    char ifr[8];
    integer jtab[90]	/* was [30][3] */, jump, nbtr1, ntime;

    /* Fortran I/O blocks */
    static cilist io___270 = { 0, 6, 0, fmt_1050, 0 };
    static cilist io___273 = { 0, 6, 0, fmt_1020, 0 };
    static cilist io___274 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___275 = { 0, 6, 0, fmt_1001, 0 };
    static cilist io___283 = { 0, 6, 0, fmt_1002, 0 };
    static cilist io___284 = { 0, 6, 0, fmt_1003, 0 };
    static cilist io___285 = { 0, 6, 0, fmt_1004, 0 };
    static cilist io___287 = { 0, 6, 0, fmt_1005, 0 };
    static cilist io___288 = { 0, 6, 0, fmt_1006, 0 };
    static cilist io___290 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___291 = { 0, 6, 0, fmt_1012, 0 };
    static cilist io___292 = { 0, 6, 0, fmt_1013, 0 };
    static cilist io___293 = { 0, 6, 0, fmt_1014, 0 };
    static cilist io___295 = { 0, 6, 0, fmt_1015, 0 };
    static cilist io___296 = { 0, 6, 0, fmt_1030, 0 };
    static cilist io___297 = { 0, 6, 0, fmt_1040, 0 };
    static cilist io___298 = { 0, 6, 0, fmt_1020, 0 };
    static cilist io___299 = { 0, 6, 0, fmt_1020, 0 };
    static cilist io___300 = { 0, 6, 0, fmt_1020, 0 };
    static cilist io___301 = { 0, 6, 0, fmt_1020, 0 };
    static cilist io___302 = { 0, 6, 0, fmt_1021, 0 };
    static cilist io___303 = { 0, 6, 0, fmt_1022, 0 };



#define j23_ref(a_1,a_2) j23[(a_2)*24 + a_1 - 25]
#define kw_ref(a_1,a_2) kw[(a_2)*6 + a_1 - 7]
#define arr_ref(a_1,a_2) arr[(a_2)*48 + a_1 - 49]
#define tab1_ref(a_1,a_2) tab1[(a_2)*60 + a_1 - 61]
#define jtab_ref(a_1,a_2) jtab[(a_2)*30 + a_1 - 31]
#define ldel_ref(a_1,a_2) ldel[(a_2)*20 + a_1 - 21]
#define line_ref(a_1,a_2) line[(a_2)*60 + a_1 - 61]
#define lcol_ref(a_1,a_2) lcol[(a_2)*60 + a_1 - 61]
#define jdiag_ref(a_1,a_2) jdiag[(a_2)*48 + a_1 - 49]
#define arrow_ref(a_1,a_2) arrow[(a_2)*24 + a_1 - 25]

/*     ***************************                                       ABBY2270 */
/*                                                                       ABBY2271 */
/*  ***THIS SUBROUTINE PRINTS INTERMEDIATE RESULTS IN STANDARD FORM FROM ABBY2272 */
/*  ***WHEREVER IT IS CALLED.                                            ABBY2273 */
/*                                                                       ABBY2274 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY2275 */
/*                                                                       ABBY2276 */
/*                                                                       ABBY2280 */
/*                                                                       ABBY2282 */
/*                                                                       ABBY2284 */
/*                                                                       ABBY2289 */
/*                                                                       ABBY2291 */
/*                                                                       ABBY2293 */
/*                                                                       ABBY2296 */
/*                                                                       ABBY2299 */
/*                                                                       ABBY2304 */
/*                                                                       ABBY2308 */
/*                                                                       ABBY2312 */
/* L1000: */
/* L1001: */
/* L1002: */
/* L1003: */
/* L1004: */
/* L1005: */
/* L1006: */
/* L1010: */
/* L1012: */
/* L1013: */
/* L1014: */
/* L1015: */
/* L1020: */
/* L1021: */
/* L1022: */
/* L1030: */
/* L1040: */
/* L1050: */
/*                                                                       ABBY2333 */
/*                                                                       ABBY2334 */
/*                                                                       ABBY2335 */
    if (ibug3 != 1) {
	return 0;
    }
    s_wsfe(&io___270);
    do_fio(&c__1, names, (ftnlen)6);
    e_wsfe();
/*                                                                       ABBY2338 */
/*                                                                       ABBY2339 */
    jump = *jp;
    if (jump == 0) {
/*                                                                       ABBY2342 */
	for (i__ = 1; i__ <= 7; ++i__) {
	    ix[i__ - 1] = 1;
/* L9: */
	}
/*                                                                       ABBY2346 */
	s_wsfe(&io___273);
	do_fio(&c__1, ij1, (ftnlen)4);
	i__1 = m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&j1[i__ - 1], (ftnlen)sizeof(integer));
	}
	e_wsfe();
    }
/*                                                                       ABBY2349 */
    if (jump < 8) {
	goto L20;
    }
    s_wsfe(&io___274);
    do_fio(&c__1, (char *)&nbnode, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nbtr, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nfin, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ifirst, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ilast, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nfree, (ftnlen)sizeof(integer));
    e_wsfe();
    jump += -8;
    s_wsfe(&io___275);
    e_wsfe();
    k = 0;
/*                                                                       ABBY2355 */
    i__1 = nbnode;
    for (i__ = 1; i__ <= i__1; ++i__) {
	it = ih[i__ - 1];
	s_copy(ifr, iblank, (ftnlen)8, (ftnlen)8);
	jt = jdiag_ref(it, 1);
/*                                                                       ABBY2360 */
	if (tab1_ref(jt, 2) != it || jt == jdiag_ref(ifirst, 3)) {
	    ++k;
	    jtab_ref(k, 1) = jt;
	    jtab_ref(k, 2) = tab1_ref(jt, 1);
	    jtab_ref(k, 3) = tab1_ref(jt, 2);
	}
/*                                                                       ABBY2367 */
	if (tab1_ref(jt, 2) > ilast) {
	    s_copy(ifr, ifree, (ftnlen)8, (ftnlen)8);
	}
/*                                                                       ABBY2369 */
	for (j = 1; j <= 3; ++j) {
	    *(unsigned char *)&is[j - 1] = *(unsigned char *)&ip[0];
	    if (arr_ref(it, j) < 1) {
		*(unsigned char *)&is[j - 1] = *(unsigned char *)&im[0];
	    }
/* L2: */
	}
/*                                                                       ABBY2374 */
	s_wsfe(&io___283);
	for (j = 1; j <= 3; ++j) {
	    do_fio(&c__1, is + (j - 1), (ftnlen)1);
	}
	e_wsfe();
	s_wsfe(&io___284);
	do_fio(&c__1, (char *)&il[it - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&it, (ftnlen)sizeof(integer));
	do_fio(&c__1, ifr, (ftnlen)8);
	for (j = 1; j <= 3; ++j) {
	    do_fio(&c__1, (char *)&jdiag_ref(it, j), (ftnlen)sizeof(integer));
	}
	e_wsfe();
/*                                                                       ABBY2377 */
/* L1: */
    }
/*                                                                       ABBY2379 */
    s_wsfe(&io___285);
    e_wsfe();
    ntime = 0;
    jt = jdiag_ref(ifirst, 3);
    if (jt != jdiag_ref(ilast, 2)) {
	if (tab1_ref(jt, 2) < 1000) {
	    goto L5;
	}
    }
L4:
    ++k;
    jtab_ref(k, 1) = jt;
    jtab_ref(k, 2) = tab1_ref(jt, 1);
    jtab_ref(k, 3) = tab1_ref(jt, 2);
L5:
    ++ntime;
/*                                                                       ABBY2391 */
    if (ntime != 2) {
	jt = jdiag_ref(ilast, 2);
	if (tab1_ref(jt, 2) == 1000) {
	    goto L4;
	}
    }
/*                                                                       ABBY2396 */
    s_wsfe(&io___287);
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    do_fio(&c__1, (char *)&jtab_ref(i__, j), (ftnlen)sizeof(integer));
	}
    }
    e_wsfe();
    s_wsfe(&io___288);
    i__1 = mp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&sumvar[i__ - 1], (ftnlen)sizeof(logical));
    }
    e_wsfe();
L20:
    if (jump < 4) {
	goto L30;
    }
    jump += -4;
    nbtr1 = (n << 1) - 2;
    s_wsfe(&io___290);
    do_fio(&c__1, (char *)&nbtr1, (ftnlen)sizeof(integer));
    e_wsfe();
    k = 0;
/*                                                                       ABBY2404 */
    i__1 = nbtr1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (tabs[i__ - 1]) {
	    goto L11;
	}
	++k;
/*                                                                       ABBY2408 */
	for (j = 1; j <= 3; ++j) {
	    *(unsigned char *)&is[j - 1] = *(unsigned char *)&ip[0];
	    if (arrow_ref(i__, j) < 1) {
		*(unsigned char *)&is[j - 1] = *(unsigned char *)&im[0];
	    }
/* L12: */
	}
/*                                                                       ABBY2413 */
	s_wsfe(&io___291);
	for (j = 1; j <= 3; ++j) {
	    do_fio(&c__1, is + (j - 1), (ftnlen)1);
	}
	e_wsfe();
	s_wsfe(&io___292);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	for (j = 1; j <= 3; ++j) {
	    do_fio(&c__1, (char *)&j23_ref(i__, j), (ftnlen)sizeof(integer));
	}
	e_wsfe();
/*                                                                       ABBY2416 */
L11:
	;
    }
/*                                                                       ABBY2418 */
    s_wsfe(&io___293);
    e_wsfe();
    mm = m;
    if (s_cmp(names, nsettb, (ftnlen)6, (ftnlen)6) != 0) {
	mm = m - 1;
    }
    s_wsfe(&io___295);
    i__1 = mm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	for (j = 1; j <= 2; ++j) {
	    do_fio(&c__1, (char *)&line_ref(i__, j), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&lcol_ref(i__, j), (ftnlen)sizeof(integer));
	}
    }
    e_wsfe();
/*                                                                       ABBY2423 */
L30:
    if (jump >= 2) {
	jump += -2;
	s_wsfe(&io___296);
	do_fio(&c__1, (char *)&nc, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&npart, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ipartl, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&iparts, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&icross, (ftnlen)sizeof(integer));
	i__1 = nc;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&npoint[i__ - 1], (ftnlen)sizeof(integer));
	}
	e_wsfe();
    }
/*                                                                       ABBY2428 */
    if (jump >= 1) {
	s_wsfe(&io___297);
	do_fio(&c__1, (char *)&nzero, (ftnlen)sizeof(integer));
	i__1 = nzero;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&jzero[i__ - 1], (ftnlen)sizeof(integer));
	}
	e_wsfe();
    }
    if (j6c >= i6c) {
	s_wsfe(&io___298);
	do_fio(&c__1, i6, (ftnlen)4);
	i__1 = j6c;
	for (i__ = i6c; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&j6[i__ - 1], (ftnlen)sizeof(integer));
	}
	e_wsfe();
    }
    if (j7c >= i7c) {
	s_wsfe(&io___299);
	do_fio(&c__1, i7, (ftnlen)4);
	i__1 = j7c;
	for (i__ = i7c; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&j7[i__ - 1], (ftnlen)sizeof(integer));
	}
	e_wsfe();
    }
    if (j8c >= i8c) {
	s_wsfe(&io___300);
	do_fio(&c__1, i8, (ftnlen)4);
	i__1 = j8c;
	for (i__ = i8c; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&j8[i__ - 1], (ftnlen)sizeof(integer));
	}
	e_wsfe();
    }
    if (j9c >= i9c) {
	s_wsfe(&io___301);
	do_fio(&c__1, i9, (ftnlen)4);
	i__1 = j9c;
	for (i__ = i9c; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&j9[i__ - 1], (ftnlen)sizeof(integer));
	}
	e_wsfe();
    }
    if (jdel >= idel) {
	s_wsfe(&io___302);
	i__1 = jdel;
	for (i__ = idel; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 2; ++j) {
		do_fio(&c__1, (char *)&ldel_ref(i__, j), (ftnlen)sizeof(
			integer));
	    }
	}
	e_wsfe();
    }
    if (jwc >= iwc) {
	s_wsfe(&io___303);
	i__1 = jwc;
	for (i__ = iwc; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 6; ++j) {
		do_fio(&c__1, (char *)&kw_ref(j, i__), (ftnlen)sizeof(integer)
			);
	    }
	}
	e_wsfe();
    }
    i6c = j6c + 1;
    i7c = j7c + 1;
    i8c = j8c + 1;
    i9c = j9c + 1;
    idel = jdel + 1;
    iwc = jwc + 1;
    return 0;
} /* printj_ */

#undef j23
#undef arrow
#undef line
#undef lcol
#undef tabs
#undef nbtr
#undef jdiag
#undef arr
#undef tab1
#undef il
#undef ih
#undef npoint
#undef nbnode
#undef ifirst
#undef ilast
#undef iparts
#undef ipartl
#undef npart
#undef icross
#undef nfree
#undef nfin
#undef nc
#undef m
#undef n
#undef j1
#undef j6c
#undef j7c
#undef j8c
#undef j9c
#undef jwc
#undef j6
#undef j7
#undef j8
#undef j9
#undef kw
#undef jdel
#undef ldel
#undef sumvar
#undef mp
#undef i6c
#undef i7c
#undef i8c
#undef i9c
#undef idel
#undef iwc
#undef ibug3
#undef nzero
#undef jzero


#undef arrow_ref
#undef jdiag_ref
#undef lcol_ref
#undef line_ref
#undef ldel_ref
#undef jtab_ref
#undef tab1_ref
#undef arr_ref
#undef kw_ref
#undef j23_ref
#undef ix

#define jdiag (graph_1.jdiag)
#define tab1 (graph_1.tab1)
#define il (graph_1.il)
#define ih (graph_1.ih)
#define npoint (graph_1.npoint)
#define nbnode (graph_1.nbnode)
#define ifirst (graph_1.ifirst)
#define ilast (graph_1.ilast)
#define iparts (graph_1.iparts)
#define ipartl (graph_1.ipartl)
#define npart (graph_1.npart)
#define icross (graph_1.icross)
#define nc (graph_1.nc)

/*                                                                       ABBY2444 */
/*                                                                       ABBY2445 */
/*     ***********************                                           ABBY2446 */
/* Subroutine */ int search_(logical *find)
{
    /* Initialized data */

    static char name__[6] = "SEARCH";

    /* Format strings */
    static char fmt_1000[] = "(\002 ERROR IN SEARCH.I,I1,I2,I3,NPART,IPART,N"
	    "C=\002,7i5)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer i_sign(integer *, integer *), s_wsfe(cilist *), do_fio(integer *, 
	    char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer i__, i1, i2, i3, k2, k3, k4, i4, ja, ic, jc, i20, i21, i31, ii, 
	    jb, jd, in, i30, it, jt, nc2, ii2, ii1, ncm, nbn, ips, jps, ncm1;
    extern /* Subroutine */ int phase_(integer *, integer *, integer *);
    integer idist, ipart;
    extern /* Subroutine */ int change_(integer *, integer *), printj_(char *,
	     integer *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___335 = { 0, 6, 0, fmt_1000, 0 };



#define tab1_ref(a_1,a_2) tab1[(a_2)*60 + a_1 - 61]
#define jdiag_ref(a_1,a_2) jdiag[(a_2)*48 + a_1 - 49]

/*     ***********************                                           ABBY2448 */
/*                                                                       ABBY2449 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY2450 */
/*                                                                       ABBY2451 */
/*                                                                       ABBY2454 */
/*                                                                       ABBY2458 */
/*                                                                       ABBY2462 */
/*                                                                       ABBY2464 */
/* L1000: */
/*                                                                       ABBY2466 */
/*  ***THIS SUBROUTINE LOCATES CIRCUITS OR LOOPS OF ORDER NC.NPOINT(NC)  ABBY2467 */
/*  ***ARE THE INDICES OF THE POINTS(TRIADS) PERTAINING TO THE FIRST     ABBY2468 */
/*  ***SUCH LOOP FOUND.                                                  ABBY2469 */
/*  ***NPART IS THE NUMBER OF SEPARATE PARTS(GROUPS OF CONTIGUOUS POINTS)ABBY2470 */
/*  ***ON THE AXIS OF THE FLAT GRAPH.IPARTS IS THE NUMBER OF POINTS IN   ABBY2471 */
/*  ***THE SMALLEST PART.IPARTL IS THE NUMBER OF POINTS IN THE LARGEST   ABBY2472 */
/*  ***PART.                                                             ABBY2473 */
/*  ***THIS SUBROUTINE FINDS ALL THE POSSIBLE LOOPS OF ORDER 3 AND 4.FOR ABBY2474 */
/*  ***NC.GE.5,IT LOOKS FOR ONLY THOSE WHO ARE PARTITIONNED IN NPART.LE.2ABBY2475 */
/*  ***WHICH CAN EVENTUALLY REDUCE TP A LOOP OF ORDER 4 WITHOUT BREAKING ABBY2476 */
/*  ***THE BASIC STRUCTURE OF THE FLAT GRAPH. ICROSS=-1,IF LINES CROSS   ABBY2477 */
/* --------------------------------------------------------------------   ABBY2478 */
/*                                                                       ABBY2479 */
/*                                                                       ABBY2480 */
/*  ***INITIALIZATION                                                    ABBY2481 */
/*                                                                       ABBY2482 */
    *find = FALSE_;
    ncm1 = nc - 1;
    ncm = nc - 2;
    icross = 0;
/*                                                                       ABBY2487 */
/*  ***FIRST ARE TREATED TWO CASES THAT DO NOT INVOLVE DO LOOPS          ABBY2488 */
/*  ***1.ONE ISOLATED POINT,EITHER THE FIRST OR THE LAST                 ABBY2489 */
/*                                                                       ABBY2490 */
    npart = 1;
    ipartl = nc - 1;
    iparts = 1;
/*                                                                       ABBY2494 */
/*  ***A.FIRST                                                           ABBY2495 */
/*                                                                       ABBY2496 */
    i1 = ifirst;
    k3 = 3;
    k2 = 2;
L200:
    ja = jdiag_ref(i1, 1);
    jc = jdiag_ref(i1, k3);
/*                                                                       ABBY2502 */
    if (ja == jc) {
	if (nc > 1) {
	    goto L800;
	}
	npoint[0] = i1;
	goto L900;
    }
/*                                                                       ABBY2508 */
    i2 = tab1_ref(ja, k2);
    i3 = tab1_ref(jc, k2);
/*                                                                       ABBY2511 */
    if ((i__1 = il[i3 - 1] - il[i2 - 1], abs(i__1)) - ncm < 0) {
	goto L800;
    }
/*                                                                       ABBY2513 */
    if ((i__1 = il[i3 - 1] - il[i2 - 1], abs(i__1)) - ncm > 0) {
/*                                                                       ABBY2515 */
/*  ***B.LAST                                                            ABBY2516 */
/*                                                                       ABBY2517 */
	if (i1 != ifirst) {
	    goto L250;
	}
	i1 = ilast;
	k3 = 2;
	k2 = 1;
	goto L200;
    }
/*                                                                       ABBY2524 */
    ic = 1;
    npoint[ic - 1] = i1;
    i20 = min(i2,i3);
    i21 = il[i20 - 1];
    i31 = i21 + ncm1;
/*                                                                       ABBY2530 */
    i__1 = i31;
    for (ii = i21; ii <= i__1; ++ii) {
	++ic;
	npoint[ic - 1] = ih[ii - 1];
/* L203: */
    }
/*                                                                       ABBY2535 */
    if (nc <= 2) {
	if (jdiag_ref(ifirst, 1) != jdiag_ref(ilast, 1)) {
	    phase_(&i1, jdiag, &c__48);
	}
	goto L900;
    }
/*                                                                       ABBY2540 */
    if (i1 != ilast) {
	it = i2;
	jt = jdiag_ref(ilast, 2);
	k4 = 2;
	i4 = ilast;
    } else {
	it = i3;
	jt = jdiag_ref(ifirst, 3);
	k4 = 3;
	i4 = ifirst;
    }
/*                                                                       ABBY2552 */
    if (it == i20) {
	phase_(&i1, jdiag, &c__48);
    }
    if (jt == ja || jt == jc) {
	change_(&i4, &k4);
    }
    goto L900;
/*                                                                       ABBY2556 */
/*  ***2.TWO ISOLATED POINTS,FIRST AND LAST.                             ABBY2557 */
/*                                                                       ABBY2558 */
L250:
    if (nc == 1) {
	return 0;
    }
    if (nc <= 3) {
	goto L100;
    }
    ipartl = nc - 2;
    iparts = 1;
    i1 = ifirst;
    i2 = ilast;
    ja = jdiag_ref(i1, 1);
    jb = jdiag_ref(i1, 3);
/*                                                                       ABBY2567 */
    if (tab1_ref(ja, 2) != i2) {
	ja = jdiag_ref(i1, 3);
	jb = jdiag_ref(i1, 1);
	if (tab1_ref(ja, 2) != i2) {
	    goto L100;
	}
    }
/*                                                                       ABBY2573 */
    if (ja == jdiag_ref(i2, 1)) {
	jc = jdiag_ref(i2, 2);
    } else {
	jc = jdiag_ref(ilast, 1);
    }
/*                                                                       ABBY2579 */
    i3 = tab1_ref(jb, 2);
    i4 = tab1_ref(jc, 1);
    idist = il[i4 - 1] - il[i3 - 1];
/*                                                                       ABBY2583 */
    if (abs(idist) - (ncm - 1) < 0) {
	goto L800;
    }
    if (abs(idist) - (ncm - 1) == 0) {
	npoint[0] = ilast;
	npoint[1] = ifirst;
	icross = i_sign(&c__1, &idist);
	ic = 2;
	i20 = min(i3,i4);
	i21 = il[i20 - 1];
	i31 = i21 + ncm;
/*                                                                       ABBY2593 */
	i__1 = i31;
	for (ii = i21; ii <= i__1; ++ii) {
	    ++ic;
	    npoint[ic - 1] = ih[ii - 1];
/* L261: */
	}
/*                                                                       ABBY2598 */
	if (ja == jdiag_ref(ifirst, 1)) {
	    change_(&ifirst, &c__3);
	}
	if (ja == jdiag_ref(ilast, 1)) {
	    change_(&ilast, &c__2);
	}
	goto L900;
    }
/*                                                                       ABBY2603 */
/*  ***FIRST GENERAL CASE:ALL POINTS IN ONE GROUP                        ABBY2604 */
/*                                                                       ABBY2605 */
L100:
    npart = 1;
    iparts = 0;
    ipartl = nc;
    k3 = 1;
/*                                                                       ABBY2610 */
    i__1 = nbnode;
    for (in = 1; in <= i__1; ++in) {
	i__ = ih[in - 1];
L108:
	ja = jdiag_ref(i__, k3);
	if (i__ != tab1_ref(ja, 2)) {
	    i2 = tab1_ref(ja, 2);
/*                                                                       ABBY2616 */
	    if (il[i2 - 1] - in - ncm1 < 0) {
		goto L800;
	    }
	    if (il[i2 - 1] - in - ncm1 == 0) {
		i21 = il[i2 - 1];
		ic = 0;
/*                                                                       ABBY2621 */
		i__2 = i21;
		for (ii = in; ii <= i__2; ++ii) {
		    ++ic;
		    npoint[ic - 1] = ih[ii - 1];
/* L103: */
		}
/*                                                                       ABBY2626 */
		if (ja == jdiag_ref(ifirst, 3)) {
		    change_(&ifirst, &c__3);
		}
		if (ja == jdiag_ref(ilast, 2)) {
		    change_(&ilast, &c__2);
		}
		goto L900;
	    }
	}
/*                                                                       ABBY2632 */
	if (in == 1) {
	    if (k3 != 3) {
		k3 = 3;
		goto L108;
	    } else {
		k3 = 1;
	    }
	}
/*                                                                       ABBY2641 */
/* L101: */
    }
/*                                                                       ABBY2643 */
/*  ***SEARCH DID NOT FIND LOOP NC.LE.3                                  ABBY2644 */
/*                                                                       ABBY2645 */
    if (nc <= 3) {
	return 0;
    }
/*                                                                       ABBY2647 */
/*  ***GENERAL CASE OF LOOP PARTITIONNED IN 2 GROUPS.DO LOOP             ABBY2648 */
/*  ***ON IPARTS                                                         ABBY2649 */
/*                                                                       ABBY2650 */
    npart = 2;
    nc2 = nc / 2;
    k3 = 1;
    k2 = 1;
/*                                                                       ABBY2655 */
    i__1 = nc2;
    for (ips = 2; ips <= i__1; ++ips) {
	jps = ips - 1;
	nbn = nbnode - jps;
/*                                                                       ABBY2659 */
	i__2 = nbn;
	for (i1 = 1; i1 <= i__2; ++i1) {
	    i__ = ih[i1 - 1];
	    i2 = ih[i1 + jps - 1];
L302:
	    ja = jdiag_ref(i__, k3);
	    jd = jdiag_ref(i2, k2);
/*                                                                       ABBY2665 */
	    if (i__ == tab1_ref(ja, 1)) {
		ii2 = tab1_ref(jd, 2);
		ii1 = tab1_ref(ja, 2);
	    } else {
		ii1 = tab1_ref(ja, 1);
		ii2 = tab1_ref(jd, 1);
	    }
/*                                                                       ABBY2673 */
	    idist = il[ii1 - 1] - il[ii2 - 1];
/*                                                                       ABBY2675 */
	    if (abs(idist) - (ncm - jps) < 0) {
		goto L800;
	    }
	    if (abs(idist) - (ncm - jps) > 0) {
		goto L320;
	    }
/* L306: */
	    icross = i_sign(&c__1, &idist);
	    ic = 0;
	    i21 = il[i2 - 1];
/*                                                                       ABBY2681 */
	    i__3 = i21;
	    for (ii = i1; ii <= i__3; ++ii) {
		++ic;
		npoint[ic - 1] = ih[ii - 1];
/* L310: */
	    }
/*                                                                       ABBY2686 */
	    i20 = min(ii1,ii2);
	    i30 = max(ii1,ii2);
	    i21 = il[i20 - 1];
	    i31 = il[i30 - 1];
/*                                                                       ABBY2691 */
	    i__3 = i31;
	    for (ii = i21; ii <= i__3; ++ii) {
		++ic;
		npoint[ic - 1] = ih[ii - 1];
/* L311: */
	    }
/*                                                                       ABBY2696 */
	    iparts = ips;
	    ipartl = nc - ips;
	    if (jdiag_ref(ifirst, 3) == ja || jdiag_ref(ifirst, 3) == jd) {
		change_(&ifirst, &c__3);
	    }
	    if (jdiag_ref(ilast, 2) == ja || jdiag_ref(ilast, 2) == jd) {
		change_(&ilast, &c__2);
	    }
	    goto L900;
/*                                                                       ABBY2704 */
L320:
	    if (i1 == 1) {
		if (k3 == 3) {
		    k3 = 1;
		    goto L301;
		} else {
		    k3 = 3;
		    goto L302;
		}
	    }
/*                                                                       ABBY2714 */
	    if (i2 == ilast) {
		if (k2 != 2) {
		    k2 = 2;
		    goto L302;
		}
	    }
/*                                                                       ABBY2721 */
L301:
	    ;
	}
/* L400: */
    }
/*                                                                       ABBY2724 */
/*  ***SEARCH DID NOT FIND CIRCUIT OF ORDER NC                           ABBY2725 */
/*                                                                       ABBY2726 */
    return 0;
/*                                                                       ABBY2728 */
/*  ***LOOP FOUND                                                        ABBY2729 */
/*                                                                       ABBY2730 */
L900:
    *find = TRUE_;
    printj_(name__, &c__10, (ftnlen)6);
/*                                                                       ABBY2733 */
    return 0;
/*                                                                       ABBY2735 */
/*  ***ERROR PRINTOUT                                                    ABBY2736 */
/*                                                                       ABBY2737 */
L800:
    s_wsfe(&io___335);
    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&i3, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&npart, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ipart, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nc, (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
/*                                                                       ABBY2740 */
    return 0;
} /* search_ */

#undef jdiag
#undef tab1
#undef il
#undef ih
#undef npoint
#undef nbnode
#undef ifirst
#undef ilast
#undef iparts
#undef ipartl
#undef npart
#undef icross
#undef nc


#undef jdiag_ref
#undef tab1_ref

#define j6c (argu_1.j6c)
#define j7c (argu_1.j7c)
#define j8c (argu_1.j8c)
#define j9c (argu_1.j9c)
#define jwc (argu_1.jwc)
#define jdel (argu_1.jdel)
#define j6cc (dim_1.j6cc)
#define j7cc (dim_1.j7cc)
#define j8cc (dim_1.j8cc)
#define j9cc (dim_1.j9cc)
#define jwcc (dim_1.jwcc)
#define jdelc (dim_1.jdelc)

/*                                                                       ABBY2742 */
/*                                                                       ABBY2743 */
/*     *****************                                                 ABBY2744 */
/* Subroutine */ int setdim_(void)
{
/*     *****************                                                 ABBY2746 */
/*                                                                       ABBY2747 */
/*  ***SET DIMENSIONS OF ARRAYS.                                         ABBY2748 */
/*                                                                       ABBY2749 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY2750 */
/*                                                                       ABBY2751 */
/*                                                                       ABBY2754 */
/*                                                                       ABBY2756 */
/*                                                                       ABBY2759 */
/*                                                                       ABBY2761 */
/*                                                                       ABBY2762 */
/*                                                                       ABBY2763 */
    jwcc = jwc;
    jdelc = jdel;
    j6cc = j6c;
    j7cc = j7c;
    j8cc = j8c;
    j9cc = j9c;
/*                                                                       ABBY2770 */
    return 0;
} /* setdim_ */

#undef j6c
#undef j7c
#undef j8c
#undef j9c
#undef jwc
#undef jdel
#undef j6cc
#undef j7cc
#undef j8cc
#undef j9cc
#undef jwcc
#undef jdelc

#define j23 (tree_1.j23)
#define arrow (tree_1.arrow)
#define line (tree_1.line)
#define lcol (tree_1.lcol)
#define tabs (tree_1.tabs)
#define nbtr (tree_1.nbtr)
#define m (couple_1.m)
#define n (couple_1.n)
#define j2 (couple_1.j2)
#define j3 (couple_1.j3)
#define j6c (argu_1.j6c)
#define j8c (argu_1.j8c)
#define j9c (argu_1.j9c)
#define j6 (argu_1.j6)
#define j8 (argu_1.j8)
#define j9 (argu_1.j9)
#define sumvar (argu_1.sumvar)
#define ial (build_1.ial)

/*                                                                       ABBY2773 */
/*                                                                       ABBY2774 */
/*     ***********************                                           ABBY2775 */
/* Subroutine */ int settab_(logical *fail)
{
    /* Initialized data */

    static char name__[6] = "SETTAB";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, j, k, l, kc, lc, ii, ji, it, jt, nb1, jt1, ipr, ipr1;
    extern /* Subroutine */ int delta_(integer *, integer *, logical *), 
	    otherj_(integer *, integer *, integer *, integer *, integer *), 
	    printj_(char *, integer *, ftnlen);


#define j2_ref(a_1,a_2) j2[(a_2)*12 + a_1 - 13]
#define j3_ref(a_1,a_2) j3[(a_2)*12 + a_1 - 13]
#define j23_ref(a_1,a_2) j23[(a_2)*24 + a_1 - 25]
#define line_ref(a_1,a_2) line[(a_2)*60 + a_1 - 61]
#define lcol_ref(a_1,a_2) lcol[(a_2)*60 + a_1 - 61]
#define arrow_ref(a_1,a_2) arrow[(a_2)*24 + a_1 - 25]

/*     ***********************                                           ABBY2777 */
/*                                                                       ABBY2778 */
/*  ***BUILDS UP THE UNSTRUCTURED GRAPH                                  ABBY2779 */
/*  ***SETS THE ARRAY J23,CONTAINING THE TWO LISTS OF ORIGINAL TRIADS    ABBY2780 */
/*  ***J2 AND J3,AND THE CORRESPONDING ARROWS ON THE ANGULAR MOMENTA     ABBY2781 */
/*  ***LINES.ALSO ESTABLISHES THE NUMERICAL AND PHASE FACTORS CONNECTING ABBY2782 */
/*  ***RECOUPLING COEFFICIENT AND GRAPHS,ACCORDING TO YUTSIS,LEVINSON ANDABBY2783 */
/*  ***VANAGAS.FOR THIS PURPOSE DETERMINES THE TOTAL J                   ABBY2784 */
/*                                                                       ABBY2785 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY2786 */
/*                                                                       ABBY2787 */
/*                                                                       ABBY2790 */
/*                                                                       ABBY2792 */
/*                                                                       ABBY2794 */
/*                                                                       ABBY2796 */
/*                                                                       ABBY2799 */
/*                                                                       ABBY2802 */
/*                                                                       ABBY2805 */
/*                                                                       ABBY2807 */
/*                                                                       ABBY2808 */
/*                                                                       ABBY2809 */
/*                                                                       ABBY2810 */
    ipr = n - 1;
    nbtr = ipr + ipr;
/*                                                                       ABBY2813 */
    i__1 = ipr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 2; ++j) {
	    j23_ref(i__, j) = j2_ref(i__, j);
	    arrow_ref(i__, j) = 1;
/* L5: */
	}
	tabs[i__ - 1] = FALSE_;
	j23_ref(i__, 3) = j2_ref(i__, 3);
	arrow_ref(i__, 3) = -1;
/* L4: */
    }
/*                                                                       ABBY2823 */
    ipr1 = ipr + 1;
/*                                                                       ABBY2825 */
    i__1 = nbtr;
    for (i__ = ipr1; i__ <= i__1; ++i__) {
	ii = i__ - ipr;
	for (j = 1; j <= 2; ++j) {
	    j23_ref(i__, j) = j3_ref(ii, j);
	    arrow_ref(i__, j) = -1;
/* L6: */
	}
	tabs[i__ - 1] = FALSE_;
	j23_ref(i__, 3) = j3_ref(ii, 3);
	arrow_ref(i__, 3) = 1;
/* L7: */
    }
/*                                                                       ABBY2836 */
    i__1 = nbtr;
    for (j = 1; j <= i__1; ++j) {
	j8[j - 1] = j23_ref(j, 1);
/* L11: */
    }
/*                                                                       ABBY2840 */
    j8c = nbtr + ipr;
    nb1 = nbtr + 1;
/*                                                                       ABBY2843 */
    i__1 = j8c;
    for (j = nb1; j <= i__1; ++j) {
	i__ = j - ipr;
	j8[j - 1] = j23_ref(i__, 3);
/* L12: */
    }
/*                                                                       ABBY2848 */
    j6c = nbtr;
/*                                                                       ABBY2850 */
    i__1 = j6c;
    for (j = 1; j <= i__1; ++j) {
	j6[j - 1] = j23_ref(j, 3);
/* L13: */
    }
/*                                                                       ABBY2854 */
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sumvar[i__ - 1] = FALSE_;
	ial[i__ - 1] = 1;
/* L10: */
    }
/*                                                                       ABBY2859 */
    i__1 = nbtr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    ji = j23_ref(i__, j);
	    k = ial[ji - 1];
	    line_ref(ji, k) = i__;
	    lcol_ref(ji, k) = j;
	    ial[ji - 1] = k + 1;
/* L8: */
	}
/* L9: */
    }
/*                                                                       ABBY2869 */
    it = 0;
/*                                                                       ABBY2871 */
    i__1 = nbtr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	jt = j23_ref(i__, 3);
/*                                                                       ABBY2874 */
	if (ial[jt - 1] == 3) {
	    otherj_(&i__, &jt, &l, &lc, &k);
	    if (lc == 3) {
		goto L19;
	    }
	    goto L18;
	}
/*                                                                       ABBY2880 */
	if (it == 1) {
	    delta_(&jt1, &jt, fail);
	    if (*fail) {
		goto L20;
	    }
	    k = line_ref(jt, 1);
	    kc = lcol_ref(jt, 1);
	    line_ref(jt1, 2) = k;
	    lcol_ref(jt1, 2) = kc;
	    line_ref(jt, 2) = line_ref(jt1, 1);
	    lcol_ref(jt, 2) = lcol_ref(jt1, 1);
	    j23_ref(k, kc) = jt1;
	    ial[jt - 1] = 1;
	    goto L19;
	}
/*                                                                       ABBY2894 */
	jt1 = jt;
	it = 1;
/*                                                                       ABBY2897 */
L18:
	;
    }
/*                                                                       ABBY2899 */
L19:
    j9[j9c] = jt;
    j9c += 2;
    j9[j9c - 1] = jt;
/*                                                                       ABBY2903 */
L20:
    printj_(name__, &c__4, (ftnlen)6);
/*                                                                       ABBY2905 */
    return 0;
} /* settab_ */

#undef j23
#undef arrow
#undef line
#undef lcol
#undef tabs
#undef nbtr
#undef m
#undef n
#undef j2
#undef j3
#undef j6c
#undef j8c
#undef j9c
#undef j6
#undef j8
#undef j9
#undef sumvar
#undef ial


#undef arrow_ref
#undef lcol_ref
#undef line_ref
#undef j23_ref
#undef j3_ref
#undef j2_ref

#define sum6j (tree_2.sum6j)
#define t6j (tree_2.t6j)
#define jt (tree_2.jt)
#define js (tree_2.js)
#define inver (tree_2.inver)
#define jnsum (tree_2.jnsum)
#define jinv (tree_2.jinv)
#define n6jn (tree_2.n6jn)
#define in6j (tree_2.in6j)
#define jsumt (tree_2.jsumt)
#define cut (cutdig_1.cut)
#define jtem4 (graph_3.jtem4)
#define jtem5 (graph_3.jtem5)
#define jtem6 (graph_3.jtem6)
#define nsum6j (graph_3.nsum6j)
#define j6sum (graph_3.j6sum)
#define j6c (argu_2.j6c)
#define j7c (argu_2.j7c)
#define j8c (argu_2.j8c)
#define j9c (argu_2.j9c)
#define jwc (argu_2.jwc)
#define j6 (argu_2.j6)
#define j7 (argu_2.j7)
#define j8 (argu_2.j8)
#define j9 (argu_2.j9)
#define jw (argu_2.jw)
#define sumvar (argu_2.sumvar)
#define mp (argu_2.mp)
#define j6cc (dim_1.j6cc)
#define j7cc (dim_1.j7cc)
#define j8cc (dim_1.j8cc)
#define j9cc (dim_1.j9cc)
#define j6p (sumarg_1.j6p)
#define j7p (sumarg_1.j7p)
#define j8p (sumarg_1.j8p)
#define j9p (sumarg_1.j9p)
#define jword (sumarg_1.jword)
#define nlsum (sumarg_1.nlsum)
#define nbj (sumarg_1.nbj)
#define nb6j (sumarg_1.nb6j)
#define k6cp (sumarg_1.k6cp)
#define k7cp (sumarg_1.k7cp)
#define k8cp (sumarg_1.k8cp)
#define k9cp (sumarg_1.k9cp)
#define jsum6 (sumarg_1.jsum6)
#define jsum4 (sumarg_1.jsum4)
#define jsum5 (sumarg_1.jsum5)
#define inv6j (sumarg_1.inv6j)

/*                                                                       ABBY2908 */
/*                                                                       ABBY2909 */
/*     ********************                                              ABBY2910 */
/* Subroutine */ int sprate_(integer *m)
{
    /* Initialized data */

    static integer kflw = 20;
    static integer kfl6 = 180;
    static integer kfl7 = 180;
    static integer kfl8 = 180;
    static integer kfl9 = 40;
    static integer kfls = 12;

    /* Format strings */
    static char fmt_1000[] = "(2x,\002DIMENSION ERROR FOR  \002,a4,i5,\002 I"
	    "S OUT OF ALLOWED RAN   GE\002,i3)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer i__, j, k, i1, i2, i3, j1, m1, ik, jk, jj, nj, kt, k6c, k7c, k8c, 
	    k9c, i6j, j6j, n6j, js6, isk;
    extern /* Subroutine */ int var_(integer *, integer *, integer *, integer 
	    *, integer *, logical *, integer *, integer *, integer *);
    integer isu, nmx, npx, j6cp, j7cp, j8cp, j9cp;
    char name__[4];
    integer nsum;
    extern /* Subroutine */ int chvar_(integer *, integer *, integer *, 
	    logical *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___361 = { 0, 6, 0, fmt_1000, 0 };



#define jw_ref(a_1,a_2) jw[(a_2)*6 + a_1 - 7]
#define jtem4_ref(a_1,a_2) jtem4[(a_2)*12 + a_1 - 13]
#define jtem5_ref(a_1,a_2) jtem5[(a_2)*12 + a_1 - 13]
#define jsum4_ref(a_1,a_2) jsum4[(a_2)*12 + a_1 - 13]
#define jsum5_ref(a_1,a_2) jsum5[(a_2)*12 + a_1 - 13]
#define jword_ref(a_1,a_2) jword[(a_2)*6 + a_1 - 7]
#define jsumt_ref(a_1,a_2) jsumt[(a_2)*20 + a_1 - 21]

/*     ********************                                              ABBY2912 */
/*                                                                       ABBY2913 */
/*  ***THIS SUBROUTINE PREPARES THE INFORMATION TO BE TRANSFERED TO      ABBY2914 */
/*  ***GENSUM FOR NUMERICAL EVALUATION.THE COMMON BLOCKS /GRAPH/ AND     ABBY2915 */
/*  **+/TREE/ ARE USED AS WORKING MEMORY,AND THEIR PREVIOUS CONTENT      ABBY2916 */
/*  ***IS DESTROYED.                                                     ABBY2917 */
/*                                                                       ABBY2918 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY2919 */
/*                                                                       ABBY2920 */
/*                                                                       ABBY2923 */
/*                                                                       ABBY2925 */
/*                                                                       ABBY2927 */
/*                                                                       ABBY2930 */
/*                                                                       ABBY2934 */
/*                                                                       ABBY2937 */
/*                                                                       ABBY2940 */
/*                                                                       ABBY2945 */
/*                                                                       ABBY2946 */
/*                                                                       ABBY2949 */
/* L1000: */
/*                                                                       ABBY2952 */
/*                                                                       ABBY2953 */
/*                                                                       ABBY2954 */
/*  ***TEST THAT ARRAY DIMENSIONS HAVE NOT BEEN EXCEEDED.                ABBY2955 */
/*                                                                       ABBY2956 */
    if (jwc > kflw) {
	nmx = kflw;
	npx = jwc;
	s_copy(name__, "KFLW", (ftnlen)4, (ftnlen)4);
    } else {
	if (j6c > kfl6) {
	    nmx = kfl6;
	    npx = j6c;
	    s_copy(name__, "KFL6", (ftnlen)4, (ftnlen)4);
	} else {
	    if (j7c > kfl7) {
		nmx = kfl7;
		npx = j7c;
		s_copy(name__, "KFL7", (ftnlen)4, (ftnlen)4);
	    } else {
		if (j8c > kfl8) {
		    nmx = kfl8;
		    npx = j8c;
		    s_copy(name__, "KFL8", (ftnlen)4, (ftnlen)4);
		} else {
		    if (j9c <= kfl9) {
			goto L54;
		    }
		    nmx = kfl9;
		    npx = j9c;
		    s_copy(name__, "KFL9", (ftnlen)4, (ftnlen)4);
		}
	    }
	}
    }
/*                                                                       ABBY2985 */
L60:
    s_wsfe(&io___361);
    do_fio(&c__1, name__, (ftnlen)4);
    do_fio(&c__1, (char *)&npx, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nmx, (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", (ftnlen)0);
/*                                                                       ABBY2988 */
/*  ***DETERMINATION OF EFFECTIVE SUMMATION VARIABLES AND THEIR          ABBY2989 */
/*  ***RELATIONSHIPS WITH 6J COEFFICIENTS.                               ABBY2990 */
/*                                                                       ABBY2991 */
/*                                                                       ABBY2992 */
L54:
    i__1 = jwc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	inv6j[i__ - 1] = 0;
	sum6j[i__ - 1] = FALSE_;
/* L2: */
    }
/*                                                                       ABBY2997 */
    nsum = 0;
    nlsum = 0;
    if (mp == *m) {
	return 0;
    }
    m1 = *m + 1;
/*                                                                       ABBY3002 */
    i__1 = mp;
    for (i__ = m1; i__ <= i__1; ++i__) {
	if (sumvar[i__ - 1]) {
	    ++nsum;
	    jsum6[nsum - 1] = 0;
	    inver[i__ - 1] = nsum;
	}
/* L1: */
    }
/*                                                                       ABBY3010 */
    if (nsum == 0) {
	return 0;
    }
/*                                                                       ABBY3012 */
    if (nsum > kfls) {
	nmx = kfls;
	npx = nsum;
	s_copy(name__, "NSUM", (ftnlen)4, (ftnlen)4);
	goto L60;
    }
/*                                                                       ABBY3019 */
    kt = 0;
/*                                                                       ABBY3021 */
    i__1 = jwc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 6; ++j) {
	    ik = jw_ref(j, i__);
	    if (! sumvar[ik - 1]) {
		goto L5;
	    }
/*                                                                       ABBY3026 */
	    if (! sum6j[i__ - 1]) {
		sum6j[i__ - 1] = TRUE_;
		++kt;
		j6sum[kt - 1] = 0;
		nsum6j[kt - 1] = i__;
		inv6j[i__ - 1] = kt;
	    }
/*                                                                       ABBY3034 */
	    isk = inver[ik - 1];
	    i2 = jsum6[isk - 1] + 1;
	    jsum6[isk - 1] = i2;
	    jsum4_ref(isk, i2) = j;
	    jsum5_ref(isk, i2) = kt;
	    i3 = j6sum[kt - 1] + 1;
	    j6sum[kt - 1] = i3;
	    jsumt_ref(kt, i3) = isk;
L5:
	    ;
	}
/* L4: */
    }
/*                                                                       ABBY3045 */
    var_(j6, j6p, &j6c, &j6cp, &j6cc, sumvar, &mp, m, inver);
    var_(j7, j7p, &j7c, &j7cp, &j7cc, sumvar, &mp, m, inver);
    var_(j8, j8p, &j8c, &j8cp, &j8cc, sumvar, &mp, m, inver);
    var_(j9, j9p, &j9c, &j9cp, &j9cc, sumvar, &mp, m, inver);
/*                                                                       ABBY3050 */
    if (! cut) {
	nlsum = 1;
	nbj[0] = nsum;
	nb6j[0] = kt;
	k6cp[0] = j6cp;
	k7cp[0] = j7cp;
	k8cp[0] = j8cp;
	k9cp[0] = j9cp;
/*                                                                       ABBY3059 */
	i__1 = kt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i1 = nsum6j[i__ - 1];
	    for (j = 1; j <= 6; ++j) {
		jword_ref(j, i__) = jw_ref(j, i1);
/* L22: */
	    }
/* L21: */
	}
	i__1 = nsum;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    isu = jsum6[i__ - 1];
	    i__2 = isu;
	    for (j = 1; j <= i__2; ++j) {
		i1 = jsum5_ref(i__, j);
		j1 = jsum4_ref(i__, j);
		jword_ref(j1, i1) = mp + i__;
/* L81: */
	    }
/* L80: */
	}
/*                                                                       ABBY3075 */
	return 0;
    }
/*                                                                       ABBY3078 */
/*  ***SEPARATION OF VARIABLES AND SUMS IN CASE A CUT WAS DETECTED.      ABBY3079 */
/*                                                                       ABBY3080 */
    k6c = 0;
    k7c = 0;
    k8c = 0;
    k9c = 0;
    nj = 0;
    n6j = 0;
/*                                                                       ABBY3087 */
    i__1 = kt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t6j[i__ - 1] = FALSE_;
/* L9: */
    }
/*                                                                       ABBY3091 */
    i__1 = nsum;
    for (i__ = 1; i__ <= i__1; ++i__) {
	jt[i__ - 1] = FALSE_;
	js[i__ - 1] = FALSE_;
/* L7: */
    }
/*                                                                       ABBY3096 */
    j = 1;
/*                                                                       ABBY3098 */
L10:
    ++nj;
    jnsum[nj - 1] = j;
    jinv[j - 1] = nj;
    jt[j - 1] = TRUE_;
L18:
    js[j - 1] = TRUE_;
    js6 = jsum6[j - 1];
/*                                                                       ABBY3105 */
    i__1 = js6;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i6j = jsum5_ref(j, i__);
/*                                                                       ABBY3108 */
	if (! t6j[i6j - 1]) {
	    t6j[i6j - 1] = TRUE_;
	    ++n6j;
	    n6jn[n6j - 1] = nsum6j[i6j - 1];
	    in6j[i6j - 1] = n6j;
	}
/*                                                                       ABBY3115 */
	j6j = j6sum[i6j - 1];
/*                                                                       ABBY3117 */
	i__2 = j6j;
	for (k = 1; k <= i__2; ++k) {
	    jk = jsumt_ref(i6j, k);
	    if (! jt[jk - 1]) {
		++nj;
		jnsum[nj - 1] = jk;
		jinv[jk - 1] = nj;
		jt[jk - 1] = TRUE_;
	    }
/* L12: */
	}
/*                                                                       ABBY3127 */
/* L11: */
    }
/*                                                                       ABBY3129 */
    i__1 = nsum;
    for (jj = 1; jj <= i__1; ++jj) {
	j = jj;
	if (! js[jj - 1] && jt[jj - 1]) {
	    goto L18;
	}
/* L13: */
    }
/*                                                                       ABBY3134 */
    ++nlsum;
    nbj[nlsum - 1] = nj;
    nb6j[nlsum - 1] = n6j;
/*                                                                       ABBY3138 */
    if (j6cp != 0) {
	chvar_(j6p, &j6cp, &k6c, jt, jinv, &nsum);
    }
    k6cp[nlsum - 1] = k6c;
    if (j7cp != 0) {
	chvar_(j7p, &j7cp, &k7c, jt, jinv, &nsum);
    }
    k7cp[nlsum - 1] = k7c;
    if (j8cp != 0) {
	chvar_(j8p, &j8cp, &k8c, jt, jinv, &nsum);
    }
    k8cp[nlsum - 1] = k8c;
    if (j9cp != 0) {
	chvar_(j9p, &j9cp, &k9c, jt, jinv, &nsum);
    }
    k9cp[nlsum - 1] = k9c;
/*                                                                       ABBY3147 */
    if (nj != nsum) {
	i__1 = nsum;
	for (jj = 1; jj <= i__1; ++jj) {
	    j = jj;
	    if (! jt[jj - 1]) {
		goto L10;
	    }
/* L16: */
	}
    }
/*                                                                       ABBY3154 */
    i__1 = kt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = n6jn[i__ - 1];
	for (j = 1; j <= 6; ++j) {
	    jword_ref(j, i__) = jw_ref(j, i1);
/* L27: */
	}
/* L26: */
    }
/*                                                                       ABBY3161 */
    i__1 = nsum;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ik = jnsum[i__ - 1];
	i2 = jsum6[ik - 1];
	jtem6[i__ - 1] = i2;
	i__2 = i2;
	for (j = 1; j <= i__2; ++j) {
	    jtem4_ref(i__, j) = jsum4_ref(ik, j);
	    k = jsum5_ref(ik, j);
	    jtem5_ref(i__, j) = in6j[k - 1];
/* L29: */
	}
/* L28: */
    }
/*                                                                       ABBY3172 */
    i__1 = nsum;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i2 = jtem6[i__ - 1];
	jsum6[i__ - 1] = i2;
	i__2 = i2;
	for (j = 1; j <= i__2; ++j) {
	    i1 = jtem5_ref(i__, j);
	    j1 = jtem4_ref(i__, j);
	    jsum4_ref(i__, j) = j1;
	    jsum5_ref(i__, j) = i1;
	    jword_ref(j1, i1) = i__ + mp;
/* L41: */
	}
/* L40: */
    }
/*                                                                       ABBY3184 */
    return 0;
} /* sprate_ */

#undef sum6j
#undef t6j
#undef jt
#undef js
#undef inver
#undef jnsum
#undef jinv
#undef n6jn
#undef in6j
#undef jsumt
#undef cut
#undef jtem4
#undef jtem5
#undef jtem6
#undef nsum6j
#undef j6sum
#undef j6c
#undef j7c
#undef j8c
#undef j9c
#undef jwc
#undef j6
#undef j7
#undef j8
#undef j9
#undef jw
#undef sumvar
#undef mp
#undef j6cc
#undef j7cc
#undef j8cc
#undef j9cc
#undef j6p
#undef j7p
#undef j8p
#undef j9p
#undef jword
#undef nlsum
#undef nbj
#undef nb6j
#undef k6cp
#undef k7cp
#undef k8cp
#undef k9cp
#undef jsum6
#undef jsum4
#undef jsum5
#undef inv6j


#undef jsumt_ref
#undef jword_ref
#undef jsum5_ref
#undef jsum4_ref
#undef jtem5_ref
#undef jtem4_ref
#undef jw_ref

#define namsub (nam_1.namsub)
#define jdiag (graph_1.jdiag)
#define arr (graph_1.arr)
#define tab1 (graph_1.tab1)
#define il (graph_1.il)
#define ih (graph_1.ih)
#define npoint (graph_1.npoint)
#define nbnode (graph_1.nbnode)
#define ilast (graph_1.ilast)
#define npart (graph_1.npart)
#define icross (graph_1.icross)
#define j6c (argu_1.j6c)
#define j7c (argu_1.j7c)
#define jwc (argu_1.jwc)
#define j6 (argu_1.j6)
#define j7 (argu_1.j7)
#define kw (argu_1.kw)
#define sumvar (argu_1.sumvar)
#define mp (argu_1.mp)

/*                                                                       ABBY3187 */
/*                                                                       ABBY3188 */
/*     *****************                                                 ABBY3189 */
/* Subroutine */ int square_(void)
{
    /* Initialized data */

    static char name__[6] = "SQUARE";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    integer i__, k, l2, l4, k23, k32, it, jj1, jj3, jj2, it1, it2, it3, it4, 
	    ith, ilp, itl, ithl, itll, itmn, itmx, itmin, itmax;
    extern /* Subroutine */ int phase2_(integer *);


#define kw_ref(a_1,a_2) kw[(a_2)*6 + a_1 - 7]
#define arr_ref(a_1,a_2) arr[(a_2)*48 + a_1 - 49]
#define tab1_ref(a_1,a_2) tab1[(a_2)*60 + a_1 - 61]
#define jdiag_ref(a_1,a_2) jdiag[(a_2)*48 + a_1 - 49]

/*     *****************                                                 ABBY3191 */
/*                                                                       ABBY3192 */
/*  ***REDUCES A CIRCUIT OF ORDER 4 IN THE TWO CASES WHICH ARE LEFT      ABBY3193 */
/*  ***OVER BY POLYGN,NAMELY TWO DISCONNECTED GROUPS OF TWO POINTS       ABBY3194 */
/*  ***AND ONE GROUP OF TWO POINTS PLUS THE TWO ENDS OF THE AXIS.IN      ABBY3195 */
/*  ***THE LATTER, THE END OF THE AXIS IS TRANSFERRED TO THE BEGINNING.  ABBY3196 */
/*  ***IN THIS PROCESS,ONE SUMMATION VARIABLE AND TWO 6J SYMBOLS ARE     ABBY3197 */
/*  ***INTRODUCED.                                                       ABBY3198 */
/*                                                                       ABBY3199 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY3200 */
/*                                                                       ABBY3201 */
/*                                                                       ABBY3204 */
/*                                                                       ABBY3208 */
/*                                                                       ABBY3212 */
/*                                                                       ABBY3215 */
/*                                                                       ABBY3217 */
/*                                                                       ABBY3219 */
/*                                                                       ABBY3220 */
/*                                                                       ABBY3221 */
    s_copy(namsub, name__, (ftnlen)6, (ftnlen)6);
    ++mp;
    sumvar[mp - 1] = TRUE_;
    k = 1;
    it1 = npoint[0];
    it2 = npoint[1];
/*                                                                       ABBY3228 */
    if (icross == 1) {
	it3 = npoint[2];
	it4 = npoint[3];
	k23 = 3;
	k32 = 2;
    } else {
	it3 = npoint[3];
	it4 = npoint[2];
	k23 = 2;
	k32 = 3;
    }
/*                                                                       ABBY3240 */
    l4 = jdiag_ref(it2, 1);
/*                                                                       ABBY3242 */
    if (arr_ref(it2, 1) <= 0) {
	phase2_(&l4);
	arr_ref(it2, 1) = 1;
	arr_ref(it3, 1) = -1;
    }
/*                                                                       ABBY3248 */
    l2 = jdiag_ref(it1, 1);
    if (arr_ref(it1, 1) > 0) {
	phase2_(&l2);
    }
    ++jwc;
    kw_ref(1, jwc) = l4;
    kw_ref(2, jwc) = l2;
    kw_ref(3, jwc) = jdiag_ref(it2, 2);
    jj1 = jdiag_ref(it1, 3);
    kw_ref(4, jwc) = jj1;
    kw_ref(5, jwc) = mp;
    kw_ref(6, jwc) = jdiag_ref(it1, 2);
    if (arr_ref(it1, 2) < 0) {
	phase2_(&jdiag_ref(it1, 2));
    }
    ++jwc;
    kw_ref(1, jwc) = l4;
    kw_ref(2, jwc) = l2;
    jj3 = jdiag_ref(it3, k23);
    jj2 = jdiag_ref(it4, k32);
    kw_ref(3, jwc) = jj3;
    kw_ref(4, jwc) = jj2;
    kw_ref(5, jwc) = mp;
    kw_ref(6, jwc) = jdiag_ref(it3, k32);
    if (arr_ref(it3, k32) < 0) {
	phase2_(&jdiag_ref(it3, k32));
    }
    j6[j6c] = mp;
    j6c += 2;
    j6[j6c - 1] = mp;
/*                                                                       ABBY3273 */
    if (npart == 1) {
	itmin = it2;
	itmax = it3;
	itl = max(it3,it4);
	ith = ilast;
    } else {
	itmin = min(it2,it3);
	itmax = max(it2,it3);
	itmn = min(it1,it4);
	itmx = max(it1,it4);
	itl = max(itmin,itmn);
	ith = min(itmax,itmx);
    }
/*                                                                       ABBY3287 */
    tab1_ref(mp, 1) = itmin;
    tab1_ref(mp, 2) = itmax;
    jdiag_ref(it2, 1) = mp;
    jdiag_ref(it3, 1) = mp;
    jdiag_ref(it2, 3) = jj1;
    arr_ref(it2, 3) = arr_ref(it1, 3);
    jdiag_ref(it3, k32) = jj2;
    arr_ref(it3, k32) = arr_ref(it4, k32);
/*                                                                       ABBY3296 */
    if (icross == 1) {
	j7[j7c] = l2;
	j7[j7c + 1] = l4;
	phase2_(&l4);
	j7c += 3;
	j7[j7c - 1] = mp;
    } else {
	phase2_(&jj2);
    }
/*                                                                       ABBY3306 */
    itll = il[itl - 1];
    if (npart == 1 && icross == 1) {
	++itll;
    }
    ithl = il[ith - 1];
    if (npart == 2 && icross != 1) {
	--ithl;
    }
/*                                                                       ABBY3311 */
L5:
    i__1 = ithl;
    for (i__ = itll; i__ <= i__1; ++i__) {
	it = ih[i__ - 1];
	ilp = i__ - k;
	il[it - 1] = ilp;
	ih[ilp - 1] = it;
/* L6: */
    }
/*                                                                       ABBY3318 */
    if (ithl != nbnode) {
	itll = ithl + 2;
	if (itll <= nbnode) {
	    ithl = nbnode;
	    k = 2;
	    goto L5;
	}
    }
/*                                                                       ABBY3327 */
    if (npart != 2) {
	tab1_ref(jj1, 1) = ih[0];
	tab1_ref(jj1, 2) = ih[nbnode - 3];
    }
/*                                                                       ABBY3332 */
    return 0;
} /* square_ */

#undef namsub
#undef jdiag
#undef arr
#undef tab1
#undef il
#undef ih
#undef npoint
#undef nbnode
#undef ilast
#undef npart
#undef icross
#undef j6c
#undef j7c
#undef jwc
#undef j6
#undef j7
#undef kw
#undef sumvar
#undef mp


#undef jdiag_ref
#undef tab1_ref
#undef arr_ref
#undef kw_ref

#define cut (cutdig_1.cut)
#define j1 (couple_1.j1)
#define free (couple_1.free)
#define sumvar (argu_1.sumvar)

/*                                                                       ABBY3335 */
/*                                                                       ABBY3336 */
/*     **************************************                            ABBY3337 */
/* Subroutine */ int trdel_(integer *jj1, integer *jj2, integer *jj3, integer 
	*nbn, logical *fail)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i1, i2, i3;

/*     **************************************                            ABBY3339 */
/*                                                                       ABBY3340 */
/*  ***TEST FOR TRIANGULAR DELTA.IF NOT SATISFIED FAIL=.TRUE.            ABBY3341 */
/*                                                                       ABBY3342 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY3343 */
/*                                                                       ABBY3344 */
/*                                                                       ABBY3347 */
/*                                                                       ABBY3349 */
/*                                                                       ABBY3352 */
/*                                                                       ABBY3355 */
/*                                                                       ABBY3356 */
/*                                                                       ABBY3357 */
    if (sumvar[*jj1 - 1] || sumvar[*jj2 - 1] || sumvar[*jj3 - 1]) {
	return 0;
    }
    if (*nbn > 4) {
	cut = TRUE_;
    }
    if (! free[*jj1 - 1] && ! free[*jj2 - 1] && ! free[*jj3 - 1]) {
	i1 = j1[*jj1 - 1];
	i2 = j1[*jj2 - 1];
	i3 = j1[*jj3 - 1];
	if (i1 < (i__1 = i2 - i3, abs(i__1)) + 1 || i1 > i2 + i3 - 1) {
	    *fail = TRUE_;
	}
    }
/*                                                                       ABBY3366 */
    return 0;
} /* trdel_ */

#undef cut
#undef j1
#undef free
#undef sumvar

#define namsub (nam_1.namsub)
#define jdiag (graph_1.jdiag)
#define arr (graph_1.arr)
#define tab1 (graph_1.tab1)
#define il (graph_1.il)
#define ih (graph_1.ih)
#define npoint (graph_1.npoint)
#define nbnode (graph_1.nbnode)
#define ifirst (graph_1.ifirst)
#define ilast (graph_1.ilast)
#define jwc (argu_1.jwc)
#define kw (argu_1.kw)

/*                                                                       ABBY3369 */
/*                                                                       ABBY3370 */
/*     ***********************                                           ABBY3371 */
/* Subroutine */ int triang_(logical *fail)
{
    /* Initialized data */

    static char name__[6] = "TRIANG";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    integer i__, k12, k23, it, il2, il3, it1, it2, it3, jt1, ilp;
    extern /* Subroutine */ int trdel_(integer *, integer *, integer *, 
	    integer *, logical *), phase2_(integer *);


#define kw_ref(a_1,a_2) kw[(a_2)*6 + a_1 - 7]
#define arr_ref(a_1,a_2) arr[(a_2)*48 + a_1 - 49]
#define tab1_ref(a_1,a_2) tab1[(a_2)*60 + a_1 - 61]
#define jdiag_ref(a_1,a_2) jdiag[(a_2)*48 + a_1 - 49]

/*     ***********************                                           ABBY3373 */
/*                                                                       ABBY3374 */
/*  ***REDUCES A TRIANGLE HAVING ONE APEX AT EITHER END OF THE AXIS OF   ABBY3375 */
/*  ***THE FLAT DIAGRAM.                                                 ABBY3376 */
/*  ***THIS INTRODUCES ONE 6J SYMBOL AND SOME PHASE FACTORS .            ABBY3377 */
/*                                                                       ABBY3378 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY3379 */
/*                                                                       ABBY3380 */
/*                                                                       ABBY3383 */
/*                                                                       ABBY3387 */
/*                                                                       ABBY3391 */
/*                                                                       ABBY3394 */
/*                                                                       ABBY3396 */
/*                                                                       ABBY3398 */
/*                                                                       ABBY3399 */
/*                                                                       ABBY3400 */
    s_copy(namsub, name__, (ftnlen)6, (ftnlen)6);
    it1 = npoint[0];
    it2 = npoint[1];
    it3 = npoint[2];
    ++jwc;
    kw_ref(1, jwc) = jdiag_ref(it3, 2);
    kw_ref(2, jwc) = jdiag_ref(it2, 3);
    kw_ref(3, jwc) = jdiag_ref(it3, 1);
    if (arr_ref(it3, 1) > 0) {
	phase2_(&kw_ref(3, jwc));
    }
    kw_ref(4, jwc) = jdiag_ref(it2, 1);
    if (arr_ref(it2, 1) < 0) {
	phase2_(&kw_ref(4, jwc));
    }
    k23 = 3;
    if (it1 == ifirst) {
	k23 = 2;
    }
    kw_ref(5, jwc) = jdiag_ref(it1, k23);
    kw_ref(6, jwc) = jdiag_ref(it3, 3);
    trdel_(&kw_ref(1, jwc), &kw_ref(2, jwc), &kw_ref(5, jwc), &nbnode, fail);
    if (*fail) {
	goto L15;
    }
    if (arr_ref(it3, 3) > 0) {
	phase2_(&kw_ref(6, jwc));
    }
    jt1 = kw_ref(5, jwc);
    jdiag_ref(it3, 1) = jt1;
    jdiag_ref(it3, 3) = kw_ref(2, jwc);
    arr_ref(it3, 1) = arr_ref(it1, k23);
    arr_ref(it3, 3) = arr_ref(it2, 3);
/*                                                                       ABBY3424 */
    if (it1 != ifirst) {
	tab1_ref(jt1, 1) = it3;
	tab1_ref(jt1, 2) = ih[nbnode - 2];
	k12 = 1;
    } else {
	tab1_ref(jt1, 1) = ih[1];
	tab1_ref(jt1, 2) = it3;
	k12 = 2;
    }
/*                                                                       ABBY3434 */
    il3 = il[it3 - 1];
/*                                                                       ABBY3436 */
    if (it1 != ilast) {
	il2 = il[it2 - 1] - 1;
/*                                                                       ABBY3439 */
	i__1 = il2;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    it = ih[i__ - 1];
	    ilp = i__ - 1;
	    il[it - 1] = ilp;
	    ih[ilp - 1] = it;
/* L2: */
	}
    }
/*                                                                       ABBY3447 */
    i__1 = nbnode;
    for (i__ = il3; i__ <= i__1; ++i__) {
	it = ih[i__ - 1];
	ilp = i__ - k12;
	il[it - 1] = ilp;
	ih[ilp - 1] = it;
/* L1: */
    }
/*                                                                       ABBY3454 */
L15:
    return 0;
} /* triang_ */

#undef namsub
#undef jdiag
#undef arr
#undef tab1
#undef il
#undef ih
#undef npoint
#undef nbnode
#undef ifirst
#undef ilast
#undef jwc
#undef kw


#undef jdiag_ref
#undef tab1_ref
#undef arr_ref
#undef kw_ref


/*                                                                       ABBY3457 */
/*                                                                       ABBY3458 */
/*     *****************************************************             ABBY3459 */
/* Subroutine */ int var_(integer *jn, integer *jns, integer *jnc, integer *
	jnsc, integer *jbc, logical *sumvar, integer *mp, integer *m, integer 
	*inver)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, j, i1, jbbc;

/*     *****************************************************             ABBY3461 */
/*                                                                       ABBY3462 */
/*  ***TEST FOR VARIABLE CHARACTER AND PUT IN JNS IF YES,AND JN NOW      ABBY3463 */
/*  ***CONTAINS 0.                                                       ABBY3464 */
/*                                                                       ABBY3465 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY3466 */
/*                                                                       ABBY3467 */
/*                                                                       ABBY3470 */
/*                                                                       ABBY3473 */
/*                                                                       ABBY3474 */
/*                                                                       ABBY3475 */
    /* Parameter adjustments */
    --jns;
    --jn;
    --inver;
    --sumvar;

    /* Function Body */
    *jnsc = 0;
    if (*jbc != *jnc) {
	jbbc = *jbc + 1;
/*                                                                       ABBY3479 */
	i__1 = *jnc;
	for (i__ = jbbc; i__ <= i__1; ++i__) {
	    i1 = jn[i__];
	    if (sumvar[i1]) {
		++(*jnsc);
		j = inver[i1];
		jns[*jnsc] = j;
		jn[i__] = *m;
	    }
/* L1: */
	}
    }
/*                                                                       ABBY3490 */
    return 0;
} /* var_ */
#define j23 (tree_1.j23)
#define ial (build_1.ial)
#define if1 (build_1.if1)
#define if2 (build_1.if2)

/*                                                                       ABBY3493 */
/*                                                                       ABBY3494 */
/*     ******************************                                    ABBY3495 */
/* Subroutine */ int way_(integer *l, integer *ka, integer *kb, integer *ich, 
	integer *nb)
{
    integer i1, i2, k1, k2, l1, l2, i3, i4, ia, ib, la, lb, nb1, lc1, lc2, 
	    ji1, ji2, ji3, ji4, nbm, nbp;
    extern /* Subroutine */ int neibor_(integer *, integer *, integer *), 
	    otherj_(integer *, integer *, integer *, integer *, integer *);


#define j23_ref(a_1,a_2) j23[(a_2)*24 + a_1 - 25]

/*     ******************************                                    ABBY3497 */
/*                                                                       ABBY3498 */
/*  ***TESTS ONE STEP FORWARD  IF THE WAY IS FREE.FIRST AND SECOND       ABBY3499 */
/*  ***ARGUMENTS ARE INTERCHANGED OR NOT ACCORDING TO ICH=-1,OR +1       ABBY3500 */
/*                                                                       ABBY3501 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY3502 */
/*                                                                       ABBY3503 */
/*                                                                       ABBY3506 */
/*                                                                       ABBY3509 */
/*                                                                       ABBY3512 */
/*                                                                       ABBY3514 */
/*                                                                       ABBY3515 */
    k1 = j23_ref(*l, *ka);
    k2 = j23_ref(*l, *kb);
    *nb = ial[k1 - 1] + ial[k2 - 1] - 1;
    if (*nb < 0) {
	goto L3;
    } else if (*nb == 0) {
	goto L2;
    } else {
	goto L8;
    }
L2:
    nb1 = ial[k1 - 1] - ial[k2 - 1];
    if (nb1 >= 0) {
	goto L8;
    } else {
	goto L9;
    }
L3:
    otherj_(l, &k1, &l1, &lc1, &la);
    otherj_(l, &k2, &l2, &lc2, &lb);
    neibor_(&lc1, &i1, &i2);
    neibor_(&lc2, &i3, &i4);
    ji1 = j23_ref(l1, i1);
    ji2 = j23_ref(l1, i2);
    ji3 = j23_ref(l2, i3);
    ji4 = j23_ref(l2, i4);
    ia = ial[ji1 - 1] + ial[ji2 - 1];
    ib = ial[ji3 - 1] + ial[ji4 - 1];
    nbp = ib + ia + 1;
    nbm = ib - ia;
    switch (nbp) {
	case 1:  goto L8;
	case 2:  goto L4;
	case 3:  goto L5;
	case 4:  goto L4;
	case 5:  goto L6;
    }
L4:
    if (nbm >= 0) {
	goto L8;
    } else {
	goto L9;
    }
L5:
    if (nbm < 0) {
	goto L9;
    } else if (nbm == 0) {
	goto L6;
    } else {
	goto L8;
    }
L6:
    if (ji3 == if1 || ji3 == if2 || ji4 == if1 || ji4 == if2) {
	goto L9;
    }
L8:
    *ich = 1;
    goto L10;
L9:
    *ich = -1;
L10:
    return 0;
} /* way_ */

#undef j23
#undef ial
#undef if1
#undef if2


#undef j23_ref

#define j23 (tree_1.j23)
#define arrow (tree_1.arrow)
#define line (tree_1.line)
#define lcol (tree_1.lcol)
#define tabs (tree_1.tabs)
#define nbtr (tree_1.nbtr)
#define cut (cutdig_1.cut)
#define m (couple_1.m)
#define j1 (couple_1.j1)
#define free (couple_1.free)
#define j6c (argu_1.j6c)
#define j9c (argu_1.j9c)
#define j6 (argu_1.j6)
#define j9 (argu_1.j9)
#define sumvar (argu_1.sumvar)
#define ial (build_1.ial)
#define nzero (zer_1.nzero)
#define jzero (zer_1.jzero)

/*                                                                       ABBY3543 */
/*                                                                       ABBY3544 */
/*     **************************                                        ABBY3545 */
/* Subroutine */ int zero_(integer *j, integer *jz, logical *fail)
{
    /* Initialized data */

    static char name__[6] = "ZERO  ";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, k, l, k1, l1, l2, k2, lc, jt, jj1, jj2, lo1, lo2, lin, jtf, 
	    jjx, jjz, lco1, lco2;
    extern /* Subroutine */ int delta_(integer *, integer *, logical *);
    logical nocut;
    extern /* Subroutine */ int phase2_(integer *), neibor_(integer *, 
	    integer *, integer *), otherj_(integer *, integer *, integer *, 
	    integer *, integer *), printj_(char *, integer *, ftnlen);


#define j23_ref(a_1,a_2) j23[(a_2)*24 + a_1 - 25]
#define line_ref(a_1,a_2) line[(a_2)*60 + a_1 - 61]
#define lcol_ref(a_1,a_2) lcol[(a_2)*60 + a_1 - 61]
#define arrow_ref(a_1,a_2) arrow[(a_2)*24 + a_1 - 25]

/*     **************************                                        ABBY3547 */
/*                                                                       ABBY3548 */
/*  ***SUPPRESSES ONE LINE AND TWO NODES OF THE UNSTRUCTURED GRAPH       ABBY3549 */
/*  ***INTRODUCES  ZEROS IN THE TRIADS J23.AS A CONSEQUENCE THE OTHER    ABBY3550 */
/*  ***TWO ARGUMENTS OF THE TRIAD ARE PUT EQUAL.IF THERE WAS ALREADY     ABBY3551 */
/*  ***A ZERO IN THE TRIAD WHICH IS CHANGED,IT IS A SPECIAL CASE.        ABBY3552 */
/*                                                                       ABBY3553 */
/*  ***IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                ABBY3554 */
/*                                                                       ABBY3555 */
/*                                                                       ABBY3558 */
/*                                                                       ABBY3560 */
/*                                                                       ABBY3562 */
/*                                                                       ABBY3564 */
/*                                                                       ABBY3570 */
/*                                                                       ABBY3573 */
/*                                                                       ABBY3576 */
/*                                                                       ABBY3578 */
/*                                                                       ABBY3579 */
/*                                                                       ABBY3580 */
    nocut = FALSE_;
    nzero = 0;
/*                                                                       ABBY3583 */
    if (*j >= 1) {
	otherj_(&c__0, jz, &lin, &lc, &k1);
	i__ = nzero;
	goto L8;
    }
/*                                                                       ABBY3589 */
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (j1[i__ - 1] != 1 || free[i__ - 1] || ial[i__ - 1] <= 1) {
	    goto L11;
	}
	++nzero;
	jzero[nzero - 1] = i__;
L11:
	;
    }
/*                                                                       ABBY3595 */
    nocut = TRUE_;
    ++m;
    j1[m - 1] = 1;
    sumvar[m - 1] = FALSE_;
    free[m - 1] = FALSE_;
    if (nzero == 0) {
	goto L7;
    }
    printj_(name__, &c__1, (ftnlen)6);
    i__ = 0;
L1:
    ++i__;
    *jz = jzero[i__ - 1];
    *j = 0;
L13:
    ++(*j);
    lin = line_ref(*jz, *j);
    if (tabs[lin - 1]) {
	goto L2;
    }
    lc = lcol_ref(*jz, *j);
L8:
    neibor_(&lc, &l1, &l2);
    jj1 = j23_ref(lin, l1);
    jj2 = j23_ref(lin, l2);
/*                                                                       ABBY3614 */
    if (jj1 == jj2) {
	++j6c;
	j6[j6c - 1] = jj1;
	goto L10;
    }
/*                                                                       ABBY3620 */
    delta_(&jj1, &jj2, fail);
    if (*fail) {
	goto L7;
    }
/*                                                                       ABBY3623 */
    if (j1[jj1 - 1] != 1 && j1[jj2 - 1] != 1) {
	goto L15;
    }
    if (j1[jj1 - 1] < j1[jj2 - 1]) {
	goto L15;
    }
    if (j1[jj1 - 1] > j1[jj2 - 1]) {
	goto L19;
    }
/*                                                                       ABBY3627 */
    if (nzero != 0) {
	i__1 = nzero;
	for (jjx = i__; jjx <= i__1; ++jjx) {
	    jjz = jzero[jjx - 1];
	    if (jj1 == jjz) {
		goto L15;
	    }
	    if (jj2 == jjz) {
		goto L19;
	    }
/* L17: */
	}
    }
/*                                                                       ABBY3635 */
    goto L15;
L19:
    jjz = jj2;
    jj2 = jj1;
    jj1 = jjz;
/*                                                                       ABBY3640 */
L15:
    otherj_(&lin, &jj1, &lo1, &lco1, &k1);
    otherj_(&lin, &jj2, &lo2, &lco2, &k2);
    ++j9c;
    j9[j9c - 1] = jj1;
    j23_ref(lo2, lco2) = jj1;
    line_ref(jj1, k1) = lo2;
    lcol_ref(jj1, k1) = lco2;
/*                                                                       ABBY3648 */
L10:
    if (arrow_ref(lin, l1) < arrow_ref(lin, l2)) {
	phase2_(&jj1);
    } else {
	if (arrow_ref(lin, l1) == arrow_ref(lin, l2)) {
	    arrow_ref(lo1, lco1) = 1;
	    arrow_ref(lo2, lco2) = -1;
	}
    }
/*                                                                       ABBY3657 */
    tabs[lin - 1] = TRUE_;
    --nbtr;
    if (nbtr == 0) {
	goto L7;
    }
    if (lo1 != lo2) {
	goto L2;
    }
    l = 6 - lco1 - lco2;
    jt = j23_ref(lo1, l);
    if (j1[jt - 1] == 1 && ! free[jt - 1]) {
	goto L2;
    }
    delta_(&jt, &m, fail);
    if (*fail) {
	goto L7;
    }
    neibor_(&l, &l1, &l2);
    jtf = j23_ref(lo1, l1);
    if (arrow_ref(lo1, l1) < arrow_ref(lo1, l2)) {
	phase2_(&jtf);
    }
    ++j6c;
    j6[j6c - 1] = jtf;
    --nbtr;
    tabs[lo1 - 1] = TRUE_;
    otherj_(&lo1, &jt, &lin, &lc, &k);
    goto L8;
L2:
    if (*j == 1) {
	goto L13;
    }
/*                                                                       ABBY3677 */
    if (nbtr != 0) {
	if (i__ < nzero) {
	    goto L1;
	}
    }
/*                                                                       ABBY3681 */
L7:
    printj_(name__, &c__4, (ftnlen)6);
    if (nocut) {
	cut = FALSE_;
    }
/*                                                                       ABBY3684 */
    return 0;
} /* zero_ */

#undef j23
#undef arrow
#undef line
#undef lcol
#undef tabs
#undef nbtr
#undef cut
#undef m
#undef j1
#undef free
#undef j6c
#undef j9c
#undef j6
#undef j9
#undef sumvar
#undef ial
#undef nzero
#undef jzero


#undef arrow_ref
#undef lcol_ref
#undef line_ref
#undef j23_ref


int second_(real *data)
{
	return 0;
}

#define MANGM	(60)
#define MTRIAD	(12)

#define ibug3 (debug_1.ibug3)
#define ird (inform_1.ird)
#define ipd (inform_1.ipd)
#define nograf (timeg_1.nograf)
#define tgraf (timeg_1.tgraf)
#define nogen (timeg_1.nogen)
#define tgen (timeg_1.tgen)
#define norac (timeg_1.norac)
#define trac (timeg_1.trac)

#define j2_ref(a_1,a_2) couple_1.j2[(a_2)*12 + a_1 - 13]
#define j3_ref(a_1,a_2) couple_1.j3[(a_2)*12 + a_1 - 13]

double do_njgraf(int m,int n,int j1[MANGM],int j2[MTRIAD][3],int j3[MTRIAD][3],int isfree[MANGM])
{
	real recup;
	logical fail;

        factt_();

        ird = 5;
        ipd = 6;

        norac = 0;
        tgraf = 0.f;
        nograf = 0;
        trac = 0.f;
        tgen = 0.f;
        nogen = 0;
        fail = FALSE_;

	ibug3=0;

	//assert(MANGM==)
	//assert(MTRIAD==)

	couple_1.m=m;
	couple_1.n=n+1;

	for(int c=0;c<MANGM;c++)
	{
		couple_1.j1[c]=j1[c];
		couple_1.free[c]=isfree[c];
	}

	for(int c=0;c<MTRIAD;c++)
	{
		for(int d=0;d<3;d++)
		{
			j2_ref(c+1,d+1)=j2[c][d];
			j3_ref(c+1,d+1)=j3[c][d];
		}
	}

	njgraf_(&recup,&fail);

	return recup;
}
