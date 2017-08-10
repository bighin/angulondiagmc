#ifndef __NJFORMUL_H__
#define __NJFORMUL_H__

/*
 *	njformul.h
 *	----------
 */

#include <stdio.h>
#include <stdlib.h>

#define MAXJ            35
#define MAXK            10
#define MAXSIXJ         20
#define MAXFLEN         5
#define MAXS            4

#define FALSE	        0
#define TRUE	        1
#define OFF             0
#define ON              1
#define SHALLOW         1
#define DEEP            2
#define VERYDEEP        4

typedef int SET[MAXS];

typedef struct node {
	int index;
	SET leafs;
	struct node *left, *right;
} NODE;

typedef int SIXJ[7];

typedef struct formula {
        int nrjs, nrks, nrsixjs;
        int jsigns[MAXJ], jsqrts[MAXJ];
        int ksigns[MAXK], ksqrts[MAXK];
        SIXJ sixjs[MAXSIXJ];
        int js[MAXJ];
        SIXJ asixjs[MAXSIXJ];
        int kp[MAXK], kptmp[MAXK];
	int ordered;
} FORMULA;

typedef struct deltas {
        int n;
	int d[MAXJ][3];
} DELTAS;

extern int AUTONUM, SEARCH;

void read_njsym (int *pnrjs, NODE **pbra, NODE **pket, DELTAS *pd);
/* Reads a recoupling coefficient in the style of NJSYM,
   from the standard input device. */

void read_expr (int *pnrjs, NODE **pbra, NODE **pket, DELTAS *pd);
/* Reads an expression in bra-ket notation and returning coupling trees.
   If AUTONUM is ON, needs number of leaf nodes, given in *pnrjs;
   returns total number of nodes in *pnrjs.
   If an error occurs while reading: clears line and retries. */

void write_formula (FORMULA *f);
/* Writes the formula f to the standard output device.  */

FORMULA *generate_formula (int nrjs, NODE *bra, NODE *ket);
/* Transforms expression bra into ket, building formula f;
   returns list of delta functions for j coefficients in *pd,
   returns NULL if transformation is not successful.
   Possible searchlevels : SHALLOW, DEEP, VERYDEEP. */

#endif //__NJFORMUL_H__
