/* ------------------------- cut here ------------------------------------ */

/*
 *       njsummat.h
 *       ----------
 */

#include "njformul.h"

extern int count;
/* Counter for number of 6-j's computed,
   initialised to 0 by evaluate_formula, updated in each call by comp_6j.  */

void comp_fac ();
/* Build look-up table with logarithms of factorials. */

void read_jvals (FORMULA *f);
/* Reads the actual values for the j's (value read = 2j+1)
   for a given formula *f and stores them in f->js */

void copy_jvals (FORMULA *orig, FORMULA *copy);
/* Copy j values of original formula to copy formula. */

void split (DELTAS delta, int *pnrtrees, NODE *bras[], NODE *kets[]);
/* Splits a pair of bra and ket trees (given in 0th component)
   into a set of independent pairs of bra and ket trees.
   Given are deltas, discovered while reading.
   Returns also the number of pairs of trees obtained. */

double comp_6j (SIXJ jv);
/* Compute a 6-j coefficient with given parameters. */

int trian (FORMULA *f);
/* Checks the triangle conditions involving the angular momenta.
   Returns 1 if they are satisfied, 0 otherwise. */

void order (FORMULA *f);
/* Orders the 6-j's in the formula w.r.t. the summation variables */

double evaluate_formula (FORMULA *f);
/* Evaluates a given formula *f with built structures and actual j-values */

