/* -------------------- cut here -------------------------------------*/

/*
 *       njsummat.c
 *       ----------
 */

#include "njsummat.h"
#include <math.h>

int count=0;

/*
 *	BUILD LOOK_UP TABLE OF FACTORIALS
 */

#define MAXFAC	501				/* Size of look-up table */

static double fac[MAXFAC];			/* Look-up table */

void comp_fac ()
/* Build look-up table with logarithms of factorials. */
{
	int i;
        double *pfac, tmp;

	*(pfac=fac) = 0; 			/* fac[0] = 1; */
	*++pfac = 0;				/* fac[1] = 1; */
	for (i=2; i<MAXFAC; i++) {
	  tmp = *pfac + log((double) i);
	  *++pfac = tmp;			/* fac[i] = fac[i-1] * i; */
	}
}

/*
 *        COMPUTATION OF 6-J COEFFICIENTS
 */

#define PSIGN(a) ((a)%2 ? -1 : 1)
#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a > b ? a : b)
#define DELTA6(a,b,c) (fac[(a+b-c)/2] + fac[(a-b+c)/2] + fac[(-a+b+c)/2]\
                       - fac[(a+b+c)/2+1])

double comp_6j (SIXJ jv)
/* Compute a 6-j coefficient with given parameters. */
{
        int j1, j2, j3, l1, l2, l3, i;
	int kfirst, klast, k, sign, kf2, kf3, kl1, kl2, kl3, kl4, kn;
	double res, sum, term, factor;

	count++;

        j1 = jv[1]; j2 = jv[2]; j3 = jv[3];
        l1 = jv[4]; l2 = jv[5]; l3 = jv[6];

	sign = PSIGN((j1+j2+l1+l2)/2);
	factor = (DELTA6(j1,j2,j3) + DELTA6(l1,l2,j3) + DELTA6(l1,j2,l3)
	          + DELTA6(j1,l2,l3)) / 2;
	res = sign * exp(factor);

	kf2 = (j1+l1-j3-l3)/2; kf3 = (j2+l2-j3-l3)/2; kn = (j1+j2+l1+l2)/2+1;
	kl1 = (j1+j2-j3)/2; kl2 = (l1+l2-j3)/2;
	kl3 = (j1+l2-l3)/2; kl4 = (l1+j2-l3)/2;
	kfirst = MAX(MAX(0,kf2),kf3);
	klast = MIN(MIN(kl1,kl2),MIN(kl3,kl4));

	sign = PSIGN(kfirst);
	kf2 = kfirst - kf2; kf3 = kfirst - kf3; kn -= kfirst;
	kl1 -= kfirst; kl2 -= kfirst; kl3 -= kfirst; kl4 -= kfirst;
	factor = fac[kn] - fac[kfirst] - fac[kf2] - fac[kf3]
	         - fac[kl1] - fac[kl2] - fac[kl3] - fac[kl4];
	res *= sign * exp(factor);

	sum = 1; term = 1;
	for (k=kfirst+1; k<=klast; k++) {
	  term *= - (kl1-- * kl3-- * kl2-- * kl4--);
	  term /= (k * ++kf2 * ++kf3 * kn--);
	  sum += term;
	}
	res *= sum;
	return res;
}

/*
 *      MISCELLANEOUS
 */

void read_jvals (FORMULA *f)
/* Reads the actual values for the j's (value read = 2j+1)
   for a given formula *f and stores them in f->js */
{
        int i, j;
	char s[80];

	count = 0;
        printf("Give the array of %d j values\n", f->nrjs);
        for (i=1; i<=f->nrjs; i++) {
          scanf("%d", &j);
          f->js[i] = j-1;
        }
        for (i=1; i<=f->nrjs; i++)
          scanf ("%s", s);
}

void copy_jvals (FORMULA *orig, FORMULA *copy)
/* Copy j values of original formula to copy formula. */
{
        int i;

	copy->nrjs = orig->nrjs;
	for (i=1; i<=orig->nrjs; i++)
	  copy->js[i] = orig->js[i];
}

NODE *find_and_delete_node (int j, NODE *tree)
{
       NODE *current;

       if (tree->left == NULL)
	 return NULL;
       else if (j == tree->index) {
	 current = (NODE *) malloc (sizeof(NODE));
	 current->index = j;
	 copy_set (tree->leafs, current->leafs);
	 current->left = tree->left; current->right = tree->right;
	 tree->left = NULL; tree->right = NULL;
	 return current;
       }
       else if ((current = find_and_delete_node (j, tree->left)) != NULL)
	 return current;
       else
	 return find_and_delete_node (j, tree->right);
}

void split (DELTAS delta, int *pnrtrees, NODE *bras[], NODE *kets[])
/* Splits a pair of bra and ket trees (given in 0th component)
   into a set of independent pairs of bra and ket trees.
   Given are deltas, discovered while reading.
   Returns also the number of pairs of trees obtained. */
{
        int i;
	NODE *j1, *j2;

        for (i = 1; i <= delta.n; i++) {
	  j1 = find_and_delete_node (delta.d[i][2], bras[0]);
	  j2 = find_and_delete_node (delta.d[i][2], kets[0]);
	  bras[i] = j1; kets[i] = j2;
        }
	rebuild_leafs (bras[0]); rebuild_leafs (kets[0]);
	*pnrtrees = delta.n;
}

/*
 *        CHECKING TRIANGLE CONDITIONS
 */

int tr1 (s1, s2, s3)
{
      if (-s1 + s2 + s3 < 0) return 0;
      if (s1 - s2 + s3 < 0) return 0;
      if (s1 + s2 - s3 < 0) return 0;
      if ((s1 + s2 + s3) % 2) return 0;
      return 1;
}

int trian (FORMULA *f)
/* Checks the triangle conditions involving the angular momenta.
   Returns 1 if they are satisfied, 0 otherwise. */
{
         int i,*ptr;

         for (i=1; i<=f->nrsixjs; i++) {
           ptr=f->sixjs[i];
           if (ptr[1]>0 && ptr[2]>0 && ptr[3]>0)
	     if (!tr1(f->js[ptr[1]], f->js[ptr[2]], f->js[ptr[3]]))
               return 0;
           if (ptr[1]>0 && ptr[5]>0 && ptr[6]>0)
	     if (!tr1(f->js[ptr[1]], f->js[ptr[5]], f->js[ptr[6]]))
               return 0;
           if (ptr[2]>0 && ptr[4]>0 && ptr[6]>0)
	     if (!tr1(f->js[ptr[2]], f->js[ptr[4]], f->js[ptr[6]]))
               return 0;
           if (ptr[3]>0 && ptr[4]>0 && ptr[5]>0)
	     if (!tr1(f->js[ptr[3]], f->js[ptr[4]], f->js[ptr[5]]))
               return 0;
         }
         return 1;
}

/*
 *        ORDERING THE FORMULA
 */

void order(FORMULA *f)
/* Orders the 6-j's in the formula w.r.t. the summation variables */
{
        int i, j, k, k1, tmp, comp[MAXSIXJ];

        for (i=0; i<=f->nrks; i++) f->kp[i] = 0;
        for (i=1; i<=f->nrsixjs; i++) {
          comp[i]=0;
          for (j=1; j<=6; j++)
 	    if (f->sixjs[i][j]<=0)
              comp[i] = MAX(comp[i], abs(f->sixjs[i][j]));
          f->kp[comp[i]]++;
        }
        for (i=1; i<=f->nrks; i++)
          f->kp[i] = f->kp[i] + f->kp[i-1];
        for (i=1; i<=f->nrsixjs-1; i++) {
          if (comp[i+1] < comp[i])
            for (k=1; k<=i; k++)
              if (comp[k] > comp[i+1]){
                /* interchange 6j arrays k and i+1 */
                for (k1=1; k1<=6; k1++){
                   tmp = f->sixjs[k][k1];
		   f->sixjs[k][k1] = f->sixjs[i+1][k1];
                   f->sixjs[i+1][k1] = tmp;
                }
	        tmp = comp[k]; comp[k] = comp[i+1]; comp[i+1] = tmp;
              }
        }
	f->ordered = TRUE;
}

/*
 *        EVALUATION OF A SUMMATION FORMULA
 */

void range (int s1, int s2, int s3, int s_ind, int *plow, int *pup)
/* Update *plow and *pup for summation variable s_ind
   using triangle conditions s1,s2,s3  */
{
	int tmp;

        if (s1>=0 && s2>=0 && s3>=0) return;
	else if ((s1<0 && s2<0) || (s1<0 && s3<0) || (s2<0 && s3<0))
	  return;
        else if (s1<0 && s2<0 && s3<0) return;
	else if (s2<0) { tmp=s1; s1=s2; s2=tmp; }
	else if (s3<0) { tmp=s1; s1=s3; s3=tmp; }
        if (s1 == s_ind) {
	  if ((tmp = abs(s2-s3)) > *plow) *plow = tmp;
	  if ((tmp = s2 + s3) < *pup) *pup = tmp;
        }
}

double r_eval_formula (FORMULA *f, int it, int sign)
/* Recursive evaluation of a formula *f */
{
	int low, upp, nsign;
	int i, k, l, isum, *row;
	double coeff, sqrtfact, prod6j;

	if (it == f->nrks + 1) {
	  if (sign / 2 % 2 == 0) coeff = 1; else coeff = -1;
	  return coeff;
	}
	else {
	  /* Determine the range of the summation variable (it)  */
	  low = 0; upp = 2000;
          isum = f->kp[it-1] + 1;
	  for (i=isum; i<=f->nrsixjs; i++) {
             row = f->asixjs[i];
	     range(row[1], row[2], row[3], -it, &low, &upp);
	     range(row[1], row[5], row[6], -it, &low, &upp);
	     range(row[3], row[4], row[5], -it, &low, &upp);
	     range(row[2], row[4], row[6], -it, &low, &upp);
	  }
          /* Evaluation summation over variable it */
	  coeff = 0.0;
	  for (k=low; k<=upp; k+=2) {
	     nsign = k * f->ksigns[it];
	     sqrtfact = pow (sqrt((double)(k+1)), (double)f->ksqrts[it]);
	     for (i=isum; i<=f->nrsixjs; i++)
	       for (l=1; l<=6; l++)
		 if (f->sixjs[i][l] == -it) f->asixjs[i][l] = k;
             prod6j = 1.0;
             for (l=isum; l<=f->kp[it]; l++)
               prod6j *= comp_6j (f->asixjs[l]);
	     coeff += sqrtfact * prod6j *
                      r_eval_formula (f, it+1, sign+nsign);
	  }
          /* Restore old situation */
          for (i=isum; i<=f->nrsixjs; i++) {
            for (l=1; l<=6; l++) {
              if (f->sixjs[i][l] == -it) f->asixjs[i][l] = -it;
            }
          }
	  return coeff;
	}
}

double evaluate_formula (FORMULA *f)
/* Evaluates a given formula *f with built structures and actual j-values */
{
	int i, j, signfact;
	double sqrtfact, prod6j, coeff;

        /* Check triangle conditions on triads with angular momenta */
	if (trian(f) == 0) return 0.0;
        /* Order 6-j coefficients and build kp */
	if (f->ordered == FALSE) order(f);
        /* Initialise asixjs */
	for (i=1; i<=f->nrsixjs; i++) {
	  for (j=1; j<=6; j++) {
	    if (f->sixjs[i][j] < 0) {
	      f->asixjs[i][j] = f->sixjs[i][j];
	    }
	    else
	      f->asixjs[i][j] = f->js[f->sixjs[i][j]];
	  }
	}
        /* Compute sign and sqrt factors for angular momenta */
        sqrtfact = 1.0; signfact = 0;
	for (i=1; i<=f->nrjs; i++) {
	  if (f->jsqrts[i] != 0)
	    sqrtfact *= pow(sqrt((double)(f->js[i]+1)),
                            (double)f->jsqrts[i]);
          signfact += f->jsigns[i] * f->js[i];
        }
        /* Pull out 6-j's with only angular momenta  */
        prod6j = 1.0;
        for (i=1; i<=f->kp[0]; i++)
          prod6j *= comp_6j (f->asixjs[i]);
        /* Check if there is no sum */
        if (f->nrks == 0) {
          coeff = sqrtfact * prod6j;
	  if ((signfact / 2) % 2 != 0) coeff *= -1;
          return coeff;
        }
	else { /* start recursive evaluation process */
          coeff = sqrtfact * prod6j * r_eval_formula (f, 1, signfact);
          return coeff;
	}
}
