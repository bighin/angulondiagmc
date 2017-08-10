/*
 *	njformul.c
 *	----------
 *	Generates a summation formula for a general recoupling coefficient.
 *
 *	V.Fack, March-June 1994.
 */

#include "njformul.h"

int AUTONUM = ON, SEARCH = DEEP;

static int currentj = 0, left_first = TRUE;

/*
 *      Debugging
 *      ---------
 */

static FILE *debugfile;
static int debugging = OFF;

void init_debug (char *dname)
{
        debugging = ON; debugfile = fopen (dname, "w");
}

void close_debug ()
{
        fprintf (debugfile, "\n-----------------------------------\n");
        fclose (debugfile); debugging = OFF;
}

/*
 *	Basic operations on sets
 *	------------------------
 */

void clear_set (SET l)
{
	int i;

	for (i=0; i<MAXS; i++) *l++ = 0;
}

void add_element (int i, SET a)
{
	int s = 8*sizeof(int);

	a[i/s] |= 1 << (i%s);
}

void add_sets (SET a, SET b, SET sum)
{
	int i;

	for (i=0; i<MAXS; i++) { *sum++ = *a++ | *b; b++; }
}

int equal_sets (SET a, SET b)
{
	int i;

	for (i=0; i<MAXS; i++) if (*a++ != *b++) return FALSE;
	return TRUE;
}

int empty_set (SET a)
{
	int i;

	for (i=0; i<MAXS; i++) if (*a++) return FALSE;
	return TRUE;
}

void copy_set (SET a, SET b)
{
        int i;

        for (i=0; i<MAXS; i++) *b++ = *a++;
}

/*
 *      Basic operations on queues
 *      --------------------------
 */

#define QMAX    1024

typedef struct hqueue {
  int first, last;
  NODE *queue[QMAX];
} QUEUE;

QUEUE *init_queue ()
{
        QUEUE *q;

        q = (QUEUE *) malloc (sizeof (QUEUE));
        q->first = 0; q->last = -1;
        return q;
}

void add_to_queue (QUEUE *q, NODE *n)
{
        if (n->left == NULL) return;
        q->last++;
        if (q->last == QMAX) {
          printf ("ERROR : queue overflow\n"); exit(-1);
        }
        q->queue[q->last] = n;
}

void get_from_queue (QUEUE *q, NODE **pn)
{
        if (q->first <= q->last) {
          *pn = q->queue[q->first];
          q->first++;
        }
        else
          *pn = NULL;
}

void free_queue (QUEUE *q)
{
        free (q);
}

/*
 *      Working with coupling trees and recoupling coefficients
 *      -------------------------------------------------------
 */

void file_write_tree (FILE *ofile, NODE *tree)
/* Writes a coupling tree to a given file,
   which has to be opened already */
{
	if (tree->left != NULL) {
	  fprintf (ofile, "(");
	  file_write_tree (ofile, tree->left);
	  fprintf (ofile, ",");
	  file_write_tree (ofile, tree->right);
	  fprintf (ofile, ")");
	}
        fprintf (ofile, "%d", tree->index);
}

void file_write_expr (FILE *ofile, NODE *bra, NODE *ket)
/* Writes a recoupling coefficient bra-ket to a file,
   which has to be opened already */
{
        fprintf (ofile, "\n----------------------------------\n\n");
        fprintf (ofile, "recoupling coefficient = < ");
        file_write_tree (ofile, bra);
        fprintf (ofile, "\n                         | ");
        file_write_tree (ofile, ket);
        fprintf (ofile, " > \n");
}

NODE *copy_tree (NODE *tree)
/* Constructs a copy of a coupling tree */
{
        NODE *current;

        current = (NODE *) malloc (sizeof(NODE));
        current->index = tree->index;
        copy_set (tree->leafs, current->leafs);
        if (tree->left == NULL) {
          current->left = NULL; current->right = NULL;
        }
        else {
          current->left = copy_tree (tree->left);
          current->right = copy_tree (tree->right);
        }
        return current;
}

/*
 *      Building list of delta functions
 *      --------------------------------
 */

int h_build_deltas (NODE *this, NODE *tree, DELTAS *pd)
/* Checks wether the leafset of this is one of the leafsets of tree;
   if so, a new delta is added to *pd */
{
        int i;

	if (equal_sets (this->leafs, tree->leafs)) {
	  if (tree->index != this->index) {
	    pd->n++; i = pd->n;
            pd->d[i][1] = tree->index; pd->d[i][2] = this->index;
            tree->index = this->index;
	  }
	  return TRUE;
	}
	else if (tree->left == NULL) return FALSE;
	else if (h_build_deltas (this, tree->left, pd)) return TRUE;
	else return h_build_deltas (this, tree->right, pd);
}

void build_deltas (NODE *bra, NODE *ket, DELTAS *pd)
/* Builds a list of delta functions for bra en ket, returns in *pd */
{
        if (bra->left != NULL) {
          h_build_deltas (bra, ket, pd);
          build_deltas (bra->left, ket, pd);
          build_deltas (bra->right, ket, pd);
        }
}

/*
 *      Input in the style of NJSYM
 *      ---------------------------
 */

typedef int TRIADS[MAXJ][4];

int find_pos (int nrcoups, TRIADS tr, int val)
/* Finds position of given value in third position in a triad
   of given list; returns zero if not found */
{
        int i;

        for (i=1; i<=nrcoups; i++) if (val == tr[i][3]) return i;
        return 0;
}

NODE *build_tree (int nrcoups, TRIADS tr, int rootval)
/* Builds a tree from given list of triads,
   starting from given root value */
{
        int rootpos;
        NODE *c;

        c = (NODE *) malloc (sizeof(NODE));
        c->index = rootval;
        rootpos = find_pos (nrcoups, tr, rootval);
        if (rootpos == 0) {
          c->left = NULL; c->right = NULL;
          clear_set (c->leafs); add_element (c->index, c->leafs);
        }
        else {
          c->left = build_tree (nrcoups, tr, tr[rootpos][1]);
          c->right = build_tree (nrcoups, tr, tr[rootpos][2]);
          add_sets (c->left->leafs, c->right->leafs, c->leafs);
        }
        return c;
}

void read_njsym (int *pnrjs, NODE **pbra, NODE **pket, DELTAS *pd)
/* Reads a recoupling coefficient in the style of NJSYM */
{
        int nrcoups, dummy, i, rootval;
        TRIADS tr1, tr2;

        printf ("Give nrjs, nrcoups, dummy : ");
        scanf ("%d %d %d", pnrjs, &nrcoups, &dummy);
        printf ("Give list of triads for bra vector : \n");
        for (i=1; i<=nrcoups; i++)
	  scanf ("%d %d %d", &tr1[i][1], &tr1[i][2], &tr1[i][3]);
        printf ("Give list of triads for ket vector : \n");
        for (i=1; i<=nrcoups; i++)
          scanf ("%d %d %d", &tr2[i][1], &tr2[i][2], &tr2[i][3]);
	tr2[nrcoups][3] = tr1[nrcoups][3];
        rootval = tr1[nrcoups][3];
        *pbra = build_tree (nrcoups, tr1, rootval);
        *pket = build_tree (nrcoups, tr2, rootval);
	pd->n = 0;
	build_deltas (*pbra, *pket, pd);
        if (debugging == ON)
          file_write_expr (debugfile, *pbra, *pket);
}

/*
 *	Input in the form of expressions (only interactive input)
 *      ---------------------------------------------------------
 */

void read_char (char *pch)
/* Reads the next non-space character; returns it in *pch */
{
	scanf ("%1s", pch);
}

int read_int (char *pch)
/* Reads the integer starting with digit *pch;
   returns the first non-digit in *pch */
{
	int i = 0;

	while (*pch >= '0' && *pch <= '9') {
	  i = i * 10 + *pch - '0';
	  read_char (pch);
	}
	return i;
}

NODE *read_tree (char *pch)
/* Reads the expression for a coupling tree, starting with character *pch;
   returns the first character not belonging to expression in *pch;
   if an error occurs : returns NULL */
{
	NODE *current;

	current = (NODE *) malloc (sizeof(NODE));
	if (*pch == '(') {
	  read_char (pch);
	  if ((current->left = read_tree (pch)) == NULL) return NULL;
	  if (*pch != ',') return NULL;
	  read_char (pch);
	  if ((current->right = read_tree (pch)) == NULL) return NULL;
	  if (*pch != ')') return NULL;
	  add_sets (current->left->leafs, current->right->leafs,
	            current->leafs);
	  read_char (pch);
          if (AUTONUM == ON)
            current->index = ++currentj;
          else {
  	    if (*pch < '0' || *pch > '9') return NULL;
	    current->index = read_int (pch);
          }
	}
	else if (*pch >= '0' && *pch <= '9') {
	  current->left = NULL;	current->right = NULL;
	  current->index = read_int (pch);
	  clear_set (current->leafs);
	  add_element (current->index, current->leafs);
	}
	else
	  return NULL;
	return current;
}

void read_expr (int *pnrjs, NODE **pbra, NODE **pket, DELTAS *pd)
/* Reads an expression in bra-ket notation and returning coupling trees.
   If AUTONUM is ON, needs number of leaf nodes, given in *pnrjs;
   returns total number of nodes in *pnrjs.
   If an error occurs while reading: clears line and retries */
{
        int stop;
	char ch, s[80];

        if (AUTONUM == ON) {
          printf ("Number of j's given ? ");
          scanf ("%d", pnrjs); gets (s);
        }
        printf ("Give recoupling coefficient ");
        printf ("< state vector | state vector > : \n");
        stop = TRUE;
        do {
          if (stop == FALSE) printf ("  retry ! ");
          if (AUTONUM == ON) currentj = *pnrjs;
          stop = TRUE;
          read_char (&ch);
          if (ch != '<')
            stop = FALSE;
          else {
            read_char (&ch);
    	    *pbra = read_tree (&ch);
            if (ch != '|' || *pbra == NULL)
              stop = FALSE;
            else {
              read_char (&ch);
              *pket = read_tree (&ch);
              if (ch != '>' || *pket == NULL)
                stop = FALSE;
              else {
                (*pket)->index = (*pbra)->index;
		stop = TRUE;
	      }
            }
          }
          gets (s);
        } while (stop == FALSE);
        *pnrjs = currentj - 1;
	pd->n = 0;
	build_deltas (*pbra, *pket, pd);
        if (debugging == ON) file_write_expr (debugfile, *pbra, *pket);
}

/*
 *      Working with formulae
 *      ---------------------
 */

FORMULA *init_formula (int cnrjs)
/* Initialises a formula f, given the number of j's */
{
        int i;
        FORMULA *f;

        f = (FORMULA *) malloc (sizeof(FORMULA));
        f->nrjs = cnrjs; f->nrks = 0; f->nrsixjs = 0;
        for (i=0; i<MAXJ; i++) { f->jsigns[i] = 0; f->jsqrts[i] = 0; }
        for (i=0; i<MAXK; i++) { f->ksigns[i] = 0; f->ksqrts[i] = 0; }
	f->ordered = FALSE;
        return f;
}

void file_write_formula (FILE *ofile, FORMULA *f)
/* Writes the formula to a given file, which has to be opened already */
{
        int i, j;
	
        fprintf (ofile, "FORMULA : \n");
        fprintf (ofile, "nrjs = %d, nrks = %d, nrsixjs = %d\n",
                f->nrjs, f->nrks, f->nrsixjs);
        fprintf (ofile, "jsigns : ");
        for (i=1; i<=f->nrjs; i++)
          fprintf (ofile, "%d ", f->jsigns[i]%4);
        fprintf (ofile, "\n");
        fprintf (ofile, "ksigns : ");
        for (i=1; i<=f->nrks; i++)
          fprintf (ofile, "%d ", f->ksigns[i]%4);
        fprintf (ofile, "\n");
        fprintf (ofile, "jsqrts : ");
        for (i=1; i<=f->nrjs; i++)
          fprintf (ofile, "%d ", f->jsqrts[i]);
        fprintf (ofile, "\n");
        fprintf (ofile, "ksqrts : ");
        for (i=1; i<=f->nrks; i++)
          fprintf (ofile, "%d ", f->ksqrts[i]);
        fprintf (ofile, "\n");
        fprintf (ofile, "sixjs : \n");
        for (i=1; i<=f->nrsixjs; i++) {
          fprintf (ofile, "  %d : ", i);
          for (j=1; j<=6; j++)
            fprintf (ofile, "%d ", f->sixjs[i][j]);
          fprintf (ofile, "\n");
        }
}

void write_formula (FORMULA *f)
/* Writes the formula f to a standard output device */
{
        if (f == NULL) printf ("No proper formula generated\n");
        else file_write_formula (stdout, f);
}

void update_sign (FORMULA *f, int index, int val)
/* In formula f updates sign of given index (>0 : j; <0 : k),
   val == -1 or +1 */
{
        if (index > 0) {
          if (val == +1) f->jsigns[index]++; else f->jsigns[index]--;
        }
        else {
          if (val == +1) f->ksigns[-index]++; else f->ksigns[-index]--;
        }
}

void update_sqrt (FORMULA *f, int index, int val)
/* In formula f updates sqrt of given index (>0 : j; <0 : k),
   val == -1 or +1 */
{
        if (index > 0) {
          if (val == +1) f->jsqrts[index]++; else f->jsqrts[index]--;
        }
        else {
          if (val == +1) f->ksqrts[-index]++; else f->ksqrts[-index]--;
        }
}

void update_sixjs (FORMULA *f, int a, int b, int c, int d, int g, int e)
/* In formula f adds a new sixj with parameters (a,b,c,d,g,e) */
{
        int *lsixj;

        f->nrsixjs++; lsixj = f->sixjs[f->nrsixjs];
        lsixj[1] = a; lsixj[2] = b; lsixj[3] = c;
        lsixj[4] = d; lsixj[5] = g; lsixj[6] = e;
}

void substitute_k (FORMULA *f, int j)
/* Substitute k added in last step by given j in formula f */
{
        int i, lastk = f->nrks;

        f->ksqrts[lastk]--; f->jsqrts[j]++; f->nrks--;
        for (i=1; i<=6; i++)
          if (f->sixjs[f->nrsixjs][i] == -lastk)
            f->sixjs[f->nrsixjs][i] = j;
}

/*
 *      Working with leafsets
 *      ---------------------
 */

int check_leafs (NODE *this, NODE *tree)
/* Checks wether the leafset of this is one of the leafsets of tree;
   if so, clears leafset of expr */
{
	if (equal_sets (this->leafs, tree->leafs)) {
	  clear_set (tree->leafs);
	  if (tree->index != this->index) tree->index = this->index;
	  return TRUE;
	}
	else if (tree->left == NULL) return FALSE;
	else if (check_leafs (this, tree->left)) return TRUE;
	else return check_leafs (this, tree->right);
}

int check_all_leafs (NODE *bra, NODE *ket)
/* Checks if any of the leafsets of ket already occurs in bra;
   if so, appropriate leafset in ket is cleared;
   returns the number of nodes that do NOT occur in ket. */
{
        int tmp = 0;

        if (bra->left != NULL) {
          if (!check_leafs (bra, ket)) tmp++;
          tmp += check_all_leafs (bra->left, ket);
          tmp += check_all_leafs (bra->right, ket);
        }
        return tmp;
}

int test_leafs (NODE *this, NODE *tree, FORMULA *f)
/* Tests wether the leafset of this is one of the leafsets of tree;
   if so, substitutes this k by the appropriate j,
          updates index in this and clears leafset of tree */
{
	if (equal_sets (this->leafs, tree->leafs)) {
	  clear_set (tree->leafs);
          if (this->index < 0) substitute_k (f, tree->index);
	  this->index = tree->index;
	  return TRUE;
	}
	else if (tree->left == NULL) return FALSE;
	else if (test_leafs (this, tree->left, f)) return TRUE;
	else return test_leafs (this, tree->right, f);
}

int finished (NODE *from)
/* Checks if the transformation is finished :
   all leafsets (except bottom ones) empty */
{
	if (from->left == NULL) return TRUE;
	else if (!empty_set (from->leafs)) return FALSE;
	else if (!finished (from->left)) return FALSE;
	else return finished (from->right);
}

void rebuild_leafs (NODE *tree)
/* Rebuilds the leafsets of a given expression */
{
        if (tree->left == NULL) {
          clear_set (tree->leafs);
          add_element (tree->index, tree->leafs);
        }
        else {
          rebuild_leafs (tree->left); rebuild_leafs (tree->right);
          add_sets (tree->left->leafs,tree->right->leafs,tree->leafs);
        }
}

/*
 * 	Transformation of expressions
 *	-----------------------------
 */

void exchange (NODE *tree, FORMULA *f, int up)
/* Exchanges left and right subtrees of tree;
   a+b-c if up==1; -a-b+c if up==-1 */
{
	NODE *tmp;
        int a, b, c;

        a = tree->left->index; b = tree->right->index; c = tree->index;
        if (up == 1) {
          update_sign (f,a,1); update_sign (f,b,1);
          update_sign (f,c,-1);
        }
        else {
          update_sign (f,a,-1); update_sign (f,b,-1);
          update_sign (f,c,1);
        }
	tmp = tree->left; tree->left = tree->right;
        tree->right = tmp;
}

void flop_left (NODE *tree, FORMULA *f)
/* Flops a tree ((a,b)d,c)g to (a,(b,c)e)g  */
{
	NODE *nleft, *nright;
        int a, b, c, d, e, g;

        f->nrks++;
        a = tree->left->left->index; b = tree->left->right->index;
        c = tree->right->index; d = tree->left->index;
        e = - f->nrks; g = tree->index;
        update_sign (f,a,1); update_sign (f,b,1);
        update_sign (f,c,1); update_sign (f,g,1);
        update_sqrt (f,d,1); update_sqrt (f,e,1);
        update_sixjs (f, a, b, d, c, g, e);
	nleft  = tree->left->left; nright = tree->left;
	nright->left = tree->left->right; nright->right = tree->right;
	add_sets (nright->left->leafs,nright->right->leafs,
                  nright->leafs);
	tree->left = nleft; tree->right = nright;
        tree->right->index = e;
}

void unflop_left (NODE *tree, FORMULA *f)
/* Unflops a tree (a,(b,c)e)g back to ((a,b)d,c)g */
{
	NODE *nleft, *nright;
        int a, b, c, d, e, g;

        f->nrks--;
        a = tree->left->index; b = tree->right->left->index;
        c = tree->right->right->index; d = f->sixjs[f->nrsixjs][3];
        e = tree->right->index; g = tree->index;
        update_sign (f,a,-1); update_sign (f,b,-1);
        update_sign (f,c,-1); update_sign (f,g,-1);
        update_sqrt (f,d,-1); update_sqrt (f,e,-1);
        f->nrsixjs--;
	nright = tree->right->right; nleft = tree->right;
	nleft->right = tree->right->left; nleft->left = tree->left;
	add_sets (nleft->right->leafs, nleft->left->leafs,
                  nleft->leafs);
	tree->right = nright; tree->left = nleft;
        tree->left->index = d;
}

int try_allflopseq (int flen, NODE *bra, NODE *ket, FORMULA *f);

int try_flop_left (int flen, NODE *frombra, NODE *ket, FORMULA *f)
/* Tries a flop left on frombra and checks if it can be the first flop
   in a sequence of length flen creating a new node of ket. */
{
	if (frombra->left == NULL) return FALSE;
	if (frombra->left->left == NULL) return FALSE;
	flop_left (frombra, f);
        if (flen == 1) {
	  if (test_leafs (frombra->right, ket, f)) return TRUE;
        }
        else {
          if (try_allflopseq (flen-1, frombra, ket, f)) return TRUE;
        }
	unflop_left (frombra, f);
	return FALSE;
}

int try_4flops (int flen, NODE *frombra, NODE *ket, FORMULA *f)
/* Tries to find a flop sequence of length flen,
   starting with one of the 4 possible flops on frombra,
   creating a node of ket. */
{
        if (!left_first) exchange (frombra, f, -1);
	if (frombra->left != NULL) {
	  if (frombra->left->left != NULL) {
  	    if (try_flop_left (flen, frombra, ket, f)) return TRUE;
	    exchange (frombra->left, f, 1);
	    if (try_flop_left (flen, frombra, ket, f)) return TRUE;
	    exchange (frombra->left, f, -1);
	  }
        }
	exchange (frombra, f, 1);
	if (frombra->left != NULL) {
	  if (frombra->left->left != NULL) {
  	    if (try_flop_left (flen, frombra, ket, f)) return TRUE;
	    exchange (frombra->left, f, 1);
	    if (try_flop_left (flen, frombra, ket, f)) return TRUE;
	    exchange (frombra->left, f, -1);
          }
	}
	if (left_first) exchange (frombra, f, -1);
	return FALSE;
}

int try_allflopseq (int flen, NODE *bra, NODE *ket, FORMULA *f)
/* Tries to find a flop sequence of length flen (breadth first search),
   starting from bra, creating a node of ket. */
{
        QUEUE *q;
        NODE *c;
        int found;

        q = init_queue (); c = bra; found = FALSE;
        while (!found && c != NULL) {
	  if (try_4flops (flen, c, ket, f)) {
            found = TRUE;
	  }
	  else {
            if (c->left != NULL) {
              if (left_first) {
                add_to_queue (q, c->left); add_to_queue (q, c->right);
              }
              else {
                add_to_queue (q, c->right); add_to_queue (q, c->left);
              }
            }
            get_from_queue (q, &c);
          }
        }
        free_queue (q);
        return found;
}

void swap_final (NODE *frombra, NODE *fromket, FORMULA *f)
/* Final scan of expression, to swap left and right nodes if necessary;
   assumes bra expression has been transformed properly into ket one */
{
        if (frombra->left != NULL) {
          if (frombra->left->index != fromket->left->index)
            exchange (frombra, f, -1);
          swap_final (frombra->left, fromket->left, f);
          swap_final (frombra->right, fromket->right, f);
        }
}

FORMULA *h_generate_formula (int nrjs, NODE *bra, NODE *ket)
/* Transforms expression bra into ket, building formula f;
   returns NULL if transformation is not successful. */
{
        int flen, possible, n;
        NODE *cpbra;
        FORMULA *f;

        if (SEARCH == SHALLOW) cpbra = bra;
        else cpbra = copy_tree (bra);
        ket->index = cpbra->index;
        n = check_all_leafs (cpbra, ket);
        f = init_formula (nrjs);
        if (debugging == ON) {
	  fprintf (debugfile, "\n---------------------------\n");
  	  fprintf (debugfile, "\nStarting from : ");
          file_write_tree (debugfile, cpbra);
          fprintf (debugfile, "\n");
        }
	do {
	  flen = 0;
       	  do {
	    flen++;
	    possible = try_allflopseq (flen, cpbra, ket, f);
	  } while (!possible && flen < MAXFLEN);
	  if (possible) {
            n--;
            if (debugging == ON) {
	      fprintf (debugfile, "%d : ", flen);
              file_write_tree (debugfile, cpbra);
              fprintf (debugfile, "\n");
	    }
          }
	} while (n > 0 && possible);
        if (n > 0) return NULL;
        swap_final (cpbra, ket, f);
        if (debugging == ON) {
          fprintf (debugfile, "After final sweep : ");
          file_write_tree (debugfile, cpbra);
          fprintf (debugfile, "\n\n");
          file_write_formula (debugfile, f);
        }
        if (SEARCH != SHALLOW) rebuild_leafs (ket);
        return f;
}

FORMULA *generate_formula (int nrjs, NODE *bra, NODE *ket)
/* Transforms expression bra into ket, building formula f;
   returns NULL if transformation is not successful.
   Possible searchlevels : SHALLOW, DEEP, VERYDEEP. */
{
        FORMULA *f, *hf;

        if (!equal_sets (bra->leafs, ket->leafs)) return NULL;
        /* bra -> ket , right first */
        left_first = FALSE;
        f = h_generate_formula (nrjs, bra, ket);
        if (SEARCH == SHALLOW || f->nrks < 2) return f;
        if (f == NULL) return NULL;
        /* bra -> ket , left first */
        left_first = TRUE;
        hf = h_generate_formula (nrjs, bra, ket);
        if (hf->nrks < f->nrks) f = hf;
        if (SEARCH == DEEP || f->nrks < 2) return f;
        /* ket -> bra , right first */
        left_first = FALSE;
        hf = h_generate_formula (nrjs, ket, bra);
        if (hf->nrks < f->nrks) f = hf;
        if (f->nrks < 2) return f;
        /* ket -> bra , left first */
        left_first = TRUE;
        hf = h_generate_formula (nrjs, ket, bra);
        if (hf->nrks < f->nrks) f = hf;
        return f;
}
