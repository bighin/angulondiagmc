#ifndef __GRAPHS_H__
#define __GRAPHS_H__

#include "njsummat/njformul.h"
#include "njsummat/njsummat.h"
#include "diagrams.h"

void reset_formula(FORMULA *f);
void formula_to_wolfram(FORMULA f);

struct graph_t
{
	/*
		The number of vertices, lines and arcs.

		The information is redundant, since at every time
		one must have:
	
		nr_lines == nr_vertices + 1
		nr_vertices == nr_arcs * 2
	*/
	
	int nr_vertices,nr_lines,nr_arcs;

#define MAX_LINES		(1024)
#define MAX_ARCS		(512)
#define MAX_DELTAS		(512)

	/*
		lines[i] contains the label number of the momentum
		circulating on the i-th line. For instance if the
		second line has the momentum labelled j7, then one
		would have to set
	
		lines[1]=7
	
		Note that the momentum j0 is not used.

		arcs[j] is an array of three ints, indicating, respectively,
		the starting vertex, the ending vertex and the momentum
		label circulating inside the vertex.
	*/

	int lines[MAX_LINES];
	int arcs[MAX_ARCS][3];

	/*
		js holds the concrete values of angular momenta
	*/

	int js[MAX_LINES];
	int maxj,maxk;

	int delta[MAX_DELTAS][2];
	int nr_deltas;

	/*
		This structure will hold the resulting formula, in the format used
		by NJSUMMAT
	*/

	FORMULA f;
};

void emit_6j(struct graph_t *gt,int a,int b,int c,int d,int e,int f);
void emit_jhat(struct graph_t *gt,int a,int power);
void emit_delta(struct graph_t *gt,int a,int b);
void emit_summation(struct graph_t *gt,int newj);
void emit_phase(struct graph_t *gt,int a,int mul);

int find_arc_with_vertex(struct graph_t *gt,int vertex);
void show_graph(struct graph_t *gt);

#define SAVE_LEFT	(1)
#define SAVE_RIGHT	(2)

void remove_vertex(struct graph_t *gt,int position,int mode);
void remove_arc(struct graph_t *gt,int index);

bool has_nspanning_arc(struct graph_t *gt,int vertex,int span);
bool has_bubble(struct graph_t *gt,int vertex);
bool has_triangle(struct graph_t *gt,int vertex);
bool has_square(struct graph_t *gt,int vertex);

void remove_bubble(struct graph_t *gt,int vertex);
void remove_triangle(struct graph_t *gt,int vertex);
void remove_square(struct graph_t *gt,int vertex);
void exchange_lines(struct graph_t *gt,int v1,int v2);

void k_to_j(struct graph_t *gt,int kindex,int jindex);
void k_to_k(struct graph_t *gt,int kindex1,int kindex2);

double evaluate_graph(struct graph_t *gt,bool debugswap);
void diagram_to_graph(struct diagram_t *dgr,struct graph_t *gt);

struct hashentry_t
{
	int valid;

	struct graph_t gt;
	double value;
};

#define HASHTABLE_ENTRIES	(4096)

#define HASHTABLE_VALID_ENTRY	(101)
#define HASHTABLE_INVALID_ENTRY	(102)
	
void hashtable_init(void);
bool hashtable_probe(struct graph_t *gt,double *value,int *hashindex);
void hashtable_insert(struct graph_t *gt,double value,int hashindex);
void hasthtable_show_stats(void);

double diagram_m_weight(struct diagram_t *dgr,bool use_hashtable);
double diagram_m_weight_reference(struct diagram_t *dgr);

void init_njsummat(void);

#endif //__GRAPHS_H__
