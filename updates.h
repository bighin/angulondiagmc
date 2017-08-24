#ifndef __UPDATES_H__
#define __UPDATES_H__

#include <stdbool.h>

#include "diagrams.h"

double calculate_free_propagator_weight(struct diagram_t *dgr,struct g0_t *g0);
double calculate_arc_weight(struct diagram_t *dgr,struct arc_t *arc);
double calculate_vertex_weight(struct diagram_t *dgr,int index);

void diagram_remove_midpoint(struct diagram_t *dgr,int c);
void diagram_add_midpoint(struct diagram_t *dgr,int c,double tau,struct arc_t *phononline);
bool diagram_remove_phonon_line(struct diagram_t *dgr,int c);
void diagram_add_phonon_line(struct diagram_t *dgr,double tau1,double tau2,double k,int lambda,int mu);
void diagram_update_length(struct diagram_t *dgr,double newendtau);

bool diagram_add_worm(struct diagram_t *dgr,int target1,int target2,int deltalambda);
bool diagram_remove_worm(struct diagram_t *dgr,int index);

bool recouple(struct diagram_t *dgr,int lo,int hi);

struct free_propagators_ctx_t
{
	int nr_free_propagators;
	int lo,hi;
	int *js;
};

void save_free_propagators(struct diagram_t *dgr,struct free_propagators_ctx_t *fpc,int lo,int hi);
void unload_free_propagators_ctx(struct free_propagators_ctx_t *fpc);
void restore_free_propagators(struct diagram_t *dgr,struct free_propagators_ctx_t *fpc);

#endif //__UPDATES_H__
