#ifndef __UPDATES_H__
#define __UPDATES_H__

#include <stdbool.h>

#include "diagrams.h"

void diagram_remove_midpoint(struct diagram_t *dgr,int c);
void diagram_add_midpoint(struct diagram_t *dgr,int c,double tau,struct arc_t *phononline);
void diagram_remove_phonon_line(struct diagram_t *dgr,int c);
void diagram_add_phonon_line(struct diagram_t *dgr,double tau1,double tau2,int lambda,int mu,double *weight);
void diagram_update_length(struct diagram_t *dgr,double newendtau);

struct free_propagators_ctx_t
{
	int nr_free_propagators;
	int lo,hi;
	int *js;
};

void save_free_propagators(struct diagram_t *dgr,struct free_propagators_ctx_t *fpc,int lo,int hi);
void release_free_propagators_ctx(struct free_propagators_ctx_t *fpc);
void restore_free_propagators(struct diagram_t *dgr,struct free_propagators_ctx_t *fpc);

bool recouple(struct diagram_t *dgr,int lo,int hi);

#endif //__UPDATES_H__
