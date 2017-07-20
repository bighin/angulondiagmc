#ifndef __UPDATES_H__
#define __UPDATES_H__

#include <stdbool.h>

#include "diagrams.h"

double calculate_free_propagator_weight(struct diagram_t *dgr,struct g0_t *g0);
double calculate_arc_weight(struct diagram_t *dgr,struct arc_t *arc);

void diagram_remove_midpoint(struct diagram_t *dgr,int c);
void diagram_add_midpoint(struct diagram_t *dgr,int c,double tau,struct arc_t *phononline);
void diagram_remove_phonon_line(struct diagram_t *dgr,int c);
void diagram_add_phonon_line(struct diagram_t *dgr,double tau1,double tau2,double k,int lambda,int mu);
void diagram_update_length(struct diagram_t *dgr,double newendtau);

bool diagram_add_worm(struct diagram_t *dgr,int target1,int target2,int deltalambda);
bool diagram_remove_worm(struct diagram_t *dgr,int index);

struct free_propagators_ctx_t
{
	int nr_free_propagators;
	int *js;
	int *ms;
};

void save_free_propagators(struct diagram_t *dgr,struct free_propagators_ctx_t *fpc);
void restore_free_propagators(struct diagram_t *dgr,struct free_propagators_ctx_t *fpc);

bool check_couplings_ms(struct diagram_t *dgr);
bool recouple_ms(struct diagram_t *dgr);
void recouple_ms_and_assert(struct diagram_t *dgr);

#endif //__UPDATES_H__
