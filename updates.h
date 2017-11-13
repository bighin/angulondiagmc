#ifndef __UPDATES_H__
#define __UPDATES_H__

#include <stdbool.h>

#include "diagrams.h"

void diagram_remove_midpoint(struct diagram_t *dgr,int c);
void diagram_add_midpoint(struct diagram_t *dgr,int c,double tau,struct arc_t *phononline);
void diagram_remove_phonon_line(struct diagram_t *dgr,int c);
void diagram_add_phonon_line(struct diagram_t *dgr,double tau1,double tau2,int lambda,int mu);
void diagram_update_length(struct diagram_t *dgr,double newendtau);

bool diagram_add_worm(struct diagram_t *dgr,int target1,int target2,int deltalambda);
bool diagram_remove_worm(struct diagram_t *dgr,int index);

bool recouple(struct diagram_t *dgr,int lo,int hi);

#endif //__UPDATES_H__
