#ifndef __UPDATES_H__
#define __UPDATES_H__

#include "diagrams.h"

void diagram_remove_midpoint(struct diagram_t *dgr,int c);
void diagram_add_midpoint(struct diagram_t *dgr,int c,double tau,struct arc_t *phononline);
void diagram_remove_phonon_line(struct diagram_t *dgr,int c);
void diagram_add_phonon_line(struct diagram_t *dgr,double tau1,double tau2,double k,int lambda,int mu);
void diagram_update_length(struct diagram_t *dgr,double newmaxtau);

#endif //__UPDATES_H__
