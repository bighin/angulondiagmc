#ifndef __PHYSICS_H__
#define __PHYSICS_H__

#include <stdbool.h>
#include "diagrams.h"

double calculate_arc_weight(struct diagram_t *dgr,struct arc_t *arc);
double calculate_free_propagator_weight(struct diagram_t *dgr,struct g0_t *g0);
double calculate_vertex_weight(struct diagram_t *dgr,int index);

bool propagators_are_allowed(struct diagram_t *dgr);
bool configuration_is_physical(struct diagram_t *dgr);
bool angular_momentum_is_conserved(struct diagram_t *dgr);

double calculate_propagators_and_vertices(struct diagram_t *dgr,int startmidpoint,int endmidpoint);

int deltaj(struct diagram_t *dgr,int vertex);
void change_deltaj(struct diagram_t *dgr,int index,int newdeltaj);

#endif //__PHYSICS_H__
