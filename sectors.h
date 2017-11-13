#ifndef __SECTORS_H__
#define __SECTORS_H__

#include <stdbool.h>

#include "diagrams.h"

bool check_triangle_condition(struct diagram_t *dgr,struct vertex_t *thisvertex);
bool check_parity(struct diagram_t *dgr,struct vertex_t *thisvertex);

bool is_in_P(struct diagram_t *dgr);
bool is_in_E(struct diagram_t *dgr);

#endif //__SECTORS_H__
