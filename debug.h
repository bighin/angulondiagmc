#ifndef __DEBUG_H__
#define __DEBUG_H__

#include "diagrams.h"

void debug_propagators(struct diagram_t *dgr);
void debug_vertices(struct diagram_t *dgr);
void debug_vertices_ext(struct diagram_t *dgr);
void debug_weight(struct diagram_t *dgr);

int stresstest(void);

#endif //__DEBUG_H__
