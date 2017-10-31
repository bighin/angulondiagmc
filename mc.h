#ifndef __MC_H__
#define __MC_H__

#include <stdbool.h>

#include "diagrams.h"
#include "config.h"

int update_length(struct diagram_t *dgr,struct configuration_t *cfg);
int update_add_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg);
int update_remove_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg);
int update_add_worm(struct diagram_t *dgr,struct configuration_t *cfg);
int update_remove_worm(struct diagram_t *dgr,struct configuration_t *cfg);

void load_config_defaults(struct configuration_t *config);

int do_diagmc(struct configuration_t *config);

#endif //__MC_H__
