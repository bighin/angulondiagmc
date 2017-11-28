#ifndef __MC_H__
#define __MC_H__

#include <stdbool.h>

#include "diagrams.h"
#include "config.h"

bool propagators_are_physical(struct diagram_t *dgr);
bool angular_momentum_is_conserved(struct diagram_t *dgr);

int update_length(struct diagram_t *dgr,struct configuration_t *cfg);
int update_add_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg);
int update_remove_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg);
int update_add_worm(struct diagram_t *dgr,struct configuration_t *cfg);
int update_remove_worm(struct diagram_t *dgr,struct configuration_t *cfg);

void load_config_defaults(struct configuration_t *config);

struct mc_output_data_t
{
	double avglength;
	double avgorder;
};

int do_diagmc(struct configuration_t *config,struct mc_output_data_t *output);

#endif //__MC_H__
