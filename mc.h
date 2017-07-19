#ifndef __MC_H__
#define __MC_H__

#include <stdbool.h>

#include "diagrams.h"

struct configuration_t
{
	/* "general" section */

	char *prefix;
	bool seedrng;
	bool progressbar;

	/* "parameters" section */

	int j;
	int m;
	double endtau;
	double chempot;
	double maxtau;

	/* "potential" section */
	
	double c0;
	double c1;
	double c2;
	double omega0;
	double omega1;
	double omega2;

	/* "sampling" section */

	int iterations;
	int bins;
	double width;
};

#define DIAGRAM_NR_UPDATES	(5)

#define UPDATE_UNPHYSICAL	(0)
#define UPDATE_REJECTED		(1)
#define UPDATE_ACCEPTED		(2)
#define UPDATE_ERROR		(3)

int update_length(struct diagram_t *dgr,struct configuration_t *cfg);
int update_add_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg);
int update_remove_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg);
int update_add_worm(struct diagram_t *dgr,struct configuration_t *cfg);
int update_remove_worm(struct diagram_t *dgr,struct configuration_t *cfg);

int do_diagmc(char *configfile);

void load_config_defaults(struct configuration_t *config);

#endif //__MC_H__
