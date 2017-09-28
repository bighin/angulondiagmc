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
	bool animate;
	bool liveplot;
	bool use_hashtable;

	/* "parameters" section */

	int j;
	double endtau;
	double chempot;
	double maxtau;
	int maxorder;

	/* "potential" section */
	
	double n;

	/* "sampling" section */

	int iterations;
	int bins;
	double width;

	/* "parallel" section */
	
	bool parallel;
	int nthreads;
};

int update_length(struct diagram_t *dgr,struct configuration_t *cfg);
int update_add_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg);
int update_remove_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg);
int update_add_worm(struct diagram_t *dgr,struct configuration_t *cfg);
int update_remove_worm(struct diagram_t *dgr,struct configuration_t *cfg);

void load_config_defaults(struct configuration_t *config);

int do_diagmc(char *configfile);

#endif //__MC_H__
