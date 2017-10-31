#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdbool.h>

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
	double timelimit;
	int bins;
	double width;

	/* "parallel" section */
	
	bool parallel;
	int nthreads;
	
	/* The name of the file the configuration has been loaded from */
	
	char *ininame;
};

int configuration_handler(void *user,const char *section,const char *name,const char *value);
void load_config_defaults(struct configuration_t *config);
bool load_configuration(char *configfile,struct configuration_t *config);

#endif //__CONFIG_H__
