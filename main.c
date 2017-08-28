#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>

#include "diagrams.h"
#include "updates.h"
#include "stat.h"
#include "debug.h"
#include "aux.h"
#include "mc.h"
#include "graphs.h"
#include "tests.h"

void usage(char *argv0)
{
	printf("Usage: %s <inifile>\n",argv0);
	printf("       %s --testgraphs\n",argv0);
	printf("       %s --stresstest\n",argv0);


	exit(0);
}

int main(int argc,char *argv[])
{
	int c;

	init_njsummat();

	if(argc<2)
		usage(argv[0]);

	for(c=1;c<argc;c++)
	{
		if(strcmp(argv[c],"--testgraphs")==0)
		{
			test_graphical_machinery();
			continue;
		}

		if(strcmp(argv[c],"--stresstest")==0)
		{
			stresstest();
			continue;
		}

		do_diagmc(argv[c]);
	}

	return 0;
}
