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
	printf("       %s --test\n",argv0);

	exit(0);
}

int main(int argc,char *argv[])
{
	init_njsummat();

	if(argc!=2)
		usage(argv[0]);

	if(strcmp(argv[1],"--test")==0)
	{
		return test_graphical_machinery();

		//test_graph();
		//return 0;

		//stresstest();
		//return 0;
	}
	
	do_diagmc(argv[1]);
	
	return 0;
}
