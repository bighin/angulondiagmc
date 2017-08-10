#include <stdio.h>
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
#include "njgraf/njgraf.h"

int main(int argc,char *argv[])
{
	//stresstest();
	//return 0;
	
	if(argc!=2)
	{
		printf("Usage: %s <inifile>\n",argv[0]);
		return 0;
	}
	
	do_diagmc(argv[1]);
	
	return 0;
}
