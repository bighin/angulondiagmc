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
#include "phonon.h"
#include "selfenergies.h"
#include "config.h"

#define RETURN_AVGLENGTH	(22)
#define RETURN_AVGORDER		(23)

double try_diagmc_and_get_length(struct configuration_t *config,double chempot,short retval)
{
	struct mc_output_data_t output;

	config->timelimit=5*60.0f;
	config->chempot=chempot;

	printf("Trying DiagMC run for %f seconds, with chemical potential %f\n",config->timelimit,config->chempot);

	do_diagmc(config,&output);
	
	printf("Result: (chempot, avglength) = (%f %f)\n\n",config->chempot,output.avglength);
	
	if(retval==RETURN_AVGLENGTH)
		return output.avglength;

	if(retval==RETURN_AVGORDER)
		return output.avgorder;

	printf("Error in try_diagmc_and_get_length()!\n");
	exit(0);
}

void usage(char *argv0)
{
	printf("Usage: %s <inifile>\n",argv0);
	printf("       %s --tuning <inifile>\n",argv0);
	printf("       %s --testgraphs\n",argv0);

	exit(0);
}

int main(int argc,char *argv[])
{
	int c;
	bool first=true;

	init_njsummat();
	hashtable_init();

	if(argc<2)
		usage(argv[0]);

	for(c=1;c<argc;c++)
	{
		struct configuration_t config;
		struct mc_output_data_t output;

		if(strcmp(argv[c],"--testgraphs")==0)
		{
			test_graphical_machinery();
			continue;
		}

		if(strcmp(argv[c],"--testphonons")==0)
		{
			struct phonon_ctx_t *ctx=init_phonon_ctx(120.0f,128*1024,exp(-10),true);

			test_phonon_ctx(ctx);
			continue;
		}

		if(strcmp(argv[c],"--stressphonons")==0)
		{
			struct phonon_ctx_t *ctx=init_phonon_ctx(120.0f,128*1024,exp(-10),true);
			gsl_rng *rctx=gsl_rng_alloc(gsl_rng_mt19937);

			for(int d=0;d<1000000000;d++)
			{
				double res0,res1;

				res0=phonon_dist(rctx,ctx,0);
				res1=phonon_dist(rctx,ctx,1);
			
				if((res0<0)||(res1<0))
				{
					printf("Error detected!");
					exit(0);
				}
			}

			continue;
		}

		if(strcmp(argv[c],"--tuning")==0)
		{
			double targetlength;
			double lo,hi,mid;
			int d,j;
			
			if((c+1)>=argc)
				usage(argv[0]);

			c++;

			if(load_configuration(argv[c],&config)==false)
				exit(0);

			fprintf(stderr,"Diagrammatic Monte Carlo for the angulon (automatic chemical potential tuning)\n");

			targetlength=7.5f;

			j=config.j;

			lo=-15.0f+j*(j+1);
			hi=0.0f+j*(j+1);

			for(d=0;d<=10;d++)
			{
				double fmid;
				
				mid=(lo+hi)/2.0f;
				load_configuration(argv[c],&config);
				fmid=try_diagmc_and_get_length(&config,mid,RETURN_AVGLENGTH);

				if(fmid<targetlength)
				{
					lo=mid;
				}
				else
				{
					hi=mid;
				}
			}

			printf("Performing final run with full time control\n");

			load_configuration(argv[c],&config);
			config.chempot=(lo+hi)/2.0f;
			do_diagmc(&config,&output);

			printf("Final result: (chempot, avglength) = (%f %f)\n\n",config.chempot,output.avglength);

			continue;
		}

		if(strcmp(argv[c],"--alttuning")==0)
		{
			double targetorder;
			double lo,hi,mid;
			int d,j;
			
			if((c+1)>=argc)
				usage(argv[0]);

			c++;

			if(load_configuration(argv[c],&config)==false)
				exit(0);

			fprintf(stderr,"Diagrammatic Monte Carlo for the angulon (automatic chemical potential tuning, alternate version)\n");

			targetorder=30.0f;

			j=config.j;

			lo=-15.0f+j*(j+1);
			hi=0.0f+j*(j+1);

			for(d=0;d<=10;d++)
			{
				double fmid;
				
				mid=(lo+hi)/2.0f;
				load_configuration(argv[c],&config);
				fmid=try_diagmc_and_get_length(&config,mid,RETURN_AVGORDER);

				if(fmid<targetorder)
				{
					lo=mid;
				}
				else
				{
					hi=mid;
				}
			}

			printf("Performing final run with full time control\n");

			load_configuration(argv[c],&config);
			config.chempot=(lo+hi)/2.0f;
			do_diagmc(&config,&output);

			printf("Final result: (chempot, avglength) = (%f %f)\n\n",config.chempot,output.avglength);

			continue;
		}

		if(first==true)
			fprintf(stderr,"Diagrammatic Monte Carlo for the angulon\n");

		if(load_configuration(argv[c],&config)==false)
			continue;

		do_diagmc(&config,&output);
		first=false;
		
		if((c+1)!=argc)
			printf("\n");
	}

	{
		extern int hashtable_lookups;
		
		if(hashtable_lookups>0)
			hasthtable_show_stats();
	}

	return 0;
}
