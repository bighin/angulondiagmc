#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

#include "diagrams.h"
#include "updates.h"
#include "stat.h"
#include "inih/ini.h"

/*
	A pointer to a GSL random number generator object

	This is NOT thread safe, and must modified if parallel code is wanted!
*/

gsl_rng *rng_ctx;

void stresstest(void)
{
	struct diagram_t *dgr;
	int x,y,z;
	int parx,pary,parz;

	if((rng_ctx=gsl_rng_alloc(gsl_rng_mt19937))==NULL)
	{
		printf("Couldn't initialize the random number generator\n");
		return;
	}

	dgr=init_diagram(1.0f,2,1,5.0f);

	parx=64;
	pary=128;
	parz=46;

	assert(parz<pary);

	for(x=0;x<parx;x++)
	{
		for(y=0;y<pary;y++)
		{
			double lo,hi,k;
			int lambda,mu;

			lo=gsl_rng_uniform(rng_ctx);
			hi=gsl_rng_uniform(rng_ctx);
			k=gsl_rng_uniform(rng_ctx);

			lambda=gsl_rng_uniform_int(rng_ctx,8);
			mu=gsl_rng_uniform_int(rng_ctx,lambda+1);

			if(gsl_rng_uniform(rng_ctx)<0.5f)
				mu=-mu;

			if(lo>hi)
			{
				double tmp;
		
				tmp=hi;
				hi=lo;
				lo=tmp;
			}

			printf("<+ (%d) (w=%f)>\n",get_nr_phonons(dgr),diagram_weight(dgr));

			diagram_add_phonon_line(dgr,lo,hi,k,lambda,mu);
			diagram_check_consistency(dgr);
		}

		for(z=0;z<parz;z++)
		{
			if(get_nr_phonons(dgr)>0)
			{
				int target=gsl_rng_uniform_int(rng_ctx,get_nr_phonons(dgr));
				
				diagram_remove_phonon_line(dgr,target);
				diagram_check_consistency(dgr);
			
				printf("<- (%d)>\n",get_nr_phonons(dgr));
			}
		}
	}

	while(get_nr_vertices(dgr)>0)
	{
		if(get_nr_vertices(dgr)<=19)
			print_diagram(dgr);

		printf("<- (%d)>\n",get_nr_phonons(dgr));
		diagram_remove_phonon_line(dgr,gsl_rng_uniform_int(rng_ctx,get_nr_phonons(dgr)));
	
	}
}

struct configuration_t
{
	int iterations;
	char *prefix;
};

static int configuration_handler(void *user,const char *section,const char *name,const char *value)
{
	struct configuration_t *pconfig=(struct configuration_t *)(user);

#define MATCH(s,n) ((strcmp(section,s)==0)&&(strcmp(name,n)==0))

	if(MATCH("diagmc","iterations"))
	{
		pconfig->iterations=atoi(value);
	}
	else if(MATCH("diagmc","prefix"))
	{
		pconfig->prefix=strdup(value);
	}
	else
	{
		/* Unknown section/name, error */
		return 0;
	}

	return 1;
}

int do_diagmc(char *configfile)
{
	struct diagram_t *dgr;
	struct configuration_t config;
	struct samples_t *samples;

	int c;

        if(ini_parse(configfile,configuration_handler,&config)<0)
	{
		fprintf(stderr,"Couldn't read or parse '%s'\n",configfile);
		exit(0);
	}

	fprintf(stderr,"Loaded '%s'\n",configfile);
	fprintf(stderr,"Performing %d iterations\n",config.iterations);
	fprintf(stderr,"Writing results to prefix '%s'\n",config.prefix);

	if((rng_ctx=gsl_rng_alloc(gsl_rng_mt19937))==NULL)
	{
		printf("Couldn't initialize the random number generator\n");
		exit(0);
	}

	dgr=init_diagram(1.0f,2,1,5.0f);
	samples=samples_init();

	for(c=0;c<config.iterations;c++)
	{
		samples_add_entry(samples,diagram_weight(dgr));
	}

	printf("%f %f\n",samples_get_average(samples),samples_get_variance(samples));

	samples_fini(samples);
	fini_diagram(dgr);

	return 0;
}

int main(int argc,char *argv[])
{
	if(argc!=2)
	{
		printf("Usage: %s <inifile>\n",argv[0]);
		return 0;
	}
	
	do_diagmc(argv[1]);
	
	return 0;
}
