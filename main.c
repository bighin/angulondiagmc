#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

#include "diagrams.h"
#include "updates.h"
#include "stat.h"

#include "inih/ini.h"
#include "libprogressbar/progressbar.h"

/*
	A pointer to a GSL random number generator object

	This is NOT thread safe, and must modified if parallel code is wanted!
*/

gsl_rng *rng_ctx;

void seed_rng(gsl_rng *rng)
{
	FILE *dev;
        unsigned long seed;

	if(!(dev=fopen("/dev/random","r")))
	{
		fread(&seed,sizeof(unsigned long),1,dev);
		fclose(dev);
	}

	gsl_rng_set(rng,seed);
}

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

	/* "sampling" section */

	int iterations;
	int bins;
	double width;
};

static int configuration_handler(void *user,const char *section,const char *name,const char *value)
{
	struct configuration_t *pconfig=(struct configuration_t *)(user);

#define MATCH(s,n) ((strcmp(section,s)==0)&&(strcmp(name,n)==0))

	if(MATCH("general","prefix"))
	{
		pconfig->prefix=strdup(value);
	}
	else if(MATCH("general","seedrng"))
	{
		if(!strcmp(value,"true"))
			pconfig->seedrng=true;
		else if(!strcmp(value,"false"))
			pconfig->seedrng=false;
		else
			return 0;
	}
	else if(MATCH("general","progressbar"))
	{
		if(!strcmp(value,"true"))
			pconfig->progressbar=true;
		else if(!strcmp(value,"false"))
			pconfig->progressbar=false;
		else
			return 0;
	}
	else if(MATCH("parameters","j"))
	{
		pconfig->j=atoi(value);
	}
	else if(MATCH("parameters","m"))
	{
		pconfig->m=atoi(value);
	}
	else if(MATCH("parameters","endtau"))
	{
		pconfig->endtau=atof(value);
	}
	else if(MATCH("parameters","chempot"))
	{
		pconfig->chempot=atof(value);
	}
	else if(MATCH("parameters","maxtau"))
	{
		pconfig->maxtau=atof(value);
	}
	else if(MATCH("sampling","iterations"))
	{
		pconfig->iterations=atoi(value);
	}
	else if(MATCH("sampling","bins"))
	{
		pconfig->bins=atoi(value);
	}
	else if(MATCH("sampling","width"))
	{	
		pconfig->width=atof(value);
	}
	else
	{
		/* Unknown section/name, error */
		return 0;
	}

	return 1;
}

/*
	Returns a random number from a doubly truncated exponential distribution,
	i.e. a distribution with probability density

	\frac{lambda \exp(-\lambda \tau)}{\exp(-\lambda \tau_1) - \exp(-\lambda \tau_2)}

	See DTED.nb for the derivation.
*/

double doubly_truncated_exp_dist(gsl_rng *rctx,double lambda,double tau1,double tau2)
{
	double x=gsl_rng_uniform(rctx);
	double delta=-x+x*exp((tau1-tau2)*lambda);

	return (tau1-log1p(delta))/lambda;
}

int do_diagmc(char *configfile)
{
	struct diagram_t *dgr;
	struct configuration_t config;
	struct histogram_t *ht;
	
	FILE *out;
	char fname[1024];
	progressbar *progress;

	int c;

#warning FIXME Load sensible defaults in case they are not in the .ini file

        if(ini_parse(configfile,configuration_handler,&config)<0)
	{
		fprintf(stderr,"Couldn't read or parse '%s'\n",configfile);
		exit(0);
	}

	fprintf(stderr,"Diagrammatic Monte Carlo for the angulon\n");
	fprintf(stderr,"Loaded '%s'\n",configfile);
	fprintf(stderr,"Performing %d iterations\n",config.iterations);

	snprintf(fname,1024,"%s.dat",config.prefix);
	
	if(!(out=fopen(fname,"w+")))
	{
		fprintf(stderr,"Error: couldn't open %s for writing\n",fname);
		return 0;
	}

	fprintf(stderr,"Writing results to '%s'\n",fname);

	if((rng_ctx=gsl_rng_alloc(gsl_rng_mt19937))==NULL)
	{
		printf("Couldn't initialize the random number generator\n");
		exit(0);
	}

	if(config.seedrng)
		seed_rng(rng_ctx);

	dgr=init_diagram(config.endtau,config.j,config.m,config.chempot);
	ht=init_histogram(config.bins,config.width);

#define DIAGRAM_UPDATE_LENGTH	(0)
#define DIAGRAM_ADD_LINE	(1)
#define DIAGRAM_REMOVE_LINE	(2)

#warning FIXME Should be 3
#define DIAGRAM_NR_UPDATES	(1)

	if(config.progressbar)
		progress=progressbar_new("Progress",config.iterations/16384);
	else
		progress=NULL;

	for(c=0;c<config.iterations;c++)
	{
		int update_type;
		double oldweight=diagram_weight(dgr);
		bool accepted;
		
		update_type=gsl_rng_uniform_int(rng_ctx,DIAGRAM_NR_UPDATES);
			
		switch(update_type)
		{
			case DIAGRAM_UPDATE_LENGTH:
			{
				int nr_midpoints=get_nr_midpoints(dgr);
				int nr_free_propagators=get_nr_free_propagators(dgr);

				struct g0_t *lastg0;				
				double lasttau,oldendtau,newendtau,rate;

				lasttau=get_midpoint(dgr,nr_midpoints-1);
				lastg0=get_free_propagator(dgr,nr_free_propagators-1);
				
				/*
					We calculate newendtau between lasttau and maxtau
					using an exponential distribution with rate rate.
				*/

				oldendtau=dgr->endtau;
				rate=lastg0->j*(lastg0->j+1.0f)-dgr->chempot;
				newendtau=doubly_truncated_exp_dist(rng_ctx,rate,lasttau,config.maxtau);

				/*
					This update is always accepted
				*/
	
				diagram_update_length(dgr,newendtau);
				accepted=true;

				/*
					This codepath is taken only when debugging, as the
					update is always accepted
				*/

				if(accepted==false)
				{
					diagram_update_length(dgr,oldendtau);
					assert(diagram_weight(dgr)==oldweight);
				}
			}
			break;

			case DIAGRAM_ADD_LINE:
			{
				//double k=
				//int lambda=
				//int mu=
				//double tau1=
				//double tau2=
				//diagram_add_phonon_line(dgr,tau1,tau2,k,lambda,mu);
				//recouple

				// Calculate the accepted variable
				// accepted = ...

				diagram_add_phonon_line(dgr,0.5,0.6,1,1,0);

				accepted=false;

				if(accepted==false)
				{
					int lastline=get_nr_phonons(dgr)-1;
					diagram_remove_phonon_line(dgr,lastline);
					// Fix free propagators, they should have be saved!
					assert(diagram_weight(dgr)==oldweight);
				}

			}
			break;
			
			case DIAGRAM_REMOVE_LINE:
			{
				//int target=gsl_rng_uniform_int(rng_ctx,get_nr_phonon_lines(nr));
				//diagram_remove_line(dgr,target);
				//recouple
			}
			break;
			
			//default:
			//assert(false);
		}

		// Update stats here!
		// Update type
		// accepted?

		histogram_add_sample(ht,diagram_weight(dgr),dgr->endtau);

		if((config.progressbar)&&((c%16384)==0))
			progressbar_inc(progress);

		// FIXME Optimization: this weight can be used later in oldweight, no need to recalculate
	}

	if(config.progressbar)
		progressbar_finish(progress);

	fprintf(out,"# Diagrammatic Monte Carlo for the angulon\n");
	fprintf(out,"#\n");
	fprintf(out,"# Configuration loaded from '%s'\n",configfile);
	fprintf(out,"# Output file is '%s'\n",fname);
	fprintf(out,"#\n");
	fprintf(out,"# Initial and final state: (j=%d, m=%d)\n",config.j,config.m);
	fprintf(out,"# Chemical potential: %f\n",config.chempot);
	fprintf(out,"# Initial diagram length: %f\n",config.endtau);
	fprintf(out,"# Max length: %f\n",config.maxtau);
	fprintf(out,"#\n");
	fprintf(out,"# Sampled quantity: Green's function (G)\n");
	fprintf(out,"# Iterations: %d\n",config.iterations);
	fprintf(out,"# Nr of bins: %d\n",config.bins);
	fprintf(out,"# Bin width: %f\n",config.width);
	fprintf(out,"# (last bin is overflow)\n");
	fprintf(out,"#\n");

	//fprintf(out,"# Update statistics:\n");
	//fprintf(out,"# Update #1: proposed %d (%f%%), accepted %d (%f%%), rejected %d (%f%%).\n");
	//fprintf(out,"# Update #2: proposed %d (%f%%), accepted %d (%f%%), rejected %d (%f%%).\n");
	//fprintf(out,"# Update #3: proposed %d (%f%%), accepted %d (%f%%), rejected %d (%f%%).\n");
	//fprintf(out,"# Total: proposed %d (%f%%), accepted %d (%f%%), rejected %d (%f%%).\n");

	for(c=0;c<=config.bins;c++)
		fprintf(out,"%f %f %f\n",config.width*c+config.width/2.0f,histogram_get_bin_average(ht,c),histogram_get_bin_variance(ht,c));

	fini_histogram(ht);
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
