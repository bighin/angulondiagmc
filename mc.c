#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

#include "diagrams.h"
#include "updates.h"
#include "stat.h"
#include "aux.h"
#include "debug.h"

#include "inih/ini.h"
#include "libprogressbar/progressbar.h"

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

#define UPDATE_UNPHYSICAL	(0)
#define UPDATE_REJECTED		(1)
#define UPDATE_ACCEPTED		(2)
#define UPDATE_ERROR		(3)

int update_length(struct diagram_t *dgr,struct configuration_t *cfg)
{
	int nr_midpoints,nr_free_propagators;
	double lasttau,oldendtau,newendtau,rate,oldweight;
	struct g0_t *lastg0;
	bool is_accepted;

	oldweight=diagram_weight(dgr);

	nr_midpoints=get_nr_midpoints(dgr);
	nr_free_propagators=get_nr_free_propagators(dgr);

	lasttau=get_midpoint(dgr,nr_midpoints-1);
	lastg0=get_free_propagator(dgr,nr_free_propagators-1);
	
	/*
		We calculate newendtau between lasttau and maxtau
		using an exponential distribution with rate rate.
	*/

	oldendtau=dgr->endtau;
	rate=lastg0->j*(lastg0->j+1.0f)-dgr->chempot;
	newendtau=doubly_truncated_exp_dist(dgr->rng_ctx,rate,lasttau,cfg->maxtau);

	/*
		This update is always accepted
	*/

	diagram_update_length(dgr,newendtau);
	is_accepted=true;

	/*
		This codepath is taken only when debugging, as the
		update is always accepted.
	*/

	if(is_accepted==false)
	{
		diagram_update_length(dgr,oldendtau);
		assert(diagram_weight(dgr)==oldweight);
	
		return UPDATE_REJECTED;
	}
	
	return UPDATE_ACCEPTED;
}

int update_add_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg)
{
	int lambda,mu;
	double k,tau1,tau2,acceptance_ratio,oldweight;
	bool is_accepted;

	struct arc_t *thisline;
	struct vertex_info_t *v1,*v2;

	double omegas[3]={cfg->omega0,cfg->omega1,cfg->omega2};

	/*
		Do we want to limit the simulation only for first order diagrams?
	*/

#if 0
	if(get_nr_phonons(dgr)>=1)
		return UPDATE_UNPHYSICAL;
#endif

	oldweight=diagram_weight(dgr);

	/*
		The new lambda and mu are chosen using a uniform distribution
	*/
	
	lambda=gsl_rng_uniform_int(dgr->rng_ctx,3);
	mu=gsl_rng_uniform_int(dgr->rng_ctx,2*lambda+1)-lambda;

	/*
		The start time tau1 is sampled uniformly, while tau2 is sampled from
		a truncated exponential distribution.
	*/
	
	k=0.0f;
	tau1=gsl_rng_uniform(dgr->rng_ctx)*dgr->endtau;
	tau2=doubly_truncated_exp_dist(dgr->rng_ctx,omegas[lambda],tau1,dgr->endtau);

	/*
		The line is added...
	*/

	diagram_add_phonon_line(dgr,tau1,tau2,k,lambda,mu);
	
	/*
		...and then we check if the new line makes sense physically!
	*/

	thisline=get_phonon_line(dgr,get_nr_phonons(dgr)-1);
	v1=get_vertex(dgr,thisline->startmidpoint);
	v2=get_vertex(dgr,thisline->endmidpoint);

	if((recouple_ms(dgr)==false)||(check_triangle_condition(dgr,v1)==false)||(check_triangle_condition(dgr,v2)==false))
	{
		int lastline;
		bool result;
			
		lastline=get_nr_phonons(dgr)-1;
		diagram_remove_phonon_line(dgr,lastline);

		result=recouple_ms(dgr);
		assert(result==true);

		return UPDATE_UNPHYSICAL;
	}

	acceptance_ratio=diagram_weight(dgr)/oldweight;
	acceptance_ratio*=3.0f*(2*lambda+1)*dgr->endtau;
	acceptance_ratio/=get_nr_phonons(dgr)+1;
	acceptance_ratio/=doubly_truncated_exp_pdf(dgr->rng_ctx,omegas[lambda],tau1,dgr->endtau,tau2);

	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
		int lastline;
		bool result;

		lastline=get_nr_phonons(dgr)-1;
		diagram_remove_phonon_line(dgr,lastline);
		
		result=recouple_ms(dgr);
		assert(result==true);

		assert(fabs(diagram_weight(dgr)-oldweight)<1e-7*oldweight);

		return UPDATE_REJECTED;
	}

	assert(check_couplings_ms(dgr)==true);

	return UPDATE_ACCEPTED;
}

int update_remove_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg)
{
	struct arc_t *arc;

	int target,lambda,mu,nr_phonons;
	double acceptance_ratio,k,tau1,tau2,oldweight;
	bool is_accepted;

	double omegas[3]={cfg->omega0,cfg->omega1,cfg->omega2};

	oldweight=diagram_weight(dgr);
	nr_phonons=get_nr_phonons(dgr);
	
	if(nr_phonons<=0)
		return UPDATE_UNPHYSICAL;
	
	target=gsl_rng_uniform_int(dgr->rng_ctx,nr_phonons);
	arc=get_phonon_line(dgr,target);

	tau1=arc->starttau;
	tau2=arc->endtau;
	k=arc->k;
	lambda=arc->lambda;
	mu=arc->mu;

	if((get_vertex(dgr,arc->startmidpoint)->refs!=0)||(get_vertex(dgr,arc->endmidpoint)->refs!=0))
		return false;

	diagram_remove_phonon_line(dgr,target);

	if(recouple_ms(dgr)==false)
	{
		bool result;

		diagram_add_phonon_line(dgr,tau1,tau2,k,lambda,mu);

		result=recouple_ms(dgr);		

		assert(result==true);
		assert(fabs(diagram_weight(dgr)-oldweight)<1e-7*oldweight);
	
		return UPDATE_REJECTED;
	}

	acceptance_ratio=diagram_weight(dgr)/oldweight;
	acceptance_ratio/=3.0f*(2*lambda+1)*dgr->endtau;
	acceptance_ratio*=nr_phonons;
	acceptance_ratio*=doubly_truncated_exp_pdf(dgr->rng_ctx,omegas[lambda],tau1,dgr->endtau,tau2);

	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
		bool result;

		diagram_add_phonon_line(dgr,tau1,tau2,k,lambda,mu);

		result=recouple_ms(dgr);
	
		assert(result==true);
		assert(fabs(diagram_weight(dgr)-oldweight)<1e-7*oldweight);
	
		return UPDATE_REJECTED;
	}

	return UPDATE_ACCEPTED;
}

int update_add_worm(struct diagram_t *dgr,struct configuration_t *cfg)
{
	int nr_midpoints,target1,target2,lambda1,lambda2,minlambda,deltalambda;
	bool is_accepted;
	double oldweight,acceptance_ratio;

	/*
		We select two midpoints, with uniform probability.
	*/

	nr_midpoints=get_nr_midpoints(dgr);

	if(nr_midpoints<2)
		return UPDATE_UNPHYSICAL;

	target1=gsl_rng_uniform_int(dgr->rng_ctx,nr_midpoints);
	target2=gsl_rng_uniform_int(dgr->rng_ctx,nr_midpoints-1);

	if(target2>=target1)
		target2++;

	lambda1=vertex_get_lambda(get_vertex(dgr,target1));
	lambda2=vertex_get_lambda(get_vertex(dgr,target2));

	oldweight=diagram_weight(dgr);

	/*
		deltalambda is chosen with uniform probability
		between -min(lambda1,lambda2) and min(lambda1,lambda2)
	*/

	minlambda=MIN(lambda1,lambda2);
	deltalambda=gsl_rng_uniform_int(dgr->rng_ctx,2*minlambda+1)-minlambda;
			
	if(diagram_add_worm(dgr,MIN(target1,target2),MAX(target1,target2),deltalambda)==false)
	{
		assert(check_couplings_ms(dgr)==true);
		return UPDATE_UNPHYSICAL;
	}

	assert(check_couplings_ms(dgr)==true);

	acceptance_ratio=diagram_weight(dgr)/oldweight;
	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
		bool result;
		
		result=diagram_remove_worm(dgr,get_nr_worms(dgr)-1);

		if(result==false)
		{
			print_diagram(dgr,PRINT_TOPOLOGY|PRINT_PROPAGATORS);

			assert(false);
		}

		assert(check_couplings_ms(dgr)==true);

		assert(fabs(diagram_weight(dgr)-oldweight)<1e-7*oldweight);
	
		return UPDATE_REJECTED;
	}

	return UPDATE_ACCEPTED;
}

int update_remove_worm(struct diagram_t *dgr,struct configuration_t *cfg)
{
	double oldweight,acceptance_ratio;
	int nr_worms,target,v1,v2,deltalambda;
	struct worm_t *worm;
	bool is_accepted;

	nr_worms=get_nr_worms(dgr);

	if(nr_worms<=0)
		return UPDATE_UNPHYSICAL;

	oldweight=diagram_weight(dgr);
	target=gsl_rng_uniform_int(dgr->rng_ctx,nr_worms);

	worm=get_worm(dgr,target);
	v1=worm->startmidpoint;
	v2=worm->endmidpoint;
	deltalambda=worm->deltalambda;

	if(diagram_remove_worm(dgr,target)==false)
	{
		assert(check_couplings_ms(dgr)==true);
		return UPDATE_UNPHYSICAL;
	}

	assert(check_couplings_ms(dgr)==true);

	acceptance_ratio=diagram_weight(dgr)/oldweight;
	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
		bool result;
		
		result=diagram_add_worm(dgr,v1,v2,deltalambda);

		assert(result==true);
		assert(check_couplings_ms(dgr)==true);
		assert(fabs(diagram_weight(dgr)-oldweight)<1e-7*oldweight);

		return UPDATE_REJECTED;
	}

	return UPDATE_ACCEPTED;
}

int update_change_lambda(struct diagram_t *dgr,struct configuration_t *cfg)
{
	return UPDATE_UNPHYSICAL;
}

void load_config_defaults(struct configuration_t *config)
{
	config->prefix=strdup("default");
	config->seedrng=false;
	config->progressbar=true;

	config->j=1;
	config->m=0;
	config->endtau=1.0f;
	config->chempot=0.95f;
	config->maxtau=100.0f;
	
	config->c0=1.5;
	config->c1=0.75;
	config->c2=0.35;
	config->omega0=0.8;
	config->omega1=0.8;
	config->omega2=0.8;

	config->iterations=10000000;
	config->bins=50;
	config->width=0.25;
}

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
	else if(MATCH("potential","c0"))
	{
		pconfig->c0=atof(value);
	}
	else if(MATCH("potential","c1"))
	{
		pconfig->c1=atof(value);
	}
	else if(MATCH("potential","c2"))
	{
		pconfig->c2=atof(value);
	}
	else if(MATCH("potential","omega0"))
	{
		pconfig->omega0=atof(value);
	}
	else if(MATCH("potential","omega1"))
	{
		pconfig->omega1=atof(value);
	}
	else if(MATCH("potential","omega2"))
	{
		pconfig->omega2=atof(value);
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

void show_update_statistics(FILE *out,int proposed,int accepted,int rejected)
{
	double accepted_pct,rejected_pct;

	if(proposed>0)
	{
		accepted_pct=100.0f*((double)(accepted))/((double)(proposed));
		rejected_pct=100.0f*((double)(rejected))/((double)(proposed));
	}
	else
	{
		accepted_pct=rejected_pct=0.0f;
	}

	fprintf(out,"proposed %d, accepted %d (%f%%), rejected %d (%f%%).\n",proposed,accepted,accepted_pct,rejected,rejected_pct);
}

int do_diagmc(char *configfile)
{
	struct diagram_t *dgr;
	struct diagram_cfg_t dcfg;
	struct configuration_t config;
	struct histogram_t *ht;

	int c;	
	FILE *out;
	char fname[1024];
	progressbar *progress;

#define DIAGRAM_NR_UPDATES	(5)

	int (*updates[DIAGRAM_NR_UPDATES])(struct diagram_t *,struct configuration_t *);

	int proposed[DIAGRAM_NR_UPDATES],accepted[DIAGRAM_NR_UPDATES],rejected[DIAGRAM_NR_UPDATES];
	int total_proposed,total_accepted,total_rejected;
	int avgorder[2];

	updates[0]=update_length;
	updates[1]=update_add_phonon_line;
	updates[2]=update_remove_phonon_line;
	updates[3]=update_add_worm;
	updates[4]=update_remove_worm;
	//updates[5]=update_change_lambda;

	for(c=0;c<DIAGRAM_NR_UPDATES;c++)
		proposed[c]=accepted[c]=rejected[c]=0;
	
	avgorder[0]=avgorder[1]=0;
	
	/*
		We load some sensible defaults in case they are not in the .ini file
	*/

	load_config_defaults(&config);

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

	dcfg.j=config.j;
	dcfg.m=config.m;
	dcfg.endtau=config.endtau;
	dcfg.chempot=config.chempot;

	dcfg.c0=config.c0;
	dcfg.c1=config.c1;
	dcfg.c2=config.c2;
	dcfg.omega0=config.omega0;

	dgr=init_diagram(&dcfg);
	ht=init_histogram(config.bins,config.width);

	if(config.seedrng)
		seed_rng(dgr->rng_ctx);

	if(config.progressbar)
		progress=progressbar_new("Progress",config.iterations/65536);
	else
		progress=NULL;

	for(c=0;c<config.iterations;c++)
	{
		int update_type,status;

		if(check_couplings_ms(dgr)!=true)
		{
			printf("Coupling error!\n");
			printf("<<%s>>\n",recouple_ms(dgr)?"true":"false");
			printf("<<%s>>\n",check_couplings_ms(dgr)?"true":"false");
		
			assert(false);
		}

		update_type=gsl_rng_uniform_int(dgr->rng_ctx,DIAGRAM_NR_UPDATES);
		status=updates[update_type](dgr,&config);

		//print_diagram(dgr,PRINT_TOPOLOGY|PRINT_PROPAGATORS);

		if(check_couplings_ms(dgr)!=true)
		{
			printf("Coupling error!\n");
			printf("<<%s>>\n",recouple_ms(dgr)?"true":"false");
			printf("<<%s>>\n",check_couplings_ms(dgr)?"true":"false");
			printf("After update type: %d, with status: %d\n",update_type,status);

			assert(false);
		}

		switch(status)
		{
			case UPDATE_ACCEPTED:
			proposed[update_type]++;
			accepted[update_type]++;
			break;

			case UPDATE_REJECTED:
			proposed[update_type]++;
			rejected[update_type]++;
			break;

#warning CHECKME!

			case UPDATE_UNPHYSICAL:
			break;

			case UPDATE_ERROR:
			assert(false);
		}

		histogram_add_sample(ht,diagram_weight(dgr),dgr->endtau);

		if((config.progressbar)&&((c%65536)==0))
			progressbar_inc(progress);

		avgorder[0]+=get_nr_phonons(dgr);
		avgorder[1]++;

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
	fprintf(out,"# Max diagram length: %f\n",config.maxtau);
	fprintf(out,"# Potential parameters: (c0, c1, c2, omega0) = (%f, %f, %f, %f)\n",config.c0,config.c1,config.c2,config.omega0);
	fprintf(out,"#\n");
	fprintf(out,"# Sampled quantity: Green's function (G)\n");
	fprintf(out,"# Iterations: %d\n",config.iterations);
	fprintf(out,"# Nr of bins: %d\n",config.bins);
	fprintf(out,"# Bin width: %f\n",config.width);
	fprintf(out,"# (last bin is overflow)\n");
	fprintf(out,"#\n");

	total_proposed=total_accepted=total_rejected=0;

	fprintf(out,"# Update statistics:\n");

	for(c=0;c<DIAGRAM_NR_UPDATES;c++)
	{
		fprintf(out,"# Update #%d: ",c);
		show_update_statistics(out,proposed[c],accepted[c],rejected[c]);
	
		total_proposed+=proposed[c];
		total_accepted+=accepted[c];
		total_rejected+=rejected[c];
	}

	fprintf(out,"# Total: ");
	show_update_statistics(out,total_proposed,total_accepted,total_rejected);
	fprintf(out,"# Average order: %f\n",((double)(avgorder[0]))/((double)(avgorder[1])));
	fprintf(out,"#\n");

	fprintf(out,"# <Bin center> <Average> <Stddev>\n");

	for(c=0;c<=config.bins;c++)
	{
		fprintf(out,"%f %f ",config.width*c+config.width/2.0f,histogram_get_bin_average(ht,c));
		
		if(ht->sctx->smpls[c]->next>1)
			fprintf(out,"%f\n",sqrtf(histogram_get_bin_variance(ht,c)/(ht->sctx->smpls[c]->next)));
		else
			fprintf(out,"%f\n",0.0f);
	}

	fini_histogram(ht);
	fini_diagram(dgr);

	return 0;
}