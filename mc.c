#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <unistd.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "diagrams.h"
#include "updates.h"
#include "stat.h"
#include "aux.h"
#include "graphs.h"
#include "debug.h"
#include "mc.h"

#include "inih/ini.h"
#include "libprogressbar/progressbar.h"

/*
	The following headers need to be included after graphs.h, becuase it defines a 'lines' macro
	conflicting with a field in struct graph_t
*/

#include <curses.h>
#include <term.h>

static char term_buffer[2048];

void init_terminal_data(void)
{
	char *termtype=getenv("TERM");
	int success;

	if(termtype==NULL)
		fprintf(stderr,"Please specify a terminal type with 'setenv TERM <yourtype>'.\n");

	success=tgetent(term_buffer,termtype);
	
	if(success<0)
		fprintf(stderr, "Could not access the termcap database.\n");

	if(success==0)
		fprintf(stderr, "Terminal type `%s' is not defined.\n",termtype);
}

#define UPDATE_UNPHYSICAL	(0)
#define UPDATE_REJECTED		(1)
#define UPDATE_ACCEPTED		(2)
#define UPDATE_ERROR		(3)

int update_length(struct diagram_t *dgr,struct configuration_t *cfg)
{
	int nr_midpoints,nr_free_propagators;
	double lasttau,oldendtau,newendtau,rate;
	struct g0_t *lastg0;
	bool is_accepted;

#ifndef NDEBUG
	double oldweight=diagram_weight(dgr);
#endif

	nr_midpoints=get_nr_midpoints(dgr);
	nr_free_propagators=get_nr_free_propagators(dgr);

	lasttau=get_midpoint(dgr,nr_midpoints-1);
	lastg0=get_free_propagator(dgr,nr_free_propagators-1);

	/*
		We calculate newendtau between lasttau and maxtau
		using an exponential distribution with rate 'rate'.
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
	double k,tau1,tau2,acceptance_ratio,omegaeff,g0eff;
	bool is_accepted;

	struct arc_t *thisline;
	struct vertex_t *v1,*v2;
	struct diagram_t *old;

	/*
		Do we want to limit the simulation only to first order diagrams (or at a certain maximum order)?
	*/

	if(get_nr_phonons(dgr)>=cfg->maxorder)
		return UPDATE_UNPHYSICAL;

	/*
		The new lambda is chosen using a uniform distribution.
		Note that mu does not need to be initialized, we will sum over all mu's
		when evaluating the weight of a diagram.
	*/

	lambda=gsl_rng_uniform_int(dgr->rng_ctx,2);
	mu=0;

	omegaeff=(lambda==0)?(fabs(dgr->v0slope)):(fabs(dgr->v1slope));
	g0eff=(lambda==0)?(fabs(dgr->v0intercept)):(fabs(dgr->v1intercept));

	/*
		The start time tau1 is sampled uniformly, whereas tau2 is sampled from
		a truncated exponential distribution.
	*/
	
	k=0.0f;
	tau1=gsl_rng_uniform(dgr->rng_ctx)*dgr->endtau;
	tau2=doubly_truncated_exp_dist(dgr->rng_ctx,omegaeff,tau1,dgr->endtau);

	/*
		We save the current diagram and the line is added...
	*/

	old=diagram_clone(dgr);
	diagram_add_phonon_line(dgr,tau1,tau2,k,lambda,mu);

	/*
		...and then we check if the new line makes sense physically!
	*/

	thisline=get_phonon_line(dgr,get_nr_phonons(dgr)-1);
	v1=get_vertex(dgr,thisline->startmidpoint);
	v2=get_vertex(dgr,thisline->endmidpoint);

	if((check_triangle_condition_and_parity(dgr,v1)==false)||(check_triangle_condition_and_parity(dgr,v2)==false))
	{
		diagram_copy(old,dgr);
		fini_diagram(old);

		return UPDATE_UNPHYSICAL;
	}

	recouple(dgr,thisline->startmidpoint,thisline->endmidpoint);

	acceptance_ratio=(diagram_weight(dgr)*diagram_m_weight(dgr,cfg->use_hashtable))/(diagram_weight(old)*diagram_m_weight(old,cfg->use_hashtable));
	acceptance_ratio*=3.0f*dgr->endtau;
	acceptance_ratio/=get_nr_phonons(dgr)+1;
	acceptance_ratio/=doubly_truncated_exp_pdf(dgr->rng_ctx,omegaeff,tau1,dgr->endtau,tau2);
	acceptance_ratio/=g0eff;

	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
		diagram_copy(old,dgr);
		fini_diagram(old);

		return UPDATE_REJECTED;
	}

	fini_diagram(old);

	return UPDATE_ACCEPTED;
}

int update_remove_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg)
{
	struct arc_t *arc;
	struct diagram_t *old;
	
	int target,lambda,nr_phonons,startmidpoint,endmidpoint;
	double acceptance_ratio,tau1,tau2,omegaeff,g0eff;
	bool is_accepted;

	nr_phonons=get_nr_phonons(dgr);

	if(nr_phonons<=0)
		return UPDATE_UNPHYSICAL;
	
	target=gsl_rng_uniform_int(dgr->rng_ctx,nr_phonons);
	arc=get_phonon_line(dgr,target);

	tau1=arc->starttau;
	tau2=arc->endtau;
	lambda=arc->lambda;
	startmidpoint=arc->startmidpoint;
	endmidpoint=arc->endmidpoint;

	old=diagram_clone(dgr);

	/*
		A line removal can fail, since it may leave the remaining
		lines in a unphysical state.
	*/
	
	if(diagram_remove_phonon_line(dgr,target)==false)
	{	
		diagram_copy(old,dgr);
		fini_diagram(old);
	
		return UPDATE_UNPHYSICAL;
	}

	/*
		This is a bit tricky: if we just removed a bubble,
		then there are no free propagators to recouple.
	*/

	if((startmidpoint+1)!=endmidpoint)
	{
		recouple(dgr,startmidpoint,endmidpoint-2);
	}
	else
	{
		/*
			...however the bubble may be 'asymmetric': that would be an error!
		*/
	
		if(get_left_neighbour(old,startmidpoint)->j!=get_right_neighbour(old,endmidpoint)->j)
		{
			print_diagram(dgr,PRINT_TOPOLOGY|PRINT_PROPAGATORS);
		
			printf("%f\n",diagram_weight(old)*diagram_m_weight(old,cfg->use_hashtable));

			exit(0);
		}
	}

	omegaeff=(lambda==0)?(fabs(dgr->v0slope)):(fabs(dgr->v1slope));
	g0eff=(lambda==0)?(fabs(dgr->v0intercept)):(fabs(dgr->v1intercept));

	acceptance_ratio=(diagram_weight(dgr)*diagram_m_weight(dgr,cfg->use_hashtable))/(diagram_weight(old)*diagram_m_weight(old,cfg->use_hashtable));
	acceptance_ratio/=3.0f*dgr->endtau;
	acceptance_ratio*=nr_phonons;
	acceptance_ratio*=doubly_truncated_exp_pdf(dgr->rng_ctx,omegaeff,tau1,dgr->endtau,tau2);
	acceptance_ratio*=g0eff;

	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
		diagram_copy(old,dgr);
		fini_diagram(old);

		return UPDATE_REJECTED;
	}

	fini_diagram(old);

	return UPDATE_ACCEPTED;
}

void load_config_defaults(struct configuration_t *config)
{
	config->prefix=strdup("default");
	config->seedrng=false;
	config->progressbar=true;
	config->animate=false;

	config->j=1;
	config->endtau=1.0f;
	config->chempot=0.95f;
	config->maxtau=100.0f;
	config->maxorder=1024;
	
	config->n=exp(-10.0f);

	config->iterations=10000000;
	config->bins=50;
	config->width=0.25;

	config->parallel=false;
	config->nthreads=1;
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
	else if(MATCH("general","animate"))
	{
		if(!strcmp(value,"true"))
			pconfig->animate=true;
		else if(!strcmp(value,"false"))
			pconfig->animate=false;
		else
			return 0;
	}
	else if(MATCH("general","hashtable"))
	{
		if(!strcmp(value,"true"))
			pconfig->use_hashtable=true;
		else if(!strcmp(value,"false"))
			pconfig->use_hashtable=false;
		else
			return 0;
	}
	else if(MATCH("parameters","j"))
	{
		pconfig->j=atoi(value);
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
	else if(MATCH("parameters","maxorder"))
	{
		pconfig->maxorder=atoi(value);
	}
	else if(MATCH("potential","logn"))
	{
		pconfig->n=exp(atof(value));
	}
	else if(MATCH("sampling","bins"))
	{
		pconfig->bins=atoi(value);
	}
	else if(MATCH("sampling","iterations"))
	{
		pconfig->iterations=atoi(value);
	}
	else if(MATCH("sampling","width"))
	{	
		pconfig->width=atof(value);
	}
	else if(MATCH("parallel","parallel"))
	{
		if(!strcmp(value,"true"))
			pconfig->parallel=true;
		else if(!strcmp(value,"false"))
			pconfig->parallel=false;
		else
			return 0;
	}
	else if(MATCH("parallel","nthreads"))
	{		
		if(!strcmp(value,"auto"))
		{
#if defined(_OPENMP)
			pconfig->nthreads=omp_get_max_threads();
#else
			pconfig->nthreads=1;
#endif
		}
		else
			pconfig->nthreads=atoi(value);
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
	struct diagram_cfg_t dcfg;
	struct configuration_t config;
	struct histogram_t *ht;

	int d;
	FILE *out;
	char fname[1024];
	progressbar *progress;

#define DIAGRAM_NR_UPDATES	(3)

	int (*updates[DIAGRAM_NR_UPDATES])(struct diagram_t *,struct configuration_t *);
	char *update_names[DIAGRAM_NR_UPDATES];

	int proposed[DIAGRAM_NR_UPDATES],accepted[DIAGRAM_NR_UPDATES],rejected[DIAGRAM_NR_UPDATES];
	int total_proposed,total_accepted,total_rejected;
	int avgorder[2];

	updates[0]=update_length;
	updates[1]=update_add_phonon_line;
	updates[2]=update_remove_phonon_line;

	update_names[0]="UpdateLength";
	update_names[1]="AddPhononLine";
	update_names[2]="RemovePhononLine";

	for(d=0;d<DIAGRAM_NR_UPDATES;d++)
		proposed[d]=accepted[d]=rejected[d]=0;
	
	avgorder[0]=avgorder[1]=0;

	/*
		We load some sensible defaults in case they are not in the .ini file, the we try
		loading a .ini file.
	*/

	load_config_defaults(&config);

	if(ini_parse(configfile,configuration_handler,&config)<0)
	{
		fprintf(stderr,"Couldn't read or parse '%s'\n",configfile);
		return 1;
	}

	if((config.animate==true)&&(config.progressbar==true))
	{
		fprintf(stderr,"The options 'animate' and 'progressbar' cannot be set at the same time!\n");
		exit(0);
	}

	if((config.animate==true)&&(config.parallel==true))
	{
		fprintf(stderr,"The options 'animate' and 'parallel' cannot be set at the same time!\n");
		exit(0);
	}

#if !defined(_OPENMP)
	if(config.nthreads>1)
	{
		fprintf(stderr,"Warning: parallel run requested, but the current binary has been compiled without OpenMP support.\n");
		config.nthreads=1;
	}
#endif

	if(config.j*(config.j+1)<=config.chempot)
	{
		fprintf(stderr,"Warning: the chemical potential should be set below j(j+1).\n");
	}

	if((config.parallel==true)&&(config.use_hashtable==true))
	{
		fprintf(stderr,"Warning: the hashtable is incompatible with the 'parallel' option. The hashtable will be switched off.\n");
		config.use_hashtable=false;
	}

	if(config.parallel==false)
		config.nthreads=1;

	fprintf(stderr,"Loaded '%s'\n",configfile);
	fprintf(stderr,"Performing %d iterations, using %d thread(s)\n",config.iterations,config.nthreads);

	snprintf(fname,1024,"%s.dat",config.prefix);
	
	if(!(out=fopen(fname,"w+")))
	{
		fprintf(stderr,"Error: couldn't open %s for writing\n",fname);
		return 0;
	}

	fprintf(stderr,"Writing results to '%s'\n",fname);

	dcfg.j=config.j;
	dcfg.endtau=config.endtau;
	dcfg.maxtau=config.maxtau;
	dcfg.chempot=config.chempot;

	dcfg.n=config.n;

        init_terminal_data();

	if(config.progressbar)
		progress=progressbar_new("Progress",config.iterations/16384);
	else
		progress=NULL;

	ht=init_histogram(config.bins,config.width);

	/*
		This is the main DiagMC loop
	*/

#if defined(_OPENMP)
#pragma omp parallel for
#endif

	for(d=0;d<config.nthreads;d++)
	{
		struct diagram_t *dgr;
		int c;

#define LOCALSTORAGE	(1024*16)

		double localtotalweight[LOCALSTORAGE][2];
		int localavgorder[2];
		int progressbarticks;

		localavgorder[0]=localavgorder[1]=0;
		progressbarticks=0;

		/*
			We initialize some basic data structures and various other stuff, as well as the
			random number generator, in case a NON deterministic seed is requested.
		*/

		dgr=init_diagram(&dcfg);

		dgr->rng_ctx=gsl_rng_alloc(gsl_rng_mt19937);
		assert(dgr->rng_ctx!=NULL);

		if(config.seedrng)
			seed_rng(dgr->rng_ctx);

		/*
			If we are running many threads in parallel, we must initialized
			the RNG in each thread in a different way, otherwise we are sampling
			the same diagrams.
		
			In case a random seed has NOT been required (config.seedrng==false),
			we use a deterministic seed (just the thread number).
		*/

		if((config.parallel)&&(!config.seedrng))
			gsl_rng_set(dgr->rng_ctx,d);

		for(c=0;c<config.iterations/config.nthreads;c++)
		{
			int update_type,status;
			double totalweight;

			update_type=gsl_rng_uniform_int(dgr->rng_ctx,DIAGRAM_NR_UPDATES);
			status=updates[update_type](dgr,&config);

			if((config.nthreads==1)&&(config.animate)&&(status==UPDATE_ACCEPTED)&&((update_type==1)||(update_type==2)))
			{
				if((c%16)==0)
				{
					int d,plines;

					plines=print_diagram(dgr,PRINT_TOPOLOGY|PRINT_INFO0|PRINT_DRYRUN);
		
					for(d=plines;d<tgetnum("li");d++)
							printf("\n");

					print_diagram(dgr,PRINT_TOPOLOGY|PRINT_INFO0);
					usleep(5000);
				}
			}

#if defined(_OPENMP)				
#pragma omp atomic
#endif
			proposed[update_type]++;

			switch(status)
			{
				case UPDATE_ACCEPTED:

#if defined(_OPENMP)				
#pragma omp atomic
#endif
				accepted[update_type]++;

				break;

				/*
					FIXME

					Here unphysical and rejected updates are treated in exactly the same
					way. Is this OK?
				*/

				case UPDATE_UNPHYSICAL:
				case UPDATE_REJECTED:

#if defined(_OPENMP)
#pragma omp atomic
#endif
				rejected[update_type]++;

				break;

				case UPDATE_ERROR:
				assert(false);
			}

			assert(fabs(diagram_weight(dgr)-diagram_weight_non_incremental(dgr))<10e-7*diagram_weight(dgr));
			diagram_check_consistency(dgr);

			/*
				Note that by using diagram_m_weight() we are effectively summing over
				all values of m, so that we have to divide the weight by (2j + 1).
			*/

			totalweight=diagram_weight(dgr)*diagram_m_weight(dgr,config.use_hashtable)/(2.0f*dcfg.j+1.0f);
			assert(totalweight>=0.0f);

			localtotalweight[c%LOCALSTORAGE][0]=totalweight;
			localtotalweight[c%LOCALSTORAGE][1]=dgr->endtau;

			localavgorder[0]+=get_nr_phonons(dgr);
			localavgorder[1]++;

			if((config.progressbar)&&((c%16384)==0))
				progressbarticks++;

			if(((c+1)%LOCALSTORAGE)==0)
			{
#if defined(_OPENMP)
#pragma omp atomic
#endif
				avgorder[0]+=localavgorder[0];

#if defined(_OPENMP)
#pragma omp atomic
#endif

				avgorder[1]+=localavgorder[1];
				localavgorder[0]=localavgorder[1]=0;
		
#if defined(_OPENMP)
#pragma omp critical
#endif
				{
					int i;
					
					for(i=0;i<LOCALSTORAGE;i++)
						histogram_add_sample(ht,localtotalweight[i][0],localtotalweight[i][1]);
				}

				if(config.progressbar)
				{
#if defined(_OPENMP)
#pragma omp atomic
#endif
					progress->value+=progressbarticks;
				}

				progressbarticks=0;

#if defined(_OPENMP)
				if((config.progressbar)&&(omp_get_thread_num()==0))
				{
					progressbar_update(progress,progress->value);
				}
#else
				if(config.progressbar)
				{
					progressbar_update(progress,progress->value);
				}
#endif		
			}
		}

		if(dgr->rng_ctx)
			gsl_rng_free(dgr->rng_ctx);

		fini_diagram(dgr);
	}

	if(config.progressbar)
		progressbar_finish(progress);

	/*
		Now we print the statistics we collected to the output file in a nice way.
	*/

	fprintf(out,"# Diagrammatic Monte Carlo for the angulon\n");
	fprintf(out,"#\n");
	fprintf(out,"# Configuration loaded from '%s'\n",configfile);
	fprintf(out,"# Output file is '%s'\n",fname);
	fprintf(out,"#\n");
	fprintf(out,"# Initial and final state: (j=%d, sampling over all possible m values)\n",config.j);
	fprintf(out,"# Chemical potential: %f\n",config.chempot);
	fprintf(out,"# Initial diagram length: %f\n",config.endtau);
	fprintf(out,"# Max diagram length: %f\n",config.maxtau);
	fprintf(out,"# Potential parameters: (logn) = (%f)\n",log(config.n));
	fprintf(out,"#\n");

	total_proposed=total_accepted=total_rejected=0;

	fprintf(out,"# Update statistics:\n");

	for(d=0;d<DIAGRAM_NR_UPDATES;d++)
	{
		fprintf(out,"# Update #%d (%s): ",d,update_names[d]);
		show_update_statistics(out,proposed[d],accepted[d],rejected[d]);
	
		total_proposed+=proposed[d];
		total_accepted+=accepted[d];
		total_rejected+=rejected[d];
	}

	fprintf(out,"# Total: ");
	show_update_statistics(out,total_proposed,total_accepted,total_rejected);
	fprintf(out,"#\n# Average order: %f\n",((double)(avgorder[0]))/((double)(avgorder[1])));
	fprintf(out,"#\n");
	fprintf(out,"# Sampled quantity: Green's function (G)\n");
	fprintf(out,"# Iterations: %d\n",config.iterations);
	fprintf(out,"# Nr of bins: %d\n",config.bins);
	fprintf(out,"# Bin width: %f\n",config.width);
	fprintf(out,"# (last bin is overflow)\n");
	fprintf(out,"#\n");

	fprintf(out,"# <Bin center> <Average> <Stddev>\n");

	for(d=0;d<=config.bins;d++)
	{
		fprintf(out,"%f %f ",config.width*d+config.width/2.0f,histogram_get_bin_average(ht,d));
		
		if(ht->sctx->smpls[d]->next>1)
			fprintf(out,"%f\n",sqrtf(histogram_get_bin_variance(ht,d)/(ht->sctx->smpls[d]->next)));
		else
			fprintf(out,"%f\n",0.0f);
	}

	/*
		That's all folks. Cleaning up.
	*/

	if(out)
		fclose(out);

	fini_histogram(ht);

	return 0;
}
