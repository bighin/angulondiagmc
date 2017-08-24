#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <unistd.h>

#include "diagrams.h"
#include "updates.h"
#include "stat.h"
#include "aux.h"
#include "graphs.h"
#include "debug.h"

#include "inih/ini.h"
#include "libprogressbar/progressbar.h"

/*
	The following headers need to be included after graphs.h, becuase it defines a 'lines' macro
	which conflicts with a field in struct graph_t
*/

#include <curses.h>
#include <term.h>

static char term_buffer[2048];

void init_terminal_data(void)
{
	char *termtype=getenv("TERM");
	int success;

	if (termtype==NULL)
		fprintf(stderr,"Please specify a terminal type with 'setenv TERM <yourtype>'.\n");

	success=tgetent(term_buffer,termtype);
	
	if(success<0)
		fprintf(stderr, "Could not access the termcap data base.\n");

	if(success==0)
		fprintf(stderr, "Terminal type `%s' is not defined.\n",termtype);
}

struct configuration_t
{
	/* "general" section */

	char *prefix;
	bool seedrng;
	bool progressbar;
	bool animate;

	/* "parameters" section */

	int j;
	int m;
	double endtau;
	double chempot;
	double maxtau;
	int maxorder;

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

#define DO_RECOUPLE

int update_add_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg)
{
	int lambda,mu;
	double k,tau1,tau2,acceptance_ratio,oldweight;
	bool is_accepted;

	struct arc_t *thisline;
	struct vertex_t *v1,*v2;
	struct free_propagators_ctx_t fpc;

	double omegas[3]={cfg->omega0,cfg->omega1,cfg->omega2};

	/*
		Do we want to limit the simulation only to first order diagrams?
	*/

	if(get_nr_phonons(dgr)>=cfg->maxorder)
		return UPDATE_UNPHYSICAL;

	oldweight=diagram_weight(dgr);

	/*
		The new lambda is chosen using a uniform distribution
	*/

	lambda=gsl_rng_uniform_int(dgr->rng_ctx,3);
	mu=0;

	/*
		The start time tau1 is sampled uniformly, whereas tau2 is sampled from
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

	if((check_triangle_condition_and_parity(dgr,v1)==false)||(check_triangle_condition_and_parity(dgr,v2)==false))
	{
		int lastline=get_nr_phonons(dgr)-1;
		
		diagram_remove_phonon_line(dgr,lastline);

		/*
			The weight went to zero due to bad couplings, so we cannot trust
			the incremental update...
		*/

		dgr->weight=oldweight;

		return UPDATE_UNPHYSICAL;
	}

	save_free_propagators(dgr,&fpc,thisline->startmidpoint,thisline->endmidpoint);

#ifdef DO_RECOUPLE
	recouple(dgr,thisline->startmidpoint,thisline->endmidpoint);
#endif

	acceptance_ratio=diagram_weight(dgr)/oldweight;	
	acceptance_ratio*=3.0f*dgr->endtau;
	acceptance_ratio/=get_nr_phonons(dgr)+1;
	acceptance_ratio/=doubly_truncated_exp_pdf(dgr->rng_ctx,omegas[lambda],tau1,dgr->endtau,tau2);

	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
		int lastline;

		/*
			Restoring the free propagators effectively
			undoes the recouple() call.
		*/

#ifdef DO_RECOUPLE
		restore_free_propagators(dgr,&fpc);
#endif

		lastline=get_nr_phonons(dgr)-1;
		diagram_remove_phonon_line(dgr,lastline);

		assert(fabs(diagram_weight(dgr)-oldweight)<1e-7*oldweight);

		return UPDATE_REJECTED;
	}

#ifdef DO_RECOUPLE
	release_free_propagators_ctx(&fpc);
#endif

	return UPDATE_ACCEPTED;
}

int update_remove_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg)
{
	struct arc_t *arc;
	struct free_propagators_ctx_t fpc;

	int target,lambda,mu,nr_phonons,startmidpoint,endmidpoint;
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
	startmidpoint=arc->startmidpoint;
	endmidpoint=arc->endmidpoint;

	/*
		Note the order of the following two function calls and
		of the subsequent recouple() call!
	*/

	diagram_check_consistency(dgr);

	save_free_propagators(dgr,&fpc,startmidpoint,endmidpoint);

	/*
		A line removal can fail, since it may leave the remaining
		lines in a unphysical state.
	*/
	
	diagram_check_consistency(dgr);
	
	if(diagram_remove_phonon_line(dgr,target)==false)
	{	
		diagram_add_phonon_line(dgr,tau1,tau2,k,lambda,mu);
		restore_free_propagators(dgr,&fpc);

		/*
			The update failed, so this means that the incremental
			weight went to zero, and we need to fix it...
		*/

		dgr->weight=oldweight;
	
		return UPDATE_UNPHYSICAL;
	}

	/*
		This is a bit tricky: if we just removed a bubble,
		then there are no free propagators to recouple.
	*/

	diagram_check_consistency(dgr);

#ifdef DO_RECOUPLE
	if((startmidpoint+1)!=endmidpoint)
		recouple(dgr,startmidpoint,endmidpoint-2);
#endif

	diagram_check_consistency(dgr);

	acceptance_ratio=diagram_weight(dgr)/oldweight;
	acceptance_ratio/=3.0f*dgr->endtau;
	acceptance_ratio*=nr_phonons;
	acceptance_ratio*=doubly_truncated_exp_pdf(dgr->rng_ctx,omegas[lambda],tau1,dgr->endtau,tau2);

	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
		diagram_check_consistency(dgr);

		diagram_add_phonon_line(dgr,tau1,tau2,k,lambda,mu);

		diagram_check_consistency(dgr);

		/*
			Nasty, nasty, nasty: diagram_add_phonon_line() doesn't 
		*/

		restore_free_propagators(dgr,&fpc);
		
		assert(fabs(diagram_weight(dgr)-oldweight)<1e-7*oldweight);
	
		return UPDATE_REJECTED;
	}

	release_free_propagators_ctx(&fpc);

	return UPDATE_ACCEPTED;
}

void load_config_defaults(struct configuration_t *config)
{
	config->prefix=strdup("default");
	config->seedrng=false;
	config->progressbar=true;
	config->animate=false;

	config->j=1;
	config->m=0;
	config->endtau=1.0f;
	config->chempot=0.95f;
	config->maxtau=100.0f;
	config->maxorder=1024;
	
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
	else if(MATCH("general","animate"))
	{
		if(!strcmp(value,"true"))
			pconfig->animate=true;
		else if(!strcmp(value,"false"))
			pconfig->animate=false;
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
	else if(MATCH("parameters","maxorder"))
	{
		pconfig->maxorder=atoi(value);
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

	for(c=0;c<DIAGRAM_NR_UPDATES;c++)
		proposed[c]=accepted[c]=rejected[c]=0;
	
	avgorder[0]=avgorder[1]=0;

	/*
		We load some sensible defaults in case they are not in the .ini file, the we try
		loading a .ini file.
	*/

	load_config_defaults(&config);

	if(ini_parse(configfile,configuration_handler,&config)<0)
	{
		fprintf(stderr,"Couldn't read or parse '%s'\n",configfile);
		exit(0);
	}

	if((config.animate==true)&&(config.progressbar==true))
	{
		fprintf(stderr,"The options 'animate' and 'progressbar' cannot be set at the same time!\n");
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
	dcfg.omega1=config.omega1;
	dcfg.omega2=config.omega2;

	/*
		We initialize some basic data structures and various other stuff, as well as the
		random number generator, in case a NON deterministic seed is requested.
	*/

	dgr=init_diagram(&dcfg);
	ht=init_histogram(config.bins,config.width);

	if(config.seedrng)
		seed_rng(dgr->rng_ctx);

	if(config.progressbar)
		progress=progressbar_new("Progress",config.iterations/16384);
	else
		progress=NULL;

        init_terminal_data();

	/*
		This is the main DiagMC loop
	*/

	for(c=0;c<config.iterations;c++)
	{
		int update_type,status;

		update_type=gsl_rng_uniform_int(dgr->rng_ctx,DIAGRAM_NR_UPDATES);
		update_type=gsl_rng_uniform_int(dgr->rng_ctx,3);

		status=updates[update_type](dgr,&config);

		if((config.animate)&&(status==UPDATE_ACCEPTED)&&((update_type==1)||(update_type==2)))
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

		switch(status)
		{
			case UPDATE_ACCEPTED:
			proposed[update_type]++;
			accepted[update_type]++;
			break;

			/*
				FIXME

				Here unphysical and rejected updates are treated in exactly the same
				way. Is this OK?
			*/

			case UPDATE_UNPHYSICAL:
			case UPDATE_REJECTED:
			proposed[update_type]++;
			rejected[update_type]++;
			break;

			case UPDATE_ERROR:
			assert(false);
		}

		assert(fabs(diagram_weight(dgr)-diagram_weight_non_incremental(dgr))<10e-7*diagram_weight(dgr));
		diagram_check_consistency(dgr);

		histogram_add_sample(ht,diagram_weight(dgr)*diagram_m_weight(dgr),dgr->endtau);

		if((config.progressbar)&&((c%16384)==0))
			progressbar_inc(progress);

		avgorder[0]+=get_nr_phonons(dgr);
		avgorder[1]++;
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
	fprintf(out,"# Initial and final state: (j=%d, m=%d)\n",config.j,config.m);
	fprintf(out,"# Chemical potential: %f\n",config.chempot);
	fprintf(out,"# Initial diagram length: %f\n",config.endtau);
	fprintf(out,"# Max diagram length: %f\n",config.maxtau);
	fprintf(out,"# Potential parameters: (c0, c1, c2, omega0, omega1, omega2) = (%f, %f, %f, %f, %f, %f)\n",config.c0,config.c1,config.c2,config.omega0,config.omega1,config.omega2);
	fprintf(out,"#\n");

	total_proposed=total_accepted=total_rejected=0;

	fprintf(out,"# Update statistics:\n");

	for(c=0;c<DIAGRAM_NR_UPDATES;c++)
	{
		fprintf(out,"# Update #%d (%s): ",c,update_names[c]);
		show_update_statistics(out,proposed[c],accepted[c],rejected[c]);
	
		total_proposed+=proposed[c];
		total_accepted+=accepted[c];
		total_rejected+=rejected[c];
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

	for(c=0;c<=config.bins;c++)
	{
		fprintf(out,"%f %f ",config.width*c+config.width/2.0f,histogram_get_bin_average(ht,c));
		
		if(ht->sctx->smpls[c]->next>1)
			fprintf(out,"%f\n",sqrtf(histogram_get_bin_variance(ht,c)/(ht->sctx->smpls[c]->next)));
		else
			fprintf(out,"%f\n",0.0f);
	}

	/*
		That's all folks. Cleaning up.
	*/

	if(out)
		fclose(out);

	fini_histogram(ht);
	fini_diagram(dgr);

	return 0;
}
