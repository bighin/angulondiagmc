#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <unistd.h>
#include <signal.h>
#include <sys/time.h>

#include "mc.h"
#include "physics.h"
#include "diagrams.h"
#include "updates.h"
#include "stat.h"
#include "aux.h"
#include "graphs.h"
#include "phonon.h"
#include "debug.h"
#include "config.h"
#include "histograms.h"
#include "inih/ini.h"
#include "libprogressbar/progressbar.h"
#include "gnuplot_i/gnuplot_i.h"

int get_nr_available_phonons(struct diagram_t *dgr)
{
	int c,available;

	for(c=available=0;c<get_nr_phonons(dgr);c++)
	{
		int startmidpoint,endmidpoint;
		
		startmidpoint=get_phonon_line(dgr,c)->startmidpoint;
		endmidpoint=get_phonon_line(dgr,c)->endmidpoint;
		
		if(deltaj(dgr,startmidpoint)==-deltaj(dgr,endmidpoint))
			available++;
	}
	
	return available;
}

int get_random_available_phonon(struct diagram_t *dgr)
{
	int c,available,target;

	available=0;
	target=gsl_rng_uniform_int(dgr->rng_ctx,get_nr_available_phonons(dgr));

	for(c=0;c<get_nr_phonons(dgr);c++)
	{
		int startmidpoint,endmidpoint;
		
		startmidpoint=get_phonon_line(dgr,c)->startmidpoint;
		endmidpoint=get_phonon_line(dgr,c)->endmidpoint;
		
		if(deltaj(dgr,startmidpoint)==-deltaj(dgr,endmidpoint))
		{
			if(available==target)
				return c;
			
			available++;
		}
	}

	assert(false);
	
	return -1;
}

/*
	Now it's time to start with all the different updates!
*/

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
		This update is always accepted.
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

		assert(almost_same_float(diagram_weight(dgr),oldweight));
	
		return UPDATE_REJECTED;
	}
	
	return UPDATE_ACCEPTED;
}

#define MAXLAMBDA	(1)

int update_add_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg)
{
	int lambda,mu,deltaj1,deltaj2;
	double tau1,tau2,weightratio,acceptance_ratio;
	bool is_accepted;

	struct arc_t *thisline;

#ifndef NDEBUG
	double oldweight=diagram_weight(dgr);
#endif

	/*
		Do we want to limit the simulation only to first order diagrams (or at a certain maximum order)?
	*/

	if(get_nr_phonons(dgr)>=cfg->maxorder)
		return UPDATE_UNPHYSICAL;

	/*
		The new lambda and mu are chosen using a uniform distribution.
	
		Deltaj_1 and Delta_j are chosen as to respect angular momentum conservation
		and parity, i.e. among the following integer values
	
		-lambda, -lambda + 2, ..., lambda - 2, lambda
	*/

	lambda=gsl_rng_uniform_int(dgr->rng_ctx,1+MAXLAMBDA);
	mu=gsl_rng_uniform_int(dgr->rng_ctx,1+2*lambda)-lambda;
	deltaj1=2*gsl_rng_uniform_int(dgr->rng_ctx,1+lambda)-lambda;
	deltaj2=-deltaj1;

	/*
		The start time tau1 is sampled uniformly in a randomly chosen propagator,
		whereas tau2 is sampled from a truncated exponential distribution.
	*/

	tau1=dgr->endtau*gsl_rng_uniform_pos(dgr->rng_ctx);
	tau2=tau1+phonon_dist(dgr->rng_ctx,dgr->phonon_ctx,lambda);

	if(tau1>tau2)
	{
		fprintf(stderr,"Unexpected negative return value from phonon_dist(). Please debug me! (%f %f)\n",tau1,tau2);
		exit(0);
	}

	if(tau2>=dgr->endtau)
		return UPDATE_UNPHYSICAL;

	/*
		We save the current diagram and the line is added...
	*/

	weightratio=1.0f/calculate_propagators_and_vertices(dgr,0,get_nr_vertices(dgr)-1);

	diagram_add_phonon_line(dgr,tau1,tau2,lambda,mu);

	/*
		...and finally we fix the \Delta_j's at every vertex
	*/

	thisline=get_phonon_line(dgr,get_nr_phonons(dgr)-1);
	change_deltaj(dgr,thisline->endmidpoint,deltaj2);
	change_deltaj(dgr,thisline->startmidpoint,deltaj1);
	45
	if(propagators_are_allowed(dgr)==false)
	{
		change_deltaj(dgr,thisline->endmidpoint,0);
		change_deltaj(dgr,thisline->startmidpoint,0);
		diagram_remove_phonon_line(dgr,get_nr_phonons(dgr)-1);

#ifndef NDEBUG
		if(isfinite(oldweight)&&isfinite(diagram_weight(dgr)))
			assert(almost_same_float(oldweight,diagram_weight(dgr))==true);
#endif

		return UPDATE_REJECTED;
	}
	
	/*
		We calculate the weight ratio (the fundamental quantity in accepting/rejecting the update)
		and we use it to update the diagram weight, as well.
	*/

	weightratio*=calculate_arc_weight(dgr,thisline);
	weightratio*=calculate_propagators_and_vertices(dgr,0,get_nr_vertices(dgr)-1);

	if(weightratio<0.0f)
		dgr->sign*=-1;

#ifndef NDEBUG
	if((!isinf(diagram_weight(dgr)))&&(!isinf(oldweight)))
		assert(almost_same_float(weightratio,diagram_weight(dgr)/oldweight));
#endif

	/*
		Finally we calculate the acceptance ratio for the update.
	*/

	acceptance_ratio=fabs(weightratio);
	acceptance_ratio/=1.0f/(1+MAXLAMBDA);
	acceptance_ratio/=1.0f/(1+2*lambda);
	acceptance_ratio/=1.0f/(1+lambda);
	acceptance_ratio/=1.0f/dgr->endtau;
	acceptance_ratio/=phonon_pdf(dgr->phonon_ctx,lambda,tau2-tau1);
	acceptance_ratio*=1.0f/get_nr_available_phonons(dgr);

	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
		change_deltaj(dgr,thisline->endmidpoint,0);
		change_deltaj(dgr,thisline->startmidpoint,0);
		diagram_remove_phonon_line(dgr,get_nr_phonons(dgr)-1);

		if(weightratio<0.0f)
			dgr->sign*=-1;

#ifndef NDEBUG
		if(isfinite(oldweight)&&isfinite(diagram_weight(dgr)))
			assert(almost_same_float(oldweight,diagram_weight(dgr))==true);
#endif

		return UPDATE_REJECTED;
	}

	return UPDATE_ACCEPTED;
}

int update_remove_phonon_line(struct diagram_t *dgr,struct configuration_t *cfg)
{
	struct arc_t *arc;

	int target,lambda,mu,nr_available_phonons,startmidpoint,endmidpoint,deltaj1,deltaj2;
	double weightratio,targetweight,acceptance_ratio,tau1,tau2;
	bool is_accepted;

#ifndef NDEBUG
	double oldweight=diagram_weight(dgr);
#endif

	nr_available_phonons=get_nr_available_phonons(dgr);

	if(nr_available_phonons<=0)
		return UPDATE_UNPHYSICAL;

	target=get_random_available_phonon(dgr);
	arc=get_phonon_line(dgr,target);
	targetweight=calculate_arc_weight(dgr,arc);

	tau1=arc->starttau;
	tau2=arc->endtau;
	lambda=arc->lambda;
	mu=arc->mu;
	startmidpoint=arc->startmidpoint;
	endmidpoint=arc->endmidpoint;
	deltaj1=deltaj(dgr,startmidpoint);
	deltaj2=deltaj(dgr,endmidpoint);

	weightratio=1.0f/calculate_propagators_and_vertices(dgr,0,get_nr_vertices(dgr)-1);

	change_deltaj(dgr,endmidpoint,0);
	change_deltaj(dgr,startmidpoint,0);
	diagram_remove_phonon_line(dgr,target);

	if(propagators_are_allowed(dgr)==false)
	{
		diagram_add_phonon_line(dgr,tau1,tau2,lambda,mu);
		change_deltaj(dgr,endmidpoint,deltaj2);
		change_deltaj(dgr,startmidpoint,deltaj1);

#ifndef NDEBUG
		if(isfinite(oldweight)&&isfinite(diagram_weight(dgr)))
			assert(almost_same_float(oldweight,diagram_weight(dgr))==true);
#endif

		return UPDATE_REJECTED;
	}

	weightratio/=targetweight;
	weightratio*=calculate_propagators_and_vertices(dgr,0,get_nr_vertices(dgr)-1);

	if(weightratio<0.0f)
		dgr->sign*=-1;

#ifndef NDEBUG
	if((!isinf(diagram_weight(dgr)))&&(!isinf(oldweight)))
		assert(almost_same_float(weightratio,diagram_weight(dgr)/oldweight));
#endif

	/*
		Finally we calculate the acceptance ratio for the update.
	*/

	acceptance_ratio=fabs(weightratio);
	acceptance_ratio*=1.0f/(1+MAXLAMBDA);
	acceptance_ratio*=1.0f/(1+2*lambda);
	acceptance_ratio*=1.0f/(1+lambda);
	acceptance_ratio*=1.0f/dgr->endtau;
	acceptance_ratio*=phonon_pdf(dgr->phonon_ctx,lambda,tau2-tau1);
	acceptance_ratio/=1.0f/nr_available_phonons;

	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
		diagram_add_phonon_line(dgr,tau1,tau2,lambda,mu);
		change_deltaj(dgr,endmidpoint,deltaj2);
		change_deltaj(dgr,startmidpoint,deltaj1);

		if(weightratio<0.0f)
			dgr->sign*=-1;

#ifndef NDEBUG
		if(isfinite(oldweight)&&isfinite(diagram_weight(dgr)))
			assert(almost_same_float(oldweight,diagram_weight(dgr))==true);
#endif

		return UPDATE_REJECTED;
	}

	return UPDATE_ACCEPTED;
}

void identify_lambda1_clusters(struct diagram_t *dgr,int *clusters,int *iclusters,int maxnrclusters)
{
	int nr_vertices=get_nr_vertices(dgr);

	bool in_cluster=false;
	bool cluster_has_more_than_one_phonon=false;
	int start=-1,end=-1;

	*iclusters=0;

	for(int c=0;c<nr_vertices;c++)
	{
		struct vertex_t *thisvertex=get_vertex(dgr,c);
		
		/*
			We are entering a new cluster
		*/
		
		if((thisvertex->left->arcs_over_me_lambda1==0)&&(thisvertex->right->arcs_over_me_lambda1>0))
		{
			in_cluster=true;
			start=c;
		}

		/*
			We are leaving a cluster
		*/

		if((thisvertex->left->arcs_over_me_lambda1>0)&&(thisvertex->right->arcs_over_me_lambda1==0))
		{
			assert(in_cluster==true);

			end=c;

#if 0
			printf("Found cluster starting and %d and ending at %d. ",start,end);

			if(cluster_has_more_than_one_phonon==true)
				printf("The cluster has more than one phonon.\n");
			else
				printf("The cluster has one phonon.\n");
#endif

			if(cluster_has_more_than_one_phonon==true)
			{
				clusters[(*iclusters)*2]=start;
				clusters[(*iclusters)*2+1]=end;
				(*iclusters)++;
			
				assert((*iclusters)<maxnrclusters);
			}

			in_cluster=false;
			cluster_has_more_than_one_phonon=false;
		}

		/*
			The cluster we are in has more than one phonon inside.
		*/

		if((in_cluster==true)&&(thisvertex->right->arcs_over_me_lambda1>1))
			cluster_has_more_than_one_phonon=true;
	}
}

int update_shuffle(struct diagram_t *dgr,struct configuration_t *cfg)
{
	int nr_vertices;
	double weightratio,acceptance_ratio;
	bool is_accepted;

#ifndef NDEBUG
	double oldweight=diagram_weight(dgr);
#endif

#define MAX_NR_CLUSTERS	(1024)

	int clusters[2*MAX_NR_CLUSTERS];
	int iclusters;
	int target,targetstart,targetend;

	struct free_propagators_ctx_t *fpc;

	nr_vertices=get_nr_vertices(dgr);
	
	if(nr_vertices<1)
		return UPDATE_UNPHYSICAL;

	assert(MAX_NR_CLUSTERS>=cfg->maxorder);
	identify_lambda1_clusters(dgr,clusters,&iclusters,MAX_NR_CLUSTERS);

#warning Debug this function! Use the check below and other checks!

#if 0
	{
		if((nr_vertices>20)&&(iclusters>1))
		{
			print_diagram(dgr,PRINT_TOPOLOGY|PRINT_PROPAGATORS);
		
			for(int d=0;d<iclusters;d++)
				printf("%d %d\n",clusters[2*d],clusters[2*d+1]);
	
			exit(0);
		}
	}
#endif

	if(iclusters==0)
		return UPDATE_UNPHYSICAL;

	target=gsl_rng_uniform_int(dgr->rng_ctx,iclusters);
	targetstart=clusters[2*target];
	targetend=clusters[2*target+1];

	assert(targetstart>=0);
	assert(targetend>=0);
	assert(targetstart<nr_vertices);
	assert(targetend<nr_vertices);
	assert(targetend>targetstart);

	fpc=malloc(sizeof(struct free_propagators_ctx_t));
	assert(fpc!=NULL);

	save_free_propagators(dgr,fpc,targetstart,targetend);

	weightratio=1.0f;
	weightratio/=calculate_propagators_and_vertices(dgr,targetstart,targetend);

	recouple(dgr,targetstart,targetend);

	weightratio*=calculate_propagators_and_vertices(dgr,targetstart,targetend);

	acceptance_ratio=weightratio;

	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
		restore_free_propagators(dgr,fpc);

#ifndef NDEBUG
		if(isfinite(oldweight)&&isfinite(diagram_weight(dgr)))
			assert(almost_same_float(oldweight,diagram_weight(dgr))==true);
#endif
		if(fpc)
			free(fpc);

		return UPDATE_REJECTED;
	}

	release_free_propagators_ctx(fpc);

	if(fpc)
		free(fpc);

	return UPDATE_ACCEPTED;
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

void send_to_gnuplot(gnuplot_ctrl *gp,gsl_histogram *g,struct configuration_t *config)
{
	double *x,*y1,*y2;
	char desc1[1024],desc2[1024];
	int c;
	
	gnuplot_resetplot(gp);
	gnuplot_setstyle(gp,"points");

	if(!(x=malloc(sizeof(double)*config->bins)))
		return;

	if(!(y1=malloc(sizeof(double)*config->bins)))
	{
		if(x)	free(x);
		
		return;
	}

	if(!(y2=malloc(sizeof(double)*config->bins)))
	{
		if(x)	free(x);
		if(y1)	free(y1);

		return;
	}

	for(c=0;c<config->bins;c++)
	{
		double tau,Ej;
	
		tau=config->width*c+config->width/2.0f;
		Ej=config->j*(config->j+1);

		x[c]=tau;
		y1[c]=gsl_histogram_get(g,c)/gsl_histogram_get(g,0)*exp(-(Ej-config->chempot)*config->width/2.0f);
		y2[c]=exp(-(Ej-config->chempot)*tau);
	}

	for(c=0;c<config->bins;c++)
	{
		y1[c]=log(y1[c]);
		y2[c]=log(y2[c]);
	}

	snprintf(desc1,1024,"Angulon (j=%d, logn=%f)",config->j,log(config->n));
	snprintf(desc2,1024,"Free rotor (j=%d)",config->j);

	gnuplot_cmd(gp,"set xlabel \"{/Symbol t}\" font \" ,30\"");
	gnuplot_cmd(gp,"set ylabel \"log(G({/Symbol t}))\" font \" ,30\"");
	gnuplot_cmd(gp,"set xtics font \" ,22\"");
	gnuplot_cmd(gp,"set ytics font \" ,22\"");
	gnuplot_cmd(gp,"set key font \" ,22\"");
	gnuplot_cmd(gp,"set pointsize 1.5");

	gnuplot_plot_xyy(gp,x,y1,y2,config->bins,desc1,desc2);

	if(x)	free(x);
	if(y1)	free(y1);
	if(y2)	free(y2);
}

static volatile int keep_running=1;

void interrupt_handler(int dummy)
{
	keep_running=0;
}

int do_diagmc(struct configuration_t *config,struct mc_output_data_t *output)
{
	struct diagram_t *dgr;
	struct diagram_parameters_t dpars;

	gsl_histogram *g,*g0,*g1,*g2;

#define NR_BLOCKSIZES	(5)

	int blocksizes[NR_BLOCKSIZES]={500,1000,10000,25000,50000};
	struct super_sampler_t *sst;

	int c,d;
	FILE *out;
	char fname[1024],progressbar_text[1024];

	progressbar *progress;
	gnuplot_ctrl *gp;

	struct timeval starttime;

#define DIAGRAM_NR_UPDATES	(4)

	int (*updates[DIAGRAM_NR_UPDATES])(struct diagram_t *,struct configuration_t *);
	char *update_names[DIAGRAM_NR_UPDATES];

	int proposed[DIAGRAM_NR_UPDATES],accepted[DIAGRAM_NR_UPDATES],rejected[DIAGRAM_NR_UPDATES];
	int total_proposed,total_accepted,total_rejected;

	long long int avgorder[2];
	long long int physical_updates,unphysical_updates;
	long long int g0stats=0,gstats=0;
	double avglength[2];

	/*
		We set up the updates we will be using
	*/

	updates[0]=update_length;
	updates[1]=update_add_phonon_line;
	updates[2]=update_remove_phonon_line;
	updates[3]=update_shuffle;

	update_names[0]="UpdateLength";
	update_names[1]="AddPhononLine";
	update_names[2]="RemovePhononLine";
	update_names[3]="Shuffle";

	/*
		We reset the update statistics
	*/

	for(d=0;d<DIAGRAM_NR_UPDATES;d++)
		proposed[d]=accepted[d]=rejected[d]=0;
	
	physical_updates=unphysical_updates=0;
	
	avgorder[0]=avgorder[1]=0;
	avglength[0]=avglength[1]=0.0f;

	/*
		We print some informative message, and then we open the log file
	*/

	fprintf(stderr,"Performing %d iterations, %d thermalization sweeps, using %d thread(s)\n",config->iterations,config->thermalization,config->nthreads);

	snprintf(fname,1024,"%s.dat",config->prefix);
	if(!(out=fopen(fname,"w+")))
	{
		fprintf(stderr,"Error: couldn't open %s for writing\n",fname);
		return 0;
	}

	fprintf(stderr,"Writing results to '%s'\n",fname);

	/*
		The diagram parameters are loaded from the configuration, and then a new diagram is created
	*/

	dpars.j=config->j;
	dpars.endtau=config->endtau;
	dpars.maxtau=config->maxtau;
	dpars.chempot=config->chempot;
	dpars.n=config->n;

	dgr=init_diagram(&dpars,true);

	/*
		We setup an interript handler to gracefully handle a CTRL-C, and initialize a structure needed
		by the ncurses library to return info about the current terminal.
	*/

	keep_running=1;
	signal(SIGINT,interrupt_handler);
        init_terminal_data();

	if(config->progressbar==true)
		progress=progressbar_new("Progress",config->iterations/16384);
	else
		progress=NULL;

	if(config->liveplot==true)
		gp=gnuplot_init();
	else
		gp=NULL;

	/*
		We save the start time
	*/

	gettimeofday(&starttime,NULL);

	/*
		We initialize the histograms
	*/

	g=gsl_histogram_alloc(config->bins);
	g0=gsl_histogram_alloc(config->bins);
	g1=gsl_histogram_alloc(config->bins);
	g2=gsl_histogram_alloc(config->bins);

	gsl_histogram_set_ranges_uniform(g,0.0f,config->bins*config->width);
	gsl_histogram_set_ranges_uniform(g0,0.0f,config->bins*config->width);
	gsl_histogram_set_ranges_uniform(g1,0.0f,config->bins*config->width);
	gsl_histogram_set_ranges_uniform(g2,0.0f,config->bins*config->width);

	sst=init_super_sampler(blocksizes,NR_BLOCKSIZES,config->bins,config->width);

	/*
		We initialize the random number generator, and we seed it in case a NON deterministic seed is requested.
	*/

	dgr->rng_ctx=gsl_rng_alloc(gsl_rng_mt19937);
	assert(dgr->rng_ctx!=NULL);

	if(config->seedrng)
		seed_rng(dgr->rng_ctx);

	/*
		This is the main DiagMC loop
	*/

	for(c=0;(c<config->iterations)&&(keep_running==1);c++)
	{
		int update_type,status;

		update_type=gsl_rng_uniform_int(dgr->rng_ctx,DIAGRAM_NR_UPDATES);
		status=updates[update_type](dgr,config);

		if((config->nthreads==1)&&(config->animate)&&(status==UPDATE_ACCEPTED)&&((update_type==1)||(update_type==2)))
		{
			if((c%16)==0)
			{
				int d,plines,flags;

				flags=PRINT_TOPOLOGY|PRINT_PROPAGATORS|PRINT_INFO0;
				plines=print_diagram(dgr,flags|PRINT_DRYRUN);
	
				for(d=plines;d<get_terminal_nr_lines();d++)
						printf("\n");

				print_diagram(dgr,flags);
				usleep(5000);
			}
		}

		proposed[update_type]++;

		switch(status)
		{
			case UPDATE_ACCEPTED:

			accepted[update_type]++;

			break;

			/*
				FIXME

				Here unphysical and rejected updates are treated in exactly the same
				way. Is this OK?
			*/

			case UPDATE_UNPHYSICAL:
			case UPDATE_REJECTED:

			rejected[update_type]++;

			break;

			case UPDATE_ERROR:
			assert(false);
		}

#ifndef NDEBUG
		diagram_check_consistency(dgr);
		assert(angular_momentum_is_conserved(dgr)==true);
#endif

		if(configuration_is_physical(dgr)==true)
		{	
			if((c>config->thermalization)&&((physical_updates%config->decorrelation)==0))
			{
				super_sampler_add_time_sample(sst,dgr->endtau,dgr->sign);
				
				gsl_histogram_accumulate(g,dgr->endtau,dgr->sign);

				gstats++;

				if(get_nr_phonons(dgr)==0)
					g0stats++;

				if(get_nr_phonons(dgr)==0)
					gsl_histogram_accumulate(g0,dgr->endtau,dgr->sign);

				if(get_nr_phonons(dgr)==1)
					gsl_histogram_accumulate(g1,dgr->endtau,dgr->sign);

				if(get_nr_phonons(dgr)==2)
					gsl_histogram_accumulate(g2,dgr->endtau,dgr->sign);
			}

			physical_updates++;
		}
		else
		{
			unphysical_updates++;
		}

		avgorder[0]+=get_nr_phonons(dgr);
		avgorder[1]++;

		avglength[0]+=dgr->endtau;
		avglength[1]++;

		if((c%16384)==0)
		{
			struct timeval now;
			double elapsedtime;
			
			if(config->progressbar)
			{
				double avg1,avg2,ratio;

				avg1=((double)(avgorder[0]))/((double)(avgorder[1]));
				avg2=((double)(avglength[0]))/((double)(avglength[1]));

				ratio=((double)(physical_updates))/((double)(physical_updates+unphysical_updates));

				snprintf(progressbar_text,1024,"Progress (avg. order: %f, avg. length: %f, P/(P+U): %f)",avg1,avg2,ratio);

				progressbar_update_label(progress,progressbar_text);
				progressbar_inc(progress);
			}

			if(config->liveplot)
				send_to_gnuplot(gp,g,config);	
		
			gettimeofday(&now,NULL);

			elapsedtime=(now.tv_sec-starttime.tv_sec)*1000.0;
			elapsedtime+=(now.tv_usec-starttime.tv_usec)/1000.0;
			elapsedtime/=1000;
		
			if((config->timelimit>0.0f)&&(elapsedtime>config->timelimit))
				keep_running=0;
		}
	}

	if(keep_running==0)
	{
		printf("Caught SIGINT, exiting earlier.\n");
		config->iterations=c;
	}

	if(config->progressbar)
		progressbar_finish(progress);

	if(dgr->rng_ctx)
		gsl_rng_free(dgr->rng_ctx);

	fini_diagram(dgr);

	/*
		Now we print the statistics we collected to the output file in a nice way.
	*/

	fprintf(out,"# Diagrammatic Monte Carlo for the angulon\n");
	fprintf(out,"#\n");
	fprintf(out,"# Configuration loaded from '%s'\n",config->ininame);
	fprintf(out,"# Output file is '%s'\n",fname);
	fprintf(out,"#\n");
	fprintf(out,"# Initial and final state: (j=%d, sampling over all possible m values)\n",config->j);
	fprintf(out,"# Chemical potential: %f\n",config->chempot);
	fprintf(out,"# Initial diagram length: %f\n",config->endtau);
	fprintf(out,"# Max diagram length: %f\n",config->maxtau);
	fprintf(out,"# Potential parameters: (logn) = (%f)\n",log(config->n));
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
	fprintf(out,"# Average length: %f\n",((double)(avglength[0]))/((double)(avglength[1])));
	fprintf(out,"#\n");
	fprintf(out,"# Sampled quantity: Green's function (G)\n");
	fprintf(out,"# Iterations: %d\n",config->iterations);
	fprintf(out,"# Thermalization sweeps: %d\n",config->thermalization);

	/*
		FIXME Poco elegante!
	*/
	
	{
		struct timeval now;
		double elapsedtime;

		gettimeofday(&now,NULL);

		elapsedtime=(now.tv_sec-starttime.tv_sec)*1000.0;
		elapsedtime+=(now.tv_usec-starttime.tv_usec)/1000.0;
		elapsedtime/=1000;

		fprintf(out,"# Total time: %f seconds\n",elapsedtime);
	}
	
	fprintf(out,"# Nr of bins: %d\n",config->bins);
	fprintf(out,"# Bin width: %f\n",config->width);
	fprintf(out,"# (last bin is overflow)\n");
	fprintf(out,"#\n");

	fprintf(out,"# Blocksizes: ");
	
	for(int k=0;k<NR_BLOCKSIZES;k++)
	{
		fprintf(out,"%d",blocksizes[k]);
		
		if((k+1)!=NR_BLOCKSIZES)
			fprintf(out,", ");
	}
	
	fprintf(out,"\n#\n");
	fprintf(out,"# <Bin center> <G> <Sigma_G blocksize[0]> ... <Sigma_G blocksize[N]> <G0>\n");

	/*
		We normalize the histogram...
	*/

	{
		double lower,upper,bin0center,Ej,scalefactor;

		gsl_histogram_get_range(g,0,&lower,&upper);

		bin0center=(lower+upper)/2.0f;
		Ej=config->j*(config->j+1)-config->chempot;
		scalefactor=exp(-Ej*bin0center)/gsl_histogram_get(g,0);

		gsl_histogram_scale(g,scalefactor);
		gsl_histogram_scale(g0,scalefactor);
		gsl_histogram_scale(g1,scalefactor);
		gsl_histogram_scale(g2,scalefactor);
	}

	/*
		...and then we print it
	*/

	for(d=0;d<config->bins;d++)
	{
		double lower,upper,bincenter,Ej;
		double mean,sigma,I0;

		gsl_histogram_get_range(g,d,&lower,&upper);
		bincenter=(lower+upper)/2.0f;

		assert(almost_same_float(bincenter,config->width*d+config->width/2.0f)==true);

		Ej=config->j*(config->j+1)-config->chempot;
		I0=(1.0f-exp(-Ej*config->maxtau))/Ej;

		fprintf(out,"%f ",bincenter);

		for(int k=0;k<NR_BLOCKSIZES;k++)
		{
			mean=block_histogram_get_mean(sst->bhs[k],d);
			sigma=sqrtf(block_histogram_get_variance(sst->bhs[k],d));

#warning Check the following six assignments!

			mean*=I0*gstats;
			mean/=g0stats;
			mean/=(upper-lower);

			sigma*=I0*gstats;
			sigma/=g0stats;
			sigma/=(upper-lower);

			if(k==0)
				fprintf(out,"%f %f ",mean,sigma);
			else
				fprintf(out,"%f ",sigma);
		}

		fprintf(out,"%f\n",exp(-Ej*bincenter));
	}

	/*
		Finally we fill the structure output
	*/

	output->avglength=((double)(avglength[0]))/((double)(avglength[1]));
	output->avgorder=((double)(avgorder[0]))/((double)(avgorder[1]));

	/*
		That's all folks. Cleaning up.
	*/

	if(out)
		fclose(out);

	gsl_histogram_free(g);
	gsl_histogram_free(g0);
	gsl_histogram_free(g1);
	gsl_histogram_free(g2);

	fini_super_sampler(sst);

	if(config->liveplot)
		gnuplot_close(gp);

	return 0;
}
