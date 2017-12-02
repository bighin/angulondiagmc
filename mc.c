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
#include "diagrams.h"
#include "updates.h"
#include "stat.h"
#include "aux.h"
#include "graphs.h"
#include "phonon.h"
#include "debug.h"
#include "config.h"

#warning UseMe!
#if 0

double calculate_propagators_and_vertices_ratio(struct diagram_t *numerator,int startmidpoint1,int endmidpoint1,
						struct diagram_t *denominator,int startmidpoint2,int endmidpoint2)
{
	assert(endmidpoint1>=startmidpoint1);
	assert(endmidpoint2>=startmidpoint2);

	assert((endmidpoint1-startmidpoint1)==(endmidpoint2-startmidpoint2));
	x
	int length=endmidpoint1-startmidpoint1;
	double ret=0.0f;

	for(c=0;c<=length;c++)
	{

#warning This ratio can be optimized (taken in a single step)

		ret*=calculate_free_propagator_weight(numerator,get_left_neighbour(numerator,startmidpoint1+c));
		ret/=calculate_free_propagator_weight(denominator,get_left_neighbour(denominator,startmidpoint2+c));

		ret*=calculate_vertex_weight(numerator,startmidpoint1+c);
		ret/=calculate_vertex_weight(denominator,startmidpoint2+c);
	}

	assert((startmidpoint1+c)==endmidpoint1);
	assert((startmidpoint2+c)==endmidpoint2);

#warning This ratio can be optimized (taken in a single step)

	ret*=calculate_free_propagator_weight(numerator,get_right_neighbour(numerator,endmidpoint1));
	ret/=calculate_free_propagator_weight(denominator,get_right_neighbour(denominator,endmidpoint2));

	return ret;
}

#endif

#include "inih/ini.h"
#include "libprogressbar/progressbar.h"
#include "gnuplot_i/gnuplot_i.h"

/*
	This function checks if -- for every propagator -- one has |m| <= j
*/

bool propagators_are_physical(struct diagram_t *dgr)
{
	int c;

	for(c=0;c<get_nr_free_propagators(dgr);c++)
	{
		struct g0_t *g0=get_free_propagator(dgr,c);

		if(g0->j<0)
			return false;
	
		if(abs(g0->m)>g0->j)
			return false;
	}

	return true;
}

/*
	This function checks if the diagrams conserves the angular momentum,
	i.e. if the first and last free propagators have the same angular momentum.
*/

bool angular_momentum_is_conserved(struct diagram_t *dgr)
{
	struct g0_t *first,*last;
	
	first=get_free_propagator(dgr,0);
	last=get_free_propagator(dgr,get_nr_free_propagators(dgr)-1);

	return (first->j==last->j)?(true):(false);
}

/*
	The weight 'penalty' for diagrams not conserving angular momentum
*/

double unphysical_penalty(struct diagram_t *dgr)
{
	struct g0_t *first,*last;
	double firstlen,lastlen,x;
	
	first=get_free_propagator(dgr,0);
	last=get_free_propagator(dgr,get_nr_free_propagators(dgr)-1);

	firstlen=first->endtau-first->starttau;
	lastlen=last->endtau-last->starttau;

	assert(firstlen>0);
	assert(lastlen>0);

	x=pow(last->j-first->j,2.0f)*(firstlen+lastlen);
	
	return exp(-x);
}

/*
	This routine calculates the weight all the propagators between two midpoints,
	including the propagators preceding the first midpoint, and the propagators following
	the last midpoint.

	The weight of all enclosed vertices is also calculated.
*/

double calculate_propagators_and_vertices(struct diagram_t *dgr,int startmidpoint,int endmidpoint)
{
	double ret=1.0f;
	int c;

	for(c=startmidpoint;c<=endmidpoint;c++)
	{
		ret*=calculate_free_propagator_weight(dgr,get_left_neighbour(dgr,c));
		ret*=calculate_vertex_weight(dgr,c);
	}

	ret*=calculate_free_propagator_weight(dgr,get_right_neighbour(dgr,endmidpoint));

	return ret;
}

/*
	(Delta_j)_i is defined as the difference in angular momentum between
	the free propagators before and after the i-th vertex, i.e.

	(Delta_j)_i = j_{i+1} - j_i

	The following functions allow to get and set (Delta_j)_i.
*/

int deltaj(struct diagram_t *dgr,int vertex)
{
	struct g0_t *left,*right;
	
	left=get_left_neighbour(dgr,vertex);
	right=get_right_neighbour(dgr,vertex);

	return right->j-left->j;
}

void change_deltaj(struct diagram_t *dgr,int index,int newdeltaj)
{
	int c,nr_free_propagators,olddeltaj;
	
	nr_free_propagators=get_nr_free_propagators(dgr);
	olddeltaj=deltaj(dgr,index);

	for(c=index+1;c<nr_free_propagators;c++)
		get_free_propagator(dgr,c)->j+=newdeltaj-olddeltaj;
}

/*
	Now it's time to start with all the different updates!
*/

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
	double tau1,tau2,weightratio,acceptance_ratio,oldweight;
	bool is_accepted;

	struct arc_t *thisline;
	struct diagram_t *old;

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

	if(tau2>=dgr->endtau)
		return UPDATE_UNPHYSICAL;

	/*
		We save the current diagram and the line is added...
	*/

	old=diagram_clone(dgr);
	oldweight=diagram_weight(dgr);

	weightratio=1.0f/calculate_propagators_and_vertices(dgr,0,get_nr_vertices(dgr)-1);

	diagram_add_phonon_line(dgr,tau1,tau2,lambda,mu);

	/*
		...and finally we fix the \Delta_j's at every vertex
	*/

	thisline=get_phonon_line(dgr,get_nr_phonons(dgr)-1);
	weightratio*=calculate_arc_weight(dgr,thisline);
	
	change_deltaj(dgr,thisline->endmidpoint,deltaj2);
	change_deltaj(dgr,thisline->startmidpoint,deltaj1);

	if(propagators_are_physical(dgr)==false)
	{	
#if 0
		diagram_copy(old,dgr);
#else
		change_deltaj(dgr,thisline->endmidpoint,0);
		change_deltaj(dgr,thisline->startmidpoint,0);
		diagram_remove_phonon_line(dgr,get_nr_phonons(dgr)-1);

                assert(almost_same_float(oldweight,diagram_weight(dgr))==true);
#endif

		fini_diagram(old);

		return UPDATE_REJECTED;
	}

	/*
		We calculate the weight ratio (the fundamental quantity in accepting/rejecting the update)
		and we use it to update the diagram weight, as well.
	*/

	weightratio*=calculate_propagators_and_vertices(dgr,0,get_nr_vertices(dgr)-1);

	if(weightratio<0.0f)
		dgr->sign*=-1;

	if((!isinf(diagram_weight(dgr)))&&(!isinf(oldweight)))
		assert(almost_same_float(weightratio,diagram_weight(dgr)/oldweight));

	/*
		Finally we calculate the acceptance ratio for the update.
	*/

	acceptance_ratio=fabs(weightratio);
	acceptance_ratio/=1.0f/(1+MAXLAMBDA);
	acceptance_ratio/=1.0f/(1+2*lambda);	
	acceptance_ratio/=1.0f/(1+lambda);
	acceptance_ratio/=1.0f/dgr->endtau;
	acceptance_ratio/=phonon_pdf(dgr->phonon_ctx,lambda,tau2-tau1);
	acceptance_ratio*=1.0f/get_nr_phonons(dgr);

	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
#if 0
		diagram_copy(old,dgr);
#else
		change_deltaj(dgr,thisline->endmidpoint,0);
		change_deltaj(dgr,thisline->startmidpoint,0);
		diagram_remove_phonon_line(dgr,get_nr_phonons(dgr)-1);

		if(weightratio<0.0f)
			dgr->sign*=-1;

                assert(almost_same_float(oldweight,diagram_weight(dgr))==true);
#endif

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

	int target,lambda,mu,nr_available_phonons,startmidpoint,endmidpoint,deltaj1,deltaj2;
	double weightratio,targetweight,acceptance_ratio,tau1,tau2,oldweight;
	bool is_accepted;

	nr_available_phonons=get_nr_phonons(dgr);

	if(nr_available_phonons<=0)
		return UPDATE_UNPHYSICAL;

	target=gsl_rng_uniform_int(dgr->rng_ctx,get_nr_phonons(dgr));
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

	old=diagram_clone(dgr);
	oldweight=diagram_weight(dgr);

	change_deltaj(dgr,endmidpoint,0);
	change_deltaj(dgr,startmidpoint,0);

	diagram_remove_phonon_line(dgr,target);

	if(propagators_are_physical(dgr)==false)
	{
#if 0
		diagram_copy(old,dgr);
#else
		diagram_add_phonon_line(dgr,tau1,tau2,lambda,mu);
		change_deltaj(dgr,endmidpoint,deltaj2);
		change_deltaj(dgr,startmidpoint,deltaj1);

                assert(almost_same_float(oldweight,diagram_weight(dgr))==true);
#endif

		fini_diagram(old);

		return UPDATE_REJECTED;
	}

	weightratio=1.0f/targetweight;

	weightratio*=calculate_propagators_and_vertices(dgr,0,get_nr_vertices(dgr)-1);
	weightratio/=calculate_propagators_and_vertices(old,0,get_nr_vertices(old)-1);

	if(weightratio<0.0f)
		dgr->sign*=-1;

	if((!isinf(diagram_weight(dgr)))&&(!isinf(diagram_weight(old))))
		assert(almost_same_float(weightratio,diagram_weight(dgr)/diagram_weight(old)));

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
#if 0
		diagram_copy(old,dgr);
#else
		diagram_add_phonon_line(dgr,tau1,tau2,lambda,mu);
		change_deltaj(dgr,endmidpoint,deltaj2);
		change_deltaj(dgr,startmidpoint,deltaj1);

		if(weightratio<0.0f)
			dgr->sign*=-1;

                assert(almost_same_float(oldweight,diagram_weight(dgr))==true);
#endif

		fini_diagram(old);

		return UPDATE_REJECTED;
	}

	fini_diagram(old);

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

	int c,d;
	FILE *out;
	char fname[1024],progressbar_text[1024];

	progressbar *progress;
	gnuplot_ctrl *gp;

	struct timeval starttime;

#define DIAGRAM_NR_UPDATES	(3)

	int (*updates[DIAGRAM_NR_UPDATES])(struct diagram_t *,struct configuration_t *);
	char *update_names[DIAGRAM_NR_UPDATES];

	int proposed[DIAGRAM_NR_UPDATES],accepted[DIAGRAM_NR_UPDATES],rejected[DIAGRAM_NR_UPDATES];
	int total_proposed,total_accepted,total_rejected;
	long long int avgorder[2];
	long long int physical_diagrams,unphysical_diagrams;
	double avglength[2];

	/*
		We set up the updates we will be using
	*/

	updates[0]=update_length;
	updates[1]=update_add_phonon_line;
	updates[2]=update_remove_phonon_line;

	update_names[0]="UpdateLength";
	update_names[1]="AddPhononLine";
	update_names[2]="RemovePhononLine";

	/*
		We reset the update statistics
	*/

	for(d=0;d<DIAGRAM_NR_UPDATES;d++)
		proposed[d]=accepted[d]=rejected[d]=0;
	
	avgorder[0]=avgorder[1]=0;
	avglength[0]=avglength[1]=0.0f;
	
	physical_diagrams=unphysical_diagrams=0;

	/*
		We print some informative message, and then we open the log file
	*/

	fprintf(stderr,"Performing %d iterations, using %d thread(s)\n",config->iterations,config->nthreads);

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

				flags=PRINT_TOPOLOGY|PRINT_INFO0;
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
#endif

#warning Statistics here, even order by order!

		if(angular_momentum_is_conserved(dgr))
		{
			double w;
			
			/*
				The histograms are updated only if we are in the physical sector,
				where the diagram conserves angular momentum.
			*/

			w=dgr->sign;
			
			gsl_histogram_accumulate(g,dgr->endtau,w);

			if(get_nr_phonons(dgr)==0)
				gsl_histogram_accumulate(g0,dgr->endtau,w);

			if(get_nr_phonons(dgr)==1)
				gsl_histogram_accumulate(g1,dgr->endtau,w);

			if(get_nr_phonons(dgr)==2)
				gsl_histogram_accumulate(g2,dgr->endtau,w);
		
			physical_diagrams++;
		}
		else
		{
			unphysical_diagrams++;
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

				ratio=((double)(physical_diagrams))/((double)(physical_diagrams+unphysical_diagrams));

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

	fprintf(out,"# <Bin center> <Average> <Stddev> <Free rotor>\n");

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

		gsl_histogram_get_range(g,d,&lower,&upper);
		bincenter=(lower+upper)/2.0f;

		assert(almost_same_float(bincenter,config->width*d+config->width/2.0f)==true);

		Ej=config->j*(config->j+1)-config->chempot;

		fprintf(out,"%f %f %f ",bincenter,gsl_histogram_get(g,d),exp(-Ej*bincenter));
		fprintf(out,"%f %f %f\n",gsl_histogram_get(g0,d),gsl_histogram_get(g1,d),gsl_histogram_get(g2,d));
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

	if(config->liveplot)
		gnuplot_close(gp);

	return 0;
}
