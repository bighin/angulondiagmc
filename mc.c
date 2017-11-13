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

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "mc.h"
#include "diagrams.h"
#include "updates.h"
#include "stat.h"
#include "aux.h"
#include "graphs.h"
#include "phonon.h"
#include "debug.h"
#include "config.h"
#include "sectors.h"

#include "inih/ini.h"
#include "libprogressbar/progressbar.h"
#include "gnuplot_i/gnuplot_i.h"

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

int get_nr_available_phonons(struct diagram_t *dgr)
{
	int c,available;

	for(c=available=0;c<get_nr_phonons(dgr);c++)
	{
		int startmidpoint,endmidpoint;
		
		startmidpoint=get_phonon_line(dgr,c)->startmidpoint;
		endmidpoint=get_phonon_line(dgr,c)->endmidpoint;
		
		if((deltaj(dgr,startmidpoint)==0)&&(deltaj(dgr,endmidpoint)==0))
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
		
		if((deltaj(dgr,startmidpoint)==0)&&(deltaj(dgr,endmidpoint)==0))
		{
			if(available==target)
				return c;
			
			available++;
		}
	}

	assert(false);
	
	return -1;
}

struct g0_t *find_propagator_containing_time(struct diagram_t *dgr,double tau)
{
	int c;
	
	assert(tau>=dgr->mintau);
	assert(tau<=dgr->endtau);

	for(c=0;c<get_nr_free_propagators(dgr);c++)
	{
		struct g0_t *thisg0=get_free_propagator(dgr,c);
		
		if((tau>=thisg0->starttau)&&(tau<=thisg0->endtau))
			return thisg0;
	}

	assert(false);

	return NULL;
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
	int lambda,mu;
	double tau1,tau2,weightratio,acceptance_ratio;
	bool is_accepted;

	struct arc_t *thisline;
	struct diagram_t *old;

	/*
		Do we want to limit the simulation only to first order diagrams (or at a certain maximum order)?
	*/

	if(get_nr_phonons(dgr)>=cfg->maxorder)
		return UPDATE_UNPHYSICAL;

	/*
		The new lambda is chosen using a uniform distribution.
	*/

	lambda=gsl_rng_uniform_int(dgr->rng_ctx,1+MAXLAMBDA);
	mu=gsl_rng_uniform_int(dgr->rng_ctx,1+2*lambda)-lambda;

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

	diagram_add_phonon_line(dgr,tau1,tau2,lambda,mu);

	/*
		...and then we check if the new line makes sense physically!
	*/

	thisline=get_phonon_line(dgr,get_nr_phonons(dgr)-1);

	/*
		We calculate the weight ratio (the fundamental quantity in accepting/rejecting the update)
		and we use it to update the diagram weight, as well.
	*/

	weightratio=calculate_arc_weight(dgr,thisline);
	weightratio*=calculate_vertex_weight(dgr,thisline->startmidpoint);
	weightratio*=calculate_vertex_weight(dgr,thisline->endmidpoint);

	if((!isinf(diagram_weight(dgr)))&&(!isinf(diagram_weight(old))))
		assert(almost_same_float(weightratio,diagram_weight(dgr)/diagram_weight(old)));

	assert(diagram_weight(dgr)>=0);
	assert(diagram_weight(old)>=0);

	/*
		Finally we calculate the acceptance ratio for the update.
	*/

	acceptance_ratio=weightratio;
	acceptance_ratio/=1.0f/(1+MAXLAMBDA);
	acceptance_ratio/=1.0f/(1+2*lambda);
	acceptance_ratio/=1.0f/dgr->endtau;
	acceptance_ratio/=phonon_pdf(dgr->phonon_ctx,lambda,tau2-tau1);
	acceptance_ratio*=1.0f/get_nr_available_phonons(dgr);

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

	int target,lambda,nr_available_phonons,startmidpoint,endmidpoint;
	double weightratio,targetweight,acceptance_ratio,tau1,tau2;
	bool is_accepted;

	nr_available_phonons=get_nr_available_phonons(dgr);

	if(nr_available_phonons<=0)
		return UPDATE_UNPHYSICAL;

	target=get_random_available_phonon(dgr);
	arc=get_phonon_line(dgr,target);
	targetweight=calculate_arc_weight(dgr,arc);

	tau1=arc->starttau;
	tau2=arc->endtau;
	lambda=arc->lambda;
	startmidpoint=arc->startmidpoint;
	endmidpoint=arc->endmidpoint;

	old=diagram_clone(dgr);
	
	diagram_remove_phonon_line(dgr,target);

	weightratio=1.0f/targetweight;
	weightratio/=calculate_vertex_weight(old,startmidpoint);
	weightratio/=calculate_vertex_weight(old,endmidpoint);

	if((!isinf(diagram_weight(dgr)))&&(!isinf(diagram_weight(old))))
		assert(almost_same_float(weightratio,diagram_weight(dgr)/diagram_weight(old)));

	assert(diagram_weight(dgr)>=0);
	assert(diagram_weight(old)>=0);

	/*
		Finally we calculate the acceptance ratio for the update.
	*/

	acceptance_ratio=weightratio;
	acceptance_ratio*=1.0f/(1+MAXLAMBDA);
	acceptance_ratio*=1.0f/(1+2*lambda);
	acceptance_ratio*=1.0f/dgr->endtau;
	acceptance_ratio*=phonon_pdf(dgr->phonon_ctx,lambda,tau2-tau1);
	acceptance_ratio/=1.0f/nr_available_phonons;

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

int update_recouple(struct diagram_t *dgr,struct configuration_t *cfg)
{
	int nr_vertices;
	double weightratio,acceptance_ratio;
	bool is_accepted;

#ifndef NDEBUG
	bool result;
#endif

	struct diagram_t *old;

	nr_vertices=get_nr_vertices(dgr);
	
	if(nr_vertices<1)
		return UPDATE_UNPHYSICAL;

	old=diagram_clone(dgr);

#ifndef NDEBUG
	result=recouple(dgr,0,nr_vertices-1);
	assert(result==true);
#else
	recouple(dgr,0,nr_vertices-1);
#endif

	weightratio=calculate_propagators_and_vertices(dgr,0,nr_vertices-1);
	weightratio/=calculate_propagators_and_vertices(old,0,nr_vertices-1);

	acceptance_ratio=weightratio;

	if((!isinf(diagram_weight(dgr)))&&(!isinf(diagram_weight(old))))
		assert(almost_same_float(weightratio,diagram_weight(dgr)/diagram_weight(old)));

	assert(diagram_weight(dgr)>=0);
	assert(diagram_weight(old)>=0);

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

#define DIAGRAM_NR_UPDATES	(4)

	int (*updates[DIAGRAM_NR_UPDATES])(struct diagram_t *,struct configuration_t *);
	char *update_names[DIAGRAM_NR_UPDATES];

	int proposed[DIAGRAM_NR_UPDATES],accepted[DIAGRAM_NR_UPDATES],rejected[DIAGRAM_NR_UPDATES];
	int total_proposed,total_accepted,total_rejected;
	int samples_in_P_sector,samples_not_in_P_sector;
	int avgorder[2];
	double avglength[2];

	/*
		We set up the updates we will be using
	*/

	updates[0]=update_length;
	updates[1]=update_add_phonon_line;
	updates[2]=update_remove_phonon_line;
	updates[3]=update_recouple;

	update_names[0]="UpdateLength";
	update_names[1]="AddPhononLine";
	update_names[2]="RemovePhononLine";
	update_names[3]="RecoupleFull";

	/*
		We reset the update statistics
	*/

	for(d=0;d<DIAGRAM_NR_UPDATES;d++)
		proposed[d]=accepted[d]=rejected[d]=0;
	
	avgorder[0]=avgorder[1]=0;
	avglength[0]=avglength[1]=0.0f;

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

		assert(diagram_weight(dgr)>=0.0f);
		assert(is_in_E(dgr)==true);
		
		/*
			The histograms are updated only if we are in the physical sector
		*/
		
		if(is_in_P(dgr)==true)
		{
			gsl_histogram_increment(g,dgr->endtau);

			if(get_nr_phonons(dgr)==0)
				gsl_histogram_increment(g0,dgr->endtau);

			if(get_nr_phonons(dgr)==1)
				gsl_histogram_increment(g1,dgr->endtau);

			if(get_nr_phonons(dgr)==2)
				gsl_histogram_increment(g2,dgr->endtau);
		
			samples_in_P_sector++;
		}
		else
		{
			samples_not_in_P_sector++;
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
				double avg1,avg2;

				avg1=((double)(avgorder[0]))/((double)(avgorder[1]));
				avg2=((double)(avglength[0]))/((double)(avglength[1]));

				snprintf(progressbar_text,1024,"Progress (avg. order: %f, avg. length: %f)",avg1,avg2);

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
	fprintf(out,"# Samples taken in physical sector: %d\n",samples_in_P_sector);
	fprintf(out,"# Total samples: %d\n",samples_in_P_sector+samples_not_in_P_sector);
	fprintf(out,"# Percentage in physical sector: %f\n",((double)(samples_in_P_sector))/((double)(samples_in_P_sector+samples_not_in_P_sector)));

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

/*
	Things to add here:

	- CTRL-c handler
	- gnuplot interface
	- Number of updates
*/

#if 0
int do_diagmc_parallel(char *configfile)
{
	struct diagram_cfg_t dcfg;
	struct configuration_t config;
	struct histogram_t *ht;

	int d;
	FILE *out;
	char fname[1024];
	progressbar *progress;

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

			assert(fabs(diagram_weight(dgr)-diagram_weight(dgr))<10e-7*diagram_weight(dgr));
			diagram_check_consistency(dgr);

			totalweight=diagram_weight(dgr)*diagram_m_weight(dgr,config.use_hashtable);
			assert(totalweight>=0.0f);

			if(totalweight>10e4)
			{
				debug_weight(dgr);
				exit(0);
			}

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
	fprintf(out,"# Extrapolated quasiparticle weight: %f\n",calculate_qpw(&config,ht));
	fprintf(out,"# Extrapolated energy: %f\n",calculate_energy(&config,ht));
	fprintf(out,"#\n");
	fprintf(out,"# Sampled quantity: Green's function (G)\n");
	fprintf(out,"# Iterations: %d\n",config.iterations);
	fprintf(out,"# Nr of bins: %d\n",config.bins);
	fprintf(out,"# Bin width: %f\n",config.width);
	fprintf(out,"# (last bin is overflow)\n");
	fprintf(out,"#\n");

	fprintf(out,"# <Bin center> <Average> <Stddev> <Free rotor>\n");

	for(d=0;d<=config.bins;d++)
	{
		double tau,Ej;
		
		tau=config.width*d+config.width/2.0f;
		Ej=config.j*(config.j+1);
		
		fprintf(out,"%f %f ",tau,histogram_get_bin_average(ht,d));
		
		if(ht->sctx->smpls[d]->next>1)
			fprintf(out,"%f ",sqrtf(histogram_get_bin_variance(ht,d)/(ht->sctx->smpls[d]->next)));
		else
			fprintf(out,"%f ",0.0f);

		fprintf(out,"%f\n",exp(-(Ej-config.chempot)*tau));
	}

	/*
		That's all folks. Cleaning up.
	*/

	if(out)
		fclose(out);

	fini_histogram(ht);

	return 0;
}
#endif
