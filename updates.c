#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <gsl/gsl_sf.h>

#include "diagrams.h"
#include "aux.h"
#include "debug.h"
#include "updates.h"
#include "phonon.h"

/*
	This function updates all the cross-referencing inside a diagram,
	when a new midpoint and vertex are added at a certain position.
*/

void diagram_update_xrefs(struct diagram_t *dgr,int position)
{
	int c;
	
	/*
		We update all the references contained in the phonon lines,
		except for the last one which has just been added.
	*/
	
	for(c=0;c<get_nr_phonons(dgr)-1;c++)
	{
		struct arc_t *ph=get_phonon_line(dgr,c);

		if(ph->startmidpoint>=position)
			ph->startmidpoint++;

		if(ph->endmidpoint>=position)
			ph->endmidpoint++;
	}

	/*
		We update all the propagators coming after
	*/
	
	for(c=position+2;c<get_nr_free_propagators(dgr);c++)
	{
		struct g0_t *g0=get_free_propagator(dgr,c);

		g0->startmidpoint++;
		g0->endmidpoint++;
	}

	/*
		Finally we also need to update the pointers inside each vertex
	*/

	for(c=position+1;c<get_nr_vertices(dgr);c++)
	{
		struct vertex_t *vtx=get_vertex(dgr,c);

		if(c!=position+1)
			vtx->left++;

		vtx->right++;
	}
}

void diagram_add_start_midpoint(struct diagram_t *dgr,int c,double tau,struct arc_t *phononline)
{
	struct vertex_t *thisvertex,*leftvertex,*rightvertex;
	struct g0_t *leftline,*rightline;

	/*
		We insert a new midpoint and a new vertex
	*/

	vlist_add_element(dgr->midpoints,&tau,c);
	thisvertex=vlist_add_empty(dgr->vertices,c);

	leftvertex=(c>0)?(get_vertex(dgr,c-1)):(NULL);
	rightvertex=((c+1)<get_nr_vertices(dgr))?(get_vertex(dgr,c+1)):(NULL);

	/*
		We insert a new propagator, as well...
	*/

	vlist_add_empty(dgr->free_propagators,c+1);

	/*
		...and we fix all the interconnections between propagators
	*/

	leftline=get_free_propagator(dgr,c);
	rightline=get_free_propagator(dgr,c+1);

	leftline->startmidpoint=c-1;
	leftline->endmidpoint=c;

	rightline->startmidpoint=c;
	rightline->endmidpoint=c+1;

	leftline->starttau=get_midpoint(dgr,c-1);
	leftline->endtau=get_midpoint(dgr,c);

	rightline->starttau=get_midpoint(dgr,c);
	rightline->endtau=get_midpoint(dgr,c+1);

	/*
		We now setup the newly-created vertex and the connections between vertices.
	*/

	thisvertex->left=leftline;
	thisvertex->right=rightline;
	thisvertex->phononline=phononline;

	if(leftvertex)
		leftvertex->right=leftline;

	if(rightvertex)
		rightvertex->left=rightline;

	/*
		Finally the angular momentum from the "oldest" line is
		given to the newly-created one, as well.
	*/

	rightline->j=leftline->j;
	rightline->m=leftline->m;
	rightline->arcs_over_me=leftline->arcs_over_me;
}

void diagram_add_end_midpoint(struct diagram_t *dgr,int c,double tau,struct arc_t *phononline)
{
	struct vertex_t *thisvertex,*leftvertex,*rightvertex;
	struct g0_t *leftline,*rightline;

	/*
		We insert a new midpoint and a new vertex
	*/

	vlist_add_element(dgr->midpoints,&tau,c);
	thisvertex=vlist_add_empty(dgr->vertices,c);
	
	leftvertex=(c>0)?(get_vertex(dgr,c-1)):(NULL);
	rightvertex=((c+1)<get_nr_vertices(dgr))?(get_vertex(dgr,c+1)):(NULL);

	/*
		We insert a new propagator, as well...
	*/

	vlist_add_empty(dgr->free_propagators,c);

	/*
		...and we fix all the interconnections between propagators
	*/

	leftline=get_free_propagator(dgr,c);
	rightline=get_free_propagator(dgr,c+1);

	leftline->startmidpoint=c-1;
	leftline->endmidpoint=c;

	rightline->startmidpoint=c;
	rightline->endmidpoint=c+1;

	leftline->starttau=get_midpoint(dgr,c-1);
	leftline->endtau=get_midpoint(dgr,c);

	rightline->starttau=get_midpoint(dgr,c);
	rightline->endtau=get_midpoint(dgr,c+1);

	/*
		We now setup the newly-created vertex and the connections between vertices.
	*/

	thisvertex->left=leftline;
	thisvertex->right=rightline;
	thisvertex->phononline=phononline;

	if(leftvertex)
		leftvertex->right=leftline;

	if(rightvertex)
		rightvertex->left=rightline;

	/*
		Finally the angular momentum from the "oldest" line is
		given to the newly-created one, as well.
	*/

	leftline->j=rightline->j;
	leftline->m=rightline->m;
	leftline->arcs_over_me=rightline->arcs_over_me;
}

double calculate_arc_weight(struct diagram_t *dgr,struct arc_t *arc)
{
	double timediff;
	double ret=1.0f;

	timediff=arc->endtau-arc->starttau;

	assert(timediff>=0);

	ret*=chi_lambda(dgr->phonon_ctx,arc->lambda,timediff);
	ret*=pow(-1.0f,arc->lambda);

	return ret;
}

void diagram_add_phonon_line(struct diagram_t *dgr,double tau1,double tau2,double k,int lambda,int mu)
{
	struct arc_t *arc;
	int c,lo,hi;

	/*
		We start by appending a new phonon arc, and setting it up.
	*/

	arc=vlist_append_empty(dgr->phonons);

	arc->k=k;
	arc->lambda=lambda;
	arc->mu=mu;

	arc->starttau=tau1;
	arc->endtau=tau2;

	assert(tau1<dgr->endtau);
	assert(tau1>0.0f);

	/*
		Next we locate the insertion point...
	*/

	lo=get_nr_midpoints(dgr);
	for(c=0;c<get_nr_midpoints(dgr);c++)
	{		
		if(tau1<get_midpoint(dgr,c))
		{
			lo=c;
			break;
		}
	}

	/*
		...and we attach a new midpoint and vertex there!
	*/

	diagram_add_start_midpoint(dgr,lo,tau1,arc);
	diagram_update_xrefs(dgr,lo);
	arc->startmidpoint=lo;

	assert(tau2>tau1);
	assert(tau2<dgr->endtau);
	assert(tau2>0.0f);

	/*
		The same procedure is repeated to attach the end of the
		phonon arc to the main diagram.
	*/

	hi=get_nr_midpoints(dgr);
	for(c=0;c<get_nr_midpoints(dgr);c++)
	{
		if(tau2<get_midpoint(dgr,c))
		{
			hi=c;
			break;
		}
	}

	diagram_add_end_midpoint(dgr,hi,tau2,arc);
	diagram_update_xrefs(dgr,hi);
	arc->endmidpoint=hi;

	for(c=arc->startmidpoint;c<arc->endmidpoint;c++)
		get_right_neighbour(dgr,c)->arcs_over_me++;

	assert(hi>lo);
}

void diagram_remove_start_midpoint(struct diagram_t *dgr,int c)
{
	struct g0_t *left,*right;
	struct vertex_t *leftvertex,*rightvertex,*thisvertex;

	assert(c>=0);
	assert(c<get_nr_vertices(dgr));

	/*
		The right propagator is removed, and connections between propagators are fixed...
	*/

	delete_right_neighbour(dgr,c);

	vlist_remove_element(dgr->midpoints,c);
	vlist_remove_element(dgr->vertices,c);

	left=get_left_neighbour(dgr,c);
	right=get_right_neighbour(dgr,c);

	left->endtau=get_midpoint(dgr,c);
	right->starttau=get_midpoint(dgr,c);

	/*
		...and also the connection between vertices are fixed.
	*/

	assert(get_nr_vertices(dgr)>0);

	thisvertex=get_vertex(dgr,c);

	leftvertex=(c>0)?(get_vertex(dgr,c-1)):(NULL);
	rightvertex=((c+1)<get_nr_vertices(dgr))?(get_vertex(dgr,c+1)):(NULL);

	thisvertex->left=left;
	thisvertex->right=right;

	if(leftvertex)
		leftvertex->right=left;

	if(rightvertex)
		rightvertex->left=right;
}

void diagram_remove_end_midpoint(struct diagram_t *dgr,int c)
{
	struct g0_t *left,*right;
	struct vertex_t *leftvertex,*rightvertex,*thisvertex;

	assert(c>=0);
	assert(c<get_nr_midpoints(dgr));
	assert((c+1)<get_nr_free_propagators(dgr));

	/*
		The left propagator is removed, and connections between propagators are fixed.
	*/

	delete_left_neighbour(dgr,c);

	vlist_remove_element(dgr->midpoints,c);
	vlist_remove_element(dgr->vertices,c);
	
	/*
		The value of c is decreased because a propagator
		has been deleted to the left...
	*/

	c--;

	/*
		...but then we have to check everywhere for negative values of c!
	*/

	if(c>=0)
	{
		left=get_left_neighbour(dgr,c);
		left->endtau=get_midpoint(dgr,c);
	}


	right=get_right_neighbour(dgr,c);
	right->starttau=get_midpoint(dgr,c);

	/*
		Finally the connections between vertices are fixed.
	*/

	thisvertex=(c>=0)?(get_vertex(dgr,c)):(NULL);

	leftvertex=(c>0)?(get_vertex(dgr,c-1)):(NULL);
	rightvertex=((c+1)<get_nr_vertices(dgr))?(get_vertex(dgr,c+1)):(NULL);

	if(thisvertex)
	{
		thisvertex->left=left;
		thisvertex->right=right;
	}

	if(leftvertex)
		leftvertex->right=left;

	if(rightvertex)
		rightvertex->left=right;
}

bool diagram_remove_phonon_line(struct diagram_t *dgr,int position)
{
	int c,startmidpoint,endmidpoint;
	struct arc_t *arc;
	bool result=true;

	assert(position<get_nr_phonons(dgr));

	/*
		At first we find the arc we are about to remove and fetch some basic data
	*/

	arc=vlist_get_element(dgr->phonons,position);

	startmidpoint=arc->startmidpoint;
	endmidpoint=arc->endmidpoint;

	for(c=arc->startmidpoint;c<arc->endmidpoint;c++)
		get_right_neighbour(dgr,c)->arcs_over_me--;

	/*
		Now the hardcore part starts: we update all the field in a diagram_t struct.
	*/

	assert(startmidpoint<endmidpoint);

	for(c=0;c<get_nr_vertices(dgr);c++)
	{
		struct vertex_t *vtx=get_vertex(dgr,c);
		
		if(vtx->phononline>arc)
			vtx->phononline--;
	}

	vlist_remove_element(dgr->phonons,position);

	/*
		Before removing a midpoint we have to update all the references:
		the index of every midpoint after the removed one will be decreased
		by one, and so we have to update all the structures pointing to a midpoint.
	*/

	for(c=0;c<get_nr_phonons(dgr);c++)
	{
		struct arc_t *ph=get_phonon_line(dgr,c);
	
		if(ph->startmidpoint>=startmidpoint)
			ph->startmidpoint--;

		if(ph->endmidpoint>=startmidpoint)
			ph->endmidpoint--;
	}

	for(c=startmidpoint+1;c<get_nr_free_propagators(dgr);c++)
	{
		struct g0_t *g0=get_free_propagator(dgr,c);
		
		g0->startmidpoint--;
		g0->endmidpoint--;
	}

	for(c=startmidpoint+1;c<get_nr_vertices(dgr);c++)
	{
		struct vertex_t *vtx=get_vertex(dgr,c);

		if(c>startmidpoint+1)
			vtx->left--;
		
		vtx->right--;
	}

	diagram_remove_start_midpoint(dgr,startmidpoint);

	endmidpoint--;

	for(c=0;c<get_nr_phonons(dgr);c++)
	{
		struct arc_t *ph=get_phonon_line(dgr,c);
	
		if(ph->startmidpoint>=endmidpoint)
			ph->startmidpoint--;

		if(ph->endmidpoint>=endmidpoint)
			ph->endmidpoint--;
	}

	for(c=endmidpoint+1;c<get_nr_free_propagators(dgr);c++)
	{
		struct g0_t *g0=get_free_propagator(dgr,c);
		
		g0->startmidpoint--;
		g0->endmidpoint--;
	}

	for(c=endmidpoint+1;c<get_nr_vertices(dgr);c++)
	{
		struct vertex_t *vtx=get_vertex(dgr,c);

		if(c>endmidpoint+1)
			vtx->left--;
		
		vtx->right--;
	}

	diagram_remove_end_midpoint(dgr,endmidpoint);

	return result;
}

bool is_last_propagator(struct diagram_t *dgr,struct g0_t *g0)
{
	bool result;
	int nr_propagators;
	
	nr_propagators=get_nr_free_propagators(dgr);
	result=((g0==get_free_propagator(dgr,nr_propagators-1))?(true):(false));

	assert(result==((g0->endtau==dgr->endtau)?(true):(false)));

	return result;
}

double calculate_free_propagator_weight(struct diagram_t *dgr,struct g0_t *g0)
{
	double en,tau,phase;
	int j;
	
	j=g0->j;
	en=j*(j+1)-dgr->chempot;
	tau=g0->endtau-g0->starttau;

	phase=1.0f;

	if(is_last_propagator(dgr,g0)==false)
		phase=pow(-1.0f,j);

	return exp(-en*tau)*phase;
}

double calculate_vertex_weight(struct diagram_t *dgr,int index)
{
	int j1,j2,j3;
	double coupling;

	struct vertex_t *thisvertex=get_vertex(dgr,index);

	j1=thisvertex->left->j;
	j2=thisvertex->right->j;
	j3=thisvertex->phononline->lambda;

	coupling=gsl_sf_coupling_3j(2*j1,2*j2,2*j3,0,0,0)*
	         sqrtf((2.0f*j1+1)*(2.0f*j2+1)*(2.0f*j3+1)/(4.0f*M_PI));

	return coupling;
}

void diagram_update_length(struct diagram_t *dgr,double newendtau)
{	
	int nr_free_propagators=get_nr_free_propagators(dgr);

	/*
		The following variable is used only in assert() statements,
		so we need it only if the NDEBUG macro is *not* defined.
	*/

#ifndef NDEBUG
	int nr_midpoints=get_nr_midpoints(dgr);
#endif

	assert(dgr->endtau==get_midpoint(dgr,nr_midpoints));
	assert(newendtau>get_midpoint(dgr,nr_midpoints-1));

	dgr->endtau=newendtau;
	get_free_propagator(dgr,nr_free_propagators-1)->endtau=newendtau;
}

void save_free_propagators(struct diagram_t *dgr,struct free_propagators_ctx_t *fpc,int lo,int hi)
{
	int c;
	
	assert(fpc);
	assert(hi>=lo);

#ifdef DEBUG7
	printf("I was asked to save the propagators between v1=%d and v2=%d\n",lo,hi);
#endif

	fpc->nr_free_propagators=hi-lo+1;
	fpc->lo=lo;
	fpc->hi=hi;

	fpc->js=malloc(sizeof(int)*fpc->nr_free_propagators);

	for(c=lo;c<hi;c++)
	{
		fpc->js[c-lo]=get_right_neighbour(dgr,c)->j;

#ifdef DEBUG7
		printf("Saving %d in position %d\n",fpc->js[c-lo],c);
#endif
	}

#ifdef DEBUG7	
	printf("Done!\n\n");
#endif
}

void release_free_propagators_ctx(struct free_propagators_ctx_t *fpc)
{
	if(fpc)
	{
		if(fpc->js)
			free(fpc->js);
	}
}

void restore_free_propagators(struct diagram_t *dgr,struct free_propagators_ctx_t *fpc)
{
	int c;

	assert(fpc);
	assert(fpc->js);

	for(c=fpc->lo;c<=fpc->hi;c++)
	{
		dgr->weight/=calculate_vertex_weight(dgr,c);
	
		if(c!=fpc->hi)
			dgr->weight/=calculate_free_propagator_weight(dgr,get_right_neighbour(dgr,c));
	}

	for(c=fpc->lo;c<fpc->hi;c++)
	{
		get_right_neighbour(dgr,c)->j=fpc->js[c-fpc->lo];

#ifdef DEBUG7
		printf("Restoring %d in position %d\n",fpc->js[c-fpc->lo],c);
#endif
	}

	for(c=fpc->lo;c<=fpc->hi;c++)
	{
		dgr->weight*=calculate_vertex_weight(dgr,c);
	
		if(c!=fpc->hi)
			dgr->weight*=calculate_free_propagator_weight(dgr,get_right_neighbour(dgr,c));
	}

#ifdef DEBUG7
	printf("Done!\n\n");
#endif

	release_free_propagators_ctx(fpc);
}

/*
	Given a diagram and a starting and ending vertex (lo and hi),
	this function recouples momenta on free propagators between lo and hi included.

	Great care must be taken in extracting *fairly* one possible coupling
	among all of them.
*/

bool recouple(struct diagram_t *dgr,int lo,int hi)
{
	int c,d,length,targettotal;
	int *deltalist,*jlist;
	
	assert(hi>=lo);
	
	length=hi-lo+1;
	deltalist=malloc(sizeof(int)*length);
	jlist=malloc(sizeof(int)*length);

#define MAX_TRIES	(1024)

	targettotal=get_right_neighbour(dgr,hi)->j-get_left_neighbour(dgr,lo)->j;

	for(d=0;d<=MAX_TRIES;d++)
	{
		int total=0;
		bool failedtc;
		
		/*
			This is just to silence CLang's static analyzer. Actually deltalist[0]
			is correctly inizitialed in the subsequent 'for' loop.
		
			Probably... I should check this better!
		*/

		deltalist[0]=0;

		/*
			At each vertex one can define a delta value, defined as the difference
			in angular momenta between the left propagator and the right propagator.

			Clearly, instead of choosing one configuration between all angular
			momentum values, one can equivalently choose a delta configuration.

			Which values can delta take at each vertex? Given lambda, the phonon
			propagator's angular momentum, the triangular condition at each vertex
			dictates that the following delta values are allowed:
		
			lambda = 0 --> delta = 0
			lambda = 1 --> delta = -1, 0, 1
			lambda = 2 --> delta = -2, -1, 0, 1, 2
		
			However, it turns out that a more restrictive condition holds and that
			some of this integer values are not allowed: at every vertex,
			according to the diagrammatic rules, one has a 3j symbol with all m's
			set to zero. This vertex imposes that the sums of angular momenta at a
			vertex should be an even integer.
		
			Some of the delta values above would change the parity, resuling in an
			odd total sum of angular momenta at a vertex.

			Thus, the correct list is:

			lambda = 0 --> delta = 0
			lambda = 1 --> delta = -1, 1
			lambda = 2 --> delta = -2, 0, 2
			and so on...

			and we select one random value among these deltas, for each vertex.
		*/

		for(c=lo;c<=hi;c++)
		{
			struct vertex_t *vtx;
			int lambda;
			
			vtx=get_vertex(dgr,c);
			lambda=vtx->phononline->lambda;

			deltalist[c-lo]=2*gsl_rng_uniform_int(dgr->rng_ctx,1+lambda)-lambda;
		}

		for(c=lo;c<=hi;c++)
		{
			total+=deltalist[c-lo];
		}

		/*
			jlist[i] referes to the propagator *after* the i-th vertex
		*/

		jlist[0]=get_left_neighbour(dgr,lo)->j+deltalist[0];
		
		failedtc=false;
		for(c=lo+1;c<=hi;c++)
		{
			int lambda=get_vertex(dgr,c)->phononline->lambda;
			
			jlist[c-lo]=jlist[c-lo-1]+deltalist[c-lo];
	
			if(check_triangle_condition_and_parity_from_js(jlist[c-lo],jlist[c-lo-1],lambda)==false)
				failedtc=true;
		}

		if((total==targettotal)&&(failedtc==false))
			break;
	}

	/*
		We have tried to get a deltalist summing to 0 for MAX_TRIES times.

		If we didn't managed to do so, it is time to throw an error...
	*/

	if(d==MAX_TRIES)
	{
		printf("Failed to find an acceptable recoupling after %d tries, shamefully exiting. Please debug this.\n",MAX_TRIES);
		exit(0);
	}

	/*
		We have a working deltalist, let's use it to rewrite the free propagators.
	*/

	for(c=lo;c<hi;c++)
	{
		struct vertex_t *vtx=get_vertex(dgr,c);

		vtx->right->j=vtx->left->j+deltalist[c-lo];
	}

	/*
		As a final check, we verify the consistency of the last coupling.
	*/

#ifndef NDEBUG
	{
		struct vertex_t *vtx=get_vertex(dgr,hi);
		assert(vtx->right->j==vtx->left->j+deltalist[hi-lo]);
	}
#endif

	if(jlist)
		free(jlist);

	if(deltalist)
		free(deltalist);
	
	return true;
}
