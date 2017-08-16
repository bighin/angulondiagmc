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
		Also the worms need to be modified
	*/

	for(c=0;c<get_nr_worms(dgr);c++)
	{
		struct worm_t *worm=get_worm(dgr,c);

		if(worm->startmidpoint>=position)
			worm->startmidpoint++;

		if(worm->endmidpoint>=position)
			worm->endmidpoint++;
	}

	/*
		Finally we also need to update the pointers inside each vertex
	*/

	for(c=position+1;c<get_nr_vertices(dgr);c++)
	{
		struct vertex_info_t *vtx=get_vertex(dgr,c);

		if(c!=position+1)
			vtx->left++;

		vtx->right++;
	}
}

void diagram_add_start_midpoint(struct diagram_t *dgr,int c,double tau,struct arc_t *phononline)
{
	struct vertex_info_t *thisvertex,*leftvertex,*rightvertex;
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
	thisvertex->refs=0;

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
}

void diagram_add_end_midpoint(struct diagram_t *dgr,int c,double tau,struct arc_t *phononline)
{
	struct vertex_info_t *thisvertex,*leftvertex,*rightvertex;
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
	thisvertex->refs=0;

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
}

double calculate_arc_weight(struct diagram_t *dgr,struct arc_t *arc)
{
	int j1,j2,j3;
	struct vertex_info_t *thisvertex;
	double timediff;
	double c0,c1,c2,omega0,omega1,omega2;
	double ret=1.0f;
	
	thisvertex=get_vertex(dgr,arc->startmidpoint);

	j1=thisvertex->left->j;
	j2=thisvertex->right->j;
	j3=thisvertex->phononline->lambda;

	ret*=gsl_sf_coupling_3j(2*j1,2*j2,2*j3,0,0,0)*
	     sqrtf((2.0f*j1+1)*(2.0f*j2+1)*(2.0f*j3+1)/(4.0f*M_PI));

	thisvertex=get_vertex(dgr,arc->endmidpoint);

	j1=thisvertex->left->j;
	j2=thisvertex->right->j;
	j3=thisvertex->phononline->lambda;

	ret*=gsl_sf_coupling_3j(2*j1,2*j2,2*j3,0,0,0)*
	     sqrtf((2.0f*j1+1)*(2.0f*j2+1)*(2.0f*j3+1)/(4.0f*M_PI));

	timediff=arc->endtau-arc->starttau;

	c0=dgr->c0;
	c1=dgr->c1;
	c2=dgr->c2;
	omega0=dgr->omega0;
	omega1=dgr->omega1;
	omega2=dgr->omega2;

	assert(timediff>=0);

	switch(arc->lambda)
	{
		case 0:
		ret*=c0*exp(-timediff*omega0);
		break;

		case 1:
		ret*=c1*exp(-timediff*omega1);
		break;

		case 2:
		ret*=c2*exp(-timediff*omega2);
		break;
		
		default:
		assert(false);
	}
	
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

	dgr->weight*=calculate_arc_weight(dgr,arc);

	assert(hi>lo);
}

void diagram_remove_start_midpoint(struct diagram_t *dgr,int c)
{
	struct g0_t *left,*right;
	struct vertex_info_t *leftvertex,*rightvertex,*thisvertex;

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
	struct vertex_info_t *leftvertex,*rightvertex,*thisvertex;

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
		...but then we have to check everywhere for negative zero!
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

void diagram_remove_phonon_line(struct diagram_t *dgr,int position)
{
	int c,startmidpoint,endmidpoint;
	struct arc_t *arc;

	assert(position<get_nr_phonons(dgr));

	arc=vlist_get_element(dgr->phonons,position);

	dgr->weight/=calculate_arc_weight(dgr,arc);

	startmidpoint=arc->startmidpoint;
	endmidpoint=arc->endmidpoint;

	assert(startmidpoint<endmidpoint);

	for(c=0;c<get_nr_vertices(dgr);c++)
	{
		struct vertex_info_t *vtx=get_vertex(dgr,c);
		
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

	for(c=0;c<get_nr_worms(dgr);c++)
	{
		struct worm_t *worm=get_worm(dgr,c);
	
		if(worm->startmidpoint>=startmidpoint)
			worm->startmidpoint--;

		if(worm->endmidpoint>=startmidpoint)
			worm->endmidpoint--;
	}

	for(c=startmidpoint+1;c<get_nr_free_propagators(dgr);c++)
	{
		struct g0_t *g0=get_free_propagator(dgr,c);
		
		g0->startmidpoint--;
		g0->endmidpoint--;
	}

	for(c=startmidpoint+1;c<get_nr_vertices(dgr);c++)
	{
		struct vertex_info_t *vtx=get_vertex(dgr,c);

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

	for(c=0;c<get_nr_worms(dgr);c++)
	{
		struct worm_t *worm=get_worm(dgr,c);
	
		if(worm->startmidpoint>=endmidpoint)
			worm->startmidpoint--;

		if(worm->endmidpoint>=endmidpoint)
			worm->endmidpoint--;
	}

	for(c=endmidpoint+1;c<get_nr_free_propagators(dgr);c++)
	{
		struct g0_t *g0=get_free_propagator(dgr,c);
		
		g0->startmidpoint--;
		g0->endmidpoint--;
	}

	for(c=endmidpoint+1;c<get_nr_vertices(dgr);c++)
	{
		struct vertex_info_t *vtx=get_vertex(dgr,c);

		if(c>endmidpoint+1)
			vtx->left--;
		
		vtx->right--;
	}

	diagram_remove_end_midpoint(dgr,endmidpoint);
}

double calculate_free_propagator_weight(struct diagram_t *dgr,struct g0_t *g0)
{
	double en,tau;
	int j;
	
	j=g0->j;
	en=j*(j+1)-dgr->chempot;
	tau=g0->endtau-g0->starttau;

	return exp(-en*tau);
}

void diagram_update_length(struct diagram_t *dgr,double newendtau)
{	
	int nr_midpoints=get_nr_midpoints(dgr);
	int nr_free_propagators=get_nr_free_propagators(dgr);

	assert(dgr->endtau==get_midpoint(dgr,nr_midpoints));
	assert(newendtau>get_midpoint(dgr,nr_midpoints-1));

	dgr->weight/=calculate_free_propagator_weight(dgr,get_free_propagator(dgr,nr_free_propagators-1));

	dgr->endtau=newendtau;
	get_free_propagator(dgr,nr_free_propagators-1)->endtau=newendtau;

	dgr->weight*=calculate_free_propagator_weight(dgr,get_free_propagator(dgr,nr_free_propagators-1));
}

bool diagram_add_worm(struct diagram_t *dgr,int target1,int target2,int deltalambda)
{	
	int c,nr_worms;
	struct worm_t *worm;

	assert(target1<get_nr_midpoints(dgr));
	assert(target2<get_nr_midpoints(dgr));
	assert(target1<target2);

	nr_worms=get_nr_worms(dgr);
	for(c=0;c<nr_worms;c++)
	{
		struct worm_t *local=get_worm(dgr,c);

		if((local->startmidpoint==target1)&&(local->endmidpoint==target2))
			return false;
	}

	for(c=target1;c<target2;c++)
	{
		struct vertex_info_t *thisvertex=get_vertex(dgr,c);

		if((thisvertex->right->j+deltalambda)<0)
			return false;
	}

	for(c=target1;c<target2;c++)
	{
		struct vertex_info_t *thisvertex=get_vertex(dgr,c);

		dgr->weight/=calculate_free_propagator_weight(dgr,thisvertex->right);

		thisvertex->right->j+=deltalambda;

		dgr->weight*=calculate_free_propagator_weight(dgr,thisvertex->right);
	}

#warning FIXME: the contribution to vertices needs to be calculated as well!

	/*
		At last we keep track of the worm by adding it on the global list
	*/

	worm=vlist_append_empty(dgr->worms);
	
	worm->startmidpoint=target1;
	worm->endmidpoint=target2;
	worm->deltalambda=deltalambda;

	assert(get_worm(dgr,get_nr_worms(dgr)-1)->startmidpoint==target1);
	assert(get_worm(dgr,get_nr_worms(dgr)-1)->endmidpoint==target2);
	assert(get_worm(dgr,get_nr_worms(dgr)-1)->deltalambda==deltalambda);

	get_vertex(dgr,target1)->refs++;
	get_vertex(dgr,target2)->refs++;
	
	return true;
}

bool diagram_remove_worm(struct diagram_t *dgr,int index)
{
	int c,target1,target2,deltalambda;
	struct worm_t *worm;

	assert(index<get_nr_worms(dgr));
	
	worm=vlist_get_element(dgr->worms,index);

	target1=worm->startmidpoint;
	target2=worm->endmidpoint;
	deltalambda=worm->deltalambda;

	assert(get_vertex(dgr,target1)->refs>0);
	assert(get_vertex(dgr,target2)->refs>0);

	for(c=target1;c<target2;c++)
	{
		struct vertex_info_t *thisvertex=get_vertex(dgr,c);

		if((thisvertex->right->j-deltalambda)<0)
		{
			printf("Removal failed! Is it physical?");
			
			exit(0);
			
			return false;
		}
	}

	vlist_remove_element(dgr->worms,index);

	for(c=target1;c<target2;c++)
	{
		struct vertex_info_t *thisvertex=get_vertex(dgr,c);

		dgr->weight/=calculate_free_propagator_weight(dgr,thisvertex->right);

		thisvertex->right->j-=deltalambda;

		dgr->weight*=calculate_free_propagator_weight(dgr,thisvertex->right);
	}

#warning FIXME: the contribution to vertices needs to be calculated as well!

	get_vertex(dgr,target1)->refs--;
	get_vertex(dgr,target2)->refs--;

	return true;
}
