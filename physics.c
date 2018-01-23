#include <math.h>
#include <assert.h>
#include <gsl/gsl_sf.h>

#include "physics.h"
#include "diagrams.h"
#include "phonon.h"
#include "aux.h"

/*
	The following three functions implement the diagrammatic rules for the angulon, return:

	- the weight of a phonon arc
	- the weight of a free rotor propagator
	- the weight of a vertex
*/

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

double calculate_free_propagator_weight(struct diagram_t *dgr,struct g0_t *g0)
{
	double en,tau,phase,result;
	int j;

	/*
		Note that we may be in the unphysical sector where j < 0.
	
		We can choose the weights in the unphysical sector arbitrarily,
		the easiest way is to take the absolute value of j.
	*/

	j=abs(g0->j);
	en=j*(j+1)-dgr->chempot;
	tau=g0->endtau-g0->starttau;

	phase=1.0f;

	if(is_last_propagator(dgr,g0)==false)
		phase=pow(-1.0f,j);

	result=exp(-en*tau)*phase;

	assert((result*phase)>=0.0f);

	return result;
}

double calculate_vertex_weight(struct diagram_t *dgr,int index)
{
	int j1,j2,j3,m1,m2,m3;
	double coupling,unphysical_penalty;

	struct vertex_t *thisvertex=get_vertex(dgr,index);

	/*
		Again, in the unphysical sector one may have j < 0.

		Again, since we can define the rules in the unphysical sector arbitrarily,
		we just take the modulus.
	*/

	j1=abs(thisvertex->left->j);
	j2=abs(thisvertex->right->j);
	j3=thisvertex->phononline->lambda;


	m1=thisvertex->left->m;
	m2=thisvertex->right->m;
	m3=thisvertex->phononline->mu;

	coupling=sqrtf((2.0f*j1+1)*(2.0f*j2+1)*(2.0f*j3+1)/(4.0f*M_PI));
	coupling*=gsl_sf_coupling_3j(2*j1,2*j2,2*j3,0,0,0);
	
	/*
		The unphysical penalty is zero. This means that no unphysical states will be visited.
	*/

	unphysical_penalty=0.0f;

	if((thisvertex->left->j<0)||(thisvertex->right->j<0))
		coupling*=unphysical_penalty;

	if(index==thisvertex->phononline->startmidpoint)
	{
		if((-m1+m2+m3)==0)
			coupling*=gsl_sf_coupling_3j(2*j1,2*j2,2*j3,-2*m1,2*m2,2*m3);
		else
			coupling*=unphysical_penalty*gsl_sf_coupling_3j(2*j1,2*j2,2*j3,0,0,0);
	}
	else
	{
		if((-m1+m2-m3)==0)
			coupling*=gsl_sf_coupling_3j(2*j1,2*j2,2*j3,-2*m1,2*m2,-2*m3);
		else
			coupling*=unphysical_penalty*gsl_sf_coupling_3j(2*j1,2*j2,2*j3,0,0,0);
	}

	return coupling;
}

bool propagators_are_allowed(struct diagram_t *dgr)
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
	On the other hand, a physical configuration must have j >= 0 and |m| <= j for every rotor propagator
*/

bool configuration_is_physical(struct diagram_t *dgr)
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

	for(c=0;c<get_nr_vertices(dgr);c++)
	{
		int j1,m1,j2,m2,j3,m3,count;

		struct vertex_t *thisvertex=get_vertex(dgr,c);

		j1=thisvertex->left->j;
		m1=thisvertex->left->m;

		j2=thisvertex->right->j;
		m2=thisvertex->right->m;

		j3=thisvertex->phononline->lambda;
		m3=thisvertex->phononline->mu;

		if(c==thisvertex->phononline->startmidpoint)
			count=-m1+m2+m3;
		else
			count=-m1+m2-m3;

		if(count!=0)
			return false;
		
		if(!ISEVEN(j1+j2+j3))
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
