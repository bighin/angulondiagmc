#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>

#include "diagrams.h"
#include "updates.h"
#include "debug.h"
#include "aux.h"
#include "phonon.h"
#include "mc.h"
#include "physics.h"

void debug_propagators(struct diagram_t *dgr)
{
	int c;
	
	for(c=0;c<get_nr_free_propagators(dgr);c++)
		printf("[%d %d] ",get_free_propagator(dgr,c)->startmidpoint,get_free_propagator(dgr,c)->endmidpoint);

	printf("\n");

	for(c=0;c<get_nr_free_propagators(dgr);c++)
		printf("[%f %f] ",get_free_propagator(dgr,c)->starttau,get_free_propagator(dgr,c)->endtau);

	printf("\n");
}

void debug_vertices(struct diagram_t *dgr)
{
	int c;

	for(c=0;c<get_nr_vertices(dgr);c++)
	{
		struct vertex_t *thisvertex=get_vertex(dgr,c);
		
		if(thisvertex)
			printf("%d %d [%d %d]\n",thisvertex->left->endmidpoint,thisvertex->right->startmidpoint,thisvertex->phononline->startmidpoint,thisvertex->phononline->endmidpoint);
		else
			printf("(null)\n");
	}
}

void debug_vertices_ext(struct diagram_t *dgr)
{
	int c,d;

	debug_propagators(dgr);

	for(c=0;c<get_nr_vertices(dgr);c++)
	{
		struct vertex_t *thisvertex=get_vertex(dgr,c);

		if(!thisvertex)
		{
			printf("(null, null)\n");
			continue;
		}

		for(d=0;d<get_nr_free_propagators(dgr);d++)
		{
			if(thisvertex->left==get_free_propagator(dgr,d))
				printf("(%d, ",d);
		}

		for(d=0;d<get_nr_free_propagators(dgr);d++)
		{
			if(thisvertex->right==get_free_propagator(dgr,d))
				printf("%d)",d);
		}
	}
}

void debug_weight_old(struct diagram_t *dgr)
{
	int c;
	double ret=1.0f;
	
	/*
		We calculate the weight associated to each free rotor line...
	*/
	
	for(c=0;c<get_nr_free_propagators(dgr);c++)
	{
		struct g0_t *g0=get_free_propagator(dgr,c);
		double en,tau;
		int j;
		
		j=g0->j;
		en=j*(j+1)-dgr->chempot;
		tau=g0->endtau-g0->starttau;

		ret*=exp(-en*tau);
		
		printf("Adding free propagator (j=%d): exp(%f) = %f\n",j,-en*tau,exp(-en*tau));
	}

	/*
		...then we add the phonon arcs...
	*/

	for(c=0;c<get_nr_phonons(dgr);c++)
	{
		struct arc_t *arc;
		double factor,timediff;

		arc=get_phonon_line(dgr,c);
		timediff=arc->endtau-arc->starttau;

		assert(timediff>=0);

		factor=chi_lambda(dgr->phonon_ctx,arc->lambda,timediff);

		ret*=factor;
		printf("Adding %f (vertex from precalculated tables, lambda=%d, length=%f)\n",factor,arc->lambda,timediff);
	}

	/*
		...and finally we consider the vertices.
	*/

	for(c=0;c<get_nr_vertices(dgr);c++)
	{
		int j1,j2,j3;
		double coupling;

		struct vertex_t *thisvertex=get_vertex(dgr,c);

		j1=thisvertex->left->j;
		j2=thisvertex->right->j;
		j3=thisvertex->phononline->lambda;

		coupling=gsl_sf_coupling_3j(2*j1,2*j2,2*j3,0,0,0)*
		         sqrtf((2.0f*j1+1)*(2.0f*j2+1)*(2.0f*j3+1)/(4.0f*M_PI));
	
		ret*=coupling;

		printf("Adding %f (coupling of %d %d %d)\n",coupling,j1,j2,j3);
	}

	printf("Final result: %f\n",ret);
}

void debug_weight(struct diagram_t *dgr)
{
	int c;
	double ret=1.0f;
	
	/*
		We calculate the weight associated to each free rotor line...
	*/
	
	for(c=0;c<get_nr_free_propagators(dgr);c++)
	{
		double localret;
		
		localret=calculate_free_propagator_weight(dgr,get_free_propagator(dgr,c));
		printf("Free propagator: %f\n",localret);

		ret*=localret;
	}

	/*
		...then we add the phonon arcs...
	*/

	for(c=0;c<get_nr_phonons(dgr);c++)
	{
		double localret;
		
		localret=calculate_arc_weight(dgr,get_phonon_line(dgr,c));
		printf("Arc: %f\n",localret);
	
		ret*=localret;
	}

	/*
		...and finally we consider the vertices.
	*/

	for(c=0;c<get_nr_vertices(dgr);c++)
	{
		double localret;
		
		localret=calculate_vertex_weight(dgr,c);
		printf("Vertex: %f\n",localret);
	
		ret*=localret;
	}

	printf("Final result: %f\n",ret);
}
