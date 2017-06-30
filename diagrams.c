#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_sf.h>

#include "diagrams.h"
#include "aux.h"
#include "debug.h"

int vertex_get_j1(struct vertex_info_t *vif)
{
	return vif->left->j;
}

int vertex_get_m1(struct vertex_info_t *vif)
{
	return vif->left->m;
}

int vertex_get_j2(struct vertex_info_t *vif)
{
	return vif->right->j;
}

int vertex_get_m2(struct vertex_info_t *vif)
{
	return vif->right->m;
}

int vertex_get_lambda(struct vertex_info_t *vif)
{
	return vif->phononline->lambda;
}

int vertex_get_mu(struct vertex_info_t *vif)
{
	return vif->phononline->mu;
}

struct diagram_t *init_diagram(double maxtau,int j,int m,double chempot)
{
	struct diagram_t *ret;
	struct g0_t *g0;

	if(!(ret=malloc(sizeof(struct diagram_t))))
		return NULL;

	ret->mintau=0.0f;
	ret->maxtau=maxtau;
	ret->chempot=chempot;
	
	ret->phonons=init_vlist(sizeof(struct arc_t),16*1024);
	ret->midpoints=init_vlist(sizeof(double),16*1024);
	ret->free_propagators=init_vlist(sizeof(struct g0_t),1+32*1024);
	ret->vertices=init_vlist(sizeof(struct vertex_info_t),32*1024);

	/*
		Finally we add the initial, lone propagator.
	*/

	g0=vlist_append_empty(ret->free_propagators);

	g0->j=j;
	g0->m=m;

	g0->startmidpoint=-1;
	g0->endmidpoint=0;

	g0->starttau=0.0f;
	g0->endtau=maxtau;

	return ret;
}

void fini_diagram(struct diagram_t *dgr)
{
	if(dgr)
	{
		fini_vlist(dgr->phonons);
		fini_vlist(dgr->midpoints);
		fini_vlist(dgr->free_propagators);
		fini_vlist(dgr->vertices);
		
		free(dgr);
	}
}

struct arc_t *get_phonon_line(struct diagram_t *dgr,int c)
{
	return vlist_get_element(dgr->phonons,c);
}

double get_midpoint(struct diagram_t *dgr,int c)
{
	double *data;
	
	if(c==dgr->midpoints->nelements)
		return dgr->maxtau;

	if(c==-1)
		return 0.0f;

	data=vlist_get_element(dgr->midpoints,c);

	return *data;
}

struct g0_t *get_free_propagator(struct diagram_t *dgr,int c)
{
	return vlist_get_element(dgr->free_propagators,c);
}

struct vertex_info_t *get_vertex(struct diagram_t *dgr,int c)
{
	return vlist_get_element(dgr->vertices,c);
}

int get_nr_phonons(struct diagram_t *dgr)
{
	return vlist_get_nr_elements(dgr->phonons);
}

int get_nr_midpoints(struct diagram_t *dgr)
{
	return vlist_get_nr_elements(dgr->midpoints);
}

int get_nr_free_propagators(struct diagram_t *dgr)
{
	return vlist_get_nr_elements(dgr->free_propagators);
}

int get_nr_vertices(struct diagram_t *dgr)
{
	return vlist_get_nr_elements(dgr->vertices);
}

struct g0_t *get_left_neighbour(struct diagram_t *dgr,int midpoint)
{
	return vlist_get_element(dgr->free_propagators,midpoint);
}

struct g0_t *get_right_neighbour(struct diagram_t *dgr,int midpoint)
{
	return vlist_get_element(dgr->free_propagators,midpoint+1);
}

void delete_left_neighbour(struct diagram_t *dgr,int midpoint)
{
	vlist_remove_element(dgr->free_propagators,midpoint);
}

void delete_right_neighbour(struct diagram_t *dgr,int midpoint)
{
	vlist_remove_element(dgr->free_propagators,midpoint+1);
}

struct arc_t *get_phonon_line_after_propagator(struct diagram_t *dgr,int c)
{
	assert(c<get_nr_vertices(dgr));
	
	return get_vertex(dgr,c)->phononline;
}

void diagram_check_consistency_of_times(struct diagram_t *dgr,double tau,int c)
{
	if(c==-1)
		assert(tau==0.0f);

	if(c==get_nr_midpoints(dgr))
		assert(tau==dgr->maxtau);

	assert(get_midpoint(dgr,c)==tau);
}

void diagram_check_consistency(struct diagram_t *dgr)
{
	int c;

	/*
		Check that taus are correctly ordered
	*/

	assert(get_free_propagator(dgr,0)->starttau==0.0f);
	assert(get_free_propagator(dgr,get_nr_free_propagators(dgr)-1)->endtau==dgr->maxtau);

	for(c=0;c<get_nr_midpoints(dgr)-1;c++)
		assert(get_midpoint(dgr,c)<get_midpoint(dgr,c+1));

	/*
		Check the consistency of time information for free lines and phonons
	*/

	for(c=0;c<get_nr_free_propagators(dgr);c++)
	{
		if(c!=0)
		{
			assert(get_free_propagator(dgr,c)->startmidpoint==get_free_propagator(dgr,c-1)->endmidpoint);
			assert(get_free_propagator(dgr,c)->starttau==get_free_propagator(dgr,c-1)->endtau);
		}

		if(c!=get_nr_free_propagators(dgr)-1)
		{
			assert(get_free_propagator(dgr,c)->endmidpoint==get_free_propagator(dgr,c+1)->startmidpoint);
			assert(get_free_propagator(dgr,c)->endtau==get_free_propagator(dgr,c+1)->starttau);
		}
	}

	/*
		Check the consistency of vertices and phonon lines attached to them
	*/

	for(c=0;c<get_nr_vertices(dgr);c++)
	{
		struct vertex_info_t *thisvertex=get_vertex(dgr,c);
		int d=0;
		
		assert(thisvertex->left->endmidpoint==thisvertex->right->startmidpoint);
		assert(thisvertex->left->endtau==thisvertex->right->starttau);

		if(thisvertex->phononline->startmidpoint==thisvertex->left->endmidpoint)
			d++;

		if(thisvertex->phononline->endmidpoint==thisvertex->left->endmidpoint)
			d++;
	
		assert(d>0);

		d=0;

		if(thisvertex->phononline->starttau==thisvertex->left->endtau)
			d++;

		if(thisvertex->phononline->endtau==thisvertex->left->endtau)
			d++;
		
		assert(d>0);
	}

#if 0
	/*
		Check for consistency of angular momentum couplings
	*/

	for(c=0;c<get_nr_vertices(dgr);c++)
	{
		int j1,m1,j2,m2,j3,m3;
		double coupling,epsilon;

		struct vertex_info_t *thisvertex=get_vertex(dgr,c);

		j1=thisvertex->left->j;
		m1=thisvertex->left->m;

		j2=thisvertex->right->j;
		m2=thisvertex->right->m;

		j3=thisvertex->phononline->lambda;
		m3=thisvertex->phononline->mu;

		coupling=gsl_sf_coupling_3j(2*j1,2*j2,2*j3,2*m1,2*m2,2*m3);		
		epsilon=1e-10;

		if(fabs(coupling)<=epsilon)
		{
			printf("Wrong coupling: (%d, %d) (%d, %d) (%d, %d)",j1,m1,j2,m2,j3,m3);
			print_diagram(dgr);
		}

		assert(fabs(coupling)>epsilon);
	}
#endif

	/*
		More checks: again we check propagators...
	*/

	for(c=0;c<get_nr_free_propagators(dgr);c++)
	{
		struct g0_t *thisg0=get_free_propagator(dgr,c);

		diagram_check_consistency_of_times(dgr,thisg0->starttau,thisg0->startmidpoint);
		diagram_check_consistency_of_times(dgr,thisg0->endtau,thisg0->endmidpoint);
	}

	/*
		...and phonon lines
	*/

	for(c=0;c<get_nr_phonons(dgr);c++)
	{
		struct arc_t *thisline=get_phonon_line(dgr,c);
	
		diagram_check_consistency_of_times(dgr,thisline->starttau,thisline->startmidpoint);
		diagram_check_consistency_of_times(dgr,thisline->endtau,thisline->endmidpoint);
	}
}

double diagram_weight(struct diagram_t *dgr)
{
	int c;
	double ret=1.0f;
	
	/*
		We calculate the weight associated to each free rotor line...
	*/
	
#warning FIXME Probabilmente mancano dei fattori (-1)^\mu che pero non dovrebbero contare
	
	for(c=0;c<get_nr_free_propagators(dgr);c++)
	{
		struct g0_t *g0=get_free_propagator(dgr,c);

		ret*=exp(-g0->j*(g0->j+1)*(g0->endtau-g0->starttau));
	}

	/*
		...then we add the phonon arcs...
	*/

	for(c=0;c<get_nr_phonons(dgr);c++)
	{
		struct arc_t *arc;
		double timediff;
		
		arc=get_phonon_line(dgr,c);
		timediff=arc->endtau-arc->starttau;

#warning FIXME

		// FIXME : we need to decide which interaction potential to use!
		// Maybe it would be better to precalculate it...
		// \chi_\lambda(\Delta \tau) = \sum_{k} |U_\lambda (k)|^2 \exp(- \Delta \tau \omega_k)
		
		ret*=1.0f;
	}

	/*
		...and finally we consider the vertices.
	*/

	for(c=0;c<get_nr_vertices(dgr);c++)
	{
		int j1,m1,j2,m2,j3,m3;
		double coupling;

		struct vertex_info_t *thisvertex=get_vertex(dgr,c);

		j1=thisvertex->left->j;
		m1=thisvertex->left->m;

		j2=thisvertex->right->j;
		m2=thisvertex->right->m;

		j3=thisvertex->phononline->lambda;
		m3=thisvertex->phononline->mu;

		coupling=gsl_sf_coupling_3j(2*j1,2*j2,2*j3,2*m1,2*m2,2*m3)*
			 gsl_sf_coupling_3j(2*j1,2*j2,2*j3,0,0,0)*
		         sqrtf((2.0f*j1+1)*(2.0f*j2+1)*(2.0f*j3+1)/(4.0f*M_PI));
	
		ret*=coupling;
	}

	return ret;
}

void print_chars(char ch,int n)
{
	int c;
	
	for(c=0;c<n;c++)
		printf("%c",ch);
}

/*
	Nice printout of a diagram
*/

void print_diagram(struct diagram_t *dgr)
{
	int c,dashing=4;

	char *pattern;
	int pattern_len;
	
	pattern_len=(dashing+4)*(get_nr_free_propagators(dgr));
	pattern=malloc(sizeof(char)*pattern_len);

	for(c=0;c<pattern_len-1;c++)
		pattern[c]=' ';

	pattern[pattern_len-1]='\0';

	for(c=0;c<get_nr_phonons(dgr);c++)
	{
		struct arc_t *thisline=vlist_get_element(dgr->phonons,c);
		
		print_chars(' ',2+(dashing+4)*(1+thisline->startmidpoint));
		print_chars('_',-1+(dashing+4)*(thisline->endmidpoint-thisline->startmidpoint));

		pattern[1+(dashing+4)*(1+thisline->startmidpoint)]='|';
		pattern[1+(dashing+4)*(1+thisline->endmidpoint)]='|';

		printf("\n%s\n",pattern);
	}

	printf("| %d|",0);

	for(c=0;c<get_nr_free_propagators(dgr);c++)
	{
		print_chars('-',dashing);
		
		if(c+1<10)
			printf("| %d|",c+1);
		else
			printf("|%d|",c+1);
	}

	printf("\n\n");

	for(c=0;c<get_nr_free_propagators(dgr);c++)
	{
		if(c==0)
			printf("|%d| %f\n",c,get_midpoint(dgr,c-1));
		else
			printf("|%d| %f ---> (lambda=%d, mu=%d)\n",c,get_midpoint(dgr,c-1),get_vertex(dgr,c-1)->phononline->lambda,get_vertex(dgr,c-1)->phononline->mu);

		printf(" |\n");
		printf(" | G0(j=%d, m=%d)\n",get_free_propagator(dgr,c)->j,get_free_propagator(dgr,c)->m);
		printf(" |\n");
	}

	printf("|%d| %f\n",1+get_nr_midpoints(dgr),dgr->maxtau);
	
	if(pattern)
		free(pattern);
}
