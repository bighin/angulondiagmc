#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>

#include "diagrams.h"
#include "aux.h"
#include "debug.h"
#include "selfenergies.h"
#include "stat.h"
#include "phonon.h"
#include "mc.h"
#include "physics.h"

int vertex_get_j1(struct vertex_t *vif)
{
	return vif->left->j;
}

int vertex_get_m1(struct vertex_t *vif)
{
	return vif->left->m;
}

int vertex_get_j2(struct vertex_t *vif)
{
	return vif->right->j;
}

int vertex_get_m2(struct vertex_t *vif)
{
	return vif->right->m;
}

int vertex_get_lambda(struct vertex_t *vif)
{
	return vif->phononline->lambda;
}

int vertex_get_mu(struct vertex_t *vif)
{
	return vif->phononline->mu;
}

struct diagram_t *init_diagram(struct diagram_parameters_t *dpars,bool verbose)
{
	struct diagram_t *ret;
	struct g0_t *g0;

	if(!(ret=malloc(sizeof(struct diagram_t))))
		return NULL;

	ret->mintau=0.0f;
	ret->endtau=dpars->endtau;
	ret->chempot=dpars->chempot;

	ret->n=dpars->n;

	ret->phonons=init_vlist(sizeof(struct arc_t),16*1024);
	ret->midpoints=init_vlist(sizeof(double),16*1024);
	ret->free_propagators=init_vlist(sizeof(struct g0_t),1+32*1024);
	ret->vertices=init_vlist(sizeof(struct vertex_t),32*1024);

	/*
		We add the initial, lone propagator.
	*/

	g0=vlist_append_empty(ret->free_propagators);

	g0->j=dpars->j;

	g0->startmidpoint=-1;
	g0->endmidpoint=0;

	g0->starttau=0.0f;
	g0->endtau=dpars->endtau;
	
	g0->arcs_over_me=0;
	g0->arcs_over_me_lambda1=0;

	ret->sign=1;
	
	/*
		Finally we initialize the phonon context
	*/

	ret->phonon_ctx=init_phonon_ctx(30.0f,128*1024,dpars->n,verbose);

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
		
		dgr->phonon_ctx->refs--;

		if(dgr->phonon_ctx->refs<0)
			fini_phonon_ctx(dgr->phonon_ctx);

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
		return dgr->endtau;

	if(c==-1)
		return 0.0f;

	data=vlist_get_element(dgr->midpoints,c);

	return *data;
}

struct g0_t *get_free_propagator(struct diagram_t *dgr,int c)
{
	return vlist_get_element(dgr->free_propagators,c);
}

struct vertex_t *get_vertex(struct diagram_t *dgr,int c)
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
		assert(tau==dgr->endtau);

	assert(get_midpoint(dgr,c)==tau);
}

void diagram_check_consistency(struct diagram_t *dgr)
{
	int c;

	/*
		Check that taus are correctly ordered
	*/

	assert(get_free_propagator(dgr,0)->starttau==0.0f);
	assert(get_free_propagator(dgr,get_nr_free_propagators(dgr)-1)->endtau==dgr->endtau);

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
		struct vertex_t *thisvertex=get_vertex(dgr,c);
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
		assert((c==thisvertex->phononline->startmidpoint)||(c==thisvertex->phononline->endmidpoint));
	}

	/*
		Check for consistency of angular momentum couplings,
		if we are in the physical sector
	*/

	if(configuration_is_physical(dgr)==true)
	{
		for(c=0;c<get_nr_vertices(dgr);c++)
		{
			int j1,m1,j2,m2,j3,m3;
			double coupling,epsilon;

			struct vertex_t *thisvertex=get_vertex(dgr,c);

			j1=thisvertex->left->j;
			m1=thisvertex->left->m;

			j2=thisvertex->right->j;
			m2=thisvertex->right->m;

			j3=thisvertex->phononline->lambda;
			m3=thisvertex->phononline->mu;

			coupling=gsl_sf_coupling_3j(2*j1,2*j2,2*j3,0,0,0);		
			epsilon=1e-10;

			if(fabs(coupling)<=epsilon)
			{
				printf("Wrong coupling: (%d, %d) (%d, %d) (%d, %d)\n",j1,m1,j2,m2,j3,m3);
			
				if(get_nr_vertices(dgr)<10)
					print_diagram(dgr,PRINT_TOPOLOGY|PRINT_PROPAGATORS);
			}

			assert(fabs(coupling)>epsilon);
		}
	}

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
		...and phonon lines...
	*/

	for(c=0;c<get_nr_phonons(dgr);c++)
	{
		struct arc_t *thisline=get_phonon_line(dgr,c);
	
		diagram_check_consistency_of_times(dgr,thisline->starttau,thisline->startmidpoint);
		diagram_check_consistency_of_times(dgr,thisline->endtau,thisline->endmidpoint);
	}

	/*
		Now we check the consistency of the arcs_over_me field
	*/
	
	{
		int count=0;

		for(c=0;c<get_nr_vertices(dgr);c++)
		{
			struct arc_t *thisline=get_vertex(dgr,c)->phononline;
			
			assert((c==thisline->startmidpoint)||(c==thisline->endmidpoint));

			if(thisline->startmidpoint==c)
				count++;
			else
				count--;
			
			assert(get_right_neighbour(dgr,c)->arcs_over_me==count);
		}

		assert(count==0);

		for(c=0;c<get_nr_vertices(dgr);c++)
		{
			struct arc_t *thisline=get_vertex(dgr,c)->phononline;
			
			assert((c==thisline->startmidpoint)||(c==thisline->endmidpoint));

			if(thisline->lambda>=1)
			{

				if(thisline->startmidpoint==c)
					count++;
				else
					count--;
			}
			
			assert(get_right_neighbour(dgr,c)->arcs_over_me_lambda1==count);
		}

		assert(count==0);
	}

	/*
		Conservation of the m quantum number at each vertex and parity
	*/

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

		/*
			The following silences unused variable warnings
		*/

		(void)(j1);
		(void)(j2);
		(void)(j3);
		(void)(count);

		if(configuration_is_physical(dgr)==true)
		{
			assert(count==0);

			assert(ISEVEN(j1+j2+j3));
		}
	}
	
	/*
		We check the sign
	*/

	if(diagram_weight(dgr)>0.0f)
		assert(dgr->sign==+1);
	else
		assert(dgr->sign==-1);
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

double diagram_weight(struct diagram_t *dgr)
{
	int c;
	double ret=1.0f;
	
	/*
		We calculate the weight associated to each free rotor line...
	*/
	
	for(c=0;c<get_nr_free_propagators(dgr);c++)
		ret*=calculate_free_propagator_weight(dgr,get_free_propagator(dgr,c));

	/*
		...then we add the phonon arcs...
	*/

	for(c=0;c<get_nr_phonons(dgr);c++)
		ret*=calculate_arc_weight(dgr,get_phonon_line(dgr,c));

	/*
		...and finally we consider the vertices.
	*/

	for(c=0;c<get_nr_vertices(dgr);c++)
		ret*=calculate_vertex_weight(dgr,c);

	return ret;
}

void print_chars(char ch,int n)
{
	int c;
	
	for(c=0;c<n;c++)
		printf("%c",ch);
}

/*
	Nice printout of a diagram, returns the number of printed lines.
*/

int print_diagram(struct diagram_t *dgr,int flags)
{
	int c,dashing=4;

	char *pattern;
	int pattern_len,plines;
	
	plines=0;
	pattern_len=(dashing+4)*(get_nr_free_propagators(dgr));
	pattern=malloc(sizeof(char)*pattern_len);

	for(c=0;c<pattern_len-1;c++)
		pattern[c]=' ';

	pattern[pattern_len-1]='\0';

	if(flags&PRINT_TOPOLOGY)
	{
		for(c=0;c<get_nr_phonons(dgr);c++)
		{
			struct arc_t *thisline=vlist_get_element(dgr->phonons,c);
			
			if(!(flags&PRINT_DRYRUN))
			{
				print_chars(' ',2+(dashing+4)*(1+thisline->startmidpoint));
				print_chars('_',-1+(dashing+4)*(thisline->endmidpoint-thisline->startmidpoint));
			}

			pattern[1+(dashing+4)*(1+thisline->startmidpoint)]='|';
			pattern[1+(dashing+4)*(1+thisline->endmidpoint)]='|';

			if(!(flags&PRINT_DRYRUN))
				printf("\n%s\n",pattern);

			plines+=2;
		}

		if(!(flags&PRINT_DRYRUN))
			printf("| %d|",0);

		for(c=0;c<get_nr_free_propagators(dgr);c++)
		{
			if(!(flags&PRINT_DRYRUN))
				print_chars('-',dashing);
		
			if(c+1<10)
			{
				if(!(flags&PRINT_DRYRUN))
					printf("| %d|",c+1);
			}
			else
			{
				if(!(flags&PRINT_DRYRUN))
					printf("|%d|",c+1);
			}
		}

		if(!(flags&PRINT_DRYRUN))
			printf("\n");

		plines+=1;

		if(flags&PRINT_ARCS_OVER_ME)
		{
			if(!(flags&PRINT_DRYRUN))
				printf("|  |");

			for(c=0;c<get_nr_free_propagators(dgr);c++)
			{
				char ch='0'+(get_left_neighbour(dgr,c)->arcs_over_me%10);
				
				if(!(flags&PRINT_DRYRUN))
					print_chars(ch,dashing);
		
				if(!(flags&PRINT_DRYRUN))
					printf("|  |");
			}
	
			if(!(flags&PRINT_DRYRUN))
				printf("\n");

			plines+=1;
		}

		if(!(flags&PRINT_DRYRUN))
			printf("\n\n");

		plines+=2;
	}

	if(flags&PRINT_INFO0)
	{	
		if(!(flags&PRINT_DRYRUN))
			printf("\n(diagram order: %d, local weight: %f, length: %f)\n",get_nr_phonons(dgr),diagram_weight(dgr),dgr->endtau);
		
		plines+=2;
	}

	if(flags&PRINT_PROPAGATORS)
	{
		for(c=0;c<get_nr_free_propagators(dgr);c++)
		{
			if(!(flags&PRINT_DRYRUN))
			{
				if(c==0)
					printf("|%d| %f\n",c,get_midpoint(dgr,c-1));
				else
					printf("|%d| %f ---> (lambda=%d, mu=%d)\n",c,get_midpoint(dgr,c-1),get_vertex(dgr,c-1)->phononline->lambda,get_vertex(dgr,c-1)->phononline->mu);

				printf(" |\n");
				printf(" | G0(j=%d, m=%d)\n",get_free_propagator(dgr,c)->j,get_free_propagator(dgr,c)->m);
				printf(" |\n");
			}

			plines+=4;
		}

		if(!(flags&PRINT_DRYRUN))
			printf("|%d| %f\n",1+get_nr_midpoints(dgr),dgr->endtau);

		plines++;
	}
	
	if(pattern)
		free(pattern);

	return plines;
}

void diagram_copy(struct diagram_t *src,struct diagram_t *dst)
{
	int c;

	dst->mintau=src->mintau;
	dst->endtau=src->endtau;
	dst->chempot=src->chempot;
	dst->sign=src->sign;

	dst->n=src->n;

	dst->rng_ctx=src->rng_ctx;

	vlist_copy(src->phonons,dst->phonons);
	vlist_copy(src->midpoints,dst->midpoints);
	vlist_copy(src->free_propagators,dst->free_propagators);
	vlist_copy(src->vertices,dst->vertices);
	
	for(c=0;c<get_nr_vertices(dst);c++)
	{
		struct vertex_t *thisvertex=get_vertex(dst,c);
	
		thisvertex->left=get_left_neighbour(dst,c);
		thisvertex->right=get_right_neighbour(dst,c);
	}

	for(c=0;c<get_nr_phonons(dst);c++)
	{
		struct arc_t *arc=get_phonon_line(dst,c);
	
		get_vertex(dst,arc->startmidpoint)->phononline=arc;
		get_vertex(dst,arc->endmidpoint)->phononline=arc;
	}

	dst->phonon_ctx=src->phonon_ctx;
	dst->phonon_ctx->refs++;
}

struct diagram_t *diagram_clone(struct diagram_t *src)
{
	struct diagram_t *ret;

	if(!(ret=malloc(sizeof(struct diagram_t))))
		return NULL;

	ret->phonons=init_vlist(sizeof(struct arc_t),16*1024);
	ret->midpoints=init_vlist(sizeof(double),16*1024);
	ret->free_propagators=init_vlist(sizeof(struct g0_t),1+32*1024);
	ret->vertices=init_vlist(sizeof(struct vertex_t),32*1024);

	diagram_copy(src,ret);

	assert(get_nr_free_propagators(src)==get_nr_free_propagators(ret));
	assert(get_nr_phonons(src)==get_nr_phonons(ret));
	assert(get_nr_vertices(src)==get_nr_vertices(ret));

	diagram_check_consistency(ret);

	return ret;
}
