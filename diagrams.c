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

struct diagram_t *init_diagram(struct diagram_cfg_t *cfg)
{
	struct diagram_t *ret;
	struct g0_t *g0;

	if(!(ret=malloc(sizeof(struct diagram_t))))
		return NULL;

	ret->mintau=0.0f;
	ret->endtau=cfg->endtau;
	ret->chempot=cfg->chempot;

	ret->n=cfg->n;

	ret->phonons=init_vlist(sizeof(struct arc_t),16*1024);
	ret->midpoints=init_vlist(sizeof(double),16*1024);
	ret->free_propagators=init_vlist(sizeof(struct g0_t),1+32*1024);
	ret->vertices=init_vlist(sizeof(struct vertex_t),32*1024);

	init_vertex_tables(ret,cfg->maxtau,1024);

	/*
		We add the initial, lone propagator.
	*/

	g0=vlist_append_empty(ret->free_propagators);

	g0->j=cfg->j;

	g0->startmidpoint=-1;
	g0->endmidpoint=0;

	g0->starttau=0.0f;
	g0->endtau=cfg->endtau;
	
	g0->arcs_over_me=0;

	/*
		We set the initial value for the diagram weight, it will
		be updated incrementally from now on.
	*/

	ret->weight=diagram_weight_non_incremental(ret);

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

		fini_vertex_tables(dgr);
	
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
	}

	/*
		Check for consistency of angular momentum couplings
	*/

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
		Finally, the weight calculate incrementally should match the non-incremental
		reference version
	*/

	{
		double w1,w2;
		
		w1=diagram_weight(dgr);
		w2=diagram_weight_non_incremental(dgr);
	
		if(!(fabs(w1-w2)<10e-7*w1))
		{
			printf("Consistency check for the incremental weight failed! (%f, %f)\n",w1,w2);
			debug_weight(dgr);
			exit(0);
		}
	
		assert(fabs(w1-w2)<10e-7*w1);
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
	}
}

double diagram_weight(struct diagram_t *dgr)
{
	return dgr->weight;
}

double diagram_weight_non_incremental(struct diagram_t *dgr)
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

		assert(timediff>=0);

		ret*=chi_lambda(dgr,arc->lambda,timediff);
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

	dst->n=src->n;

	dst->v0intercept=src->v0intercept;
	dst->v0slope=src->v0slope;
	dst->v1intercept=src->v1intercept;
	dst->v1slope=src->v1slope;

	copy_interpolation(dst->v0table,src->v0table);
	copy_interpolation(dst->v1table,src->v1table);

	dst->weight=src->weight;
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

	ret->v0table=malloc(sizeof(struct interpolation_t));
	ret->v1table=malloc(sizeof(struct interpolation_t));

	ret->v0table->x=malloc(sizeof(double)*src->v0table->n);
	ret->v0table->y=malloc(sizeof(double)*src->v0table->n);
	ret->v0table->y2=malloc(sizeof(double)*src->v0table->n);

	ret->v1table->x=malloc(sizeof(double)*src->v1table->n);
	ret->v1table->y=malloc(sizeof(double)*src->v1table->n);
	ret->v1table->y2=malloc(sizeof(double)*src->v1table->n);

	diagram_copy(src,ret);

	assert(get_nr_free_propagators(src)==get_nr_free_propagators(ret));
	assert(get_nr_phonons(src)==get_nr_phonons(ret));
	assert(get_nr_vertices(src)==get_nr_vertices(ret));

	diagram_check_consistency(ret);

	return ret;
}

bool check_triangle_condition_and_parity(struct diagram_t *dgr,struct vertex_t *thisvertex)
{
	int j1,j2,j3;
	bool result;

	j1=thisvertex->left->j;
	j2=thisvertex->right->j;
	j3=thisvertex->phononline->lambda;

	result=true;

	if(j1+j2<j3)
		result=false;

	if(j2+j3<j1)
		result=false;

	if(j3+j1<j2)
		result=false;

	/*
		This additional condition is not stricly the triangular condition,
		however we must check them because of the 3j symbols with all the
		magnetic quantum numbers to zero, that appear at every vertex.
	*/

	if(!ISEVEN(j1+j2+j3))
		result=false;

	return result;
}

bool check_triangle_condition_and_parity_from_js(int j1,int j2,int j3)
{
	bool result=true;

	if(j1+j2<j3)
		result=false;

	if(j2+j3<j1)
		result=false;

	if(j3+j1<j2)
		result=false;

	if(!ISEVEN(j1+j2+j3))
		result=false;

	return result;
}

/*
	The values of corresponding to the k integral for each line
	joining two vertices are precalculated and an interpolation is
	created.
*/

bool init_vertex_tables(struct diagram_t *dgr,double maxtau,int steps)
{
	int c;
	double n,step;
	double *x,*y0,*y1;

	n=dgr->n;
	step=maxtau/steps;

	if(!(x=malloc(sizeof(double)*steps)))
		return false;

	if(!(y0=malloc(sizeof(double)*steps)))
	{
		if(x)
			free(x);
	
		return false;
	}

	if(!(y1=malloc(sizeof(double)*steps)))
	{
		if(x)
			free(x);

		if(y0)
			free(y0);

		return false;
	}

	for(c=0;c<steps;c++)
	{
		double deltatau=step*c;

		x[c]=deltatau;
		y0[c]=chi(0,deltatau,n);
		y1[c]=chi(1,deltatau,n);
	}

	dgr->v0table=init_interpolation(x,y0,steps);
	dgr->v1table=init_interpolation(x,y1,steps);

	/*
		Now we use the same data to get an exponential fit of the function,
		which will be used when adding/removing a phonon line...
	*/

	{
		struct linreg_ctx_t *lct0,*lct1;
		
		lct0=init_linreg_ctx();
		lct1=init_linreg_ctx();

		for(c=0;c<steps;c++)
		{
			if((step*c)>3.5f)
				break;
			
			linreg_add_entry(lct0,x[c],log(y0[c]));
			linreg_add_entry(lct1,x[c],log(y1[c]));
		}

		dgr->v0slope=y0[0];
		dgr->v0intercept=slope(lct0);

		dgr->v1slope=y1[0];
		dgr->v1intercept=slope(lct1);
	
		fini_linreg_ctx(lct0);
		fini_linreg_ctx(lct1);
	}

	if(x) free(x);
	if(y0) free(y0);
	if(y1) free(y1);
	
	return true;
}

void fini_vertex_tables(struct diagram_t *dgr)
{
	fini_interpolation(dgr->v0table);
	fini_interpolation(dgr->v1table);
}

/*
	This function evaluates the function \chi_{\lambda}, as defined in my paper.

	The values are from an interpolation, which is calculated at the creation of a new diagram,
	i.e. by init_diagram().
*/

double chi_lambda(struct diagram_t *dgr,int lambda,double timediff)
{
	double ret=1.0f;
	
	switch(lambda)
	{
		case 0:
		ret*=get_point(dgr->v0table,timediff);
		break;

		case 1:
		ret*=get_point(dgr->v1table,timediff);
		break;
		
		default:
		assert(false);
	}
	
	return ret;
}

/*
	The weight of a phonon arc (i.e. chi_lambda as a function of t) does not
	have a closed form, being the result of an integral.

	However, we can approximate it as alphaeff*exp(-omega0eff*t)
*/

double get_alphaeff(struct diagram_t *dgr,int lambda)
{
	switch(lambda)
	{
		case 0:
		return fabs(dgr->v0slope);
		break;

		case 1:
		return fabs(dgr->v1slope);
		break;
	}
	
	assert(false);

	return 0.0;
}

double get_omega0eff(struct diagram_t *dgr,int lambda)
{
	switch(lambda)
	{
		case 0:
		return fabs(dgr->v0intercept);
		break;

		case 1:
		return fabs(dgr->v1intercept);
		break;
	}

	assert(false);

	return 0.0;
}
