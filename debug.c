#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>

#include "diagrams.h"
#include "updates.h"
#include "debug.h"

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
		struct vertex_info_t *thisvertex=get_vertex(dgr,c);
		
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
		struct vertex_info_t *thisvertex=get_vertex(dgr,c);

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

void debug_weight(struct diagram_t *dgr)
{
	int c;
	
	/*
		We calculate the weight associated to each free rotor line...
	*/
	
	printf("<");
	
	for(c=0;c<get_nr_free_propagators(dgr);c++)
	{
		struct g0_t *g0=get_free_propagator(dgr,c);
		double en,tau;
		
		en=g0->j*(g0->j+1.0f)-dgr->chempot;
		tau=g0->endtau-g0->starttau;

		printf("%f (%f,%f)",-en*tau,en,tau);
	
		if(c!=get_nr_free_propagators(dgr)-1)
			printf(", ");
	}

	printf("> <");

	/*
		...then we add the phonon arcs...
	*/

	for(c=0;c<get_nr_phonons(dgr);c++)
	{
		struct arc_t *arc;
		double timediff,omega0;

		omega0=dgr->omega0;

		arc=get_phonon_line(dgr,c);
		timediff=arc->endtau-arc->starttau;

		assert(timediff>=0);

		printf("%f",-timediff*omega0);

		if(c!=get_nr_phonons(dgr)-1)
			printf(", ");
	}

	printf("> <");

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
	
		printf("%f ((%d,%d),(%d,%d),(%d,%d))",coupling,j1,m1,j2,m2,j3,m3);

		if(c!=get_nr_vertices(dgr)-1)
			printf(", ");
	}

	printf("> Total: <%f>\n",diagram_weight(dgr));
}

void stresstest(void)
{
	struct diagram_t *dgr;
	struct diagram_cfg_t cfg;

	gsl_rng *rctx;

	int x,y,z;
	int parx,pary,parz;

	if((rctx=gsl_rng_alloc(gsl_rng_mt19937))==NULL)
	{
		printf("Couldn't initialize the random number generator\n");
		return;
	}

	cfg.endtau=1.0f;
	cfg.j=2;
	cfg.m=1;
	cfg.chempot=5.0f;

	dgr=init_diagram(&cfg);

	parx=64;
	pary=128;
	parz=46;

	assert(parz<pary);

	for(x=0;x<parx;x++)
	{
		for(y=0;y<pary;y++)
		{
			double lo,hi,k;
			int lambda,mu;

			lo=gsl_rng_uniform(rctx);
			hi=gsl_rng_uniform(rctx);
			k=gsl_rng_uniform(rctx);

			lambda=gsl_rng_uniform_int(rctx,8);
			mu=gsl_rng_uniform_int(rctx,lambda+1);

			if(gsl_rng_uniform(rctx)<0.5f)
				mu=-mu;

			if(lo>hi)
			{
				double tmp;
		
				tmp=hi;
				hi=lo;
				lo=tmp;
			}

			printf("<+ (%d) (w=%f)>\n",get_nr_phonons(dgr),diagram_weight(dgr));

			diagram_add_phonon_line(dgr,lo,hi,k,lambda,mu);

			/*
				The newly added line may violate the triangle condition, in this
				case it is removed!
			*/

			{
				struct arc_t *thisline=get_phonon_line(dgr,get_nr_phonons(dgr)-1);
				struct vertex_info_t *v1,*v2;
			
				v1=get_vertex(dgr,thisline->startmidpoint);
				v2=get_vertex(dgr,thisline->startmidpoint);
			
				if((check_triangle_condition(dgr,v1)==false)||(check_triangle_condition(dgr,v2)==false))
				{
					int lastline=get_nr_phonons(dgr)-1;
					diagram_remove_phonon_line(dgr,lastline);
				}
			}

			diagram_check_consistency(dgr);
			print_diagram(dgr);
		}

		for(z=0;z<parz;z++)
		{
			if(get_nr_phonons(dgr)>0)
			{
				int target=gsl_rng_uniform_int(rctx,get_nr_phonons(dgr));
				
				diagram_remove_phonon_line(dgr,target);
				diagram_check_consistency(dgr);
			
				printf("<- (%d)>\n",get_nr_phonons(dgr));
			}
		}
	}

	while(get_nr_vertices(dgr)>0)
	{
		if(get_nr_vertices(dgr)<=19)
			print_diagram(dgr);

		printf("<- (%d)>\n",get_nr_phonons(dgr));
		diagram_remove_phonon_line(dgr,gsl_rng_uniform_int(rctx,get_nr_phonons(dgr)));
	
	}
}
