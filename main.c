#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

#include "diagrams.h"
#include "updates.h"

/*
	A pointer to a GSL random number generator object

	This is NOT thread safe, and must modified if parallel code is wanted!
*/

gsl_rng *rng_ctx;

void stresstest(void)
{
	struct diagram_t *dgr;
	int x,y,z;
	int parx,pary,parz;

	if((rng_ctx=gsl_rng_alloc(gsl_rng_mt19937))==NULL)
	{
		printf("Couldn't initialize the random number generator\n");
		return;
	}

	dgr=init_diagram(1.0f,7,6,0.5f);

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

			lo=gsl_rng_uniform(rng_ctx);
			hi=gsl_rng_uniform(rng_ctx);
			k=gsl_rng_uniform(rng_ctx);

			lambda=gsl_rng_uniform_int(rng_ctx,8);
			mu=gsl_rng_uniform_int(rng_ctx,lambda+1);

			if(gsl_rng_uniform(rng_ctx)<0.5f)
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
			diagram_check_consistency(dgr);
		}

		for(z=0;z<parz;z++)
		{
			if(get_nr_phonons(dgr)>0)
			{
				int target=gsl_rng_uniform_int(rng_ctx,get_nr_phonons(dgr));
				
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
		diagram_remove_phonon_line(dgr,gsl_rng_uniform_int(rng_ctx,get_nr_phonons(dgr)));
	
	}
}

int main(void)
{
	stresstest();
	
	return 0;
}

int main_old(void)
{
	struct diagram_t *dgr;

	if((rng_ctx=gsl_rng_alloc(gsl_rng_mt19937))==NULL)
	{
		printf("Couldn't initialize the random number generator\n");
		return 0;
	}
		
	dgr=init_diagram(1.0f,7,6,0.5f);
	print_diagram(dgr);

	diagram_add_phonon_line(dgr,0.5,0.6,1,3,2);
	print_diagram(dgr);
	diagram_check_consistency(dgr);

	diagram_add_phonon_line(dgr,0.3,0.9,1,4,-1);
	print_diagram(dgr);
	diagram_check_consistency(dgr);

	diagram_add_phonon_line(dgr,0.05,0.95,1,2,0);
	print_diagram(dgr);
	diagram_check_consistency(dgr);

	diagram_remove_phonon_line(dgr,0);
	print_diagram(dgr);
	diagram_check_consistency(dgr);

	diagram_remove_phonon_line(dgr,0);
	print_diagram(dgr);
	diagram_check_consistency(dgr);

	diagram_update_length(dgr,1.9f);
	print_diagram(dgr);
	diagram_check_consistency(dgr);

	diagram_remove_phonon_line(dgr,0);
	print_diagram(dgr);
	diagram_check_consistency(dgr);

	fini_diagram(dgr);
	
	gsl_rng_free(rng_ctx);

	return 0;
}
