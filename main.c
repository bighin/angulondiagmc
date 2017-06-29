#include <stdio.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>

#include "diagrams.h"
#include "updates.h"

/*
	A pointer to a GSL random number generator object

	This is NOT thread safe, and must modified if parallel code is wanted!
*/

gsl_rng *rng_ctx;

int main(void)
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
