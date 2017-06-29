#include <stdio.h>
#include <stdbool.h>

#include "diagrams.h"
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
