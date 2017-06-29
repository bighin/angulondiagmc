#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "diagrams.h"
#include "aux.h"

#define MAX(a,b)	(((a)>(b))?(a):(b))
#define MIN(a,b)	(((a)>(b))?(b):(a))

int positive_part(int x)
{
	return (x>=0)?(x):(0);
}

bool diagram_calculate_couplings(struct diagram_t *dgr,int lo,int hi)
{
	int c;
	int localminj,localmaxj;
	int *minj,*maxj;
	
	minj=malloc(sizeof(int)*get_nr_free_propagators(dgr));
	maxj=malloc(sizeof(int)*get_nr_free_propagators(dgr));

	assert(lo>0);
	assert(hi>0);
	assert(lo<get_nr_free_propagators(dgr)-1);
	assert(hi<get_nr_free_propagators(dgr)-1);
	assert(lo<=hi);

	assert(get_free_propagator(dgr,lo)->j==-1);
	assert(get_free_propagator(dgr,lo-1)->j!=-1);

	assert(get_free_propagator(dgr,hi)->j==-1);
	assert(get_free_propagator(dgr,hi+1)->j!=-1);

	// Forward propagation of bounds
	localminj=localmaxj=get_free_propagator(dgr,lo-1)->j;
	for(c=lo;c<=hi;c++)
	{
		struct arc_t *phonon_line;
		int lambda,mu;

		phonon_line=get_vertex(dgr,c-1)->phononline;
		lambda=phonon_line->lambda;
		mu=phonon_line->mu;

		assert(lambda>=0);
		assert(abs(mu)<=lambda);

#ifdef DEBUG8
		printf("Propagating (%d, %d) to ",localminj,localmaxj);
#endif
		
		localminj=positive_part(localminj-lambda);
		localmaxj+=lambda;

		minj[c]=localminj;
		maxj[c]=localmaxj;

#ifdef DEBUG8
		printf("(%d, %d) due to lambda=%d\n",minj[c],maxj[c],lambda);
#endif
	}

	// Back propagation of bounds
	localminj=localmaxj=get_free_propagator(dgr,hi+1)->j;
	for(c=hi;c>=lo;c--)
	{
		struct arc_t *phonon_line;
		int lambda,mu;
	
		phonon_line=get_vertex(dgr,c)->phononline;
		lambda=phonon_line->lambda;
		mu=phonon_line->mu;

		assert(lambda>=0);
		assert(abs(mu)<=lambda);

#ifdef DEBUG8
		printf("Backpropagating (%d, %d) to ",localminj,localmaxj);
#endif

		localminj=positive_part(localminj-lambda);
		localmaxj=localmaxj+lambda;

		minj[c]=MAX(minj[c],localminj);
		maxj[c]=MIN(maxj[c],localmaxj);

#ifdef DEBUG8
		printf("(%d, %d) due to lambda=%d\n",minj[c],maxj[c],lambda);
#endif
	}

	// Selection rules
	for(c=lo;c<=hi;c++)
	{
		int j1,m1,j2,m2,j3,m3;
		struct randomized_list_t *lst=init_rlist();

		j1=get_free_propagator(dgr,c-1)->j;
		m1=get_free_propagator(dgr,c-1)->m;

		j2=get_phonon_line_after_propagator(dgr,c-1)->lambda;
		m2=get_phonon_line_after_propagator(dgr,c-1)->mu;

#ifdef DEBUG7
		printf("Coupling: (%d %d) and (%d %d)",j1,m1,j2,m2);
#endif

		m3=-m1-m2;

		for(j3=minj[c];j3<=maxj[c];j3++)
		{
			/*
				m3 cannot exceed j3
			*/

			if(abs(m3)>j3)
				continue;

			/*
				If all ms are zero, then the sum of all js must be even.
			*/

			if((m1==m2)&&(m2==m3)&&(m3==0))
				if(((j1+j2+j3)%2)==1)
					continue;

			/*
				Check if the angular momenta satisfy the triangular inequality
			*/

			if(abs(j1-j2)>j3)
				continue;

			if(abs(j2-j3)>j1)
				continue;

			if(abs(j3-j1)>j2)
				continue;

			rlist_add_item(lst,j3);
		}

		if(lst->nitems==0)
		{
			/*
				FIXME
			
				Probably for very big diagrams is possible that
				the algorithm is not able to find a set of couplings...
			
				There should be some fallback to c--, etc...
			*/
			
			print_diagram(dgr);
			printf("Couldn't find an allowed coupling, please debug this!\n");
			return false;
		}

		get_free_propagator(dgr,c)->j=rlist_get_random_item(lst);
		get_free_propagator(dgr,c)->m=m3;

#ifdef DEBUG7
		printf(" to (%d %d)\n",get_free_propagator(dgr,c)->j,get_free_propagator(dgr,c)->m);
#endif

		fini_rlist(lst);
	}

	if(maxj)
		free(maxj);

	if(minj)
		free(minj);
	
	return true;
}

void diagram_identify_invalid_couplings(struct diagram_t *dgr,int *first,int *last)
{
	int c;

	*first=-1;
	for(c=0;c<get_nr_free_propagators(dgr);c++)
	{
		struct g0_t *thisg0=get_free_propagator(dgr,c);

		if((thisg0->j==-1)&&(thisg0->m==-1))
		{
			*first=c;
			break;
		}
	}

	*last=-1;
	for(c=get_nr_free_propagators(dgr)-1;c>=0;c--)
	{
		struct g0_t *thisg0=get_free_propagator(dgr,c);
		
		if((thisg0->j==-1)&&(thisg0->m==-1))
		{
			*last=c;
			break;
		}
	}

	/*
		Nothing to mark as invalid, happens only when
		we are back to just one propagator.
	*/

	if((*first==-1)&&(*last==-1))
	{
		assert(get_nr_free_propagators(dgr)==1);
		return;
	}

	assert(*first!=-1);	
	assert(*last!=-1);
	assert(*first<=*last);

	for(c=*first;c<=*last;c++)
	{
		struct g0_t *thisg0=get_free_propagator(dgr,c);
		
		thisg0->j=-1;
		thisg0->m=-1;
	}
}
