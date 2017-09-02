#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>

#include "diagrams.h"
#include "updates.h"
#include "graphs.h"
#include "tests.h"

double relative_distance(double w1,double w2)
{
	return fabs(fabs(w1)-fabs(w2))/fabs(w1);
}

int sign(double x)
{
	if(fabs(x)<10e-6)
		return 0;

	if(x>0.0f)
		return 1;

	return -1;
}

/*
	The same function as double diagram_m_weight(struct diagram_t *dgr)
	implemented in graphs.c, with the hashtable logic removed.

	However, the debugswap flag is enabled when calling evaluate_graph(), such that
	a line exchange is performed, to verify the equivalency of a diagram
	after such operation is performed.
*/

double diagram_m_weight_with_swap(struct diagram_t *dgr)
{
	struct graph_t gt;
	double value;

	memset(&gt,0,sizeof(struct graph_t));

	diagram_to_graph(dgr,&gt);

	value=evaluate_graph(&gt,true);

	formula_to_wolfram(gt.f);
	
	return value;
}

#define TEST_GRAPHS		(1)
#define TEST_REFERENCE_ALGO	(2)
#define TEST_SWAPPED_GRAPHS	(4)

int test_one_diagram(int index,int j1,int j2,int j3,int j4,int flags)
{
	struct diagram_cfg_t dcfg;
	struct diagram_t *dgr;
	double w1,w2,w3;
	int mode;
	
	dcfg.j=2;

	dcfg.endtau=1.0f;
	dcfg.chempot=5.75f;

	dcfg.c0=0.05f;
	dcfg.c1=0.05f;
	dcfg.c2=0.05f;

	dcfg.omega0=16.0f;
	dcfg.omega1=16.0f;
	dcfg.omega2=16.0f;

	mode=PRINT_TOPOLOGY;//|PRINT_PROPAGATORS;
	dgr=init_diagram(&dcfg);

	dgr->rng_ctx=gsl_rng_alloc(gsl_rng_mt19937);
	assert(dgr->rng_ctx!=NULL);

	switch(index)
	{
		/*
			One loop: bubble
		*/

		case 0:

		diagram_add_phonon_line(dgr,0.05,0.10,1.0,j1,0);

		break;

		/*
			Two loops: two subsequent bubbles
		*/

		case 1:

		diagram_add_phonon_line(dgr,0.05,0.10,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.35,0.65,1.0,j2,0);

		break;
		
		/*
			Three loops.
		*/

		case 2:

		diagram_add_phonon_line(dgr,0.05,0.10,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.35,0.65,1.0,j2,0);
		diagram_add_phonon_line(dgr,0.45,0.75,1.0,j3,0);

		break;

		/*
			Two interleaved loops
		*/

		case 3:

		diagram_add_phonon_line(dgr,0.05,0.50,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.45,0.55,1.0,j2,0);
		
		break;

		/*
			Three loops, sunset diagram, 10 momenta.
		*/

		case 4:

		diagram_add_phonon_line(dgr,0.05,0.95,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.25,0.75,1.0,j2,0);
		diagram_add_phonon_line(dgr,0.45,0.55,1.0,j3,0);

		break;

		/*
			Three loops, 10 momenta.
		*/

		case 5:

		diagram_add_phonon_line(dgr,0.05,0.45,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.35,0.75,1.0,j2,0);
		diagram_add_phonon_line(dgr,0.65,0.95,1.0,j3,0);
		
		break;
	
		/*
			Three loops, 10 momenta, contains a square.
		*/

		case 6:

		diagram_add_phonon_line(dgr,0.05,0.65,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.25,0.85,1.0,j2,0);
		diagram_add_phonon_line(dgr,0.45,0.95,1.0,j3,0);
	
		break;
	
		/*
			Four loops, 13 momenta.
		*/

		case 7:

		diagram_add_phonon_line(dgr,0.05,0.95,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.15,0.85,1.0,j2,0);
		diagram_add_phonon_line(dgr,0.25,0.75,1.0,j3,0);
		diagram_add_phonon_line(dgr,0.45,0.55,1.0,j4,0);
		
		break;

		/*
			Four loops, 13 momenta, contains a triangle.
		*/

		case 8:

		diagram_add_phonon_line(dgr,0.05,0.95,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.15,0.85,1.0,j2,0);
		diagram_add_phonon_line(dgr,0.25,0.75,1.0,j3,0);
		diagram_add_phonon_line(dgr,0.025,0.55,1.0,j4,0);
		
		break;

		/*
			Four loops, 13 momenta,  contains a square.
		*/

		case 9:

		diagram_add_phonon_line(dgr,0.05,0.95,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.15,0.85,1.0,j2,0);
		diagram_add_phonon_line(dgr,0.25,0.975,1.0,j3,0);
		diagram_add_phonon_line(dgr,0.025,0.55,1.0,j4,0);
		
		break;
	}

	print_diagram(dgr,mode);

	w1=w2=w3=0.0f; // To silence GCC (wrong) warnings...

	if(flags&TEST_GRAPHS)		w1=diagram_m_weight(dgr,false);
	if(flags&TEST_REFERENCE_ALGO)	w2=diagram_m_weight_reference(dgr);
	if(flags&TEST_SWAPPED_GRAPHS)	w3=diagram_m_weight_with_swap(dgr);

	if((flags&TEST_GRAPHS)&&(flags&TEST_REFERENCE_ALGO)&&(flags&TEST_SWAPPED_GRAPHS))
		printf("<mWeight (using graphs, reference, additional swap)>: %f %f %f\n",w1,w2,w3);

	if((flags&TEST_GRAPHS)&&(flags&TEST_REFERENCE_ALGO))
	{
		if(!(flags&TEST_SWAPPED_GRAPHS))
			printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);

		assert(relative_distance(w1,w2)<10e-6);
		assert(sign(w1)==sign(w2));
	}

	if((flags&TEST_GRAPHS)&&(flags&TEST_SWAPPED_GRAPHS))
	{
		if(!(flags&TEST_REFERENCE_ALGO))
			printf("<mWeight (using graphs, additional swap)>: %f %f\n",w1,w3);

		assert(relative_distance(w1,w3)<10e-6);
		assert(sign(w1)==sign(w3));
	}

	if((flags&TEST_SWAPPED_GRAPHS)&&(flags&TEST_REFERENCE_ALGO))
	{
		if(!(flags&TEST_GRAPHS))
			printf("<mWeight (reference, additional swap)>: %f %f\n",w2,w3);

		assert(relative_distance(w2,w3)<10e-6);
		assert(sign(w2)==sign(w3));
	}

	if(dgr->rng_ctx)
		gsl_rng_free(dgr->rng_ctx);

	fini_diagram(dgr);

	return 0;
}

int test_graphical_machinery(void)
{
	int c,flags;

	flags=TEST_GRAPHS|TEST_REFERENCE_ALGO|TEST_SWAPPED_GRAPHS;

	for(c=0;c<=9;c++)
		test_one_diagram(c,1,2,1,2,flags);

	for(c=0;c<=9;c++)
		test_one_diagram(c,2,2,2,2,flags);

	for(c=0;c<=9;c++)
		test_one_diagram(c,1,1,1,1,flags);

	return 0;
}
