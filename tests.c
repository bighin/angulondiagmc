#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>

#include "diagrams.h"
#include "updates.h"
#include "graphs.h"
#include "tests.h"

int test_graphical_machinery(void)
{
	struct diagram_cfg_t dcfg;
	struct diagram_t *dgr;
	double w1,w2;
	
	dcfg.j=2;
	dcfg.m=0;

	dcfg.endtau=1.0f;
	dcfg.chempot=5.75f;

	dcfg.c0=0.05f;
	dcfg.c1=0.05f;
	dcfg.c2=0.05f;

	dcfg.omega0=16.0f;
	dcfg.omega1=16.0f;
	dcfg.omega2=16.0f;

	dgr=init_diagram(&dcfg);

	//goto onlysquare;

	/*
		One loop: bubble
	*/

	diagram_add_phonon_line(dgr,0.05,0.10,1.0,2,0);

	print_diagram(dgr,PRINT_TOPOLOGY);
	
	w1=diagram_m_weight(dgr);
	w2=diagram_m_weight_reference(dgr);
	printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
	assert((fabs(w1-w2)/w1)<10e-6);

	/*
		Two loops: two subsequent bubbles
	*/

	diagram_add_phonon_line(dgr,0.35,0.65,1.0,2,0);

	print_diagram(dgr,PRINT_TOPOLOGY);

	w1=diagram_m_weight(dgr);
	w2=diagram_m_weight_reference(dgr);
	printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
	assert((fabs(w1-w2)/w1)<10e-6);
	/*
		Three loops.
	*/

	diagram_add_phonon_line(dgr,0.45,0.75,1.0,2,0);

	print_diagram(dgr,PRINT_TOPOLOGY);

	w1=diagram_m_weight(dgr);
	w2=diagram_m_weight_reference(dgr);
	printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
	assert((fabs(w1-w2)/w1)<10e-6);

	/*
		Two interleaved loops
	*/

	while(get_nr_phonons(dgr)>0)
		diagram_remove_phonon_line(dgr,0);

	diagram_add_phonon_line(dgr,0.05,0.50,1.0,2,0);
	diagram_add_phonon_line(dgr,0.45,0.55,1.0,2,0);

	print_diagram(dgr,PRINT_TOPOLOGY);

	w1=diagram_m_weight(dgr);
	w2=diagram_m_weight_reference(dgr);
	printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
	assert((fabs(w1-w2)/w1)<10e-6);

	/*
		Three loops, sunset diagram, 10 momenta.
	*/

	while(get_nr_phonons(dgr)>0)
		diagram_remove_phonon_line(dgr,0);

	diagram_add_phonon_line(dgr,0.05,0.95,1.0,2,0);
	diagram_add_phonon_line(dgr,0.25,0.75,1.0,2,0);
	diagram_add_phonon_line(dgr,0.45,0.55,1.0,2,0);

	print_diagram(dgr,PRINT_TOPOLOGY);

	w1=diagram_m_weight(dgr);
	w2=diagram_m_weight_reference(dgr);
	printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
	assert((fabs(w1-w2)/w1)<10e-6);

	/*
		Three loops, 10 momenta.
	*/

	while(get_nr_phonons(dgr)>0)
		diagram_remove_phonon_line(dgr,0);

	diagram_add_phonon_line(dgr,0.05,0.45,1.0,2,0);
	diagram_add_phonon_line(dgr,0.35,0.75,1.0,2,0);
	diagram_add_phonon_line(dgr,0.65,0.95,1.0,2,0);

	print_diagram(dgr,PRINT_TOPOLOGY);

	w1=diagram_m_weight(dgr);
	w2=diagram_m_weight_reference(dgr);
	printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
	assert((fabs(w1-w2)/w1)<10e-6);
	
	/*
		Three loops, 10 momenta, contains a square.
	*/

	while(get_nr_phonons(dgr)>0)
		diagram_remove_phonon_line(dgr,0);

	diagram_add_phonon_line(dgr,0.05,0.65,1.0,2,0);
	diagram_add_phonon_line(dgr,0.25,0.85,1.0,2,0);
	diagram_add_phonon_line(dgr,0.45,0.95,1.0,2,0);

	print_diagram(dgr,PRINT_TOPOLOGY);

	w1=diagram_m_weight(dgr);
	w2=diagram_m_weight_reference(dgr);
	printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
	assert((fabs(w1-w2)/w1)<10e-6);
	
	/*
		Four loops, 13 momenta.
	*/

	while(get_nr_phonons(dgr)>0)
		diagram_remove_phonon_line(dgr,0);

	diagram_add_phonon_line(dgr,0.05,0.95,1.0,2,0);
	diagram_add_phonon_line(dgr,0.15,0.85,1.0,2,0);
	diagram_add_phonon_line(dgr,0.25,0.75,1.0,2,0);
	diagram_add_phonon_line(dgr,0.45,0.55,1.0,2,0);

	print_diagram(dgr,PRINT_TOPOLOGY);

	w1=diagram_m_weight(dgr);
	w2=diagram_m_weight_reference(dgr);
	printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
	assert((fabs(w1-w2)/w1)<10e-6);

	/*
		Four loops, 13 momenta, contains a triangle.
	*/

	while(get_nr_phonons(dgr)>0)
		diagram_remove_phonon_line(dgr,0);

	diagram_add_phonon_line(dgr,0.05,0.95,1.0,2,0);
	diagram_add_phonon_line(dgr,0.15,0.85,1.0,2,0);
	diagram_add_phonon_line(dgr,0.25,0.75,1.0,2,0);
	diagram_add_phonon_line(dgr,0.025,0.55,1.0,2,0);

	print_diagram(dgr,PRINT_TOPOLOGY);

	w1=diagram_m_weight(dgr);
	w2=diagram_m_weight_reference(dgr);
	printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
	assert((fabs(w1-w2)/w1)<10e-6);

	/*
		Four loops, 13 momenta,  contains a square.
	*/

	onlysquare:

	print_diagram(dgr,PRINT_TOPOLOGY);

	diagram_add_phonon_line(dgr,0.05,0.95,1.0,2,0);
	diagram_add_phonon_line(dgr,0.15,0.85,1.0,2,0);
	diagram_add_phonon_line(dgr,0.25,0.975,1.0,2,0);
	diagram_add_phonon_line(dgr,0.025,0.55,1.0,2,0);

	diagram_m_weight(dgr);

	print_diagram(dgr,PRINT_TOPOLOGY);

	w1=diagram_m_weight(dgr);
	w2=diagram_m_weight_reference(dgr);
	printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
	assert((fabs(w1-w2)/w1)<10e-6);

	return 0;
}
