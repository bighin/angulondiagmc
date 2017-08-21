#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>

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

int test_one_diagram(int index,int j1,int j2,int j3,int j4)
{
	struct diagram_cfg_t dcfg;
	struct diagram_t *dgr;
	double w1,w2;
	int mode;
	
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

	mode=PRINT_TOPOLOGY;//|PRINT_PROPAGATORS;
	dgr=init_diagram(&dcfg);

	switch(index)
	{
		/*
			One loop: bubble
		*/

		case 0:

		diagram_add_phonon_line(dgr,0.05,0.10,1.0,j1,0);

		print_diagram(dgr,mode);
	
		w1=diagram_m_weight(dgr);
		w2=diagram_m_weight_reference(dgr);
		printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
		assert(relative_distance(w1,w2)<10e-6);
		assert(sign(w1)==sign(w2));

		break;

		/*
			Two loops: two subsequent bubbles
		*/

		case 1:

		diagram_add_phonon_line(dgr,0.05,0.10,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.35,0.65,1.0,j2,0);

		print_diagram(dgr,mode);

		w1=diagram_m_weight(dgr);
		w2=diagram_m_weight_reference(dgr);
		printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
		assert(relative_distance(w1,w2)<10e-6);
		assert(sign(w1)==sign(w2));

		break;
		
		/*
			Three loops.
		*/

		case 2:

		diagram_add_phonon_line(dgr,0.05,0.10,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.35,0.65,1.0,j2,0);
		diagram_add_phonon_line(dgr,0.45,0.75,1.0,j3,0);

		print_diagram(dgr,mode);

		w1=diagram_m_weight(dgr);
		w2=diagram_m_weight_reference(dgr);
		printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
		assert(relative_distance(w1,w2)<10e-6);
		assert(sign(w1)==sign(w2));

		break;

		/*
			Two interleaved loops
		*/

		case 3:

		diagram_add_phonon_line(dgr,0.05,0.50,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.45,0.55,1.0,j2,0);

		print_diagram(dgr,mode);

		w1=diagram_m_weight(dgr);
		w2=diagram_m_weight_reference(dgr);
		printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
		assert(relative_distance(w1,w2)<10e-6);
		assert(sign(w1)==sign(w2));
		
		break;

		/*
			Three loops, sunset diagram, 10 momenta.
		*/

		case 4:

		diagram_add_phonon_line(dgr,0.05,0.95,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.25,0.75,1.0,j2,0);
		diagram_add_phonon_line(dgr,0.45,0.55,1.0,j3,0);

		print_diagram(dgr,mode);

		w1=diagram_m_weight(dgr);
		w2=diagram_m_weight_reference(dgr);
		printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
		assert(relative_distance(w1,w2)<10e-6);
		assert(sign(w1)==sign(w2));

		break;

		/*
			Three loops, 10 momenta.
		*/

		case 5:

		diagram_add_phonon_line(dgr,0.05,0.45,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.35,0.75,1.0,j2,0);
		diagram_add_phonon_line(dgr,0.65,0.95,1.0,j3,0);

		print_diagram(dgr,mode);

		w1=diagram_m_weight(dgr);
		w2=diagram_m_weight_reference(dgr);
		printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
		assert(relative_distance(w1,w2)<10e-6);
		assert(sign(w1)==sign(w2));
		
		break;
	
		/*
			Three loops, 10 momenta, contains a square.
		*/

		case 6:

		diagram_add_phonon_line(dgr,0.05,0.65,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.25,0.85,1.0,j2,0);
		diagram_add_phonon_line(dgr,0.45,0.95,1.0,j3,0);

		print_diagram(dgr,mode);

		w1=diagram_m_weight(dgr);
		w2=diagram_m_weight_reference(dgr);
		printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
		assert(relative_distance(w1,w2)<10e-6);
		assert(sign(w1)==sign(w2));
	
		break;
	
		/*
			Four loops, 13 momenta.
		*/

		case 7:

		diagram_add_phonon_line(dgr,0.05,0.95,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.15,0.85,1.0,j2,0);
		diagram_add_phonon_line(dgr,0.25,0.75,1.0,j3,0);
		diagram_add_phonon_line(dgr,0.45,0.55,1.0,j4,0);

		print_diagram(dgr,mode);

		w1=diagram_m_weight(dgr);
		w2=diagram_m_weight_reference(dgr);
		printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
		assert(relative_distance(w1,w2)<10e-6);
		assert(sign(w1)==sign(w2));
		
		break;

		/*
			Four loops, 13 momenta, contains a triangle.
		*/

		case 8:

		diagram_add_phonon_line(dgr,0.05,0.95,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.15,0.85,1.0,j2,0);
		diagram_add_phonon_line(dgr,0.25,0.75,1.0,j3,0);
		diagram_add_phonon_line(dgr,0.025,0.55,1.0,j4,0);

		print_diagram(dgr,mode);

		w1=diagram_m_weight(dgr);
		w2=diagram_m_weight_reference(dgr);
		printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
		assert(relative_distance(w1,w2)<10e-6);
		assert(sign(w1)==sign(w2));
		
		break;

		/*
			Four loops, 13 momenta,  contains a square.
		*/

		case 9:

		print_diagram(dgr,mode);

		diagram_add_phonon_line(dgr,0.05,0.95,1.0,j1,0);
		diagram_add_phonon_line(dgr,0.15,0.85,1.0,j2,0);
		diagram_add_phonon_line(dgr,0.25,0.975,1.0,j3,0);
		diagram_add_phonon_line(dgr,0.025,0.55,1.0,j4,0);

		diagram_m_weight(dgr);

		print_diagram(dgr,mode);

		w1=diagram_m_weight(dgr);
		w2=diagram_m_weight_reference(dgr);
		printf("<mWeight (using graphs, reference)>: %f %f\n",w1,w2);
		assert(relative_distance(w1,w2)<10e-6);
		assert(sign(w1)==sign(w2));
		
		break;
	}

	return 0;
}

int test_graphical_machinery(void)
{
	for(int c=0;c<=9;c++)
		test_one_diagram(c,1,2,1,2);

	for(int c=0;c<=9;c++)
		test_one_diagram(c,2,2,2,2);

	for(int c=0;c<=9;c++)
		test_one_diagram(c,1,1,1,1);

	return 0;
}
