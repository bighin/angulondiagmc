#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>

#include "diagrams.h"
#include "updates.h"
#include "stat.h"
#include "debug.h"
#include "aux.h"
#include "mc.h"
#include "graphs.h"

int main(int argc,char *argv[])
{
	init_njsummat();
	hashtable_init();

	/*
		This can go in graph.c as a unit test for the graph technology!
	*/

	{
		struct diagram_cfg_t dcfg;
		struct diagram_t *dgr;
		
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
	
		diagram_add_phonon_line(dgr,0.05,0.10,1.0,2,0);

		print_diagram(dgr,PRINT_TOPOLOGY);
		printf("<mWeight (using graphs, reference)>: %f %f\n",diagram_m_weight(dgr),diagram_m_weight_reference(dgr));

		diagram_add_phonon_line(dgr,0.35,0.65,1.0,2,0);

		print_diagram(dgr,PRINT_TOPOLOGY);
		printf("<mWeight (using graphs, reference)>: %f %f\n",diagram_m_weight(dgr),diagram_m_weight_reference(dgr));

		diagram_add_phonon_line(dgr,0.45,0.75,1.0,2,0);

		print_diagram(dgr,PRINT_TOPOLOGY);
		printf("<mWeight (using graphs, reference)>: %f %f\n",diagram_m_weight(dgr),diagram_m_weight_reference(dgr));

		while(get_nr_phonons(dgr)>0)
			diagram_remove_phonon_line(dgr,0);

		diagram_add_phonon_line(dgr,0.05,0.50,1.0,2,0);
		diagram_add_phonon_line(dgr,0.45,0.55,1.0,2,0);

		print_diagram(dgr,PRINT_TOPOLOGY);
		printf("<mWeight (using graphs, reference)>: %f %f\n",diagram_m_weight(dgr),diagram_m_weight_reference(dgr));

		while(get_nr_phonons(dgr)>0)
			diagram_remove_phonon_line(dgr,0);

		diagram_add_phonon_line(dgr,0.05,0.95,1.0,2,0);
		diagram_add_phonon_line(dgr,0.25,0.75,1.0,2,0);
		diagram_add_phonon_line(dgr,0.45,0.55,1.0,2,0);

		print_diagram(dgr,PRINT_TOPOLOGY);
		printf("<mWeight (using graphs, reference)>: %f %f\n",diagram_m_weight(dgr),diagram_m_weight_reference(dgr));

		while(get_nr_phonons(dgr)>0)
			diagram_remove_phonon_line(dgr,0);

		diagram_add_phonon_line(dgr,0.05,0.45,1.0,2,0);
		diagram_add_phonon_line(dgr,0.35,0.75,1.0,2,0);
		diagram_add_phonon_line(dgr,0.65,0.95,1.0,2,0);

		print_diagram(dgr,PRINT_TOPOLOGY);
		printf("<mWeight (using graphs, reference)>: %f %f\n",diagram_m_weight(dgr),diagram_m_weight_reference(dgr));

		while(get_nr_phonons(dgr)>0)
			diagram_remove_phonon_line(dgr,0);

		diagram_add_phonon_line(dgr,0.05,0.95,1.0,2,0);
		diagram_add_phonon_line(dgr,0.15,0.85,1.0,2,0);
		diagram_add_phonon_line(dgr,0.25,0.75,1.0,2,0);
		diagram_add_phonon_line(dgr,0.45,0.55,1.0,2,0);

		print_diagram(dgr,PRINT_TOPOLOGY);
		printf("<mWeight (using graphs, reference)>: %f %f\n",diagram_m_weight(dgr),diagram_m_weight_reference(dgr));

		while(get_nr_phonons(dgr)>0)
			diagram_remove_phonon_line(dgr,0);

		diagram_add_phonon_line(dgr,0.05,0.95,1.0,2,0);
		diagram_add_phonon_line(dgr,0.15,0.85,1.0,2,0);
		diagram_add_phonon_line(dgr,0.25,0.75,1.0,2,0);
		diagram_add_phonon_line(dgr,0.025,0.55,1.0,2,0);

		print_diagram(dgr,PRINT_TOPOLOGY);
		printf("<mWeight (using graphs, reference)>: %f %f\n",diagram_m_weight(dgr),diagram_m_weight_reference(dgr));

		return 0;
	}

	//test_graph();
	//return 0;

	//stresstest();
	//return 0;
	
	if(argc!=2)
	{
		printf("Usage: %s <inifile>\n",argv[0]);
		return 0;
	}
	
	do_diagmc(argv[1]);
	
	return 0;
}
