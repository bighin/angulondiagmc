#include <stdlib.h>

#include "sectors.h"
#include "aux.h"

bool check_triangle_condition(struct diagram_t *dgr,struct vertex_t *thisvertex)	
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

	return result;
}

bool check_parity(struct diagram_t *dgr,struct vertex_t *thisvertex)
{
	int j1,j2,j3;

	j1=thisvertex->left->j;
	j2=thisvertex->right->j;
	j3=thisvertex->phononline->lambda;

	if(!ISEVEN(j1+j2+j3))
		return false;

	return true;
}

bool is_in_P(struct diagram_t *dgr)
{
	int c;

	for(c=0;c<get_nr_vertices(dgr);c++)
	{
		if(check_parity(dgr,get_vertex(dgr,c))==false)
			return false;
		
		if(check_triangle_condition(dgr,get_vertex(dgr,c))==false)
			return false;
	}

	return true;
}

bool is_in_E(struct diagram_t *dgr)
{
	int c;

	for(c=0;c<get_nr_vertices(dgr);c++)
	{
		struct vertex_t *thisvertex;
		int j1,j2,j3;

		thisvertex=get_vertex(dgr,c);

		j1=thisvertex->left->j;
		j2=thisvertex->right->j;
		j3=thisvertex->phononline->lambda;
		
		if(abs(j1-j2)>j3)
			return false;
	}

	return true;
}
