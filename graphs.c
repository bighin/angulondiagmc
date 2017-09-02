#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <gsl/gsl_sf.h>

#include "diagrams.h"
#include "graphs.h"

#include "njsummat/njformul.h"
#include "njsummat/njsummat.h"
#include "murmurhash3/murmurhash3.h"

/*
	FIXME: it would be nice to have everywhere a check, to see if we are hitting the maximum
	number of lines, arcs or deltas.
*/

void reset_formula(FORMULA *f)
{
	int c,d;
	
	/*
		Warning: FORMULA is the same datastructure as in the NJFORMULA
		and NJSUMMAT packages, and so uses a weird 1-index convention in C!
	*/
	
	f->nrjs=0;
	f->nrks=0;
	f->nrsixjs=0;

	for(c=0;c<MAXJ;c++)
		f->jsigns[c]=f->jsqrts[c]=f->js[c]=0;

	for(c=0;c<MAXK;c++)
		f->ksigns[c]=f->ksqrts[c]=f->kp[c]=f->kptmp[c]=0;

	for(c=0;c<MAXSIXJ;c++)
		for(d=1;d<=6;d++)
			f->sixjs[c][d]=0;

	f->ordered=FALSE;
}

void verbose_printf(char *fmt,...)
{
	va_list ap;
	
	va_start(ap,fmt);
	//vprintf(fmt,ap);
	va_end(ap);
}

void debug_formula(FORMULA f)
{
	int c,d;
	
	verbose_printf("#js: %d\n",f.nrjs);
	verbose_printf("#ks: %d\n",f.nrks);
	verbose_printf("#6js: %d\n",f.nrsixjs);

	for(c=1;c<=f.nrjs;c++)
		verbose_printf("j(%d) = %d\n",c,f.js[c]);

	for(c=1;c<=f.nrjs;c++)
		if(f.jsigns[c]!=0)
			verbose_printf("jsign(%d) = %d\n",c,f.jsigns[c]);

	for(c=1;c<=f.nrks;c++)
		if(f.ksigns[c]!=0)
			verbose_printf("ksign(%d) = %d\n",c,f.ksigns[c]);

	for(c=1;c<=f.nrjs;c++)
		if(f.jsqrts[c]!=0)
			verbose_printf("jsqrts(%d) = %d\n",c,f.jsqrts[c]);

	for(c=1;c<=f.nrks;c++)
		if(f.ksqrts[c]!=0)
			verbose_printf("ksqrts(%d) = %d\n",c,f.ksqrts[c]);

	for(c=1;c<=f.nrsixjs;c++)
	{
		verbose_printf("sixj(");
		
		for(d=1;d<=6;d++)
		{
			verbose_printf("%d",f.sixjs[c][d]);
		
			if(d!=6)
				verbose_printf(",");
		}
		
		verbose_printf(")\n");
	}

	verbose_printf("Formula is %sordered.\n",(f.ordered==FALSE)?("NOT "):(""));
}

#define J_OR_K(a)	((a>0)?('j'):('k'))

void formula_to_wolfram(FORMULA f)
{
	int c,d;

	verbose_printf("Sum[(-1)^(");

	for(c=1;c<=f.nrjs;c++)
		if(f.jsigns[c]!=0)
			verbose_printf("%+d*j%d ",f.jsigns[c],c);

	for(c=1;c<=f.nrks;c++)
		if(f.ksigns[c]!=0)
			verbose_printf("%+d*k%d ",f.ksigns[c],c);

	verbose_printf(") ");

	for(c=1;c<=f.nrjs;c++)
		if(f.jsqrts[c]!=0)
			verbose_printf("(2 j%d + 1)^(%d/2) ",c,f.jsqrts[c]);

	for(c=1;c<=f.nrks;c++)
		if(f.ksqrts[c]!=0)
			verbose_printf("(2 k%d + 1)^(%d/2) ",c,f.ksqrts[c]);

	for(c=1;c<=f.nrsixjs;c++)
	{
		verbose_printf("SixJSymbol[{");
		
		for(d=1;d<=3;d++)
		{
			verbose_printf("%c%d",J_OR_K(f.sixjs[c][d]),abs(f.sixjs[c][d]));
		
			if(d!=3)
				verbose_printf(",");
		}

		verbose_printf("},{");

		for(d=4;d<=6;d++)
		{
			verbose_printf("%c%d",J_OR_K(f.sixjs[c][d]),abs(f.sixjs[c][d]));
		
			if(d!=6)
				verbose_printf(",");
		}

		verbose_printf("}] ");
	}

	verbose_printf("/. {");

	for(c=1;c<=f.nrjs;c++)
	{
		if((f.js[c]%2)==0)
			verbose_printf("j%d -> %d",c,f.js[c]/2);
		else
			verbose_printf("j%d -> %d.5",c,(f.js[c]-1)/2);

		if(c!=f.nrjs)
			verbose_printf(", ");
	}

	verbose_printf("}, ");

	for(c=1;c<=f.nrks;c++)
	{
		verbose_printf("{k%d, 0, 10}",c);
	
		if(c!=f.nrks)
			verbose_printf(", ");
	}

	verbose_printf("]\n");
}

void emit_6j(struct graph_t *gt,int a,int b,int c,int d,int e,int f)
{
	verbose_printf("SIXJ(");
	verbose_printf("%c%d,",J_OR_K(a),abs(a));
	verbose_printf("%c%d,",J_OR_K(b),abs(b));
	verbose_printf("%c%d,",J_OR_K(c),abs(c));
	verbose_printf("%c%d,",J_OR_K(d),abs(d));
	verbose_printf("%c%d,",J_OR_K(e),abs(e));
	verbose_printf("%c%d)\n",J_OR_K(f),abs(f));

	gt->f.nrsixjs++;
	gt->f.sixjs[gt->f.nrsixjs][1]=a;
	gt->f.sixjs[gt->f.nrsixjs][2]=b;
	gt->f.sixjs[gt->f.nrsixjs][3]=c;
	gt->f.sixjs[gt->f.nrsixjs][4]=d;
	gt->f.sixjs[gt->f.nrsixjs][5]=e;
	gt->f.sixjs[gt->f.nrsixjs][6]=f;
}

/*
	jhat = sqrt(2j+1)
*/

void emit_jhat(struct graph_t *gt,int a,int power)
{
	if((power%2)==0)
		verbose_printf("JHAT: (2*%c%d + 1)^(%d)\n",J_OR_K(a),abs(a),power/2);
	else
		verbose_printf("JHAT: (2*%c%d + 1)^(%d/2)\n",J_OR_K(a),abs(a),power);

	if(a>0)
		gt->f.jsqrts[a]+=power;
	else
		gt->f.ksqrts[-a]+=power;
}

void emit_delta(struct graph_t *gt,int a,int b)
{
	verbose_printf("DELTA(%c%d,%c%d)\n",J_OR_K(a),abs(a),J_OR_K(b),abs(b));

	gt->delta[gt->nr_deltas][0]=a;
	gt->delta[gt->nr_deltas][1]=b;
	gt->nr_deltas++;
}

void emit_summation(struct graph_t *gt,int newk)
{
	assert(newk<0);

	verbose_printf("SUM_{%c%d}\n",J_OR_K(newk),abs(newk));

	gt->f.nrks++;
}

void emit_phase(struct graph_t *gt,int a,int mul)
{
	verbose_printf("PHASE: (-1)^(%d %c%d)\n",mul,J_OR_K(a),abs(a));

	if(a>0)
		gt->f.jsigns[a]+=mul;
	else
		gt->f.ksigns[-a]+=mul;
}

int find_arc_with_vertex(struct graph_t *gt,int vertex)
{
	int c;
	
	for(c=0;c<gt->nr_arcs;c++)
		if((gt->arcs[c][0]==vertex)||(gt->arcs[c][1]==vertex))
			return c;
	
	return -1;
}

void show_graph(struct graph_t *gt)
{
	int c;
	
	verbose_printf("-\n");
	
	for(c=0;c<gt->nr_vertices;c++)
	{
		int index=find_arc_with_vertex(gt,c);

		verbose_printf("|\n|\n| %c%d\n|\n|______ %c%d\n",J_OR_K(gt->lines[c]),abs(gt->lines[c]),J_OR_K(gt->arcs[index][2]),abs(gt->arcs[index][2]));
	}

	verbose_printf("|\n|\n| j%d\n|\n|\n-\n\n",gt->lines[c]);
}

void remove_vertex(struct graph_t *gt,int position,int mode)
{	
	/*
		SAVE_LEFT means: overwrite the right neighbour.
	*/
	
	if(mode==SAVE_LEFT)
		memmove(&gt->lines[position+1],&gt->lines[position+2],sizeof(int)*(MAX_LINES-position-2));

	/*
		SAVE_RIGHT means: overwrite the left neighbour.
	*/

	else if(mode==SAVE_RIGHT)
		memmove(&gt->lines[position],&gt->lines[position+1],sizeof(int)*(MAX_LINES-position-1));
}

void remove_arc(struct graph_t *gt,int index)
{
	int c,d,start,end;

	start=gt->arcs[index][0];
	end=gt->arcs[index][1];

	remove_vertex(gt,end,SAVE_RIGHT);	
	remove_vertex(gt,start,SAVE_LEFT);

	memmove(&gt->arcs[index],&gt->arcs[index+1],sizeof(int)*3*(MAX_ARCS-index-1));

	gt->nr_vertices-=2;
	gt->nr_lines-=2;
	gt->nr_arcs--;

	for(c=0;c<gt->nr_arcs;c++)
	{
		for(d=0;d<=1;d++)
		{
			if(gt->arcs[c][d]>end)
				gt->arcs[c][d]--;

			if(gt->arcs[c][d]>start)
				gt->arcs[c][d]--;
		}
	}
}

bool has_nspanning_arc(struct graph_t *gt,int vertex,int span)
{
	int c;

	for(c=0;c<gt->nr_arcs;c++)
		if((gt->arcs[c][0]==vertex)&&(gt->arcs[c][1]==(vertex+span)))
			return true;

	return false;
}

bool has_bubble(struct graph_t *gt,int vertex)
{
	return has_nspanning_arc(gt,vertex,1);
}

bool has_triangle(struct graph_t *gt,int vertex)
{
	return has_nspanning_arc(gt,vertex,2);
}

bool has_square(struct graph_t *gt,int vertex)
{
	return has_nspanning_arc(gt,vertex,3);
}

void remove_bubble(struct graph_t *gt,int vertex)
{
	int index,arcstart,arcend,before,after;

	index=find_arc_with_vertex(gt,vertex);

	arcstart=gt->arcs[index][0];
	arcend=gt->arcs[index][1];

	before=arcstart;
	after=arcend+1;

	emit_phase(gt,gt->lines[before],1);
	emit_phase(gt,gt->lines[before+1],1);
	emit_phase(gt,gt->arcs[index][2],1);

	emit_delta(gt,gt->lines[before],gt->lines[after]);
	emit_jhat(gt,gt->lines[before],-2);

	remove_arc(gt,index);
}

void remove_triangle(struct graph_t *gt,int vertex)
{
	int arcstart,arcend,before,after,firstinside,secondinside,freeline,closedline;

	int index=find_arc_with_vertex(gt,vertex);

	arcstart=gt->arcs[index][0];
	arcend=gt->arcs[index][1];

	before=gt->lines[arcstart];
	after=gt->lines[arcend+1];

	firstinside=gt->lines[arcstart+1];
	secondinside=gt->lines[arcstart+2];
	
	freeline=gt->arcs[find_arc_with_vertex(gt,vertex+1)][2];
	closedline=gt->arcs[index][2];

	emit_phase(gt,firstinside,1);
	emit_phase(gt,secondinside,1);
	emit_phase(gt,freeline,1);

	emit_6j(gt,before,after,freeline,
	        secondinside,firstinside,closedline);

	remove_arc(gt,index);
}

void remove_square(struct graph_t *gt,int vertex)
{
	int arcstart,arcend,before,after,firstinside,secondinside,thirdinside,freeline1,freeline2,closedline,newmomentum;

	int index=find_arc_with_vertex(gt,vertex);

	arcstart=gt->arcs[index][0];
	arcend=gt->arcs[index][1];

	before=gt->lines[arcstart];
	after=gt->lines[arcend+1];

	firstinside=gt->lines[arcstart+1];
	secondinside=gt->lines[arcstart+2];
	thirdinside=gt->lines[arcstart+3];

	freeline1=gt->arcs[find_arc_with_vertex(gt,vertex+1)][2];
	freeline2=gt->arcs[find_arc_with_vertex(gt,vertex+2)][2];
	closedline=gt->arcs[index][2];

	gt->maxk++;
	newmomentum=-gt->maxk;

	emit_6j(gt,freeline1,before,newmomentum,
	        closedline,secondinside,firstinside);

	emit_6j(gt,freeline2,after,newmomentum,
	        closedline,secondinside,thirdinside);


	emit_phase(gt,newmomentum,1);
	emit_phase(gt,closedline,1);
	emit_phase(gt,firstinside,1);
	emit_phase(gt,secondinside,1);
	emit_phase(gt,thirdinside,1);
	emit_phase(gt,freeline1,1);
	emit_phase(gt,freeline2,1);

	emit_jhat(gt,newmomentum,2);

	emit_summation(gt,newmomentum);

	gt->lines[arcstart+2]=newmomentum;

	remove_arc(gt,index);
}

void exchange_lines(struct graph_t *gt,int v1,int v2)
{	
	int a1,a2,before,middle,after,line1,line2,newmomentum;
	
	assert(v2==(v1+1));

	a1=find_arc_with_vertex(gt,v1);
	a2=find_arc_with_vertex(gt,v2);

	assert(a1!=a2);

	if(gt->arcs[a1][0]==v1)
		gt->arcs[a1][0]=v2;

	if(gt->arcs[a1][1]==v1)
		gt->arcs[a1][1]=v2;

	if(gt->arcs[a2][0]==v2)
		gt->arcs[a2][0]=v1;

	if(gt->arcs[a2][1]==v2)
		gt->arcs[a2][1]=v1;

	before=gt->lines[v1];
	middle=gt->lines[v1+1];
	after=gt->lines[v1+2];
	
	line1=gt->arcs[a1][2];
	line2=gt->arcs[a2][2];

	gt->maxk++;
	newmomentum=-gt->maxk;

	emit_phase(gt,line1,1);
	emit_phase(gt,line2,1);
	emit_phase(gt,middle,1);
	emit_phase(gt,newmomentum,1);

	emit_jhat(gt,newmomentum,2);

	emit_6j(gt,after,line2,newmomentum,
	        before,line1,middle);

	emit_summation(gt,newmomentum);

	gt->lines[v1+1]=newmomentum;
}

void k_to_j(struct graph_t *gt,int kindex,int jindex)
{
	int c,d;
	
	assert((kindex<0)&&(jindex>0));
	
	for(c=0;c<gt->nr_lines;c++)
	{
		if(gt->lines[c]==kindex)
			gt->lines[c]=jindex;
	}

	for(c=0;c<gt->nr_arcs;c++)
	{
		if(gt->arcs[c][2]==kindex)
			gt->arcs[c][2]=jindex;
	}

	for(c=0;c<gt->nr_deltas;c++)
	{
		if(gt->delta[c][0]==kindex)
			gt->delta[c][0]=jindex;

		if(gt->delta[c][1]==kindex)
			gt->delta[c][1]=jindex;
	}

	for(c=1;c<=gt->f.nrsixjs;c++)
	{
		for(d=1;d<=6;d++)
		{
			if(gt->f.sixjs[c][d]==kindex)
				gt->f.sixjs[c][d]=jindex;
		}
	}

	gt->f.jsigns[jindex]+=gt->f.ksigns[-kindex];
	gt->f.jsqrts[jindex]+=gt->f.ksqrts[-kindex];

	for(c=-kindex;c<=gt->f.nrks;c++)
	{
		gt->f.ksigns[c]=gt->f.ksigns[c+1];
		gt->f.ksqrts[c]=gt->f.ksqrts[c+1];
	}

	gt->f.nrks--;

	assert(gt->f.nrks>=0);
}

void k_to_k(struct graph_t *gt,int kindex1,int kindex2)
{
	int c,d;

	assert((kindex1<0)&&(kindex2<0));
	assert(kindex1!=kindex2);

	if(abs(kindex1)>abs(kindex2))
	{
		k_to_k(gt,kindex2,kindex1);
		return;
	}
	
	for(c=0;c<gt->nr_lines;c++)
	{
		if(gt->lines[c]==kindex2)
			gt->lines[c]=kindex1;
	}

	for(c=0;c<gt->nr_arcs;c++)
	{
		if(gt->arcs[c][2]==kindex2)
			gt->arcs[c][2]=kindex1;
	}

	for(c=0;c<gt->nr_deltas;c++)
	{
		if(gt->delta[c][0]==kindex2)
			gt->delta[c][0]=kindex1;

		if(gt->delta[c][1]==kindex2)
			gt->delta[c][1]=kindex1;
	}

	for(c=1;c<=gt->f.nrsixjs;c++)
	{
		for(d=1;d<=6;d++)
		{
			if(gt->f.sixjs[c][d]==kindex1)
				gt->f.sixjs[c][d]=kindex2;
		}
	}

	gt->f.ksigns[-kindex1]+=gt->f.ksigns[-kindex2];
	gt->f.ksqrts[-kindex1]+=gt->f.ksqrts[-kindex2];

	for(c=-kindex2;c<=gt->f.nrks;c++)
	{
		gt->f.ksigns[c]=gt->f.ksigns[c+1];
		gt->f.ksqrts[c]=gt->f.ksqrts[c+1];
	}

	gt->f.nrks--;

	assert(gt->f.nrks>=0);
}

double evaluate_graph(struct graph_t *gt,bool debugswap)
{
	int c;
	
	/*
		Here we could take care of zerolines immediately.
		This would likely result in simpler expressions, however is not strictly necessary.
	*/

	reset_formula(&gt->f);
	
	/*
		If the debugswap flag is enabled a swap is performed between the first two
		lines. The operation should have no effect, and thus can be used for debug purposes.
	*/
	
	if(debugswap==true)
	{
		if(gt->nr_vertices>=2)
		{
			if(find_arc_with_vertex(gt,0)!=find_arc_with_vertex(gt,1))
			{
					exchange_lines(gt,0,1);
			
					verbose_printf("evaluate_graph(): performed line exchange (debugswap flag is enabled).\n");
			}
		}
	}

	verbose_printf("Starting (initial) formula debug\n");
	debug_formula(gt->f);
	verbose_printf("Ending (initial) formula debug\n\n");

	startover:

	verbose_printf("Starting (over?) simplification loop:\n");
	show_graph(gt);

	verbose_printf("Looking for bubbles:\n");

	for(c=0;c<(gt->nr_vertices-1);c++)
	{
		if(has_bubble(gt,c)==true)
		{
			verbose_printf("Found a bubble @ %d, removing it\n",c);
			
			remove_bubble(gt,c);
			goto startover;
		}
	}

	verbose_printf("Looking for triangles:\n");

	for(c=0;c<(gt->nr_vertices-2);c++)
	{
		if(has_triangle(gt,c)==true)
		{
			verbose_printf("Found a triangle @ %d, removing it\n",c);

			remove_triangle(gt,c);
			goto startover;
		}
	}

	/*
		Here it would be nice to implement an overlap euristic
		to look for the best 4-cycle to remove!
	*/

	verbose_printf("Looking for squares:\n");

	for(c=0;c<(gt->nr_vertices-3);c++)
	{
		if(has_square(gt,c)==true)
		{
			verbose_printf("Found a square @ %d, removing it\n",c);

			remove_square(gt,c);
			goto startover;
		}
	}

	if(gt->nr_arcs>0)
	{
		int c,n;
		
		/*
			Look for a 5-loop, and simplify it by exchange...

			Alternatively, we look for a 6-loop, 7-loop, 8-loop and so on!
		*/

		for(n=4;n<gt->nr_vertices;n++)
		{
			verbose_printf("Looking for (%d)-loops:\n",n+1);
			
			for(c=0;c<(gt->nr_vertices-n);c++)
			{
				if(has_nspanning_arc(gt,c,n)==true)
				{
					verbose_printf("Found a (%d)-loop @ %d, converting it to a (%d)-loop by exchange\n",n+1,c,n);

					exchange_lines(gt,c,c+1);
					goto startover;
				}
			}
		}
	}

	assert(gt->nr_arcs==0);
	
	/*
		Finally we have a closed loop, whose value is (2j + 1).
	*/
	
	emit_jhat(gt,gt->lines[0],2);

	/*
		We are done simplifying, let's carry out the sums, starting from the deltas!
	*/
	
	for(c=0;c<gt->nr_deltas;c++)
	{
		int first,second;
		
		first=gt->delta[c][0];
		second=gt->delta[c][1];
		
		if(first==second)
			continue;
		
		if((first>0)&&(second>0))
		{
			if(gt->js[first]!=gt->js[second])
				return 0.0f;
		}

		if((first>0)&&(second<0))
		{
			k_to_j(gt,second,first);
		}

		if((first<0)&&(second>0))
		{
			k_to_j(gt,first,second);
		}

		if((first<0)&&(second<0))
		{
			k_to_k(gt,first,second);
		}
	}
	
	/*
		After taking care of the deltas, we are just left with a sums
		over a number of 6j symbol, will will give the desidred coefficient,
		and will be calculated by means of NJSUMMAT.

		The formula f has already been setup, we just need to fill the js,
		remembering that they should be multiplied by 2, as per NJSUMMAT documentation.
	*/

	gt->f.nrjs=gt->maxj;

	for(c=1;c<=gt->f.nrjs;c++)
		gt->f.js[c]=2*gt->js[c];

	verbose_printf("Starting (final) formula debug\n");
	debug_formula(gt->f);
	verbose_printf("Ending (final) formula debug\n\n");

	return evaluate_formula(&gt->f);
}

void diagram_to_graph(struct diagram_t *dgr,struct graph_t *gt)
{
	int c;
	
	gt->nr_vertices=get_nr_free_propagators(dgr)-1;
	gt->nr_lines=get_nr_free_propagators(dgr);
	gt->nr_arcs=get_nr_phonons(dgr);
	gt->maxj=gt->maxk=gt->nr_deltas=0;

	for(c=0;c<get_nr_free_propagators(dgr);c++)
	{
		gt->maxj++;
		gt->lines[c]=gt->maxj;
		gt->js[gt->maxj]=get_free_propagator(dgr,c)->j;
	}

	for(c=0;c<get_nr_phonons(dgr);c++)
	{
		gt->arcs[c][0]=get_phonon_line(dgr,c)->startmidpoint;
		gt->arcs[c][1]=get_phonon_line(dgr,c)->endmidpoint;

		gt->maxj++;
		gt->arcs[c][2]=gt->maxj;
		gt->js[gt->maxj]=get_phonon_line(dgr,c)->lambda;
	}
}

/*
	The following code implements a very simple hashtable, to avoid
	calculating the value of the same graph over and over...
*/

struct hashentry_t hashtable[HASHTABLE_ENTRIES];
int hashtable_lookups,hashtable_matches;

void hashtable_init(void)
{
	int c;
	
	for(c=0;c<HASHTABLE_ENTRIES;c++)
		hashtable[c].valid=HASHTABLE_INVALID_ENTRY;
	
	printf("Initalized graph hashtable, with %d entries, %ld KB\n",HASHTABLE_ENTRIES,HASHTABLE_ENTRIES*sizeof(struct graph_t)/1024);

	hashtable_lookups=hashtable_matches=0;
}

void hasthtable_show_stats(void)
{
	printf("Hashtable lookups: %d (%f%% matches)\n",hashtable_lookups,100.f*((float)(hashtable_matches))/((float)(hashtable_lookups)));
}

bool graphs_are_equivalent(struct graph_t *a,struct graph_t *b)
{
	if(a->nr_lines!=b->nr_lines)
		return false;
	
	if(a->nr_arcs!=b->nr_arcs)
		return false;

	if(a->maxj!=b->maxj)
		return false;
	
	if(memcmp(a->lines,b->lines,sizeof(int)*a->nr_lines)!=0)
		return false;

	if(memcmp(a->arcs,b->arcs,3*sizeof(int)*a->nr_arcs)!=0)
		return false;

	if(memcmp(a->js,b->js,sizeof(int)*(a->nr_lines+a->nr_arcs))!=0)
		return false;
	
	return true;
}

bool hashtable_probe(struct graph_t *gt,double *value,int *hashindex)
{
	unsigned int h;

	h=0xCAFEBABE;
	h^=murmur3_hash((void *)(gt->lines),sizeof(int)*gt->nr_lines);	
	h^=murmur3_hash((void *)(gt->arcs),3*sizeof(int)*gt->nr_arcs);
	h^=murmur3_hash((void *)(gt->js),sizeof(int)*(gt->nr_lines+gt->nr_arcs));
	h=h%HASHTABLE_ENTRIES;

	*hashindex=h;

	hashtable_lookups++;

	if(hashtable[*hashindex].valid==HASHTABLE_VALID_ENTRY)
	{
		if(graphs_are_equivalent(gt,&hashtable[*hashindex].gt)==true)
		{
			*value=hashtable[*hashindex].value;
			hashtable_matches++;

			return true;
		}
	}
	
	return false;
}

void hashtable_insert(struct graph_t *gt,double value,int hashindex)
{
	struct graph_t *dst=&hashtable[hashindex].gt;

	dst->nr_lines=gt->nr_lines;
	dst->nr_arcs=gt->nr_arcs;
	dst->maxj=gt->maxj;
		
	memcpy(dst->lines,gt->lines,sizeof(int)*gt->nr_lines);
	memcpy(dst->arcs,gt->arcs,3*sizeof(int)*gt->nr_arcs);
	memcpy(dst->js,gt->js,sizeof(int)*(gt->nr_lines+gt->nr_arcs));

	hashtable[hashindex].valid=HASHTABLE_VALID_ENTRY;
	hashtable[hashindex].value=value;
}

double diagram_m_weight(struct diagram_t *dgr,bool use_hashtable)
{
	struct graph_t gt;
	double value;

	int hashindex;

	/*
		FIXME: are these memset() call needed? Would gt and gt.f be undefined otherwise?
	*/

	memset(&gt,0,sizeof(struct graph_t));
	memset(&gt.f,0,sizeof(FORMULA));

	diagram_to_graph(dgr,&gt);

	if((use_hashtable==true)&&(hashtable_probe(&gt,&value,&hashindex)==true))
	{
		assert(value==evaluate_graph(&gt,false));

		return value;
	}

	value=evaluate_graph(&gt,false);

	if((gt.f.nrjs>MAXJ)||(gt.f.nrks>MAXK)||(gt.f.nrsixjs>MAXSIXJ))
	{
		printf("The current diagram translates to a formula exceeding the limits of NJSUMMAT.");
		printf(" Please tune the limits in njformul.h\n");
		exit(0);
	}

	//formula_to_wolfram(gt.f);

	if(use_hashtable==true)
	{
		memset(&gt.f,0,sizeof(FORMULA));

		hashtable_insert(&gt,value,hashindex);
	}

	return value;
}

double diagram_m_weight_reference(struct diagram_t *dgr)
{
	int *pjs[128],*pms[128];
	int c,i,j;
	
	double total=0.0f;

	for(c=i=0;c<get_nr_free_propagators(dgr);c++)
	{
		pjs[i]=&(get_free_propagator(dgr,c)->j);
		pms[i++]=&(get_free_propagator(dgr,c)->m);
	}

	for(c=0;c<get_nr_phonons(dgr);c++)
	{
		pjs[i]=&(get_phonon_line(dgr,c)->lambda);
		pms[i++]=&(get_phonon_line(dgr,c)->mu);
	}

	for(c=0;c<i;c++)
		*pms[c]=-*pjs[c];

	printf("Doing the full evaluation with %d momenta\n",i);

	while(true)
	{		
		double coupling=1.0f;

		/*
			The coupling is calculated for this particular combination of ms...
		*/

#ifdef DEBUG8
		printf("\n\n\nCALCULATING COUPLING!\n");
#endif

		for(c=0;c<get_nr_vertices(dgr);c++)
		{
			int j1,m1,j2,m2,j3,m3;

			struct vertex_t *thisvertex=get_vertex(dgr,c);

			j1=thisvertex->left->j;
			m1=thisvertex->left->m;

			j2=thisvertex->right->j;
			m2=thisvertex->right->m;

			j3=thisvertex->phononline->lambda;
			m3=thisvertex->phononline->mu;

#ifdef DEBUG8
			printf("ThreeJSymbol[{%d,%d},{%d,%d},{%d,%d}] = %f\n",j1,m1,j2,m2,j3,m3,gsl_sf_coupling_3j(2*j1,2*j2,2*j3,2*m1,2*m2,2*m3));
#endif

			/*
				Note the sign convention!
			*/

			assert((c==thisvertex->phononline->startmidpoint)||(c==thisvertex->phononline->endmidpoint));

			if(c==thisvertex->phononline->startmidpoint)
				coupling*=gsl_sf_coupling_3j(2*j1,2*j2,2*j3,-2*m1,2*m2,2*m3);
			else
				coupling*=gsl_sf_coupling_3j(2*j1,2*j2,2*j3,-2*m1,2*m2,-2*m3);
		}

		/*
			And now the phases...
		
			Note: the last phase is not taken into account, because it needs to be
			identified with the first one, within this context, hence the (-1) in the sum indices.
		*/

		for(c=0;c<get_nr_free_propagators(dgr)-1;c++)
		{
			struct g0_t *g0=get_free_propagator(dgr,c);

#ifdef DEBUG8
			printf("Adding (free) phase: (-1)^(%d - %d) ==> %f\n",g0->j,g0->m,pow(-1.0f,g0->j-g0->m));
#endif

			coupling*=pow(-1.0f,g0->j-g0->m);
		}

		for(c=0;c<get_nr_phonons(dgr);c++)
		{
			struct arc_t *arc=get_phonon_line(dgr,c);

#ifdef DEBUG8
			printf("Adding (phonon) phase: (-1)^(%d - %d) ==> %f\n",arc->lambda,arc->mu,pow(-1.0f,arc->lambda-arc->mu));
#endif
	
			coupling*=pow(-1.0f,arc->lambda-arc->mu);
		}
		
		/*
			Finally the Kronecker delta between the first and last free propagator momenta,
			since they must be the same...
		*/

		{
			struct g0_t *first,*last;
			
			first=get_free_propagator(dgr,0);
			last=get_free_propagator(dgr,get_nr_free_propagators(dgr)-1);
		
			if(first->j!=last->j)
				coupling=0.0f;

			if(first->m!=last->m)
				coupling=0.0f;
		}
		
		/*
			...and the is summed to the total result.
		*/

		total+=coupling;

#ifdef DEBUG7

		if(fabs(coupling)>10e-6)
		{
			printf(">>> ");

			for(j=0;j<i;j++)
			{
				printf("%d ",*pms[j]);
			}
		
			printf(" --> %f (cumulative: %f)\n",coupling,total);
		}

#endif
		/*
			We increment one of the ms, such that we visit
			all the possible combinations...
		*/

		for(j=0;j<i;j++)
		{		
			if(*pms[j]<*pjs[j])
			{
				*pms[j]=*pms[j]+1;
				break;
			}
			else
			{
				*pms[j]=-*pjs[j];
			}
		}

		/*
			...and we exit the loop if we have reached the end.
		*/

		if(j==i)
			break;
	}

	printf("Completed a full evaluation with %d momenta, total: %f\n",i,total);

	return total;
}

void init_njsummat(void)
{
	comp_fac();
}
