#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <curses.h>
#include <term.h>

#include "aux.h"
#include "spline.h"

/*
	Variable length list
*/

struct vlist_t *init_vlist(size_t elementsize,int maxsz)
{
	struct vlist_t *ret;
	
	if(!(ret=malloc(sizeof(struct vlist_t))))
		return NULL;

	assert(elementsize>0);
	
	ret->elementsize=elementsize;
	ret->nelements=0;
	ret->nalloced=maxsz;

	if(!(ret->mem=malloc(elementsize*ret->nalloced)))
	{
		if(ret)
			free(ret);

		return NULL;
	}
	
	return ret;
}

void fini_vlist(struct vlist_t *lst)
{
	if(lst)
	{
		if(lst->mem)
			free(lst->mem);
	
		free(lst);
	}
}

#define CAST_TO_CPTR(ptr)	((char *)(ptr))

void *vlist_get_element(struct vlist_t *lst,int n)
{
	assert(lst);
	assert(n>=0);
	assert(n<lst->nelements);

	return CAST_TO_CPTR(lst->mem)+lst->elementsize*n;
}

void vlist_remove_element(struct vlist_t *lst,int position)
{
	void *base,*target;

	base=CAST_TO_CPTR(lst->mem)+lst->elementsize*position;
	target=CAST_TO_CPTR(lst->mem)+lst->elementsize*(position+1);

	assert(position>=0);

	assert(base>=lst->mem);
	assert(CAST_TO_CPTR(base)<CAST_TO_CPTR(lst->mem)+lst->elementsize*lst->nalloced);

	assert(target>lst->mem);
	assert(CAST_TO_CPTR(target)+lst->elementsize<=CAST_TO_CPTR(lst->mem)+lst->elementsize*lst->nalloced);

	memmove(base,target,lst->elementsize*(lst->nelements-position));
	lst->nelements--;
}

void *vlist_add_element(struct vlist_t *lst,void *element,int position)
{
	void *base,*target;

	assert(lst);
	assert(lst->nelements<lst->nalloced);

	base=CAST_TO_CPTR(lst->mem)+lst->elementsize*position;
	target=CAST_TO_CPTR(lst->mem)+lst->elementsize*(position+1);

	assert(position>=0);
	assert(position<=lst->nelements);

	assert(base>=lst->mem);
	assert(CAST_TO_CPTR(base)<CAST_TO_CPTR(lst->mem)+lst->elementsize*lst->nalloced);

	assert(target>=lst->mem);
	assert(CAST_TO_CPTR(target)+lst->elementsize<=CAST_TO_CPTR(lst->mem)+lst->elementsize*lst->nalloced);
	
	memmove(target,base,lst->elementsize*(lst->nelements-position));	
	memcpy(base,element,lst->elementsize);
	lst->nelements++;
	
	return base;
}

void *vlist_add_empty(struct vlist_t *lst,int position)
{
	void *data,*ret;

	assert(lst->elementsize>0);
	assert(lst->nelements<lst->nalloced);
	
	data=malloc(lst->elementsize);
	assert(data);
	memset(data,0,lst->elementsize);

	ret=vlist_add_element(lst,data,position);
	
	if(data)
		free(data);
	
	return ret;
}

void *vlist_append(struct vlist_t *lst,void *element)
{
	void *end;

	assert(lst);
	assert(lst->nelements<lst->nalloced);

	end=CAST_TO_CPTR(lst->mem)+(lst->elementsize*lst->nelements);

	memcpy(end,element,lst->elementsize);
	lst->nelements++;
	
	return end;
}

void *vlist_append_empty(struct vlist_t *lst)
{
	void *data,*ret;

	assert(lst->elementsize>0);
	
	data=malloc(lst->elementsize);
	assert(data);
	memset(data,0,lst->elementsize);

	ret=vlist_append(lst,data);
	
	if(data)
		free(data);
	
	return ret;
}

int vlist_get_nr_elements(struct vlist_t *lst)
{
	return lst->nelements;
}

struct vlist_t *vlist_clone(struct vlist_t *lst)
{
	struct vlist_t *ret;

	assert(vlist_get_nr_elements(lst)<lst->nalloced);

	if(!(ret=init_vlist(lst->elementsize,lst->nalloced)))
		return NULL;

	memcpy(ret->mem,lst->mem,lst->elementsize*vlist_get_nr_elements(lst));

	ret->nelements=vlist_get_nr_elements(lst);

	return ret;
}

void vlist_copy(struct vlist_t *src,struct vlist_t *dst)
{
	assert((src!=NULL)&&(dst!=NULL));
	assert(src->elementsize==dst->elementsize);
	assert(src->nalloced==dst->nalloced);

	memcpy(dst->mem,src->mem,src->elementsize*vlist_get_nr_elements(src));

	dst->nelements=src->nelements;
}

/*
	Spline interpolation, given a grid of points, essentially a nice wrapper
	over the code contained in Numerical Recipes (NR).

	Note that the code in NR is __not__ thread safe, the static keyword has to
	be removed in the splint() function.
*/

struct interpolation_t *init_interpolation(double *x,double *y,int n)
{
	struct interpolation_t *ret;
	
	if(!(ret=malloc(sizeof(struct interpolation_t))))
		return NULL;

	if(!(ret->x=malloc(sizeof(double)*n)))
	{
		if(ret)
			free(ret);
	
		return NULL;
	}

	if(!(ret->y=malloc(sizeof(double)*n)))
	{
		if(ret)
		{
			if(ret->x)
				free(ret->x);

			free(ret);
		}
		
		return NULL;
	}

	if(!(ret->y2=malloc(sizeof(double)*n)))
	{
		if(ret)
		{
			if(ret->x)
				free(ret->x);

			if(ret->y)
				free(ret->y);
		
			free(ret);
		}

		return NULL;
	}

	memcpy(ret->x,x,sizeof(double)*n);
	memcpy(ret->y,y,sizeof(double)*n);
	ret->n=n;
	
	spline(ret->x,ret->y,ret->n,0.0f,0.0f,ret->y2);

	return ret;
}

void fini_interpolation(struct interpolation_t *it)
{
	if(it)
	{
		if(it->x)
			free(it->x);

		if(it->y)
			free(it->y);

		if(it->y2)
			free(it->y2);
	
		free(it);
	}
}

void copy_interpolation(struct interpolation_t *dst,struct interpolation_t *src)
{
	assert(dst);
	assert(src);

	dst->n=src->n;

	memcpy(dst->x,src->x,sizeof(double)*dst->n);
	memcpy(dst->y,src->y,sizeof(double)*dst->n);
	memcpy(dst->y2,src->y2,sizeof(double)*dst->n);
}

double get_point(struct interpolation_t *it,double x)
{
	double y;
	
	splint(it->x,it->y,it->y2,it->n,x,&y);
	
	return y;
}

/*
	Routine to seed GSL's random number generator
*/

void seed_rng(gsl_rng *rng)
{
	char *devname="/dev/random";
	FILE *dev;
        unsigned long seed;

	if((dev=fopen(devname,"r"))!=NULL)
	{
		fread(&seed,sizeof(unsigned long),1,dev);
		fclose(dev);

		gsl_rng_set(rng,seed);
	}
	else
	{
		printf("Warning: couldn't read from %s to seed the RNG.\n",devname);
	}
}

/*
	Quick routine for double comparison
*/

bool almost_same_float(double a,double b)
{
	if((fabs(a)<10e-7)&&(fabs(b)<10e-7))
		return true;

	if((fabs(a-b)/fabs(a))<10e-7)
		return true;

	return false;
}

static char term_buffer[2048];

void init_terminal_data(void)
{
	char *termtype=getenv("TERM");
	int success;

	if(termtype==NULL)
		fprintf(stderr,"Please specify a terminal type with 'setenv TERM <yourtype>'.\n");

	success=tgetent(term_buffer,termtype);
	
	if(success<0)
		fprintf(stderr, "Could not access the termcap database.\n");

	if(success==0)
		fprintf(stderr, "Terminal type `%s' is not defined.\n",termtype);
}

int get_terminal_nr_lines(void)
{
	return tgetnum("li");
}

int isign(int x)
{
	if (x>0)
		return 1;
	
	if(x<0)
		return -1;

	return 0;
}
