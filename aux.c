#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

#include "aux.h"

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
