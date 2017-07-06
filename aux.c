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
		return NULL;
	
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

void *vlist_get_element(struct vlist_t *lst,int n)
{
	assert(lst);
	assert(n>=0);
	assert(n<lst->nelements);
	
	return lst->mem+lst->elementsize*n;
}

void vlist_remove_element(struct vlist_t *lst,int position)
{
	void *base,*target;

	base=lst->mem+lst->elementsize*position;
	target=lst->mem+lst->elementsize*(position+1);

	assert(position>=0);

	assert(base>=lst->mem);
	assert(base<lst->mem+lst->elementsize*lst->nalloced);

	assert(target>lst->mem);
	assert(target+lst->elementsize<=lst->mem+lst->elementsize*lst->nalloced);

	memmove(base,target,lst->elementsize*(lst->nelements-position));
	lst->nelements--;
}

void *vlist_add_element(struct vlist_t *lst,void *element,int position)
{
	void *base,*target;

	assert(lst);
	assert(lst->nelements<lst->nalloced);

	base=lst->mem+lst->elementsize*position;
	target=lst->mem+lst->elementsize*(position+1);

	assert(position>=0);
	assert(position<=lst->nelements);

	assert(base>=lst->mem);
	assert(base<lst->mem+lst->elementsize*lst->nalloced);

	assert(target>=lst->mem);
	assert(target+lst->elementsize<=lst->mem+lst->elementsize*lst->nalloced);
	
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

	end=lst->mem+(lst->elementsize*lst->nelements);

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

/*
	A simpler list, with the possibility of requiring a random element.
*/

struct randomized_list_t *init_rlist(void)
{
	struct randomized_list_t *ret;
	
	if(!(ret=malloc(sizeof(struct randomized_list_t))))
		return NULL;
	
	ret->nitems=0;
	ret->nalloced=ret->blocksize=16*1024;

	ret->items=malloc(sizeof(int)*ret->nalloced);
	
	return ret;
}

void fini_rlist(struct randomized_list_t *lst)
{
	if(lst)
		free(lst);
}

void rlist_add_item(struct randomized_list_t *lst,int item)
{
	if(lst->nitems==lst->nalloced)
	{
		lst->nalloced+=lst->blocksize;
		lst->items=realloc(lst->items,sizeof(int)*lst->nalloced);
	}

	lst->items[lst->nitems]=item;
	lst->nitems++;
}

int rlist_get_random_item(struct randomized_list_t *lst)
{
	extern gsl_rng *rng_ctx;
	
	assert(lst->nitems>0);
	
	return lst->items[gsl_rng_uniform_int(rng_ctx,lst->nitems)];
}

int rlist_get_elements(struct randomized_list_t *lst)
{
	return lst->nitems;
}

void rlist_remove_element(struct randomized_list_t *lst,int position)
{
	void *base,*target;

	base=lst->items+sizeof(int)*position;
	target=lst->items+sizeof(int)*(position+1);

	memmove(base,target,sizeof(int)*(lst->nitems-position));
	lst->nitems--;
}

#warning TESTME

int rlist_pop_random_item(struct randomized_list_t *lst)
{
	extern gsl_rng *rng_ctx;
	int target,result;
	
	assert(lst->nitems>0);
	
	target=gsl_rng_uniform_int(rng_ctx,lst->nitems);
	result=lst->items[target];

	rlist_remove_element(lst,target);

	return result;
}
