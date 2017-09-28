#ifndef __AUX_H__
#define __AUX_H__

#include <stdbool.h>
#include <gsl/gsl_rng.h>

#define MIN(a,b)	(((a)<(b))?(a):(b))
#define MAX(a,b)	(((a)>(b))?(a):(b))

#define ISEVEN(x)	(((x)%2)==0)
#define ISODD(x)	(((x)%2)==1)

struct vlist_t
{
	void *mem;

	size_t elementsize;
	int nelements,nalloced;
};

struct vlist_t *init_vlist(size_t elementsize,int maxsz);
void fini_vlist(struct vlist_t *lst);

void *vlist_get_element(struct vlist_t *lst,int n);
void vlist_remove_element(struct vlist_t *lst,int position);
void *vlist_add_element(struct vlist_t *lst,void *element,int position);
void *vlist_add_empty(struct vlist_t *lst,int position);
void *vlist_append(struct vlist_t *lst,void *element);
void *vlist_append_empty(struct vlist_t *lst);
int vlist_get_nr_elements(struct vlist_t *lst);
struct vlist_t *vlist_clone(struct vlist_t *lst);
void vlist_copy(struct vlist_t *src,struct vlist_t *dst);

struct interpolation_t
{
	double *x;
	double *y;

	int n;
	
	double *y2;
};

struct interpolation_t *init_interpolation(double *x,double *y,int n);
void fini_interpolation(struct interpolation_t *it);
void copy_interpolation(struct interpolation_t *dst,struct interpolation_t *src);
double get_point(struct interpolation_t *it,double x);

void seed_rng(gsl_rng *rng);

bool almost_same_float(double a,double b);

#endif //__AUX_H__
