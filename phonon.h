#ifndef __PHONON_H__
#define __PHONON_H__

#include "aux.h"

struct phonon_ctx_t
{
	int nsteps;
	
	/*
		Boson density
	*/

	double n;

	/*
		Tables containing, respectively, the time, the
		phonon distribution \chi_\lambda (\Delta t)
		for \lambda=0 and \lambda=1 and the normalized
		cumulative phonon distribution.
	*/

	double *x,*y0,*y1;
	double *ncy0,*ncy1;

	/*
		The normalization constants
	*/

	double c0,c1;

	/*
		The interpolations
	*/
	
	struct interpolation_t *y0int,*y1int;
	struct interpolation_t *invncy0int,*invncy1int;
	
	/*
		Every time a diagram using this context is
		created (chiefly by copying or cloning) refs
		is increased.
	
		Every time one of such diagrams is destroyed,
		refs is decreased. When it reaches zero, it
		means that there are no diagrams using this
		context, so we can free it.
	*/

	int refs;
};

struct phonon_ctx_t *init_phonon_ctx(double maxtau,int nsteps,double n,bool verbose);
void fini_phonon_ctx(struct phonon_ctx_t *ctx);

double chi_lambda(struct phonon_ctx_t *ctx,int lambda,double deltatau);

double phonon_dist(gsl_rng *rctx,struct phonon_ctx_t *ctx,int lambda);
double phonon_pdf(struct phonon_ctx_t *ctx,int lambda,double deltatau);

#endif //__PHONON_H__
