#ifndef __STAT_H__
#define __STAT_H__

#include <gsl/gsl_histogram.h>

#include "mc.h"

struct linreg_ctx_t
{
	int n;
	double sx,sy,sxx,syy,sxy;
};

struct linreg_ctx_t *init_linreg_ctx(void);
void fini_linreg_ctx(struct linreg_ctx_t *lct);
void linreg_add_entry(struct linreg_ctx_t *lct,double x,double y);
double slope(struct linreg_ctx_t *lct);
double intercept(struct linreg_ctx_t *lct);
double slope_error(struct linreg_ctx_t *lct);

double doubly_truncated_exp_dist(gsl_rng *rctx,double lambda,double tau1,double tau2);
double doubly_truncated_exp_pdf(gsl_rng *rctx,double lambda,double tau1,double tau2,double tau);

#endif
