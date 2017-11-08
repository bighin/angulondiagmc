/*
	This struct helps in collecting and (very simply) analyzing Monte Carlo samples
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

#include "stat.h"

/*
	Linear regression
*/

struct linreg_ctx_t *init_linreg_ctx(void)
{
	struct linreg_ctx_t *lct=malloc(sizeof(struct linreg_ctx_t));
	
	if(!lct)
		return NULL;
	
	lct->n=0;
	lct->sx=lct->sy=lct->sxx=lct->syy=lct->sxy=0.0f;

	return lct;
}

void fini_linreg_ctx(struct linreg_ctx_t *lct)
{
	if(lct)
		free(lct);
}

void linreg_add_entry(struct linreg_ctx_t *lct,double x,double y)
{
	lct->sx+=x;
	lct->sy+=y;
	lct->sxx+=x*x;
	lct->syy+=y*y;
	lct->sxy+=x*y;
	lct->n++;
}

double slope(struct linreg_ctx_t *lct)
{
	double a,b,xmean,ymean;
	
	xmean=lct->sx/((double)(lct->n));
	ymean=lct->sy/((double)(lct->n));

	a=lct->sxy-ymean*lct->sx-xmean*lct->sy+lct->n*xmean*ymean;
	b=lct->sxx+lct->n*xmean*xmean-2.0f*xmean*lct->sx;

	return a/b;
}

double intercept(struct linreg_ctx_t *lct)
{
	double xmean,ymean;
	
	xmean=lct->sx/((double)(lct->n));
	ymean=lct->sy/((double)(lct->n));	

	return ymean-slope(lct)*xmean;
}

double slope_error(struct linreg_ctx_t *lct)
{
	double a,b,ret=0.0f;
	
	a=intercept(lct);
	b=slope(lct);

	ret+=lct->syy+lct->n*a*a+b*b*lct->sxx;
	ret+=2.0f*(-lct->sy*a+lct->sxy*b-a*b*lct->sx);

	return sqrt(ret/((double)(lct->n)));
}

/*
	Returns a random number from a doubly truncated exponential distribution,
	i.e. a distribution with probability density

	\frac{lambda \exp(-\lambda \tau)}{\exp(-\lambda \tau_1) - \exp(-\lambda \tau_2)}

	See DTED.nb for the derivation.
*/

double doubly_truncated_exp_dist(gsl_rng *rctx,double lambda,double tau1,double tau2)
{
	double x=gsl_rng_uniform(rctx);
	double delta=-x*(1.0f-exp(-(tau2-tau1)*lambda));

	return tau1-log1p(delta)/lambda;
}

double doubly_truncated_exp_pdf(gsl_rng *rctx,double lambda,double tau1,double tau2,double tau)
{
	return lambda*exp(-lambda*tau)/(exp(-lambda*tau1)-exp(-lambda*tau2));
}
