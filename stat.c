/*
	This struct helps in collecting and (very simply) analyzing Monte Carlo samples
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

#include "stat.h"

struct samples_t *samples_init(void)
{
	struct samples_t *ret;
	
	if(!(ret=malloc(sizeof(struct samples_t))))
		return NULL;

	ret->nalloced=1024;
	ret->next=0;
	ret->data=malloc(sizeof(double)*ret->nalloced);
	
	if(!(ret->data))
	{
		if(ret)
			free(ret);
		
		return NULL;
	}
	
	return ret;
}

void samples_fini(struct samples_t *smpls)
{
	if(smpls)
	{
		if(smpls->data)
			free(smpls->data);
	
		free(smpls);
	}
}

void samples_add_entry(struct samples_t *smpls,double x)
{
	assert(smpls);

	assert(!(smpls->next>smpls->nalloced));

	if(smpls->next==smpls->nalloced)
	{
		smpls->nalloced+=1024;
		smpls->data=realloc(smpls->data,sizeof(double)*smpls->nalloced);
	}

	smpls->data[smpls->next]=x;
	smpls->next++;
}

double samples_get_average(struct samples_t *smpls)
{
	int c;
	double total=0.0f;

	assert(!(smpls->next>smpls->nalloced));

	for(c=0;c<smpls->next;c++)
		total+=smpls->data[c];
	
	return total/((double)(smpls->next));
}

double samples_get_variance(struct samples_t *smpls)
{
	int c;
	double total,average,variance,n;

	assert(!(smpls->next>smpls->nalloced));

	total=variance=0.0f;
	n=smpls->next;

	for(c=0;c<smpls->next;c++)
		total+=smpls->data[c];

	average=total/n;

	for(c=0;c<smpls->next;c++)
		variance+=pow(smpls->data[c]-average,2.0f);

	variance*=(1.0f/(n-1));

	return variance;
}

/*
	Sampling context, basic a bunch of samples_t: they can be updated
	one by one, their averages and standard deviations can be dumped
	to a tuple
*/

struct sampling_ctx_t *sampling_ctx_init(int channels)
{
	struct sampling_ctx_t *ret;
	int c;
	
	if(!(ret=malloc(sizeof(struct sampling_ctx_t))))
		return NULL;
	
	ret->smpls=malloc(sizeof(struct samples_t *)*channels);
	ret->channels=channels;
	
	for(c=0;c<ret->channels;c++)
		ret->smpls[c]=samples_init();
	
	return ret;
}

void sampling_ctx_fini(struct sampling_ctx_t *sctx)
{
	if(sctx)
	{
		int c;

		for(c=0;c<sctx->channels;c++)
			samples_fini(sctx->smpls[c]);
	
		free(sctx);
	}
}

void sampling_ctx_add_entry_to_channel(struct sampling_ctx_t *sctx,int channel,double entry)
{
	assert(channel>=0);
	assert(channel<sctx->channels);
	
	samples_add_entry(sctx->smpls[channel],entry);
}

void sampling_ctx_to_tuple(struct sampling_ctx_t *sctx,double *average,double *variance)
{
	int c;
	
	if(average)
	{
		for(c=0;c<sctx->channels;c++)
			average[c]=samples_get_average(sctx->smpls[c]);
	}

	if(variance)
	{
		for(c=0;c<sctx->channels;c++)
			variance[c]=samples_get_variance(sctx->smpls[c]);
	}
}

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
	A simple histogram
*/

struct histogram_t *init_histogram(int nbins,double width)
{
	struct histogram_t *ret=malloc(sizeof(struct histogram_t));
	
	assert(ret);
	assert(nbins>0);
	assert(width>0.0f);

	ret->sctx=sampling_ctx_init(nbins+1);
	ret->nbins=nbins;
	ret->width=width;

	assert(ret->sctx);

	return ret;
}

void histogram_add_sample(struct histogram_t *htt,double sample,double time)
{
	int targetbin;
	
	targetbin=((int)(time/htt->width));

	/* Overflow bin, the two conditions should be equivalent. */
	if(targetbin>htt->nbins)
		targetbin=htt->nbins;

	if(sample>(htt->width*htt->nbins))
		targetbin=htt->nbins;

	assert(targetbin<=htt->nbins);

	sampling_ctx_add_entry_to_channel(htt->sctx,targetbin,sample);
}

double histogram_get_bin_average(struct histogram_t *htt,int bin)
{
	assert(bin>=0);
	assert(bin<=htt->nbins);

	if(htt->sctx->smpls[bin]->next==0)
		return 0.0f;

	return samples_get_average(htt->sctx->smpls[bin]);
}

double histogram_get_bin_variance(struct histogram_t *htt,int bin)
{
	assert(bin>=0);
	assert(bin<=htt->nbins);

	if(htt->sctx->smpls[bin]->next==0)
		return 0.0f;

	return samples_get_variance(htt->sctx->smpls[bin]);
}

void fini_histogram(struct histogram_t *htt)
{
	if(htt)
	{
		sampling_ctx_fini(htt->sctx);
		
		free(htt);
	}
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
	double delta=-x+x*exp((tau1-tau2)*lambda);

	return (lambda*tau1-log1p(delta))/lambda;
}

double doubly_truncated_exp_pdf(gsl_rng *rctx,double lambda,double tau1,double tau2,double tau)
{
	return lambda*exp(-lambda*tau)/(exp(-lambda*tau1)-exp(-lambda*tau2));
}
