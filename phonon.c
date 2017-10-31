#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <assert.h>

#include "phonon.h"
#include "aux.h"
#include "stat.h"
#include "selfenergies.h"

void test_phonon_ctx(struct phonon_ctx_t *ctx);

/*
	The values of corresponding to the k integral for each line
	joining two vertices are precalculated and an interpolation is
	created.
*/

struct phonon_ctx_t *init_phonon_ctx(double maxtau,int nsteps,double n,bool verbose)
{
	struct phonon_ctx_t *ret;

	int c,nsteps0,nsteps1;;
	double step;

	assert(nsteps>0);

	if(!(ret=malloc(sizeof(struct phonon_ctx_t))))
		return NULL;

	ret->n=n;
	ret->nsteps=nsteps;
	step=maxtau/nsteps;

	if(verbose)
	{
		printf("Precalculating phonon distributions...");
		fflush(stdout);
	}

	ret->x=malloc(sizeof(double)*nsteps);
	ret->y0=malloc(sizeof(double)*nsteps);
	ret->y1=malloc(sizeof(double)*nsteps);
	ret->ncy0=malloc(sizeof(double)*nsteps);
	ret->ncy1=malloc(sizeof(double)*nsteps);

	if((!ret->x)||(!ret->y0)||(!ret->y1)||(!ret->ncy0)||(!ret->ncy1))
	{
		if(ret)
		{
			if(ret->x)	free(ret->x);
			if(ret->y0)	free(ret->y0);
			if(ret->y1)	free(ret->y1);
			if(ret->ncy0)	free(ret->ncy0);
			if(ret->ncy1)	free(ret->ncy1);

			free(ret);
		}
		
		return NULL;
	}

	/*
		We fill the tables...
	*/

	for(c=0;c<nsteps;c++)
	{
		double deltatau=step*c;

		ret->x[c]=deltatau;
		ret->y0[c]=calculate_chi_lambda(0,deltatau,n);
		ret->y1[c]=calculate_chi_lambda(1,deltatau,n);
	}

	/*
		...then we transform them to a cumulative distribution...
	*/

	ret->ncy0[0]=ret->ncy1[0]=0.0f;

	for(c=1;c<nsteps;c++)
	{
		ret->ncy0[c]=ret->ncy0[c-1]+0.5f*step*(ret->y0[c-1]+ret->y0[c]);
		ret->ncy1[c]=ret->ncy1[c-1]+0.5f*step*(ret->y1[c-1]+ret->y1[c]);
	}

	/*
		...and finally we normalize them!
	*/

	ret->c0=ret->ncy0[nsteps-1];
	ret->c1=ret->ncy1[nsteps-1];

	for(c=0;c<nsteps;c++)
	{
		ret->ncy0[c]/=ret->ncy0[nsteps-1];
		ret->ncy1[c]/=ret->ncy1[nsteps-1];
	}

	/*
		At last we create the interpolations...
	*/

	ret->y0int=init_interpolation(ret->x,ret->y0,nsteps);
	ret->y1int=init_interpolation(ret->x,ret->y1,nsteps);

	/*
		...and the interpolations for the inverse of the normalized, cumulative distribution.
		
		Since most of datapoints will be 1.0f, we can use a smaller number of steps, saving memory and time.
	*/

	nsteps0=nsteps1=c;

	for(c=0;c<nsteps;c++)
	{
		if(ret->ncy0[c]>=1.0f)
		{
			nsteps0=c;
			break;
		}
	}

	for(c=0;c<nsteps;c++)
	{
		if(ret->ncy1[c]>=1.0f)
		{
			nsteps1=c;
			break;
		}
	}

	ret->invncy0int=init_interpolation(ret->ncy0,ret->x,nsteps0);
	ret->invncy1int=init_interpolation(ret->ncy1,ret->x,nsteps1);

	ret->refs=0;

	if(verbose)
		printf(" Done!\n");

	return ret;
}

void fini_phonon_ctx(struct phonon_ctx_t *ctx)
{
	if(ctx)
	{
		if(ctx->x) free(ctx->x);
		if(ctx->y0) free(ctx->y0);
		if(ctx->y1) free(ctx->y1);
		if(ctx->ncy0) free(ctx->ncy0);
		if(ctx->ncy1) free(ctx->ncy1);

		fini_interpolation(ctx->y0int);
		fini_interpolation(ctx->y1int);

		fini_interpolation(ctx->invncy0int);
		fini_interpolation(ctx->invncy1int);

		free(ctx);
	}
}

/*
	This function evaluates the function \chi_{\lambda}, as defined in my paper.

	The values are from a precalculated interpolation.
*/

double chi_lambda(struct phonon_ctx_t *ctx,int lambda,double deltatau)
{
	double ret=1.0f;
	
	switch(lambda)
	{
		case 0:
		ret*=get_point(ctx->y0int,deltatau);
		break;

		case 1:
		ret*=get_point(ctx->y1int,deltatau);
		break;
		
		default:
		assert(false);
	}
	
	return ret;
}

double phonon_dist(gsl_rng *rctx,struct phonon_ctx_t *ctx,int lambda)
{
        double x=gsl_rng_uniform_pos(rctx);
	struct interpolation_t *it;

	assert((lambda==0)||(lambda==1));

	switch(lambda)
	{
		case 0:
		it=ctx->invncy0int;
		break;

		case 1:
		it=ctx->invncy1int;
		break;
		
		/*
			To silence silly compiler warnings.
		*/
		
		default:
		it=NULL;
		break;
	}

	return get_point(it,x);
}

double phonon_pdf(struct phonon_ctx_t *ctx,int lambda,double deltatau)
{
	double c;

	assert((lambda==0)||(lambda==1));

	switch(lambda)
	{
		case 0:
		c=ctx->c0;
		break;

		case 1:
		c=ctx->c1;
		break;

		/*
			To silence silly compiler warnings.
		*/
	
		default:
		c=1.0f;
		break;
	}

	return chi_lambda(ctx,lambda,deltatau)/c;
}

void test_phonon_ctx(struct phonon_ctx_t *ctx)
{
	FILE *out;
	gsl_rng *rctx=gsl_rng_alloc(gsl_rng_taus);
	int lambda;
	
	/*
		Test #1: histogram generation vs. normalized distribution
	*/
	
	for(lambda=0;lambda<=1;lambda++)
	{
		int d;
		double deltatau;
		char fname[1024];
				
		snprintf(fname,1024,"histogram.%d.test.dat",lambda);
		out=fopen(fname,"w+");

		for(d=0;d<2000000;d++)
		{	
			fprintf(out,"%f\n",phonon_dist(rctx,ctx,lambda));
		}

		if(out)
			fclose(out);

		snprintf(fname,1024,"phononpdf.%d.test.dat",lambda);
		out=fopen(fname,"w+");

		for(deltatau=0.0f;deltatau<=10;deltatau+=0.00125f)
			fprintf(out,"%f %f\n",deltatau,phonon_pdf(ctx,lambda,deltatau));

		if(out)
			fclose(out);
	}
}