#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include "histograms.h"

struct simple_histogram_t *init_simple_histogram(int nbins,double binwidth)
{
	struct simple_histogram_t *ret;

	assert(nbins>0);
	assert(binwidth>0);

	if(!(ret=malloc(sizeof(struct simple_histogram_t))))
		return NULL;

	ret->means=malloc(sizeof(double)*nbins);
	ret->m2s=malloc(sizeof(double)*nbins);
	ret->ns=malloc(sizeof(int)*nbins);

	if((!ret->means)||(!ret->m2s)||(!ret->ns))
	{
		fini_simple_histogram(ret);
		return NULL;
	}

	for(int i=0;i<nbins;i++)
	{
		ret->means[i]=ret->m2s[i]=0.0f;
		ret->ns[i]=0;
	}

	ret->nbins=nbins;
	ret->binwidth=binwidth;

	return ret;
}

void fini_simple_histogram(struct simple_histogram_t *sh)
{
	if(sh)
	{
		if(sh->means)
			free(sh->means);

		if(sh->m2s)
			free(sh->m2s);

		if(sh->ns)
			free(sh->ns);
	}
}

int simple_histogram_get_bin(struct simple_histogram_t *sh,double t)
{
	int index=t/sh->binwidth;
	
	if(index<0)
		return 0;
	
	if(index>=sh->nbins)
		return sh->nbins-1;
	
	return index;
}
	
void simple_histogram_add_entry(struct simple_histogram_t *sh,int index,double value)
{
	double delta1,delta2;
	
	sh->ns[index]++;
	delta1=value-sh->means[index];
	sh->means[index]+=delta1/sh->ns[index];
	delta2=value-sh->means[index];
	sh->m2s[index]+=delta1*delta2;
}

void simple_histogram_add_time_sample(struct simple_histogram_t *sh,double t,double value)
{
	int index=simple_histogram_get_bin(sh,t);
	
	for(int c=0;c<sh->nbins;c++)
		simple_histogram_add_entry(sh,c,(c==index)?(value):(0.0f));
}

double simple_histogram_get_mean(struct simple_histogram_t *sh,int index)
{
	return sh->means[index];
}

double simple_histogram_get_variance(struct simple_histogram_t *sh,int index)
{
	return sh->m2s[index]/(sh->ns[index])/(sh->ns[index]-1);
}

struct block_histogram_t *init_block_histogram(int nbins,double binwidth,int blocksize)
{
	struct block_histogram_t *ret;

	assert(nbins>0);
	assert(binwidth>0);

	if(!(ret=malloc(sizeof(struct block_histogram_t))))
		return NULL;

	ret->means=malloc(sizeof(double)*nbins);
	ret->variances=malloc(sizeof(double)*nbins);
	ret->ns=malloc(sizeof(int)*nbins);

	if((!ret->means)||(!ret->variances)||(!ret->ns))
	{
		fini_block_histogram(ret);
		return NULL;
	}

	for(int i=0;i<nbins;i++)
	{
		ret->means[i]=ret->variances[i]=0.0f;
		ret->ns[i]=0;
	}

	ret->blocks=malloc(sizeof(double *)*nbins);
	for(int i=0;i<nbins;i++)
	{
		if(!(ret->blocks[i]=malloc(sizeof(double)*blocksize)))
		{
			fini_block_histogram(ret);
			return NULL;	
		}
	}

	if(!(ret->iblock=malloc(sizeof(int)*nbins)))
	{
		fini_block_histogram(ret);
		return NULL;
	}

	for(int i=0;i<nbins;i++)
		ret->iblock[i]=0;

	ret->blocksize=blocksize;
	
	ret->nbins=nbins;
	ret->binwidth=binwidth;

	return ret;
}

void fini_block_histogram(struct block_histogram_t *bh)
{
	if(bh)
	{
		if(bh->means)
			free(bh->means);

		if(bh->variances)
			free(bh->variances);

		if(bh->ns)
			free(bh->ns);

		if(bh->blocks)
		{
			for(int i=0;i<bh->nbins;i++)
				if(bh->blocks[i])
					free(bh->blocks[i]);

			free(bh->blocks);
		}
	
		if(bh->iblock)
			free(bh->iblock);
	}
}

int block_histogram_get_bin(struct block_histogram_t *bh,double t)
{
	int index=t/bh->binwidth;
	
	if(index<0)
		return 0;
	
	if(index>=bh->nbins)
		return bh->nbins-1;
	
	return index;
}

void block_histogram_add_entry(struct block_histogram_t *bh,int index,double value)
{
	assert(index>=0);
	assert(index<bh->nbins);
	assert(bh->iblock[index]>=0);
	assert(bh->iblock[index]<bh->blocksize);

	bh->blocks[index][bh->iblock[index]]=value;
	bh->iblock[index]++;

	if(bh->iblock[index]==bh->blocksize)
	{
		int c;
		double mean,delta1,delta2;

		bh->iblock[index]=0;

		mean=0.0f;

		for(c=0;c<bh->blocksize;c++)
			mean+=bh->blocks[index][c];

		mean/=bh->blocksize;

		/*
			Online algorithm for online mean and variance, taken from https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
		*/


		bh->ns[index]++;

		delta1=mean-bh->means[index];
		bh->means[index]+=delta1/bh->ns[index];
		delta2=mean-bh->means[index];

		bh->variances[index]+=delta1*delta2;
	}
}

void block_histogram_add_time_sample(struct block_histogram_t *bh,double t,double value)
{
	int index=block_histogram_get_bin(bh,t);
	
	for(int c=0;c<bh->nbins;c++)
		block_histogram_add_entry(bh,c,(c==index)?(value):(0.0f));
}

double block_histogram_get_mean(struct block_histogram_t *bh,int index)
{	
	return bh->means[index];
}

double block_histogram_get_variance(struct block_histogram_t *bh,int index)
{
	assert(bh->blocksize>=1);

	return bh->variances[index]/(bh->ns[index])/(bh->ns[index]-1);
}

struct super_sampler_t *init_super_sampler(int *blocksizes,int nbhs,int nbins,double binwidth)
{
	assert(nbhs>0);
	
	struct super_sampler_t *ret=malloc(sizeof(struct super_sampler_t));

	if(!ret)
		return NULL;

	ret->nbhs=nbhs;

	if(!(ret->bhs=malloc(sizeof(struct block_histogram_t)*nbhs)))
	{
		if(ret)
			free(ret);
		
		return NULL;
	}
	
	for(int c=0;c<nbhs;c++)
	{
		assert(blocksizes[c]>0);
		
		ret->bhs[c]=init_block_histogram(nbins,binwidth,blocksizes[c]);
	}

	ret->sh=init_simple_histogram(nbins,binwidth);

	return ret;
}

void fini_super_sampler(struct super_sampler_t *sst)
{
	if(sst)
	{
		if(sst->bhs)
		{
			for(int c=0;c<sst->nbhs;c++)
				fini_block_histogram(sst->bhs[c]);

			fini_simple_histogram(sst->sh);

			free(sst->bhs);
		}

		free(sst);
	}
}

void super_sampler_add_time_sample(struct super_sampler_t *sst,double t,double value)
{
	simple_histogram_add_time_sample(sst->sh,t,value);

	for(int c=0;c<sst->nbhs;c++)
		block_histogram_add_time_sample(sst->bhs[c],t,value);
}

double super_sampler_get_tau(struct super_sampler_t *sst,int bszindex,int binindex)
{
	double ratio,sigma2;

	sigma2=sst->sh->m2s[binindex]/(sst->sh->ns[binindex]-1);
	ratio=sst->sh->ns[binindex]*block_histogram_get_variance(sst->bhs[bszindex],binindex)/sigma2;

	return (ratio-1.0f)/2.0f;
}
