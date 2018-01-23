#ifndef __HISTOGRAMS_H__
#define __HISTOGRAMS_H__

struct simple_histogram_t
{
	double *means;
	double *m2s;
	int *ns;

	int nbins;
	double binwidth;
};

struct simple_histogram_t *init_simple_histogram(int nbins,double binwidth);
void fini_simple_histogram(struct simple_histogram_t *sh);

int simple_histogram_get_bin(struct simple_histogram_t *sh,double t);
void simple_histogram_add_entry(struct simple_histogram_t *sh,int index,double value);
void simple_histogram_add_time_sample(struct simple_histogram_t *sh,double t,double value);
double simple_histogram_get_mean(struct simple_histogram_t *sh,int index);
double simple_histogram_get_variance(struct simple_histogram_t *sh,int index);

struct block_histogram_t
{
	double *means;
	double *variances;
	int *ns;

	double **blocks;
	int *iblock;

	int blocksize;

	int nbins;
	double binwidth;
};

struct block_histogram_t *init_block_histogram(int nbins,double binwidth,int blocksize);
void fini_block_histogram(struct block_histogram_t *bh);

int block_histogram_get_bin(struct block_histogram_t *bh,double t);
void block_histogram_add_entry(struct block_histogram_t *bh,int index,double value);
void block_histogram_add_time_sample(struct block_histogram_t *bh,double t,double value);
double block_histogram_get_mean(struct block_histogram_t *bh,int index);
double block_histogram_get_variance(struct block_histogram_t *bh,int index);

struct super_sampler_t
{
	struct simple_histogram_t *sh;
	struct block_histogram_t **bhs;
	int nbhs;
};

struct super_sampler_t *init_super_sampler(int *blocksizes,int nbhs,int nbins,double binwidth);
void fini_super_sampler(struct super_sampler_t *sst);

void super_sampler_add_time_sample(struct super_sampler_t *sst,double t,double value);
double super_sampler_get_tau(struct super_sampler_t *sst,int bszindex,int binindex);

#endif //__HISTOGRAMS_H__