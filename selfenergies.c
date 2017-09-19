#include "selfenergies.h"
#include <gsl/gsl_sf.h>

/*
	We use the Cubature library for numerical integration,
	see http://ab-initio.mit.edu/wiki/index.php/Cubature
*/

#include "cubature/cubature.h"

/*
        Nice progress bar
*/

#include "libprogressbar/progressbar.h"

/*
    GCC does not define M_PI in C11 mode
*/

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

/*
	Tuning paramaters for numerical integration.
*/

size_t maxEval=1024*1024;
double relError=1e-5;

/*
	Physical parameters, as in R. Schmidt and M. Lemeshko, Phys. Rev. Lett. 114, 203001 (2015).
*/

double u0=218;
double u1=218.0f/1.75f;

double r0=1.5;
double r1=1.5;

double abb=3.3;

/*
        Epsilon, also defined in spectral.c

        This regulates the broadness of the peaks in the spectral function,
        they will be tending to a delta function when epsilon tends to zero.
        It is good then to keep epsilon reasonably small but not too small,
        so that the peaks are still visibile.

        Also smaller values give rise to numerical instabilities for higher omegas
        in the second order self-energies.
*/

double epsilon=1e-2;

/*
	Dawson's integral, from: http://www.ebyte.it/library/codesnippets/DawsonIntegralApproximations.html
        (maximum relative error: 10 ppm)
*/

double DawsonF(double x)
{
	double y,p,q;

	y = x*x;

	p = 1.0 + y*(0.1049934947 + y*(0.0424060604
                + y*(0.0072644182 + y*(0.0005064034
                + y*(0.0001789971)))));

	q = 1.0 + y*(0.7715471019 + y*(0.2909738639
                + y*(0.0694555761 + y*(0.0140005442
                + y*(0.0008327945 + 2*0.0001789971*y)))));

	return x*(p/q);
}

/*
	Wigner's 3j symbol from the GNU Scientific Library and the F function
	we get at every vertex in terms of the 3j symbol
*/

double wigner3j(int l1,int l2,int l3,int m1,int m2,int m3)
{
	return gsl_sf_coupling_3j(2*l1,2*l2,2*l3,2*m1,2*m2,2*m3);
}

double F(int l1,int l2,int l3,int m1,int m2,int m3)
{
	return sqrt(((2*l1+1)*(2*l2+1)*(2*l3+1))/(4*M_PI))*wigner3j(l1,l2,l3,0,0,0)*wigner3j(l1,l2,l3,m1,m2,m3);
}

/*
	Dispersion relations for free particles and Bogoliubov quasiparticles
*/

double ek(double k)
{
	return k*k/2.0f;
}

double omegak(double k,double n)
{
	return sqrt(ek(k)*(ek(k)+8.0f*M_PI*abb*n));
}

/*
	Potentials
*/

double U0(double n,double k)
{
	double a,b,c;

	a=sqrt((8.0f*n*k*k*ek(k))/(omegak(k,n)));
	b=exp(-0.5f*k*k*r0*r0);
	c=4*M_PI*pow(r0,-3.0f);

	return u0*a*b/c;
}

double U1(double n,double k)
{
	double a,b,c;

	a=sqrt((8.0f*n*k*k*ek(k))/(3.0f*omegak(k,n)));
	b=r1*(-sqrt(2.0f)*k*r1+2.0f*(1+k*k*r1*r1)*DawsonF(k*r1/sqrt(2.0f)));
	c=4*pow(M_PI,3.0f/2.0f)*k*k;
	
	return u1*a*b/c;
}

double U(short lambda,double n,double k)
{
	switch(lambda)
	{
		case 0:
		return U0(n,k);

		case 1:
		return U1(n,k);
	
		default:
		fprintf(stderr,"U() called with wrong argument lambda=%d\n",lambda);
		exit(0);
	}
}

/*
	Self-energies, first order
*/

struct sigma_ctx_t
{
	double omega,n,cg2;

	short L,lambda,j;
};

int fSigma1(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
	/* The context from which we read the global variables */
	
	struct sigma_ctx_t *ctx=(struct sigma_ctx_t *)(fdata);

	/* Global variables */

	short lambda,j;
	double omega,n,cg2;

	/* The integration variables */

	double t=x[0];
	double k=t/(1.0f-t);

	/* Auxiliary variables */
	
	double complex a,b,c,z;

	/* We load the global variables */

	omega=ctx->omega;
	n=ctx->n;
	cg2=ctx->cg2;

	lambda=ctx->lambda;
	j=ctx->j;

	/* Finally, we compute the result and return it */

	a=(2.0f*lambda+1)/(4.0f*M_PI);
	b=pow(U(lambda,n,k),2.0f)*cg2;
	c=j*(j+1)-(omega+I*epsilon)+omegak(k,n);

	z=(-1.0f)*a*b/c*pow(1-t,-2.0f);
	
	fval[0]=creal(z);
	fval[1]=cimag(z);

	return 0;
}

double complex Sigma1(short L,double omega,double n)
{
	double xmin[1]={0.0f};
	double xmax[1]={1.0f};
	double res[2],err[2];

	double complex total;

	struct sigma_ctx_t ctx;

	ctx.L=L;
	ctx.n=n;
	ctx.omega=omega;

	total=0.0f;

	/*
		The selection rules from L=0 to L=3, along with the Clebsch-Gordan weights
		are calculated in Mathematica and hardcoded here... It could be improved...
	*/

	switch(L)
	{
		case 0:
		
		ctx.lambda=0;
		ctx.j=0;
		ctx.cg2=1.0f;

		hcubature(2,fSigma1,&ctx,1,xmin,xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);
		total+=res[0]+I*res[1];
		
		ctx.lambda=1;
		ctx.j=1;
		ctx.cg2=1.0f;

		hcubature(2,fSigma1,&ctx,1,xmin,xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);
		total+=res[0]+I*res[1];

		break;

		case 1:
		
		ctx.lambda=0;
		ctx.j=1;
		ctx.cg2=1.0f;

		hcubature(2,fSigma1,&ctx,1,xmin,xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);
		total+=res[0]+I*res[1];

		ctx.lambda=1;
		ctx.j=0;
		ctx.cg2=1.0f/3.0f;

		hcubature(2,fSigma1,&ctx,1,xmin,xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);
		total+=res[0]+I*res[1];

		ctx.lambda=1;
		ctx.j=2;
		ctx.cg2=2.0f/3.0f;

		hcubature(2,fSigma1,&ctx,1,xmin,xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);
		total+=res[0]+I*res[1];

		break;

		case 2:
		
		ctx.lambda=0;
		ctx.j=2;
		ctx.cg2=1.0f;

		hcubature(2,fSigma1,&ctx,1,xmin,xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);
		total+=res[0]+I*res[1];

		ctx.lambda=1;
		ctx.j=1;
		ctx.cg2=2.0f/5.0f;

		hcubature(2,fSigma1,&ctx,1,xmin,xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);
		total+=res[0]+I*res[1];

		ctx.lambda=1;
		ctx.j=3;
		ctx.cg2=3.0f/5.0f;

		hcubature(2,fSigma1,&ctx,1,xmin,xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);
		total+=res[0]+I*res[1];

		break;

		case 3:
		
		ctx.lambda=0;
		ctx.j=3;
		ctx.cg2=1.0f;

		hcubature(2,fSigma1,&ctx,1,xmin,xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);
		total+=res[0]+I*res[1];

		ctx.lambda=1;
		ctx.j=2;
		ctx.cg2=3.0f/7.0f;

		hcubature(2,fSigma1,&ctx,1,xmin,xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);
		total+=res[0]+I*res[1];

		ctx.lambda=1;
		ctx.j=4;
		ctx.cg2=4.0f/7.0f;

		hcubature(2,fSigma1,&ctx,1,xmin,xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);
		total+=res[0]+I*res[1];

		break;
		
		default:

		fprintf(stderr,"Not implemented! Please, implement me... :)\n");
		exit(0);
	}

	return total;
}

/*
	Second order, things are getting trickier!
*/

struct sigma2_ctx_t
{
	double omega,n;

	short L,l1,l2,l3,l4,l5;
};

int fSigma2A(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
	/* The context from which we read the global variables */

	struct sigma2_ctx_t *ctx=(struct sigma2_ctx_t *)(fdata);

	/* Global variables */

	short l1,l2,l3,l4,l5;
	double omega,n;

	/* The integration variables */

	double t1=x[0];
	double k1=t1/(1.0f-t1);

	double t2=x[1];
	double k2=t2/(1.0f-t2);

	/* Auxiliary variables */
	
	double complex a,b,c,d,f,z,El2,El4,El5;

	/* We load the global variables */

	omega=ctx->omega;
	n=ctx->n;

	l1=ctx->l1;
	l2=ctx->l2;
	l3=ctx->l3;
	l4=ctx->l4;
	l5=ctx->l5;

	/* Finally, we compute the result and return it */

	El2=l2*(l2+1);
	El4=l4*(l4+1);
	El5=l5*(l5+1);

	a=pow(U(l1,n,k1),2.0f);
	b=pow(U(l3,n,k2),2.0f);
	c=(omega+I*epsilon)-omegak(k1,n)-El2;
	d=(omega+I*epsilon)-omegak(k2,n)-El5;
	f=(-omega-I*epsilon)+El4+omegak(k1,n)+omegak(k2,n);

	z=(-1.0f)*a*b/(c*d*f)*pow(1-t1,-2.0f)*pow(1-t2,-2.0f);

	fval[0]=creal(z);
	fval[1]=cimag(z);

	return 0;
}

int fSigma2B(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
	/* The context from which we read the global variables */
	
	struct sigma2_ctx_t *ctx=(struct sigma2_ctx_t *)(fdata);

	/* Global variables */

	short l1,l2,l3,l4,l5;
	double omega,n;

	/* The integration variables */

	double t1=x[0];
	double k1=t1/(1.0f-t1);

	double t2=x[1];
	double k2=t2/(1.0f-t2);

	/* Auxiliary variables */
	
	double complex a,b,c,d,f,z,El2,El4,El5;

	/* We load the global variables */

	omega=ctx->omega;
	n=ctx->n;

	l1=ctx->l1;
	l2=ctx->l2;
	l3=ctx->l3;
	l4=ctx->l4;
	l5=ctx->l5;

	/* Finally, we compute the result and return it */

	El2=l2*(l2+1);
	El4=l4*(l4+1);
	El5=l5*(l5+1);

	a=pow(U(l1,n,k1),2.0f);
	b=pow(U(l3,n,k2),2.0f);
	c=(omega+I*epsilon)-omegak(k1,n)-El2;
	d=(omega+I*epsilon)-omegak(k1,n)-El5;
	f=(-omega-I*epsilon)+El4+omegak(k1,n)+omegak(k2,n);

	z=(-1.0f)*a*b/(c*d*f)*pow(1-t1,-2.0f)*pow(1-t2,-2.0f);

	fval[0]=creal(z);
	fval[1]=cimag(z);

	return 0;
}

/*
	Selection rules for 2nd order diagrams
*/

double complex Sigma2_from_selection_rule(struct selection_rule_t rule,char type,double omega,double n)
{
	double xmin[2]={0.0f,0.0f};
	double xmax[2]={1.0f,1.0f};
	double res[2],err[2];

	double complex total;

	struct sigma2_ctx_t ctx;

	ctx.L=rule.ls[0];
	ctx.n=n;
	ctx.omega=omega;

	ctx.l1=rule.ls[1];
	ctx.l2=rule.ls[2];
	ctx.l3=rule.ls[3];
	ctx.l4=rule.ls[4];
	ctx.l5=rule.ls[5];

	total=0.0f;

	switch(type)
	{
		case 'A':
		
		hcubature(2,fSigma2A,&ctx,2,xmin,xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);
		total+=res[0]+I*res[1];

		break;
		
		case 'B':

		hcubature(2,fSigma2B,&ctx,2,xmin,xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);
		total+=res[0]+I*res[1];

		break;
		
		default:
		
		fprintf(stderr,"Wrong 2nd-order diagram type!\n");
		exit(0);
	}

	return rule.weight*total;
}

struct selection_rule_t *arules[SIGMA_HIGHEST_L_ORDER+1],*brules[SIGMA_HIGHEST_L_ORDER+1];
short iarules[SIGMA_HIGHEST_L_ORDER+1],ibrules[SIGMA_HIGHEST_L_ORDER+1];

double complex Sigma2(short L,double omega,double n)
{
	double complex res=0.0f;
	int c;

	if(L>SIGMA_HIGHEST_L_ORDER)
	{
		fprintf(stderr,"Sigma order too high! (L=%d, Lmax=%d)\n",L,SIGMA_HIGHEST_L_ORDER);
		exit(0);
	}

	for(c=0;c<iarules[L];c++)
		res+=Sigma2_from_selection_rule(arules[L][c],'A',omega,n);

	for(c=0;c<ibrules[L];c++)
		res+=Sigma2_from_selection_rule(brules[L][c],'B',omega,n);
	
	return res;
}

double complex Sigma2A(short L,double omega,double n)
{
	double complex res=0.0f;
	int c;

	if(L>SIGMA_HIGHEST_L_ORDER)
	{
		fprintf(stderr,"Sigma order too high! (L=%d, Lmax=%d)\n",L,SIGMA_HIGHEST_L_ORDER);
		exit(0);
	}

	for(c=0;c<iarules[L];c++)
		res+=Sigma2_from_selection_rule(arules[L][c],'A',omega,n);
	
	return res;
}

double complex Sigma2B(short L,double omega,double n)
{
	double complex res=0.0f;
	int c;

	if(L>SIGMA_HIGHEST_L_ORDER)
	{
		fprintf(stderr,"Sigma order too high! (L=%d, Lmax=%d)\n",L,SIGMA_HIGHEST_L_ORDER);
		exit(0);
	}

	for(c=0;c<ibrules[L];c++)
		res+=Sigma2_from_selection_rule(brules[L][c],'B',omega,n);

	return res;
}

/*
	Calculates all the allowed internal momenta on a 2nd order diagram,
	given the initial and final state (L,M).

	A list of selection rules is created, to use at a later time.
*/

#define MAXL		(8)
#define MAX_RULES	(1024)

struct selection_rule_t *gen_selection_rules(short L,short M,char type,short *nr,progressbar *progress)
{
	struct selection_rule_t *rules;
	short irules;

	short ls[6];
	short ms[6];

	irules=0;
	rules=malloc(sizeof(struct selection_rule_t)*MAX_RULES);

#define LOOP(i)		for(ls[i]=0;ls[i]<=MAXL;ls[i]++) for(ms[i]=-ls[i];ms[i]<=ls[i];ms[i]++)
#define REDUCED_LOOP(i)	for(ls[i]=0;ls[i]<=1;ls[i]++) for(ms[i]=-ls[i];ms[i]<=ls[i];ms[i]++)

	/*
		The interaction potential is defined only on the l=0 and l=1 channels. As a consequence
		if a quantum number l_i appears in a |U_{l_i} (k)|^2 term, we can restrict the sum to
		just l=0 and l=1.
	
		Hence the different summation in REDUCED_LOOP, as opposed to LOOP
	*/

	ls[0]=L;
	ms[0]=M;

	LOOP(2)
	LOOP(4)
	LOOP(5)
	REDUCED_LOOP(1)
	REDUCED_LOOP(3)
	{
		short L=ls[0];
		short M=ms[0];

		double weight;

		switch(type)
		{
			case 'A':
			
			weight=F(L,ls[1],ls[2],M,-ms[1],-ms[2])*
			       F(ls[2],ls[3],ls[4],ms[2],-ms[3],-ms[4])*
			       F(ls[4],ls[1],ls[5],ms[4],ms[1],-ms[5])*
			       F(ls[5],ls[3],L,ms[5],ms[3],-M)*
			       pow(-1.0f,ms[1]+ms[2]+ms[3]+ms[4]+ms[5]);

			break;
			
			case 'B':

			weight=F(L,ls[1],ls[2],M,-ms[1],-ms[2])*
			       F(ls[2],ls[3],ls[4],ms[2],-ms[3],-ms[4])*
			       F(ls[3],ls[4],ls[5],ms[3],ms[4],-ms[5])*
			       F(ls[1],ls[5],L,ms[1],ms[5],-M)*
			       pow(-1.0f,ms[1]+ms[2]+ms[3]+ms[4]+ms[5]);

			break;
			
			default:
			
			fprintf(stderr,"Wrong 2nd-order diagram type!\n");
			exit(0);
		}

		if(weight!=0.0f)
		{
			if(irules==MAX_RULES)
			{
				fprintf(stderr,"Too many selection rules!\n");
				exit(0);
			}

			memcpy(rules[irules].ls,ls,sizeof(short)*6);
			memcpy(rules[irules].ms,ms,sizeof(short)*6);
			rules[irules].weight=weight;
	
			irules++;
		}

		if((ls[4]==0)&&(ms[4]==0)&&(ls[5]==0)&&(ms[5]==0)&&(ls[1]==0)&&(ms[1]==0)&&(ls[3]==0)&&(ms[3]==0))
		{
#if defined(_OPENMP)
#pragma omp critical
#endif

			{
				progressbar_inc(progress);
			}
		}
	}

	*nr=irules;
	
	return rules;
}

int rule_compar(const void *x,const void *y)
{
 	const struct selection_rule_t *a,*b;
	int c;

        a=x;
        b=y;

	for(c=1;c<=5;c++)
	{
		if(a->ls[c]!=b->ls[c])
			return a->ls[c]-b->ls[c];
	}

        return 0;
}

void init_selection_rules(void)
{
	int c,d;

	progressbar *progress;

	progress=progressbar_new("Precalculating selection rules...",2*(1+MAXL)*(1+MAXL)*(SIGMA_HIGHEST_L_ORDER+1));

#if defined(_OPENMP)
#pragma omp parallel for
#endif

	for(c=0;c<=SIGMA_HIGHEST_L_ORDER;c++)
	{	
		/*
			We generate the rules
		*/

		arules[c]=gen_selection_rules(c,0,'A',&iarules[c],progress);
		brules[c]=gen_selection_rules(c,0,'B',&ibrules[c],progress);

		/*
			We sort the rules with respect to the l quantum numbers
		*/
	
		qsort(arules[c],iarules[c],sizeof(struct selection_rule_t),rule_compar);
		qsort(brules[c],ibrules[c],sizeof(struct selection_rule_t),rule_compar);
		
		/*
			We remove duplicate rules (i.e. rules with the same l numbers, differing
			just as far as the m numbers are concerned), creating a new rule with the sum of
			the weights as its new weight.
		
			This is equivalent to carrying out the sum over m_1 ... m_5 using the orthogonality relations.
		*/
	
		for(d=0;d<(iarules[c]-1);d++)
		{
			if(rule_compar(&arules[c][d],&arules[c][d+1])==0)
			{
				arules[c][d+1].weight+=arules[c][d].weight;
				arules[c][d].weight=0.0f;
			}
		}

		for(d=0;d<(iarules[c]-1);d++)
		{
			if(fabs(arules[c][d].weight)<10e-10)
			{
				memmove(&arules[c][d],&arules[c][d+1],sizeof(struct selection_rule_t)*(iarules[c]-d-1));
				iarules[c]--;
				d--;
			}
		}

		for(d=0;d<(ibrules[c]-1);d++)
		{
			if(rule_compar(&brules[c][d],&brules[c][d+1])==0)
			{
				brules[c][d+1].weight+=brules[c][d].weight;
				brules[c][d].weight=0.0f;
			}
		}

		for(d=0;d<(ibrules[c]-1);d++)
		{
			if(fabs(brules[c][d].weight)<10e-10)
			{
				memmove(&brules[c][d],&brules[c][d+1],sizeof(struct selection_rule_t)*(ibrules[c]-d-1));
				ibrules[c]--;
				d--;
			}
		}
	}

	progressbar_finish(progress);
}

void fini_selection_rules(void)
{
	int c;

	for(c=0;c<=SIGMA_HIGHEST_L_ORDER;c++)
	{
		if(arules[c])
			free(arules[c]);

		if(brules[c])
			free(brules[c]);
	}
}

/*
	The chi function used in DiagMC, representing the weight of a phonon arc including the vertex part.
*/

struct chi_ctx_t
{
	double deltat,n;

	short lambda;
};

int fchi(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
	/* The context from which we read the global variables */
	
	struct chi_ctx_t *ctx=(struct chi_ctx_t *)(fdata);

	/* Global variables */

	short lambda;
	double deltat,n;

	/* The integration variables */

	double t=x[0];
	double k=t/(1.0f-t);

	/* Auxiliary variables */

	double a,b,z;

	/* We load the global variables */

	deltat=ctx->deltat;
	n=ctx->n;

	lambda=ctx->lambda;

	/* Finally, we compute the result and return it */

	a=pow(U(lambda,n,k),2.0f);
	b=exp(-deltat*omegak(k,n));

	z=a*b*pow(1-t,-2.0f);
	
	fval[0]=z;

	return 0;
}

double chi(short lambda,double deltat,double n)
{
	double xmin[1]={0.0f};
	double xmax[1]={1.0f};
	double res,err;

	struct chi_ctx_t ctx;

	ctx.lambda=lambda;
	ctx.deltat=deltat;
	ctx.n=n;

	hcubature(1,fchi,&ctx,1,xmin,xmax,maxEval,0,relError,ERROR_INDIVIDUAL,&res,&err);

	return res;
}
