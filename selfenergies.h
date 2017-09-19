#ifndef __SELF_ENERGIES_H__
#define __SELF_ENERGIES_H__

#include <math.h>
#include <complex.h>

double U(short lambda,double n,double k);
double omegak(double k,double n);

double complex Sigma1(short L,double omega,double n);
double complex Sigma2(short L,double omega,double n);

double complex Sigma2A(short L,double omega,double n);
double complex Sigma2B(short L,double omega,double n);

void init_selection_rules(void);
void fini_selection_rules(void);

struct selection_rule_t
{
	short ls[6];
	short ms[6];
	
	double weight;
};

#define SIGMA_HIGHEST_L_ORDER	(3)

extern struct selection_rule_t *arules[SIGMA_HIGHEST_L_ORDER+1],*brules[SIGMA_HIGHEST_L_ORDER+1];
extern short iarules[SIGMA_HIGHEST_L_ORDER+1],ibrules[SIGMA_HIGHEST_L_ORDER+1];

double chi(short lambda,double deltat,double n);

#endif //__SELF_ENERGIES_H__

