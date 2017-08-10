#ifndef __NJGRAF_H__
#define __NJGRAF_H__

#define MANGM	(60)
#define MTRIAD	(12)

void init_njgraf(void);
double do_njgraf(int m,int n,int j1[MANGM],int j2[MTRIAD][3],int j3[MTRIAD][3],int isfree[MANGM]);

#endif
