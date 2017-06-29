#ifndef __COUPLINGS_H__
#define __COUPLINGS_H__

bool diagram_calculate_couplings(struct diagram_t *dgr,int lo,int hi);
void diagram_identify_invalid_couplings(struct diagram_t *dgr,int *first,int *last);

#endif //__COUPLINGS_H__
