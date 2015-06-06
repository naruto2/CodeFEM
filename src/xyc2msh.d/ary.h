#ifndef _ESTIVA_ARY_H_
#define _ESTIVA_ARY_H_

extern "C" {
void    estiva_ary1(void** v,long n_1, size_t o);
long    estiva_dim1(void* v);
}

extern void arytovec_nde(nde *N, vector<nde> &Nv);
extern void vectoary_nde(vector<nde> &Nv, nde *N);

#define ary1(b,n_1)        estiva_ary1((void *)&b,n_1,sizeof(*b))
#define dim1(v)            estiva_dim1(v)


#endif
