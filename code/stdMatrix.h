#ifndef STDMATRIX_H_INCLUDED
#define STDMATRIX_H_INCLUDED

#include <stdint.h>

typedef struct stdmatrix StdMatrix;


int stdmatrix_create(StdMatrix **m, uint32_t ni, uint32_t nj);
int stdmatrix_destroy(StdMatrix *m);
int stdmatrix_copy(const StdMatrix *m, StdMatrix *r);

int stdmatrix_multiply(const StdMatrix *m, const StdMatrix *n, StdMatrix *r);
int stdmatrix_multiplymTn(const StdMatrix *m, const StdMatrix *n, StdMatrix *r);
int stdmatrix_multiplymnT(const StdMatrix *m, const StdMatrix *n, StdMatrix *r);
double stdmatrix_add(const StdMatrix *m, const StdMatrix *n, StdMatrix *r);
double stdmatrix_subtract(const StdMatrix *m, const StdMatrix *n, StdMatrix *r);
int stdmatrix_multiplyScalar(StdMatrix *m, double scalar);
int stdmatrix_addmultiplyscalar(const StdMatrix *m, const StdMatrix *n, double scalar, StdMatrix *r);
int stdmatrix_subtractmultiplyscalar(double scalar, const StdMatrix *m, const StdMatrix *n, StdMatrix *r);
double stdmatrix_norm(const StdMatrix *m);
double stdmatrix_voigtStressNormSqr(const StdMatrix *m);
double stdmatrix_voigtTrace(const StdMatrix *m);
int stdmatrix_invert(const StdMatrix *m, StdMatrix *mInv);

double stdmatrix_getelem(const StdMatrix *m, uint32_t i, uint32_t j);
int stdmatrix_setelem(StdMatrix *m, uint32_t i, uint32_t j, double val);

int stdmatrix_print(const StdMatrix *m);

uint32_t stdmatrix_ni(const StdMatrix *m);
uint32_t stdmatrix_nj(const StdMatrix *m);
double *stdmatrix_getptr(const StdMatrix *m);

#endif
