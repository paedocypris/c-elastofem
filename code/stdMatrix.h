#ifndef STDMATRIX_H_INCLUDED
#define STDMATRIX_H_INCLUDED

#include <stdint.h>

typedef struct stdmatrix StdMatrix;

int stdmatrix_create(StdMatrix **m, uint32_t ni, uint32_t nj);
int stdmatrix_createzero(StdMatrix **m, uint32_t ni, uint32_t nj);
int stdmatrix_destroy(StdMatrix *m);
int stdmatrix_multiply(const StdMatrix *m, const StdMatrix *n, StdMatrix *r);
int stdmatrix_multiplymTn(const StdMatrix *m, const StdMatrix *n, StdMatrix *r);
double stdmatrix_getelem(const StdMatrix *m, uint32_t i, uint32_t j);
int stdmatrix_setelem(StdMatrix *m, uint32_t i, uint32_t j, double val);

#endif
