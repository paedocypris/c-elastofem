#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include "stdMatrix.h"

struct stdmatrix {
  uint32_t ni;
  uint32_t nj;
  double *val;
};

inline void sval(StdMatrix *m, uint32_t i, uint32_t j, double val);
inline double gval(const StdMatrix *m, uint32_t i, uint32_t j);
inline uint32_t gi(const StdMatrix *m, uint32_t i, uint32_t j);

int stdmatrix_create(StdMatrix **m, uint32_t ni, uint32_t nj)
{
  StdMatrix *newStdMatrix = malloc(sizeof(StdMatrix));
  newStdMatrix->ni = ni;
  newStdMatrix->nj = nj;
  newStdMatrix->val = malloc(ni*nj*sizeof(double));
  *m = newStdMatrix;
  return 0;
}

int stdmatrix_createzero(StdMatrix **m, uint32_t ni, uint32_t nj)
{
  StdMatrix *newStdMatrix = malloc(sizeof(StdMatrix));
  newStdMatrix->ni = ni;
  newStdMatrix->nj = nj;
  newStdMatrix->val = calloc(ni*nj, sizeof(double));
  *m = newStdMatrix;
  return 0;
}

int stdmatrix_destroy(StdMatrix *m)
{
  free(m->val);
  free(m);
  return 0;
}

int stdmatrix_multiply(const StdMatrix *m, const StdMatrix *n, StdMatrix *r)
{
  if (m == NULL || n == NULL || r == NULL)
  {
    fprintf(stderr, "stdmatrix_multiply: Matriz não definida.\n");
    return 1;
  }
  
  if (m->nj != n->ni || m->ni != r->ni || n->nj != r->nj)
  {
    fprintf(stderr, "stdmatrix_multiply: Multiplicação de matrizes com dimensões incorretas, m: %" PRIu32 " %" PRIu32 ", n: %" PRIu32 " %" PRIu32 ", r: %" PRIu32 " %" PRIu32 ".\n", m->ni, m->nj, n->ni, n->nj, r->ni, r->nj);
    return 1;
  }
  
  uint32_t i, j, k;
  for (i = 0; i < r->ni; i++)
  {
    for (j = 0; j < r->nj; j++)
    {
      double sum = 0;
      for (k = 0; k < m->nj; k++)
      {
        sum += gval(m, i, k) * gval(n, k, j);
      }
      sval(r, i, j, sum);
    }
  }
  
  return 0;
}

int stdmatrix_multiplymTn(const StdMatrix *m, const StdMatrix *n, StdMatrix *r)
{
  if (m == NULL || n == NULL || r == NULL)
  {
    fprintf(stderr, "stdmatrix_multiply: Matriz não definida.\n");
    return 1;
  }
  
  if (m->ni != n->ni || m->nj != r->ni || n->nj != r->nj)
  {
    fprintf(stderr, "stdmatrix_multiply: Multiplicação de matrizes com dimensões incorretas, m: %" PRIu32 " %" PRIu32 ", n: %" PRIu32 " %" PRIu32 ", r: %" PRIu32 " %" PRIu32 ".\n", m->ni, m->nj, n->ni, n->nj, r->ni, r->nj);
    return 1;
  }
  
  uint32_t i, j, k;
  for (i = 0; i < r->ni; i++)
  {
    for (j = 0; j < r->nj; j++)
    {
      double sum = 0;
      for (k = 0; k < m->ni; k++)
      {
        sum += gval(m, k, i) * gval(n, k, j);
      }
      sval(r, i, j, sum);
    }
  }
  
  return 0;
}

double stdmatrix_getelem(const StdMatrix *m, uint32_t i, uint32_t j)
{
  return gval(m, i, j);
}

int stdmatrix_setelem(StdMatrix *m, uint32_t i, uint32_t j, double val)
{
  if (i >= m->ni || j >= m->nj)
  {
    fprintf(stderr, "stdmatrix_setelem: Elemento fora da matriz.\n");
    return 1;
  }
  
  sval(m, i, j, val);
  return 0;
}

inline void sval(StdMatrix *m, uint32_t i, uint32_t j, double val)
{
  m->val[gi(m, i, j)] = val;
}

inline double gval(const StdMatrix *m, uint32_t i, uint32_t j)
{
  return m->val[gi(m, i, j)];
}

inline uint32_t gi(const StdMatrix *m, uint32_t i, uint32_t j)
{
  return i*m->ni + j;
}