#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include "helper.h"
#include "stdMatrix.h"

struct stdmatrix {
  uint32_t ni;
  uint32_t nj;
  double *val;
};

inline void sval(StdMatrix *m, uint32_t i, uint32_t j, double val);
inline double gval(const StdMatrix *m, uint32_t i, uint32_t j);
inline uint32_t gi(const StdMatrix *m, uint32_t i, uint32_t j);

static void luDec(const StdMatrix *m, StdMatrix *lu);
static void luSolv(const StdMatrix *lu, StdMatrix *b);

int stdmatrix_create(StdMatrix **m, uint32_t ni, uint32_t nj)
{
  StdMatrix *newStdMatrix = sMalloc(sizeof(StdMatrix));
  newStdMatrix->ni = ni;
  newStdMatrix->nj = nj;
  newStdMatrix->val = sCalloc(ni*nj, sizeof(double));
  *m = newStdMatrix;
  return 0;
}

int stdmatrix_destroy(StdMatrix *m)
{
  free(m->val);
  free(m);
  return 0;
}

int stdmatrix_copy(const StdMatrix *m, StdMatrix *r)
{
  if (m == NULL || r == NULL)
  {
    fprintf(stderr, "stdmatrix_copy: Matrix não definida.\n");
    return 1;
  }
  if (m->ni != r->ni || m->nj != r->nj)
  {
    fprintf(stderr, "stdmatrix_copy: Matrizes de tamanhos diferentes, m: %" PRIu32 " %" PRIu32 ", r: %" PRIu32 " %" PRIu32 ".\n", m->ni, m->nj, r->ni, r->nj);
    return 1;
  }

  for (uint32_t i = 0; i < m->ni; i++)
  {
    for (uint32_t j = 0; j < m->nj; j++)
    {
      sval(r, i, j, gval(m, i, j));
    }
  }
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
    fprintf(stderr, "stdmatrix_multiplymTn: Matriz não definida.\n");
    exit(-1);
  }

  if (m->ni != n->ni || m->nj != r->ni || n->nj != r->nj)
  {
    fprintf(stderr, "stdmatrix_multiplymTn: Multiplicação de matrizes com dimensões incorretas, m: %" PRIu32 " %" PRIu32 ", n: %" PRIu32 " %" PRIu32 ", r: %" PRIu32 " %" PRIu32 ".\n", m->ni, m->nj, n->ni, n->nj, r->ni, r->nj);
    exit(-1);
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

int stdmatrix_multiplymnT(const StdMatrix *m, const StdMatrix *n, StdMatrix *r)
{
  if (m == NULL || n == NULL || r == NULL)
  {
    fprintf(stderr, "stdmatrix_multiplymnT: Matriz não definida.\n");
    exit(-1);
  }

  if (m->nj != n->nj || m->nj != r->ni || n->ni != r->nj)
  {
    fprintf(stderr, "stdmatrix_multiplymnT: Multiplicação de matrizes com dimensões incorretas, m: %" PRIu32 " %" PRIu32 ", n: %" PRIu32 " %" PRIu32 ", r: %" PRIu32 " %" PRIu32 ".\n", m->ni, m->nj, n->ni, n->nj, r->ni, r->nj);
    exit(-1);
  }

  for (uint32_t i = 0; i < r->ni; i++)
  {
    for (uint32_t j = 0; j < r->nj; j++)
    {
      double sum = 0;
      for (uint32_t k = 0; k < m->nj; k++)
      {
        sum += gval(m, i, k) * gval(n, j, k);
      }
      sval(r, i, j, sum);
    }
  }
  return 0;
}

double stdmatrix_add(const StdMatrix *m, const StdMatrix *n, StdMatrix *r)
{
  if (m == NULL || n == NULL || r == NULL)
  {
    fprintf(stderr, "stdmatrix_add: Matriz não definida.\n");
    exit(-1);
  }

  if (m->ni != n->ni || m->nj != n->nj
    || m->ni != r->ni || m->nj != r->nj)
  {
    fprintf(stderr, "stdmatrix_add: Subtração de matrizes com dimensões incorretas, m: %" PRIu32 " %" PRIu32 ", n: %" PRIu32 " %" PRIu32 ", r: %" PRIu32 " %" PRIu32 ".\n", m->ni, m->nj, n->ni, n->nj, r->ni, r->nj);
    exit(-1);
  }

  for (uint32_t i = 0; i < m->ni; i++)
  {
    for (uint32_t j = 0; j < m->nj; j++)
    {
      sval(r, i, j, gval(m, i, j) + gval(n, i, j));
    }
  }
  return 0;
}

double stdmatrix_subtract(const StdMatrix *m, const StdMatrix *n, StdMatrix *r)
{
  if (m == NULL || n == NULL || r == NULL)
  {
    fprintf(stderr, "stdmatrix_subtract: Matriz não definida.\n");
    exit(-1);
  }

  if (m->ni != n->ni || m->nj != n->nj
    || m->ni != r->ni || m->nj != r->nj)
  {
    fprintf(stderr, "stdmatrix_subtract: Subtração de matrizes com dimensões incorretas, m: %" PRIu32 " %" PRIu32 ", n: %" PRIu32 " %" PRIu32 ", r: %" PRIu32 " %" PRIu32 ".\n", m->ni, m->nj, n->ni, n->nj, r->ni, r->nj);
    exit(-1);
  }

  for (uint32_t i = 0; i < m->ni; i++)
  {
    for (uint32_t j = 0; j < m->nj; j++)
    {
      sval(r, i, j, gval(m, i, j) - gval(n, i, j));
    }
  }
  return 0;
}

int stdmatrix_multiplyScalar(StdMatrix *m, double scalar)
{
  if (m == NULL)
  {
    fprintf(stderr, "stdmatrix_multiplyScalar: Matriz não definida.\n");
    exit(-1);
  }

  for (uint32_t i = 0; i < m->ni*m->nj; i++)
  {
    m->val[i] = m->val[i] * scalar;
  }
  return 0;
}

int stdmatrix_addmultiplyscalar(const StdMatrix *m, const StdMatrix *n, double scalar, StdMatrix *r)
{
  if (m == NULL || n == NULL || r == NULL)
  {
    fprintf(stderr, "stdmatrix_addmultiplyscalar: Matriz não definida.\n");
    exit(-1);
  }

  if (m->ni != n->ni || m->nj != n->nj
    || m->ni != r->ni || m->nj != r->nj)
  {
    fprintf(stderr, "stdmatrix_addmultiplyscalar: Subtração de matrizes com dimensões incorretas, m: %" PRIu32 " %" PRIu32 ", n: %" PRIu32 " %" PRIu32 ", r: %" PRIu32 " %" PRIu32 ".\n", m->ni, m->nj, n->ni, n->nj, r->ni, r->nj);
    exit(-1);
  }

  for (uint32_t i = 0; i < m->ni; i++)
  {
    for (uint32_t j = 0; j < m->nj; j++)
    {
      sval(r, i, j, gval(m, i, j) + scalar * gval(n, i, j));
    }
  }
  return 0;
}

int stdmatrix_subtractmultiplyscalar(double scalar, const StdMatrix *m, const StdMatrix *n, StdMatrix *r)
{
  if (m == NULL || n == NULL || r == NULL)
  {
    fprintf(stderr, "stdmatrix_subtractmultiplyscalar: Matriz não definida.\n");
    exit(-1);
  }

  if (m->ni != n->ni || m->nj != n->nj
    || m->ni != r->ni || m->nj != r->nj)
  {
    fprintf(stderr, "stdmatrix_subtractmultiplyscalar: Subtração de matrizes com dimensões incorretas, m: %" PRIu32 " %" PRIu32 ", n: %" PRIu32 " %" PRIu32 ", r: %" PRIu32 " %" PRIu32 ".\n", m->ni, m->nj, n->ni, n->nj, r->ni, r->nj);
    exit(-1);
  }

  for (uint32_t i = 0; i < m->ni; i++)
  {
    for (uint32_t j = 0; j < m->nj; j++)
    {
      sval(r, i, j, scalar * gval(m, i, j) - gval(n, i, j));
    }
  }
  return 0;
}

double stdmatrix_norm(const StdMatrix *m)
{
  if (m == NULL)
  {
    fprintf(stderr, "stdmatrix_norm: Matriz não definida.\n");
    exit(-1);
  }

  double normSquare = 0;
  double val;
  for (uint32_t i = 0; i < m->ni; i++)
  {
    for (uint32_t j = 0; j < m->nj; j++)
    {
      val = gval(m, i, j);
      normSquare += val*val;
    }
  }
  return sqrt(normSquare);
}

double stdmatrix_voigtStressNormSqr(const StdMatrix *m)
{
  if (m == NULL || m->nj != 1)
  {
    fprintf(stderr, "stdmatrix_voigtStressNorm: Matriz não definida.\n");
    exit(-1);
  }

  double norm = 0;
  for (uint32_t i = 0; i < 6; i++)
  {
    norm += m->val[i] * m->val[i];
  }
  return norm;
}

double stdmatrix_voigtTrace(const StdMatrix *m)
{
  if (m == NULL || m->nj != 1)
  {
    fprintf(stderr, "stdmatrix_voigtTrace: Matriz não definida.\n");
    exit(-1);
  }

  return m->val[0] + m->val[1] + m->val[2];
}

int stdmatrix_invert(const StdMatrix *m, StdMatrix *mInv)
{
  /* LU Decomposition */
  if (m == NULL || mInv == NULL)
  {
    fprintf(stderr, "stdmatrix_invert: Matriz não definida.\n");
    exit(-1);
  }

  if (m->ni != mInv->ni || m->ni != m->nj || mInv->ni != mInv->nj)
  {
    fprintf(stderr, "stdmatrix_invert: Matrizes com tamanhos incorretos.\n");
    exit(-1);
  }

  uint32_t n = m->ni;

  StdMatrix *lu;
  stdmatrix_create(&lu, n, n);

  StdMatrix *col;
  stdmatrix_create(&col, n, 1);

  luDec(m, lu);
  for (uint32_t j = 0; j < n; j++)
  {
    for (uint32_t i = 0; i < n; i++)
    {
      sval(col, i, 0, 0);
    }
    sval(col, j, 0, 1);
    luSolv(lu, col);
    for (uint32_t i = 0; i < n; i++)
    {
      sval(mInv, i, j, gval(col, i, 0));
    }
  }

  stdmatrix_destroy(lu);
  stdmatrix_destroy(col);
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

int stdmatrix_print(const StdMatrix *m)
{
  if (m == NULL)
  {
    fprintf(stderr, "stdmatrix_print: Invalid matrix.\n");
    exit(-1);
  }

  printf("%" PRIu32 " %" PRIu32 "\n", m->ni, m->nj);

  for (uint32_t i = 0; i < m->ni; i++)
  {
    for (uint32_t j = 0; j < m->nj; j++)
    {
      printf("%7.3g ", gval(m, i, j));
    }
    printf("\n");
  }
  printf("\n");
  return 0;
}

uint32_t stdmatrix_ni(const StdMatrix *m)
{
  return m->ni;
}

uint32_t stdmatrix_nj(const StdMatrix *m)
{
  return m->nj;
}

double *stdmatrix_getptr(const StdMatrix *m)
{
  return m->val;
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
  return i*m->nj + j;
}

static void luDec(const StdMatrix *m, StdMatrix *lu)
{
  uint32_t n = m->ni;

  for (uint32_t k = 0; k < n; k++)
  {
    for (uint32_t i = k; i < n; i++)
    {
      double sum = 0;
      for (uint32_t p = 0; p < k; p++)
      {
        sum += gval(lu, i, p) * gval(lu, p, k);
      }
      sval(lu, i, k, gval(m, i, k) - sum);
    }

    for (uint32_t j = k + 1; j < n; j++)
    {
      double sum = 0;
      for (uint32_t p = 0; p < k; p++)
      {
        sum += gval(lu, k, p) * gval(lu, p, j);
      }
      sval(lu, k, j, (gval(m, k, j) - sum) / gval(lu, k, k));
    }
  }
}

static void luSolv(const StdMatrix *lu, StdMatrix *b)
{
  uint32_t n = lu->ni;

  double *y = sMalloc(n * sizeof(double));
  for (uint32_t i = 0; i < n; i++)
  {
    double sum = 0;
    for (uint32_t k = 0; k < i; k++)
    {
      sum += gval(lu, i, k) * y[k];
    }
    y[i] = (gval(b, i, 0) - sum) / gval(lu, i, i);
  }

  for (int64_t i = n - 1; i >= 0; --i)
  {
    double sum = 0;
    for (uint32_t k = (uint32_t)i+1; k < n; k++)
    {
      sum += gval(lu, (uint32_t)i, k)*gval(b, k, 0);
    }
    sval(b, (uint32_t)i, 0, y[i] - sum);
  }
  free(y);
}