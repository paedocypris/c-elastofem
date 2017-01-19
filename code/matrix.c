#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <float.h>
#include <math.h>
#include "matrix.h"

#define internal static

typedef struct matrixElem {
  double val;
  uint32_t i;
  uint32_t j;
  struct matrixElem *next;
} MatrixElem;

struct matrix {
  uint32_t ni;
  uint32_t nj;
  struct matrixElem **lines;
};

struct matrixcrs {
  double *sa;
  uint32_t *ija;
};

struct vector {
  double *val;
  uint32_t len;
};

internal inline MatrixElem* createNewNode(uint32_t i, uint32_t j, double val);
internal inline void removeNode(Matrix *m, MatrixElem *prevNode, MatrixElem *curNode);
internal inline void insertNewNode(MatrixElem *headNode, MatrixElem *prevNode, MatrixElem *curNode, MatrixElem *newElem);
internal MatrixElem* findNextNode(const Matrix* m, const MatrixElem* node);

internal int floatIsZero(double x);
internal inline void *sMalloc(size_t memsize);

int matrix_create(Matrix** m, uint32_t ni, uint32_t nj)
{
  Matrix* newMatrix = sMalloc(sizeof(Matrix));
  newMatrix->ni = ni;
  newMatrix->nj = nj;
  newMatrix->lines = sMalloc(ni * sizeof(MatrixElem*));  
  for(uint32_t i = 0; i < ni; i++)
  {
    newMatrix->lines[i] = NULL;
  }
  
  *m = newMatrix;
  return 0;
}

int matrix_destroy(Matrix* m)
{
  MatrixElem* curNode, *nextNode;
  for (uint32_t i = 0; i < m->ni; i++)
  {
    curNode = m->lines[i];
    while (curNode != NULL)
    {
      nextNode = curNode->next;
      free(curNode);
      curNode = nextNode;
    }
  }
  free(m->lines);
  free(m);

  return 0;
}

int matrix_print(const Matrix* m)
{
  printf("%" PRIu32 " %" PRIu32 "\n", m->ni, m->nj);

  for (MatrixElem *curNode = findNextNode(m, NULL); curNode != NULL; curNode = findNextNode(m, curNode))
  {
    printf("%" PRIu32 " %" PRIu32 " %.4g\n", curNode->i, curNode->j, curNode->val);
  }

  return 0;
}

int matrix_getelem( const Matrix* m, uint32_t i, uint32_t j, double *elem)
{
  if (i >= m->ni || j >= m->nj)
  {
    fprintf(stderr, "matrix_getelem: Tentando obter elemento em matriz de tamanho incorreto.\n");
    return 1;
  }
  
  MatrixElem *curNode = m->lines[i];
  while (curNode != NULL  && curNode->j < j)
  {
    curNode = curNode->next;
  }

  if (curNode != NULL && curNode->j == j)
  {
    *elem = curNode->val;
  }
  else
  {
    *elem = 0;
  }
  return 0;
}

int matrix_setelem(Matrix* m, uint32_t i, uint32_t j, double elem)
{
  if (i >= m->ni || j >= m->nj)
  {
    fprintf(stderr, "matrix_setelem: Tentando setar elemento em matriz de tamanho incorreto.\n");
    return 1;
  }

  MatrixElem *prevNode = NULL;
  MatrixElem *curNode = m->lines[i];
  while (curNode != NULL && curNode->j < j)
  {
    prevNode = curNode;
    curNode = curNode->next;
  }

  MatrixElem *newElem = createNewNode(i, j, elem);
  insertNewNode(m->lines[i], prevNode, curNode, newElem);
  
  return 0;
}

int matrix_sumelem( Matrix* m, uint32_t i, uint32_t j, double elem)
{
  if (i >= m->ni || j >= m->nj)
  {
    fprintf(stderr, "matrix_sumelem: Tentando setar elemento em matriz de tamanho incorreto.\n");
    return 1;
  }

  MatrixElem *prevNode = NULL;
  MatrixElem *curNode = m->lines[i];
  while (curNode != NULL && curNode->j < j)
  {
    prevNode = curNode;
    curNode = curNode->next;
  }

  /* node already exists*/
  if (curNode != NULL && curNode->j == j)
  {
    double v = curNode->val + elem;
    if (!floatIsZero(v))
    {
      curNode->val = v;
    }
    else
    {
      /* removes the element */
      prevNode->next = curNode->next;
      free(curNode);
    }
  }
  else
  {
    /* Node doesn't exist. If value is zero, do nothing. Else, insert the new node */
    if (!floatIsZero(elem))
    {
      MatrixElem *newElem = createNewNode(i, j, elem);
      insertNewNode(m->lines[i], prevNode, curNode, newElem);
    }
  }

  return 0;
}

int matrix_applyDirBC(Matrix *m, Vector *v, uint32_t idx, double val)
{
  MatrixElem *prevNode;
  MatrixElem *curNode;
  MatrixElem *nextNode;
  curNode = m->lines[idx];

  /* remove row from matrix */
  while (curNode != NULL)
  {
    nextNode = curNode->next;
    if (curNode->j == idx)
    {
      m->lines[idx] = curNode;
      curNode->next = NULL;
    }
    else
    {
      free(curNode);
    }
    curNode = nextNode;
  }

  /* searches each line for elements of the desired column, and pass them to the right. */
  for (uint32_t i = 0; i < m->ni; i++)
  {
    curNode = m->lines[i];
    prevNode = NULL;
    while (curNode != NULL && curNode->j < idx)
    {
      prevNode = curNode;
      curNode = curNode->next;
    }

    if (curNode != NULL && curNode->j == idx)
    {
      /* if its the node at the line, sets the right value*/
      if (i == idx)
      {
        v->val[i] = curNode->val * val;
      }
      /* if its not, multiplies it by the inputed value and removes it from the left */
      else
      {
        v->val[i] -= curNode->val * val;
        removeNode(m, prevNode, curNode);
      }
    }
  }

  return 0;
}

uint32_t matrix_countnnz(const Matrix* m)
{
  uint32_t nnz = 0;
  MatrixElem *curNode = findNextNode(m, NULL);
  while (curNode != NULL)
  {
    nnz++;
    curNode = findNextNode(m, curNode);
  }
  return nnz;
}

int matrixcrs_create(const Matrix* m, MatrixCRS **newCRSMatrix)
{
  /* only works for square matrix */
  if (m->ni != m->nj)
  {
    fprintf(stderr, "matrixcrs_create: Tentando converter para matriz não quadrada.\n");
    return 1;
  }
  uint32_t n = m->ni;
  uint32_t nnz = matrix_countnnz(m);
  
  /* allocates the required memory */
  MatrixCRS *newMatrix = sMalloc(sizeof(MatrixCRS));
  newMatrix->sa = sMalloc((nnz + 1) * sizeof(double));
  newMatrix->ija = sMalloc((nnz + 1) * sizeof(uint32_t));
  
  /* empties the firsts nodes */
  uint32_t k;
  for (k = 0; k < n; k++)
  {
    newMatrix->sa[k] = 0;
  }
  
  /* fill the matrix */
  /* Index to 1st row off-diagonal element, if any. */
  newMatrix->ija[0] = n+1;
  
  MatrixElem *curNode;
  uint32_t i = 0;
  k = n;
  curNode = findNextNode(m, NULL);
  while (curNode != NULL)
  {
    while (curNode->i > i)
    {
      i++;
      newMatrix->ija[i] = k + 1;
    }

    /* checks if the element is diagonal */
    if (curNode->i == curNode->j)
    {
      newMatrix->sa[curNode->i] = curNode->val;
    }
    else
    {
    /* element is not at diagonal */
      if (++k > nnz)
      {
        for (uint32_t juju = 0; juju < nnz; juju++)
        {
          printf("%"PRIu32" %"PRIu32" %lf\n", juju, newMatrix->ija[juju], newMatrix->sa[juju]);
        }
        fprintf(stderr, "matrixcrs_create: tamanho insuficiente");
        return 1;
      };
      newMatrix->sa[k] = curNode->val;
      newMatrix->ija[k] = curNode->j;
    }

    curNode = findNextNode(m, curNode);
  }
    
  newMatrix->ija[i+1] = k+1;
  
  *newCRSMatrix = newMatrix;
  return 0;
}

int matrixcrs_destroy(MatrixCRS* m)
{
  free(m->sa);
  free(m->ija);
  free(m);
  return 0;
}

int matrixcrs_multiplyVector(const MatrixCRS* m, const Vector *v, Vector *r)
{
  uint32_t n = m->ija[0] - 1;
  
  if (n != v->len || n != r->len)
  {
    fprintf(stderr, "matrixcrs_multiplyVector: Matriz e vetor de tamanhos diferentes: m = %" PRIu32 "x%" PRIu32 ", v = %" PRIu32 ", r = %" PRIu32".\n", n, n, v->len, r->len);
    return 1;
  }
  
  uint32_t i, k;
  for (i = 0; i < n; i++)
  {
    r->val[i] = m->sa[i] * v->val[i]; /* diagonal term */
    for(k = m->ija[i]; k < m->ija[i+1]; k++)
    {
      r->val[i] += m->sa[k] * v->val[m->ija[k]];
    }
  }
  return 0;
}

int vector_create(Vector **v, uint32_t n)
{
  Vector *newVector = sMalloc(sizeof(Vector));
  newVector->len = n;
  newVector->val = sMalloc(n * sizeof(double));
  *v = newVector;
  return 0;
}

int vector_createzero(Vector **v, uint32_t n)
{
  Vector *newVector = sMalloc(sizeof(Vector));
  newVector->len = n;
  newVector->val = calloc(n, sizeof(double));
  
  *v = newVector;
  return 0;
}

int vector_destroy(Vector *v)
{
  free(v->val);
  free(v);
  return 0;
}

double vector_internalProduct(const Vector *v, const Vector *w)
{
  if (v->len != w->len)
  {
    fprintf(stderr, "vector_internalProduct: vector sizes inconsistencies.\n");
    return DBL_MAX;
  }
  
  uint32_t n = v->len;
  double sum = 0;
  uint32_t i;
  for (i = 0; i < n; i++)
  {
    sum += v->val[i] * w->val[i];
  }
  return sum;
}

int vector_copy(const Vector *v, Vector *r)
{
  if (v->len != r->len)
  {
    fprintf(stderr, "vector_copy: vector sizes inconsistencies.\n");
    return 1;
  }
  uint32_t n = v->len;
  
  uint32_t i;
  for (i = 0; i < n; i++)
  {
    r->val[i] = v->val[i];
  }
  return 0;
}

int vector_add (const Vector *v, const Vector *w, Vector *r)
{
  uint32_t n = v->len;
  if (n != w->len || n != r->len)
  {
    fprintf(stderr, "addVector: vetores de tamanhos diferentes: v = %" PRIu32 ", w = %" PRIu32 ", r = %" PRIu32".\n", v->len, w->len, r->len);
    return 1;
  }
  
  uint32_t i;
  for (i = 0; i < n; i++)
  {
    r->val[i] = v->val[i] + w->val[i];
  }
  return 0;
}

int vector_subtract (const Vector *v, const Vector *w, Vector *r)
{
  uint32_t n = v->len;
  if (n != w->len || n != r->len)
  {
    fprintf(stderr, "addVector: vetores de tamanhos diferentes: v = %" PRIu32 ", w = %" PRIu32 ", r = %" PRIu32".\n", v->len, w->len, r->len);
    return 1;
  }
  
  uint32_t i;
  for (i = 0; i < n; i++)
  {
    r->val[i] = v->val[i] - w->val[i];
  }
  return 0;
}

int vector_addVectorScalar (const Vector *v, const Vector *w, double alfa, Vector *r)
{
  uint32_t n = v->len;
  if (n != w->len || n != r->len)
  {
    fprintf(stderr, "addVector: vetores de tamanhos diferentes: v = %" PRIu32 ", w = %" PRIu32 ", r = %" PRIu32".\n", v->len, w->len, r->len);
    return 1;
  }
  
  uint32_t i;
  for (i = 0; i < n; i++)
  {
    r->val[i] = v->val[i] + (alfa * w->val[i]);
  }
  return 0;
}

uint32_t vector_size(const Vector *v)
{
  return v->len;
}

int vector_setelem(Vector *v, uint32_t i, double elem)
{
  if (i >= v->len)
  {
    fprintf(stderr, "vector_setelem: elemento fora do vetor.\n");
    return 1;
  }
  
  v->val[i] = elem;
  return 0;
}

int vector_sumelem(Vector *v, uint32_t i, double elem)
{
  if (i >= v->len)
  {
    fprintf(stderr, "vector_setelem: elemento fora do vetor.\n");
    return 1;
  }
  
  v->val[i] += elem;
  return 0;
}

double vector_getelem(const Vector *v, uint32_t i)
{
  if (i >= v->len)
  {
    fprintf(stderr, "vector_getelem: elemento fora do vetor.\n");
    return 1;
  }
  return v->val[i];
}

int vector_print(const Vector *v)
{
  uint32_t i;
  for (i = 0; i < v->len; i++)
  {
    printf("%" PRIu32 ": %#.3g\n", i, v->val[i]);
  }
  return 0;
}

internal inline MatrixElem* createNewNode(uint32_t i, uint32_t j, double val)
{
  MatrixElem* newElem = sMalloc(sizeof(MatrixElem));
  newElem->i = i;
  newElem->j = j;
  newElem->val = val;
  return newElem;
}

internal inline void removeNode(Matrix *m, MatrixElem *prevNode, MatrixElem *curNode)
{
  if (prevNode == NULL)
  {
    m->lines[curNode->i] = curNode->next;
  }
  else
  {
    prevNode->next = curNode->next;
  }
  free(curNode);
}

internal inline void insertNewNode(MatrixElem *headNode, MatrixElem *prevNode, MatrixElem *curNode, MatrixElem *newElem)
{
  if (prevNode == NULL)
  {
    newElem->next = curNode;
    headNode = newElem;
  }
  else
  {
    newElem->next = prevNode->next;
    prevNode->next = newElem;
  }
}

internal MatrixElem* findNextNode(const Matrix* m, const MatrixElem* node)
{
  uint32_t i;
  if (node == NULL)
  {
    i = 0;
  }
  else
  {
    /* searches the line. If it's the last row, continues from the next*/
    MatrixElem* curNode = node->next;
    if (curNode != NULL)
    {
      return curNode;
    }
    i = node->i + 1;
  }

  while (i < m->ni)
  {
    if (m->lines[i] != NULL)
    {
      return m->lines[i];
    }
    i++;
  }

  return NULL;
}


internal int floatIsZero(double x)
{
  if (fabs(x) < DBL_EPSILON) return 1;
  return 0;
}

internal inline void *sMalloc(size_t memsize)
{
    void *allocMem = malloc(memsize);
    if(!allocMem && memsize)
    {
        fprintf(stderr, "Could not allocate memory!\n");
        exit(-1);
    }
    return allocMem;
}