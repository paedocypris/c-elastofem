#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <stdint.h>
#include "stdMatrix.h"

typedef struct matrix Matrix;
typedef struct matrixcrs MatrixCRS;
typedef struct vector Vector;

int matrix_create(Matrix** m, uint32_t ni, uint32_t nj);
int matrix_destroy(Matrix* m);
int matrix_print(const Matrix* m);

int matrix_issimetric(const Matrix *m);
int matrix_getelem( const Matrix* m, uint32_t i, uint32_t j, double *elem);
int matrix_setelem( Matrix* m, uint32_t i, uint32_t j, double elem);
int matrix_sumelem( Matrix* m, uint32_t i, uint32_t j, double elem);
int matrix_applyDirBC(Matrix *m, Vector *v, uint32_t idx, double val);
uint32_t matrix_countnnz(const Matrix* m);

int matrixcrs_create(const Matrix* m, MatrixCRS **newCRSMatrix);
int matrixcrs_destroy(MatrixCRS* m);
int matrixcrs_multiplyVector(const MatrixCRS* m, const Vector *v, Vector *r);

int vector_create(Vector **v, uint32_t n);
int vector_createzero(Vector **v, uint32_t n);
int vector_destroy(Vector *v);
double vector_internalProduct(const Vector *v, const Vector *w);
int vector_copy(const Vector *v, Vector *r);
int vector_add (const Vector *v, const Vector *w, Vector *r);
int vector_subtract (const Vector *v, const Vector *w, Vector *r);
int vector_addVectorScalar (const Vector *v, const Vector *w, double alfa, Vector *r);
uint32_t vector_size(const Vector *v);
int vector_setelem(Vector *v, uint32_t i, double elem);
int vector_sumelem(Vector *v, uint32_t i, double elem);
double vector_getelem(const Vector *v, uint32_t i);
int vector_print(const Vector *v);

#endif
