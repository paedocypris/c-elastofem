#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <stdint.h>

typedef struct matrix Matrix;
typedef struct matrixcrs MatrixCRS;
typedef struct vector Vector;

int matrix_create(Matrix** m, int64_t ni, int64_t nj);
int matrix_destroy(Matrix* m);
int matrix_print(const Matrix* m);
int matrix_add(const Matrix* m, const Matrix* n, Matrix** r);
int matrix_addE(const Matrix* m, const Matrix* n, Matrix* r);
int matrix_subtract(const Matrix *m, const Matrix* n, Matrix** r);
int matrix_subtractE(const Matrix *m, const Matrix* n, Matrix* r);
int matrix_multiply(const Matrix* m, const Matrix* n, Matrix** r);
int matrix_multiplyE(const Matrix* m, const Matrix* n, Matrix* r);
int matrix_transpose(const Matrix* m, Matrix** r);
int matrix_transposeE(const Matrix* m, Matrix* r);
int matrix_getelem( const Matrix* m, int64_t x, int64_t y, double *elem);
int matrix_setelem( Matrix* m, int64_t x, int64_t y, double elem);
int matrix_sumelem( Matrix* m, int64_t x, int64_t y, double elem);
int matrix_scalarmult( double alfa, const Matrix* m, Matrix** r);
double matrix_internalProduct(const Matrix* m, const Matrix* n);
int matrix_copy(const Matrix* m, Matrix** r);

int64_t matrix_ni(const Matrix* m);
int64_t matrix_nj(const Matrix* m);

int matrixcrs_create(const Matrix* m, MatrixCRS **newCRSMatrix);
int matrixcrs_multiplyVector(const MatrixCRS* m, const Vector *v, Vector *r);

int vector_create(Vector **v, uint32_t n);
int vector_destroy(Vector *v);
double vector_internalProduct(const Vector *v, const Vector *w);
int vector_copy(const Vector *v, Vector *r);
int vector_add (const Vector *v, const Vector *w, Vector *r);
int vector_subtract (const Vector *v, const Vector *w, Vector *r);
int vector_addVectorScalar (const Vector *v, const Vector *w, double alfa, Vector *r);
uint32_t vector_size(const Vector *v);
int vector_setelem(Vector *v, uint32_t i, double elem);
int vector_sumelem(Vector *v, uint32_t i, double elem);
int vector_print(const Vector *v);

#endif
