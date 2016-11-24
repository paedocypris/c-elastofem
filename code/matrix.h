#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <stdint.h>

typedef struct matrix Matrix;

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

#endif
