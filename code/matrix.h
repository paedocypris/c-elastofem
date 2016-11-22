#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <stdint.h>

typedef struct matrix Matrix;

int matrix_create(Matrix** m, int64_t ni, int64_t nj);
int matrix_destroy(Matrix* m);
int matrix_print(const Matrix* m);
int matrix_add(const Matrix* m, const Matrix* n, Matrix** r);
int matrix_subtract(const Matrix *m, const Matrix* n, Matrix** r);
int matrix_multiply(const Matrix* m, const Matrix* n, Matrix** r);
int matrix_transpose(const Matrix* m, Matrix** r);
int matrix_getelem( const Matrix* m, int64_t x, int64_t y, double *elem);
int matrix_setelem( Matrix* m, int64_t x, int64_t y, double elem);
int matrix_sumelem( Matrix* m, int64_t x, int64_t y, double elem);

#endif
