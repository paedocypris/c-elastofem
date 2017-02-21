#include <stdio.h>
#include "matrix.h"

int main( void ) {
/* Inicializacao da aplicacao ... */
  Matrix *A=NULL;
  Matrix *B=NULL;
  if( !matrix_create( &A, 3, 3 ) )
  {
    matrix_setelem( A, 0, 0, 1);
    matrix_setelem( A, 0, 1, 2);
    matrix_setelem( A, 1, 1, 5);
    matrix_print( A );
  }
  else {
    fprintf( stderr, "Erro na alocacao de A como listas encadeadas.\n" );
    return 1;
  }
  if( !matrix_create( &B, 3, 3 ) )
  {
    matrix_setelem( B, 0, 1, 2);
    matrix_setelem( B, 1, 0, 1);
    matrix_setelem( B, 1, 1, 5);
    matrix_setelem( B, 2, 1, 1);
    matrix_print( B );
  }
  else {
    fprintf( stderr, "Erro na alocacao de B como listas encadeadas.\n" );
    return 1;
  }
  
  matrix_destroy( A );
  matrix_destroy( B );

  if (!matrix_create(&A, 5, 5))
  {
	  matrix_setelem(A, 0, 0, 3);
	  matrix_setelem(A, 0, 2, 1);
	  matrix_setelem(A, 1, 1, 4);
	  matrix_setelem(A, 2, 1, 7);
	  matrix_setelem(A, 2, 2, 5);
	  matrix_setelem(A, 2, 3, 9);
	  matrix_setelem(A, 3, 4, 2);
	  matrix_setelem(A, 4, 3, 6);
	  matrix_setelem(A, 4, 4, 5);
	  matrix_print(A);
  }
  MatrixCRS *mCrs;
  matrixcrs_create(A, &mCrs);
  return 0;
}