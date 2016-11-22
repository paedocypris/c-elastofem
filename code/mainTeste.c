#include <stdio.h>
#include "matrix.h"

int main( void ) {
/* Inicializacao da aplicacao ... */
  Matrix *A=NULL;
  Matrix *B=NULL;
  Matrix *C=NULL;
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
  if ( !matrix_add( A, B, &C ) ) {
    matrix_print( C );
  }
  else
    fprintf( stderr, "Erro na soma C=A+B.\n" );
  matrix_destroy( C );
  if ( !matrix_subtract( A, B, &C ) ) {
    matrix_print( C );
  }
  else
    fprintf( stderr, "Erro na subtração C=A-B.\n" );
  matrix_destroy( C );
  
  matrix_destroy( C );
  matrix_destroy( A );
  matrix_destroy( B );
  return 0;
}