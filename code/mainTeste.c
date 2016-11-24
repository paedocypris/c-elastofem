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
  
  if ( !matrix_scalarmult( 1.4, A, &C))
  {
    matrix_print( C );
  }
  else
    fprintf( stderr, "Erro na multiplicação por escalar alfa*A.\n" );
  matrix_destroy( C );
  
  /* create vectors, test internal product */
  Matrix *va = NULL;
  Matrix *vb = NULL;
  if( !matrix_create( &va, 4, 1 ) )
  {
    matrix_setelem( va, 0, 0, 1);
    matrix_setelem( va, 1, 0, 1.5);
    matrix_setelem( va, 2, 0, -2);
    matrix_print( va );
  }
  else {
    fprintf( stderr, "Erro na alocacao de va como listas encadeadas.\n" );
    return 1;
  }
  
  if( !matrix_create( &vb, 4, 1 ) )
  {
    matrix_setelem( vb, 0, 0, 2);
    matrix_setelem( vb, 1, 0, -1.3);
    matrix_setelem( vb, 3, 0, 4);
    matrix_print( vb );
  }
  else {
    fprintf( stderr, "Erro na alocacao de va como listas encadeadas.\n" );
    return 1;
  }
  double ip = matrix_internalProduct(va, vb);
  if (ip - 0.05 >= 1e-5)
  {
    fprintf (stderr, "Erro no cálculo do produto interno");
  }
  ip = matrix_internalProduct(vb, va);
  if (ip - 0.05 >= 1e-5)
  {
    fprintf (stderr, "Erro no cálculo do produto interno");
  }
  
  matrix_destroy(va);
  matrix_destroy(vb);
  
  matrix_destroy( A );
  matrix_destroy( B );
  return 0;
}