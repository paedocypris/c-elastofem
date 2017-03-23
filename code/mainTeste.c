#include <stdio.h>
#include "stdMatrix.h"

int main( void ) {
/* Inicializacao da aplicacao ... */
  StdMatrix *A=NULL;
  StdMatrix *invA=NULL;
  if( !stdmatrix_create( &A, 3, 3 ) )
  {
    stdmatrix_setelem( A, 0, 0, 1);
    stdmatrix_setelem( A, 0, 1, 2);
    stdmatrix_setelem( A, 0, 2, 1);

    stdmatrix_setelem( A, 1, 0, 0);
    stdmatrix_setelem( A, 1, 1, 3);
    stdmatrix_setelem( A, 1, 2, 5);

    stdmatrix_setelem( A, 2, 0, -2);
    stdmatrix_setelem( A, 2, 1, 4.12);
    stdmatrix_setelem( A, 2, 2, -1);
  }
  else {
    fprintf( stderr, "Erro na alocacao de A como listas encadeadas.\n" );
    return 1;
  }

  stdmatrix_create(&invA, 3, 3);
  stdmatrix_invert(A, invA);


  
  stdmatrix_destroy( A );
  stdmatrix_destroy(invA);
  return 0;
}