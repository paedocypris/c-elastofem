#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include "matrix.h"

#define internal static
#define FLT_TOLERANCE 1e-5

typedef int8_t int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;

typedef uint8_t uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;


struct matrix {
   struct matrix* right;
   struct matrix* below;
   int64 line;
   int64 column;
   double info;
};

internal Matrix* createMatrixHeaders(int64 ni, int64 nj);
internal int readLineAndUpdateMatrix(Matrix* m);
internal Matrix* insertNewEmptyNodeToMatrix();
internal Matrix* insertNewNodeToMatrix(Matrix* m, int64 i, int64 j, double info);
internal Matrix* insertNewNodeToMatrixPrevNodes(Matrix* prevRow, Matrix* prevColumn, int64 i, int64 j, double info);
internal void removeNodeFromMatrix(Matrix* node, Matrix* prevRow, Matrix* prevColumn);
internal int64 getMatrixIndex(const Matrix* m, const Matrix* node);

internal Matrix* findPrevRowNode(const Matrix* m, int64 i, int64 j);
internal Matrix* findPrevRowNodeGivenNode(const Matrix* node);
internal Matrix* findPrevColumnNode(const Matrix* m, int64 i, int64 j);
internal Matrix* findPrevColumnNodeGivenNode(const Matrix* node);
internal Matrix* findNode(const Matrix* m, int64 i, int64 j);
internal Matrix* findNextNode(const Matrix* m, const Matrix* node);

internal int floatIsZero(double x);

internal Matrix* createMatrixHeaders(int64 ni, int64 nj)
{
  Matrix* prevNode;
  Matrix* headNode;
  
  /* read the matrix's size */
  if (ni <= 0 || nj <= 0)
  {
    fprintf(stderr, "Tamanho das matrizes incorreto.\n" );
    return NULL;
  }
  
  /* given the matrix size, create the head nodes */
  /* the line number and column number indicates the matrix's size */
  headNode = insertNewEmptyNodeToMatrix();
  headNode->right = headNode;
  headNode->below = headNode;
  headNode->line = -ni;
  headNode->column = -nj;
  
  /* create line head nodes */
  prevNode = headNode;
  for(int64 i = 0; i < ni; i++)
  {
    Matrix* newRowNode = insertNewEmptyNodeToMatrix();
    newRowNode->right = newRowNode;
    newRowNode->below = headNode;
    newRowNode->line = i;
    newRowNode->column = -1;
    
    prevNode->below = newRowNode;
    prevNode = newRowNode;
  }
  
  /* create column head nodes */
  prevNode = headNode;
  for(int64 j = 0; j < nj; j++)
  {
    Matrix* newColumnNode = insertNewEmptyNodeToMatrix();
    newColumnNode->right = headNode;
    newColumnNode->below = newColumnNode;
    newColumnNode->line = -1;
    newColumnNode->column = j;
    
    prevNode->right = newColumnNode;
    prevNode = newColumnNode;
  }
  
  return headNode;
}

/* returns 1 if reached the end of the matrix representation,
           0 if succeeded,
           -1 if failed */
internal int readLineAndUpdateMatrix(Matrix* m)
{
  int64 x, y;
  double value;
  scanf("%"SCNd64, &x);
  
  /* read first item, if not -1, it's a valid row, else, the user stopped issuing values */
  if (x == -1) return 1;
  
  scanf("%" SCNd64 "%lf" , &y, &value);
  
  /* matrix representation is zero based, the user interface is one based */
  return matrix_setelem(m, x, y, value);
}

internal Matrix* insertNewEmptyNodeToMatrix()
{
  Matrix* newNode = (Matrix*) malloc(sizeof(Matrix));
  if (newNode == NULL)
  {
    fprintf(stderr, "Erro alocando memória.\n" );
    return NULL;
  }
  return newNode;
}

internal Matrix* insertNewNodeToMatrix(Matrix* m, int64 i, int64 j, double info)
{
  Matrix* prevRow = findPrevRowNode(m, i, j);
  Matrix* prevColumn = findPrevColumnNode(m, i, j);
  
  return insertNewNodeToMatrixPrevNodes(prevRow, prevColumn, i, j, info);
}

internal Matrix* insertNewNodeToMatrixPrevNodes(Matrix* prevRow, Matrix* prevColumn, int64 i, int64 j, double info)
{
  Matrix* newNode = insertNewEmptyNodeToMatrix();
  newNode->right = prevRow->right;
  newNode->below = prevColumn->below;
  newNode->line = i;
  newNode->column = j;
  newNode->info = info;
  
  prevRow->right = newNode;
  prevColumn->below = newNode;
  
  return newNode;
}

internal void removeNodeFromMatrix(Matrix* node, Matrix* prevRow, Matrix* prevColumn)
{
  prevRow->right = node->right;
  prevColumn->below = prevColumn->below;
  free(node);
}

internal Matrix* findPrevRowNode(const Matrix* m, int64 i, int64 j)
{
  Matrix* prevRowNode;
  Matrix* currentNode = m->below;
  
  /* searches the head node corresponding the line */
  while (currentNode->line != i && currentNode != m)
  {
    currentNode = currentNode->below;
  }
  
  if (currentNode == m)
  {
    fprintf(stderr, "Matriz mal formada.\n" );
    return NULL;
  }
  
  prevRowNode = currentNode;
  currentNode = currentNode->right;
  while (currentNode->column < j && currentNode->column != -1)
  {
    prevRowNode = currentNode;
    currentNode = currentNode->right;
  }
  
  return prevRowNode;
}

internal Matrix* findPrevRowNodeGivenNode(const Matrix* node)
{
  Matrix* currentNode;
  
  currentNode = (Matrix*) node;
  while (currentNode->right != node)
  {
    currentNode = currentNode->right;
  }
  return currentNode;
}

internal Matrix* findPrevColumnNode(const Matrix* m, int64 i, int64 j)
{
  Matrix* prevColumnNode;
  Matrix* currentNode = m->right;
  
  /* searches the head node corresponding the line */
  while (currentNode->column != j && currentNode != m)
  {
    currentNode = currentNode->right;
  }
  
  if (currentNode == m)
  {
    fprintf(stderr, "Matriz mal formada.\n" );
    return NULL;
  }
  
  prevColumnNode = currentNode;
  currentNode = currentNode->below;
  while (currentNode->line < i && currentNode->line != -1)
  {
    prevColumnNode = currentNode;
    currentNode = currentNode->below;
  }
  
  return prevColumnNode;
}

internal Matrix* findPrevColumnNodeGivenNode(const Matrix* node)
{
  Matrix* currentNode;
  
  currentNode = (Matrix*) node;
  while (currentNode->below != node)
  {
    currentNode = currentNode->below;
  }
  return currentNode;
}

internal Matrix* findNode(const Matrix* m, int64 i, int64 j)
{
  Matrix *prevNode;
  if (i < j)
  {
    prevNode = findPrevRowNode(m, i, j);
    if (prevNode->right->column == j)
    {
      return prevNode->right;
    }
  }
  else
  {
    prevNode = findPrevColumnNode(m, i, j);
    if (prevNode->below->line == i)
    {
      return prevNode->below;
    }
  }
  return NULL;
}

/* skips head nodes */
internal Matrix* findNextNode(const Matrix* m, const Matrix* node)
{
  const Matrix* currentNode;
  
  // find the first node
  if (node == m)
  {
    currentNode = m->below;
  }
  else
  {
    currentNode = node;
  }
  
  while (currentNode != m)
  {
    currentNode = currentNode->right;
    if (currentNode->column != -1)
    {
      return (Matrix*) currentNode;
    }
    currentNode = currentNode->below;
  }
    /* didn't find any nodes */
  return NULL;
}

/* if m is NULL, returns -1 */
internal int64 getMatrixIndex(const Matrix* m, const Matrix* node)
{
  if (node == NULL) return -1;
  return node->line*(-m->column) + node->column;
}

internal int floatIsZero(double x)
{
  if (-FLT_TOLERANCE < x && x < FLT_TOLERANCE) return 1;
  return 0;
}

int matrix_create(Matrix** m, int64 ni, int64 nj)
{
  Matrix* headNode = createMatrixHeaders(ni, nj);
  if (headNode == NULL)
  {
    return -1;
  }
  *m = headNode;
  return 0;
}

int matrix_destroy(Matrix* m)
{
  Matrix* currentNode;
  Matrix* currentHeadNode;
  Matrix* nextNode;
  int matrixNotDestroyed;
  
  if (m == NULL)
  {
    fprintf(stderr, "Tentando destruir matriz não definidas.\n");
    return -1;
  }
  
  matrixNotDestroyed = 1;
  currentNode = m->below;
  while (matrixNotDestroyed)
  {
    currentHeadNode = currentNode;
    currentNode = currentNode->right;
    while (currentHeadNode != currentNode)
    {
      nextNode = currentNode->right;
      free(currentNode);
      currentNode = nextNode;
    }
    currentNode = currentNode->below;
    if (currentHeadNode == m)
    {
      matrixNotDestroyed = 0;
    }
    free(currentHeadNode);
  }
  return 0;
}

int matrix_print(const Matrix* m)
{
  const Matrix* currentNode;
  
  if (m == NULL)
  {
    fprintf(stderr, "Tentando imprimir matriz não definidas.\n");
    return -1;
  }
  
  currentNode = m;
  
  fprintf( stdout, "%" PRId64 " %" PRId64 "\n", -currentNode->line, -currentNode->column);
  
  /* searches for elements row by row */
  for (currentNode = findNextNode(m, currentNode); currentNode != NULL; currentNode = findNextNode(m, currentNode))
  {
    fprintf( stdout, "%" PRId64 " %" PRId64 " %f\n", currentNode->line, currentNode->column, currentNode->info);
  }
  
  fprintf( stdout, "0\n\n");
  
  return 0;
}

int matrix_add(const Matrix* m, const Matrix* n, Matrix** r)
{
  const Matrix *currentLowerNode, *currentHigherNode, *lowerHeadNode, *higherHeadNode, *proxyNode;
  Matrix* headNodeR;
  int64 ni, nj;
  int64 lowerNodeIndex, higherNodeIndex;
  
  if (m == NULL || n == NULL)
  {
    fprintf(stderr, "Tentando somar uma ou duas matrizes não definidas.\n");
    return -1;
  }
  
  if (m->line != n->line || m->column != n->column)
  {
    fprintf( stderr, "Soma não definida: Matriz M: %" PRId64 "x%" PRId64 "; Matriz N: %" PRId64 "x%" PRId64 "\n", -m->line, -m->column, -n->line, -n->column);
    return -1;
  }
  ni = -m->line;
  nj = -m->column;
  
  headNodeR = createMatrixHeaders(ni, nj);
  
  currentLowerNode = findNextNode(m, m);
  lowerHeadNode = m;
  currentHigherNode = findNextNode(n, n);
  higherHeadNode = n;
  
  while (currentLowerNode != NULL && currentHigherNode != NULL)
  {
    lowerNodeIndex = getMatrixIndex(lowerHeadNode, currentLowerNode);
    higherNodeIndex = getMatrixIndex(higherHeadNode, currentHigherNode);
    
    while (currentHigherNode == NULL || (currentLowerNode != NULL && lowerNodeIndex < higherNodeIndex)) 
    {
      matrix_setelem(headNodeR, currentLowerNode->line, currentLowerNode->column, currentLowerNode->info);
      currentLowerNode = findNextNode(lowerHeadNode, currentLowerNode);
      lowerNodeIndex = getMatrixIndex(lowerHeadNode, currentLowerNode);
    }
    if (lowerNodeIndex == higherNodeIndex)
    {
      matrix_setelem(headNodeR, currentLowerNode->line, currentLowerNode->column, currentLowerNode->info + currentHigherNode->info);
      currentLowerNode = findNextNode(lowerHeadNode, currentLowerNode);
      currentHigherNode = findNextNode(higherHeadNode, currentHigherNode);
    }
    else
    {
      /* switch the lower and higher nodes */
      proxyNode = currentLowerNode;
      currentLowerNode = currentHigherNode;
      currentHigherNode = proxyNode;
      
      proxyNode = lowerHeadNode;
      lowerHeadNode = higherHeadNode;
      higherHeadNode = proxyNode;
    }
  }
  
  *r = headNodeR;
  return 0;
}

int matrix_subtract(const Matrix* m, const Matrix* n, Matrix** r)
{
  const Matrix *currentLowerNode, *currentHigherNode, *lowerHeadNode, *higherHeadNode, *proxyNode;
  Matrix* headNodeR;
  int64 ni, nj;
  int64 lowerNodeIndex, higherNodeIndex;
  
  if (m == NULL || n == NULL)
  {
    fprintf(stderr, "Tentando subtrair uma ou duas matrizes não definidas.\n");
    return -1;
  }
  
  if (m->line != n->line || m->column != n->column)
  {
    fprintf( stderr, "Subtração não definida: Matriz M: %" PRId64 "x%" PRId64 "; Matriz N: %" PRId64 "x%" PRId64 "\n", -m->line, -m->column, -n->line, -n->column);
    return -1;
  }
  ni = -m->line;
  nj = -m->column;
  
  headNodeR = createMatrixHeaders(ni, nj);
  
  currentLowerNode = findNextNode(m, m);
  lowerHeadNode = m;
  currentHigherNode = findNextNode(n, n);
  higherHeadNode = n;
  
  /* shows if the current lower node belongs to the matrix m
     mIsLower = 1 if positive
     mIsLower = -1 if negative */
  int mIsLower = 1;
  
  while (currentLowerNode != NULL && currentHigherNode != NULL)
  {
    lowerNodeIndex = getMatrixIndex(lowerHeadNode, currentLowerNode);
    higherNodeIndex = getMatrixIndex(higherHeadNode, currentHigherNode);
    
    while (currentHigherNode == NULL || (currentLowerNode != NULL && lowerNodeIndex < higherNodeIndex)) 
    {
      matrix_setelem(headNodeR, currentLowerNode->line, currentLowerNode->column, mIsLower*currentLowerNode->info);
      currentLowerNode = findNextNode(lowerHeadNode, currentLowerNode);
      lowerNodeIndex = getMatrixIndex(lowerHeadNode, currentLowerNode);
    }
    if (lowerNodeIndex == higherNodeIndex)
    {
      matrix_setelem(headNodeR, currentLowerNode->line, currentLowerNode->column, mIsLower*(currentLowerNode->info - currentHigherNode->info));
      currentLowerNode = findNextNode(lowerHeadNode, currentLowerNode);
      currentHigherNode = findNextNode(higherHeadNode, currentHigherNode);
    }
    else
    {
      /* switch the lower and higher nodes */
      proxyNode = currentLowerNode;
      currentLowerNode = currentHigherNode;
      currentHigherNode = proxyNode;
      
      proxyNode = lowerHeadNode;
      lowerHeadNode = higherHeadNode;
      higherHeadNode = proxyNode;
      
      mIsLower = mIsLower*-1;
    }
  }
  
  *r = headNodeR;
  return 0;
}

int matrix_multiply(const Matrix* m, const Matrix* n, Matrix** r)
{
  const Matrix *currentMNode, *currentNNode;
  Matrix *currentRNode;
  Matrix *headNodeR;
  
  if (m == NULL || n == NULL)
  {
    fprintf(stderr, "Tentando multiplicar uma ou duas matrizes não definidas.\n");
    return -1;
  }

  if (m->column != n->line)
  {
    fprintf( stderr, "Multiplicação de matrizes com dimensões incorretas, m: %" PRId64 " %" PRId64 ", n: %" PRId64 " %" PRId64 ".\n", -m->line, -m->column, -n->line, -n->column);
    return -1;
  }
  
  headNodeR = createMatrixHeaders(-m->column, -n->line);
  
  /* assumes the matrix have about the same number of non-zero elements */
  currentNNode = findNextNode(n, n);
  while (currentNNode != NULL)
  {
    currentMNode = findPrevColumnNode(m, 0, currentNNode->line)->below;
    while (currentMNode->line != -1)
    {
      currentRNode = findNode(headNodeR, currentMNode->line, currentNNode->column);
      if (currentRNode == NULL)
      {
        matrix_setelem(headNodeR, currentMNode->line, currentNNode->column, currentMNode->info * currentNNode->info);
      }
      else
      {
        matrix_setelem(headNodeR, currentMNode->line, currentNNode->column, currentRNode->info + currentMNode->info * currentNNode->info);
      }
      currentMNode = currentMNode->below;
    }
    currentNNode = findNextNode(n, currentNNode);
  }
  *r = headNodeR;
  return 0;
}

int matrix_transpose(const Matrix* m, Matrix** r)
{
  const Matrix *currentNode;
  Matrix *headNodeR;
  int64 ni, nj;
  
  if (m == NULL)
  {
    fprintf(stderr, "Tentando transpor matriz não definidas.\n");
    return -1;
  }
  
  ni = -m->line;
  nj = -m->column;
  
  headNodeR = createMatrixHeaders(nj, ni);
  
  currentNode = findNextNode(m, m);
  while (currentNode != NULL)
  {
    matrix_setelem(headNodeR, currentNode->column, currentNode->line, currentNode->info);
    currentNode = findNextNode(m, currentNode);
  }
  *r = headNodeR;
  return 0;
}

int matrix_getelem( const Matrix* m, int64 x, int64 y, double *elem)
{
  int64 i, j;
  const Matrix* desiredNode;
  
  if (m == NULL)
  {
    fprintf(stderr, "Tentando obter elemento de matriz não definida.\n");
    return -1;
  }
  
  i = x;
  j = y;
  
  if (x > -m->line || y > -m->column || x < 1 || y < 1)
  {
    fprintf( stderr, "Elemento %" PRId64 "x%" PRId64 " fora da matriz.\n", x, y);
    return -1;
  }
  
  desiredNode = findNode(m, i, j);
  
  if (desiredNode != NULL)
  {
    *elem = desiredNode->info;
    return 0;
  }
  else
  {
    *elem = 0.0;
    return 0;
  }
}

/* returns 0 if succeeded, -1 if failed */
int matrix_setelem( Matrix* m, int64 x, int64 y, double elem)
{
  Matrix *prevRowNode;
  Matrix *prevColumnNode;
  int64 i, j;
  
  if (m == NULL)
  {
    fprintf(stderr, "Tentando atualizar o valor dentro de uma matriz não definida.\n");
    return -1;
  }
  
  i = x; j = y;
  
  /* check if the matriz is in range */
  if (i < 0 || i >= -m->line || j < 0 || j >= -m->column)
  {
    fprintf(stderr, "Elemento fora da matriz.\n" );
    return -1;
  }
  
  prevRowNode = findPrevRowNode(m, i, j);
  if (floatIsZero(elem))
  {
    if (prevRowNode->right->column != j)
    {
      /* element doesn't exist, do nothing */
      return 0;
    }
    else
    {
      /* element exists, remove it from the matrix */
      prevColumnNode = findPrevColumnNodeGivenNode(prevRowNode->right);
      removeNodeFromMatrix(prevColumnNode->below, prevRowNode, prevColumnNode);
      return 0;
    }
  }
  else
  {
    if (prevRowNode->right->column != j)
    {
      /* element doesn't exist, insert a new node */
      prevColumnNode = findPrevColumnNode(m, i, j);
      if (insertNewNodeToMatrixPrevNodes(prevRowNode, prevColumnNode, i, j, elem) == NULL)
      {
        return -1;
      }
      
      return 0;
    }
    else
    {
      /* element exists, update its info */
      prevRowNode->right->info = elem;
      return 0;
    }
  }
}

int matrix_sumelem( Matrix* m, int64 x, int64 y, double elem)
{
  Matrix *prevRowNode;
  Matrix *prevColumnNode;
  int64 i, j;

  if (m == NULL)
  {
    fprintf(stderr, "Tentando atualizar o valor dentro de uma matriz não definida.\n");
    return -1;
  }
  
  i = x; j = y;
  
  /* check if the matriz is in range */
  if (i < 0 || i >= -m->line || j < 0 || j >= -m->column)
  {
    fprintf(stderr, "Elemento fora da matriz.\n" );
    return -1;
  }
  
  if (floatIsZero(elem))
  {
    /* doesn't change the element */
    return 0;
  }
  
  prevRowNode = findPrevRowNode(m, i, j);
  if (prevRowNode->right->column != j)
  {
    /* element doesn't exist, insert a new node */
    prevColumnNode = findPrevColumnNode(m, i, j);
    if (insertNewNodeToMatrixPrevNodes(prevRowNode, prevColumnNode, i, j, elem) == NULL)
    {
      return -1;
    }
    
    return 0;
  }
  else
  {
    double temp;
    
    temp = prevRowNode->right->info + elem;
    if (floatIsZero(temp))
    {
      /* element exists and the sum is zero: removes the element */
      prevColumnNode = findPrevColumnNodeGivenNode(prevRowNode->right);
      removeNodeFromMatrix(prevColumnNode->below, prevRowNode, prevColumnNode);
      return 0;
    }
    
    /* element exists, updates the matrix with the new value */
    prevRowNode->right->info = temp;
    return 0;
  }
}