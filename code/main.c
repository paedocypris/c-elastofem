#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "matrix.h"

typedef int8_t int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;

typedef uint8_t uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

struct node {
  int64 index;
  double x;
  double y;
};

typedef struct node Node;

struct element {
  int64 index;
  int64 n0;
  int64 n1;
  int64 n2;
};

typedef struct element Element;

int createNodes(int64 *nNode, Node **nodes);
int createElements(int64 *nElem, Element **elements);
int computeGlobalMatrices(const int64 nElem, const Element *elementsArray, const int64 nNode, const Node *nodesArray);
int computeElementStiffness(Matrix *stifMatrix, const Element *element, const Node *nodes);

int conjugateGradientMethod(int64 ni, int64 nj, Matrix *a, Matrix *x, Matrix *b);


int main( void ) {
  int64 nNode;
  Node *nodes;
  
  int64 nElement;
  Element *elements;
  /*
  unsigned nDOF;
  */
  
  /* (b) Representation of the triangulation Th. */
  createNodes(&nNode, &nodes);
  createElements(&nElement, &elements);
  
  /* (c) Computation of the element stiffness matrices aK and element loads bK. */
  /* (d) Assembly of the global stiffness matrix A and load vector b. */
  computeGlobalMatrices( nElement, elements, nNode, nodes);
  
  /* (e) Solution of the system of equations Ax = b. */
  /* TODO: calc (f,b) */
  Matrix *b;
  matrix_create(&b, nNode, 1);
  Matrix *x;
  matrix_create(&x, nNode, 1);
  
  /* (f) Presentation of result. */

  return 0;
}

int createNodes(int64 *nNode, Node **nodes)
{
  int64 n;

  FILE *fp;
  fp = fopen("bar.node", "r");
  fscanf(fp, "%" SCNd64 " %*[^\n]", &n);
  if (n < 1)
  {
    fprintf(stderr, "No nodes to be read\n");
    fclose(fp);
    return 1;
  }
  
  /* There are nodes, allocate memory and reads them */
  Node *newNodes = malloc(n * sizeof(Node));
  int64 i;
  double x, y;
  int64 index;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNd64" %lf %lf %*[^\n]\n", &index, &x, &y);
    (*(newNodes + i)).index = index;
    (*(newNodes + i)).x = x;
    (*(newNodes + i)).y = y;
  }
  
  *nodes = newNodes;
  *nNode = n;
  
  fclose(fp);
  return 0;
}

int createElements(int64 *nElem, Element **elements)
{
  int64 n;

  FILE *fp;
  fp = fopen("bar.ele", "r");
  fscanf(fp, "%" SCNd64 " %*[^\n]", &n);
  if (n < 1)
  {
    fprintf(stderr, "No elements to be read\n");
    fclose(fp);
    return 1;
  }
  
  /* There are nodes, allocate memory and reads them */
  Element *newElements = malloc(n * sizeof(Element));
  int64 i;
  int64 index, n0, n1, n2;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNd64 " %" SCNd64 " %" SCNd64 " %" SCNd64 "%*[^\n]\n", &index, &n0, &n1, &n2);
    (*(newElements + i)).index = index;
    (*(newElements + i)).n0 = n0;
    (*(newElements + i)).n1 = n1;
    (*(newElements + i)).n2 = n2;
  }
  
  *elements = newElements;
  *nElem = n;
  
  fclose(fp);
  return 0;
}

int computeGlobalMatrices(const int64 nElem, const Element *elementsArray, const int64 nNode, const Node *nodesArray)
{
  Matrix *stifMatrix;
  matrix_create(&stifMatrix, nNode, nNode);
  
  int i;
  for (i = 0; i < nElem; i++)
  {
    computeElementStiffness(stifMatrix, &elementsArray[i],nodesArray);
    matrix_print(stifMatrix);
  }
  return 0;
}

int computeElementStiffness(Matrix *stifMatrix, const Element *element, const Node *nodes)
{
  double a0, a1, a2;
  double b0, b1, b2;
  double det;
  
  double k00, k01, k02, k11, k12, k22;
  
  const Node *n0, *n1, *n2;
  n0 = &nodes[element->n0];
  n1 = &nodes[element->n1];
  n2 = &nodes[element->n2];
  
  /* calc auxiliar parameters */
  det = n0->x*n1->y + n0->y*n2->x + n1->x*n2->y - n1->y*n2->x - n2->y*n0->x - n0->y*n1->x;
  a0 = n1->y - n2->y;
  a1 = n2->y - n0->y;
  a2 = n0->y - n1->y;
  b0 = n2->x - n1->x;
  b1 = n0->x - n2->x;
  b2 = n1->x - n0->x;
  
  double den = det*2;
  k00 = (a0*a0 + b0*b0)/den;
  k01 = (a0*a1 + b0*b1)/den;
  k02 = (a0*a2 + b0*b2)/den;
  k11 = (a1*a1 + b1*b1)/den;
  k12 = (a1*a2 + b1*b2)/den;
  k22 = (a2*a2 + b2*b2)/den;
  
  matrix_sumelem(stifMatrix, element->n0, element->n0, k00); /* K00 */
  matrix_sumelem(stifMatrix, element->n0, element->n1, k01); /* K01 */
  matrix_sumelem(stifMatrix, element->n0, element->n2, k02); /* K02 */
  matrix_sumelem(stifMatrix, element->n1, element->n0, k01); /* K10 */
  matrix_sumelem(stifMatrix, element->n1, element->n1, k11); /* K11 */
  matrix_sumelem(stifMatrix, element->n1, element->n2, k12); /* K12 */
  matrix_sumelem(stifMatrix, element->n2, element->n0, k02); /* K20 */
  matrix_sumelem(stifMatrix, element->n2, element->n1, k12); /* K21 */
  matrix_sumelem(stifMatrix, element->n2, element->n2, k22); /* K22 */
  
  return 0;
}

int conjugateGradientMethod(int64 ni, int64 nj, Matrix *a, Matrix *x, Matrix *b)
{
  // double c, t, d;
  Matrix *p, *r, *temp;
  matrix_create(&p, nj, 1);
  matrix_create(&r, nj, 1);
  matrix_create(&temp, nj, 1);
  
  /* r = b - Ax;
     p = r;
     c = (r,r);
     for (k = 0 to M) do
       if (p,p)^0.5 < tolA then exit loop;
       z = Ap;
       t = c / (p,z);
       x = x + tp;
       r = r - tz;
       d = (r,r)
       if (d < tolB) then exit loop;
       p = r + (d/c)p;
       c = d;
     end do */
     
   /* r = b - Ax 
  matrix_multiply(a, x, temp);
  matrix_subtract(b, temp, r);
  
  multMatrix (ni, nj, a, nj, 1, x, temp);
  subtractVector(nj, b, temp, r);
  
  /* p = r 
  copyVector(nj, r, p);
  
  /* c = (r,r) 
  c = internalProduct(nj, r, r);
  
  int i;
  for (i = 0; i < M; i++)
  {
    if (internalProduct(nj, p, p) < TOL) break;
    
    /* z = Ap 
    multMatrix(ni, nj, a, nj, 1, p, temp);
    
    /* t = c / (p,z) 
    t = c / internalProduct(nj, p, temp);
    
    /* x = x + tp; 
    addVectorScalar(nj, x, p, t, x);
    
    /* r = r - tz; 
    addVectorScalar(nj, r, temp, -t, r);
    
    /* d = (r,r) 
    d = internalProduct(nj, r, r);
    
    /*   if (d < tolB) then exit loop; 
    if (d < TOL) break;
    
    /* p = r + (d/c)p; 
    addVectorScalar(nj, r, p, d/c, p);
    
    /* c = d; 
    c = d;
  }
  TODO: DEstroy this   */
  
  matrix_destroy(p);
  matrix_destroy(r);
  matrix_destroy(temp);
  
  return 0;
}