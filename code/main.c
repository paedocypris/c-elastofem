#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "matrix.h"

#define TOL 1e-6

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
  double u;
};
typedef struct node Node;

struct element {
  int64 index;
  int64 material;
  int64 n0;
  int64 n1;
  int64 n2;
};
typedef struct element Element;

struct material {
  int64 index;
  double e;
  double f;
};
typedef struct material Material;

struct dirichletBC {
  int64 nodeIndex;
  double value;
};
typedef struct dirichletBC DirichletBC;

struct neumannBC {
  int64 elementIndex;
  int64 nSide;
  double value;
};
typedef struct neumannBC NeumannBC;

int readNodes(int64 *nNode, Node **nodes);
int readElements(int64 *nElem, Element **elements);
int readMaterials(int64 *nMaterial, Material **materials);
int readBC(int64 *nDirich, DirichletBC **dirBound, int64 *nNeumann, NeumannBC **neuBound);

int computeGlobalMatrices(const int64 nElem, const Element *elementsArray, const int64 nNode, const Node *nodesArray, const int64 nDirichlet, const DirichletBC *dirBound, const int64 nEquations, const int64 *nodeEquation, const Material *elemMaterial);
int eleN3P3(Matrix *stifMatrix, Matrix *b, const Element *element, const Node *nodes, const int64 nDirichlet, const DirichletBC *dirBound, const int64 *nodeEquation, const Material *elemMaterial);

int conjugateGradientMethod(const Matrix *a, const Matrix *x, const Matrix *b, Matrix **result);

static const DirichletBC* findElementInArray(const int64 node, const int64 nDirichlet, const DirichletBC *dirBound);


int main( void ) {
  int64 nNode;
  Node *nodes;
  
  int64 nElement;
  Element *elements;
  
  int64 nMaterial;
  Material *materials;
  
  int64 nDirichlet;
  DirichletBC *dirBound;
  
  int64 nNeumann;
  NeumannBC *neuBound;
  
  /*
  unsigned nDOF;
  */
  
  /* (a) Input of data */
  /* (b) Representation of the triangulation Th. */
  readNodes(&nNode, &nodes);
  readElements(&nElement, &elements);
  readBC(&nDirichlet, &dirBound, &nNeumann, &neuBound);
  readMaterials(&nMaterial, &materials);
  
  /* find numbers of equations and builds the transformation array */
  int64 nEquations;
  int64 *nodeEquation; /* in each position will be the position of the corresponding equation */
  nEquations = nNode - nDirichlet;
  nodeEquation = malloc(nNode*sizeof(int64));
  
  /* (c) Computation of the element stiffness matrices aK and element loads bK. */
  /* (d) Assembly of the global stiffness matrix A and load vector b. */
  computeGlobalMatrices( nElement, elements, nNode, nodes, nDirichlet, dirBound, nEquations, nodeEquation, materials);
  
  /* (e) Solution of the system of equations Ax = b. */
  /* TODO: calc (f,b) */
  
  Matrix *x;
  matrix_create(&x, nNode, 1);
  
  /* (f) Presentation of result. */

  return 0;
}

int readNodes(int64 *nNode, Node **nodes)
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

int readElements(int64 *nElem, Element **elements)
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

int readMaterials(int64 *nMaterial, Material **materials)
{
  int64 n;

  FILE *fp;
  fp = fopen("bar.mat", "r");
  fscanf(fp, "%" SCNd64 "%*[^\n]", &n);
  if (n < 1)
  {
    fprintf(stderr, "No materials to be read\n");
    fclose(fp);
    return 1;
  }
  
  /* There are materials, allocate memory and reads them */
  Material *newMaterials = malloc(n * sizeof(Material));
  int64 i;
  double e, density;
  int64 index;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNd64" %lf %lf*\n", &index, &e, &density);
    (*(newMaterials + i)).index = index;
    (*(newMaterials + i)).e = e;
    (*(newMaterials + i)).f = density * 9.80665;
  }
  
  *nMaterial = n;
  *materials = newMaterials;
  
  fclose(fp);
  
  return 0;
}

int readBC(int64 *nDirich, DirichletBC **dirBound, int64 *nNeumann, NeumannBC **neuBound)
{
  int64 n;

  FILE *fp;
  fp = fopen("bar.bc", "r");
  fscanf(fp, "%" SCNd64 "*\n", &n);
  if (n < 1)
  {
    fprintf(stderr, "No dirichlet boundary conditions to be read\n");
    fclose(fp);
    return 1;
  }
  
  /* There are bcs, allocate memory and reads them */
  DirichletBC *dirBC = malloc(n * sizeof(DirichletBC));
  
  int64 i;
  int64 nodeIndex;
  double value;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNd64" %lf*\n", &nodeIndex, &value);
    (*(dirBC + i)).nodeIndex = nodeIndex;
    (*(dirBC + i)).value = value;
  }
  
  *nDirich = n;
  *dirBound = dirBC;
  
  *nNeumann = 0;
  *neuBound = NULL;
  
  fclose(fp);
  
  return 0;
}

int computeGlobalMatrices(const int64 nElem, const Element *elementsArray, const int64 nNode, const Node *nodesArray, const int64 nDirichlet, const DirichletBC *dirBound, const int64 nEquations, const int64 *nodeEquation, const Material *elemMaterial)
{
  Matrix *stifMatrix;
  matrix_create(&stifMatrix, nEquations, nEquations);
  
  Matrix *bMatrix;
  matrix_create(&bMatrix, nEquations, 1);
  
  int i;
  for (i = 0; i < nElem; i++)
  {
    eleN3P3(stifMatrix, bMatrix, &elementsArray[i], nodesArray, nDirichlet, dirBound, nodeEquation, &elemMaterial[elementsArray[i].material]);
    
    matrix_print(stifMatrix);
  }
  return 0;
}

int eleN3P3(Matrix *stifMatrix, Matrix *b, const Element *element, const Node *nodes, const int64 nDirichlet, const DirichletBC *dirBound, const int64 *nodeEquation, const Material *elemMaterial)
{
  double a0, a1, a2;
  double b0, b1, b2;
  double det, den, area;

/*
  int64 nNeumann;
  NeumannBC *neuBound;
*/

  double k00, k01, k02, k11, k12, k22;
  
  const Node *n0, *n1, *n2;
  n0 = &nodes[element->n0];
  n1 = &nodes[element->n1];
  n2 = &nodes[element->n2];
  
  const int64 eqNode0 = nodeEquation[element->n0];
  const int64 eqNode1 = nodeEquation[element->n1];
  const int64 eqNode2 = nodeEquation[element->n2];
  
  /* calc auxiliar parameters */
  det = n0->x*n1->y + n0->y*n2->x + n1->x*n2->y - n1->y*n2->x - n2->y*n0->x - n0->y*n1->x;
  den = det*2;
  area = det/2;
  a0 = n1->y - n2->y;
  a1 = n2->y - n0->y;
  a2 = n0->y - n1->y;
  b0 = n2->x - n1->x;
  b1 = n0->x - n2->x;
  b2 = n1->x - n0->x;
  
  k00 = (a0*a0 + b0*b0)/den;
  k01 = (a0*a1 + b0*b1)/den;
  k02 = (a0*a2 + b0*b2)/den;
  k11 = (a1*a1 + b1*b1)/den;
  k12 = (a1*a2 + b1*b2)/den;
  k22 = (a2*a2 + b2*b2)/den;
  
  const DirichletBC *dirNode0, *dirNode1, *dirNode2;
  dirNode0 = findElementInArray(element->n0, nDirichlet, dirBound);
  dirNode1 = findElementInArray(element->n1, nDirichlet, dirBound);
  dirNode2 = findElementInArray(element->n2, nDirichlet, dirBound);
  
  if (dirNode0 == NULL)
  {
    matrix_sumelem(stifMatrix, eqNode0, eqNode0, k00); /* K00 */
    matrix_sumelem(stifMatrix, eqNode0, eqNode1, k01); /* K01 */
    matrix_sumelem(stifMatrix, eqNode0, eqNode2, k02); /* K02 */
    matrix_sumelem(b, eqNode0, 0, elemMaterial->f/elemMaterial->e *area);
  }
  if (dirNode1 == NULL)
  {
    matrix_sumelem(stifMatrix, eqNode1, eqNode0, k01); /* K10 */
    matrix_sumelem(stifMatrix, eqNode1, eqNode1, k11); /* K11 */
    matrix_sumelem(stifMatrix, eqNode1, eqNode2, k12); /* K12 */
    matrix_sumelem(b, eqNode1, 0, elemMaterial->f/elemMaterial->e *area);
  }
  if (dirNode2 == NULL)
  {
    matrix_sumelem(stifMatrix, eqNode2, eqNode0, k02); /* K20 */
    matrix_sumelem(stifMatrix, eqNode2, eqNode1, k12); /* K21 */
    matrix_sumelem(stifMatrix, eqNode2, eqNode2, k22); /* K22 */
    matrix_sumelem(b, eqNode2, 0, elemMaterial->f/elemMaterial->e *area);
  }
  
  int conditionNumber = 0;
  if (dirNode0 != NULL) conditionNumber += 1; /* node 0 has prescribed value */
  if (dirNode1 != NULL) conditionNumber += 2; /* node 1 has prescribed value */
  if (dirNode2 != NULL) conditionNumber += 4; /* node 2 has prescribed value */
  
  switch (conditionNumber)
  {
    case 0: /* no dirichlet nodes */
      break;
    case 1: /* node 0 is a dirichlet node */
      matrix_sumelem(b, eqNode1, 0, -k01*dirNode1->value);
      matrix_sumelem(b, eqNode2, 0, -k02*dirNode2->value);
      break;
    case 2: /* node 1 is a dirichlet node */
      matrix_sumelem(b, eqNode0, 0, -k01*dirNode0->value);
      matrix_sumelem(b, eqNode2, 0, -k12*dirNode2->value);
      break;
    case 3: /* nodes 0 and 1 are dirichlet nodes */
      matrix_sumelem(b, eqNode2, 0, -k02*dirNode0->value - k12*dirNode1->value);
      break;
    case 4: /* node 2 is a dirichlet node */
      matrix_sumelem(b, eqNode0, 0, -k02*dirNode0->value);
      matrix_sumelem(b, eqNode1, 0, -k12*dirNode2->value);
      break;
    case 5: /* nodes 0 and 2 are dirichlet nodes */
      matrix_sumelem(b, eqNode1, 0, -k01*dirNode0->value - k12*dirNode2->value);
      break;
    case 6: /* nodes 1 and 2 are dirichlet nodes */
      matrix_sumelem(b, eqNode0, 0, -k01*dirNode1->value - k02*dirNode2->value);
      break;
    case 7: /* all nodes are dirichlet nodes */
      break;
    default:
      break;
  }
  
  return 0;
}

int conjugateGradientMethod(const Matrix *a, const Matrix *x, const Matrix *b, Matrix **result)
{
  double c, t, d;
  Matrix  *p = NULL, 
          *r = NULL, 
          *temp = NULL, 
          *z = NULL, 
          *nextVector = NULL, 
          *prevVector = NULL;
  
  Matrix *curX;
  matrix_copy(x, &curX);
  
  /* r = b - Ax;
     p = r;
     c = (r,r);
     for (i = 0 to M) do
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
     
   /* r = b - Ax */
  matrix_multiply(a, curX, &temp);
  matrix_subtract(b, temp, &r);
  
  matrix_destroy(temp);
  temp = NULL;
  
  /* p = r */
  matrix_copy(r, &p);
  
  /* c = (r,r) */
  c = matrix_internalProduct(r, r);
  
  /* for (i = 0 to M) do */
  int i;
  for (i = 0; i < matrix_ni(b); i++)
  {
    /* if (p,p)^0.5 < tolA then exit loop; */
    if (matrix_internalProduct(p, p) < TOL) break;
    
    /* z = Ap */
    matrix_multiply(a, p, &z);
    
    /* t = c / (p,z) */
    t = c / matrix_internalProduct(p, z);
    
    /* x = x + tp */
    matrix_scalarmult(t, p, &temp);
    matrix_add(curX, temp, &nextVector);
    prevVector = curX;
    curX = nextVector;
    nextVector=NULL;
    
    matrix_destroy(prevVector);
    prevVector=NULL;
    matrix_destroy(temp);
    temp = NULL;
    
    /* r = r - tz */
    matrix_scalarmult(t,z, &temp);
    matrix_subtract(r, temp, &nextVector);
    prevVector = r;
    r = nextVector;
    nextVector=NULL;
    
    matrix_destroy(prevVector);
    prevVector=NULL;
    matrix_destroy(temp);
    temp = NULL;
    matrix_destroy(z);
    z = NULL;
    
    /* d = (r,r) */
    d = matrix_internalProduct(r, r);
    
    /*   if (d < tolB) then exit loop */
    if (d < TOL) break;
    
    /* p = r + (d/c)p */
    matrix_scalarmult(d/c, p, &temp);
    matrix_add(r, temp, &nextVector);
    prevVector = p;
    p = nextVector;
    nextVector = NULL;
    
    matrix_destroy(prevVector);
    prevVector=NULL;
    matrix_destroy(temp);
    temp = NULL;
    
    /* c = d */
    c = d;
  }
  
  if (p != NULL || r != NULL || 
        temp != NULL || z != NULL || 
        prevVector != NULL || nextVector != NULL)
  {
    fprintf(stderr, "Memory Leak: conjugateGradientMethod.\n");
    return -1;
  }
  *result = curX;
  
  return 0;
}

static const DirichletBC* findElementInArray(const int64 node, const int64 nDirichlet, const DirichletBC *dirBound)
{
  int i;
  for (i = 0; i < nDirichlet; i++)
  {
    if (dirBound[i].nodeIndex == node)
      return dirBound+i;
  }
  return NULL;
}