#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "matrix.h"
#include "math.h"

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
  uint32 index;
  double x;
  double y;
  double u;
};
typedef struct node Node;

struct nodeArray {
  uint32 len;
  Node *ptr;
};
typedef struct nodeArray NodeArray;

struct element {
  uint32 index;
  uint32 material;
  uint32 n0;
  uint32 n1;
  uint32 n2;
};
typedef struct element Element;

struct elementArray{
  uint32 len;
  Element *ptr;
};
typedef struct elementArray ElementArray;

struct material {
  uint32 index;
  double e;
  double f;
};
typedef struct material Material;

struct dirichletBC {
  uint32 nodeIndex;
  double value;
};
typedef struct dirichletBC DirichletBC;

struct neumannBC {
  uint32 elementIndex;
  int64 nSide;
  double value;
};
typedef struct neumannBC NeumannBC;

int readNodes(NodeArray *nodes, char *fileName);
int readElements(ElementArray *elements, char *fileName);
int readMaterials(uint32 *nMaterial, Material **materials, char *fileName);
int readBC(uint32 *nDirich, DirichletBC **dirBound, uint32 *nNeumann, NeumannBC **neuBound, char *fileName);

int computeGlobalMatrices(Matrix **a, Vector **b, const ElementArray *elements, const NodeArray *nodes, const uint32 nDirichlet, const DirichletBC *dirBound, const uint32 nEquations, const int64 *nodeEquation, const Material *elemMaterial);

int eleN3P3(Matrix *stifMatrix, Vector *b, const Element *element, const NodeArray *nodes, const uint32 nDirichlet, const DirichletBC *dirBound, const int64 *nodeEquation, const Material *elemMaterial);


int conjugateGradientMethod(const MatrixCRS *a, Vector *x, const Vector *b);

static DirichletBC* findElementInArray(const int64 nodeIndex, const int64 nDirichlet, const DirichletBC *dirBound);


int main( int argc, char *argv[] ) 
{
  NodeArray nodes;
  ElementArray elements;
  
  uint32 nMaterial;
  Material *materials;
  
  uint32 nDirichlet;
  DirichletBC *dirBound;
  
  uint32 nNeumann;
  NeumannBC *neuBound;
  
  /* read arguments */
  char *fileName;
  if( argc == 2 ) 
  {
    fileName = argv[1];
  }
  else if( argc > 2 ) {
    fprintf(stderr, "Muitos argumentos submetidos.\n");
    return 1;
  }
  else {
    fprintf(stderr, "Passar o nome do arquivo.\n");
    return 1;
  }
  
  /*
  unsigned nDOF;
  */
  
  /* (a) Input of data */
  /* (b) Representation of the triangulation Th. */
  readNodes(&nodes, fileName);
  readElements(&elements, fileName);
  readBC(&nDirichlet, &dirBound, &nNeumann, &neuBound, fileName);
  readMaterials(&nMaterial, &materials, fileName);
  
  /* find numbers of equations and builds the transformation array */
  uint32 nEquations = (uint32)(nodes.len - nDirichlet);
  int64 *nodeEquation = malloc(nodes.len * sizeof(int64)); /* in each position will be the position of the corresponding equation */
  
  uint32 i = 0;
  uint32 nodeIdx = 0;
  while(i < nodes.len)
  {
    const DirichletBC *found = findElementInArray(i, nDirichlet, dirBound);
    if (found == NULL)
    {
      nodeEquation[i] = nodeIdx;
      nodeIdx++;
    }
    else
    {
      nodeEquation[i] = -1;
    }
    i++;
  }
  
  /* (c) Computation of the element stiffness matrices aK and element loads bK. */
  /* (d) Assembly of the global stiffness matrix A and load vector b. */
  Matrix *a;
  Vector *b;
  computeGlobalMatrices(&a, &b, &elements, &nodes, nDirichlet, dirBound, nEquations, nodeEquation, materials);
  
  /* (e) Solution of the system of equations Ax = b. */
  /* TODO: calc (f,b) */
  MatrixCRS *aCrs; 
  matrixcrs_create(a, &aCrs);
  Vector *x;
  vector_create(&x, nEquations);
  conjugateGradientMethod(aCrs, x, b);
  
  /* (f) Presentation of result. */
  vector_print(x);

  return 0;
}

int readNodes(NodeArray *nodes, char *fileName)
{
  uint32 n;

  FILE *fp;
  char fullPath[80];
  snprintf(fullPath, sizeof(fullPath), "%s.node", fileName);
  fp = fopen(fullPath, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "Cannot open file %s.\n", fullPath);
    return 1;
  }
  
  fscanf(fp, "%" SCNu32 " %*[^\n]", &n);
  if (n < 1)
  {
    fprintf(stderr, "No nodes to be read\n");
    fclose(fp);
    return 1;
  }
  
  /* There are nodes, allocate memory and reads them */
  nodes->ptr = malloc(n * sizeof(Node));
  uint32 i;
  double x, y;
  uint32 index;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNu32" %lf %lf %*[^\n]\n", &index, &x, &y);
    (nodes->ptr + i)->index = index;
    (nodes->ptr + i)->x = x;
    (nodes->ptr + i)->y= y;
  }
  nodes->len = n;
  
  fclose(fp);
  return 0;
}

int readElements(ElementArray *elements, char *fileName)
{
  uint32 n;

  FILE *fp;
  char fullPath[80];
  snprintf(fullPath, sizeof(fullPath), "%s.ele", fileName);
  fp = fopen(fullPath, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "Cannot open file %s.\n", fullPath);
    return 1;
  }
  
  fscanf(fp, "%" SCNu32 " %*[^\n]", &n);
  if (n < 1)
  {
    fprintf(stderr, "No elements to be read\n");
    fclose(fp);
    return 1;
  }
  
  /* There are nodes, allocate memory and reads them */
  elements->ptr = malloc(n * sizeof(Element));
  uint32 i;
  uint32 index, n0, n1, n2;
  uint32 matIdx;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNu32 " %" SCNu32 " %" SCNu32 " %" SCNu32 " %" SCNu32 "%*[^\n]\n", &index, &n0, &n1, &n2, &matIdx);
    (elements->ptr + i)->index = index;
    (elements->ptr + i)->n0 = n0;
    (elements->ptr + i)->n1 = n1;
    (elements->ptr + i)->n2 = n2;
    (elements->ptr + i)->material = matIdx;
  }
  
  elements->len = n;
  
  fclose(fp);
  return 0;
}

int readMaterials(uint32 *nMaterial, Material **materials, char *fileName)
{
  uint32 n;

  FILE *fp;
  char fullPath[80];
  snprintf(fullPath, sizeof(fullPath), "%s.mat", fileName);
  fp = fopen(fullPath, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "Cannot open file %s.\n", fullPath);
    return 1;
  }
  
  fscanf(fp, "%" SCNu32 "%*[^\n]", &n);
  if (n < 1)
  {
    fprintf(stderr, "No materials to be read\n");
    fclose(fp);
    return 1;
  }
  
  /* There are materials, allocate memory and reads them */
  Material *newMaterials = malloc(n * sizeof(Material));
  uint32 i;
  double e, density;
  uint32 index;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNu32" %lf %lf*\n", &index, &e, &density);
    (*(newMaterials + i)).index = index;
    (*(newMaterials + i)).e = e;
    (*(newMaterials + i)).f = density * 9.80665;
  }
  
  *nMaterial = n;
  *materials = newMaterials;
  
  fclose(fp);
  
  return 0;
}

int readBC(uint32 *nDirich, DirichletBC **dirBound, uint32 *nNeumann, NeumannBC **neuBound, char *fileName)
{
  uint32 n;

  FILE *fp;
  char fullPath[80];
  snprintf(fullPath, sizeof(fullPath), "%s.bc", fileName);
  fp = fopen(fullPath, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "Cannot open file %s.\n", fullPath);
    return 1;
  }
  
  fscanf(fp, "%" SCNu32 "*\n", &n);
  if (n < 1)
  {
    fprintf(stderr, "No dirichlet boundary conditions to be read\n");
    fclose(fp);
    return 1;
  }
  
  /* There are bcs, allocate memory and reads them */
  DirichletBC *dirBC = malloc(n * sizeof(DirichletBC));
  
  uint32 i;
  uint32 nodeIndex;
  double value;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNu32" %lf*\n", &nodeIndex, &value);
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

int computeGlobalMatrices(Matrix **a, Vector **b, const ElementArray *elements, const NodeArray *nodes, const uint32 nDirichlet, const DirichletBC *dirBound, const uint32 nEquations, const int64 *nodeEquation, const Material *elemMaterial)
{
  Matrix *newStif;
  matrix_create(&newStif, nEquations, nEquations);
  Vector *newB;
  vector_create(&newB, nEquations);
  
  uint32 i;
  for (i = 0; i < elements->len; i++)
  {
    eleN3P3(newStif, newB, &elements->ptr[i], nodes, nDirichlet, dirBound, nodeEquation, &elemMaterial[elements->ptr[i].material]);
  }
  
  *a = newStif;
  *b = newB;
  return 0;
}

int eleN3P3(Matrix *stifMatrix, Vector *b, const Element *element, const NodeArray *nodes, const uint32 nDirichlet, const DirichletBC *dirBound, const int64 *nodeEquation, const Material *elemMaterial)
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
  n0 = &nodes->ptr[element->n0];
  n1 = &nodes->ptr[element->n1];
  n2 = &nodes->ptr[element->n2];
  
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
  
  int conditionNumber = 0;
  if (eqNode0 == -1) conditionNumber += 1; /* node 0 has prescribed value */
  if (eqNode1 == -1) conditionNumber += 2; /* node 1 has prescribed value */
  if (eqNode2 == -1) conditionNumber += 4; /* node 2 has prescribed value */
  
  DirichletBC *dirNode0, *dirNode1, *dirNode2;
  double fb = elemMaterial->f/elemMaterial->e * area;
  switch (conditionNumber)
  {
    case 0: /* no dirichlet nodes */
      matrix_sumelem(stifMatrix, eqNode0, eqNode0, k00); /* K00 */
      matrix_sumelem(stifMatrix, eqNode0, eqNode1, k01); /* K01 */
      matrix_sumelem(stifMatrix, eqNode0, eqNode2, k02); /* K02 */
      matrix_sumelem(stifMatrix, eqNode1, eqNode0, k01); /* K10 */
      matrix_sumelem(stifMatrix, eqNode1, eqNode1, k11); /* K11 */
      matrix_sumelem(stifMatrix, eqNode1, eqNode2, k12); /* K12 */
      matrix_sumelem(stifMatrix, eqNode2, eqNode0, k02); /* K20 */
      matrix_sumelem(stifMatrix, eqNode2, eqNode1, k12); /* K21 */
      matrix_sumelem(stifMatrix, eqNode2, eqNode2, k22); /* K22 */
      vector_sumelem(b, (const uint32)eqNode0, fb);
      vector_sumelem(b, (const uint32)eqNode1, fb);
      vector_sumelem(b, (const uint32)eqNode2, fb);
      break;
      
    case 1: /* node 0 is a dirichlet node */
      dirNode0 = findElementInArray(element->n0, nDirichlet, dirBound);
    
      matrix_sumelem(stifMatrix, eqNode1, eqNode1, k11); /* K11 */
      matrix_sumelem(stifMatrix, eqNode1, eqNode2, k12); /* K12 */
      matrix_sumelem(stifMatrix, eqNode2, eqNode1, k12); /* K21 */
      matrix_sumelem(stifMatrix, eqNode2, eqNode2, k22); /* K22 */
      
      vector_sumelem(b, (const uint32)eqNode1, fb - k01*dirNode0->value);
      vector_sumelem(b, (const uint32)eqNode2, fb - k02*dirNode0->value);
      break;
      
    case 2: /* node 1 is a dirichlet node */
      dirNode1 = findElementInArray(element->n1, nDirichlet, dirBound);
    
      matrix_sumelem(stifMatrix, eqNode0, eqNode0, k00); /* K00 */
      matrix_sumelem(stifMatrix, eqNode0, eqNode2, k02); /* K02 */
      matrix_sumelem(stifMatrix, eqNode2, eqNode0, k02); /* K20 */
      matrix_sumelem(stifMatrix, eqNode2, eqNode2, k22); /* K22 */
      vector_sumelem(b, (const uint32)eqNode0, fb - k01*dirNode1->value);
      vector_sumelem(b, (const uint32)eqNode2, fb - k12*dirNode1->value);
      break;
      
    case 3: /* nodes 0 and 1 are dirichlet nodes */
      dirNode0 = findElementInArray(element->n0, nDirichlet, dirBound);
      dirNode1 = findElementInArray(element->n1, nDirichlet, dirBound);
    
      matrix_sumelem(stifMatrix, eqNode2, eqNode2, k22); /* K22 */
      vector_sumelem(b, (const uint32)eqNode2, fb - k02*dirNode0->value - k12*dirNode1->value);
      break;
      
    case 4: /* node 2 is a dirichlet node */
      dirNode2 = findElementInArray(element->n2, nDirichlet, dirBound);
      
      matrix_sumelem(stifMatrix, eqNode0, eqNode0, k00); /* K00 */
      matrix_sumelem(stifMatrix, eqNode0, eqNode1, k01); /* K01 */
      matrix_sumelem(stifMatrix, eqNode1, eqNode0, k01); /* K10 */
      matrix_sumelem(stifMatrix, eqNode1, eqNode1, k11); /* K11 */
      vector_sumelem(b, (const uint32)eqNode0, fb - k02*dirNode2->value);
      vector_sumelem(b, (const uint32)eqNode1, fb - k12*dirNode2->value);
      break;
      
    case 5: /* nodes 0 and 2 are dirichlet nodes */
      dirNode0 = findElementInArray(element->n0, nDirichlet, dirBound);
      dirNode2 = findElementInArray(element->n2, nDirichlet, dirBound);
    
      matrix_sumelem(stifMatrix, eqNode1, eqNode1, k11); /* K11 */
      vector_sumelem(b, (const uint32)eqNode1, fb - k01*dirNode0->value - k12*dirNode2->value);
      break;
      
    case 6: /* nodes 1 and 2 are dirichlet nodes */
      dirNode1 = findElementInArray(element->n1, nDirichlet, dirBound);
      dirNode2 = findElementInArray(element->n2, nDirichlet, dirBound);
    
      matrix_sumelem(stifMatrix, eqNode0, eqNode0, k00); /* K00 */
      vector_sumelem(b, (const uint32)eqNode0, fb - k01*dirNode1->value - k02*dirNode2->value);
      break;
      
    case 7: /* all nodes are dirichlet nodes */
      /* do nothing */
      break;
    default:
      break;
  }
  
  return 0;
}

int conjugateGradientMethod(const MatrixCRS *a, Vector *x, const Vector *b)
{
  uint32 n = vector_size(x);

  double c, t, d;
  Vector *p, *r, *temp;
  
  vector_create(&p, n);
  vector_create(&r, n);
  vector_create(&temp, n);
  
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
  matrixcrs_multiplyVector(a, x, temp);
  vector_subtract(b, temp, r);
  
  /* p = r */
  vector_copy(r, p);
  
  /* c = (r,r) */
  c = vector_internalProduct(r, r);
  
  /* for (i = 0 to M) do */
  uint32 i;
  for (i = 0; i < n; i++)
  {
    /* if (p,p)^0.5 < tolA then exit loop; */
    if (sqrt(vector_internalProduct(p, p)) < TOL) break;
    
    /* z = Ap */
    matrixcrs_multiplyVector(a, p, temp);
    
    /* t = c / (p,z) */
    t = c / vector_internalProduct(p, temp);
    
    /* x = x + tp */
    vector_addVectorScalar(x, p, t, x);
    
    /* r = r - tz */
    vector_addVectorScalar(r, temp, -t, r);
    
    /* d = (r,r) */
    d = vector_internalProduct(r, r);
    
    /*   if (d < tolB) then exit loop */
    if (sqrt(d) < TOL) break;
    
    /* p = r + (d/c)p */
    vector_addVectorScalar(r, p, d/c, p);
    
    /* c = d */
    c = d;
  }
  
  vector_destroy(p);
  vector_destroy(r);
  vector_destroy(temp);
  
  return 0;
}

static DirichletBC* findElementInArray(const int64 nodeIndex, const int64 nDirichlet, const DirichletBC *dirBound)
{
  int i;
  for (i = 0; i < nDirichlet; i++)
  {
    if (dirBound[i].nodeIndex == nodeIndex)
      return (DirichletBC*)(dirBound+i);
  }
  return NULL;
}