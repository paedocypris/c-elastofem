#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include "matrix.h"
#include "stdMatrix.h"

#define TOL 1e-6

struct node {
  uint32_t index;
  double x;
  double y;
  double u;
  double ui[2];
};
typedef struct node Node;

struct nodeArray {
  uint32_t len;
  Node *ptr;
};
typedef struct nodeArray NodeArray;

struct material {
  uint32_t index;
  double k;
  double fp;
  double l;
  double mu;
  double fg;
};
typedef struct material Material;

struct element {
  uint32_t index;
  Material *material;
  Node *n[3];
  double vd;
  double sEf[2][2];
};
typedef struct element Element;

struct elementArray{
  uint32_t len;
  Element *ptr;
};
typedef struct elementArray ElementArray;

struct dirichletBC {
  uint32_t nodeIndex;
  uint32_t dir;
  double val;
};
typedef struct dirichletBC DirichletBC;

struct neumannBC {
  uint32_t n0;
  uint32_t n1;
  double val[2];
};
typedef struct neumannBC NeumannBC;

int readMaterials(uint32_t *nMaterial, Material **materials, char *fileName);
int readNodes(NodeArray *nodes, char *fileName);
int readElements(ElementArray *elements, const NodeArray *nodes, const Material *materials, char *fileName);
int readBCP(uint32_t *nDirich, DirichletBC **dirBound, uint32_t *nNeumann, NeumannBC **neuBound, char *fileName);
int readBCE(uint32_t *nDirich, DirichletBC **dirBound, uint32_t *nNeumann, NeumannBC **neuBound, char *fileName);

int computeGlobalMatrices(Matrix **a, Vector **b, const ElementArray *elements, const NodeArray *nodes, const uint32_t nDirichlet, const DirichletBC *dirBound, const uint32_t nNeumann, const NeumannBC *neuBound, const uint32_t ndof);

int eleP1Poisson(Matrix *stifMatrix, Vector *b, const Element *element);
int eleP1Elast(Matrix *stifMatrix, Vector *b, const Element *element);
int computeStresses(const ElementArray *elements);

int loadPhiLFunctions(StdMatrix *phiL, Element *element);

int applyNeuBoundary1(Vector *b, const uint32_t nNeumann, const NeumannBC *neuBound, const NodeArray *nodes);
int applyNeuBoundary2(Vector *b, const uint32_t nNeumann, const NeumannBC *neuBound, const NodeArray *nodes);
int applyDirBoundary(Matrix *stifMatrix, Vector *b, const uint32_t nDirichlet, const DirichletBC *dirBound, const uint32_t ndof);

int conjugateGradientMethod(const MatrixCRS *a, Vector *x, const Vector *b);

int loadResultIntoNodePoisson(NodeArray *nodes, Vector *x);
int loadResultIntoNodeElast(NodeArray *nodes, Vector *x);

int printVTK(NodeArray *nodes, ElementArray *elements, char *fileName);

inline uint32_t gIndex(uint32_t nj, uint32_t i, uint32_t j);
inline double getSideLength(Node n0, Node n1);
inline uint32_t assemblyGlobalNumber(const Node *n[3], const uint32_t i, const uint32_t ndof);

double exactx(double x, double y)
{
  return 2*y*y;
}

double exacty(double x, double y)
{
  return 3*x*x;
}

int main( int argc, char *argv[] ) 
{
  NodeArray nodes;
  ElementArray elements;
  
  uint32_t nMaterial;
  Material *materials;
  
  uint32_t nDirichlet;
  DirichletBC *dirBound;
  
  uint32_t nNeumann;
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
  
  /* Poisson problem */
  
  /* (a) Input of data */
  /* (b) Representation of the triangulation Th. */
  readMaterials(&nMaterial, &materials, fileName);
  readNodes(&nodes, fileName);
  readElements(&elements, &nodes, materials, fileName);
  readBCP(&nDirichlet, &dirBound, &nNeumann, &neuBound, fileName);
  
  /* (c) Computation of the element stiffness matrices aK and element loads bK. */
  /* (d) Assembly of the global stiffness matrix A and load vector b. */
  Matrix *a;
  Vector *b;
  computeGlobalMatrices(&a, &b, &elements, &nodes, nDirichlet, 
                            dirBound, nNeumann, neuBound, 1);
  
  /* (e) Solution of the system of equations Ax = b. */
  MatrixCRS *aCrs; 
  matrixcrs_create(a, &aCrs);
  Vector *x;
  vector_createzero(&x, nodes.len);
  conjugateGradientMethod(aCrs, x, b);
  loadResultIntoNodePoisson(&nodes, x);
  
  /* (g) Cleanup. */
  matrix_destroy(a);
  vector_destroy(b);
  free(neuBound);
  free(dirBound);
  matrixcrs_destroy(aCrs);
  vector_destroy(x);
  
  /* Now run the elastic problem */
  
  /* (a) Input of data */
  /* (b) Representation of the triangulation Th. */
  readBCE(&nDirichlet, &dirBound, &nNeumann, &neuBound, fileName);
  
  /* (c) Computation of the element stiffness matrices aK and element loads bK. */
  /* (d) Assembly of the global stiffness matrix A and load vector b. */
  computeGlobalMatrices(&a, &b, &elements, &nodes, nDirichlet, dirBound, nNeumann, neuBound, 2);
  
  /* (e) Solution of the system of equations Ax = b. */
  matrixcrs_create(a, &aCrs);
  vector_createzero(&x, nodes.len*2);
  conjugateGradientMethod(aCrs, x, b);
  loadResultIntoNodeElast(&nodes, x);
  
  /* (f) Calculation of the stresses at elements */
  computeStresses(&elements);
  
  /* (f) Presentation of result. */
  printVTK(&nodes, &elements, fileName);
  /* printSol(&nodes, fileName); */
  
  return 0;
}

int printVTK(NodeArray *nodes, ElementArray *elements, char *fileName)
{
  FILE *fp;
  char fullPath[80];
  snprintf(fullPath, sizeof(fullPath), "%s.vtk", fileName);
  fp = fopen(fullPath, "w");
  if (fp == NULL)
  {
    fprintf(stderr, "Cannot open file %s.\n", fullPath);
    return 1;
  }

  uint32_t i;
  
  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "2D scalar data\n");
  fprintf(fp, "ASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp, "POINTS %"PRIu32" float\n", nodes->len);
  for(i = 0; i < nodes->len; i++)
  {
    fprintf(fp, "%lf %lf 0\n", nodes->ptr[i].x, nodes->ptr[i].y);
  }
  fprintf(fp, "\n");
  fprintf(fp, "CELLS %"PRIu32" %"PRIu32"\n",elements->len, elements->len*4);
  for(i = 0; i < elements->len; i++)
  {
    fprintf(fp, "3 %"PRIu32" %"PRIu32" %"PRIu32"\n", elements->ptr[i].n[0]->index, elements->ptr[i].n[1]->index, elements->ptr[i].n[2]->index);
  }
  fprintf(fp, "\n");
  fprintf(fp, "CELL_TYPES %"PRIu32"\n", elements->len);
  for(i = 0; i < elements->len; i++)
  {
    fprintf(fp, "5\n");
  }
  fprintf(fp, "\n");
  fprintf(fp, "POINT_DATA %"PRIu32"\n", nodes->len);
  
  fprintf(fp, "SCALARS pressure float 1\n");
  fprintf(fp, "LOOKUP_TABLE default\n");
  for(i = 0; i < nodes->len; i++)
  {
    fprintf(fp, "%.6g\n", nodes->ptr[i].u);
  }
  
  fprintf(fp, "VECTORS displacement float\n");
  for(i = 0; i < nodes->len; i++)
  {
    fprintf(fp, "%.6g %.6g 0\n", nodes->ptr[i].ui[0], nodes->ptr[i].ui[1]);
  }
  fprintf(fp, "VECTORS displacementExact float\n");
  for(i = 0; i < nodes->len; i++)
  {
    fprintf(fp, "%.6g %.6g 0\n", exactx(nodes->ptr[i].x,nodes->ptr[i].y), exacty(nodes->ptr[i].x,nodes->ptr[i].y));
  }
  fprintf(fp, "VECTORS error float\n");
  double errorSum = 0;
  for(i = 0; i < nodes->len; i++)
  {
    double errorx = nodes->ptr[i].ui[0] - exactx(nodes->ptr[i].x,nodes->ptr[i].y);
    double errory = nodes->ptr[i].ui[1] - exacty(nodes->ptr[i].x,nodes->ptr[i].y);
    fprintf(fp, "%.6g %.6g 0\n", errorx, errory);
    errorSum += errorx*errorx + errory*errory;
  }
  printf("erro malha %s: %.6g", fileName, sqrt(errorSum/nodes->len/2));
  fprintf(fp, "\n");
  
  fprintf(fp, "CELL_DATA %"PRIu32"\n", elements->len);
  
  fprintf(fp, "TENSORS efStress float\n");
  for(i = 0; i < elements->len; i++)
  {
    fprintf(fp, "%.6g %.6g 0\n", elements->ptr[i].sEf[0][0], elements->ptr[i].sEf[0][1]);
    fprintf(fp, "%.6g %.6g 0\n", elements->ptr[i].sEf[1][0], elements->ptr[i].sEf[1][1]);
    fprintf(fp, "0    0    0\n");
    fprintf(fp, "\n");
  }
  
  fclose(fp);
  return 0;
}

int readNodes(NodeArray *nodes, char *fileName)
{
  uint32_t n;

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
  uint32_t i;
  double x, y;
  uint32_t index;
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

int readElements(ElementArray *elements, const NodeArray *nodes, const Material *materials, char *fileName)
{
  uint32_t n;

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
  uint32_t i;
  uint32_t index, n0, n1, n2;
  uint32_t matIdx;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNu32 " %" SCNu32 " %" SCNu32 " %" SCNu32 " %" SCNu32 "%*[^\n]\n", &index, &n0, &n1, &n2, &matIdx);
    (elements->ptr + i)->index = index;
    (elements->ptr + i)->n[0] = &nodes->ptr[n0];
    (elements->ptr + i)->n[1] = &nodes->ptr[n1];
    (elements->ptr + i)->n[2] = &nodes->ptr[n2];
    (elements->ptr + i)->material = (Material*)&materials[matIdx];
  }
  
  elements->len = n;
  
  fclose(fp);
  return 0;
}

int readMaterials(uint32_t *nMaterial, Material **materials, char *fileName)
{
  uint32_t n;

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
  uint32_t i;
  double k, fPoisson;
  double e, nu, density;
  uint32_t index;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNu32" %lf %lf %lf %lf %lf*\n", &index, &k, &fPoisson, &e, &nu, &density);
    (*(newMaterials + i)).index = index;
    (*(newMaterials + i)).k = k;
    (*(newMaterials + i)).fp = fPoisson;
    (*(newMaterials + i)).l = -e*nu/(nu*(2*nu+1) - 1);
    (*(newMaterials + i)).mu = e/(2*nu + 2);
    (*(newMaterials + i)).fg = density * 9.80665;
  }
  
  *nMaterial = n;
  *materials = newMaterials;
  
  fclose(fp);
  
  return 0;
}

int readBCP(uint32_t *nDirich, DirichletBC **dirBound, uint32_t *nNeumann, NeumannBC **neuBound, char *fileName)
{
  uint32_t n;

  FILE *fp;
  char fullPath[80];
  snprintf(fullPath, sizeof(fullPath), "%s.bcp", fileName);
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
  
  uint32_t i;
  uint32_t bIdx;
  double value;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNu32" %lf*\n", &bIdx, &value);
    (*(dirBC + i)).nodeIndex = bIdx;
    (*(dirBC + i)).dir = 0;
    (*(dirBC + i)).val = value;
  }
  
  *nDirich = n;
  *dirBound = dirBC;
  
  fscanf(fp, " %" SCNu32 "*\n", &n);
  *nNeumann = n;
  if (n > 0)
  {
    NeumannBC *neuBC = malloc(n * sizeof(NeumannBC));
    uint32_t n0, n1;
    for (i = 0; i < n; i++)
    {
      fscanf(fp, "%" SCNu32" %"SCNu32" %lf*\n", &n0, &n1, &value);
      
      (*(neuBC + i)).n0 = n0;
      (*(neuBC + i)).n1 = n1;
      (*(neuBC + i)).val[0] = value;
    }
    *neuBound = neuBC;
  }
  else
  {
    *neuBound = NULL;
  }
  
  fclose(fp);
  
  return 0;
}

int readBCE(uint32_t *nDirich, DirichletBC **dirBound, uint32_t *nNeumann, NeumannBC **neuBound, char *fileName)
{
  uint32_t n;

  FILE *fp;
  char fullPath[80];
  snprintf(fullPath, sizeof(fullPath), "%s.bce", fileName);
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
  
  uint32_t i;
  uint32_t bIdx;
  uint32_t dir;
  double value[2];
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNu32" %"SCNu32" %lf*\n", &bIdx, &dir, &value[0]);
    (*(dirBC + i)).nodeIndex = bIdx;
    (*(dirBC + i)).dir = dir;
    (*(dirBC + i)).val = value[0];
  }
  
  *nDirich = n;
  *dirBound = dirBC;
  
  fscanf(fp, " %" SCNu32 "*\n", &n);
  NeumannBC *neuBC = malloc(n * sizeof(NeumannBC));
  uint32_t n0, n1;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNu32" %"SCNu32" %lf %lf*\n", &n0, &n1, &value[0], &value[1]);
    
    (*(neuBC + i)).n0 = n0;
    (*(neuBC + i)).n1 = n1;
    (*(neuBC + i)).val[0] = value[0];
    (*(neuBC + i)).val[1] = value[1];
  }
  *nNeumann = n;
  *neuBound = neuBC;
  
  fclose(fp);
  
  return 0;
}

int computeGlobalMatrices(Matrix **a, Vector **b, const ElementArray *elements, const NodeArray *nodes, const uint32_t nDirichlet, const DirichletBC *dirBound, const uint32_t nNeumann, const NeumannBC *neuBound, const uint32_t ndof)
{
  Matrix *newStif;
  matrix_create(&newStif, nodes->len*ndof, nodes->len*ndof);
  Vector *newB;
  vector_createzero(&newB, nodes->len*ndof);
  
  uint32_t i;
  for (i = 0; i < elements->len; i++)
  {
    if (ndof == 1)
    {
      eleP1Poisson(newStif, newB, &elements->ptr[i]);
    }
    else
    {
      eleP1Elast(newStif, newB, &elements->ptr[i]);
    }
  }
  
  /* apply neumann boundary */
  if (ndof == 1)
  {
    applyNeuBoundary1(newB, nNeumann, neuBound, nodes);
  }
  else
  {
    applyNeuBoundary2(newB, nNeumann, neuBound, nodes);
  }
  
  /* apply dir boundary */
  applyDirBoundary(newStif, newB, nDirichlet, dirBound, ndof);
  
  *a = newStif;
  *b = newB;
  return 0;
}

int eleP1Poisson(Matrix *stifMatrix, Vector *b, const Element *element)
{ 
  uint32_t i, j;
  Material *elMat = element->material;

  const Node *n[3];
  n[0] = element->n[0];
  n[1] = element->n[1];
  n[2] = element->n[2];
  
  /* calc auxiliar parameters */
  double det, den, area;
  det = n[0]->x*n[1]->y + n[0]->y*n[2]->x + n[1]->x*n[2]->y - n[1]->y*n[2]->x - n[2]->y*n[0]->x - n[0]->y*n[1]->x;
  den = det*2;
  area = det/2;
  
  double phiL[2][3];
  phiL[0][0] = n[1]->y - n[2]->y;  
  phiL[0][1] = n[2]->y - n[0]->y;
  phiL[0][2] = n[0]->y - n[1]->y;
  phiL[1][0] = n[2]->x - n[1]->x;
  phiL[1][1] = n[0]->x - n[2]->x;
  phiL[1][2] = n[1]->x - n[0]->x;  
  
  double k[3][3];
  for (i = 0; i < 3; i++)
  {
    for(j = 0; j < 3; j++)
    {
      k[i][j] = elMat->k * (phiL[0][i]*phiL[0][j] + phiL[1][i]*phiL[1][j])/den;
    }
  }
 
  double fb = elMat->fp * area / 3;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      matrix_sumelem(stifMatrix, n[i]->index, n[j]->index, k[i][j]);
    }
    vector_sumelem(b, n[i]->index, fb);
  }
  
  return 0;
}

int eleP1Elast(Matrix *stifMatrix, Vector *b, const Element *element)
{
  Material *elMat = element->material;

  const Node *n[3];
  n[0] = element->n[0];
  n[1] = element->n[1];
  n[2] = element->n[2];
  
  /* calc auxiliar parameters */
  double det = n[0]->x*n[1]->y + n[0]->y*n[2]->x + n[1]->x*n[2]->y - n[1]->y*n[2]->x - n[2]->y*n[0]->x - n[0]->y*n[1]->x;
  double area = det/2;
  
  double phiL[2][3];
  phiL[0][0] = n[1]->y - n[2]->y;  
  phiL[0][1] = n[2]->y - n[0]->y;
  phiL[0][2] = n[0]->y - n[1]->y;
  phiL[1][0] = n[2]->x - n[1]->x;
  phiL[1][1] = n[0]->x - n[2]->x;
  phiL[1][2] = n[1]->x - n[0]->x;
  
  /* 
  [ a0  0   a1  0   a2  0  ]
  [ 0   b0  0   b1  0   b2 ]
  [ b0  a0  b1  a1  b2  a2 ]
  */
  StdMatrix *bMatrix;
  stdmatrix_create(&bMatrix, 3, 6);
  stdmatrix_setelem(bMatrix, 0, 0, phiL[0][0]);
  stdmatrix_setelem(bMatrix, 0, 2, phiL[0][1]);
  stdmatrix_setelem(bMatrix, 0, 4, phiL[0][2]);
  stdmatrix_setelem(bMatrix, 1, 1, phiL[1][0]);
  stdmatrix_setelem(bMatrix, 1, 3, phiL[1][1]);
  stdmatrix_setelem(bMatrix, 1, 5, phiL[1][2]);
  stdmatrix_setelem(bMatrix, 2, 0, phiL[1][0]);
  stdmatrix_setelem(bMatrix, 2, 1, phiL[0][0]);
  stdmatrix_setelem(bMatrix, 2, 2, phiL[1][1]);
  stdmatrix_setelem(bMatrix, 2, 3, phiL[0][1]);
  stdmatrix_setelem(bMatrix, 2, 4, phiL[1][2]);
  stdmatrix_setelem(bMatrix, 2, 5, phiL[0][2]);
  
  StdMatrix *cMatrix;
  stdmatrix_create(&cMatrix, 3, 3);
  /* plane strain */
  /*
  [ l+2mu  l       0 ]
  [  l   l + 2mu   0 ]
  [ 0      0      mu ]
  */
  
  double c = elMat->l+2*elMat->mu;
  stdmatrix_setelem(cMatrix, 0, 0, c);
  stdmatrix_setelem(cMatrix, 0, 1, elMat->l);
  stdmatrix_setelem(cMatrix, 1, 0, elMat->l);
  stdmatrix_setelem(cMatrix, 1, 1, c);
  stdmatrix_setelem(cMatrix, 2, 2, elMat->mu);
  
  /* plane stress */
  /*
            [1    v      0  ]
  E/(1-v^2) [v    1      0  ]
            [0    0  (1-v)/2]
  */
  /*
  double constant = 2600/(1-0.3*0.3);
  stdmatrix_setelem(cMatrix, 0, 0, constant);
  stdmatrix_setelem(cMatrix, 0, 1, constant * 0.3);
  stdmatrix_setelem(cMatrix, 1, 0, constant * 0.3);
  stdmatrix_setelem(cMatrix, 1, 1, constant);
  stdmatrix_setelem(cMatrix, 2, 2, constant * (1-0.3/2));  */
  
  StdMatrix *temp;
  stdmatrix_create(&temp, 6, 3);
  StdMatrix *kMatrix;
  stdmatrix_create(&kMatrix, 6, 6);
  
  stdmatrix_multiplymTn(bMatrix, cMatrix, temp);
  stdmatrix_multiply(temp, bMatrix, kMatrix);
  
  /* cleanup matrixes */
  stdmatrix_destroy(temp);
  stdmatrix_destroy(bMatrix);
  stdmatrix_destroy(cMatrix);
  
  double fbVal = elMat->fg * area / 3;
  for (uint32_t i = 0; i < 6; i++)
  {
    uint32_t idx0 = assemblyGlobalNumber(n, i, 2);
    for (uint32_t j = 0; j < 6; j++)
    {
      uint32_t idx1 = assemblyGlobalNumber(n, j, 2);
      
      matrix_sumelem(stifMatrix, idx0, idx1, stdmatrix_getelem(kMatrix, i, j) /2/det);
    }
    
    double fb, fp;
    double pInt = (n[0]->u + n[1]->u + n[2]->u)/3 * area;
    if (i % 2 == 0)
    {
      fb = 0;
      fb = -4*elMat->mu*area/3;
      fp = phiL[0][i/2]/det * pInt;
    }
    else
    {
      fb = fbVal;
      fb = -6*elMat->mu*area/3;
      fp = phiL[1][i/2]/det * pInt;
    }
    vector_sumelem(b, idx0, fb + fp);
  }
  
  stdmatrix_destroy(kMatrix);
  
  return 0;
}

int computeStresses(const ElementArray *elements)
{
  /* sEf = \lambda div(u)I + 2\mu \epsilon(u) */
  for (uint32_t i = 0; i < elements->len; i++)
  { 
    Element *curElement = &elements->ptr[i];
    Material *elMat = curElement->material;
    
    const Node *n[3];
    n[0] = curElement->n[0];
    n[1] = curElement->n[1];
    n[2] = curElement->n[2];
    
    double det = n[0]->x*n[1]->y + n[0]->y*n[2]->x + n[1]->x*n[2]->y - n[1]->y*n[2]->x - n[2]->y*n[0]->x - n[0]->y*n[1]->x;
    
    double phiL[2][3];
    phiL[0][0] = (n[1]->y - n[2]->y)/det;
    phiL[0][1] = (n[2]->y - n[0]->y)/det;
    phiL[0][2] = (n[0]->y - n[1]->y)/det;
    phiL[1][0] = (n[2]->x - n[1]->x)/det;
    phiL[1][1] = (n[0]->x - n[2]->x)/det;
    phiL[1][2] = (n[1]->x - n[0]->x)/det;
    
    /* calc div */
    double divu = 0;
    for (uint32_t nIndex = 0; nIndex < 3; nIndex++)
    {
      for (uint32_t jDir = 0; jDir < 2; jDir++)
      {
        divu += n[nIndex]->ui[jDir] * phiL[jDir][nIndex];
      }
    }
    
    /* calc \epsilon */
    double ep[2][2] = {0};
    for (uint32_t nIndex = 0; nIndex < 3; nIndex++)
    {
      ep[0][0] += n[nIndex]->ui[0]*phiL[0][nIndex];
      ep[0][1] += n[nIndex]->ui[0]*phiL[1][nIndex] + n[nIndex]->ui[1]*phiL[0][nIndex];
      ep[1][1] += n[nIndex]->ui[1]*phiL[1][nIndex];
    }
    ep[0][1] /= 2;
    ep[1][0] = ep[0][1];
    
    /* fill stress tensor */
    double lambda = elMat->l; double mu = elMat->mu;
    curElement->sEf[0][0] = lambda*divu + 2*mu*ep[0][0];
    curElement->sEf[0][1] =               2*mu*ep[0][1];
    curElement->sEf[1][0] =               2*mu*ep[1][0];
    curElement->sEf[1][1] = lambda*divu + 2*mu*ep[1][1];
  }
  return 0;
}

int loadResultIntoNodePoisson(NodeArray *nodes, Vector *x)
{
  uint32_t i;
  for (i = 0; i < nodes->len; i++)
  {
    nodes->ptr[i].u = vector_getelem(x, i);
  }
  return 0;
}


int loadResultIntoNodeElast(NodeArray *nodes, Vector *x)
{
  uint32_t i;
  for (i = 0; i < nodes->len; i++)
  {
    nodes->ptr[i].ui[0] = vector_getelem(x, 2*i);
    nodes->ptr[i].ui[1] = vector_getelem(x, 2*i+1);
  }
  return 0;
}

int applyNeuBoundary1(Vector *b, const uint32_t nNeumann, const NeumannBC *neuBound, const NodeArray *nodes)
{
  Node n[2];
  for (uint32_t i = 0; i < nNeumann; i++)
  {
    n[0] = nodes->ptr[neuBound[i].n0];
    n[1] = nodes->ptr[neuBound[i].n1];
    
    double sLength = getSideLength(n[0], n[1]);
    double val = neuBound->val[0] * sLength /2;
    
    vector_sumelem(b, n[0].index, val);
    vector_sumelem(b, n[1].index, val);
  }
  
  return 0;
}

int applyNeuBoundary2(Vector *b, const uint32_t nNeumann, const NeumannBC *neuBound, const NodeArray *nodes)
{
  for (uint32_t i = 0; i < nNeumann; i++)
  {
    Node n[2];
    n[0] = nodes->ptr[neuBound[i].n0];
    n[1] = nodes->ptr[neuBound[i].n1];
    double sLength = getSideLength(n[0], n[1]);
    
    vector_sumelem(b, n[0].index*2, neuBound[i].val[0]*sLength/2);
    vector_sumelem(b, n[0].index*2 + 1, neuBound[i].val[1]*sLength/2);
    vector_sumelem(b, n[1].index*2, neuBound[i].val[0]*sLength/2);
    vector_sumelem(b, n[1].index*2 + 1, neuBound[i].val[1]*sLength/2);
  }
  
  return 0;
}

int applyDirBoundary(Matrix *stifMatrix, Vector *b, const uint32_t nDirichlet, const DirichletBC *dirBound, const uint32_t ndof)
{
  uint32_t nodeIdx, dirIdx;
  double val;
  uint32_t i;
  for (i = 0; i < nDirichlet; i++)
  {
    nodeIdx = dirBound[i].nodeIndex;
    dirIdx = dirBound[i].dir;
    val = dirBound[i].val;
    matrix_applyDirBC(stifMatrix, b,  gIndex(ndof, nodeIdx, dirIdx), val);
  }
  
  return 0;
}


int conjugateGradientMethod(const MatrixCRS *a, Vector *x, const Vector *b)
{
  uint32_t n = vector_size(x);

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
  uint32_t i;
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

inline uint32_t assemblyGlobalNumber(const Node *n[3], const uint32_t i, const uint32_t ndof)
{
  return n[i/ndof]->index*ndof + i%ndof;
}

inline uint32_t gIndex(uint32_t nj, uint32_t i, uint32_t j)
{
  return i*nj + j;
}

inline double getSideLength(Node n0, Node n1)
{
  return sqrt(pow(n1.x - n0.x,2) + pow(n1.y - n0.y,2));
}