#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include "matrix.h"
#include "helper.h"
#include "stdMatrix.h"

#define EPS 1e-8
#define TOL1 1e-6
#define TOL2 1e-6
#define MAXITER 100
#define LOADSTEPS 1

struct node {
  uint32_t index;
  double x;
  double y;
  double u;
  double ui[2];  
  uint32_t bm;

  double DUI[2];
  double dui[2];
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
  double yieldStrengh;
};
typedef struct material Material;

struct element {
  uint32_t index;
  Material *material;
  Node *n[3];
  double vd;
  StdMatrix *stress;
  StdMatrix *eps;
  StdMatrix *epsP;
  StdMatrix *Cep;
};
typedef struct element Element;

struct elementArray{
  uint32_t len;
  Element *ptr;
};
typedef struct elementArray ElementArray;

struct dirichletBCFile {
  uint32_t bm;
  uint32_t dir;
  double val;
};
typedef struct dirichletBCFile DirichletBCFile;

struct dirichletBC {
  uint32_t nodeIndex;
  uint32_t dir;
  double val;

  /* plast */
  double fk;
  double Dfk;
};
typedef struct dirichletBC DirichletBC;

struct dirichletBCArray {
  uint32_t len;
  DirichletBC *ptr;
};
typedef struct dirichletBCArray DirichletBCArray;

struct neumannBCFile {
  uint32_t bm;
  double val[3];
};
typedef struct neumannBCFile NeumannBCFile;

struct neumannBC {
  uint32_t n0;
  uint32_t n1;
  double val[3];

  /* plast */
  double fk[3];
  double Dfk[3];
};
typedef struct neumannBC NeumannBC;

struct neumannBCArray {
  uint32_t len;
  NeumannBC *ptr;
};
typedef struct neumannBCArray NeumannBCArray;

int readMaterials(uint32_t *nMaterial, Material **materials, char *fileName);
int readNodes(NodeArray *nodes, char *fileName);
int readElements(ElementArray *elements, const NodeArray *nodes, const Material *materials, char *fileName);
int readBC(DirichletBCArray *dirBCArray, NeumannBCArray *neuBCArray, const NodeArray *nodes, char *fullPath, char *ext);
int readBCE(DirichletBCArray *dirBCArray, NeumannBCArray *neuBCArray, const NodeArray *nodes, char *fileName);
int readBCP(DirichletBCArray *dirBCArray, NeumannBCArray *neuBCArray, const NodeArray *nodes, char *fileName);

int solveFEM(ElementArray *elements, NodeArray *nodes, const DirichletBCArray *dirBCArray, const NeumannBCArray *neuBCArray, const uint32_t mode);

int computeGlobalMatrices(Matrix **a, Vector **b, const ElementArray *elements, const NodeArray *nodes, const DirichletBCArray *dirBCArray, const NeumannBCArray *neuBCArray, const uint32_t mode);

int eleP1Poisson(Matrix *stifMatrix, Vector *b, const Element *element);
int eleP1Elast(Matrix *stifMatrix, Vector *b, const Element *element);
int eleP1ElastPlast(Matrix *stifMatrix, Vector *b, const Element *element);

int closestPointProjection(Element *element);

int computeStress(const StdMatrix *extCMatrix, const StdMatrix *totalStrain, const StdMatrix *plasticStrain, StdMatrix *stress);
int computeStresses(const ElementArray *elements);

double calcNMatrix(StdMatrix *N, const Element *element);
int calcBMatrix(StdMatrix *B, const StdMatrix *phiL);
int calcCMatrix(StdMatrix *C, const Element *element);
int calcPMatrix(StdMatrix *P);
int calcExtendedCMatrix(StdMatrix *extC, const Element *element);
int calcXIMatrix(const StdMatrix *stress, const StdMatrix *extCInv, double delGamma, StdMatrix *fll, StdMatrix *xiMatrix);

int updateIncrement(DirichletBCArray *dirBCArray, NeumannBCArray *neuBCArray, uint32_t curStep, uint32_t numLoadSteps);
int initializePlasticVariables(ElementArray *elements, NodeArray *nodes);

int applyNeuBoundary1(Vector *b, const NeumannBCArray *neuBCArray, const NodeArray *nodes);
int applyNeuBoundary2(Vector *b, const NeumannBCArray *neuBCArray, const NodeArray *nodes);
int applyDirBoundary(Matrix *stifMatrix, Vector *b, const DirichletBCArray *dirBCArray, const uint32_t ndof);

int conjugateGradientMethod(const MatrixCRS *a, Vector *x, const Vector *b);

int loadResultIntoNodePoisson(NodeArray *nodes, Vector *x);
int loadResultIntoNodeElast(NodeArray *nodes, Vector *x);
int incrementResultIntoNodeElast(NodeArray *nodes, Vector *x);

double plasticYieldF(const StdMatrix *stressVector, double flowStress);
int calcPlasticYieldFl(const StdMatrix *stressVector, StdMatrix *r);
double calcPlasticYieldFll(const StdMatrix *stressVector, StdMatrix *fll);

int calcStressDeviator(const StdMatrix *stressVector, StdMatrix *r);

int printVTK(NodeArray *nodes, ElementArray *elements, char *fileName);

inline uint32_t gIndex(uint32_t nj, uint32_t i, uint32_t j);
inline double getSideLength(Node n0, Node n1);
inline uint32_t assemblyGlobalNumber(const Node *n[3], const uint32_t i, const uint32_t ndof);

int change3to6Vector(const StdMatrix *m, StdMatrix *r);
int change3to6Vector(const StdMatrix *m, StdMatrix *r);
int change6to3Matrix(const StdMatrix *m, StdMatrix *r);

int main( int argc, char *argv[] ) 
{
  NodeArray nodes;
  ElementArray elements;
  DirichletBCArray dirBCArray;
  NeumannBCArray neuBCArray;
  
  uint32_t nMaterial;
  Material *materials;
  
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
  
  /* (a) Input of data */
  /* (b) Representation of the triangulation Th. */
  readMaterials(&nMaterial, &materials, fileName);
  readNodes(&nodes, fileName);
  readElements(&elements, &nodes, materials, fileName);
  readBCE(&dirBCArray, &neuBCArray, &nodes, fileName);

  for (uint32_t k = 0; k <= LOADSTEPS; k++)
  {
    updateIncrement(&dirBCArray, &neuBCArray, 1, LOADSTEPS);
    /* main Newton method */

    /* initialize plastic load step */
    initializePlasticVariables(&elements, &nodes);
    for (uint32_t j = 0; j < MAXITER; j++)
    {
      /*1. Compute elastic predictor */
      solveFEM(&elements, &nodes, &dirBCArray, &neuBCArray, 2);

      int isPlastic = 0;
      for (uint32_t i = 0; i < elements.len; i++)
      {
        /* for each gauss point */
        if (closestPointProjection(&elements.ptr[i]) == 0)
        {
          isPlastic = 1;
        }
      }
      if (isPlastic == 0) break;
    }


  }
  
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
  fprintf(fp, "\n");
  fprintf(fp, "CELL_DATA %"PRIu32"\n", elements->len);
  
  fprintf(fp, "TENSORS efStress float\n");
  for(i = 0; i < elements->len; i++)
  {
    fprintf(fp, "%.6g %.6g 0\n", stdmatrix_getelem(elements->ptr[i].stress, 0, 0), stdmatrix_getelem(elements->ptr[i].stress, 3, 0));
    fprintf(fp, "%.6g %.6g 0\n", stdmatrix_getelem(elements->ptr[i].stress, 3, 0), stdmatrix_getelem(elements->ptr[i].stress, 1, 0));
    fprintf(fp, "0    0    %.6g\n", stdmatrix_getelem(elements->ptr[i].stress, 2, 0));
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
  nodes->ptr = sMalloc(n * sizeof(Node));
  uint32_t i;
  double x, y;
  uint32_t index;
  uint32_t bm;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNu32" %lf %lf %" SCNu32"\n", &index, &x, &y, &bm);
    (nodes->ptr + i)->index = index;
    (nodes->ptr + i)->x = x;
    (nodes->ptr + i)->y= y;
    (nodes->ptr + i)->bm = bm;
    (nodes->ptr + i)->u = 0;
    (nodes->ptr + i)->ui[0] = 0;
    (nodes->ptr + i)->ui[1] = 0;
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
  elements->ptr = sMalloc(n * sizeof(Element));
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
    stdmatrix_create(&elements->ptr[i].eps, 6, 1);
    stdmatrix_create(&elements->ptr[i].epsP, 6, 1);
    stdmatrix_create(&elements->ptr[i].stress, 6, 1);
    stdmatrix_create(&elements->ptr[i].Cep, 6, 6);
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
  Material *newMaterials = sMalloc(n * sizeof(Material));
  uint32_t i;
  double k, fPoisson;
  double e, nu, density;
  uint32_t index;
  double yieldStrengh;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNu32" %lf %lf %lf %lf %lf %lf*\n", &index, &k, &fPoisson, &e, &nu, &density, &yieldStrengh);
    (*(newMaterials + i)).index = index;
    (*(newMaterials + i)).k = k;
    (*(newMaterials + i)).fp = fPoisson;
    (*(newMaterials + i)).l = -e*nu/(nu*(2*nu+1) - 1);
    (*(newMaterials + i)).mu = e/(2*nu + 2);
    (*(newMaterials + i)).fg = density * 9.80665;
    (*(newMaterials + i)).yieldStrengh = yieldStrengh;
  }
  
  *nMaterial = n;
  *materials = newMaterials;
  
  fclose(fp);
  
  return 0;
}

int readBCP(DirichletBCArray *dirBCArray, NeumannBCArray *neuBCArray, const NodeArray *nodes, char *fileName)
{
  return readBC(dirBCArray, neuBCArray, nodes, fileName, ".bcp");
}

int readBCE(DirichletBCArray *dirBCArray, NeumannBCArray *neuBCArray, const NodeArray *nodes, char *fileName)
{
  return readBC(dirBCArray, neuBCArray, nodes, fileName, ".bce");
}

int readBC(DirichletBCArray *dirBCArray, NeumannBCArray *neuBCArray, const NodeArray *nodes, char *fileName, char *ext)
{
  uint32_t n;

  char fullPath[80];
  snprintf(fullPath, sizeof(fullPath), "%s%s", fileName, ext);
  FILE *fp;
  fp = fopen(fullPath, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "readBC: Cannot open file %s.\n", fullPath);
    return 1;
  }
  
  fscanf(fp, "%" SCNu32 "*\n", &n);
  if (n < 1)
  {
    fprintf(stderr, "readBC: No dirichlet boundary conditions to be read\n");
    fclose(fp);
    return 1;
  }
  
  /* There are bcs, allocate memory and reads them */
  DirichletBCFile *dirBCFile = sMalloc(n * sizeof(DirichletBCFile));
  
  uint32_t bm;
  uint32_t dir;
  double value[3];
  for (uint32_t i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNu32" %"SCNu32" %lf*\n", &bm, &dir, &value[0]);
    (*(dirBCFile + i)).bm = bm;
    (*(dirBCFile + i)).dir = dir;
    (*(dirBCFile + i)).val = value[0];
  }

  /* loops through all the nodes, searching for the boundary markers */
  uint32_t k = 0;
  dirBCArray->ptr = sMalloc(nodes->len * sizeof(DirichletBC));
  for (uint32_t i = 0; i < nodes->len; i++)
  {
    if (nodes->ptr[i].bm > 0)
    {
      /* has a boundary marker: check which one, and create a boundary */
      for (uint32_t bci = 0; bci < n; bci++)
      {
        if (dirBCFile[bci].bm == nodes->ptr[i].bm)
        {
          dirBCArray->ptr[k].nodeIndex = nodes->ptr[i].index;
          dirBCArray->ptr[k].dir = dirBCFile[bci].dir;
          dirBCArray->ptr[k].val = dirBCFile[bci].val;
          k++;
        }
      }
    }
  }

  if (k == 0)
  {
    fprintf(stderr, "readBC: No dirichlet boundary conditions found on nodes.\n");
    fclose(fp);
    return 1;
  }
  /* shrinks the array to the correct size */
  dirBCArray->ptr = sRealloc(dirBCArray->ptr, k * sizeof(DirichletBC));

  /* excludes the dirichletBCFiles, it is not needed anymore. */
  free(dirBCFile);

  /* stores the number of dirichlet nodes */
  dirBCArray->len = k;


  /* reads the neumann part of the file */
  fscanf(fp, " %" SCNu32 "*\n", &n);
  if (n == 0)
  {
    /* doesn't have any neumann condition*/
    neuBCArray->len = 0;
    fclose(fp);
    return 0;
  }

  NeumannBCFile *neuBCFile = sMalloc(n * sizeof(NeumannBCFile));
  for (uint32_t i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNu32" %lf %lf %lf\n", &bm, &value[0], &value[1], &value[2]);
    neuBCFile[i].bm = bm;
    neuBCFile[i].val[0] = value[0];
    neuBCFile[i].val[1] = value[1];
    neuBCFile[i].val[2] = value[2];
  }

  /* loops through all the edges, searching for the boundary markers */
  FILE *fpEdge;
  snprintf(fullPath, sizeof(fullPath), "%s.edge", fileName);
  fpEdge = fopen(fullPath, "r");
  if (fpEdge == NULL)
  {
    fprintf(stderr, "readBC: Cannot open file %s.\n", fullPath);
    exit(-1);
  }

  uint32_t nEdge;
  fscanf(fpEdge, "%" SCNu32 " %*[^\n]\n", &nEdge);
  if (nEdge < 1)
  {
    fprintf(stderr, "readBC: No edges to be read\n");
    exit(-1);
  }

  k = 0;
  neuBCArray->ptr = sMalloc(nEdge * sizeof(NeumannBC));
  for (uint32_t i = 0; i < nEdge; i++)
  {
    uint32_t n0, n1;
    fscanf(fpEdge, "%*" SCNu32 " %" SCNu32 " %" SCNu32 " %" SCNu32 "\n", &n0, &n1, &bm);
    if (bm > 0)
    {
      /* has a boundary marker: checks which one, and create a boundary */
      for (uint32_t bci = 0; bci < n; bci++)
      {
        if (neuBCFile[bci].bm == bm)
        {
          neuBCArray->ptr[k].n0 = n0;
          neuBCArray->ptr[k].n1 = n1;
          neuBCArray->ptr[k].val[0] = neuBCFile[bci].val[0];
          neuBCArray->ptr[k].val[1] = neuBCFile[bci].val[1];
          neuBCArray->ptr[k].val[2] = neuBCFile[bci].val[2];
          k++;
          break;
        }
      }
    }
  }
  fclose(fpEdge);

  if (k > 0)
  {
    /* found Neumann edges */
    /* shrinks the array to the correct size */
    neuBCArray->ptr = sRealloc(neuBCArray->ptr, k * sizeof(NeumannBC));
  }
  else
  {
    /* doesn't have neumann edges, frees the memory. */
    free(neuBCArray->ptr);
  }

  /* excludes the neumannBCFiles, it is not needed anymore. */
  free(neuBCFile);

  /* stores the number of neumann edges */
  neuBCArray->len = k;

  fclose(fp);
  
  return 0;
}

/*  mode 0 : poisson 
    mode 1 : elasticity
    mode 2 : incremental elastoplasticity

*/
int solveFEM(ElementArray *elements, NodeArray *nodes, const DirichletBCArray *dirBCArray, const NeumannBCArray *neuBCArray, const uint32_t mode)
{
  uint32_t ndof;
  switch (mode)
  {
  case 0:
    ndof = 1;
    break;
  case 1:
  case 2:
    ndof = 2;
    break;
  default:
    fprintf(stderr, "solveFEM: inexistant mode chosen.\n");
    exit(-1);
    break;
  }

  /* (c) Computation of the element stiffness matrices aK and element loads bK. */
  /* (d) Assembly of the global stiffness matrix A and load vector b. */
  Matrix *a;
  Vector *b;
  computeGlobalMatrices(&a, &b, elements, nodes, dirBCArray, neuBCArray, mode);

  /* (e) Solution of the system of equations Ax = b. */

  if (!matrix_issimetric(a))
  {
    fprintf(stderr, "main: matriz não simétrica.\n");
    exit(-1);
  }
  MatrixCRS *aCrs;
  matrixcrs_create(a, &aCrs);
  Vector *x;
  vector_createzero(&x, nodes->len * ndof);
  conjugateGradientMethod(aCrs, x, b);

  switch (mode)
  {
  case 0:
    loadResultIntoNodePoisson(nodes, x);
    break;
  case 1:
    loadResultIntoNodeElast(nodes, x);
    break;
  case 2:
    incrementResultIntoNodeElast(nodes, x);
    break;
  default:
    break;
  }

  /* (g) Cleanup. */
  matrix_destroy(a);
  vector_destroy(b);
  matrixcrs_destroy(aCrs);
  vector_destroy(x);

  return 0;
}

int computeGlobalMatrices(Matrix **a, Vector **b, const ElementArray *elements, const NodeArray *nodes, const DirichletBCArray *dirBCArray, const NeumannBCArray *neuBCArray, const uint32_t mode)
{
  uint32_t ndof;
  switch (mode)
  {
  case 0:
    ndof = 1;
    break;
  case 1:
  case 2:
    ndof = 2;
    break;
  default:
    fprintf(stderr, "computeGlobalMatrices: inexistant mode chosen.\n");
    exit(-1);
    break;
  }

  Matrix *newStif;
  matrix_create(&newStif, nodes->len*ndof, nodes->len*ndof);
  Vector *newB;
  vector_createzero(&newB, nodes->len*ndof);
  
  uint32_t i;
  for (i = 0; i < elements->len; i++)
  {
    switch (mode)
    {
    case 0:
      eleP1Poisson(newStif, newB, &elements->ptr[i]);
      break;
    case 1:
      eleP1Elast(newStif, newB, &elements->ptr[i]);
      break;
    case 2:
      eleP1ElastPlast(newStif, newB, &elements->ptr[i]);
      break;
    default:
      fprintf(stderr, "computeGlobalMatrices: inexistant mode chosen.\n");
      exit(-1);
      break;
    }
  }
  
  /* apply neumann boundary */
  if (ndof == 1)
  {
    applyNeuBoundary1(newB, neuBCArray, nodes);
  }
  else
  {
    applyNeuBoundary2(newB, neuBCArray, nodes);
  }
  
  /* apply dir boundary */
  applyDirBoundary(newStif, newB, dirBCArray, ndof);
  
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
  
  /* calc shape function */
  StdMatrix *N;
  stdmatrix_create(&N, 2, 3);
  double det = calcNMatrix(N, element);
  
  StdMatrix *bMatrix;
  stdmatrix_create(&bMatrix, 3, 6);
  calcBMatrix(bMatrix, N);
  
  StdMatrix *cMatrix;
  stdmatrix_create(&cMatrix, 3, 3);
  calcCMatrix(cMatrix, element);
  
  StdMatrix *temp;
  stdmatrix_create(&temp, 6, 3);
  StdMatrix *kMatrix;
  stdmatrix_create(&kMatrix, 6, 6);
  
  stdmatrix_multiplymTn(bMatrix, cMatrix, temp);
  stdmatrix_multiply(temp, bMatrix, kMatrix);
  
  for (uint32_t i = 0; i < 6; i++)
  {
    uint32_t idx0 = assemblyGlobalNumber(n, i, 2);
    for (uint32_t j = 0; j < 6; j++)
    {
      uint32_t idx1 = assemblyGlobalNumber(n, j, 2);
      
      matrix_sumelem(stifMatrix, idx0, idx1, stdmatrix_getelem(kMatrix, i, j) * det / 2);
    }
    
    double fb = 0;
    double fp = 0;
    /* double pInt = (n[0]->u + n[1]->u + n[2]->u)/3 * area; */
    if (i % 2 == 0)
    {
      fb = 0;
      /* fp = phiL[0][i/2]/det * pInt; */
    }
    else
    {
      fb = -elMat->fg * det / 6;
      /* fp = phiL[1][i/2]/det * pInt; */
    }
    vector_sumelem(b, idx0, fb + fp);
  }
  
  /* cleanup matrixes */
  stdmatrix_destroy(temp);
  stdmatrix_destroy(bMatrix);
  stdmatrix_destroy(cMatrix);
  stdmatrix_destroy(kMatrix);
  stdmatrix_destroy(N);
  
  return 0;
}

int eleP1ElastPlast(Matrix *stifMatrix, Vector *b, const Element *element)
{
  Material *elMat = element->material;

  const Node *n[3];
  n[0] = element->n[0];
  n[1] = element->n[1];
  n[2] = element->n[2];

  /* calc shape function */
  StdMatrix *N;
  stdmatrix_create(&N, 2, 3);
  double det = calcNMatrix(N, element);

  StdMatrix *bMatrix;
  stdmatrix_create(&bMatrix, 3, 6);
  calcBMatrix(bMatrix, N);

  StdMatrix *temp;
  stdmatrix_create(&temp, 6, 3);
  StdMatrix *kMatrix;
  stdmatrix_create(&kMatrix, 6, 6);

  StdMatrix *Cep33;
  stdmatrix_create(&Cep33, 3, 3);

  StdMatrix *fVector;
  stdmatrix_create(&fVector, 6, 1);

  StdMatrix *ifVector;
  stdmatrix_create(&ifVector, 6, 1);

  StdMatrix *temp61;
  stdmatrix_create(&temp61, 6, 1);

  StdMatrix *stress31;
  stdmatrix_create(&stress31, 3, 1);

  /* Ctm is 6x6. Must reduce to 3x3 */
  change6to3Matrix(element->Cep, Cep33);

  /* K Matrix*/
  stdmatrix_multiplymTn(bMatrix, Cep33, temp);
  stdmatrix_multiply(temp, bMatrix, kMatrix);

  stdmatrix_multiplyScalar(kMatrix, det / 2); /* integrate */

  /* Force vector */
  stdmatrix_setelem(fVector, 0, 1, -elMat->fg);
  stdmatrix_setelem(fVector, 0, 3, -elMat->fg);
  stdmatrix_setelem(fVector, 0, 5, -elMat->fg);

  /* det/2 = area,   x/3 \basis function */
  stdmatrix_multiplyScalar(fVector, det / 2 / 3); /* integrate */

  /* additional internal force */
  change3to6Vector(element->stress, stress31);
  stdmatrix_multiplymTn(bMatrix, stress31, ifVector);

  stdmatrix_multiplyScalar(ifVector, det / 2 / 3);  /* integrate */

  /* calc resisual */
  stdmatrix_subtract(fVector, ifVector, temp61);

  /* assembly */
  for (uint32_t i = 0; i < 6; i++)
  {
    uint32_t idx0 = assemblyGlobalNumber(n, i, 2);
    for (uint32_t j = 0; j < 6; j++)
    {
      uint32_t idx1 = assemblyGlobalNumber(n, j, 2);

      matrix_sumelem(stifMatrix, idx0, idx1, stdmatrix_getelem(kMatrix, i, j));
    }

    vector_sumelem(b, idx0, stdmatrix_getelem(temp61, i, 0));
  }

  /* cleanup matrixes */
  stdmatrix_destroy(temp);
  stdmatrix_destroy(bMatrix);
  stdmatrix_destroy(kMatrix);
  stdmatrix_destroy(N);
  stdmatrix_destroy(fVector);
  stdmatrix_destroy(ifVector);
  stdmatrix_destroy(stress31);
  stdmatrix_destroy(temp61);

  return 0;
}

double calcNMatrix(StdMatrix *N, const Element *element)
{
  if (N == NULL || stdmatrix_ni(N) != 2 || stdmatrix_nj(N) != 3)
  {
    fprintf(stderr, "calcNMatrix: mismatched matrix size.\n");
    exit(-1);
  }

  double det = element->n[0]->x*element->n[1]->y + 
    element->n[0]->y*element->n[2]->x + 
    element->n[1]->x*element->n[2]->y - 
    element->n[1]->y*element->n[2]->x - 
    element->n[2]->y*element->n[0]->x - 
    element->n[0]->y*element->n[1]->x;

  stdmatrix_setelem(N, 0, 0, (element->n[1]->y - element->n[2]->y) / det);
  stdmatrix_setelem(N, 0, 1, (element->n[2]->y - element->n[0]->y) / det);
  stdmatrix_setelem(N, 0, 2, (element->n[0]->y - element->n[1]->y) / det);
  stdmatrix_setelem(N, 1, 0, (element->n[2]->x - element->n[1]->x) / det);
  stdmatrix_setelem(N, 1, 1, (element->n[0]->x - element->n[2]->x) / det);
  stdmatrix_setelem(N, 1, 2, (element->n[1]->x - element->n[0]->x) / det);
  return det;
}

int calcBMatrix(StdMatrix *bMatrix, const StdMatrix *phiL)
{
  if (bMatrix == NULL || stdmatrix_ni(bMatrix) != 3 || stdmatrix_nj(bMatrix) != 6)
  {
    fprintf(stderr, "calcBMatrix: mismatched matrix size.\n");
    exit(-1);
  }

  uint32_t defRow;
  defRow = stdmatrix_ni(bMatrix) == 3 ? 2 : 3;


  /*
  ni = 3 
  [ a0  0   a1  0   a2  0  ]
  [ 0   b0  0   b1  0   b2 ]
  [ b0  a0  b1  a1  b2  a2 ]

  ni = 6
  [ a0  0   a1  0   a2  0  ]
  [ 0   b0  0   b1  0   b2 ]
  [ 0   0   0   0   0   0  ]
  [ b0  a0  b1  a1  b2  a2 ]
  [ 0   0   0   0   0   0  ]
  [ 0   0   0   0   0   0  ]
  */
  stdmatrix_setelem(bMatrix, 0, 0, stdmatrix_getelem(phiL, 0, 0));
  stdmatrix_setelem(bMatrix, 0, 2, stdmatrix_getelem(phiL, 0, 1));
  stdmatrix_setelem(bMatrix, 0, 4, stdmatrix_getelem(phiL, 0, 2));
  stdmatrix_setelem(bMatrix, 1, 1, stdmatrix_getelem(phiL, 1, 0));
  stdmatrix_setelem(bMatrix, 1, 3, stdmatrix_getelem(phiL, 1, 1));
  stdmatrix_setelem(bMatrix, 1, 5, stdmatrix_getelem(phiL, 1, 2));
  stdmatrix_setelem(bMatrix, defRow, 0, stdmatrix_getelem(phiL, 1, 0));
  stdmatrix_setelem(bMatrix, defRow, 1, stdmatrix_getelem(phiL, 0, 0));
  stdmatrix_setelem(bMatrix, defRow, 2, stdmatrix_getelem(phiL, 1, 1));
  stdmatrix_setelem(bMatrix, defRow, 3, stdmatrix_getelem(phiL, 0, 1));
  stdmatrix_setelem(bMatrix, defRow, 4, stdmatrix_getelem(phiL, 1, 2));
  stdmatrix_setelem(bMatrix, defRow, 5, stdmatrix_getelem(phiL, 0, 2));
  return 0;
}

int calcCMatrix(StdMatrix *C, const Element *element)
{
  if (C == NULL || stdmatrix_ni(C) != 3 || stdmatrix_nj(C) != 3)
  {
    fprintf(stderr, "calcCMatrix: mismatched matrix size.\n");
    exit(-1);
  }

  /* plane strain */
  /*
  [ l+2mu  l       0 ]
  [  l   l + 2mu   0 ]
  [ 0      0      mu ]
  */
  double c = element->material->l + 2 * element->material->mu;
  stdmatrix_setelem(C, 0, 0, c);
  stdmatrix_setelem(C, 0, 1, element->material->l);
  stdmatrix_setelem(C, 1, 0, element->material->l);
  stdmatrix_setelem(C, 1, 1, c);
  stdmatrix_setelem(C, 2, 2, element->material->mu);
  return 0;
}

/* Psigma = dev[sigma]*/
int calcPMatrix(StdMatrix *P)
{
  if (P == NULL || stdmatrix_ni(P) != 6 || stdmatrix_nj(P) != 6)
  {
    fprintf(stderr, "calcPMatrix: mismatched matrix size.\n");
    exit(-1);
  }
  /*
      [ 2/3  -1/3  -1/3   0   0   0]
      [-1/3   2/3  -1/3   0   0   0]
      [-1/3  -1/3   2/3   0   0   0]
      [ 0     0     0     1   0   0]
      [ 0     0     0     0   1   0]
      [ 0     0     0     0   0   1]
  */
  stdmatrix_setelem(P, 0, 0, 2 / 3);
  stdmatrix_setelem(P, 1, 1, 2 / 3);
  stdmatrix_setelem(P, 2, 2, 2 / 3);
  stdmatrix_setelem(P, 3, 3, 1);
  stdmatrix_setelem(P, 4, 4, 1);
  stdmatrix_setelem(P, 5, 5, 1);

  stdmatrix_setelem(P, 0, 1, -1 / 3);
  stdmatrix_setelem(P, 0, 2, -1 / 3);
  stdmatrix_setelem(P, 1, 2, -1 / 3);
  stdmatrix_setelem(P, 1, 0, -1 / 3);
  stdmatrix_setelem(P, 2, 0, -1 / 3);
  stdmatrix_setelem(P, 2, 1, -1 / 3);
  return 0;
}

int calcXIMatrix(const StdMatrix *stress, const StdMatrix *extCInv, double delGamma, StdMatrix *fll, StdMatrix *xiMatrix)
{
  StdMatrix *temp66;
  stdmatrix_create(&temp66, 6, 6);

  calcPlasticYieldFll(stress, fll);
  stdmatrix_addmultiplyscalar(extCInv, fll, delGamma, temp66);
  stdmatrix_invert(temp66, xiMatrix);

  stdmatrix_destroy(temp66);
  return 0;
}

int calcExtendedCMatrix(StdMatrix *extC, const Element *element)
{
  if (extC == NULL || stdmatrix_ni(extC) != 6 || stdmatrix_nj(extC) != 6)
  {
    fprintf(stderr, "calcExtendedCMatrix: mismatched matrix size.\n");
    exit(-1);
  }

  /*
  [ l+2mu    l       l      0   0   0  ]
  [  l     l + 2mu   l      0   0   0  ]
  [  l       l     l + 2mu  0   0   0  ]
  [  0       0       0      mu  0   0  ]
  [  0       0       0      0   mu  0  ]
  [  0       0       0      0   0   mu ]
  */
  double c = element->material->l + 2 * element->material->mu;
  stdmatrix_setelem(extC, 0, 0, c);
  stdmatrix_setelem(extC, 1, 1, c);
  stdmatrix_setelem(extC, 2, 2, c);
  stdmatrix_setelem(extC, 3, 3, element->material->mu);
  stdmatrix_setelem(extC, 4, 4, element->material->mu);
  stdmatrix_setelem(extC, 5, 5, element->material->mu);
  stdmatrix_setelem(extC, 0, 1, element->material->l);
  stdmatrix_setelem(extC, 1, 0, element->material->l);
  stdmatrix_setelem(extC, 0, 2, element->material->l);
  stdmatrix_setelem(extC, 2, 0, element->material->l);
  stdmatrix_setelem(extC, 1, 2, element->material->l);
  stdmatrix_setelem(extC, 2, 1, element->material->l);
  return 0;
}

int computeStress(const StdMatrix *extCMatrix, const StdMatrix *totalStrain, const StdMatrix *plasticStrain, StdMatrix *stress)
{
  StdMatrix *elasticStrain;
  stdmatrix_create(&elasticStrain, 6, 1);

  stdmatrix_subtract(totalStrain, plasticStrain, elasticStrain);

  stdmatrix_multiply(extCMatrix, elasticStrain, stress);

  stdmatrix_destroy(elasticStrain);
  return 0;
}

int computeStresses(const ElementArray *elements)
{
  /* \sigma = C\epsilon(u) */
  StdMatrix *N;
  stdmatrix_create(&N, 2, 3);

  StdMatrix *uColumnVector;
  stdmatrix_create(&uColumnVector, 6, 1);

  StdMatrix *bMatrix;
  stdmatrix_create(&bMatrix, 3, 6);

  StdMatrix *extCMatrix;
  stdmatrix_create(&extCMatrix, 6, 6);

  StdMatrix *epsilon;
  stdmatrix_create(&epsilon, 3, 1);

  for (uint32_t i = 0; i < elements->len; i++)
  {
    Element *curElement = &elements->ptr[i];

    const Node *n[3];
    n[0] = curElement->n[0];
    n[1] = curElement->n[1];
    n[2] = curElement->n[2];

    calcNMatrix(N, curElement);

    /* calc \epsilon */
    /* \epsilon(u) = [B]{u} */
    calcBMatrix(bMatrix, N);

    for (uint32_t nIndex = 0; nIndex < 3; nIndex++)
    {
      stdmatrix_setelem(uColumnVector, nIndex * 2, 0, n[nIndex]->ui[0]);
      stdmatrix_setelem(uColumnVector, nIndex * 2 + 1, 0, n[nIndex]->ui[1]);
    }

    stdmatrix_multiply(bMatrix, uColumnVector, epsilon);

    change3to6Vector(epsilon, curElement->eps);

    /* calc stress */
    calcExtendedCMatrix(extCMatrix, curElement);

    computeStress(extCMatrix, curElement->eps, curElement->epsP, curElement->stress);
  }

  /* clean matrices */
  stdmatrix_destroy(N);
  stdmatrix_destroy(uColumnVector);
  stdmatrix_destroy(bMatrix);
  stdmatrix_destroy(extCMatrix);
  stdmatrix_destroy(epsilon);

  return 0;
}

int closestPointProjection(Element *element)
{
  /* closest point projection, computational inelasticity */
  /* BOX 4.1, pag 173 */

  if (plasticYieldF(element->stress, element->material->yieldStrengh) <= 0)
  {
    /* elastic gauss point */
    return 1;
  }
  double fKn1;

  StdMatrix *extC;
  stdmatrix_create(&extC, 6, 6);
  calcExtendedCMatrix(extC, element);

  StdMatrix *extCInv;
  stdmatrix_create(&extCInv, 6, 6);
  stdmatrix_invert(extC, extCInv);

  StdMatrix *stressKn1;
  stdmatrix_create(&stressKn1, 6, 1);

  StdMatrix *strainPKn1;
  stdmatrix_create(&strainPKn1, 6, 1);

  StdMatrix *rKn1;
  stdmatrix_create(&rKn1, 6, 1);

  StdMatrix *fl;
  stdmatrix_create(&fl, 6, 1);

  StdMatrix *fll;
  stdmatrix_create(&fll, 6, 6);

  StdMatrix *xiMatrix;
  stdmatrix_create(&xiMatrix, 6, 6);

  StdMatrix *temp11;
  stdmatrix_create(&temp11, 1, 1);

  StdMatrix *temp61;
  stdmatrix_create(&temp61, 6, 1);

  StdMatrix *temp16;
  stdmatrix_create(&temp16, 1, 6);

  StdMatrix *temp66;
  stdmatrix_create(&temp66, 6, 6);

  StdMatrix *delEpsPk;
  stdmatrix_create(&delEpsPk, 6, 1);

  double delGamma = 0;
  stdmatrix_copy(element->epsP, strainPKn1);
  for (uint32_t k = 0; k < MAXITER; k++)
  {
    /* 2a. Compute residuals */
    /* sn+1 := \gradW[\gradsu - \epsp] */
    computeStress(extC, element->eps, strainPKn1, stressKn1);

    /* fn+1 = f[sn+1] */
    fKn1 = plasticYieldF(stressKn1, element->material->yieldStrengh);
    
    /* rn+1 = epn+1k - enp - \delta\gamma * \gradf[sn+1] */
    calcPlasticYieldFl(stressKn1, fl);
    stdmatrix_addmultiplyscalar(element->epsP, fl, delGamma, temp61);
    stdmatrix_subtract(strainPKn1, temp61, rKn1);
   
    /* 2b. Check convergence */
    if (stdmatrix_norm(rKn1) < TOL1 && fKn1 < TOL2) /* strain norm */
    {
      stdmatrix_copy(strainPKn1, element->epsP);

      /* calcs Nn1 and XI */
      /* XI = [Ckn1^-1 + \delta\gammak\grad2f(\sigmakn1)]^-1 */
      /* Nn1 = (XI*\gradf)/(sqrt(\gradfn1^T*XI*\gradfn1)) */
      calcPlasticYieldFl(stressKn1, fl);
      calcXIMatrix(stressKn1, extCInv, delGamma, fll, xiMatrix);
      
      stdmatrix_multiply(xiMatrix, fl, temp61);
      stdmatrix_multiplymTn(fl, temp61, temp11);
      stdmatrix_multiplyScalar(temp61, sqrt(stdmatrix_getelem(temp11, 0, 0))); /* temp61 == NN1 */
      stdmatrix_multiplymnT(temp61, temp61, temp66);
      
      /* ctm = XI - NN1*NN1T */
      stdmatrix_subtract(xiMatrix, temp66, element->Cep);

      break;
    }

    /* 2c. Compute consistent (algorithmic) tangent moduli */
    /* Ckn1 = grad2W[\gradsu - \epsp] */
    /* isotropic homogeneous.. C is constant */
    /* calcExtendedCMatrix(extC, element); */

    /* XI = [Ckn1^-1 + \delta\gammak\grad2f(\sigmakn1)]^-1 */
    calcXIMatrix(stressKn1, extCInv, delGamma, fll, xiMatrix);

    /* 2d. Compute kth increments */
    /* delta^2\gammak = (fkn1 + [\grad(fkn1)]^T * XIkn1*rkn1)/([\grad(fkn1)]^T * XIkn1 * \grad(fkn1)) */
    stdmatrix_multiplymTn(fl, xiMatrix, temp16);
    stdmatrix_multiply(temp16, rKn1, temp11);
    double num = stdmatrix_getelem(temp11, 0, 0) + fKn1;
    stdmatrix_multiply(temp16, fl, temp11);
    double den = stdmatrix_getelem(temp11, 0, 0);
    double del2Gamma = num / den;

    /* delta \epsilonPKn1 = Ck-1n1 * Xi */
    stdmatrix_subtractmultiplyscalar(del2Gamma, fl, rKn1, temp61);
    stdmatrix_multiply(extCInv, xiMatrix, temp66);
    stdmatrix_multiply(temp66, temp61, delEpsPk);

    /* 2e. Update plastic strains and consistency parameter */
    /* delta \gammak1 = delta\gammak + delta^2\gammak */
    delGamma = delGamma + del2Gamma;

    /* \epsPk1n1 = \epspkn1 + deltaEpspkn1 */
    stdmatrix_add(strainPKn1, delEpsPk, temp61);
    stdmatrix_copy(temp61, strainPKn1);

    /* k = k+1, go to 2a.*/
  }

  stdmatrix_copy(strainPKn1, element->epsP);

  stdmatrix_destroy(extC);
  stdmatrix_destroy(extCInv);
  stdmatrix_destroy(stressKn1);
  stdmatrix_destroy(strainPKn1);
  stdmatrix_destroy(rKn1);
  stdmatrix_destroy(temp11);
  stdmatrix_destroy(temp16);
  stdmatrix_destroy(temp61);
  stdmatrix_destroy(temp66);
  stdmatrix_destroy(fl);
  stdmatrix_destroy(fll);
  stdmatrix_destroy(xiMatrix);
  stdmatrix_destroy(delEpsPk);
  
  return 0;
}

int calcStressDeviator(const StdMatrix *stressVector, StdMatrix *r)
{
  if (stressVector == NULL || r == NULL)
  {
    fprintf(stderr, "calcStressDeviator: matrices not defined.\n");
    exit(-1);
  }

  StdMatrix *P;
  stdmatrix_create(&P, 6, 6);
  calcPMatrix(P);

  stdmatrix_multiply(P, stressVector, r);

  stdmatrix_destroy(P);

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

int incrementResultIntoNodeElast(NodeArray *nodes, Vector *x)
{
  uint32_t i;
  for (i = 0; i < nodes->len; i++)
  {
    nodes->ptr[i].ui[0] += vector_getelem(x, 2 * i);
    nodes->ptr[i].ui[1] += vector_getelem(x, 2 * i + 1);
  }
  return 0;
}

double plasticYieldF(const StdMatrix *stressVector, double flowStress)
{
  return sqrt(stdmatrix_voigtStressNormSqr(stressVector) - 1 / 3 * pow(stdmatrix_voigtTrace(stressVector), 2)) - flowStress;
}

int calcPlasticYieldFl(const StdMatrix *stressVector, StdMatrix *n)
{
  if (stressVector == NULL || n == NULL)
  {
    fprintf(stderr, "plasticYieldFl: matrices not defined.\n");
    exit(-1);
  }
  if (stdmatrix_ni(stressVector) != 6 || stdmatrix_nj(stressVector) != 1)
  {
    fprintf(stderr, "plasticYieldFl: incorrect matrix size.\n");
    exit(-1);
  }

  calcStressDeviator(stressVector, n);
  double norm = sqrt(stdmatrix_voigtStressNormSqr(n));
 
  double *ptr = stdmatrix_getptr(n);
  for (uint32_t i = 0; i < 6; i++)
  {
    ptr[i] = ptr[i] / norm;
  }

  return 0;
}

/* calcs numerically the second derivative of f in relation to the stress */
double calcPlasticYieldFll(const StdMatrix *stressVector, StdMatrix *fll)
{
  if (stressVector == NULL || fll == NULL)
  {
    fprintf(stderr, "plasticYieldFll: matrices not defined.\n");
    exit(-1);
  }
  if (stdmatrix_ni(fll) != 6 || stdmatrix_nj(fll) != 6)
  {
    fprintf(stderr, "plasticYieldFll: incorrect matrix size.\n");
    exit(-1);
  }

  StdMatrix *fl0;
  stdmatrix_create(&fl0, 6, 1);

  StdMatrix *fl;
  stdmatrix_create(&fl, 6, 1);

  StdMatrix *perturbedStress;
  stdmatrix_create(&perturbedStress, 6, 1);
  double *pStress = stdmatrix_getptr(perturbedStress);
  double *stress = stdmatrix_getptr(stressVector);
  

  calcPlasticYieldFl(stressVector, fl0);
  for (uint32_t j = 0; j < 6; j++)
  {
    for (uint32_t i = 0; i < 6; i++)
    {
      /* f(sig0, ... sig6) */
      pStress[i] = stress[i];
    }
    pStress[j] += EPS; /* f(sig0, ..., sigj+eps, ..., sig6) */
    calcPlasticYieldFl(perturbedStress, fl); /* fll(stress(sig0, ... , sigj+eps, ... sig6)) */

    for (uint32_t i = 0; i < 6; i++)
    {
      stdmatrix_setelem(fll, i, j, (stdmatrix_getelem(fl, i, 0) - stdmatrix_getelem(fl0, i, 0)) / EPS);
    }
  }
  stdmatrix_destroy(fl0);
  stdmatrix_destroy(fl);
  stdmatrix_destroy(perturbedStress);
  return 0;
}

int applyNeuBoundary1(Vector *b, const NeumannBCArray *neuBCArray, const NodeArray *nodes)
{
  Node n[2];
  for (uint32_t i = 0; i < neuBCArray->len; i++)
  {
    n[0] = nodes->ptr[neuBCArray->ptr[i].n0];
    n[1] = nodes->ptr[neuBCArray->ptr[i].n1];
    
    double sLength = getSideLength(n[0], n[1]);
    double val = neuBCArray->ptr[i].val[0] * sLength /2;
    
    vector_sumelem(b, n[0].index, val);
    vector_sumelem(b, n[1].index, val);
  }
  
  return 0;
}

int applyNeuBoundary2(Vector *b, const NeumannBCArray *neuBCArray, const NodeArray *nodes)
{
  double egdeNormal[2];
  for (uint32_t i = 0; i < neuBCArray->len; i++)
  {
    Node n[2];
    n[0] = nodes->ptr[neuBCArray->ptr[i].n0];
    n[1] = nodes->ptr[neuBCArray->ptr[i].n1];

    egdeNormal[0] = n[1].y - n[0].y;
    egdeNormal[1] = -(n[1].x - n[0].x);

    double xForce = neuBCArray->ptr[i].val[0] * egdeNormal[0] + neuBCArray->ptr[i].val[2] * egdeNormal[1];
    double yForce = neuBCArray->ptr[i].val[2] * egdeNormal[0] + neuBCArray->ptr[i].val[1] * egdeNormal[1];
    
    vector_sumelem(b, n[0].index * 2, xForce / 2);
    vector_sumelem(b, n[0].index * 2 + 1, yForce / 2);
    vector_sumelem(b, n[1].index * 2, xForce / 2);
    vector_sumelem(b, n[1].index * 2 + 1, yForce / 2);
  }
  
  return 0;
}

int applyDirBoundary(Matrix *stifMatrix, Vector *b, const DirichletBCArray *dirBCArray, const uint32_t ndof)
{
  uint32_t nodeIdx, dirIdx;
  double val;
  uint32_t i;
  for (i = 0; i < dirBCArray->len; i++)
  {
    nodeIdx = dirBCArray->ptr[i].nodeIndex;
    dirIdx = dirBCArray->ptr[i].dir;
    val = dirBCArray->ptr[i].val;
    matrix_applyDirBC(stifMatrix, b,  gIndex(ndof, nodeIdx, dirIdx), val);
  }
  
  return 0;
}


int conjugateGradientMethod(const MatrixCRS *a, Vector *x, const Vector *b)
{
  uint32_t n = vector_size(x);

  double alfa, rlOld, rlNew;
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
  
  /* rlOld = (r,r) */
  rlOld = vector_internalProduct(r, r);
  
  /* for (i = 0 to M) do */
  uint32_t i;
  for (i = 0; i < n; i++)
  { 
    /* z = Ap */
    matrixcrs_multiplyVector(a, p, temp);
    
    /* alfa = rlOld / (p,z) */
    alfa = rlOld / vector_internalProduct(p, temp);
    
    /* x = x + alfa p */
    vector_addVectorScalar(x, p, alfa, x);
    
    /* r = r - alfa z */
    vector_addVectorScalar(r, temp, -alfa, r);
    
    /* rlOld = (r,r) */
    rlNew = vector_internalProduct(r, r);
    
    /*   if (d < tolB) then exit loop */
    if (sqrt(rlNew) < TOL1) break;
    
    /* p = r + (rlNew/rlOld)p */
    vector_addVectorScalar(r, p, rlNew / rlOld, p);
    
    /* rlOld = rlNew */
    rlOld = rlNew;
  }
  
  vector_destroy(p);
  vector_destroy(r);
  vector_destroy(temp);
  
  return 0;
}

int change3to6Vector(const StdMatrix *m, StdMatrix *r)
{
  if (m == NULL || r == NULL)
  {
    fprintf(stderr, "change3to6vector: undefined matrices.\n");
    exit(-1);
  }
  if (stdmatrix_ni(m) != 3 || stdmatrix_ni(r) != 6 || stdmatrix_nj(m) != 1 || stdmatrix_nj(r) != 1)
  {
    fprintf(stderr, "change3to6vector: wrong sized matrices.\n");
    exit(-1);
  }
  stdmatrix_setelem(r, 0, 0, stdmatrix_getelem(m, 0, 0));
  stdmatrix_setelem(r, 1, 0, stdmatrix_getelem(m, 1, 0));
  stdmatrix_setelem(r, 3, 0, stdmatrix_getelem(m, 2, 0));
  return 0;
}

int change6to3Vector(const StdMatrix *m, StdMatrix *r)
{
  if (m == NULL || r == NULL)
  {
    fprintf(stderr, "change6to3Vector: undefined matrices.\n");
    exit(-1);
  }
  if (stdmatrix_ni(m) != 6 || stdmatrix_ni(r) != 3 || stdmatrix_nj(m) != 1 || stdmatrix_nj(r) != 1)
  {
    fprintf(stderr, "change6to3Vector: wrong sized matrices.\n");
    exit(-1);
  }
  stdmatrix_setelem(r, 0, 0, stdmatrix_getelem(m, 0, 0));
  stdmatrix_setelem(r, 1, 0, stdmatrix_getelem(m, 1, 0));
  stdmatrix_setelem(r, 2, 0, stdmatrix_getelem(m, 3, 0));
  return 0;
}

int change6to3Matrix(const StdMatrix *m, StdMatrix *r)
{
  if (m == NULL || r == NULL)
  {
    fprintf(stderr, "change6to3Matrix: undefined matrices.\n");
    exit(-1);
  }
  if (stdmatrix_ni(m) != 6 || stdmatrix_nj(m) != 6 ||
      stdmatrix_ni(r) != 3 || stdmatrix_nj(r) != 3)
  {
    fprintf(stderr, "change6to3Matrix: wrong sized matrices.\n");
    exit(-1);
  }

  stdmatrix_setelem(r, 0, 0, stdmatrix_getelem(m, 0, 0));
  stdmatrix_setelem(r, 0, 1, stdmatrix_getelem(m, 0, 1));
  stdmatrix_setelem(r, 0, 2, stdmatrix_getelem(m, 0, 3));
  stdmatrix_setelem(r, 1, 0, stdmatrix_getelem(m, 1, 0));
  stdmatrix_setelem(r, 1, 1, stdmatrix_getelem(m, 1, 1));
  stdmatrix_setelem(r, 1, 2, stdmatrix_getelem(m, 1, 3));
  stdmatrix_setelem(r, 2, 0, stdmatrix_getelem(m, 3, 0));
  stdmatrix_setelem(r, 2, 1, stdmatrix_getelem(m, 3, 1));
  stdmatrix_setelem(r, 2, 2, stdmatrix_getelem(m, 3, 3));
  return 0;
}

int updateIncrement(DirichletBCArray *dirBCArray, NeumannBCArray *neuBCArray, uint32_t curStep, uint32_t numLoadSteps)
{
  for (uint32_t i = 0; i < dirBCArray->len; i++)
  {
    dirBCArray->ptr[i].Dfk = dirBCArray->ptr[i].val / numLoadSteps;
    dirBCArray->ptr[i].fk = dirBCArray->ptr[i].Dfk * curStep;
  }
  for (uint32_t i = 0; i < neuBCArray->len; i++)
  {
    for (uint32_t k = 0; k < 3; k++)
    {
      neuBCArray->ptr[i].Dfk[k] = neuBCArray->ptr[i].val[k] / numLoadSteps;
      neuBCArray->ptr[i].fk[k] = neuBCArray->ptr[i].Dfk[k] * curStep;
    }
  }

  return 0;
}

int initializePlasticVariables(ElementArray *elements, NodeArray *nodes)
{
  for (uint32_t i = 0; i < nodes->len; i++)
  {
    nodes->ptr[i].DUI[0] = 0;
    nodes->ptr[i].DUI[1] = 0;
  }
  for (uint32_t i = 0; i < elements->len; i++)
  {
    Element *curElement = &elements->ptr[i];
    calcExtendedCMatrix(curElement->Cep, curElement);
  }
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