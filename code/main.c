#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include "matrix.h"
#include "helper.h"
#include "stdMatrix.h"

#define TOL 1e-6

struct node {
  uint32_t index;
  double x;
  double y;
  double u;
  double ui[2];
  uint32_t bm;
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

int computeGlobalMatrices(Matrix **a, Vector **b, const ElementArray *elements, const NodeArray *nodes, const DirichletBCArray *dirBCArray, const NeumannBCArray *neuBCArray, const uint32_t ndof);

int eleP1Poisson(Matrix *stifMatrix, Vector *b, const Element *element);
int eleP1Elast(Matrix *stifMatrix, Vector *b, const Element *element);
int computeStresses(const ElementArray *elements);

int loadPhiLFunctions(StdMatrix *phiL, Element *element);

int applyNeuBoundary1(Vector *b, const NeumannBCArray *neuBCArray, const NodeArray *nodes);
int applyNeuBoundary2(Vector *b, const NeumannBCArray *neuBCArray, const NodeArray *nodes);
int applyDirBoundary(Matrix *stifMatrix, Vector *b, const DirichletBCArray *dirBCArray, const uint32_t ndof);

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
  
  /* Poisson problem */
  
  /* (a) Input of data */
  /* (b) Representation of the triangulation Th. */
  readMaterials(&nMaterial, &materials, fileName);
  readNodes(&nodes, fileName);
  readElements(&elements, &nodes, materials, fileName);
  readBCP(&dirBCArray, &neuBCArray, &nodes, fileName);
  
  /* (c) Computation of the element stiffness matrices aK and element loads bK. */
  /* (d) Assembly of the global stiffness matrix A and load vector b. */
  Matrix *a;
  Vector *b;
  computeGlobalMatrices(&a, &b, &elements, &nodes, &dirBCArray, &neuBCArray, 1);
  
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
  free(dirBCArray.ptr);
  if (neuBCArray.len > 0)
  {
    free(neuBCArray.ptr);
  }
  matrixcrs_destroy(aCrs);
  vector_destroy(x);
  
  /* Now run the elastic problem */
  
  /* (a) Input of data */
  /* (b) Representation of the triangulation Th. */
  readBCE(&dirBCArray, &neuBCArray, &nodes, fileName);
  
  /* (c) Computation of the element stiffness matrices aK and element loads bK. */
  /* (d) Assembly of the global stiffness matrix A and load vector b. */
  computeGlobalMatrices(&a, &b, &elements, &nodes, &dirBCArray, &neuBCArray, 2);
  
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

int computeGlobalMatrices(Matrix **a, Vector **b, const ElementArray *elements, const NodeArray *nodes, const DirichletBCArray *dirBCArray, const NeumannBCArray *neuBCArray, const uint32_t ndof)
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
      fb = fbVal;
      /* fp = phiL[1][i/2]/det * pInt; */
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
  for (uint32_t i = 0; i < neuBCArray->len; i++)
  {
    Node n[2];
    n[0] = nodes->ptr[neuBCArray->ptr[i].n0];
    n[1] = nodes->ptr[neuBCArray->ptr[i].n1];
    double sLength = getSideLength(n[0], n[1]);
    
    vector_sumelem(b, n[0].index*2, neuBCArray->ptr[i].val[0]*sLength/2);
    vector_sumelem(b, n[0].index*2 + 1, neuBCArray->ptr[i].val[1]*sLength/2);
    vector_sumelem(b, n[1].index*2, neuBCArray->ptr[i].val[0]*sLength/2);
    vector_sumelem(b, n[1].index*2 + 1, neuBCArray->ptr[i].val[1]*sLength/2);
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