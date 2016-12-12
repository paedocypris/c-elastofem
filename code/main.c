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
};
typedef struct node Node;

struct nodeArray {
  uint32_t len;
  Node *ptr;
};
typedef struct nodeArray NodeArray;

struct element {
  uint32_t index;
  uint32_t material;
  uint32_t n0;
  uint32_t n1;
  uint32_t n2;
};
typedef struct element Element;

struct elementArray{
  uint32_t len;
  Element *ptr;
};
typedef struct elementArray ElementArray;

struct material {
  uint32_t index;
  double k;
  double fp;
  double l;
  double mu;
  double fg;
};
typedef struct material Material;

struct dirichletBC {
  uint32_t nodeIndex;
  uint32_t dir;
  double val;
};
typedef struct dirichletBC DirichletBC;

struct neumannBC {
  uint32_t elementIndex;
  uint32_t edge;
  double val;
};
typedef struct neumannBC NeumannBC;

int readNodes(NodeArray *nodes, char *fileName);
int readElements(ElementArray *elements, char *fileName);
int readMaterials(uint32_t *nMaterial, Material **materials, char *fileName);
int readBCP(uint32_t *nDirich, DirichletBC **dirBound, uint32_t *nNeumann, NeumannBC **neuBound, char *fileName);

int computeGlobalMatrices(Matrix **a, Vector **b, const ElementArray *elements, const NodeArray *nodes, const uint32_t nDirichlet, const DirichletBC *dirBound, const uint32_t nNeumann, const NeumannBC *neuBound, const Material *materials);

int eleP1Poisson(Matrix *stifMatrix, Vector *b, const Element *element, const NodeArray *nodes, const uint32_t nDirichlet, const DirichletBC *dirBound, const uint32_t nNeumann, const NeumannBC *neuBound, const Material *elemMaterial);
int eleP1Elast(Matrix *stifMatrix, Vector *b, const Element *element, const NodeArray *nodes, const uint32_t nDirichlet, const DirichletBC *dirBound, const uint32_t nNeumann, const NeumannBC *neuBound, const Material *elemMaterial);

int applyDirBoundary(Matrix *stifMatrix, Vector *b, const uint32_t nDirichlet, const DirichletBC *dirBound, const uint32_t ndof);

int conjugateGradientMethod(const MatrixCRS *a, Vector *x, const Vector *b);

int loadResultIntoNode(NodeArray *nodes, Vector *x);
static int64_t findElementInNeuArray(const uint32_t elementIndex, const uint32_t nNeumann, const NeumannBC *neuBound);

int printVTK(NodeArray *nodes, ElementArray *elements, char *fileName);

inline uint32_t gi(uint32_t nj, uint32_t i, uint32_t j);


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
  
  /* (a) Input of data */
  /* (b) Representation of the triangulation Th. */
  readNodes(&nodes, fileName);
  readElements(&elements, fileName);
  readBCP(&nDirichlet, &dirBound, &nNeumann, &neuBound, fileName);
  readMaterials(&nMaterial, &materials, fileName);
  
  /* (c) Computation of the element stiffness matrices aK and element loads bK. */
  /* (d) Assembly of the global stiffness matrix A and load vector b. */
  Matrix *a;
  Vector *b;
  computeGlobalMatrices(&a, &b, &elements, &nodes, nDirichlet, dirBound, nNeumann, neuBound, materials);
  matrix_print(a);
  vector_print(b);
  
  /* (e) Solution of the system of equations Ax = b. */
  MatrixCRS *aCrs; 
  matrixcrs_create(a, &aCrs);
  Vector *x;
  vector_createzero(&x, nodes.len);
  conjugateGradientMethod(aCrs, x, b);
  loadResultIntoNode(&nodes, x);
  
  /* (f) Presentation of result. */
  printVTK(&nodes, &elements, fileName);
  
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
    fprintf(fp, "3 %"PRIu32" %"PRIu32" %"PRIu32"\n", elements->ptr[i].n0, elements->ptr[i].n1, elements->ptr[i].n2);
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
    fprintf(fp, "%.4g\n", nodes->ptr[i].u);
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

int readElements(ElementArray *elements, char *fileName)
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
    (elements->ptr + i)->n0 = n0;
    (elements->ptr + i)->n1 = n1;
    (elements->ptr + i)->n2 = n2;
    (elements->ptr + i)->material = matIdx;
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
  NeumannBC *neuBC = malloc(n * sizeof(NeumannBC));
  uint32_t edgeNumber;
  for (i = 0; i < n; i++)
  {
    fscanf(fp, "%" SCNu32" %"SCNu32" %lf*\n", &bIdx, &edgeNumber, &value);
    (*(neuBC + i)).elementIndex = bIdx;
    (*(neuBC + i)).edge = edgeNumber;
    (*(neuBC + i)).val = value;
  }
  *nNeumann = n;
  *neuBound = neuBC;
  
  fclose(fp);
  
  return 0;
}

int computeGlobalMatrices(Matrix **a, Vector **b, const ElementArray *elements, const NodeArray *nodes, const uint32_t nDirichlet, const DirichletBC *dirBound, const uint32_t nNeumann, const NeumannBC *neuBound, const Material *materials)
{
  Matrix *newStif;
  matrix_create(&newStif, nodes->len, nodes->len);
  Vector *newB;
  vector_createzero(&newB, nodes->len);
  
  uint32_t i;
  for (i = 0; i < elements->len; i++)
  {
    eleP1Poisson(newStif, newB, &elements->ptr[i], nodes, nDirichlet, dirBound, nNeumann, neuBound, &materials[elements->ptr[i].material]);
  }
  
  /* apply dir boundary */
  applyDirBoundary(newStif, newB, nDirichlet, dirBound, 1);
  
  *a = newStif;
  *b = newB;
  return 0;
}

int eleP1Poisson(Matrix *stifMatrix, Vector *b, const Element *element, const NodeArray *nodes, const uint32_t nDirichlet, const DirichletBC *dirBound, const uint32_t nNeumann, const NeumannBC *neuBound, const Material *elemMaterial)
{ 
  uint32_t i, j;

  const Node *n[3];
  n[0] = &nodes->ptr[element->n0];
  n[1] = &nodes->ptr[element->n1];
  n[2] = &nodes->ptr[element->n2];
  
  /* calc auxiliar parameters */
  double det, den, area;
  det = n[0]->x*n[1]->y + n[0]->y*n[2]->x + n[1]->x*n[2]->y - n[1]->y*n[2]->x - n[2]->y*n[0]->x - n[0]->y*n[1]->x;
  den = det*2;
  area = det/2;
  
  double ax[3], bx[3];
  ax[0] = n[1]->y - n[2]->y;
  ax[1] = n[2]->y - n[0]->y;
  ax[2] = n[0]->y - n[1]->y;
  bx[0] = n[2]->x - n[1]->x;
  bx[1] = n[0]->x - n[2]->x;
  bx[2] = n[1]->x - n[0]->x;
  
  
  double k[3][3];
  for (i = 0; i < 3; i++)
  {
    for(j = 0; j < 3; j++)
    {
      k[i][j] = elemMaterial->k * (ax[i]*ax[j] + bx[i]*bx
      [j])/den;
    }
  }
 
  double fb = elemMaterial->fp * area;
  /* find neumann boundary value */
  double gb[3] = {0};
  int64_t nIdx = findElementInNeuArray(element->index, nNeumann, neuBound);
  if (nIdx != -1)
  {
    switch (neuBound[nIdx].edge)
    {
      case 0:
        gb[1] = gb[2] = sqrt(pow(n[2]->x - n[1]->x,2) + pow(n[2]->y - n[1]->y,2)) * neuBound->val/2;
        break;
      case 1:
        gb[0] = gb[2] = sqrt(pow(n[0]->x - n[2]->x,2) + pow(n[0]->y - n[2]->y,2)) * neuBound->val/2;
        break;
      case 2:
        gb[1] = gb[2] = sqrt(pow(n[1]->x - n[0]->x,2) + pow(n[1]->y - n[0]->y,2)) * neuBound->val/2;
        break;
      default:
        fprintf(stderr, "eleP1Poisson: lado incorreto: %"PRIu32".\n", neuBound[nIdx].edge);
        break;
    }
  }
  
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      matrix_sumelem(stifMatrix, n[i]->index, n[j]->index, k[i][j]);
    }
    vector_sumelem(b, n[i]->index, fb + gb[0]);
  }
  
  return 0;
}

int eleP1Elast(Matrix *stifMatrix, Vector *b, const Element *element, const NodeArray *nodes, const uint32_t nDirichlet, const DirichletBC *dirBound, const uint32_t nNeumann, const NeumannBC *neuBound, const Material *elemMaterial)
{
  double a0, a1, a2;
  double b0, b1, b2;
  double c, d;
  double det, den, area;
  
  const Node *n[3];
  n[0] = &nodes->ptr[element->n0];
  n[1] = &nodes->ptr[element->n1];
  n[2] = &nodes->ptr[element->n2];
  
  /* calc auxiliar parameters */
  det = n[0]->x*n[1]->y + n[0]->y*n[2]->x + n[1]->x*n[2]->y - n[1]->y*n[2]->x - n[2]->y*n[0]->x - n[0]->y*n[1]->x;
  den = det*2;
  area = det/2;
  a0 = n[1]->y - n[2]->y;
  a1 = n[2]->y - n[0]->y;
  a2 = n[0]->y - n[1]->y;
  b0 = n[2]->x - n[1]->x;
  b1 = n[0]->x - n[2]->x;
  b2 = n[1]->x - n[0]->x;
  c = elemMaterial->l+2*elemMaterial->mu;
  d = elemMaterial->l;
  
  /* 
  [ a0  0   a1  0   a2  0  ]
  [ 0   b0  0   b1  0   b2 ]
  [ b0  a0  b1  a1  b2  a2 ]
  */
  StdMatrix *bMatrix;
  stdmatrix_createzero(&bMatrix, 3, 6);
  stdmatrix_setelem(bMatrix, 0, 0, a0);
  stdmatrix_setelem(bMatrix, 0, 2, a1);
  stdmatrix_setelem(bMatrix, 0, 4, a2);
  stdmatrix_setelem(bMatrix, 1, 1, b0);
  stdmatrix_setelem(bMatrix, 1, 3, b1);
  stdmatrix_setelem(bMatrix, 1, 5, b2);
  stdmatrix_setelem(bMatrix, 2, 0, b0);
  stdmatrix_setelem(bMatrix, 2, 1, a0);
  stdmatrix_setelem(bMatrix, 2, 2, b1);
  stdmatrix_setelem(bMatrix, 2, 3, a1);
  stdmatrix_setelem(bMatrix, 2, 4, b2);
  stdmatrix_setelem(bMatrix, 2, 5, a2);
  
  /*
  [ c  d  0 ]
  [ d  c  0 ]
  [ 0  0  d ]
  */
  StdMatrix *cMatrix;
  stdmatrix_createzero(&cMatrix, 3, 3);
  stdmatrix_setelem(cMatrix, 0, 0, c);
  stdmatrix_setelem(cMatrix, 0, 1, d);
  stdmatrix_setelem(cMatrix, 1, 0, d);
  stdmatrix_setelem(cMatrix, 1, 1, c);
  stdmatrix_setelem(cMatrix, 2, 2, d);
  
  StdMatrix *temp;
  stdmatrix_create(&temp, 6, 3);
  StdMatrix *kMatrix;
  stdmatrix_create(&kMatrix, 6, 6);
  
  stdmatrix_multiplymTn(bMatrix, cMatrix, temp);
  stdmatrix_multiply(temp, bMatrix, kMatrix);
  stdmatrix_destroy(temp);
  
  double fb = elemMaterial->fg * area;
  
  uint32_t i,j;
  for (i = 0; i < 6; i++)
  {
    for (j = 0; j < 6; j++)
    {
      matrix_sumelem(stifMatrix, n[i]->index, n[j]->index, stdmatrix_getelem(kMatrix, i, j));
    }
    vector_sumelem(b, n[i]->index, fb);
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

int loadResultIntoNode(NodeArray *nodes, Vector *x)
{
  uint32_t i;
  for (i = 0; i < nodes->len; i++)
  {
    nodes->ptr[i].u = vector_getelem(x, i);
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
    matrix_applyDirBC(stifMatrix, b, gi(ndof, nodeIdx, dirIdx), val);
  }
  return 0;
}

static int64_t findElementInNeuArray(const uint32_t elementIndex, const uint32_t nNeumann, const NeumannBC *neuBound)
{
  int64_t i;
  for (i = 0; i < nNeumann; i++)
  {
    if (neuBound[i].elementIndex == elementIndex)
      return i;
  }
  return -1;
}

inline uint32_t gi(uint32_t nj, uint32_t i, uint32_t j)
{
  return i*nj + j;
}