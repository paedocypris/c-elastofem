#ifndef HELPER_H_INCLUDED
#define HELPER_H_INCLUDED

#include <float.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

inline int floatIsZero(double x)
{
  if (fabs(x) < DBL_EPSILON) return 1;
  return 0;
}

inline void *sMalloc(size_t memsize)
{
  void *allocMem = malloc(memsize);
  if (!allocMem && memsize)
  {
    fprintf(stderr, "sMalloc: Could not allocate memory!\n");
    exit(-1);
  }
  return allocMem;
}

inline void *sRealloc(void *block, size_t memsize)
{
  void *rBuf = realloc(block, memsize);
  if (rBuf == NULL && memsize)
  {
    fprintf(stderr, "sMalloc: Could not (re)allocate memory!\n");
    exit(-1);
  }
  return rBuf;
}

#endif
