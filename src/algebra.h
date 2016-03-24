#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <stdio.h>
#include <stdlib.h>


/*#ifndef WIN32
#define WIN32
  typedef unsigned int uint;
#endif*/

#ifdef WIN32
#define uint unsigned int
#endif

typedef struct{
  double **data;
  uint row, col;
}matrix;

typedef struct{
  double *data;
  uint size;
} dvector;

typedef struct{
  uint *data;
  uint size;
} uivector;

void NewMatrix(uint, uint, matrix**);
void DelMatrix(matrix**);

void MatrixVectorDotProduct(matrix *m, dvector *v, dvector **r);
void MatrixDotProduct(matrix *m_t, matrix *m, matrix **r);
void MatrixTranspose(matrix *m, matrix *r);
void MatrixInversion(matrix *A, matrix **Y);
/*void MatrixInversion(double **A, int order, double **Y); fa cagare
int GetMinor(double **src, double **dest, int row, int col, int order);
double CalcDeterminant( double **mat, int order); 
*/
void IdentityMatrix(matrix **);




void NewVector(uint, dvector**);
void DelVector(dvector**);
void VectorAppend(dvector**, double);

void NewUintVector(uint, uivector**);
void DelUintVector(uivector**);
void VectorUintAppend(uivector**, uint);

void SolveEquation(matrix *, dvector **result);

uint Factorial(uint x);
#endif