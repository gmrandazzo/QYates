#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "algebra.h"

//memory wrapper
static inline void *xmalloc(size_t size)
{
  void *ptr = malloc(size);
  if (ptr == NULL) {
    fprintf(stderr, "Memory Exhausted!\n");
    abort();
  }
  return ptr;
}

// matrix inversioon
// the result is put in Y
void NewMatrix(uint row_ , uint col_, matrix **m )
{
  (*m) = xmalloc(sizeof(matrix));
  int i, j;
  (*m)->row = row_;
  (*m)->col = col_;
  (*m)->data = xmalloc(sizeof(double)*row_);
  for(i = 0; i < row_; i++){
    (*m)->data[i] = xmalloc(sizeof(double)*col_);
    for(j = 0; j < col_; j++)
      (*m)->data[i][j] = 0.f;
  }
}

void DelMatrix(matrix **m)
{
  int i;
  for(i = 0; i < (*m)->row; i++)
    free((*m)->data[i]);
  free((*m)->data);
  free((*m));
}

void NewVector(uint size_, dvector **v )
{
  int i;
  (*v) = xmalloc(sizeof(dvector));
  (*v)->size = size_; 
  (*v)->data = xmalloc(sizeof(double)*size_);
  for(i = 0; i < (*v)->size; i++)
    (*v)->data[i] = 0.f;
}

void DelVector(dvector **v )
{
  free((*v)->data);
  free((*v));
}

void VectorAppend(dvector** v, double val)
{
  uint i;
  dvector *tmp;
  NewVector((*v)->size+1, &tmp);
  
  for(i = 0; i < (*v)->size; i++)
    tmp->data[i] = (*v)->data[i];
  tmp->data[tmp->size-1] = val;

  DelVector(v);
  NewVector(tmp->size, v);
  for(i = 0; i < tmp->size; i++)
    (*v)->data[i] = tmp->data[i];
  
  DelVector(&tmp);
}


void NewUintVector(uint size_, uivector **v)
{
  int i;
  (*v) = xmalloc(sizeof(uivector));
  (*v)->size = size_; 
  (*v)->data = xmalloc(sizeof(uint)*size_);
  for(i = 0; i < (*v)->size; i++)
    (*v)->data[i] = 0;
}

void DelUintVector(uivector **v)
{
  free((*v)->data);
  free((*v));
}

void VectorUintAppend(uivector** v, uint val )
{
  uint i;
  uivector *tmp;
  NewUintVector((*v)->size+1, &tmp);
  
  for(i = 0; i < (*v)->size; i++)
    tmp->data[i] = (*v)->data[i];
  tmp->data[tmp->size-1] = val;

  DelUintVector(v);
  NewUintVector(tmp->size, v);
  for(i = 0; i < tmp->size; i++)
    (*v)->data[i] = tmp->data[i];
  
  DelUintVector(&tmp);
}


/*
 * v is a vector multiplied for the matrix m that must have the same column as the
 * lenght of the vector.
 * r is the result vector cfrom
 */

void MatrixVectorDotProduct(matrix *m, dvector *v, dvector **r)
{
  int i, j;
 /* RM = numero di righe della matrice */
   /* CM = numero di colonne della matrice (uguale al numero di righe del vettore) */
   /* M = matrice [RM] × [CM] */
   /* V = vettore [CM]        */
   /* il vettore risultato sarà VR [RM] con stesso numero di righe della matrice */
   /* poniamo VR [RM] sia inizialmente posto con tutti i valori a zero. */
   for (i = 0; i < m->row; i++)       /* scandisco le righe con l'indice i */
       for (j = 0; j < m->col; j++)     /* e le colonne con j  */
           (*r)->data[i] = (*r)->data[i] + m->data[i][j] * v->data[j];
}

/*
 * m_t is the transposed matrix of m. The result is a square matrix named r.
 */
void MatrixDotProduct(matrix *m_t, matrix *m, matrix **r)
{
  int i, j, k, l, q;
//   (*r).row =(*r).col = m.col; // square matrix
  for(i = 0; i < (*r)->row; i++){
    for(j = 0; j < m_t->row; j++){
      for(k = 0; k < m_t->col; k++){
        (*r)->data[j][i] += m_t->data[j][k] * m->data[k][i];
      }
    }
  }
}

void MatrixTranspose(matrix *m, matrix *r){
  int i, j;
  for(i = 0; i < m->row; i++)
    for(j = 0; j < m->col; j++)
      (*r).data[j][i] = m->data[i][j];
}

/*Gauss-Jordan
 *
 * 1) If the first row have the first element 0, then excange this with an other row that have the first element != 0. If all the row have the first element 0 go to the step 3.
 * 2) For eac row A_i with the first element != 0, except the first row considered, multiply the first row for a coefficient c that must 
 * Per ogni riga Ai con primo elemento non nullo, eccetto la prima (i > 1), moltiplica la prima riga per un coefficiente scelto in maniera tale che la somma tra la prima riga e Ai abbia il primo elemento nullo (quindi coefficiente = − Ai1 / A11). Sostituisci Ai con la somma appena ricavata.
 */
void MatrixInversion(matrix *A, matrix **Y)
{
  matrix *A_;
  memcpy(&A_, &A, sizeof(matrix));
  uint i, j, k, l;
  IdentityMatrix(Y); // fill the identity matrix
  
  /*for(i = 0; i < (*Y)->row; i++){
    for(j = 0; j < (*Y)->col; j++)
      printf("%f\t", (*Y)->data[i][j]);
    printf("\n");
  }
  
  for(i = 0; i < A_->row; i++){
    for(j = 0; j < A_->col; j++)
      printf("%f\t", A_->data[i][j]);
    printf("\n");
  }
  */
  double c; // this is the multiplier constant for make the step 2
  for(i = 0 ; i < A_->row; i++){
    if(A_->data[i][i] != 0 ){
      /*printf("pivot %f\n", A_->data[i][i]);*/ 
      for(k = 0; k < A_->row; k++){
        if(k != i){
          c = -1 * A_->data[k][i]/A_->data[i][i];

          dvector *tmprow;
          NewVector(A_->col*2, &tmprow);
          for(j = 0; j < A_->col; j++){
            tmprow->data[j] = (c*A_->data[i][j]) + A_->data[k][j];
          }
          
          for(j = A_->col; j < 2*A_->col; j++){
            tmprow->data[j] = (c*(*Y)->data[i][j-A_->col]) + (*Y)->data[k][j-A_->col];
          }
          
          for(j = 0; j < A_->col; j++){
            A_->data[k][j] = tmprow->data[j];
          }
          
          for(j = A_->col; j < 2*A_->col; j++){
            (*Y)->data[k][j-A_->col] = tmprow->data[j];
          }

          DelVector(&tmprow);
        }
        /*
        printf(">>>>>>>>>> STEP %d.%d\n",i, k);
        for(l = 0; l < A_->row; l++){
          for(j = 0; j < A_->col; j++)
            printf("%f\t", A_->data[l][j]);
          printf("\n");
        }
        
        printf("\n IDENTITY \n");
        for(l = 0; l < (*Y)->row; l++){
          for(j = 0; j < (*Y)->col; j++)
            printf("%f\t", (*Y)->data[l][j]);
          printf("\n");
        }
        */
      }
    }
  }
  
  for(i = 0; i < (*Y)->row; i++){
    for(j = 0; j < (*Y)->row; j++){
      (*Y)->data[i][j] = (*Y)->data[i][j] / A->data[i][i];
    }
  }
  
  /*printf("Inverted Matrix\n");
  for(i = 0; i < (*Y)->row; i++){
    for(j = 0; j < (*Y)->row; j++)
      printf("%f\t", (*Y)->data[i][j]);
    printf("\n");
  }*/ 
  
}

void IdentityMatrix(matrix **m)
{
  uint i;
  if((*m)->row == (*m)->col)
    for(i=0; i < (*m)->row; i++)
      (*m)->data[i][i]=1;
}

/*
void MatrixInversion(double **A, int order, double **Y)
{
  int i, j;
  // get the determinant of a
  double det = 1.0/CalcDeterminant(A,order);

  // memory allocation
  double *temp = malloc(sizeof(double)*((order-1)*(order-1)));
  double **minor = malloc(sizeof(double)*(order-1));
  for(i=0;i<order-1;i++)
    minor[i] = temp+(i*(order-1));

    for(j=0;j<order;j++){
      for(i=0;i<order;i++){
      // get the co-factor (matrix) of A(j,i)
        GetMinor(A,minor,j,i,order);
        Y[i][j] = det*CalcDeterminant(minor,order-1);
        if( (i+j)%2 == 1)
          Y[i][j] = -Y[i][j];
      }
    }
    // release memory
    free(temp);
    free(minor);
}

// calculate the cofactor of element (row,col)
int GetMinor(double **src, double **dest, int row, int col, int order)
{
  int i, j;
  // indicate which col and row is being copied to dest
  int colCount=0,rowCount=0;

  for(i = 0; i < order; i++ ){
    if( i != row ){
      colCount = 0;
      for(j = 0; j < order; j++ ){
      // when j is not the element
        if( j != col ){
          dest[rowCount][colCount] = src[i][j];
          colCount++;
        }
      }
      rowCount++;
    }
  }
  return 1;
}

// Calculate the determinant recursively.
double CalcDeterminant( double **mat, int order)
{
  int i;
  // order must be >= 0
  // stop the recursion when matrix is a single element
  if( order == 1 )
    return mat[0][0];

  // the determinant value
  double det = 0;

  // allocate the cofactor matrix
  double **minor;
  minor = malloc(sizeof(double*)*(order-1));
  for(i=0;i<order-1;i++)
    minor[i] = malloc(sizeof(double)*(order-1));

  for(i = 0; i < order; i++ ){
    // get minor of element (0,i)
    GetMinor( mat, minor, 0, i , order);
    // the recusion is here!

    det += (i%2==1?-1.0:1.0) * mat[0][i] * CalcDeterminant(minor,order-1);
    //det += pow( -1.0, i ) * mat[0][i] * CalcDeterminant( minor,order-1 );
  }

  // release memory
  for(i=0;i<order-1;i++)
    free(minor[i]);
  free(minor);
  return det;
}
*/

/* limitato a 3 equazioni */
void SolveEquation(matrix *equation, dvector **solution)
{
  uint i, j, eq;
  double a;
  /*
   Problem: Solve the following system:
      x + y + z  = 4
      x - 2y - z = 1
      2x - y - 2z = -1

    Start out by multiplying the
    first equation by -1 and add
    it to the second equation to 
    eliminate x from the second
    (*equation).

    -x  - y - z = -4
      x - 2y - z = 1
    ----------------
        -3y - 2z = -3

    Now eliminate x from the third
    equation by multiplying the first
    equation by -2 and add it to
    the third (*equation).

    -2x - 2y - 2z = -8
      2x -  y - 2z = -1
    ------------------
          -3y - 4z = -9

    Next, eliminate y from the third
    equation by multiplying the second
    equation by -1 and adding it to
    the third (*equation).

      3y +  2z = 3 
    -3y -  4z = -9
    --------------
          -2z = -6

    Solve the third equation for z.
          
    -2z = -6
      z = 3

    Substitute 3 for z in the
    second equation and solve for y.

    -3y - 2z = -3
    -3y - 2(3) = -3
    -3y - 6 = -3
    -3y = 3
    y = -1

    Lastly, substitute -1 for y and
    3 for z in the first equation
    and solve for x.

    x + (-1) + 3 = 4
    x + 2 = 4
    x = 2

    The answer is (2, -1, 3).

  */
  
  for(eq = 0; eq < (*equation).row; eq++){
    for(i = eq+1; i < 3; i++){
      a = (*equation).data[i][eq] / (*equation).data[eq][eq];
      for(j = 0; j < 4; j++){
        (*equation).data[i][j] = (*equation).data[eq][j] *( -a) + (*equation).data[i][j];
      }
    }    
  }

  for(eq = (*equation).row-1; eq < -1; eq--){
    if(eq == (*equation).row-1){
      (*solution)->data[eq] = (*equation).data[eq][(*equation).col-1] / (*equation).data[eq][eq];
    }
    else{
      double k = 0.f;
      for(i = eq+1; i < (*equation).col-1; i++){
        if((*solution)->data[i] != 0){
          k += (*solution)->data[i]*(*equation).data[eq][i];
        }
      }
      (*solution)->data[eq] = ((*equation).data[eq][3] - k) / (*equation).data[eq][eq];
    } 
  }
  
  /* the solution for 3 equation in 3 unknowns variable
  (*solution)->data[2] = (*equation).data[2][3] / (*equation).data[2][2];
  (*solution)->data[1] = ((*equation).data[1][3] - ((*equation).data[1][2] * (*solution)->data[2])) / (*equation).data[1][1];
  (*solution)->data[0] = ((*equation).data[0][3] - (*equation).data[0][2]*(*solution)->data[2] - (*equation).data[0][1]*(*solution)->data[1]) / (*equation).data[0][0];
  */
}

uint Factorial(uint x)
{
  uint f;
  uint i;
  f = x;
  i = x;

  while(i > 0){
    f *= i;
    i--;
  }
  return f;
}

