#include "algebra.h"

void test1()
{
  int i, j;
  double s[16][10] = {{1, -1, -1, -1, 1, 1, 1, 1, 1, 1},
{1, 1, -1, -1, 1, 1, 1, -1, -1, 1},
{1, -1, 1, -1, 1, 1, 1, -1, 1, -1},
{1, 1, 1, -1, 1, 1, 1, 1, -1, -1},
{1, -1, -1, 1, 1, 1, 1, 1, -1, -1},
{1, 1, -1, 1, 1, 1, 1, -1, 1, -1},
{1, -1, 1, 1, 1, 1, 1, -1, -1, 1},
{1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
{1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, -2, 0, 0, 4, 0, 0, 0, 0, 0},
{1, 2, 0, 0, 4, 0, 0, 0, 0, 0},
{1, 0, -2, 0, 0, 4, 0, 0, 0, 0}, 
{1, 0, 2, 0, 0, 4, 0, 0, 0, 0},
{1, 0, 0, -2, 0, 0, 4, 0, 0, 0},
{1, 0, 0, 2, 0, 0, 4, 0, 0, 0}};

  double v_[16] = {25.74,
                  48.98,
                  42.78,
                  35.94,
                  41.50,
                  50.10,
                  46.06,
                  27.70,
                  57.52,
                  59.68,
                  35.50,
                  44.18,
                  35.58,
                  28.46,
                  33.50,
                  42.02};
  
 
  dvector *v;
  NewVector(16, &v);
  
  for(i = 0; i < 16; i++)
    v->data[i] = v_[i];
  
  
  matrix *m;
  matrix *r;
  NewMatrix(16,10, &m);
  for(i = 0; i < 16; i++)
    for(j = 0; j < 10; j++)
      m->data[i][j] = s[i][j];

  NewMatrix(m->col, m->row, &r);
  MatrixTranspose(m, r);
  
 /* printf("transposed matrix\n");
  for(i = 0; i<r->row; i++){
    for(j = 0; j < r->col; j++)
      printf("%f  ", r->data[i][j]);
    printf("\n");
  }
  
  */
  matrix *dprod;
  NewMatrix(10,10,&dprod);
  MatrixDotProduct(r, m, &dprod);
  
  /*
  printf("Dot Product\n");
  for(i = 0; i<dprod->row; i++){
    for(j = 0; j < dprod->col; j++)
      printf("%f  ", dprod->data[i][j]);
    printf("\n");
  }
  */
  
  matrix *inv;
  NewMatrix(dprod->row, dprod->col, &inv);
  MatrixInversion(dprod, &inv);
  
  dvector *d;
  NewVector(10, &d);
  MatrixVectorDotProduct(r, v, &d);

  
  dvector *f;
  NewVector(10, &f);
  MatrixVectorDotProduct(inv, d, &f);
  
  
//   printf("Vettore Z'y\n");
//   for(i = 0; i < 10; i++)
//     printf("%f\t", f->data[i]);
//   printf("\n");
  
  printf("final\n");
  for(i = 0; i < 10; i++)
    printf("%f\n", f->data[i]);

  
  DelMatrix(&inv);
  DelMatrix(&dprod);
  DelMatrix(&r);
  DelMatrix(&m);
  
  DelVector(&v);
  DelVector(&d);
  DelVector(&f);
}

void test2(){
  uint i, j;
  matrix *m1;
  matrix *m2;
  NewMatrix(2,3, &m1);
  NewMatrix(3,2, &m2);
  
  m1->data[0][0] = 2; m1->data[0][1] = 1;  m1->data[0][2] = 3; 
  m1->data[1][0] = -2; m1->data[1][1] = 2; m1->data[1][2] = 1; 
  
  m2->data[0][0] = 2; m2->data[0][1] = 1;
  m2->data[1][0] = 3; m2->data[1][1] = 2;
  m2->data[2][0] = -2; m2->data[2][1] = 2;
  
  matrix *r;
  
  NewMatrix(m1->row, m1->row, &r);
  MatrixDotProduct(m1, m2, &r);
  
  for(i = 0; i < r->row; i++){
    for(j = 0; j < r->col; j++)
      printf("%f\t",r->data[i][j]);
    printf("\n");
  }
  DelMatrix(&m1);
  DelMatrix(&m2);
  DelMatrix(&r);
}

void test3()
{
  int i;
  matrix *m;
  dvector *v;
  dvector *r;
  NewMatrix(2, 3, &m);
  NewVector(3, &v);
  NewVector(2, &r);
  
  m->data[0][0] = 1; m->data[0][1] = 2; m->data[0][2]= 0;
  m->data[1][0] = 3; m->data[1][1] = -1; m->data[1][2]= 4;
  
  v->data[0] = 1; v->data[1] = 0; v->data[2] = -1; 
  
  MatrixVectorDotProduct(m, v, &r);

  for(i = 0; i < 2; i++)
    printf("%f\n", r->data[i]);
  
  DelVector(&r);
  DelVector(&v);
  DelMatrix(&m);
}

void test4()
{
  printf("Test 4\n");
  int i;
  dvector *v;
  matrix *equation;
  
  NewMatrix(3,4, &equation);
  NewVector(3,&v);
  
  /*
  * results of this equation is x = 2, y = -1 , z = 3
  */
  equation->data[0][0] = 1; equation->data[0][1] = 1; equation->data[0][2] = 1; equation->data[0][3] = 4;
  equation->data[1][0] = 1; equation->data[1][1] = -2; equation->data[1][2] = -1; equation->data[1][3] = 1;
  equation->data[2][0] = 2; equation->data[2][1] = -1; equation->data[2][2] = -2; equation->data[2][3] = -1;
  
  SolveEquation(equation, &v);
  
  for(i = 0; i < v->size; i++)
    printf("%f\n", v->data[i]);
  
  DelVector(&v);
  DelMatrix(&equation);
}

void test5()
{
  printf("Test 5\n");
  int i;
  dvector *v;
  matrix *equation;
  
  NewMatrix(3,4, &equation);
  NewVector(3,&v);

  /*
  * results of this equation is x = 1, y = -2 , z = -2
  */
  equation->data[0][0] = 3; equation->data[0][1] = 2; equation->data[0][2] = -1; equation->data[0][3] = 1;
  equation->data[1][0] = 2; equation->data[1][1] = -2; equation->data[1][2] = 4; equation->data[1][3] = -2;
  equation->data[2][0] = -1; equation->data[2][1] = 0.5; equation->data[2][2] = -1; equation->data[2][3] = -0;
  
  SolveEquation(equation, &v);
  
  for(i = 0; i < v->size; i++)
    printf("%f\n", v->data[i]);
  
  DelVector(&v);
  DelMatrix(&equation);
}


void test6()
{
  printf("Test 6\n");
  int i;
  dvector *v;
  matrix *equation;
  
  NewMatrix(3,4, &equation);
  NewVector(3,&v);
  
  equation->data[0][0] = 9.38; equation->data[0][1] = 7.13; equation->data[0][2] = 3.27; equation->data[0][3] = 1.50;
  equation->data[1][0] = 7.13; equation->data[1][1] = 12.54; equation->data[1][2] = 2.73; equation->data[1][3] = -2.13;
  equation->data[2][0] = 3.27; equation->data[2][1] = 2.73; equation->data[2][2] = 10.42; equation->data[2][3] = 1.81;
  
  SolveEquation(equation, &v);
  
  for(i = 0; i < v->size; i++)
    printf("%f\n", v->data[i]);
  
  DelVector(&v);
  DelMatrix(&equation);
}

void test7()
{
  uint i, j;
  printf("Test 7\n");
  matrix *m1;
  matrix *m2;
    
  NewMatrix(3, 3, &m1);
  NewMatrix(3, 3, &m2);
  
  m1->data[0][0] = 1; m1->data[0][1] = 0; m1->data[0][2]= 3;
  m1->data[1][0] = 4; m1->data[1][1] = 2; m1->data[1][2]= -1;
  m1->data[2][0] = -2; m1->data[2][1] = 0; m1->data[2][2]= 2;
  
  MatrixInversion(m1, &m2);
  
  printf("Inverted Matrix\n");
  for(i = 0; i < m2->row; i++){
    for(j = 0; j < m2->col; j++)
      printf("%f\t", m2->data[i][j]);
    printf("\n");
  }
  
  DelMatrix(&m1);
  DelMatrix(&m2);

}

void test8()
{
  uint i, j;
  printf("Test 8\n");
  matrix *m1;
  matrix *m2;
    
  NewMatrix(3, 3, &m1);
  NewMatrix(3, 3, &m2);
  
  m1->data[0][0] = 1; m1->data[0][1] = 4; m1->data[0][2]= 8;
  m1->data[1][0] = 2; m1->data[1][1] = 1; m1->data[1][2]= 3;
  m1->data[2][0] = 6; m1->data[2][1] = 9; m1->data[2][2]= 4;
  
  MatrixInversion(m1, &m2);
  
  printf("Inverted Matrix\n");
  for(i = 0; i < m2->row; i++){
    for(j = 0; j < m2->col; j++)
      printf("%f\t", m2->data[i][j]);
    printf("\n");
  }
  
  DelMatrix(&m1);
  DelMatrix(&m2);
  
  uivector *v;
  NewUintVector(6, &v);
  for(i = 0; i < 6; i++)
    v->data[i] = i;

  VectorUintAppend(&v, 10);
  for(i = 0; i < v->size; i++)
    printf("%d\n", v->data[i]);

  DelUintVector(&v);
}

void test10()
{
  uint i;
  uivector *v;
  NewUintVector(6, &v);
  for(i = 0; i < 6; i++)
    v->data[i] = i;

  VectorUintAppend(&v, 10);
  for(i = 0; i < v->size; i++)
    printf("%d\n", v->data[i]);

  DelUintVector(&v);
}

int main(void){
//   test1();
//  test2();
//   test3();
//   test4();
//   test5();
//   test6();
  test7();
  test8();
//   test9();
  return 0;
}