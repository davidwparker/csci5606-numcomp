#include <stdio.h> /* printing */
#include <float.h> /* machine epsilon */
#include <math.h>  /* sqrt */

#define NUM2 2
#define NUM3 3
#define PRINT 0

/*
 * Structs (for ease of use)
 */
typedef struct matrix2 {
  double v[NUM2][NUM2];
} matrix2;

typedef struct matrix3 {
  double v[NUM3][NUM3];
} matrix3;

typedef struct matrix {
  int size;
  matrix2 m2;
  matrix3 m3;
} matrix;

typedef struct vector2 {
  double v[NUM2];
} vector2;

typedef struct vector3 {
  double v[NUM3];
} vector3;

typedef struct vector {
  int size;
  vector2 v2;
  vector3 v3;
} vector;

/*
 * Prototypes
 */
/*
 * Algorithms
 */

/*
 * Printing
 */
void printVector(struct vector,char* arrayName);
void printMatrix(struct matrix);

/*
 * Homework problems
 */
void problem46_1(void);
void problem46_1a(void);
void problem46_1b(void);
void problem46_2(void);

/*
 * Methods
 */
/*
 * Algorithms
 */
/*
 * Printing convenience
 */
void printMatrix(struct matrix a)
{
  int k,i;
  for (k = 0; k < a.size; k++) {
    for (i = 0; i < a.size; i++) {
      // size check
      if (a.size == NUM2) 
	printf("%f  ",a.m2.v[k][i]);
      else if (a.size == NUM3)
	printf("%f  ",a.m3.v[k][i]);
    }
    printf("\n");
  }
}

// size
void printVector(struct vector a, char* s)
{
  int k;
  for (k = 0; k < NUM2; k++) {
    printf("%s[%d] = %f  \n",s,k,a.v2.v[k]);
  }
  printf("\n");
}

/*
 * Problems...
 * ------
 * Computer problem 4.6 #1
 */
void problem46_1()
{
  // print
  printf("Problem 4.6 #1\n");
  printf("Program the Gauss-Seidel method and test it on these examples.\n");

  // do the work
  printf("a.\n");
  printf("3x + 1y + 1z =  5\n");
  printf("1x + 3y - 1z =  3\n");
  printf("3x + 1y - 5z = -1\n");
  problem46_1a();
  printf("\n");
  printf("b.\n");
  printf("3x + 1y + 1z =  5\n");
  printf("3x + 1y - 5z = -1\n");
  printf("1x + 3y - 1z =  3\n");
  problem46_1b();
  printf("\n");
  return;
}

void problem46_1a()
{
  // variables
  struct matrix a;
  struct matrix3 eq = {{
      {3,1,1},{1,3,-1},{3,1,-5},
    }};
  struct vector b,x;
  struct vector3 ans = {
    {5,3,-1}
  };
  struct vector3 start = {
    {0,0,0}
  };
  int max = 50;
  int k,i,j,z;

  // setup
  a.size = NUM3;
  a.m3 = eq;
  b.size = NUM3;
  b.v3 = ans;
  x.v3 = start;

  // TODO: refactor this into own method
  // solve gauss-seidel iteration
  for (k = 0; k < max; k++) {
    for (i = 0; i < a.size; i++) {
      double sum = 0;
      for (j = 0; j < a.size; j++) {
	if (i != j) {
	  sum += a.m3.v[i][j] * x.v3.v[j];
	  if (PRINT) printf("k = %d, j = %d, i = %d",k,j,i);
	  if (PRINT) printf(", a = %f, x = %f, sum = %f\n",a.m3.v[i][j],x.v3.v[j],sum);
	}
      }
      x.v3.v[i] = (b.v3.v[i] - sum)/a.m3.v[i][i];
      if (PRINT) printf("%f\n",x.v3.v[i]);
    }
  }
  
  for (z = 0; z < a.size; z++) {
    printf("%f\n",x.v3.v[z]);
  }
}

void problem46_1b()
{
  // variables
  struct matrix a;
  struct matrix3 eq = {{
      {3,1,1},{3,1,-5},{1,3,-1},
    }};
  struct vector b,x;
  struct vector3 ans = {
    {5,3,-1}
  };
  struct vector3 start = {
    {0,0,0}
  };
  int max = 50;
  int k,i,j,z;

  // setup
  a.size = NUM3;
  a.m3 = eq;
  b.size = NUM3;
  b.v3 = ans;
  x.v3 = start;

  // solve gauss-seidel iteration
  for (k = 0; k < max; k++) {
    for (i = 0; i < a.size; i++) {
      double sum = 0;
      for (j = 0; j < a.size; j++) {
	if (i != j) {
	  sum += a.m3.v[i][j] * x.v3.v[j];
	  if (PRINT) printf("k = %d, j = %d, i = %d",k,j,i);
	  if (PRINT) printf(", a = %f, x = %f, sum = %f\n",a.m3.v[i][j],x.v3.v[j],sum);
	}
      }
      x.v3.v[i] = (b.v3.v[i] - sum)/a.m3.v[i][i];
      if (PRINT) printf("%f\n",x.v3.v[i]);
    }
  }
  
  for (z = 0; z < a.size; z++) {
    printf("%f\n",x.v3.v[z]);
  }
}

/*
 * Computer problem 4.6 #2
 */
void problem46_2()
{
  // variables
  struct matrix a;
  struct matrix2 eq = {{
      {0.96326,0.81321},
      {0.81321,0.68654},
    }};
  struct vector b,x;
  struct vector2 ans = {
    {0.88824,0.74988},
  };
  struct vector2 start = {
    {0.33116,0.70000},
  };
  int max = 50000;
  int k,i,j,z;

  // setup
  a.size = NUM2;
  a.m2 = eq;
  b.size = NUM2;
  b.v2 = ans;
  x.v2 = start;

  // print
  printf("Problem 4.6 #2\n");
  printf("Apply Gauss-Seidel iteration to the system:\n");
  printf("0.96326x1 + 0.81321x2 = 0.88824\n");
  printf("0.81321x1 + 0.68654x2 = 0.74988\n");
  printf("Using (0.33116, 0.70000)^T as the starting point\n");

  // solve gauss-seidel iteration
  for (k = 0; k < max; k++) {
    for (i = 0; i < a.size; i++) {
      double sum = 0;
      for (j = 0; j < a.size; j++) {
	if (i != j) {
	  sum += a.m2.v[i][j] * x.v2.v[j];
	  if (PRINT) printf("k = %d, j = %d, i = %d",k,j,i);
	  if (PRINT) printf(", a = %f, x = %.10f, sum = %.10f\n",a.m2.v[i][j],x.v2.v[j],sum);
	}
      }
      x.v2.v[i] = (b.v2.v[i] - sum)/a.m2.v[i][i];
      if (PRINT) printf("%.10f\n",x.v2.v[i]);
    }
    if (PRINT) printf("%.10f\n",x.v2.v[0]);
  }

  // print out results
    for (z = 0; z < a.size; z++) {
    printf("%.10f\n",x.v2.v[z]);
  }
  printf("\n");
  return;
}

int main(int argc,char* argv[])
{
  /* class print header */
  printf("\nPrinciples of Numerical Computation\n");
  printf("David Parker\n");
  printf("Homework 4, October 8, 2011\n\n");

  problem46_1();
  problem46_2();

  return 0;
}
