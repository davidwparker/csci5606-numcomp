#include <stdio.h> /* printing */
#include <float.h> /* machine epsilon */
#include <math.h>  /* sqrt */

#define NUM4 4
#define NUM10 10
#define PRINT 0

/*
 * Structs (for ease of use)
 */
typedef struct matrix10 {
  double v[NUM10][NUM10];
} matrix10;

typedef struct lu10 {
  double l[NUM10][NUM10];
  double u[NUM10][NUM10];
} lu10;

typedef struct matrix4 {
  double v[NUM4][NUM4];
} matrix4;

typedef struct matrix {
  int size;
  matrix4  m4;
  matrix10 m10;
} matrix;

typedef struct lu {
  lu10 lu10;
} lu;

typedef struct array {
  double v[NUM4];
} array;

/*
 * Prototypes
 */
/*
 * Algorithms
 */
struct matrix cholesky(struct matrix);
struct matrix transposeMatrix(struct matrix);
struct array forwardSubstitution(struct matrix, struct array);
struct array backSubstitution(struct matrix, struct array);

struct lu luDecomposition(struct matrix);
struct matrix10 evaluateLu(struct lu10 lu, struct matrix10 im);
struct matrix multiplyMatrix(struct matrix matrix1, struct matrix matrix2);

/*
 * Printing
 */
void printArray(struct array,char* arrayName);
void printMatrix(struct matrix);
void printLu10(struct lu10);

/*
 * Other methods
 */
struct matrix10 loadMatrix42_1(void);
struct matrix loadIdentityMatrix(int size);

/*
 * Homework problems
 */
void problem42_1(void);
void problem42_2(void);

/*
 * Methods
 */
/*
 * Algorithms
 */
/*
 * Cholesky
 * NOTE: this is setup for a NUMxNUM matrix only, and not setup for the general case
 *       in this case, NUM = 4
 */
struct matrix cholesky(struct matrix eq)
{
  struct matrix a;
  int i,k,s;
  a.size = eq.size;
  // being 0-based rather than 1-based
  for (k = 0; k < NUM4; k++) {
    // lkk = (akk - sum...)
    a.m4.v[k][k] = eq.m4.v[k][k];
    for (s = 0; s < k; s++) {
      if (PRINT) printf("a[k][s] = a[%d][%d] = %f  \n",k,s, a.m4.v[k][s]);
      a.m4.v[k][k] -= a.m4.v[k][s]*a.m4.v[k][s];
      if (PRINT) printf("a[k][k] = a[%d][%d] = %f  \n",k,k, a.m4.v[k][k]);
    }
    // sanity check
    if (a.m4.v[k][k] <= DBL_EPSILON) {
      printf("a[k][k] cannot be less than or equal to machine epsilon\n");
      printf("a[%d][%d] = %f  \n",k,k, a.m4.v[k][k]);
      return a;
    }
    a.m4.v[k][k] = sqrt(a.m4.v[k][k]);
    if (PRINT) printf("a[k][k] = a[%d][%d] = %f  \n\n",k,k, a.m4.v[k][k]);

    // next loop
    for (i = k+1; i < NUM4; i++) {
      a.m4.v[i][k] = eq.m4.v[i][k];
      for (s = 0; s < k; s++) {
	if (PRINT) {
	  printf("a[i][s] = a[%d][%d] = %f  \n",i,s, a.m4.v[i][s]);
	  printf("a[k][s] = a[%d][%d] = %f  \n",k,s, a.m4.v[k][s]);
	  printf("a[i][s]*a[k][s] = %f  \n", a.m4.v[i][s]*a.m4.v[k][s]);
	  printf("a[i][k] = a[%d][%d] = %f  \n",i,k, a.m4.v[i][k]);
	}
	a.m4.v[i][k] -= a.m4.v[i][s]*a.m4.v[k][s];
	if (PRINT) printf("a[i][k] = a[%d][%d] = %f  \n",i,k, a.m4.v[i][k]);
      }
      
      a.m4.v[i][k] /= a.m4.v[k][k];
    }
  }

  return a;
}

struct matrix transposeMatrix(struct matrix a)
{
  struct matrix at;
  int i,k;
  at.size = a.size;
  for (i = 0; i < a.size; i++) {
    for (k = 0; k < a.size; k++) {
      // check here
      at.m4.v[i][k] = a.m4.v[k][i];
      if (PRINT) printf("%f  ",at.m4.v[i][k]);
    }
    if (PRINT) printf("\n");
  }
  return at;
}

struct array forwardSubstitution(struct matrix a, struct array b)
{
  struct array y;
  int k,i;
  for (k = 0; k < a.size; k++) {
    y.v[k] = b.v[k];
    if (PRINT) printf("y[k] = y[%d] = %f  \n",k,y.v[k]);
    for (i = 0; i < k; i++) {
      // check here
      y.v[k] -= a.m4.v[k][i]*y.v[i];
      if (PRINT) printf("y[k]    = y[%d]    = %f  \n",k,y.v[k]);
    }
    // check here
    y.v[k] /= a.m4.v[k][k];
  }
  return y;
}

struct array backSubstitution(struct matrix at, struct array y)
{
  struct array x;
  int k,i;
  for (k = at.size-1; k >= 0; k--) {
    x.v[k] = y.v[k];
    if (PRINT) printf("x[k] = x[%d] = %f  \n",k,x.v[k]);
    for (i = k+1; i < at.size; i++) {
      // check here
      x.v[k] -= at.m4.v[k][i]*x.v[i];
      if (PRINT) printf("x[k] = x[%d] = %f  \n",k,x.v[k]);
    }
    // check here
    x.v[k] /= at.m4.v[k][k];
  }
  return x;
}

/*
 * nxn matrix; using doolittle algorithm
 */
struct lu luDecomposition(struct matrix a)
{
  struct lu lu;
  int k,j,i,s;

  for (k = 0; k < a.size; k++) {
    // check here
    lu.lu10.l[k][k] = 1;
    for (j = k; j < a.size; j++) {
      double sum = 0;
      for (s = 0; s < k-2; s++) {
	//check here
	sum += lu.lu10.l[k][s]*lu.lu10.u[s][j];
      }
      // check here
      lu.lu10.u[k][j] = a.m10.v[k][j] - sum;
    }
    for (i = k+1; i < a.size; i++) {
      double sum2 = 0;
      for (s = 0; s < k-2; s++) {
	// check here
	sum2 += lu.lu10.l[i][s]*lu.lu10.u[s][k];
      }
      // check here
      lu.lu10.l[i][k] = a.m10.v[i][k] - sum2;
      lu.lu10.l[i][k] /= lu.lu10.u[k][k];
    }
  }
  return lu;
}

struct matrix10 evaluateLu(struct lu10 lu, struct matrix10 im)
{
  struct matrix10 ans;
  int i,j,k;
  double b[NUM10];
  double y[NUM10];
  double x[NUM10];

  for (k = 0; k < NUM10; k++) {
    // load b
    for (i = 0; i < NUM10; i++) {
      b[i] = im.v[k][i];
    }

    // Forward solve Ly = b
    for (i = 0; i < NUM10; i++) {
      y[i] = b[i];
      for (j = 0; j < i; j++) {
	y[i] -= lu.l[i][j] * y[j];
      }
      y[i] /= lu.l[i][i];
    }
    // Backward solve Ux = y
    for (i = NUM10 - 1; i >= 0; i--) {
      x[i] = y[i];
      for (j = i + 1; j < NUM10; j++) {
	x[i] -= lu.u[i][j] * x[j];
      }
      x[i] /= lu.u[i][i];
      ans.v[i][k] = x[i];
    }
  }
  return ans;
}

struct matrix multiplyMatrix(struct matrix a, struct matrix b)
{
  struct matrix c;
  int k,j,s;
  c.size = a.size;

  if (a.size == NUM10) {
    // row of a
    for (k = 0; k < a.size; k++) {
      // column of b
      for (j = 0; j < a.size; j++) {
	// set to 0
	c.m10.v[k][j]=0;
	for (s = 0; s < a.size; s++) {
	  c.m10.v[k][j] += a.m10.v[k][s]*b.m10.v[s][j];
	}
      }
    }    
  }
  return c;
}

/*
 * Other methods
 */
struct matrix10 loadMatrix42_1()
{
  struct matrix10 a;
  int k,i;
  for (k = 0; k < NUM10; k++) {
    for (i = 0; i <= k; i++) {
      // +1 due to change from 0-based to 1-based
      a.v[k][i] = pow(((i+1)+(k+1)),2);
    }
  }
  return a;
}

struct matrix loadIdentityMatrix(int size)
{
  struct matrix im;
  im.size = size;
  int k,i;
  for (k = 0; k < size; k++) {
    for (i = 0; i < size; i++) {
      if (i == k) {
	if (size == NUM4)
	  im.m4.v[k][i] = 1;
	if (size == NUM10)
	  im.m10.v[k][i] = 1;
      } else {
	if (size == NUM4)
	  im.m4.v[k][i] = 0;
	if (size == NUM10)
	  im.m10.v[k][i] = 0;
      }
    }
  }
  return im;
}

/*
 * Printing convenience
 */
void printMatrix(struct matrix a)
{
  int k,i;
  for (k = 0; k < a.size; k++) {
    for (i = 0; i < a.size; i++) {
      // size check
      if (a.size == 4) 
	printf("%f  ",a.m4.v[k][i]);
      else if (a.size == 10)
	printf("%f  ",a.m10.v[k][i]);
    }
    printf("\n");
  }
}

void printArray(struct array a, char* s)
{
  int k;
  for (k = 0; k < NUM4; k++) {
    printf("%s[%d] = %f  \n",s,k,a.v[k]);
  }
  printf("\n");
}

void printLu10(struct lu10 lu)
{
  int k,i;
  printf("l:\n");
  for (k = 0; k < NUM10; k++) {
    for (i = 0; i < NUM10; i++) {
      printf("%.6f     ",lu.l[k][i]);
    }
    printf("\n");
  }

  printf("\nu:\n");
  for (k = 0; k < NUM10; k++) {
    for (i = 0; i < NUM10; i++) {
      printf("%.6f     ",lu.u[k][i]);
    }
    printf("\n");
  }
}


/*
 * Problems...
 * ------
 * Computer problem 4.2 #1
 */
void problem42_1()
{
  // variables
  struct matrix a,ai,aai,im;
  struct lu lu;

  // setup
  a.size = NUM10;
  ai.size = NUM10;
  aai.size = NUM10;
  im.size = NUM10;

  // print
  printf("Problem 4.2 #1\n");
  printf("Algorithm for inverting an nxn lower triangular matrix A.\n");
  printf("Testing on matrix whose elements are aij = (i+j)^2 when i >= j.\n");
  printf("Use n = 10. Form AA^-1 as test.\n");

  // do the work
  a.m10 = loadMatrix42_1();
  im = loadIdentityMatrix(NUM10);
  lu = luDecomposition(a);
  ai.m10 = evaluateLu(lu.lu10,im.m10);
  aai = multiplyMatrix(a,ai);

  // Original matrix
  printf("\nOriginal Matrix: A\n");
  printMatrix(a);
  // Identity matrix
  printf("\nIdentity Matrix: I\n");
  printMatrix(im);
  // LU decomposition
  printf("\nLU decomposition:\n");
  printLu10(lu.lu10);
  // Inverted matrix
  printf("\nInverted Matrix: A^-1\n");
  printMatrix(ai);
  // Product of A and A^-1
  printf("\nProduct of AA^-1 Matrix\n");
  printMatrix(aai);
  printf("\n");
  return;
}

/*
 * Computer problem 4.2 #2
 */
void problem42_2()
{
  // variables
  struct matrix mat,chol,cholT;
  struct matrix4 eq = {{
      {0.05,0.07,0.06,0.05},
      {0.07,0.10,0.08,0.07},
      {0.06,0.08,0.10,0.09},
      {0.05,0.07,0.09,0.10}
    }};
  struct array ans = {{0.23,0.32,0.33,0.31}};
  struct array y,x;

  // setup
  mat.m4 = eq;
  mat.size = NUM4;

  // print
  printf("Problem 4.2 #2\n");
  printf("Solve this system by Cholesky method:\n");
  printf("0.05x1 + 0.07x2 + 0.06x3 + 0.05x4 = 0.23\n");
  printf("0.07x1 + 0.10x2 + 0.08x3 + 0.07x4 = 0.32\n");
  printf("0.06x1 + 0.08x2 + 0.10x3 + 0.09x4 = 0.33\n");
  printf("0.05x1 + 0.07x2 + 0.09x3 + 0.10x4 = 0.31\n");

  // do the work
  chol = cholesky(mat);
  // transpose
  cholT = transposeMatrix(chol);
  // solve Ly = b (for y)
  y = forwardSubstitution(chol,ans);
  // solve L^Tx = y (for x)
  x = backSubstitution(cholT,y);

  // print out results
  printf("\nCholesky\n");
  printMatrix(chol);
  printf("\nTranspose \n");
  printMatrix(cholT);
  printf("\nLy = b (for y)\n");
  printArray(y,"y");
  printf("L^Tx = y (for x)\n");
  printArray(x,"x");
  printf("\n");
  return;
}

int main(int argc,char* argv[])
{
  /* class print header */
  printf("\nPrinciples of Numerical Computation\n");
  printf("David Parker\n");
  printf("Homework 3, October 6, 2011\n\n");

  /* For whatever reason, I'm getting some global state issues.
     Need to only run one problem at a time */
  problem42_1();
  //problem42_2();

  return 0;
}
