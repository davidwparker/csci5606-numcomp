#include <stdio.h> /* printing */
#include <float.h> /* machine epsilon */
#include <math.h>  /* pow */

/* Poor man's approximation of PI */
#define PI 3.1415926535898

/*
 * Prototypes
 */
/*
 * Algorithms, etc
 */
void newton(double (*f)(double),double (*fp)(double),double,int);

/*
 * Homework problems
 */
double function_p93_14a(double);
double function_p93_14aprime(double);
double function_p93_14b(double);
double function_p93_14bprime(double);
void problem93_14();

/*
 * newton
 * ------
 * implementation of the newton algorithm
 * takes the function, the fprime, starting x0, and max as input
 */
void newton(double (*f)(double), double (*fp)(double), double x0, int max) {
  double v = (*f)(x0);
  double delta = DBL_EPSILON;
  double eps = DBL_EPSILON;
  double x1;
  int i;

  if (fabs(v) < eps) {
    printf("below eps\n\n");
    return;
  }
  printf("i              x                  f(x)=v\n");
  printf("    %20.20g     %20.20g\n",x0,v);
  for (i = 0; i < max; i++) {
    double v1 = (*fp)(x0);
    if (v1 == 0) {
      printf("Cannot divide by zero\n");
      return;
    }
    x1 = x0 - v/v1;
    v = (*f)(x1);
    printf("%d  %.20g    %.20g\n",i,x1,v);
    /* sanity checks */
    if (fabs(v) < eps) {
      printf("below eps\n\n");
      return;
    }
    if (fabs(x1-x0) < delta) {
      printf("below delta\n\n");
      return;
    }
    x0 = x1;
  }
}

/*
 * Function wrappers
 */

/*
 * bisect_example_func
 * ------
 * e^x - sin(x)
 */
double bisect_example_func(double x) {
  return (exp(x) - sin(x));
}

/*
 * newton_example_func and newton_example_func_prime
 * ------
 * e^x - 1.5 - tan^-1(x)
 * e^x - (1+x^2)^-1
 */
double newton_example_func(double x) {
  return (exp(x) - 1.5 - atan(x));
}

double newton_example_func_prime(double x) {
  return (exp(x) - pow((1+pow(x,2)),-1));
}


/*
 * functions
 */
double function_p93_14a(double x) {
  //          4y2 + 4y + 52x = 19
  //169x2 + 3y2 + 111x - 10y = 10
  return 0.0;
}

double function_p93_14aprime(double x) {
  return 0.0;
}

double function_p93_14b(double x) {
  //          x + e-1x + y3 = 0
  //x2 + 2xy - y2 +  tan(x) = 0
  return 0.0;
}

double function_p93_14bprime(double x) {
  return 0.0;
}

// problems
void problem32_4() {
  printf("Problem 3.2 #4\n");
  printf("f(x)  = x^3 - 5x^2 + 3x - 7 with x0 = 5\n");
  printf("f'(x) = 3x^2 - 10x + 3\n");
  //  newton(function_p32_4,function_p32_4prime,5,10);
  return;
}

/*
 * Final - Page 93, computer problem 14
 */
void problem93_14() {
  printf("Problem 14\n");
  printf("a.\n");
  printf("          4y2 + 4y + 52x = 19\n");
  printf("169x2 + 3y2 + 111x - 10y = 10\n");
  //  newton(function_p93_14a,function_p93_14aprime,19,MAX);

  printf("b.\n");
  printf("          x + e-1x + y3 = 0\n");
  printf("x2 + 2xy - y2 +  tan(x) = 0\n");
  //  newton(function_p93_14b,function_p93_14bprime,19,MAX);
  return;
}

int main(int argc,char* argv[])
{
  /* class print header */
  printf("\nPrinciples of Numerical Computation\n");
  printf("David Parker\n");
  printf("Final Project, November 29, 2011\n\n");

  problem93_14();
  
  return 0;
}
