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
int isPositive(double);
int isSignSame(double, double);
void bisect(double (*f)(double),double,double,int);
void bisectRecursive(double (*f)(double),double,double,int,int);
void newton(double (*f)(double),double (*fp)(double),double,int);

/*
 * Homework problems
 */
double bisect_example_func(double);
double newton_example_func(double);
double newton_example_func_prime(double);
double function_p31_1a(double);
double function_p31_1b(double);
double function_p31_1c(double);
double function_p31_1d(double);
double function_p31_4a(double);
double function_p31_4b(double);
double function_p32_1(double);
double function_p32_1prime(double);
double function_p32_4(double);
double function_p32_4prime(double);
void problem_bisect_example();
void problem_newton_example();
void problem31_1();
void problem31_4();
void problem32_1();
void problem32_4();

/*
 * isPositive
 * ------
 * determines whether a number is positive or not
 */
int isPositive(double s) {
  /* if number is greater or equal to zero, then yes, it is positive */
  if (s >= 0)
    return 1;
  return 0;
}

/*
 * isSignSame
 * ------
 * determines whether the sign between two numbers is the same or not
 */
int isSignSame(double s, double t) {
  /* if sign is same */
  if (isPositive(s) == isPositive(t))
    return 1;
  return 0;
}

/*
 * bisect
 * ------
 * implementation of the bisection algorithm
 * takes the function, two endpoints, and max as input
 */
void bisect(double (*f)(double), double a, double b, int max) {
  /* midpoint */
  double c;
  /* function values */
  double u = (*f)(a);
  double v = (*f)(b);
  double w;       /* function value of midpoint */

  double e = b - a; /* error */
  double delta = DBL_EPSILON; 
  double eps = DBL_EPSILON;
  int i;          /* for loop counter */

  if (isSignSame(u,v)) {
    printf("Sign doesn't change. No roots to find.\n\n");
    return;
  }
  printf("i          c                       w\n");
  for (i=0; i <= max; i++) {
    e = e/2.0;
    c = a + e;
    w = (*f)(c);
    printf("%d  %20.20g     %20.20g\n",i,c,w);

    if (w == 0) {
      printf("Perfect root found!\n\n");
      return;
    }

    /* sanity checks */
    if (fabs(e) < eps) {
      printf("below eps\n\n");
      return;
    }
    if (fabs(w) < delta) {
      printf("below delta\n\n");
      return;
    }

    /* reassignment */
    if (!isSignSame(w,u)) {
      b = c;
      v = w;
    }
    else {
      a = c;
      u = w;
    }
  }
  printf("\n");
}

void bisectRecursive(double (*f)(double), double a, double b, int max, int i) {
  /* midpoint */
  double c;
  /* function values */
  double u = (*f)(a);
  double v = (*f)(b);
  double w;       /* function value of midpoint */

  double e = b - a; /* error */
  double delta = DBL_EPSILON; 
  double eps = DBL_EPSILON;

  if (isSignSame(u,v)) {
    printf("Sign doesn't change. No roots to find.\n\n");
    return;
  }
  if (i <= max) {
    e = e/2.0;
    c = a + e;
    w = (*f)(c);
    printf("%d  %20.20g     %20.20g\n",i,c,w);

    if (w == 0) {
      printf("Perfect root found!\n\n");
      return;
    }

    /* sanity checks */
    if (fabs(e) < eps) {
      printf("below eps\n\n");
      return;
    }
    if (fabs(w) < delta) {
      printf("below delta\n\n");
      return;
    }

    /* reassignment */
    i++;
    if (!isSignSame(w,u)) {
      bisectRecursive((*f),a,c,max,i);
    }
    else {
      bisectRecursive((*f),c,b,max,i);
    }
  }
}

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
 * function_p31_1a
 * ------
 * x^-1 - tan(x) on [0, pi/2]
 */
double function_p31_1a(double x) {
  return (pow(x,-1.0) - tan(x));
}

/*
 * function_p31_1b
 * ------
 * x^-1 - 2^x on [0,1]
 */
double function_p31_1b(double x) {
  return (pow(x,-1.0) - pow(2,x));
}

/*
 * function_p31_1c
 * ------
 * 2^-x + e^x + 2cos(x) - 6 on [1,3]
 */
double function_p31_1c(double x) {
  return (pow(2,-x) + exp(x) + 2*cos(x) - 6);
}

/*
 * function_p31_1d
 * ------
 * (x^3 + 4x^2 + 3x + 5)/(2x^3 - 9x^2 + 18x - 2) on [0,4]
 */
double function_p31_1d(double x) {
  return ((pow(x,3) + 4*pow(x,2) + 3*x + 5)/(2*pow(x,3) - 9*pow(x,2) + 18*x - 2));
}

/*
 * function_p31_4a
 * ------
 * x^8 - 36x^7 + 546x^6 - 4536x^5 + 22449x^4 - 67284x^3 + 118124x^2 - 109584x + 40320 on [5.5,6.5]
 */
double function_p31_4a(double x) {
  return (pow(x,8) - 36*pow(x,7) + 546*pow(x,6) - 4536*pow(x,5) + 22449*pow(x,4) - 67284*pow(x,3) + 118124*pow(x,2) - 109584*x + 40320);
}

/*
 * function_p31_4b
 * ------
 * x^8 - 36.001x^7 + 546x^6 - 4536x^5 + 22449x^4 - 67284x^3 + 118124x^2 - 109584x + 40320 on [5.5,6.5]
 */
double function_p31_4b(double x) {
  return (pow(x,8) - 36.001*pow(x,7) + 546*pow(x,6) - 4536*pow(x,5) + 22449*pow(x,4) - 67284*pow(x,3) + 118124*pow(x,2) - 109584*x + 40320);
}

/*
 * function_p32_1 and function_p32_1prime
 * ------
 * x - tan(x)
 * 1 - sec^2(x)
 * --TODO:is this secant okay?
 */
double function_p32_1(double x) {
  return (x - tan(x));
}

double function_p32_1prime(double x) {
  return (1 - pow((1/cos(x)),2));
}

/*
 * function_p32_4 and function_p32_4
 * ------
 * x^3 - 5x^2 + 3x - 7
 * 3x^2 - 10x + 3
 */
double function_p32_4(double x) {
  return (pow(x,3) - 5*pow(x,2) + 3*x - 7);
}

double function_p32_4prime(double x) {
  return (3*pow(x,2) - 10*x + 3);
}

/*
 * Example Problems (to test my coded algorithms)
 */
void problem_bisect_example() {
  printf("Bisect example problem\n");
  printf("e^x - sin x on [-4,-3]\n");
  bisect(bisect_example_func,-4.0,-3.0,20);
  printf("Bisect example problem - recursive\n");
  printf("i          c                       w\n");
  bisectRecursive(bisect_example_func,-4.0,-3.0,20,0);
  printf("\n");
  return;
}

void problem_newton_example() {
  printf("Newton example problem\n");
  printf("f(x)  = e^x - 1.5 - tan^-1(x)\n");
  printf("f'(x) = e^x - (1 + x^2)^-1\n");
  newton(newton_example_func,newton_example_func_prime,-7.0,8);
  return;
}

/*
 * Problems...
 * ------
 * Computer problem 3.1 #1
 */
void problem31_1() {
  /* need to move it from undefined to a region that is computable */
  double x1 = (PI/2.0)-DBL_EPSILON*20;
  printf("Problem 3.1 #1\n");
  printf("a. x^-1 - tan x on [0,pi/2]\n");
  bisect(function_p31_1a,DBL_EPSILON,x1,10000);
  printf("b. x^-1 - 2^x on [0,1]\n");
  bisect(function_p31_1b,0.0,1.0,10000);
  printf("c. 2^-x + e^x + 2cos(x) - 6 on [1,3]\n");
  bisect(function_p31_1c,1.0,3.0,10000);
  printf("d. (x^3 + 4x^2 + 3x + 5)/(2x^3 - 9x^2 + 18x - 2) on [0,4]\n");
  bisect(function_p31_1d,0.0,4.0,10000);
  return;
}

/*
 * Computer problem 3.1 #4
 */
void problem31_4() {
  printf("Problem 3.1 #4\n");
  printf("a. x^8 - 36x^7 + 546x^6 - 4536x^5 + 22449x^4 - 67284x^3 + 118124x^2 - 109584x + 40320 on [5.5,6.5]\n");
  bisect(function_p31_4a,5.5,6.5,20);
  printf("b. x^8 - 36.001x^7 + 546x^6 - 4536x^5 + 22449x^4 - 67284x^3 + 118124x^2 - 109584x + 40320 on [5.5,6.5]\n");
  bisect(function_p31_4b,5.5,6.5,20);
  return;
}

/*
 * Computer problem 3.2 #1
 */
void problem32_1() {
  printf("Problem 3.2 #1\n");
  printf("a.\n");
  printf("f(x)  = x - tan(x) with x0 = 4.5\n");
  printf("f'(x) = 1 - sec^2(x)\n");
  newton(function_p32_1,function_p32_1prime,4.5,20);
  printf("b.\n");
  printf("f(x)  = x - tan(x) with x0 = 7.7\n");
  printf("f'(x) = 1 - sec^2(x)\n");
  newton(function_p32_1,function_p32_1prime,7.7,20);
  return;
}

/*
 * Computer problem 3.2 #4
 */
void problem32_4() {
  printf("Problem 3.2 #4\n");
  printf("f(x)  = x^3 - 5x^2 + 3x - 7 with x0 = 5\n");
  printf("f'(x) = 3x^2 - 10x + 3\n");
  newton(function_p32_4,function_p32_4prime,5,10);
  return;
}

int main(int argc,char* argv[])
{
  /* turn on/off function check */
  int checkFunctions = 1; 

  /* class print header */
  printf("\nPrinciples of Numerical Computation\n");
  printf("David Parker\n");
  printf("Homework 2, September 17, 2011\n\n");

  if (checkFunctions) {
  printf("C-defined single precision  = %g\n",FLT_EPSILON);
  printf("C-defined double precision  = %g\n\n",DBL_EPSILON);
  printf("Checking functions\n");
  printf("-0.01 positive? = %s\n",isPositive(-0.01)?"true":"false");
  printf("+0.01 positive? = %s\n",isPositive(+0.01)?"true":"false");
  printf("-0.01*-0.01 positive? = %s\n",isSignSame(-0.01,-0.01)?"true":"false");
  printf("+0.01*+0.01 positive? = %s\n",isSignSame(+0.01,+0.01)?"true":"false");
  printf("-0.01*+0.01 positive? = %s\n",isSignSame(-0.01,+0.01)?"true":"false");
  printf("+0.01*-0.01 positive? = %s\n\n",isSignSame(+0.01,-0.01)?"true":"false");
  problem_bisect_example();  
  problem_newton_example();
  printf("End checking functions\n\n");
  }
  
  problem31_1();
  problem31_4();
  problem32_1();
  problem32_4();
  
  return 0;
}
