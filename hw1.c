#include <stdio.h> /* printing */
#include <float.h> /* check against system answers */
#include <math.h>  /* pow */

int main(int argc,char* argv[])
{
  /* variables for problem */
  int iter1=0,iter2=0,iter3=0,iter4=0;
  float single_epsilon = 1.0f;
  double double_epsilon = 1.0;
  int s_exp = -1;
  int d_exp = -1;

  /* class print header */
  printf("\nPrinciples of Numerical Computation\n");
  printf("David Parker\n");
  printf("Homework 1, September 2, 2011\n\n");

  /* Traditional way of solving the Machine Epsilon problem 
     (Note: this is for me) */
  /* Single Precision */
  do {
    //    printf("Iteration %d = %g\n", iter1, single_epsilon);
    /* reassign single_epsilon by dividing by 2 */
    single_epsilon /= 2.0f;
    /* iterate iterator */
    iter1 += 1; 
    /* 
       check if single_epsilon + 1.0 is greater than 1.0.
       if it is, then reloop, otherwise we found our machine epsilon.
       Note that this must be cast as a float otherwise C returns the
       incorrect answer 
    */
  } while ((float)(single_epsilon + 1.0) > 1.0);
  single_epsilon *= 2.0f; /* return to correct answer*/
  printf("\nCalculated single precision = %g\n",single_epsilon);
  printf("C-defined single precision  = %g\n\n",FLT_EPSILON);

  /* Double Precision */
  do {
    //    printf("Iteration %d = %g\n", iter2, double_epsilon);
    /* reassign double_epsilon by dividing by 2 */
    double_epsilon /= 2.0f;
    /* iterate iterator */
    iter2 += 1;
    /* 
       check if double_epsilon + 1.0 is greater than 1.0.
       if it is, then reloop, otherwise we found our machine epsilon 
    */
  } while ((double_epsilon + 1.0) > 1.0);
  double_epsilon *= 2.0f; /* return to correct answer */
  printf("\nCalculated double precision = %g\n",double_epsilon);
  printf("C defined double precision = %g\n\n",DBL_EPSILON);


  /* 
     PROBLEM 1 
  */
  printf("******************************\n");
  printf("*          PROBLEM 1         *\n");
  printf("******************************\n");
  /* 1) Single Precision */
  do {
    //    printf("Iteration %d = %g\n", iter3, pow(2.0,s_exp));
    /* reassign s_exp by dividing by 2 */
    s_exp -= 1;
    /* iterate iterator */
    iter3 += 1;
    /* 
       Check if (2^s_exp)+1 is greater than 1. If it is, then reloop, 
       as we haven't found our machine epsilon.
       Otherwise, we found our machine epsilon.
       Note that this must be cast as a float, otherwise C computes
       it as a double.
    */
  } while ((float)(pow(2.0,s_exp) + 1.0) > 1.0);
  s_exp += 1; /* return to correct answer as we went one past it */
  printf("\nCalculated lowest value of exponent = %d\n", s_exp);
  printf("Calculated single precision = %g\n",pow(2.0,s_exp));
  printf("C-defined single precision  = %g\n\n",FLT_EPSILON);

  /* 2) Double Precision */
  do {
    //    printf("Iteration %d = %g\n", iter4, pow(2.0,d_exp));
    /* reassign d_exp by dividing by 2 */
    d_exp -= 1;
    /* iterate iterator */
    iter4 += 1;
    /* 
       Check if (2^d_exp)+1 is greater than 1. If it is, then reloop, 
       as we haven't found our machine epsilon.
       Otherwise, we found our machine epsilon.
    */
  } while ((pow(2.0,d_exp) + 1.0) > 1.0);
  d_exp += 1; /* return to correct answer as we went one past it */
  printf("\nCalculated lowest value of exponent = %d\n", d_exp);
  printf("Calculated double precision = %g\n",pow(2.0,d_exp));
  printf("C-defined double precision  = %g\n\n",DBL_EPSILON);

  return 0;
}
