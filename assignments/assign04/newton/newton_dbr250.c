/*
  ------------------------------------------------------------------------------
  Author: David Ruffner
  Date: 10/21/09

  Takes values of P and gives values of T, given the equation P = a*T + b*T^4. 
  To do this it uses Newton-Rapheson method to solve for each T.

*/

#include <stdio.h>              // for printf, fprintf
#include <stdlib.h>             // for calloc, free
#include <math.h>               // for sqrt, exp
#include "config.h"             // for get_double_param, get_int_param


// Emulate numpy.linspace behavior
// -----------------------------------------------------------------------
void linspace(double *vals, double a, double b, int N, double *dx)
{
  *dx = (b-a) / (N-1.0);
  int i;
  for (i=0; i<N; ++i)
    {
      vals[i] = a + i * (*dx);
    }
}

// Apply a function to a whole array
// -----------------------------------------------------------------------
void apply_over_array(double (*f)(double), double *y, double *x, int N)
{
  int i;
  for (i=0; i<N; ++i)
    {
      y[i] = f(x[i]);
    }
}

// Output for plotting
// -----------------------------------------------------------------------
void write_ascii_table(const char *outname,
		       double *x, double *vals, int N)
{
  FILE *outfile = fopen(outname, "w");
  int i;
  for (i=0; i<N; ++i)
    {
      fprintf(outfile, "%f %f\n", x[i], vals[i]);
    }
  fclose(outfile);
}



double func(double T, double P, double a, double b)
{
  return P - a*T - b*pow(T,4);
}

double funcPrime(double T, double a, double b)
{
  return -a - 4*b*pow(T,3);
}
    
double newtonStep(double T, double P, double a, double b)
{
  return T - func(T,P,a,b)/funcPrime(T,a,b);
}

// Initial P values
// -----------------------------------------------------------------------


int main()
{
  int numSteps;
  
  double dP;
  double Tcur;
  double Tnew;
  double diff;
  int    NP    =  get_int_param("param.cfg", "NP", 10);
  double P0    =  get_double_param("param.cfg", "P0", 0.1);
  double P1    =  get_double_param("param.cfg", "P1", 1.0);
  double a    =  get_double_param("param.cfg", "a", 0.1);
  double b    =  get_double_param("param.cfg", "b", 0.1);




  double T0    =  get_double_param("param.cfg", "T0", 1.0);

  double tol   =  get_double_param("param.cfg", "tol", .01);

  // Initial P values
  // ---------------------------------------------------------------------
  double *P = (double*) calloc(NP, sizeof(double));
  double *T = (double*) calloc(NP, sizeof(double));
  
  linspace(P, P0, P1, NP, &dP);            // like numpy.linspace
  
  //Run Newton Rapheson for each pressure
  //-----------------------------------------------------------------------
  int i;
  for(i=0; i<NP; ++i)
    {
      numSteps = 0;
      diff = tol + 1;//to make it so the while loop starts
      Tcur = T0;
      while( fabs(diff)>= tol)
	{
	  Tnew = newtonStep(Tcur,P[i],a,b);
	  diff = Tnew - Tcur;
	  Tcur = Tnew;
	  numSteps++;
	}
      T[i] = Tcur;
      printf("number of iterations: %d \n", numSteps);
    }
  write_ascii_table("TvsP.dat", P, T, NP);

  free(P);
  free(T);

  return 0;
}
