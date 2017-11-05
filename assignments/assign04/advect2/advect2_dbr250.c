
/*
  ------------------------------------------------------------------------------
  Author: Jonathan Zrake, NYU CCPP
  Date: 10/9/09

  Demonstration of simple scalar advection in 1d, using the C language. As in
  the C++ version, I have tried to draw parallels to the Python way of doing
  things. Concepts to be familiar with are dynamic memory allocation, using
  calloc(), passing pointers to functions to modify the contents of an array,
  and passing functions as arguments to functions.

*/

#include <stdio.h>              // for printf, fprintf
#include <stdlib.h>             // for calloc, free
#include <math.h>               // for sqrt, exp
#include "config.h"             // for get_double_param, get_int_param

// Emulate numpy.linspace behavior
// -----------------------------------------------------------------------
void linspace(double *vals, double a, double b, int N, double *dx)
{
  *dx = (b-a) / (N);//Important, I found that if I used the factor N instead
  //of N-1, then I wouldn't get an over count as my analytical solution went
  // around the boundary
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


// Initial condition functions
// -----------------------------------------------------------------------
// Gaussian
double initial1(double x)
{
  return exp(-x*x / 0.1);
}
// Square
double initial2(double x)
{
  if(x<= 0.1 && x>= -0.1)
    {
      return 1;
    }
  else
    {
      return 0;
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

// Advection schemes
// -----------------------------------------------------------------------
//Lax Friedrichs
//--------------
void evolve_lax_friedrichs(double *U, double a, double lam, int N)
{
  double *Uip1 = (double*) calloc(N, sizeof(double));
  double *Uim1 = (double*) calloc(N, sizeof(double));
  int i;
  for (i=0; i<N; ++i)
    {
      Uip1[i] = U[(i+1+N) % N]; // this acts like the numpy.roll function
      Uim1[i] = U[(i-1+N) % N];
    }
  for (i=0; i<N; ++i)
    {
      U[i] = 0.5*(Uip1[i] + Uim1[i]) - 0.5*a*lam * (Uip1[i] - Uim1[i]);
    }
  free(Uip1);
  free(Uim1);
}

//Lax Wendroff
//--------------
void evolve_lax_Wendroff(double *U, double a, double lam, int N)
{
  double *Uip1 = (double*) calloc(N, sizeof(double));
  double *Uim1 = (double*) calloc(N, sizeof(double));
  int i;
  for (i=0; i<N; ++i)
    {
      Uip1[i] = U[(i+1+N) % N]; // this acts like the numpy.roll function
      Uim1[i] = U[(i-1+N) % N]; // the % N does the boundary conditions
    }
  for (i=0; i<N; ++i)
    {
      U[i] = U[i] - 0.5*a*lam * (Uip1[i] - Uim1[i]) +
	.5*a*a*lam*lam*(Uip1[i] -2*U[i] + Uim1[i]);
    }
  free(Uip1);
  free(Uim1);
}


//method for calculationg the L1_Norm (outputs L1_norm to a file L1_norm.dat)
//------------------------------------------------------------------------ 
void L1norm(double *U_anal, double *U, double step, int N,const char *outname)
{
  double sum;
  double diff;
  double norm;
  int i;

  sum = 0;
  norm = 0;
  diff = 0;
  for (i=0; i<N; ++i)
    {
      diff = U_anal[i] - U[i];
      sum = sum + fabs(diff);
    }
  norm = step*sum;
  
  FILE *outfile = fopen(outname, "w");
  fprintf(outfile, "%lf\n", norm);
  fprintf(outfile, "%d\n", N);
  fprintf(outfile, "%lf\n", sum);
  fprintf(outfile, "%lf\n", diff);
  fprintf(outfile, "%lf\n", step);
  fclose(outfile);
}





int main()
{
  // Parameters for this run
  // ---------------------------------------------------------------------

  int    Nx    =  get_int_param("param.cfg", "Nx", 10);
  double x0    =  get_double_param("param.cfg", "x0",-1.0);
  double x1    =  get_double_param("param.cfg", "x1", 1.0);
  double CFL   =  get_double_param("param.cfg", "CFL", 0.5);
  double a     =  get_double_param("param.cfg", "a", 1.0);
  double t_max =  get_double_param("param.cfg", "t_max", 1.0);

  //which initial condition (1. is Gaussian, 2. is Square wave)
  int init = get_int_param("param.cfg", "init", 1 );

  //which method to integrate (1. is Lax Friedrichs, 2. is Lax Wendroff)
  int meth = get_int_param("param.cfg", "meth", 1);

  



  // Initial conditions
  // ---------------------------------------------------------------------
  
  double dx, dt;
  double *U = (double*) calloc(Nx, sizeof(double));
  double *x = (double*) calloc(Nx, sizeof(double));
  double *U_anal = (double*) calloc(Nx, sizeof(double));
  double *x_t = (double*) calloc(Nx, sizeof(double));

  linspace(x, x0, x1, Nx, &dx);            // like numpy.linspace
  
  if(init == 1)
    {
      apply_over_array(initial1, U, x, Nx);// like U = initial(x)
    }
  else if(init == 2)
    {
      apply_over_array(initial2, U, x, Nx);// like U = initial(x)
    }
  else
    {
      apply_over_array(initial1, U, x, Nx);// like U = initial(x)
    }

  int i;
  for(i=0;i<Nx;++i)
    {
      x_t[i]=x[i];
    }
  
 
  
  write_ascii_table("advect_ini.dat", x, U, Nx);
  
  // Time evolution
  // ---------------------------------------------------------------------
  double time    = 0.0;
 

  int Nt;
  int tstep;
  int ncycle =0;

  dt = CFL * dx / fabs(a);
  Nt = (t_max)*fabs(a)*Nx/(CFL*(x1-x0));
  printf("num time steps %d \n", Nt);
  printf("dt %f \n", dt);
  
  //for(tstep = 0; tstep<=Nt; ++tstep)
  while(time<= t_max) //
    {
      if(meth == 1)
	{
	  evolve_lax_friedrichs(U, a, dt/dx, Nx);
	}
      else if(meth == 2)
	{
	  evolve_lax_Wendroff(U, a, dt/dx, Nx);
	}
      else
	{
	  evolve_lax_friedrichs(U, a, dt/dx, Nx);
	}
      time += dt;
      ncycle++;
    }
  printf("The number of cycles: %d \n", ncycle);
  printf("time %f \n", time);
  //Calculate the analytical array at this time
  //--------------------------------------------

  
  for(i=0;i<Nx;++i)
    {
      x_t[i]=x[i]-a*(time);
      while(x_t[i] < x0 || x_t[i] > x1)
	{
	  if(x_t[i] < x0)
	    {
	      x_t[i] = x_t[i]-(x0-x1);
	    }
	  else
	    {
	      x_t[i] = x_t[i]+(x0-x1);
	    }
	}
    }

  write_ascii_table("x_t.dat", x, x_t, Nx);

  if(init == 1)
    {
      apply_over_array(initial1, U_anal, x_t, Nx);// like U = initial(x)
    }
  else if(init == 2)
    {
      apply_over_array(initial2, U_anal, x_t, Nx);// like U = initial(x)
    }
  else
    {
      apply_over_array(initial1, U_anal, x_t, Nx);// like U = initial(x)
    }
  
  write_ascii_table("advect_ana.dat", x, U_anal, Nx);

  
  
  L1norm(U_anal, U, dx, Nx,"norm.dat");
  //L1normSimp(U_anal, U, "normSimp.dat");
  write_ascii_table("advect_fin.dat", x, U, Nx);

  
  
  free(x);
  free(U);

  free(x_t);
  free(U_anal);

  return 0;
}
