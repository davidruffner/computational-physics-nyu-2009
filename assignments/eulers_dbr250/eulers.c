//David Ruffner
//Computational Physics
//11-10-09
//1D Euler solver
//---------------------------------------------------------------------------
//Used some template functions from Jonathan Zrake's advection code

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h> 

#include "euler.h"
#include "config.h"

#define STRING_LENGTH 128
#define gamma 1.5

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

// Apply a function to a whole array of Conserved Structs
// -----------------------------------------------------------------------
void apply_over_arrayCons2(struct Conserved (*f)(double x),
			   struct Conserved *y, double *x, int N)
{
  int i;
  for (i=0; i<N; ++i)
    {
      y[i].rho = f(x[i]).rho;
      y[i].E = f(x[i]).E;
      y[i].px = f(x[i]).px;
    }
}


// Initial condition functions
// -----------------------------------------------------------------------
//Shock Tube with interface at x=0
struct Conserved shock(double x, struct Conserved Ur, struct Conserved Ul)
{
  if ( x < 0)
    {
      return Ul;
    }
  else
    {
      return Ur;
    }
}

//sets values of U array for a shock tube
//---------------------------------------------------------------------------
void shockInitial(struct Conserved *y, struct Conserved Ur, 
		  struct Conserved Ul, double *x, int N)
{
  int i;
  for (i=0; i<N; ++i)
    {
      y[i] = shock(x[i],Ur, Ul);
    }
}


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
void write_ascii_tableCons(const char *outname,
		       double *x, struct Conserved *vals, int N)
{
  FILE *outfile = fopen(outname, "w");
  int i;
  for (i=0; i<N; ++i)
    {
      fprintf(outfile, "%f %f %f %f\n", x[i], vals[i].rho, vals[i].E, 
	      vals[i].px);
    }
  fclose(outfile);
}

// -----------------------------------------------------------------------
void write_ascii_tablePrim(const char *outname,
		       double *x, struct Primitive *vals, int N)
{
  FILE *outfile = fopen(outname, "w");
  int i;
  for (i=0; i<N; ++i)
    {
      fprintf(outfile, "%f %f %f %f\n", x[i], vals[i].pre, vals[i].e, 
	      vals[i].vx);
    }
  fclose(outfile);
}

//Numerical Approx of dx*dU/dt
//------------------------------------------------------------
void L(struct Conserved (*flux)(struct Conserved, 
	   struct Conserved) ,struct Conserved *U, struct Conserved *LU, int N)
{
  struct Conserved *Fiph = (struct Conserved*) calloc(N+1, 
						      sizeof(struct Conserved));
  //N+1 Fiph cells because flux is needed on both sides of each U
  
  int i;
  for(i = 0; i< N+1; ++i)//Check about these edges
    {
      Fiph[i] = flux(U[i+1],U[i]);
      //printf("Fiph[i].px = %f \n", Fiph[i].px);
    }

  int j;
  for(j=1; j < N+1; ++j)
    {
      LU[j-1].rho = - Fiph[j].rho + Fiph[j-1].rho;
      LU[j-1].px = - Fiph[j].px + Fiph[j-1].px;
      LU[j-1].E = - Fiph[j].E + Fiph[j-1].E;
      //printf("LU.px = %f \n", LU[j-1].px);
    } 
      
  free(Fiph);
  
}

//Time step of Euler Equations first order, using F_HLL
//----------------------------------------------------------------------
double evolveF_HLL_1stO(struct Conserved *U, int N, double dx, double CFL)
{
  struct Conserved *LU = (struct Conserved*) calloc(N, 
			    sizeof(struct Conserved));
  
  double lamdaMax;
  double dt;

  lamdaMax = MaxEigenVal(U, N);

  dt = CFL*dx/lamdaMax;//Satisfies the Courant Condition

  L(F_HLL, U, LU, N);//This creates the dx*dU/dx for all Ui

  int j;
  for(j=1; j < N+1; ++j)
    {
      
      U[j].rho = U[j].rho + (dt/dx)*LU[j-1].rho;
      U[j].px = U[j].px + (dt/dx)*LU[j-1].px;
      U[j].E = U[j].E + (dt/dx)*LU[j-1].E;
      
    } 
  //Outflow BCs
  U[0].rho = U[1].rho;
  U[0].px = U[1].px;
  U[0].E = U[1].E;

  U[N+1].rho = U[N].rho;
  U[N+1].px = U[N].px;
  U[N+1].E = U[N].E;
  free(LU);
  return dt;
}

int main()
{

  struct Conserved Ur;
  struct Conserved Ul;
  //put in initial values of Ur Ul
  int    Nx    =  get_int_param("param.cfg", "Nx", 10);
  
  double Mrpre = get_double_param("param.cfg", "Mr.pre", 0);
  double Mrrho = get_double_param("param.cfg", "Mr.rho", 0);
  double Mrv   = get_double_param("param.cfg", "Mr.v", 0);
  
  double Mlpre = get_double_param("param.cfg", "Ml.pre", 0);
  double Mlrho = get_double_param("param.cfg", "Ml.rho", 0);
  double Mlv   = get_double_param("param.cfg", "Ml.v", 0);

  double runTime   = get_double_param("param.cfg", "runTime", 0);

  Ur = mixed_to_cons(Mrpre, Mrrho, Mrv );
  Ul = mixed_to_cons(Mlpre, Mlrho, Mlv );

  double x0    = -3.0;
  double x1    = 5.0; 
  double dx;
  
  double *x = (double*) calloc(Nx+2, sizeof(double));
  
  linspace(x, x0, x1, Nx+2, &dx); 
     
  
  struct Conserved *Uinitial = (struct Conserved*) 
    calloc(Nx+2, sizeof(struct Conserved));//Size N+2 for the two ghost cells
  struct Primitive *Pinitial = (struct Primitive*) 
    calloc(Nx+2, sizeof(struct Primitive));//Size N+2 for the two ghost cells
  
  shockInitial(Uinitial,Ur,Ul,x,Nx+2);
  
  write_ascii_tableCons("testInitial.cfg",x, Uinitial,Nx+2);
  
  double CFL = 0.5;
  double time = 0;
  
  /*
  time = time+ evolveF_HLL_1stO(Uinitial, Nx, dx, CFL);
  
  write_ascii_tableCons("testFirstStep.cfg",x, Uinitial ,Nx +2);

  time = time + evolveF_HLL_1stO(Uinitial, Nx, dx, CFL);
  
  write_ascii_tableCons("testSecondStep.cfg",x, Uinitial,Nx+2);
  
  time = time + evolveF_HLL_1stO(Uinitial, Nx, dx, CFL);
  
  write_ascii_tableCons("testThirdStep.cfg",x, Uinitial,Nx+2);
  */

  int j=0;
  while(time<= runTime)
    {
      j++;
      time = time + evolveF_HLL_1stO(Uinitial, Nx, dx, CFL);
    }
   printf("the time after the %d th  step is: %f \n" ,j, time);
   
   int i;
   
   for(i=1;i<Nx+2;++i)
     {
       Pinitial[i] = cons_to_prim(Uinitial[i]);
     }
   write_ascii_tableCons("testFinalCons.cfg",x, Uinitial,Nx+2);
   write_ascii_tablePrim("testFinalPrim.cfg",x, Pinitial,Nx+2);
  
  /*  Parameters for this run
  // ---------------------------------------------------------------------

  //int    Nx    =  get_int_param("param.cfg", "Nx", 10);
  //double x0    =  get_double_param("param.cfg", "x0",-1.0);
  //double x1    =  get_double_param("param.cfg", "x1", 1.0);
  double CFL   =  get_double_param("param.cfg", "CFL", 0.5);
  //double a     =  get_double_param("param.cfg", "a", 1.0);
  double t_max =  get_double_param("param.cfg", "t_max", 1.0);

  //which initial condition (1. is Gaussian, 2. is Square wave)
  //int init = get_int_param("param.cfg", "init", 1 );

  //which method to integrate (1. is Lax Friedrichs, 2. is Lax Wendroff)
  //int meth = get_int_param("param.cfg", "meth", 1);

  



  // Initial conditions
  // ---------------------------------------------------------------------
  
  //double dx, dt;
  //struct Conserved *U = (struct Conserved*) 
  //calloc(Nx, sizeof(struct Conserved));
    //double *x = (double*) calloc(Nx, sizeof(double));
  //struct Conserved *U_anal = (struct Conserved*)
  //calloc(Nx, sizeof(struct Conserved));
  //double *x_t = (double*) calloc(Nx, sizeof(double));

  //linspace(x, x0, x1, Nx, &dx);            // like numpy.linspace
  
   if(init == 1)
    {
      // apply_over_array(initial1, U, x, Nx);// like U = initial(x)
    }
  else if(init == 2)
    {
      // apply_over_array(initial2, U, x, Nx);// like U = initial(x)
    }
  else
    {
      //apply_over_array(initial1, U, x, Nx);// like U = initial(x)
    }
  
  //int i;
  for(i=0;i<Nx;++i)
    {
      x_t[i]=x[i];
    }
  
 
  
  //write_ascii_table("advect_ini.dat", x, U, Nx);
  
  // Time evolution
  // ---------------------------------------------------------------------
  double time    = 0.0;
 

  int Nt;
  //int tstep;
  int ncycle =0;

  //dt = CFL * dx / fabs(a);
  Nt = (t_max)*fabs(a)*Nx/(CFL*(x1-x0));
  printf("num time steps %d \n", Nt);
  printf("dt %f \n", dt);
  
  //for(tstep = 0; tstep<=Nt; ++tstep)
  while(time<= t_max) //
    {
      if(meth == 1)
	{
	  //evolve_lax_friedrichs(U, a, dt/dx, Nx);
	}
      else if(meth == 2)
	{
	  //evolve_lax_Wendroff(U, a, dt/dx, Nx);
	}
      else if(meth == 3)
	{
	  //evolve_ConsrvForm(godunov_NonLinflux, U, a, dt/dx, Nx);
	}
      else
	{
          //evolve_lax_friedrichs(U, a, dt/dx, Nx);
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

  
  
  //L1norm(U_anal, U, dx, Nx,"norm.dat");
  //L1normSimp(U_anal, U, "normSimp.dat");
  //write_ascii_table("advect_fin.dat", x, U, Nx);

  
  
  free(x);
  free(U);
  */
  //free(x_t);
  free(Pinitial);
  free(Uinitial);
  free(x);
  return 0;
}
