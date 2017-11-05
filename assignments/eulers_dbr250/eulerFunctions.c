
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h> 

#define STRING_LENGTH 128
#define gamma 1.5

#include "euler.h"

//Functions for dealing with conserved and primitive quantities
//---------------------------------------------------------------------
/*struct Conserved
{
  double rho, E, px;
};

struct Primitive
{
  double pre, e, vx;
}; */

//Converts from input for shock tube to Conserved
//---------------------------------------------------------------
struct Conserved mixed_to_cons(double pre, double rho, double vx)
{
  struct Conserved U;
  U.rho = rho;
  U.px  = vx*rho;
  U.E   = 0.5*rho*pow(vx,2) + pre/(gamma-1);
  return U;
}

//---------------------------------------------------------------
struct Conserved cons_to_mixed(double pre, double rho, double vx)
{
  struct Conserved U;
  U.rho = rho;
  U.px  = vx*rho;
  U.E   = 0.5*rho*pow(vx,2) + pre/(gamma-1);
  return U;
}

//--------------------------------------------------------
struct Conserved prim_to_cons(struct Primitive inprim)
{
  double P;
  double e;
  double vx;
  //double gamma;

  
  struct Conserved outcons;

  P = inprim.pre;
  e = inprim.e;
  vx = inprim.vx;
   
  
  outcons.E = P*(0.5*pow(vx,2)+e)/((gamma-1)*e);
  outcons.px = P*vx/((gamma-1)*e);
  outcons.rho = P/((gamma-1)*e);

  return outcons;
}

//---------------------------------------------------
struct Primitive cons_to_prim(struct Conserved incons)
{
  double rho;
  double px;
  double E;
  //double gamma;

  

  struct Primitive outprim;

  rho = incons.rho;
  px = incons.px;
  E = incons.E;
   
  
  outprim.pre = (gamma-1)*(E-0.5*pow(px,2)/rho);
  outprim.vx = px/rho;
  outprim.e = E/rho - 0.5*pow(px,2)/pow(rho,2);

  return outprim;
}

//Flux at for a given value of U
//--------------------------------------------
struct Conserved flux(struct Primitive prim)
{
  double P;
  double e;
  double vx;

  struct Conserved cons;
  struct Conserved outFlux;
  double rho;
  double px;
  double E;
  //double gamma;

  //gamma = 1.5;//for a 3/2 gas

  P = prim.pre;
  e = prim.e;
  vx = prim.vx;
  
  cons = prim_to_cons(prim);

  rho = cons.rho;
  E   = cons.E;
  px  = cons.px;

  //printf("%f %f %f \n" , P, E, vx);

  outFlux.rho = rho*vx;
  outFlux.px   = rho*pow(vx,2) + P;
  outFlux.E  = (E+P)*vx;

  return outFlux;
}

//Eigenvalues of the flux jacobian, effective wave speeds
//---------------------------------------------------------------
void Eigenvalues(struct Primitive inprim, double* lamdap, double* lamdam)
{
  double vx;
  double P;
  double rho;
  double cs;

  vx = inprim.vx;
  P  = inprim.pre;
  rho = prim_to_cons(inprim).rho;

  //double gamma;
  //gamma = 1.5;//for a 3/2 gas
 
  cs = pow(gamma*P/rho,0.5);

  *lamdap = vx + cs;
  *lamdam = vx - cs;
}

// Minimum of two numbers
//---------------------------------
double Min2(double num1, double num2)
{
  if(num1>num2)
    {
      return num2;
    }
  else
    {
      return num1;
    }
}
// Max of two numbers
//---------------------------------
double Max2(double num1, double num2)
{
  if(num1<num2)
    {
      return num2;
    }
  else
    {
      return num1;
    }
}

//Finds Alphas needed for HLL function
//-------------------------------------------
double Alphavalue(double lamdal, double lamdar)
{
  double alpha;
  alpha = 0;
  alpha = Max2(0 , Max2(lamdal,lamdar) );
  return alpha; 
}

double MaxEigenVal(struct Conserved *U, int N)
{
  struct Primitive *Uprim = (struct Primitive*) calloc(N, 
			    sizeof(struct Primitive));
  double lamda1 = 0;
  double lamda2 = 0;
  double lamdaMax = 0;
  double tempMax = 0;

  int i;

  for(i=0; i<N-1; ++i)
    {
      Uprim[i] = cons_to_prim(U[i]);
      Eigenvalues(Uprim[i], &lamda1, &lamda2);
      tempMax = Max2(lamda1, lamda2);//maybe abs value
      lamdaMax = Max2(tempMax, lamdaMax);
    }
  free(Uprim);
  return lamdaMax;
}
  

//The HLL flux
//------------------------------------------------
struct Conserved F_HLL(struct Conserved Ur, struct Conserved Ul)
{
  struct Conserved fluxHLL;

  struct Conserved fluxr;
  struct Conserved fluxl;
  struct Primitive Pr;
  struct Primitive Pl;

  double lamdarP = 0;
  double lamdarM = 0;
  double lamdalP = 0;
  double lamdalM = 0;

  double alphaP = 0;
  double alphaM = 0;
  
  Pr = cons_to_prim(Ur);
  Pl = cons_to_prim(Ul);

  fluxr = flux(Pr);
  fluxl = flux(Pl);

  Eigenvalues(Pr, &lamdarP, &lamdarM);
  Eigenvalues(Pl, &lamdalP, &lamdalM); 

  alphaP = Alphavalue(lamdalP, lamdarP);
  alphaM = Alphavalue(-lamdalM, -lamdarM);

  //printf("fluxr.px = %f \n" , fluxr.px);
  //printf("fluxl.px = %f \n" , fluxl.px);

  fluxHLL.rho = (alphaP*fluxl.rho + alphaM*fluxr.rho - 
		 alphaP*alphaM*(Ur.rho-Ul.rho) )/(alphaP+alphaM);

  fluxHLL.E = (alphaP*fluxl.E + alphaM*fluxr.E - 
		 alphaP*alphaM*(Ur.E-Ul.E) )/(alphaP+alphaM);

  fluxHLL.px = (alphaP*fluxl.px + alphaM*fluxr.px - 
		 alphaP*alphaM*(Ur.px-Ul.px) )/(alphaP+alphaM);

  return fluxHLL;
}

