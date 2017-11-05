#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h> 

#include "euler.h"

#define STRING_LENGTH 128
#define gamma 1.5

int main()
{
  double preL;
  double rhoL;
  double vL;
  double preR;
  double rhoR;
  double vR;
 
  preL = 1.0;
  rhoL = 10.0;
  vL   = 0.0;

  preR = 0.125;
  rhoR = 0.10;
  vR   = 0.0;

  struct Conserved Ur;
  struct Conserved Ul;

  Ul = mixed_to_cons(preL, rhoL, vL);
  Ur = mixed_to_cons(preR, rhoR, vR);

  printf("Ur: rho = %f, E= %f, px = %f \n"
	 , Ur.rho, Ur.E, Ur.px); 
  printf("Ul: rho = %f, E= %f, px = %f \n"
	 , Ul.rho, Ul.E, Ul.px); 
  /*
  //Test fluxHLL
  //------------------------------------
  struct Conserved Ur;
  struct Conserved Ul;
  struct Conserved fluxHLLvalue;

  struct Conserved fluxR;
  struct Conserved fluxL;

  struct Conserved U;

  U = mixed_to_cons(0.125,0.1,0.0);
  printf("U given Pre = 0.125, rho = .1, vL=0.0 : rho = %f, E= %f, px = %f \n"
	 , U.rho, U.E, U.px); 


  Ur.rho = 1.0;
  Ur.E = 6.5;
  Ur.px = -3.0;

  Ul.rho = 1.0;
  Ul.E = 6.5;
  Ul.px = 3.0;

  fluxR = flux( cons_to_prim(Ur) );
  fluxL = flux( cons_to_prim(Ul) );

  fluxHLLvalue = F_HLL(Ur,Ul);
  

  printf("Ur is rho = %f, E = %f, px = %f \n" , Ur.rho, Ur.E, Ur.px); 
  printf("Ul is rho = %f, E = %f, px = %f \n" , Ul.rho, Ul.E, Ul.px);
  printf("The flux(Ur) of rho = %f, of E = %f, of px =%f \n" 
	 , fluxR.rho, fluxR.E, fluxR.px);
  printf("The flux(UL) of rho = %f, of E = %f, of px =%f \n" 
	 , fluxL.rho, fluxL.E, fluxL.px);
  printf("The flux HLL of rho = %f,of  E = %f,of  px = %f \n"
	 , fluxHLLvalue.rho, fluxHLLvalue.E, fluxHLLvalue.px); 
  
  */


  //prim_to_cons Test  and flux function test
  //-------------------------
  /*struct Primitive primtest;
  struct Conserved constest;
  
  primtest.pre = 1.0;
  primtest.e   = 2.0;
  primtest.vx  = 3.0;

  double lamdap;
  double lamdam;

  lamdap = 0;
  lamdam = 0;

  Eigenvalues(primtest, &lamdap, &lamdam);
  printf("here are the eigenvalues: %f %f \n" , lamdap, lamdam);

  printf("here is the alpha(not right): %f \n" , Alphavalue(lamdap, lamdam) );

  double vals[3];

  vals[0] = primtest.pre;
  vals[1] = primtest.e;
  vals[2] = primtest.vx;

  printf("here are the primitive values: %f, %f, %f \n",vals[0],vals[1],vals[2] );
  
  struct Conserved fluxtest;
  fluxtest = flux(primtest);

  vals[0] = fluxtest.E;
  vals[1] = fluxtest.rho;
  vals[2] = fluxtest.px;

  printf("here are the corresponding flux values: \n %f, %f, %f \n",
                         vals[0],vals[1],vals[2] );
  
  constest = prim_to_cons(primtest);

  vals[0] = constest.E;
  vals[1] = constest.rho;
  vals[2] = constest.px;

 
  printf("here are the corresponding conserved values: \n %f, %f, %f \n",
                         vals[0],vals[1],vals[2] );

  primtest = cons_to_prim(constest);
  vals[0] = primtest.pre;
  vals[1] = primtest.e;
  vals[2] = primtest.vx;

  printf("here are the primitive values again: %f, %f, %f \n",vals[0],vals[1],vals[2] );
  
  */
  return 0;
}
