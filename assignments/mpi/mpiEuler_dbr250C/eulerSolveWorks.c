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
#include "mpi.h"
#include "euler.h"
#include "config.h"

#define STRING_LENGTH 128
#define gamma 1.5


// Emulate numpy.linspace behavior
// -----------------------------------------------------------------------
void linspace2nd(double *vals, double a, double b, int N, double *dx)
{
  *dx = (b-a) / (N);//Important, I found that if I used the factor N instead
  //of N-1, then I wouldn't get an over count as my analytical solution went
  // around the boundary
  int i;
  for (i=0; i<(N+4); ++i)//For the two ghost cell on each side
    {
      vals[i] = a + (i-2) * (*dx);
    }
}
// First Order Emulate numpy.linspace behavior
// -----------------------------------------------------------------------
void linspace(double *vals, double a, double b, int N, double *dx)
{
  *dx = (b-a) / (N);//Important, I found that if I used the factor N instead
  //of N-1, then I wouldn't get an over count as my analytical solution went
  // around the boundary
  int i;
  for (i=0; i<(N+2); ++i)//For the two ghost cell on each side
    {
      vals[i] = a + (i-1) * (*dx);
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
      //printf("%f \n", y[i].rho);
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
			   double *x, struct Conserved *vals, int N, char* perm)
{
  FILE *outfile = fopen(outname, perm);
  int i;
  for (i=2; i<N+2; ++i)//for 2nd euler
    {
      fprintf(outfile, "%f %f %f %f\n", x[i], vals[i].rho, vals[i].E, 
	      vals[i].px);
    }
  fclose(outfile);
}

// -----------------------------------------------------------------------
void write_ascii_tablePrim(const char *outname,
			   double *x, struct Primitive *vals, int N, char* perm)
{
  FILE *outfile = fopen(outname, perm);
  int i;
  for (i=2; i<N+2; ++i)//for 2nd euler
    {
      fprintf(outfile, "%f %f %f %f\n", x[i], vals[i].pre, vals[i].e, 
	      vals[i].vx);
    }
  fclose(outfile);
}

// -----------------------------------------------------------------------
void read_ascii_tableCons(const char *inname, double *x, 
			  struct Conserved *vals, int N)
//This N is the actual N it reads
{
  
  FILE *infile = fopen(inname, "r");
  int i;

  
  for (i=0; i<N; ++i)
    {
      fscanf(infile, "%lf %lf %lf %lf\n", &x[i], &vals[i].rho, &vals[i].E, 
	      &vals[i].px);
    }
  fclose(infile);
}

// -----------------------------------------------------------------------
void read_ascii_tablePrim(const char *inname, double *x, 
			  struct Primitive *vals, int N)
//This N is the actual N it reads
{
  
  FILE *infile = fopen(inname, "r");
  int i;

  
  for (i=0; i<N; ++i)
    {
      fscanf(infile, "%lf %lf %lf %lf\n", &x[i], &vals[i].pre, &vals[i].e, 
	      &vals[i].vx);
    }
  fclose(infile);
}

//Sign function
//--------------------
double sgn(double x)
{
  if(x==0)
    {
      return 1.0;
    }
  else
    {
      return x/fabs(x);
    }
}

//Min Mod Function
//-----------------------------------------------------------------
double minmod(double x, double y, double z)
{
  return .25*fabs( sgn(x)+sgn(y) )*(sgn(x)+sgn(y))*Min2(fabs(x),Min2(fabs(y),fabs(z) ) );
}

//PLM interpolation
//------------------------------------------------------------------

//Reconstruct (Assumes two ghost points at the end of the array)
//-------------------------------------------------------------
void reconstruct(struct Conserved *U, struct Conserved *dUdx, int N)
{
  double theta = 1.0;
  int i;
  for(i=0; i< N+2 ; ++i)
    {
      dUdx[i].rho = minmod(theta*(U[i+1].rho-U[i].rho),
			   .5*(U[i+2].rho-U[i].rho),
			   theta*(U[i+2].rho-U[i+1].rho) );
      dUdx[i].px = minmod(theta*(U[i+1].px-U[i].px),
			   .5*(U[i+2].px-U[i].px),
			  theta*(U[i+2].px-U[i+1].px) );
      dUdx[i].E = minmod(theta*(U[i+1].E-U[i].E),
			   .5*(U[i+2].E-U[i].E),
			 theta*(U[i+2].E-U[i+1].E) );
    }
}

//Flux at Fiph in secondOrder
//----------------------------------------------------------------------
void Flux_at_iph2nd(struct Conserved (*flux)(struct Conserved, 
	  struct Conserved), 
	  struct Conserved *U, struct Conserved *Fiph, double N)
{
  struct Conserved *dUdx = (struct Conserved*) calloc(N+2, 
						     sizeof(struct Conserved));
  reconstruct(U,dUdx, N);
  
  int i;
  struct Conserved Ul;
  struct Conserved Ur;
  for(i=0; i<N+1; ++i)
    {
      Ul.rho = U[i+1].rho + .5*dUdx[i].rho;
      Ul.px = U[i+1].px + .5*dUdx[i].px;
      Ul.E = U[i+1].E + .5*dUdx[i].E;

      Ur.rho = U[i+2].rho - .5*dUdx[i+1].rho;
      Ur.px = U[i+2].px - .5*dUdx[i+1].px;
      Ur.E = U[i+2].E - .5*dUdx[i+1].E;

      Fiph[i] = flux(Ur, Ul);
    }
}
      
//Flux at Fiph
//----------------------------------------------------------------
void Flux_at_iph(struct Conserved (*flux)(struct Conserved, 
	  struct Conserved) ,struct Conserved *U, struct Conserved *Fiph, int N)
{
  int i;
  for(i = 0; i< N+1; ++i)//Check about these edges
    {
      Fiph[i] = flux(U[i+1],U[i]);
      //printf("Fiph[i].px = %f \n", Fiph[i].px);
    }
}
//Numerical Approx of dx*dU/dt  2nd Order
//----------------------------------------------------------------
void L2nd(struct Conserved *U, struct Conserved *LU, int N)
{
  struct Conserved *Fiph = (struct Conserved*) calloc(N+1, 
						      sizeof(struct Conserved));
  //N+1 Fiph cells because flux is needed on both sides of each U
  
  Flux_at_iph2nd(F_HLL, U, Fiph, N);

  int j;
  for(j=1; j < N+1; ++j)
    {
      LU[j-1].rho = - Fiph[j].rho + Fiph[j-1].rho;
      LU[j-1].px = - Fiph[j].px + Fiph[j-1].px;
      LU[j-1].E = - Fiph[j].E + Fiph[j-1].E;
      //printf("LU.px = %f \n", LU[j-1].px);
    } 
  //printf("Fiph[0].px = %f, Fiph[1].px = %f \n", Fiph[0].px, Fiph[1].px);   
  free(Fiph);
  
}

//Numerical Approx of dx*dU/dt
//------------------------------------------------------------
void L(struct Conserved *U, struct Conserved *LU, int N)
{
  struct Conserved *Fiph = (struct Conserved*) calloc(N+1, 
						      sizeof(struct Conserved));
  //N+1 Fiph cells because flux is needed on both sides of each U
  
  Flux_at_iph(F_HLL, U, Fiph, N);

  int j;
  for(j=1; j < N+1; ++j)
    {
      LU[j-1].rho = - Fiph[j].rho + Fiph[j-1].rho;
      LU[j-1].px = - Fiph[j].px + Fiph[j-1].px;
      LU[j-1].E = - Fiph[j].E + Fiph[j-1].E;
      //printf("LU.px = %f \n", LU[j-1].px);
    } 
  //printf("Fiph[0].px = %f, Fiph[1].px = %f \n", Fiph[0].px, Fiph[1].px);   
  free(Fiph);
  
}

//Boundary Conditions
//----------------------------------------------------------------------
//First Order:assuming U is N+2 long

//Outflow
void BC_1stO_out(struct Conserved *U, int N)
{
  U[0].rho = U[1].rho;
  U[0].px  = U[1].px;
  U[0].E   = U[1].E;

  U[N+1].rho = U[N].rho;
  U[N+1].px  = U[N].px;
  U[N+1].E   = U[N].E;
}
//Reflecting
void BC_1stO_rflt(struct Conserved *U, int N)
{
  U[0].rho = U[1].rho;
  U[0].px  = -U[1].px;
  U[0].E   = U[1].E;

  U[N+1].rho = U[N].rho;
  U[N+1].px  = -U[N].px;
  U[N+1].E   = U[N].E;
}

//Reflecting left out right
void BC_1stO_Rout_Lrflt(struct Conserved *U, int N)
{
  U[0].rho = U[1].rho;
  U[0].px  = -U[1].px;
  U[0].E   = U[1].E;

  U[N+1].rho = U[N].rho;
  U[N+1].px  = U[N].px;
  U[N+1].E   = U[N].E;
}

//Second Order Outflow
//-----------------------------------------
void build_CONSERVED(struct Conserved incons, MPI_Datatype* message_type_ptr);

void BC_2ndO_out(struct Conserved *U, int N, int my_rank, int np)
{
  MPI_Status status;
  MPI_Datatype CONSERVED;
  struct Conserved incons;
  build_CONSERVED(incons, &CONSERVED);
  printf("befoer BC, rank: %d \n", my_rank);
  
  struct Conserved GhostsR[2];//from the right of the sender
  struct Conserved GhostsL[2];
  struct Conserved IR[2];//from the right of the sender
  struct Conserved IL[2];

  MPI_Barrier(MPI_COMM_WORLD);
  if(my_rank == 0)
    {
     
      U[0] = U[2];
      U[1] = U[2];

            
      IR[0] = U[N];
      IR[1] = U[N+1];

      MPI_Send(IR, 2, CONSERVED, my_rank+1, 50, MPI_COMM_WORLD);

      MPI_Recv(GhostsR, 2, CONSERVED, my_rank+1, 50, MPI_COMM_WORLD, &status);
      
      U[N+3] = GhostsR[1];
      U[N+2] = GhostsR[0];

      
    }
  else if(my_rank == np-1)
    {
      U[N+3] = U[N+1];
      U[N+2] = U[N+1];
          
        
      IL[0] = U[2];
      IL[1] = U[3];

      MPI_Send(IL, 2, CONSERVED, my_rank-1, 50, MPI_COMM_WORLD);
      MPI_Recv(GhostsL, 2, CONSERVED, my_rank-1, 50, MPI_COMM_WORLD, &status);
      
      U[0] = GhostsL[0];
      U[1] = GhostsL[1];
      
          
    }
  else
    {
       
      
      IL[0] = U[2];
      IL[1] = U[3];
      IR[0] = U[N];
      IR[1] = U[N+1];

      printf("U[0].rho %f \n", U[0].rho);
      printf("GhostsL[0].rho %f \n", GhostsL[0].rho);
      
      MPI_Send(IL, 2, CONSERVED, my_rank-1, 50, MPI_COMM_WORLD);
      MPI_Recv(GhostsL, 2, CONSERVED, my_rank-1, 50, MPI_COMM_WORLD, &status);      
      
      MPI_Send(IR, 2, CONSERVED, my_rank+1, 50, MPI_COMM_WORLD);       
      MPI_Recv(GhostsR, 2, CONSERVED, my_rank+1, 50, MPI_COMM_WORLD, &status);
      printf("After received GhostsL[0].rho %f \n", GhostsL[0].rho);
      
      U[0] = GhostsL[0];
      U[1] = GhostsL[1];

      U[N+2] = GhostsR[0];
      U[N+3] = GhostsR[1];

           
      
    }
  printf("after BC, rank: %d \n", my_rank);
  MPI_Barrier(MPI_COMM_WORLD);  
}

void Utemp_BC(struct Conserved *Utemp, struct Conserved *U, int N)
{
  Utemp[0].rho = U[2].rho;
  Utemp[0].px  = U[2].px;
  Utemp[0].E   = U[2].E;

  Utemp[1].rho = U[2].rho;
  Utemp[1].px  = U[2].px;
  Utemp[1].E   = U[2].E;
  
  Utemp[N+3].rho = U[N+1].rho;
  Utemp[N+3].px  = U[N+1].px;
  Utemp[N+3].E   = U[N+1].E;

  Utemp[N+2].rho = U[N+1].rho;
  Utemp[N+2].px  = U[N+1].px;
  Utemp[N+2].E   = U[N+1].E;
}

//Applys the appropriate BCs based on BCnum
//------------------------------------------
void BC(int BCnum, struct Conserved *U, int N, int my_rank, int np)
{
  BC_2ndO_out(U,N, my_rank, np);
   
}
  



//Time step of Euler Equations 2nd order, using F_HLL
//----------------------------------------------------------------------
double evolveF_HLL_2ndO(struct Conserved *U, int N, double dx, double CFL,
			double BCnum, int my_rank, int np)
{
  struct Conserved *LU = (struct Conserved*) calloc(N, 
			    sizeof(struct Conserved));
  struct Conserved *Utemp = (struct Conserved*) calloc(N+4, 
			    sizeof(struct Conserved));
  double lamdaMax;
  double lamdaMaxGlobal;
  double dt;

  lamdaMax = MaxEigenVal(U, N);
  

  MPI_Allreduce(&lamdaMax, &lamdaMaxGlobal, 1, MPI_DOUBLE,MPI_MAX, MPI_COMM_WORLD);
 

  dt = CFL*dx/lamdaMaxGlobal;//Satisfies the Courant Condition
  
  L2nd(U, LU, N);//This creates the dx*dU/dx for all Ui
 
  //third order Runge-Kutta Time Step
  int j;
  for(j=2; j < N+2; ++j)
    {
      Utemp[j].rho = U[j].rho + (dt/dx)*LU[j-2].rho;
      Utemp[j].px = U[j].px + (dt/dx)*LU[j-2].px;
      Utemp[j].E = U[j].E + (dt/dx)*LU[j-2].E;
    }
  BC(BCnum, Utemp, N, my_rank, np);  

  if(my_rank==0)
    {
      int count;
      for(count=0; count<N+4;count++)
	{
	  printf("Rank 0 Utemp[%d].rho: %f \n", count, Utemp[count].rho);
	}
    }
  if(my_rank==1)
    {
      int count;
      for(count=0; count<N+4;count++)
	{
	  printf("Rank 1 Utemp[%d].rho: %f \n", count, Utemp[count].rho);
	}
    }
  if(my_rank==2)
    {
      int count;
      for(count=0; count<N+4;count++)
	{
	  printf("Rank 2 Utemp[%d].rho: %f \n", count, Utemp[count].rho);
	}
    }
  
  

  L2nd(Utemp, LU, N);
  
  
  for(j=2; j < N+2; ++j)
    {
      Utemp[j].rho = .75*U[j].rho + 0.25*Utemp[j].rho+ 0.25*(dt/dx)*LU[j-2].rho;
      Utemp[j].px  = .75*U[j].px  + 0.25*Utemp[j].px + 0.25*(dt/dx)*LU[j-2].px;
      Utemp[j].E   = .75*U[j].E   + 0.25*Utemp[j].E  + 0.25*(dt/dx)*LU[j-2].E;
    }
  BC(BCnum, Utemp, N, my_rank, np);
  

  L2nd(Utemp, LU, N);
  
  
  for(j=2; j < N+2; ++j)
    {
      U[j].rho = (1/3.0)*U[j].rho + (2/3.0)*Utemp[j].rho+ (2/3.0)*(dt/dx)*LU[j-2].rho;
      U[j].px  = (1/3.0)*U[j].px  + (2/3.0)*Utemp[j].px + (2/3.0)*(dt/dx)*LU[j-2].px;
      U[j].E   = (1/3.0)*U[j].E   + (2/3.0)*Utemp[j].E  + (2/3.0)*(dt/dx)*LU[j-2].E;
    }
  BC(BCnum, U, N,my_rank, np);
  
  free(Utemp);
  free(LU);
  
  return dt;
}

// MPI: this makes a MPI data type for structs
//---------------------------------------------------------------------------
void build_CONSERVED(struct Conserved incons, MPI_Datatype* message_type_ptr)
  {
    int block_lengths[3];
    MPI_Aint displacements[3];
    MPI_Aint addresses[4];
    MPI_Datatype typelist[3];
  
  

    typelist[0] = MPI_DOUBLE;
    typelist[1] = MPI_DOUBLE;
    typelist[2] = MPI_DOUBLE;

  

    block_lengths[0] = block_lengths[1] = block_lengths[2] = 1;

  
  
    MPI_Address(&incons, &addresses[0]);
    MPI_Address(&(incons.rho), &addresses[1]);
    MPI_Address(&(incons.E), &addresses[2]);
    MPI_Address(&(incons.px), &addresses[3]);
    displacements[0] = addresses[1]-addresses[0];
    displacements[1] = addresses[2]-addresses[0];
    displacements[2] = addresses[3]-addresses[0];

    
  
    MPI_Type_struct(3, block_lengths, displacements, typelist, message_type_ptr);
    MPI_Type_commit(message_type_ptr);
  }



int main(int argc, char **argv)
{
  printf("start \n \n \n \n \n");
  //start up MPI
  //------------------------------------------
  int my_rank;
  int p;       //number of processes
  //int source;  //rank of sender
  //int dest;    //rank of receiver
  //int tag = 50;//tag of messages
  //char message[100]; //storage for the message
  //MPI_Status status; //Return status for receive

  MPI_Init(&argc, &argv); 
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  printf("initialized MPI really \n");

  
  //--------------------------------------------

  
  

  

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
  
  int    BCnum2 = get_int_param("param.cfg", "BCnum2", 4);
  int    InitialType = get_int_param("param.cfg", "InitialType", 1);

  double runTime   = get_double_param("param.cfg", "runTime", 0);

  Ur = mixed_to_cons(Mrpre, Mrrho, Mrv );
  Ul = mixed_to_cons(Mlpre, Mlrho, Mlv );

  double x0Global    = -3.0;
  double x1Global    = 5.0; 
  double x0, x1;
  int N;

  double dx2ndGlobal;
  double *x2ndGlobal = (double*) calloc(Nx+4, sizeof(double));
  double dx2nd;
  double *x2nd = (double*) calloc(Nx+4, sizeof(double));//shouldn't this by N not Nx?
  
  linspace2nd(x2ndGlobal, x0Global, x1Global, Nx, &dx2ndGlobal); 
  //Gets the dx
  
  N = Nx/p;//needs to be integer division
  printf("Nx = %d\n", Nx);
  printf("p = %d\n", p);
  printf("N = %d\n", N);
  if(my_rank < (p-1))
    {
      
      if(N==0)
	{
	  printf("Too small an array\n");
	}
      x0 = x0Global + dx2ndGlobal*my_rank*(N);
      x1 = x0Global + dx2ndGlobal*(my_rank*N+N);
      linspace2nd(x2nd, x0, x1, N, &dx2nd); 
      printf("dx2nd = %f \n", dx2nd);
      printf("x1 = %f in proc %d \n", x1, my_rank);
      printf("x0 = %f  in proc %d \n", x0, my_rank);
    }
  else 
    {
      x0 = x0Global + dx2ndGlobal*my_rank*(N);
      N = Nx - (p-1)*(N);
      x1 = x1Global;
      linspace2nd(x2nd, x0, x1, N, &dx2nd); 
      printf("x1 = %f in proc %d \n", x0, my_rank);
      printf("x0 = %f  in proc %d \n", x1, my_rank);
      printf("dx2nd = %f \n", dx2nd);
    }
  
  //Second Order
  struct Conserved *Uinitial2nd = (struct Conserved*) 
    calloc(N+4, sizeof(struct Conserved));//Size N+4 for the four ghost cells
  struct Primitive *Pinitial2nd = (struct Primitive*) 
    calloc(N+4, sizeof(struct Primitive));//Size N+4 for the four ghost cells
  
  
  
  if(InitialType == 1)
    {
      
      shockInitial(Uinitial2nd,Ur,Ul,x2nd,N+4);
      int count;
      for(count=0; count<N+4;count++)
	{
	  printf("the initial value of U[%d].rho: %f \n", count, Uinitial2nd[count].rho);
	}
     
      
      if(my_rank == 0)
	{
	printf("here\n");
	write_ascii_tableCons("run/initialIn2nd.cfg",x2nd, 
			    Uinitial2nd ,N, "w");
	} 
      printf("At barrier proc: %d \n", my_rank);
      MPI_Barrier(MPI_COMM_WORLD);

      int proc; 
      for(proc = 1; proc< p; proc++)
	{
	  printf("here proc: %d \n", my_rank);
	  if(proc == my_rank)
	    {
	      write_ascii_tableCons("run/initialIn2nd.cfg",x2nd, 
				Uinitial2nd , N , "a");
	    }
	  printf("here proc: %d \n", my_rank);
	  MPI_Barrier(MPI_COMM_WORLD);
	}
    }
  else if(InitialType == 2)
    {
      
      read_ascii_tableCons("run/testFinalCons2nd.cfg",x2nd, Uinitial2nd,Nx+4);
    }
  else if(InitialType == 3)
    {
      
      read_ascii_tablePrim("run/initialIn2.cfg",x2nd, Pinitial2nd,Nx+4);
      int i;
      
      for(i=0;i<Nx+4;++i)
	{
	  Uinitial2nd[i] = prim_to_cons(Pinitial2nd[i]);
	}
      
    }
  else if(InitialType == 4)
    {
      
      read_ascii_tableCons("run/initialIn2.cfg",x2nd, Uinitial2nd,Nx+4);
      
    }

  double CFL = 0.5;
  double time = 0;
  
    
  int j=0;
     
  time = 0;
  j=0;
  while(time<= runTime)
    {
      j++;
      time = time + evolveF_HLL_2ndO(Uinitial2nd, N, dx2nd, CFL, BCnum2, my_rank, p);
    }
  printf("For Second Order the time after the %d th  step is: %f \n" ,j, time);
  int i;
   
  
  for(i=0;i<N+4;++i)
     {
       Pinitial2nd[i] = cons_to_prim(Uinitial2nd[i]);
     }
 




  //output results
  //-------------------------------------------------------------
  if(my_rank == 0){
  printf("here\n");
  write_ascii_tableCons("run/testFinalCons2nd.cfg",x2nd, Uinitial2nd,N, "w");
  write_ascii_tablePrim("run/testFinalPrim2nd.cfg",x2nd, Pinitial2nd,N, "w");
  } /*end of if for rank == 0*/

  MPI_Barrier(MPI_COMM_WORLD);

  int procNum; 
  for(procNum = 1; procNum< p; procNum++)
    {
      if(procNum == my_rank)
	{
	  write_ascii_tableCons("run/testFinalCons2nd.cfg",x2nd, 
				Uinitial2nd,N, "a");
	  write_ascii_tablePrim("run/testFinalPrim2nd.cfg",x2nd, 
				Pinitial2nd,N, "a");
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }
  //--------------------------------------------------------   


  free(x2nd);
  free(x2ndGlobal);
  free(Pinitial2nd);
  free(Uinitial2nd);
  
  MPI_Finalize();
  return 0;
}
