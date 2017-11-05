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
  printf("in write file\n");
  FILE *outfile = fopen(outname, "w");
  int i;
  for (i=0; i<N; ++i)
    {
      printf("%f %f %f %f\n", x[i], vals[i].rho, vals[i].E, 
	      vals[i].px);
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
void BC_2ndO_out(struct Conserved *U, int N)
{
  U[0].rho = U[2].rho;
  U[0].px  = U[2].px;
  U[0].E   = U[2].E;

  U[1].rho = U[2].rho;
  U[1].px  = U[2].px;
  U[1].E   = U[2].E;
  
  U[N+3].rho = U[N+1].rho;
  U[N+3].px  = U[N+1].px;
  U[N+3].E   = U[N+1].E;

  U[N+2].rho = U[N+1].rho;
  U[N+2].px  = U[N+1].px;
  U[N+2].E   = U[N+1].E;
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
void BC(int BCnum, struct Conserved *U, int N)
{
  if(BCnum == 1)
    {
      BC_1stO_out(U,N);
    }
  else if(BCnum == 2)
    {
      BC_1stO_rflt(U,N);
    }
  else if(BCnum ==3)
    {
      BC_1stO_Rout_Lrflt(U,N);
    }
  else if(BCnum == 4)
    {
      BC_2ndO_out(U,N);
    }
}
  
//Time step of Euler Equations first order, using F_HLL
//----------------------------------------------------------------------
double evolveF_HLL_1stO(struct Conserved *U, int N, double dx, double CFL,
			double BCnum)
{
  struct Conserved *LU = (struct Conserved*) calloc(N, 
			    sizeof(struct Conserved));
  
  double lamdaMax;
  double dt;

  lamdaMax = MaxEigenVal(U, N);

  dt = CFL*dx/lamdaMax;//Satisfies the Courant Condition

  L(U, LU, N);//This creates the dx*dU/dx for all Ui
  
  
  
  int j;
  for(j=1; j < N+1; ++j)
    {
      U[j].rho = U[j].rho + (dt/dx)*LU[j-1].rho;
      U[j].px = U[j].px + (dt/dx)*LU[j-1].px;
      U[j].E = U[j].E + (dt/dx)*LU[j-1].E;
    }
   
  
  BC(BCnum, U, N);
  
  
  free(LU);
  return dt;
}


//Time step of Euler Equations 2nd order, using F_HLL
//----------------------------------------------------------------------
double evolveF_HLL_2ndO(struct Conserved *U, int N, double dx, double CFL,
			double BCnum)
{
  struct Conserved *LU = (struct Conserved*) calloc(N, 
			    sizeof(struct Conserved));
  struct Conserved *Utemp = (struct Conserved*) calloc(N+4, 
			    sizeof(struct Conserved));
  double lamdaMax;
  double dt;

  lamdaMax = MaxEigenVal(U, N);

  dt = CFL*dx/lamdaMax;//Satisfies the Courant Condition

  L2nd(U, LU, N);//This creates the dx*dU/dx for all Ui
 
  //third order Runge-Kutta Time Step
  int j;
  for(j=2; j < N+2; ++j)
    {
      Utemp[j].rho = U[j].rho + (dt/dx)*LU[j-2].rho;
      Utemp[j].px = U[j].px + (dt/dx)*LU[j-2].px;
      Utemp[j].E = U[j].E + (dt/dx)*LU[j-2].E;
    }
  BC(BCnum, Utemp, N);    
  

  L2nd(Utemp, LU, N);
  
  
  for(j=2; j < N+2; ++j)
    {
      Utemp[j].rho = .75*U[j].rho + 0.25*Utemp[j].rho+ 0.25*(dt/dx)*LU[j-2].rho;
      Utemp[j].px  = .75*U[j].px  + 0.25*Utemp[j].px + 0.25*(dt/dx)*LU[j-2].px;
      Utemp[j].E   = .75*U[j].E   + 0.25*Utemp[j].E  + 0.25*(dt/dx)*LU[j-2].E;
    }
  BC(BCnum, Utemp, N);
  

  L2nd(Utemp, LU, N);
  
  
  for(j=2; j < N+2; ++j)
    {
      U[j].rho = (1/3.0)*U[j].rho + (2/3.0)*Utemp[j].rho+ (2/3.0)*(dt/dx)*LU[j-2].rho;
      U[j].px  = (1/3.0)*U[j].px  + (2/3.0)*Utemp[j].px + (2/3.0)*(dt/dx)*LU[j-2].px;
      U[j].E   = (1/3.0)*U[j].E   + (2/3.0)*Utemp[j].E  + (2/3.0)*(dt/dx)*LU[j-2].E;
    }
  BC(BCnum, U, N);
  
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
  printf("initialized MPI \n");
  MPI_Datatype CONSERVED;
  struct Conserved incons;
  build_CONSERVED(incons, &CONSERVED);
  //--------------------------------------------



  //Parameters used by all the processors
  //-----------------------------------------------------

  double runTime   = get_double_param("param.cfg", "runTime", 0);
  int    BCnum2 = get_int_param("param.cfg", "BCnum2", 4);
  int    Nx    =  get_int_param("param.cfg", "Nx", 10);
  double x0    = -3.0;
  double x1    = 5.0; 
  double dx2ndGlobal;     
  double *x2ndGlobal = (double*) calloc(Nx+4, sizeof(double));
  linspace(x2ndGlobal, x0, x1, Nx+4, &dx2ndGlobal);   

  int i;
  double dx2nd;     
  double *x2nd = (double*) calloc(Nx+4, sizeof(double));

  struct Conserved *Uinitial2ndGlobal = (struct Conserved*) 
	calloc(Nx+4, sizeof(struct Conserved));
      //Size N+4 for the four ghost cells
  struct Primitive *Pinitial2ndGlobal = (struct Primitive*) 
	calloc(Nx+4, sizeof(struct Primitive));
      //Size N+4 for the four ghost cells
  struct Conserved *Uinitial2ndmsg = (struct Conserved*) 
	calloc(Nx+4, sizeof(struct Conserved));
      //Size N+4 for the four ghost cells
  struct Primitive *Pinitial2ndmsg = (struct Primitive*) 
	calloc(Nx+4, sizeof(struct Primitive));
  //Second Order
  struct Conserved *Uinitial2nd = (struct Conserved*) 
	calloc(Nx+4, sizeof(struct Conserved));
      //Size N+4 for the four ghost cells
  struct Primitive *Pinitial2nd = (struct Primitive*) 
	calloc(Nx+4, sizeof(struct Primitive));
      //Size N+4 for the four ghost cells

  //Setting up variable and Getting initial conditions
  //processor one will get the global initial conditions and send
  // it to all the other processors the each will determine which part 
  // they have
  //-------------------------------------------------------------------
  
  if(my_rank == 0)
    {
      struct Conserved Ur;
      struct Conserved Ul;
      
      
  
      double Mrpre = get_double_param("param.cfg", "Mr.pre", 0);
      double Mrrho = get_double_param("param.cfg", "Mr.rho", 0);
      double Mrv   = get_double_param("param.cfg", "Mr.v", 0);
  
      double Mlpre = get_double_param("param.cfg", "Ml.pre", 0);
      double Mlrho = get_double_param("param.cfg", "Ml.rho", 0);
      double Mlv   = get_double_param("param.cfg", "Ml.v", 0);
      
      
      int    InitialType = get_int_param("param.cfg", "InitialType", 1);
          
      
      Ur = mixed_to_cons(Mrpre, Mrrho, Mrv );
      Ul = mixed_to_cons(Mlpre, Mlrho, Mlv );

           
      
  
      
      linspace(x2nd, x0, x1, Nx+4, &dx2nd);   
  
      
      
  
      
  

      if(InitialType == 1)
	{
	  
	  shockInitial(Uinitial2nd,Ur,Ul,x2nd,Nx+4);
	  
	  
	  write_ascii_tableCons("run/testInitial2nd.cfg",
				x2nd, Uinitial2nd ,Nx+4);
	}
      else if(InitialType == 2)
	{
	  
	  read_ascii_tableCons("run/testFinalCons2nd.cfg",
			       x2nd, Uinitial2nd,Nx+4);
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
      
      for(i=0;i<Nx+4;++i)
	{
	  //printf("in proc %d Uinitial2nd[%d].rho = %f \n", my_rank, i, Uinitial2nd[i].rho);
	  Uinitial2ndGlobal[i] = Uinitial2nd[i];
	  //printf("in proc %d Uinitial2ndGlobal[%d].rho = %f \n"
	  //, my_rank, i, Uinitial2ndGlobal[i].rho);
	  
          Pinitial2ndGlobal[i] = cons_to_prim( Uinitial2nd[i] );
	  //printf("in proc %d Pinitial2ndGlobal[%d].e = %f \n"
	  //	 , my_rank, i, Pinitial2ndGlobal[i].e);
	  x2ndGlobal[i] = x2nd[i];
	}
      
      
      

    }/*end of the if statement to have proc0 initialize things*/


  MPI_Barrier(MPI_COMM_WORLD);
  printf("Before Uinitial2ndGlobal[2].rho = %f in proc %d \n", 
	 Uinitial2ndGlobal[2].rho, my_rank);
  printf("Before Uinitial2ndGlobal[2].rho = %f in proc %d \n", 
	 Uinitial2ndGlobal[2].rho, my_rank);

  MPI_Bcast(Uinitial2ndGlobal, Nx+4, CONSERVED, 0, MPI_COMM_WORLD);
  MPI_Bcast(Pinitial2ndGlobal, Nx+4, CONSERVED, 0, MPI_COMM_WORLD);
  MPI_Bcast(x2ndGlobal, Nx+4, CONSERVED, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  printf("here in proc %d \n", my_rank);
  int k;
  for(k = 0; k<Nx+4; k++)
    {
      //printf("Uinitial2ndGlobal[%d].E = %f in proc %d \n", k,
      //     Uinitial2ndGlobal[k].E, my_rank);
    }
  printf("Uinitial2ndGlobal[2].rho = %f in proc %d \n", Uinitial2ndGlobal[2].rho, my_rank);
  printf("x2ndGlobal[2] = %f in proc %d \n", x2ndGlobal[2], my_rank);
  double CFL = 0.5;
  double time = 0;
  printf("here1 in proc %d \n", my_rank);
  
  int j=0;

  for(j = 0; j<Nx+4; j++)
    {
      Uinitial2nd[i] = Uinitial2ndGlobal[i];
	  //printf("in proc %d Uinitial2ndGlobal[%d].rho = %f \n"
	  //, my_rank, i, Uinitial2ndGlobal[i].rho);
	  
      Pinitial2nd[i] = cons_to_prim( Uinitial2ndGlobal[i] );
	  //printf("in proc %d Pinitial2ndGlobal[%d].e = %f \n"
	  //	 , my_rank, i, Pinitial2ndGlobal[i].e);
      x2nd[i] = x2ndGlobal[i];
      
    }



  if(my_rank == 0)
    {
      
      time = 0;
      j=0;
      while(time<= runTime)
	{
	  j++;
	  time = time + evolveF_HLL_2ndO(Uinitial2nd, Nx, dx2nd, CFL, BCnum2);
	}
      printf("Proc %d: For Second Order the time after the %d th  step is: %f \n"
	     , my_rank, j, time);
      int i;
   
  //for(i=0;i<Nx+2;++i)
  // {
  //   Pinitial[i] = cons_to_prim(Uinitial[i]);
  // }
      for(i=0;i<Nx+4;++i)
	{
	  printf("After iteration in proc %d Uinitial2nd[%d].rho = %f \n"
	  	 , my_rank, i, Uinitial2nd[i].rho);
	  Pinitial2nd[i] = cons_to_prim(Uinitial2nd[i]);
	  printf("After iteration in proc %d Pinitial2nd[%d].pre = %f \n"
	  	 , my_rank, i, Pinitial2nd[i].pre);
	  printf("After iteration in proc %d x2nd[%d] = %f \n"
	  	 , my_rank, i, x2nd[i]);
	}
    }
  printf("here2 in proc %d \n", my_rank);

  
  if(my_rank == 1)
   {
      
      printf("about to write to file \n");
      write_ascii_tableCons("run/testFinalCons2nd.cfg",x2nd, Uinitial2nd,Nx+4);
      //write_ascii_tablePrim("run/testFinalPrim2nd.cfg",x2nd, Pinitial2nd,Nx+4);
  
      
   }
  
  printf("here3 in proc %d \n", my_rank);

  //free(Pinitial);
  //free(Uinitial);
  free(Pinitial2ndGlobal);
  free(Uinitial2ndGlobal);
  free(Pinitial2ndmsg);
  free(Uinitial2ndmsg);
  free(Uinitial2nd);
  free(Pinitial2nd);
  free(x2nd);
  //free(x);
  free(x2ndGlobal);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
