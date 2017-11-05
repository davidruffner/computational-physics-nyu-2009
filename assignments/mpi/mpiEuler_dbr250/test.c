#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "euler.h"

//#include "/usr/lib/openmpi/1.2.7-gcc/include/mpi.h"
// this will only work on desktop after I install mpicc

// mpd& gets the message passing deamon running

// use mpiexec -np 2 to get the thing running on 2 processors.  
//int darts(int N)
//{
int darts(int N, int rank)
{
  
  double x;
  double y;
  double r;
  int hits;

  
  srand(N*rank + time(NULL));

  hits = 0;
  int count;
  for(count = 0; count< N; count++)
    {
      x = ( (double)rand()/((double)(RAND_MAX)+(double)(1)) );
      y = ( (double)rand()/((double)(RAND_MAX)+(double)(1)) );

      r = pow((x - .5),2) + pow((y-.5),2);
      
      if( r <= .25)
	{
	  hits++;
	}
    }
  return hits;
}

struct my_type
{
  float a; 
  float b;
  int n;
};


/*void Build_derived_type(INDATA_TYPE* indata, MPI_Dataype* message_type_ptr)
{
  int block_lengths[3];
  MPI_Aint displacements[3];
  MPI_Aint addresses[4];
  MPI_Datatype typelist[3];
  
  

  typelist[0] = MPI_FLOAT;
  typelist[1] = MPI_FlOAT;
  typelist[2] = MPI_INT;

  

  block_lengths[0] = block_lengths[1] = block_lengths[2] = 1;

  
  
  MPI_Address(indata, &addresses[0]);
  MPI_Address(&(indata.a), &addresses[1]);
  MPI_Address(&(indata.b), &addresses[2]);
  MPI_Address(&(indata.n), &addresses[3]);
  displacements[0] = addresses[1]-addresses[0];
  displacements[1] = addresses[2]-addresses[0];
  displacements[2] = addresses[3]-addresses[0];

  
  
  MPI_Type_struct(3, block_lengths, displacements, typelist, message_type_ptr);

  

  MPI_Type_commit(message_type_ptr);
  } */ /*Build_derived_type*/

void Get_data2(int my_rank, float* a_ptr, float* b_ptr, int* n_ptr)
{
  int root = 0; /*arguements to MPI_Bcast*/
  //int count = 1;

  if(my_rank == 0)
    {
      printf("Enter a, b, and n \n");
      scanf("%f %f %d", a_ptr, b_ptr, n_ptr);
    }
  MPI_Bcast(a_ptr, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
  MPI_Bcast(b_ptr, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
  MPI_Bcast(n_ptr, 1, MPI_FLOAT, root, MPI_COMM_WORLD);

  return;
}

/*
void Get_data3(INDATA_TYPE* indata, int my_rank)
{
  MPI_Datatype message_type; 
  int root = 0;
  int count = 1;

  if(my_rank == 0)
  {
    printf("Enter a, b, and n \n");
    scanf("%f %f %d", &(indata->a), &(indata->b), &(indata->n));
  }
  
  Build_derived_type(indata, &message_type);
  MPI_Bcast(indata, count, message_type, root, MPI_COMM_WORLD);
} 
*/

int main(int argc, char **argv)
{
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
  
  float a;
  float b;
  int n;
  Get_data2(my_rank, &a, &b, &n);
   
  printf("Now the values of a, b, n = %f %f %d in processor %d \n" 
	 , a , b, n, my_rank);

  //Making a data structure type in MPI
  //----------------------------------------------------------------
  struct my_type indata;
  indata.a = a;
  indata.b = b;
  indata.n = n;

  
  int block_lengths[3];
  MPI_Aint displacements[3];
  MPI_Aint addresses[4];
  MPI_Datatype  MY_TYPE;
  MPI_Datatype typelist[3];
  
  

  typelist[0] = MPI_DOUBLE;
  typelist[1] = MPI_DOUBLE;
  typelist[2] = MPI_INT;

  

  block_lengths[0] = block_lengths[1] = block_lengths[2] = 1;

  
  
  MPI_Address(&indata, &addresses[0]);
  MPI_Address(&(indata.a), &addresses[1]);
  MPI_Address(&(indata.b), &addresses[2]);
  MPI_Address(&(indata.n), &addresses[3]);
  displacements[0] = addresses[1]-addresses[0];
  displacements[1] = addresses[2]-addresses[0];
  displacements[2] = addresses[3]-addresses[0];

  printf("here");
  
  MPI_Type_struct(3, block_lengths, displacements, typelist, &MY_TYPE);
  MPI_Type_commit(&MY_TYPE);
  //-----------------------------------------------------------------

  MPI_Bcast(&indata, 1, MY_TYPE, 0, MPI_COMM_WORLD);

  printf("This is indata on processor %d:  %f %f %d \n", my_rank, indata.a, indata.b,
	 indata.n);
  //struct indata values;
  //struct indata outvalues;
  
  //values.a = 2.0;
  //values.b = 3.0;
  //values.n = 1;

  //G/et_data3(&values, myrank);

  //MPI_Receive(

  //MPI_Recv(&outvalues, 1, message_type, 0, tag, 
  //	   MPI_COMM_WORLD, &status);





  /*
  
  int Ndarts;
  Ndarts = 1000;
  int hits;
  int localhits;
  int messagehits;

  if(my_rank != 0)
    {
      
  
      localhits = 0;
      localhits = darts(Ndarts, my_rank);
      printf("Number of hits %d in processor %d \n", localhits, my_rank);
      printf("Estimate of Pi: %f \n" , 4.0*(double)(localhits)/(double)(Ndarts));

      dest = 0;

      messagehits = localhits;
      MPI_Send(&messagehits, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
    }
  else // so in processor 0
    {
          
  
      
  
      localhits = 0;
      localhits = darts(Ndarts, my_rank);
      printf("Number of hits %d in processor %d \n", localhits, my_rank);
      printf("Estimate of Pi: %f \n" , 4.0*(double)(localhits)/(double)(Ndarts));

      hits = localhits;

      for(source = 1; source<p; source++)
	{
	  MPI_Recv(&messagehits, 1, MPI_INT, source, tag, 
		   MPI_COMM_WORLD, &status);
	  hits = hits + messagehits;
	}
      printf("Number of hits %d \n", hits);
      printf("Estimate of Pi with %d processors: %f \n"
	     , p,  4.0*(double)(hits)/((double)(Ndarts)*p));
    }

  */
  MPI_Finalize();
  return 0;
}
