#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include <stdlib.h>
#include <math.h>

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

int main(int argc, char **argv)
{
  int my_rank;
  int p;       //number of processes
  int source;  //rank of sender
  int dest;    //rank of receiver
  int tag = 50;//tag of messages
  char message[100]; //storage for the message
  MPI_Status status; //Return status for receive

  MPI_Init(&argc, &argv); 
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);


  
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


  MPI_Finalize();
  return 0;
}
