# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <mpi.h>
# include <stddef.h>
# include "params.hpp"
# include "sga_proc.hpp"
# include "pga_proc.hpp"

using namespace std;

genotype population[POPSIZE+1];
genotype newpopulation[POPSIZE+1];


/* MPI Global variables */
int myId, senderId, receiverId, numRows, numNodes, numRowsToReceive, ierr,
        avgRowsPerNode, numRowsReceived, startRow, endRow,
        numRowsToSend;
MPI_Status status;

//*******************************************************
int main(int argc, char** argv)
{
    clock_t start, end;
    double cpu_time_used;
    start = clock();
    
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int length; // To store character length when required
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &numNodes); // How many nodes
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myId); // My Rank ID
   
  string filename;
  
    switch(FITNESS_F) {
      case F1:
          filename = "simple_ga_input.txt";
          break;
      case F8:
          filename = "simple_ga_input2.txt";
          break;
      default:
          filename = "simple_ga_input.txt";
          break;
  }
    
  int generation;
  int i;
  int seed;
  
  /* Create MPI_Type for genotype struct */
    MPI_Datatype mpi_genotype;
    create_mpi_genotype_struct(mpi_genotype);

  
   seed = 123456789;
   
   
   /* Initialise first gen population in master process */
   
    if(myId == MASTER_ID) // Master Process
    {        
       timestamp ( );
       cout << "\n";
       cout << "PARALLEL_GA:\n";
       cout << "  A simple example of a genetic algorithm.\n";

       if ( NVARS < 2 )
       {
         cout << "\n";
         cout << "  The crossover modification will not be available,\n";
         cout << "  since it requires 2 <= NVARS.\n";
       }
        
        initialize ( filename, seed );
        evaluate (0, POPSIZE );
        keep_the_best ( );
    }
   
   /* Sync message - 1 */
   MPI_Barrier(MPI_COMM_WORLD);
   
   /* Broadcast entire first generation from master process */
   MPI_Bcast(population, POPSIZE + 1, mpi_genotype,
        MASTER_ID, MPI_COMM_WORLD);
   
   /* Sync message - 2 */
   MPI_Barrier(MPI_COMM_WORLD);
   
   
    /* Calculate and store start-end index for each nodes */
    int nodes_data[numNodes][2];
    int start_i, num_rows;
    int chunk_size = (int)(ceil((double) POPSIZE / numNodes));
    
            
    for(int i = 0; i < numNodes; i++)
    {
        start_i = i * chunk_size;
        num_rows = min(chunk_size, POPSIZE - start_i);
        nodes_data[i][0] = start_i;
        nodes_data[i][1] = num_rows;
    }
    
    
    /* Each process does GA ops here on individual data parts */
    
    int my_start = nodes_data[myId][0];
    int row_count = nodes_data[myId][1];
    for ( generation = 0; generation < MAXGENS; generation++ )
    {      
        node_selector (my_start, row_count, seed );
        crossover ( seed );
        node_mutate ( my_start, row_count,seed );   
        evaluate (my_start, row_count);
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        for(int i = 0; i < numNodes; i++)
        {
            MPI_Bcast(&population[nodes_data[i][0]], nodes_data[i][1], 
                    mpi_genotype, i, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        if(myId == MASTER_ID) // Master Process
        {
           report(generation);
            elitist ( );
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
   MPI_Barrier(MPI_COMM_WORLD);
   
   /* Master process prints final report */
   
  if(myId == MASTER_ID) // Master Process
  {       
       cout << "\n";
       cout << "  Best member after " << MAXGENS << " generations:\n";
       cout << "\n";

       for ( i = 0; i < NVARS; i++ )
       {
         cout << "  var(" << i << ") = " << population[POPSIZE].gene[i] << "\n";
       }

       cout << "\n";
       cout << "  Best fitness = " << population[POPSIZE].fitness << "\n";
     
       cout << "\n";
       cout << "SIMPLE_GA:\n";
       cout << "  Normal end of execution.\n";
       cout << "\n";
       timestamp ( );
        
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        printf("%f\n",cpu_time_used);
  }
   
   /* Free MPI Datatype genotype */
   MPI_Type_free(&mpi_genotype);

   /* Terminate MPI */
   MPI_Finalize();
   return 0;
}