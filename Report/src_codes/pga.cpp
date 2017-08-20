# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <mpi.h>
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

//**********************************************
int main(int argc, char** argv)
{
    
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
  
  int seed;
  
  /* Create MPI_Type for genotype struct */
    MPI_Datatype mpi_genotype;
    create_mpi_genotype_struct(mpi_genotype);
    
   seed = 123456789;
   
    /* Create a data structure to store data partition values for each nodes
     * 
     * nodes_data[node][0]: start_index
     * nodes_data[node][1]: rowCount
     */
    int nodes_data[numNodes][2];
    /* Calculate and store start-end index for each nodes */
    int num_chunk = numNodes;
    int start_i, rowCount;
    int chunk_size = (int)(ceil((double)POPSIZE / num_chunk));

    for(int i = 0; i < num_chunk; i++)
    {
        start_i = i * chunk_size;
        rowCount = min(chunk_size, POPSIZE - start_i);
        nodes_data[i][0] = start_i;
        nodes_data[i][1] = rowCount;      
    }
   
  if(myId == MASTER_ID) // Master Process
  {
        master_process(filename, seed, mpi_genotype, nodes_data, start);
  }
  else // Slave processes
  {  
      slave_process(seed, mpi_genotype, nodes_data);
  }
   
   /* Free MPI Datatype genotype */
   MPI_Type_free(&mpi_genotype);
   
   /* Terminate MPI */
   MPI_Finalize();
   return 0;
}