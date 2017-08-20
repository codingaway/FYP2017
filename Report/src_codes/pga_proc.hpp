# include "params.hpp"
# include "sga_proc.hpp"

#ifndef PGA_PROC_HPP

#define MASTER_ID 0

/* MPI message tags */
#define SEND_DATA_TAG 2001
#define RETN_DATA_TAG 2002
#define CONT_TAG 2003
#define PART_INDEX_TAG 3001
#define BEST_WORST_TAG 4001

enum boolean {FALSE, TRUE};

typedef struct
{
  double gene[NVARS];
  double fitness;
  double upper[NVARS];
  double lower[NVARS];
  double rfitness;
  double cfitness;
} genotype;
extern genotype population[];
extern newpopulation[];

/* MPI Global variables */
extern int myId, senderId, receiverId, numRows, numNodes, numRowsToReceive, ierr,
        avgRowsPerNode, numRowsReceived, startRow, endRow,
        numRowsToSend;
extern MPI_Status status;

void master_process(string& filename, int &seed,
        MPI_Datatype &mpi_genotype, int node_data[][2],
        clock_t start);
void slave_process( int &seed, MPI_Datatype& mpi_genotype,
        int node_data[][2]);
void create_mpi_genotype_struct(MPI_Datatype& mpi_genotype);
void evaluate ( int my_start, int row_count );
void f1(int my_start, int row_count);
void f8(int my_start, int row_count);
void node_selector(int my_start, int row_count, int &seed );
void node_crossover(int my_start, int row_count, int &seed );
void node_mutate(int my_start, int row_count, int &seed );
void node_evaluate(int my_start, int row_count);
void node_Xover ( int one, int two, int &seed );
void print_population();
#endif /* PGA_PROC_HPP */