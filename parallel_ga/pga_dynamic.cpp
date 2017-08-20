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

using namespace std;

// 
//  Change any of these parameters to match your needs 
//
# define POPSIZE 5
# define MAXGENS 10
# define NVARS 3
# define PXOVER 0.8
# define PMUTATION 0.15

//  Each GENOTYPE is a member of the population, with
//  gene: a string of variables,
//  fitness: the fitness
//  upper: the variable upper bounds,
//  lower: the variable lower bounds,
//  rfitness: the relative fitness,
//  cfitness: the cumulative fitness.
//
//typedef struct genotype_s
//{
//  double* gene;
//  double fitness;
//  double* upper;
//  double* lower;
//  double rfitness;
//  double cfitness;
//} genotype;

struct genotype
{
  double fitness;
  double rfitness;
  double cfitness;
  double* gene;
};




struct genotype* population = NULL;
struct genotype* newpopulation = NULL;

int main (int, char**);
void crossover ( int &seed );
void elitist ( );
void evaluate ( );
int i4_uniform_ab ( int a, int b, int &seed );
void initialize ( string filename, int &seed,  double params_range[NVARS][2]);
void keep_the_best ( );
void mutate ( int &seed );
double r8_uniform_ab ( double a, double b, int &seed );
void report ( int generation );
void selector ( int &seed );
void timestamp ( );
void Xover ( int one, int two, int &seed );
struct genotype* initialize_poparray(int num_var, int popsize);
void copy_gene(genotype &copyTo, genotype &copyFrom, int geneCount);

void node_selector(int my_start, int row_count, int &seed );
void node_crossover(int my_start, int row_count, int &seed );
void node_mutate(int my_start, int row_count, int &seed, double params_range[NVARS][2]);
void node_evaluate(int my_start, int row_count);
void node_elitist(int my_start, int row_count, int gen_best_worst[]);
void master_elitist(int best_worst[][2], int node_count);
void node_Xover ( int one, int two, int &seed );
void print_population();


/* MPI Global variables */
int myId, senderId, receiverId, numRows, numNodes, numRowsToReceive, ierr,
        avgRowsPerNode, numRowsReceived, startRow, endRow,
        numRowsToSend;
MPI_Status status;

#define MASTER_ID 0

/* MPI message tags */
#define SEND_DATA_TAG 2001
#define RETN_DATA_TAG 2002
#define CONT_TAG 2003
#define PART_INDEX_TAG 3001
#define BEST_WORST_TAG 4001

struct nodeParts {
    int start;
    int end;
};


enum boolean {FALSE, TRUE};
boolean moreGen; // Flag used by nodes to continue another loop 

//****************************************************************************80

int main(int argc, char** argv)

//****************************************************************************80
//
//  Purpose:
//
//    MAIN supervises the genetic algorithm.
//
//  Discussion:
//
//    Each generation involves selecting the best 
//    members, performing crossover & mutation and then 
//    evaluating the resulting population, until the terminating 
//    condition is satisfied   
//
//    This is a simple genetic algorithm implementation where the 
//    evaluation function takes positive values only and the 
//    fitness of an individual is the same as the value of the 
//    objective function.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Zbigniew Michalewicz,
//    Genetic Algorithms + Data Structures = Evolution Programs,
//    Third Edition,
//    Springer, 1996,
//    ISBN: 3-540-60676-9,
//    LC: QA76.618.M53.
//
//  Parameters:
//
//    MAXGENS is the maximum number of generations.
//
//    NVARS is the number of problem variables.
//
//    PMUTATION is the probability of mutation.
//
//    POPSIZE is the population size. 
//
//    PXOVER is the probability of crossover.                          
//
{
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int length; // To store character length when required
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &numNodes); // How many nodes
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myId); // My Rank ID
    
    string filename = "simple_ga_input.txt";
    int generation;
    int i;
    int seed;
  
    int num_var = NVARS;
    int popsize = POPSIZE;
    double params_range[NVARS][2];

    population = initialize_poparray(num_var, popsize);
    newpopulation = initialize_poparray(num_var, popsize);
  
  
  /* Create MPI_Type for genotype struct */
    struct genotype a_genotpe;
    const int nitems = 4;
    int          blocklengths[4] = {1, 1, 1, NVARS};
    MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Datatype mpi_genotype;
    MPI_Aint addr1, addr2, addr3, addr4, baseaddr, offsets[4];
    
    MPI_Get_address(&a_genotpe, &baseaddr);
    MPI_Get_address(&a_genotpe.fitness,    &addr1);
//    MPI_Get_address(population->upper,      &addr3);
//    MPI_Get_address(population->lower,      &addr4);
    MPI_Get_address(&a_genotpe.rfitness,      &addr2);
    MPI_Get_address(&a_genotpe.cfitness,      &addr3);
    MPI_Get_address(&(a_genotpe.gene), &addr4);
    
    
//    offsets[0] = offsetof(genotype, gene);
//    offsets[1] = offsetof(genotype, fitness);
//    offsets[2] = offsetof(genotype, upper);
//    offsets[3] = offsetof(genotype, lower);
//    offsets[4] = offsetof(genotype, rfitness);
//    offsets[5] = offsetof(genotype, cfitness);
    
    offsets[0] = 0;
    offsets[1] = addr2 - baseaddr;
    offsets[2] = addr3 - baseaddr;
    offsets[3] = addr4 - baseaddr;
//    offsets[4] = addr5 - baseaddr;
//    offsets[5] = addr6 - baseaddr;
    
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_genotype);
    MPI_Type_commit(&mpi_genotype);

  
   seed = 123456789;
   moreGen = TRUE;
  
  if(myId == MASTER_ID) // Master Process
  {
        timestamp ( );
        cout << "\n";
        cout << "SIMPLE_GA:\n";
        cout << "  C++ version\n";
        cout << "  A simple example of a genetic algorithm.\n";

        if ( NVARS < 2 )
        {
          cout << "\n";
          cout << "  The crossover modification will not be available,\n";
          cout << "  since it requires 2 <= NVARS.\n";
        }

        initialize ( filename, seed, params_range);

        evaluate ( );
        cout << "Evaluation is done\n";
        keep_the_best ( );
        cout << "  Keep the best done.\n";
        
//        cout << " ORIGINAL POPULATION " << endl;
//        print_population();
//        cout << " =================== " << endl;
        
        int nodeId; // Variable used to address nodes
        
        /* Create a data structure to store data partition values for each nodes
         * 
         * nodes_data[node][0]: start_index
         * nodes_data[node][1]: end_index
         * nodes_data[node][2]: row_count
         */
        int nodes_data[numNodes - 1][4];
        
        /* Calculate and store start-end index for each nodes */
        int num_chunk = numNodes - 1;
        int start_i, end_i, num_rows;
        int chunk_size = (int)(ceil((double)POPSIZE / num_chunk));
        
        /* Array to store nodes best-worst index */
        int best_worst[numNodes - 1][2];
        
//        cout << "Population size: " << POPSIZE << endl;
        
        for(int i = 0; i < num_chunk; i++)
        {
            start_i = i * chunk_size;
            num_rows = min(chunk_size, POPSIZE - start_i);
            end_i = (start_i + num_rows) - 1;
            nodes_data[i][0] = start_i;
            nodes_data[i][1] = end_i;
            nodes_data[i][2] = num_rows;
            
//            cout << "Node [" << i + 1 << "] - " << " Row count: " << num_rows 
//                    << " Start: " << start_i << " End: " << end_i << endl;
        }
        
        
        /* Send each node its start and end index */
        for(nodeId = 1; nodeId < numNodes; nodeId++)
        {
            moreGen = TRUE;
            MPI_Send(&nodes_data[nodeId - 1], 3, MPI_INT, nodeId, PART_INDEX_TAG, 
                    MPI_COMM_WORLD);
        }
        
               
        /* Run GA ops for given number of generation */
        for(generation = 0; generation < MAXGENS; generation++)
        {
            cout << "Master initiating generation: " << generation << endl;
            
            /* Master controls generation cycles */
            for(nodeId = 1; nodeId < numNodes; nodeId++)
            {
                moreGen = TRUE;
                MPI_Send(&moreGen, 1, MPI_INT, nodeId, CONT_TAG, MPI_COMM_WORLD);
            }
            
            /* Update all nodes with current population */
            for(nodeId = 1; nodeId < numNodes; nodeId++)
            {
                MPI_Send(population, POPSIZE + 1, mpi_genotype, nodeId, SEND_DATA_TAG,
                        MPI_COMM_WORLD);
            }
            cout << "Population data sent." << endl;
            
            /* Collect all parts of new population from each node */
            int start_index, item_count;
            for(nodeId = 1; nodeId < numNodes; nodeId++)
            {
                start_index = nodes_data[nodeId - 1][0];
                item_count = nodes_data[nodeId - 1][2];
                MPI_Recv(&population[start_index], item_count, mpi_genotype,
                        nodeId, RETN_DATA_TAG, MPI_COMM_WORLD, &status);
            }
            
            /* Collect best and worst member index from each node and update
             * list */
            for(nodeId = 1; nodeId < numNodes; nodeId++)
            {
                MPI_Recv(&best_worst[nodeId - 1], 2, MPI_INT,
                        nodeId, BEST_WORST_TAG, MPI_COMM_WORLD, &status);
            }
            
        
            /* With update best_worst list perform elitist() on entire population */
            master_elitist(best_worst, numNodes - 1);
            
//            cout << " AFTER GENERATION: " << generation << endl;
//            print_population();
//            cout << " =================== " << endl;
            
            /* Master generate report  */
            report(generation);
        }
        
        
        /* Send GA termination signal */
        moreGen = FALSE;
        for(nodeId = 1; nodeId < numNodes; nodeId++)
        {
            MPI_Send(&moreGen, 1, MPI_INT, nodeId, CONT_TAG, MPI_COMM_WORLD);
//            cout << "Master: sending end token to node: " << nodeId << endl;
        }
        
        
        cout << "\n";
        cout << "  Best member after " << MAXGENS << " generations:\n";
        cout << "\n";

        for ( i = 0; i < NVARS; i++ )
        {
          cout << "  var(" << i << ") = " << population[POPSIZE].gene[i] << "\n";
        }

        cout << "\n";
        cout << "  Best fitness = " << population[POPSIZE].fitness << "\n";
      //
      //  Terminate.
      //
        cout << "\n";
        cout << "SIMPLE_GA:\n";
        cout << "  Normal end of execution.\n";
        cout << "\n";
        timestamp ( );
  }
  else // Slave processes
  {  
      int my_start;
      int my_end;
      int row_count;
      int gen_best_worst[2]; // To send best-worst index of this part to master
      int part_data[3]; // Holds start, end indices and row_count
      /* Get my hostname [Not required ] */
      MPI_Get_processor_name(hostname, &length);
      
      /* First get my start and end index from Master */ 
      ierr = MPI_Recv(&part_data, 3, MPI_INT, MASTER_ID, PART_INDEX_TAG,
              MPI_COMM_WORLD, &status);
      if(ierr == MPI_SUCCESS)
      {
          my_start = part_data[0];
          my_end = part_data[1];
          row_count = part_data[2];
          
//          cout << hostname << "[" << myId << "] - " << " Row count: " 
//                  << row_count << " Start: " << my_start << " End: "
//                  << my_end << endl;
      }
      
      /* Is there a generation to calculate? Check with Master */
      MPI_Recv(&moreGen, 1, MPI_INT, MASTER_ID, CONT_TAG, MPI_COMM_WORLD,
              &status);
      generation = 0; 
      while(moreGen)
      {
//          cout << hostname << " :" << myId << " processing generation: " << generation++
//                  << endl;
          /* Get updated entire population from master */
          ierr = MPI_Recv(population, POPSIZE+1, mpi_genotype, MASTER_ID, 
                  SEND_DATA_TAG, MPI_COMM_WORLD, &status);
          
          if(ierr == MPI_SUCCESS)
          {
//              cout << hostname << "[" << myId << "] - "
//                      << " Population[] updated" << endl;
                            
             /* Node relevant functions to replace above ( need to implement ) */
                node_selector (my_start, row_count, seed );
                node_crossover (my_start, row_count, seed );                
                node_mutate ( my_start, row_count,seed, params_range );
                node_evaluate (my_start, row_count);
                node_elitist (my_start, row_count, gen_best_worst);
                
            /* Send my part of new population */
            MPI_Send(&population[my_start], row_count, mpi_genotype, MASTER_ID,
                    RETN_DATA_TAG, MPI_COMM_WORLD);
          
            /* Send Best and worst members index of my part */
            MPI_Send(gen_best_worst, 2, MPI_INT, MASTER_ID,
                    BEST_WORST_TAG, MPI_COMM_WORLD);
            
          }
          
          /* Loop sentinel check from master process */
          MPI_Recv(&moreGen, 1, MPI_INT, MASTER_ID, CONT_TAG, MPI_COMM_WORLD,
              &status);
      }
//      cout << "Node: " << myId  << " ended GA processing " << endl;
  }
   
   /* Free MPI Datatype genotype */
   MPI_Type_free(&mpi_genotype);
   /* Terminate MPI */
   MPI_Finalize();
   return 0;
}
//****************************************************************************80

void crossover ( int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    CROSSOVER selects two parents for the single point crossover.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int FIRST, is a count of the number of members chosen.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  const double a = 0.0;
  const double b = 1.0;
  int mem;
  int one;
  int first = 0;
  double x;

  for ( mem = 0; mem < POPSIZE; ++mem )
  {
    x = r8_uniform_ab ( a, b, seed );

    if ( x < PXOVER )
    {
      ++first;

      if ( first % 2 == 0 )
      {
        Xover ( one, mem, seed );
      }
      else
      {
        one = mem;
      }

    }
  }
  return;
}
//****************************************************************************80

void elitist ( )

//****************************************************************************80
// 
//  Purpose:
//
//    ELITIST stores the best member of the previous generation.
//
//  Discussion:
//
//    The best member of the previous generation is stored as 
//    the last in the array. If the best member of the current 
//    generation is worse then the best member of the previous 
//    generation, the latter one would replace the worst member 
//    of the current population.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double BEST, the best fitness value.
//
//    Local, double WORST, the worst fitness value.
//
{
  int i;
  double best;
  int best_mem;
  double worst;
  int worst_mem;

  best = population[0].fitness;
  worst = population[0].fitness;

  for ( i = 0; i < POPSIZE - 1; ++i )
  {
    if ( population[i+1].fitness < population[i].fitness )
    {

      if ( best <= population[i].fitness )
      {
        best = population[i].fitness;
        best_mem = i;
      }

      if ( population[i+1].fitness <= worst )
      {
        worst = population[i+1].fitness;
        worst_mem = i + 1;
      }

    }
    else
    {

      if ( population[i].fitness <= worst )
      {
        worst = population[i].fitness;
        worst_mem = i;
      }

      if ( best <= population[i+1].fitness )
      {
        best = population[i+1].fitness;
        best_mem = i + 1;
      }

    }

  }
// 
//  If the best individual from the new population is better than 
//  the best individual from the previous population, then 
//  copy the best from the new population; else replace the 
//  worst individual from the current population with the 
//  best one from the previous generation                     
//
  if ( population[POPSIZE].fitness <= best )
  {
    for ( i = 0; i < NVARS; i++ )
    {
      population[POPSIZE].gene[i] = population[best_mem].gene[i];
    }
    population[POPSIZE].fitness = population[best_mem].fitness;
  }
  else
  {
    for ( i = 0; i < NVARS; i++ )
    {
      population[worst_mem].gene[i] = population[POPSIZE].gene[i];
    }
    population[worst_mem].fitness = population[POPSIZE].fitness;
  } 

  return;
}
//****************************************************************************80

void evaluate ( )

//****************************************************************************80
// 
//  Purpose:
//
//    EVALUATE implements the user-defined valuation function
//
//  Discussion:
//
//    Each time this is changed, the code has to be recompiled.
//    The current function is:  x[1]^2-x[1]*x[2]+x[3]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
{
  int member;
  int i;
  double x[NVARS+1];

  for ( member = 0; member < POPSIZE; member++ )
  {
    for ( i = 0; i < NVARS; i++ )
    {
      x[i+1] = population[member].gene[i];
    } 
    population[member].fitness = ( x[1] * x[1] ) - ( x[1] * x[2] ) + x[3];
  }
  return;
}
//****************************************************************************80

int i4_uniform_ab ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}
//****************************************************************************80

void initialize ( string filename, int &seed, double params_range[NVARS][2] )

//****************************************************************************80
// 
//  Purpose:
//
//    INITIALIZE initializes the genes within the variables bounds. 
//
//  Discussion:
//
//    It also initializes (to zero) all fitness values for each
//    member of the population. It reads upper and lower bounds 
//    of each variable from the input file `gadata.txt'. It 
//    randomly generates values between these bounds for each 
//    gene of each genotype in the population. The format of 
//    the input file `gadata.txt' is 
//
//      var1_lower_bound var1_upper bound
//      var2_lower_bound var2_upper bound ...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the input file.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  int i;
  ifstream input;
  int j;
  double lbound;
  double ubound;
  double value;

  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "INITIALIZE - Fatal error!\n";
    cerr << "  Cannot open the input file!\n";
    exit ( 1 );
  }
// 
//  Initialize variables within the bounds 
//
  for ( i = 0; i < NVARS; i++ )
  {
    input >> lbound >> ubound;
    params_range[i][0] = lbound;
    params_range[i][1] = ubound;
    
    for ( j = 0; j < POPSIZE; j++ )
    {
      population[j].fitness = 0;
      population[j].rfitness = 0;
      population[j].cfitness = 0;
//      population[j].lower[i] = lbound;
//      population[j].upper[i]= ubound;
      population[j].gene[i] = r8_uniform_ab ( lbound, ubound, seed );
    }
  }

  input.close ( );

  return;
}
//****************************************************************************80

void keep_the_best ( )

//****************************************************************************80
// 
//  Purpose:
//
//    KEEP_THE_BEST keeps track of the best member of the population. 
//
//  Discussion:
//
//    Note that the last entry in the array Population holds a 
//    copy of the best individual.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int CUR_BEST, the index of the best individual.
//
{     
  int cur_best;
  int mem;
  int i;

  cur_best = 0;

  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    if ( population[POPSIZE].fitness < population[mem].fitness )
    {
      cur_best = mem;
      population[POPSIZE].fitness = population[mem].fitness;
    }
  }
// 
//  Once the best member in the population is found, copy the genes.
//
//  for ( i = 0; i < NVARS; i++ )
//  {
//    population[POPSIZE].gene[i] = population[cur_best].gene[i];
//  }
  copy_gene(population[POPSIZE], population[cur_best], NVARS);
  return;
}
//****************************************************************************80

void mutate ( int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    MUTATE performs a random uniform mutation. 
//
//  Discussion:
//
//    A variable selected for mutation is replaced by a random value 
//    between the lower and upper bounds of this variable.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  double lbound;
  double ubound;
  double x;

  for ( i = 0; i < POPSIZE; i++ )
  {
    for ( j = 0; j < NVARS; j++ )
    {
      x = r8_uniform_ab ( a, b, seed );
      if ( x < PMUTATION )
      {
//        lbound = population[i].lower[j];
//        ubound = population[i].upper[j];
//        population[i].gene[j] = r8_uniform_ab ( lbound, ubound, seed );
      }
    }
  }

  return;
}
//****************************************************************************80

double r8_uniform_ab ( double a, double b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_AB returns a scaled pseudorandom R8.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_AB, a number strictly between A and B.
//
{
  int i4_huge = 2147483647;
  int k;
  double value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  value = ( double ) ( seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
}
//****************************************************************************80

void report ( int generation )

//****************************************************************************80
// 
//  Purpose:
//
//    REPORT reports progress of the simulation. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double avg, the average population fitness.
//
//    Local, best_val, the best population fitness.
//
//    Local, double square_sum, square of sum for std calc.
//
//    Local, double stddev, standard deviation of population fitness.
//
//    Local, double sum, the total population fitness.
//
//    Local, double sum_square, sum of squares for std calc.
//
{
  double avg;
  double best_val;
  int i;
  double square_sum;
  double stddev;
  double sum;
  double sum_square;

  if ( generation == 0 )
  {
    cout << "\n";
    cout << "  Generation       Best            Average       Standard \n";
    cout << "  number           value           fitness       deviation \n";
    cout << "\n";
  }

  sum = 0.0;
  sum_square = 0.0;

  for ( i = 0; i < POPSIZE; i++ )
  {
    sum = sum + population[i].fitness;
    sum_square = sum_square + population[i].fitness * population[i].fitness;
  }

  avg = sum / ( double ) POPSIZE;
  square_sum = avg * avg * POPSIZE;
  stddev = sqrt ( ( sum_square - square_sum ) / ( POPSIZE - 1 ) );
  best_val = population[POPSIZE].fitness;

  cout << "  " << setw(8) << generation 
       << "  " << setw(14) << best_val 
       << "  " << setw(14) << avg 
       << "  " << setw(14) << stddev << "\n";

  return;
}
//****************************************************************************80

void selector ( int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    SELECTOR is the selection function.
//
//  Discussion:
//
//    Standard proportional selection for maximization problems incorporating 
//    the elitist model.  This makes sure that the best member always survives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  int mem;
  double p;
  double sum;
//
//  Find the total fitness of the population.
//
  sum = 0.0;
  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    sum = sum + population[mem].fitness;
  }
//
//  Calculate the relative fitness of each member.
//
  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    population[mem].rfitness = population[mem].fitness / sum;
  }
// 
//  Calculate the cumulative fitness.
//
  population[0].cfitness = population[0].rfitness;
  for ( mem = 1; mem < POPSIZE; mem++ )
  {
    population[mem].cfitness = population[mem-1].cfitness +       
      population[mem].rfitness;
  }
// 
//  Select survivors using cumulative fitness. 
//
  for ( i = 0; i < POPSIZE; i++ )
  { 
    p = r8_uniform_ab ( a, b, seed );
    if ( p < population[0].cfitness )
    {
      newpopulation[i] = population[0];      
    }
    else
    {
      for ( j = 0; j < POPSIZE; j++ )
      { 
        if ( population[j].cfitness <= p && p < population[j+1].cfitness )
        {
          newpopulation[i] = population[j+1];
        }
      }
    }
  }
// 
//  Overwrite the old population with the new one.
//
  for ( i = 0; i < POPSIZE; i++ )
  {
    population[i] = newpopulation[i]; 
  }

  return;     
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2003
//
//  Author:
//
//    John Burkardt
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

void Xover ( int one, int two, int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    XOVER performs crossover of the two selected parents. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int point, the crossover point.
//
//  Parameters:
//
//    Input, int ONE, TWO, the indices of the two parents.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  int i;
  int point;
  double t;
// 
//  Select the crossover point.
//
  point = i4_uniform_ab ( 0, NVARS - 1, seed );
//
//  Swap genes in positions 0 through POINT-1.
//
  for ( i = 0; i < point; i++ )
  {
    t                       = population[one].gene[i];
    population[one].gene[i] = population[two].gene[i];
    population[two].gene[i] = t;
  }

  return;
}

/* Initialize Population Array */
struct genotype* initialize_poparray(int num_var, int popsize)
{        
    struct genotype *poparray = (struct genotype *)malloc(popsize * sizeof(
            struct genotype) + 1);
    for(int i = 0; i < popsize + 1; i++)
    {
        poparray[i].gene = (double*)malloc(sizeof(double) * num_var);
//        poparray[i].upper = (double*)malloc(sizeof(double) * num_var);
//        poparray[i].lower = (double*)malloc(sizeof(double) * num_var);
    }

    return poparray;
    
}

void copy_gene(genotype &copyTo, genotype &copyFrom, int geneCount){
    
    cout << "Copyto fitness: " << copyTo.gene[0] << endl;
    for(int i = 0; i < geneCount; i++)
    {
        copyTo.gene[i] = copyFrom.gene[i];
        cout << "Copying gene\n";
    }
    copyTo.fitness = copyFrom.fitness;
}

void node_crossover(int my_start, int row_count,int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    CROSSOVER selects two parents for the single point crossover.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int FIRST, is a count of the number of members chosen.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  const double a = 0.0;
  const double b = 1.0;
  int mem;
  int one;
  int first = 0;
  double x;

  for ( mem = my_start; mem < (my_start + row_count); ++mem )
  {
    x = r8_uniform_ab ( a, b, seed );

    if ( x < PXOVER )
    {
      ++first;

      if ( first % 2 == 0 )
      {
        Xover ( one, mem, seed );
      }
      else
      {
        one = mem;
      }

    }
  }
  return;
}
//****************************************************************************80

void node_elitist(int my_start, int row_count, int gen_best_worst[])

//****************************************************************************80
// 
//  Purpose:
//
//    ELITIST stores the best member of the previous generation.
//
//  Discussion:
//
//    The best member of the previous generation is stored as 
//    the last in the array. If the best member of the current 
//    generation is worse then the best member of the previous 
//    generation, the latter one would replace the worst member 
//    of the current population.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double BEST, the best fitness value.
//
//    Local, double WORST, the worst fitness value.
//  Parameters:
//
//    Input, int my_start, starting index of population array
//    Input, int row_count, Number of rows to be processed
//    Output, int gen_best_worst, index of local best[0] and worst[1] member
//
{
  int i;
  double best;
  int best_mem;
  double worst;
  int worst_mem;
  int ends = my_start + row_count;

  best = population[my_start].fitness;
  worst = population[my_start].fitness;
  best_mem = my_start;
  worst_mem = my_start;

  for ( i = my_start; i < ends - 1; ++i )
  //for ( i = 0; i < POPSIZE - 1; ++i )
  {
    if ( population[i+1].fitness < population[i].fitness )
    {

      if ( best <= population[i].fitness )
      {
        best = population[i].fitness;
        best_mem = i;
      }

      if ( population[i+1].fitness <= worst )
      {
        worst = population[i+1].fitness;
        worst_mem = i + 1;
      }

    }
    else
    {

      if ( population[i].fitness <= worst )
      {
        worst = population[i].fitness;
        worst_mem = i;
      }

      if ( best <= population[i+1].fitness )
      {
        best = population[i+1].fitness;
        best_mem = i + 1;
      }

    }

  }
  
  gen_best_worst[0] = best_mem;
  gen_best_worst[1] = worst_mem;
//  cout << endl << "Best index: " << best_mem << " : "
//          << population[best_mem].fitness
//          << " Worst: " << worst_mem << " : " 
//          << population[worst_mem].fitness<< endl;
  return;
}
//****************************************************************************80

void master_elitist(int best_worst[][2], int node_count)

//****************************************************************************80
// 
//  Purpose:
//
//    ELITIST stores the best member of the previous generation.
//
//  Discussion:
//
//    The best member of the previous generation is stored as 
//    the last in the array. If the best member of the current 
//    generation is worse then the best member of the previous 
//    generation, the latter one would replace the worst member 
//    of the current population.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double BEST, the best fitness value.
//
//    Local, double WORST, the worst fitness value.
//
// Parameters:
//  
//      Input int best_worst[][], array containing best and worst index for 
//          each node.
//      Input int node_count, number of nodes.
//
{
  int i;
  double best;
  int best_mem;
  double worst;
  int worst_mem;

  best = population[best_worst[0][0]].fitness;
  worst = population[best_worst[0][1]].fitness;
  best_mem = best_worst[0][0];
  worst_mem = best_worst[0][1];
  
  
  for(int node = 1; node < node_count; node++)
  {
      if(best < population[ best_worst[node][0] ].fitness)
      {
          best = population[ best_worst[node][0] ].fitness;
          best_mem = best_worst[node][0];
      }
      
      if(worst > population[ best_worst[node][1] ].fitness)
      {
          worst = population[ best_worst[node][1] ].fitness;
          worst_mem = best_worst[node][1];
      }
  }

// 
//  If the best individual from the new population is better than 
//  the best individual from the previous population, then 
//  copy the best from the new population; else replace the 
//  worst individual from the current population with the 
//  best one from the previous generation                     
//
  if ( population[POPSIZE].fitness <= best )
  {
    for ( i = 0; i < NVARS; i++ )
    {
      population[POPSIZE].gene[i] = population[best_mem].gene[i];
    }
    population[POPSIZE].fitness = population[best_mem].fitness;
//    cout << "New best update with index: " << best_mem << endl;
  }
  else
  {
    for ( i = 0; i < NVARS; i++ )
    {
      population[worst_mem].gene[i] = population[POPSIZE].gene[i];
    }
    population[worst_mem].fitness = population[POPSIZE].fitness;
  } 

  return;
}

void node_evaluate(int my_start, int row_count)

//****************************************************************************80
// 
//  Purpose:
//
//    EVALUATE implements the user-defined valuation function
//
//  Discussion:
//
//    Each time this is changed, the code has to be recompiled.
//    The current function is:  x[1]^2-x[1]*x[2]+x[3]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
{
  int member;
  int i;
  double x[NVARS+1];
  int ends = my_start + row_count;
  
  for ( member = my_start; member < ends; member++ )
  {
    for ( i = 0; i < NVARS; i++ )
    {
      x[i+1] = population[member].gene[i];
    } 
    population[member].fitness = ( x[1] * x[1] ) - ( x[1] * x[2] ) + x[3];
  }
  return;
}
//****************************************************************************80

void node_mutate(int my_start, int row_count, int &seed, 
        double params_range[NVARS][2] )

//****************************************************************************80
// 
//  Purpose:
//
//    MUTATE performs a random uniform mutation. 
//
//  Discussion:
//
//    A variable selected for mutation is replaced by a random value 
//    between the lower and upper bounds of this variable.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  double lbound;
  double ubound;
  double x;
  int ends = my_start + row_count;
  
  for ( i = my_start; i < ends; i++ )
  {
    for ( j = 0; j < NVARS; j++ )
    {
      x = r8_uniform_ab ( a, b, seed );
      if ( x < PMUTATION )
      {
//        lbound = population[i].lower[j];
//        ubound = population[i].upper[j];
        lbound = params_range[j][0];
        ubound = params_range[j][1]; 
        population[i].gene[j] = r8_uniform_ab ( lbound, ubound, seed );
      }
    }
  }

  return;
}
//****************************************************************************80
void node_selector (int my_start, int row_count,int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    SELECTOR is the selection function.
//
//  Discussion:
//
//    Standard proportional selection for maximization problems incorporating 
//    the elitist model.  This makes sure that the best member always survives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  int mem;
  double p;
  double sum;
//
//  Find the total fitness of the population.
//
  sum = 0.0;
  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    sum = sum + population[mem].fitness;
  }
//
//  Calculate the relative fitness of each member.
//
  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    population[mem].rfitness = population[mem].fitness / sum;
  }
// 
//  Calculate the cumulative fitness.
//
  population[0].cfitness = population[0].rfitness;
  for ( mem = 1; mem < POPSIZE; mem++ )
  {
    population[mem].cfitness = population[mem-1].cfitness +       
      population[mem].rfitness;
  }
// 
//  Select survivors using cumulative fitness. 
//
  int ends = my_start + row_count;
  for ( i = my_start; i < ends; i++ )
  { 
    p = r8_uniform_ab ( a, b, seed );
    if ( p < population[0].cfitness )
    {
      newpopulation[i] = population[0];      
    }
    else
    {
      for ( j = 0; j < POPSIZE; j++ )
      { 
        if ( population[j].cfitness <= p && p < population[j+1].cfitness )
        {
          newpopulation[i] = population[j+1];
        }
      }
    }
  }
// 
//  Overwrite the old population with the new one.
//
  for ( i = my_start; i < ends; i++ )
  {
    population[i] = newpopulation[i]; 
  }

  return;     
}
//****************************************************************************80

void node_Xover ( int one, int two, int &seed )

//****************************************************************************80
// 
//  Purpose:
//
//    XOVER performs crossover of the two selected parents. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int point, the crossover point.
//
//  Parameters:
//
//    Input, int ONE, TWO, the indices of the two parents.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
  int i;
  int point;
  double t;
// 
//  Select the crossover point.
//
  point = i4_uniform_ab ( 0, NVARS - 1, seed );
//
//  Swap genes in positions 0 through POINT-1.
//
  for ( i = 0; i < point; i++ )
  {
    t                       = population[one].gene[i];
    population[one].gene[i] = population[two].gene[i];
    population[two].gene[i] = t;
  }

  return;
}

void print_population()
{
//    cout.setf(ios::fixed, ios::floatfield);
//    cout.setf(ios::showpoint);
    for(int i = 0; i < POPSIZE + 1; i++)
    {
        cout << " | ";
        for ( int j = 0; j < NVARS; j++ )
        {
          cout << "[" << population[i].gene[j] << "] ";
        }
        cout << " | ";
    }
    cout << endl;
}
