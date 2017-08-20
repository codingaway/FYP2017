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
# include "defines.hpp"

using namespace std;

//  Each GENOTYPE is a member of the population, with
//  gene: a string of variables,
//  fitness: the fitness
//  upper: the variable upper bounds,
//  lower: the variable lower bounds,
//  rfitness: the relative fitness,
//  cfitness: the cumulative fitness.
//
typedef struct genotype_s
{
  double gene[NVARS];
  double fitness;
  double upper[NVARS];
  double lower[NVARS];
  double rfitness;
  double cfitness;
} genotype;

genotype population[POPSIZE+1];
genotype newpopulation[POPSIZE+1];

int main (int, char**);
void crossover ( int &seed );
void elitist ( );
void evaluate ( int my_start, int row_count );
int i4_uniform_ab ( int a, int b, int &seed );
void initialize ( string filename, int &seed );
void keep_the_best ( );
void mutate ( int &seed );
double r8_uniform_ab ( double a, double b, int &seed );
void report ( int generation );
void selector ( int &seed );
void timestamp ( );
void Xover ( int one, int two, int &seed );
void f1(int my_start, int row_count);
void f8(int my_start, int row_count);

/* Parallel procedures */
void node_selector(int my_start, int row_count, int &seed );
void node_crossover(int my_start, int row_count, int &seed );
void node_mutate(int my_start, int row_count, int &seed );
void node_evaluate(int my_start, int row_count);
void node_elitist(int my_start, int row_count, int gen_best_worst[]);
void master_elitist(int best_worst[][2], int node_count);
void node_Xover ( int one, int two, int &seed );
void print_population();
void create_mpi_genotype_struct(MPI_Datatype& mpi_genotype);

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
//        timestamp ( );
//        cout << "\n";
//        cout << "PARALLEL_GA:\n";
//      //  cout << "  C++ version\n";
//        cout << "  A simple example of a genetic algorithm.\n";
//
//        if ( NVARS < 2 )
//        {
//          cout << "\n";
//          cout << "  The crossover modification will not be available,\n";
//          cout << "  since it requires 2 <= NVARS.\n";
//        }

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
        node_crossover (my_start, row_count, seed );
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
//            report(generation);
            elitist ( );
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
   MPI_Barrier(MPI_COMM_WORLD);
   
   /* Master process prints final report */
   
  if(myId == MASTER_ID) // Master Process
  {       
//        cout << "\n";
//        cout << "  Best member after " << MAXGENS << " generations:\n";
//        cout << "\n";
//
//        for ( i = 0; i < NVARS; i++ )
//        {
//          cout << "  var(" << i << ") = " << population[POPSIZE].gene[i] << "\n";
//        }
//
//        cout << "\n";
//        cout << "  Best fitness = " << population[POPSIZE].fitness << "\n";
//      
//        cout << "\n";
//        cout << "PARALLEL_GA:\n";
//        cout << "  Normal end of execution.\n";
//        cout << "\n";
//        timestamp ( );
        
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        printf("%f\n",cpu_time_used);
  }
   
   /* Free MPI Datatype genotype */
   MPI_Type_free(&mpi_genotype);
   
   /* terminate log */
//   MPE_Finish_log(cpilog.log);
   /* Terminate MPI */
   MPI_Finalize();
   return 0;
}
//****************************************************************************80

void create_mpi_genotype_struct(MPI_Datatype& mpi_genotype){
    const int nitems = 6;
    int          blocklengths[6] = {NVARS, 1, NVARS, NVARS, 1, 1};
    MPI_Datatype types[6] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
    MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint     offsets[6];
    offsets[0] = offsetof(genotype, gene);
    offsets[1] = offsetof(genotype, fitness);
    offsets[2] = offsetof(genotype, upper);
    offsets[3] = offsetof(genotype, lower);
    offsets[4] = offsetof(genotype, rfitness);
    offsets[5] = offsetof(genotype, cfitness);
    
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_genotype);
    MPI_Type_commit(&mpi_genotype);
}

void crossover ( int &seed )
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

void evaluate ( int my_start, int row_count )
{
//    return FITNESS_F == F1 ? f1() : f8();
    switch(FITNESS_F) {
        case F1:
            return f1(my_start, row_count);
            break;
        case F8:
            return f8(my_start, row_count);
            break;
        default:
            break;
    }
}


void f1(int my_start, int row_count)
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

/*
 * Griewank Function
 */
void f8(int my_start, int row_count)
{
  int member;
  int i;
  double x[NVARS+1];
  double part1, part2;
  

  int ends = my_start + row_count;
  
  for ( member = my_start; member < ends; member++ )
  {
    for ( i = 0; i < NVARS; i++ )
    {
      x[i+1] = population[member].gene[i];
    }
    
    part1 = 0.0;
    for(i = 1; i < NVARS; i++)
    {
        part1 += pow(x[i], 2) / 4000;
    }
    
    part2 = 1.0;
    for(i = 1; i < NVARS; i++)
    {
        part2 *= cos( x[i] / sqrt(i) );
    }
       
    population[member].fitness = 1.0 + part1 - part2;
  }
  return;
}

//****************************************************************************80

int i4_uniform_ab ( int a, int b, int &seed )
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

void initialize ( string filename, int &seed )
{
  int i;
  ifstream input;
  int j;
  double lbound;
  double ubound;

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

    for ( j = 0; j < POPSIZE; j++ )
    {
      population[j].fitness = 0;
      population[j].rfitness = 0;
      population[j].cfitness = 0;
      population[j].lower[i] = lbound;
      population[j].upper[i]= ubound;
      population[j].gene[i] = r8_uniform_ab ( lbound, ubound, seed );
    }
  }

  input.close ( );

  return;
}
//****************************************************************************80

void keep_the_best ( )
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
  for ( i = 0; i < NVARS; i++ )
  {
    population[POPSIZE].gene[i] = population[cur_best].gene[i];
  }

  return;
}
//****************************************************************************80

void mutate ( int &seed )
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
        lbound = population[i].lower[j];
        ubound = population[i].upper[j];  
        population[i].gene[j] = r8_uniform_ab ( lbound, ubound, seed );
      }
    }
  }

  return;
}
//****************************************************************************80

double r8_uniform_ab ( double a, double b, int &seed )
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

void node_crossover(int my_start, int row_count,int &seed )
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

void node_mutate(int my_start, int row_count, int &seed )
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
        lbound = population[i].lower[j];
        ubound = population[i].upper[j];  
        population[i].gene[j] = r8_uniform_ab ( lbound, ubound, seed );
      }
    }
  }

  return;
}
//****************************************************************************80
void node_selector (int my_start, int row_count,int &seed )
{
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  int mem;
  double p;
  double sum;
  int ends = my_start + row_count;
//
//  Find the total fitness of the population.
//
  sum = 0.0;
  for ( mem = my_start; mem < ends; mem++ )
  {
    sum = sum + population[mem].fitness;
  }
//
//  Calculate the relative fitness of each member.
//
  for ( mem = my_start; mem < ends; mem++ )
  {
    population[mem].rfitness = population[mem].fitness / sum;
  }
// 
//  Calculate the cumulative fitness.
//
  population[my_start].cfitness = population[my_start].rfitness;
  for ( mem = my_start + 1; mem < ends; mem++ )
  {
    population[mem].cfitness = population[mem-1].cfitness +       
      population[mem].rfitness;
  }
// 
//  Select survivors using cumulative fitness. 
//
  for ( i = my_start; i < ends; i++ )
  { 
    p = r8_uniform_ab ( a, b, seed );
    if ( p < population[my_start].cfitness )
    {
      newpopulation[i] = population[my_start];      
    }
    else
    {
      for ( j = my_start; j < ends; j++ )
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
