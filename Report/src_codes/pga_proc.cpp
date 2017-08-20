# include <mpi.h>
# include "params.hpp"
# include "sga_proc.hpp"
# include "pga_proc.hpp"

//*******************************************************
void slave_process(int &seed, MPI_Datatype& mpi_genotype,
        int nodes_data[][2])
{
    int my_start = nodes_data[myId][0];
    int row_count = nodes_data[myId][1];
    for(int generation = 0; generation < MAXGENS; generation++)
    {
        /* Get updated subset of population from master */
        ierr = MPI_Recv(&population[my_start], row_count, mpi_genotype, MASTER_ID, 
                SEND_DATA_TAG, MPI_COMM_WORLD, &status);
        /* Evaluate */
        evaluate (my_start, row_count);
        
         /* Send evaluated population to master */
        MPI_Send(&population[my_start], row_count, mpi_genotype, MASTER_ID,
                RETN_DATA_TAG, MPI_COMM_WORLD);
    }
    
}

void master_process(string& filename, int &seed, 
        MPI_Datatype &mpi_genotype, int nodes_data[][2], clock_t start)
{
   timestamp ( );
   cout << "\n";
   cout << "PARALLEL_GA:\n";
   cout << "  C++ version\n";
   cout << "  A simple example of a parellel genetic algorithm.\n";
   
   if ( NVARS < 2 )
   {
       cout << "\n";
       cout << "  The crossover modification will not be available,\n";
       cout << "  since it requires 2 <= NVARS.\n";
   }
   
    initialize ( filename, seed );
    
    evaluate (0, POPSIZE );
    
    keep_the_best ( );
       
    /* Run GA ops for given number of generation */
    for(int generation = 0; generation < MAXGENS; generation++)
    {   
        selector ( seed );
        crossover ( seed );
        mutate ( seed );
        /* Master generate report  */
       report(generation);

        /* Update all nodes with it's part of current population */
        int start_index, item_count;
        for(int nodeId = 1; nodeId < numNodes; nodeId++)
        {
            start_index = nodes_data[nodeId][0];
            item_count = nodes_data[nodeId][1];
            MPI_Send(&population[start_index], item_count, mpi_genotype, nodeId, SEND_DATA_TAG,
                    MPI_COMM_WORLD);
        }
        
        /* Evaluate master's parts of the population */
        evaluate(nodes_data[0][0], nodes_data[0][1]);
        
        /* Collect all parts of new population from each node */
        for(int nodeId = 1; nodeId < numNodes; nodeId++)
        {
            start_index = nodes_data[nodeId][0];
            item_count = nodes_data[nodeId][1];
            MPI_Recv(&population[start_index], item_count, mpi_genotype,
                    nodeId, RETN_DATA_TAG, MPI_COMM_WORLD, &status);
        }
        
        elitist();
    }
   
   cout << "\n";
   cout << "  Best member after " << MAXGENS << " generations:\n";
   cout << "\n";
   
   for ( int i = 0; i < NVARS; i++ )
   {
       cout << "  var(" << i << ") = " << population[POPSIZE].gene[i] << "\n";
   }
   
   cout << "\n";
   cout << "  Best fitness = " << population[POPSIZE].fitness << "\n";
   //
   //  Terminate.
   //
   cout << "\n";
   cout << "PARALLEL_GA:\n";
   cout << "  Normal end of execution.\n";
   cout << "\n";
   timestamp ( );
    
}
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

//*******************************************************

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

//*******************************************************

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

//*******************************************************

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
//*******************************************************

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