FYP - 2017
==========
    Author: Abdul Halim
    Email: 13029096@studentmail.ul.ie
    FYP Supervisor: J. J. Collins
    Dept. of CSIS
    University of Limerick


Description:
===========
This project was undertaken to parallelise a Simple Genetic Algorithm using
Message Passing Interface(MPI) on Beowulf Cluster to carry on some comparative 
analysis on performance gain by varying number of nodes in a cluster and varying
population size of the Genetic Algorithm.

For this project, A Beowulf cluster was setup using Ubuntu 16.04 server edition.

Selected sequential SGA was taken from a open source implementation which was 
released under GPL LGPL license.

This code is available at: 
https://people.sc.fsu.edu/~jburkardt/cpp_src/simple_ga/simple_ga.cpp

This document is intended to give a basic run down on this project's 
dependencies and compiling and running the test scripts.

DEPENDENCIES
============
Environment:
    OS: Ubuntu 16.04
    MPI Libraries: MPICH3.2 with Hydra process manager
    Compiler: GCC & mpicc
    Language: C/C++
    Scripting language: Python 2.7

Beowulf cluster: (Consult FYP report)
    *) Dedicated MPI user
    *) Requires NFS setup
    *) Password less SSH login

IMPLEMENTATION
==============
In this attempt of parallelising there were three different versions implemented
of which source codes are included here. Full description of these 
implementation are out of the scope of this README file and are discussed on
corresponding FYP report for this project.
Actual source codes contains code fragments that are commented out which shows 
the GA output. For this project running time analysis there are codes added to
get the system time which could be commented out if not required.


COMPILATION
===========
Compilation is done using GNU Make. There is Makefile name "makefile" in the
project directory. This make file lets creating different binaries using 
different parameter values defined at compile time.

Variables:
    POPSIZE = { population size }
    MAXGENS = { number of generation }
    PGA = { Parallel GA source file name }
    SGA = { SGA source file name }
    PGA_EXEC = { name of PGA output binary }
    SGA_EXEC = { name of SGA output binary }

Executable binary files are saved in bin/ directory.


To compile run following command in project directory:
    make

TEST SCRIPTS
============
Note: the parallel GA can only be run if a Beowulf Cluster is setup and 
configured as described in the FYP report.

There are number of python test scripts that were used gather running time data
of SGA and PGA.
    
    *) test_xNodes.py -  Get running time comparison data of SGA vs. PGA by
varying number of processors for the PGA.
    
    *) test_xPopulation.py - Get running time comparison data of SGA vs. PGA by
by varying population size.

    *) get20runs_data.py - Used to get program running time of all three PGA
implementations to get AVG running time.

All of these three scripts contains variables to define, target program binary,
output directory & file.  
    

    

