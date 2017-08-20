/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   parallel_ga.hpp
 * Author: Abdul Halim <13029096@studentmail.ul.ie>
 *
 * Created on 09 February 2017, 15:28
 */
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
using namespace std;

#ifndef PARALLEL_GA_HPP
#define PARALLEL_GA_HPP

void node_selector(int my_start, int row_count, int &seed );
void node_crossover(int my_start, int row_count, int &seed );
void node_mutate(int my_start, int row_count, int &seed );
void node_evaluate(int my_start, int row_count);
void node_elitist(int my_start, int row_count, int gen_best_worst[]);
void master_elitist(int best_worst[][2], int node_count);
void node_Xover ( int one, int two, int &seed );
void print_population();


/* Temporary defines [ *** REMOVE THEM LATER *** ] */
extern genotype population[];
extern genotype newpopulation[];
extern void crossover ( int &seed );
extern void elitist ( );
extern void evaluate ( );
extern int i4_uniform_ab ( int a, int b, int &seed );
extern void initialize ( string filename, int &seed );
extern void keep_the_best ( );
extern void mutate ( int &seed );
extern double r8_uniform_ab ( double a, double b, int &seed );
extern void report ( int generation );
extern void selector ( int &seed );
extern void timestamp ( );
extern void Xover ( int one, int two, int &seed );

#endif /* PARALLEL_GA_HPP */

