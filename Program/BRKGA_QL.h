#ifndef _BRKGA_QL_H
#define _BRKGA_QL_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <cstring>
#include <string.h>
#include <algorithm>
#include <sys/time.h>
#include <utility>  // pair
#include <numeric>  // iota
#include <map>
#include <limits>
#include <random>
#include <chrono>
#include <iomanip> //graph
#include <sstream> //graph
#include <fstream> //graph

#include "Data.h"
#include "Problem.h"
#include "Method.h"


// Variables declared in main.cpp
extern char nameTable[100];                 // name of instance
extern int debug;                           // 0 - run mode      		    1 - debug mode
extern int numDecoders;                     // number of decoders
extern int numLS;                           // 0 - without local search     > k - number of local search heuristics
extern int MAXTIME;                         // maximum runtime
extern float OPTIMAL;                       // optimal solution (if it is known)
extern struct timeval Tstart, Tend, Tbest;  // computational time (unix systems)  

extern int n;                                      // size of cromossoms
extern TSol bestSolution;                          // best solution found in the A-BRKGA


//BRKGA
static int p;          	                        // size of population
static double pe;              	                // fraction of population to be the elite-set
static double pm;          	                    // fraction of population to be replaced by mutants
static double rhoe;             	                // probability that offspring inherit an allele from elite parent

static double sigma;                               // pearson correlation factor

static const double PI = 3.14159265;               // pi

static std::vector <TSol> Pop;                     // current population
static std::vector <TSol> PopInter;               	// intermediary population

// Q-Learning parameters
static double epsilon;                             // greed choice possibility
static double lf;                                  // learning factor
static double df;                                  // discount factor
static double R;                                   // reward
static double qTotal;                              // q*

// list of actions RL
static int sizeP []    = {233, 377, 610, 987, 1597, 2584};
static double Pe[]     = {0.10, 0.15, 0.20, 0.25, 0.30}; 
static double Pm[]     = {0.01, 0.02, 0.03, 0.04, 0.05}; 
static double Rhoe[]   = {0.55, 0.60, 0.65, 0.70, 0.75, 0.80}; 
 

// number of parameters in Q-table
static const int par = 4;

// actions
static int a0 = 0;                                 // p (current action)
static int a1 = 0;                                 // pe (current action)
static int a2 = 0;                                 // pm (current action)
static int a3 = 0;                                 // rhoe (current action)

static float Qmax = 0;

static std::vector <std::vector <TQ> > Q;          // Q-Table

// QL parameters (auxiliar)
static float epsilon_max = 1.0;                    // maximum epsilon 
static float epsilon_min = 0.1;                    // minimum epsilon
static int Ti = 1;                                 // number of epochs performed
static int restartEpsilon = 1;                     // number of restart epsilon



//****************************** General Functions **********************************

// /************************************************************************************
//  Method: BRKGA()
//  Description: Apply the method BRKGA to solve the problem
// *************************************************************************************/
// void BRKGA_QL();

/************************************************************************************
 Method: updateBestSolution()
 Description: Update the best solution found during the run
*************************************************************************************/
static void updateBestSolution(TSol s);

/************************************************************************************
 Method: InitiateQTable()
 Description: Initiate the Q-Table with random values
*************************************************************************************/
static void InitiateQTable();

/************************************************************************************
 Method: ChooseAction()
 Description: Choose actions and update the parameters
*************************************************************************************/
static void ChooseAction();

/************************************************************************************
 Method: UpdatePopulationSize()
 Description: Update the population size with new value of p
*************************************************************************************/
static void UpdatePopulationSize();

/************************************************************************************
 Method: UpdateQLParameters(currentTime)
 Description: Update the parameters of the Q-Learning method
*************************************************************************************/
static void SetQLParameters(float currentTime);

/************************************************************************************
 Method: UpdateQTable()
 Description: Update the values of Q-Table
*************************************************************************************/
static void UpdateQTable();
 
/************************************************************************************
 Method: CREATE INITIAL SOLUTIONS
 Description: create a initial chromossom with random keys
*************************************************************************************/
static TSol CreateInitialSolutions();

/************************************************************************************
 Method: ChaoticInd
 Description: create a solution between a mutant and a elite individual
*************************************************************************************/
static void ChaoticInd(TSol &s);

/************************************************************************************
 Method: PARAMETRICUNIFORMCROSSOVER
 Description: create a new offspring with parametric uniform crossover
*************************************************************************************/
static TSol ParametricUniformCrossover(int elitesize, int popSize);

/************************************************************************************
 Method: PEARSON CORRELATION
 Description: calculate the Pearson correlation coefficient between two chromossoms
*************************************************************************************/
static double PearsonCorrelation(std::vector <TVecSol> s1, std::vector <TVecSol> s2);

/************************************************************************************
 Metodo: IC(TSol Pop)
 Description: apply clustering method to find promising solutions in the population
*************************************************************************************/
static void IC();

/************************************************************************************
 Method: LP
 Description: Apply Label Propagation to find communities in the population
*************************************************************************************/
static void LP(std::vector<std::vector<std::pair<int, double> > > listaArestas);

/************************************************************************************
 Method: PROMISINGLP
 Description: Find the promising solutions to represent the communities
*************************************************************************************/
static void PromisingLP();

/************************************************************************************
Method: FREE MEMORY
Description: free memory of global vector
*************************************************************************************/
static void FreeMemory();

/************************************************************************************
 Method: RANDOMICO
 Description: Generate a double random number between min and max
*************************************************************************************/
static double randomico(double min, double max);

/************************************************************************************
 Method: IRANDOMICO
 Description: Generate a int random number between min and max
*************************************************************************************/
static int irandomico(int min, int max);

#endif