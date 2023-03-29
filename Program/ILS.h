#ifndef _ILS_H
#define _ILS_H

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

extern int n;                               // size of cromossoms
extern TSol bestSolution;                   // best solution found in the A-BRKGA


//ILS
static double beta = 0;                            //taxa de perturbação
static int Iter = 0;
static int IterMax = 5000;
static int IterMelhora = 0;
static int delta = 0;

static TSol s,
            sMelhor,     
            sViz,
            sMelhorViz;



//****************************** General Functions **********************************

/************************************************************************************
 Method: updateBestSolution()
 Description: Update the best solution found during the run
*************************************************************************************/
static void updateBestSolution(TSol s);

/************************************************************************************
 Method: CREATE INITIAL SOLUTIONS
 Description: create a initial chromossom with random keys
*************************************************************************************/
static TSol CreateInitialSolutions();

/************************************************************************************
 Method: PEARSON CORRELATION
 Description: calculate the Pearson correlation coefficient between two chromossoms
*************************************************************************************/
// static double PearsonCorrelation(std::vector <TVecSol> s1, std::vector <TVecSol> s2);

/************************************************************************************
 Metodo: IC(TSol Pop)
 Description: apply clustering method to find promising solutions in the population
*************************************************************************************/
// static void IC();

/************************************************************************************
 Method: LP
 Description: Apply Label Propagation to find communities in the population
*************************************************************************************/
// static void LP(std::vector<std::vector<std::pair<int, double> > > listaArestas);

/************************************************************************************
 Method: PROMISINGLP
 Description: Find the promising solutions to represent the communities
*************************************************************************************/
// static void PromisingLP();

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