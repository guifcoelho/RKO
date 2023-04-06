#ifndef _PROBLEM_H
#define _PROBLEM_H

#define INFINITO 999999999

#include <math.h>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "Data.h"

// Variables declared in main.cpp
extern int numDecoders;                     // number of decoders
extern int numLS;                           // 0 - without local search     > k - number of local search heuristics
extern int n;                               // size of cromossoms

//------ DEFINITION OF TYPES OF PROBLEM SPECIFIC --------

static std::vector <int> a;							// assigment vector
static std::vector <int> h;							// hub vector
static std::vector <int> hP;						// hub position vector 
static std::vector <std::vector <int> > Tp;			// matrix with tree
static std::vector <std::vector <double> > cT;		// matrix with cost of tree


//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  --------

static std::vector <std::vector <double> > dist;	// matrix with Euclidean distance
static std::vector <std::vector <double> > flow;	// matrix with flow

static int nPontos;
static int nHubs;
static double alfaH;


//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------


/************************************************************************************
 Method: ReadData
 Description: read input data of the problem
*************************************************************************************/
void ReadData(char nameTable[]);

/************************************************************************************
 Method: Decoder()
 Description: Convert a random key solution in a real problem solution
*************************************************************************************/
void Decoder(TSol &s);

/************************************************************************************
 Method: LocalSearch
 Description: RVND
*************************************************************************************/
void LocalSearch(TSol &s);


/************************************************************************************
 Method: CalculateFitness
 Description: calculate the fitness of a chromossom s
*************************************************************************************/
double CalculateFitness(TSol s);

/************************************************************************************
 Method: Dec1
 Description: sort decoder 
*************************************************************************************/
void Dec1(TSol &s);

/************************************************************************************
 Method: LS1
 Description: 2-Opt
*************************************************************************************/
void LS1(TSol &s);

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem();

#endif
