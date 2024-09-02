#include <sys/time.h>
#include <math.h>
#include <cstring>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string.h>
#include <algorithm>
#include <utility>  // pair
#include <numeric>  // iota
#include <map>
#include <limits>
#include <random>
#include <chrono>
#include <iomanip> //graph
#include <sstream> //graph
#include <fstream> //graph
#include <omp.h>
#include <atomic>

#include "Data.h" 
#include "GlobalVariables.h" 
#include "GlobalFunctions.h"
#include "Output.h"
#include "../Problem/Problem.h"
#include "../MH/IPR.h"
#include "../MH/Method.h"
#include "../MH/BRKGA.h"
#include "../MH/SA.h"
#include "../MH/GRASP.h"
#include "../MH/ILS.h"
#include "../MH/VNS.h"
#include "../MH/PSO.h"
#include "../MH/GA.h"
#include "../MH/BRKGA_CS.h"
#include "../MH/LNS.h"
#include "../MH/MultiStart.h"

// /************************************************************************************
// 								MAIN FUNCTION AREA
// *************************************************************************************/
int main(int argc, char *argv[ ])
{ 
    // read scenario name file 
    char nameScenario[256];
    strncpy(nameScenario,argv[1],255);
    int method = atoi(argv[2]);  
    int control = atoi(argv[3]);  
    debug = atoi(argv[4]);

    // file with test instances and input data
	FILE *arqProblems; 
    arqProblems = fopen (nameScenario, "r"); 

    if (arqProblems == NULL){
        printf("\nERROR: File %s not found\n", nameScenario);
        exit(1);
    }

    // read first line of testScenario file
    if (fgets(nameTable, sizeof(nameTable), arqProblems) == NULL) {
        printf("\nERROR: File %s not found\n", nameTable);
        exit(1);
    }
    
    // best solution that is saved in out file
    TSol sBest;
    sBest.flag = 0;
    sBest.label = 0;
    sBest.ofv = 0;
    sBest.similar = 0;
    sBest.promising = 0;
    sBest.rk.clear();

	// run the RKO for all test instances
	while (!feof(arqProblems))
	{
		// read the name of instances, debug mode, local search module, maximum time, maximum number of runs, maximum number of threads
        if (fscanf(arqProblems,"%s %d %d", nameTable, &MAXTIME, &MAXRUNS) == 0) {
            printf("\nERROR: File %s not found\n", nameTable);
            exit(1);
        }
        strcpy(instance,nameTable);
        
        double foBest = INFINITY,
               foAverage = 0;

        float timeBest = 0,
              timeTotal = 0;

        std::vector <double> ofvs;
        ofvs.clear();

        // best solutions found in MAXRUNS
        sBest.ofv = INFINITY;

		// run RKO MaxRuns for each instance
        printf("\n\nInstance: %s \nRun: ", instance);
        for (run=0; run<MAXRUNS; run++)
        {
            // current random seed
            int RSEED = 0;        

            // obtain a seed from Taillard seeds
            // int i1 = run+1;
            // i1 = i1 % 400;
            // FILE *arq = fopen ("../Instances/Seeds_Taillard.txt","r");
            // if ( arq != NULL ) for (int i = 0; i <= i1; ++i) {
            //     if (fscanf(arq, "%d", &RSEED) != 1) {
            //         printf("Error reading data from seed file.\n");
            //     }
            // }
            // else printf ( "ERROR: Unable to open file for reading.\n" );
            // fclose (arq);

            // obatin a seed of the clock
            RSEED = std::chrono::steady_clock::now().time_since_epoch().count();

            // set new seed
            if (!debug) rng.seed (RSEED);
    
            // use a fixed seed in debug mode
            if (debug) rng.seed (1234);

            // runs
            printf("%d ", run+1);
            
            // computational times
            clock_gettime(CLOCK_MONOTONIC, &Tstart);
            clock_gettime(CLOCK_MONOTONIC, &Tend);
            clock_gettime(CLOCK_MONOTONIC, &Tbest);

            // free memory with problem data
            FreeMemoryProblem();

            // read data of the instance
            ReadData(nameTable);

            // define the rk order
            RKorder.resize(n);
            for (int i = 0; i < n; i++){
                RKorder[i] = i;
            }

            // ****************** execute the RKO method ******************
            
            // define the total number of metaheuristics available
            #define NUM_MH 9        

            // array of pointers to metaheuristic functions
            void (*functions_MH[NUM_MH])(int, int) = {
                BRKGA,      // 0
                SA,         // 1
                GRASP,      // 2
                ILS,        // 3
                VNS,        // 4
                PSO,        // 5
                GA,         // 6
                LNS,        // 7
                BRKGA_CS    // 8
            };

            // declare and initialize the array of strings
            const char *algorithms[] = {
                "BRKGA",     // 0
                "SA",        // 1
                "GRASP",     // 2
                "ILS",       // 3
                "VNS",       // 4
                "PSO",       // 5
                "GA",        // 6
                "LNS",       // 7
                "BRKGA_CS",  // 8
                "RKO",       // 9
                "MultiStart" // 10
            };

            if (method <= NUM_MH) strcpy(nameMH,algorithms[method]);

            // run all metaheuristics in parallel using OpenMP
            stop_execution.store(false);
            if (method == 9){
                // create initial solutions in the pool of solutions
                CretePoolSolutions();

                // best solution found in this run
                bestSolution = pool[0];

                omp_set_num_threads(NUM_MH);
                #pragma omp parallel private(RKorder, rng) shared(bestSolution, Tbest, stop_execution)
                {
                    #pragma omp for 
                    for (int i = 0; i < NUM_MH; ++i) { 
                        // checks the cancellation point
                        #pragma omp cancellation point for

                        // get the associated function
                        void (*function_mh)(int, int) = functions_MH[i];

                        // int thread_id = omp_get_thread_num();
                        if (debug) printf("\nThread %d executing MH_%d %s.", omp_get_thread_num(), i, algorithms[i]);
                        
                        // calls the metaheuristic function
                        function_mh(i, control);

                        // cancels when a thread ends
                        stop_execution.store(true);
                        #pragma omp cancel for
                    }
                }
                stop_execution.store(false);
            }
            // run a specific metaheuristic method
            else if (method < 9){
                // create initial solutions in the pool of solutions
                CretePoolSolutions();

                // best solution found in this run
                bestSolution = pool[0];

                // get the associated function
                void (*function_mh)(int, int) = functions_MH[method];

                // calls the metaheuristic [method] function
                function_mh(method, control);
            }
            // run a multi-start method
            else{
                strcpy(nameMH,algorithms[method]);
                MultiStart(method);
            }

            // run the path-relinking
            if (0){  
                for (int i=0; i<(int)pool.size(); i++)
                {
                    TSol solPR = pool[0];
                    TSol bestPR = solPR;

                    int k1, k2;
                    do {
                        k1 = irandomico(0,pool.size()-1);
                        k2 = irandomico(0,pool.size()-1);
                    }
                    while (k1 == k2);

                    solPR = IPR(pool[k1], pool[k2], 1.0, 1);
                }
            }

            clock_gettime(CLOCK_MONOTONIC, &Tend);

            // store the best solution found in MAXRUNS
            if (bestSolution.ofv < sBest.ofv)
                sBest = bestSolution;

            // calculate best and average results
            if (bestSolution.ofv < foBest)
                foBest = bestSolution.ofv;

            foAverage += bestSolution.ofv;

            // fitness of each solution found in the runs
            ofvs.push_back(bestSolution.ofv);

            timeBest += (Tbest.tv_sec - Tstart.tv_sec) + (Tbest.tv_nsec - Tstart.tv_nsec) / 1e9;
            timeTotal += (Tend.tv_sec - Tstart.tv_sec) + (Tend.tv_nsec - Tstart.tv_nsec) / 1e9; 
        }

        // create a .csv file with average results
        foAverage = foAverage / MAXRUNS;
        timeBest = timeBest / MAXRUNS;
        timeTotal = timeTotal / MAXRUNS;

        if (!debug)
        {
        	WriteSolution(nameMH, sBest, n, timeBest, timeTotal, instance);
        	WriteResults(nameMH, foBest, foAverage, ofvs, timeBest, timeTotal, instance);
        }
        else
        {
            WriteSolutionScreen(nameMH, sBest, n, timeBest, timeTotal, instance);
            // sBest.ofv = CalculateFitness(sBest,1);
        }

        // free memory with problem data
        FreeMemoryProblem();
    }

    fclose(arqProblems);
    return 0;
}
