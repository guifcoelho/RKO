#include "Method.h"
#include "Output.h"
#include "Data.h"

#include <sys/time.h>
#include <math.h>

// /************************************************************************************
// 								GLOBAL VARIABLES
// *************************************************************************************/

// Input File
char instance[100];                         // name of instance
int debug = 1;                              // 0 - run mode      		    1 - debug mode
int numDecoders = 1;                        // number of decoders
int numLS = 1;                              // 0 - without local search     > k - number of local search heuristics
int MAXTIME = 1;                            // maximum runtime
int MAXRUNS =  1;                           // maximum number of runs of the method
unsigned MAX_THREADS = 1;            		// number of threads
float OPTIMAL = 0;                          // optimal solution (if it is known)

int n;                                      // size of the solution

// Run
struct timeval Tstart, Tend, Tbest;         // computational time (unix systems)
char nameTable[100];                        // name of instance
TSol bestSolution;                          // best solution found in the A-BRKGA

// /************************************************************************************
// 								MAIN FUNCTION AREA
// *************************************************************************************/
int main(int argc, char *argv[ ])
{ 
    // read scenario name file
    char nameScenario[100];
    strncpy(nameScenario,argv[1],100);
    int method = atoi(argv[2]);  

    // file with test instances and input data
	FILE *arqProblems; 
    arqProblems = fopen (nameScenario, "r"); 

    if (arqProblems == NULL){
        printf("\nERROR: File %s not found\n", nameScenario);
        getchar();
        exit(1);
    }

    //read first line of arqProblems file
    fgets(nameTable, sizeof(nameTable), arqProblems); 

    // best solution that is saved in out file
    TSol sBest;

	// run the BRKGA-QL for all test instances
	while (!feof(arqProblems))
	{
		// read the name of instances, debug mode, local search module, maximum time, maximum number of runs, maximum number of threads
		fscanf(arqProblems,"%s %d %d %d %d %d %d %f", nameTable, &debug, &numDecoders, &numLS, &MAXTIME, &MAXRUNS, &MAX_THREADS, &OPTIMAL);
        strcpy(instance,nameTable);
        
        double foBest = INFINITY,
               foAverage = 0;

        float timeBest = 0,
              timeTotal = 0;

        std::vector <double> ofvs;
        ofvs.clear();

        // best solutions found in MAXRUNS
        sBest.ofv = INFINITY;

		// run BRKGA-QL MaxRuns for each instance
        printf("\n\nInstance: %s \nRun: ", instance);
        for (int j=0; j<MAXRUNS; j++)
        {
            // fixed seed
            if (debug == 1)
                srand(j+1); 
            else
                srand(time(NULL));

            printf("%d ", j+1);
            
            gettimeofday(&Tstart, NULL);
            gettimeofday(&Tend, NULL);
            gettimeofday(&Tbest, NULL);

            // best solution found in this run
            bestSolution.ofv = INFINITY;

            // execute the RKO method
            switch (method)
            {
            case 1:
                MultiStart();
                break;

            case 2: 
                SA(); 
                break;

            case 3:
                BRKGA();
                break;

            case 4:
                BRKGA_QL();
                break;
            
            default:
                break;
            } 

            gettimeofday(&Tend, NULL);

            // store the best solution found in MAXRUNS
            if (bestSolution.ofv < sBest.ofv)
                sBest = bestSolution;

            // calculate best and average results
            if (bestSolution.ofv < foBest)
                foBest = bestSolution.ofv;

            foAverage += bestSolution.ofv;

            // fitness of each solution found in the runs
            ofvs.push_back(bestSolution.ofv);

            timeBest += ((Tbest.tv_sec  - Tstart.tv_sec) * 1000000u + Tbest.tv_usec - Tstart.tv_usec) / 1.e6;
            timeTotal += ((Tend.tv_sec  - Tstart.tv_sec) * 1000000u + Tend.tv_usec - Tstart.tv_usec) / 1.e6; 
        }

        // create a .csv file with average results
        foAverage = foAverage / MAXRUNS;
        timeBest = timeBest / MAXRUNS;
        timeTotal = timeTotal / MAXRUNS;

        if (!debug)
        {
        	WriteSolution(sBest, n, timeBest, timeTotal, instance);
        	WriteResults(foBest, foAverage, ofvs, timeBest, timeTotal, instance);
        }
        else
        {
            WriteSolutionScreen(sBest, n, timeBest, timeTotal, instance);
        }
    }

    fclose(arqProblems);
    return 0;
}