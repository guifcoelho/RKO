#if defined(_WIN32) || defined(_WIN64)  // Windows
    #ifndef NOMINMAX
        # define NOMINMAX
    #endif
    #include <windows.h>
#else  // Unix-like (Linux, macOS)
    #include <time.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <numeric>
#include <map>
#include <limits>
#include <random>
#include <chrono>
#include <iomanip>
// #include <sstream>
// #include <fstream>
#include <omp.h>
#include <atomic>
#include <functional>
#include <memory>

#include "Main/Data.h"
#include "Main/GlobalVariables.h"
#include "Main/GlobalFunctions.h"
// #include "Main/Output.h"
// #include "Problem/Problem.h"
#include "MH/Method.h"
#include "MH/QLearning.h"
#include "MH/BRKGA.h"
#include "MH/SA.h"
#include "MH/GRASP.h"
#include "MH/ILS.h"
#include "MH/VNS.h"
#include "MH/PSO.h"
#include "MH/GA.h"
#include "MH/LNS.h"
#include "MH/BRKGA_CS.h"
#include "MH/MultiStart.h"


std::vector<double> solve(char nameScenario[256],
                          int method,
                          bool find_best_mh_params,
                          char _mh_params_file_path[255],
                          bool _debug,
                          int size,
                          std::function<double(std::vector<double>)> _decoder){

    decoder = _decoder;
    debug = _debug;
    strcpy(mh_params_file_path, _mh_params_file_path);

    // best solution that is saved in out file
    TSol sBest;
    sBest.flag = 0;
    sBest.label = 0;
    sBest.ofv = 0;
    sBest.similar = 0;
    sBest.promising = 0;
    sBest.rk.clear();

    // current random seed
    int RSEED = 0;

    // obatin a seed of the clock
    RSEED = std::chrono::steady_clock::now().time_since_epoch().count();

    // set new seed
    if (!debug) rng.seed (RSEED);

    // use a fixed seed in debug mode
    if (debug) rng.seed (1234);

    n = size;

    // define the rk order
    RKorder.resize(n);
    for (int i = 0; i < n; i++){
        RKorder[i] = i;
    }

    // define the total number of metaheuristics available
    #define NUM_MH 9

    // array of pointers to metaheuristic functions
    void (*functions_MH[NUM_MH])(int, bool) = {
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

    if (method <= NUM_MH) strcpy(nameMH, algorithms[method]);

    // create initial solutions in the pool of solutions
    CreatePoolSolutions();

    // best solution found in this run
    bestSolution = pool[0];

    // run all metaheuristics in parallel using OpenMP
    stop_execution.store(false);
    if (method == 9){
        omp_set_num_threads(NUM_MH);
        #pragma omp parallel private(RKorder, rng) shared(bestSolution, best_time, stop_execution)
        {
            #pragma omp for
            for (int i = 0; i < NUM_MH; ++i) {
                // checks the cancellation point
                #pragma omp cancellation point for

                // get the associated function
                void (*function_mh)(int, bool) = functions_MH[i];

                // int thread_id = omp_get_thread_num();
                if (debug) printf("\nThread %d executing MH_%d %s.", omp_get_thread_num(), i, algorithms[i]);

                // calls the metaheuristic function
                function_mh(i, find_best_mh_params);

                // cancels when a thread ends
                stop_execution.store(true);
                #pragma omp cancel for
            }
        }
        stop_execution.store(false);
    }
    // run a specific metaheuristic method
    else if (method < 9){
        // get the associated function
        void (*function_mh)(int, bool) = functions_MH[method];

        // calls the metaheuristic [method] function
        function_mh(method, find_best_mh_params);
    }
    // run a multi-start method
    else{
        strcpy(nameMH, algorithms[method]);
        MultiStart(method);
    }

    std::cout<<""<<std::endl;

    //When setting the decoder with a Python function: The decoder must be cleared in order to release GIL
    decoder = nullptr;

    return pool[0].rk;
}