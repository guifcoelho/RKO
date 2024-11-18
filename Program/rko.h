// #include <sys/time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_WIN32) || defined(_WIN64)  // Windows
    #include <windows.h>
#else  // Unix-like (Linux, macOS)
    #include <time.h>
#endif

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
#include <sstream>
#include <fstream>
#include <omp.h>
#include <atomic>
#include <functional>

#include "Main/Data.h"
#include "Main/GlobalVariables.h"
#include "Main/GlobalFunctions.h"
// #include "Main/Output.h"
// #include "Problem/Problem.h"
#include "MH/Method.h"
#include "MH/QLearning.h"
// #include "MH/BRKGA.h"
#include "MH/SA.h"
// #include "MH/GRASP.h"
// #include "MH/ILS.h"
// #include "MH/VNS.h"
// #include "MH/PSO.h"
// #include "MH/GA.h"
// #include "MH/BRKGA_CS.h"
// #include "MH/LNS.h"
// #include "MH/MultiStart.h"



std::vector<double> solve(char nameScenario[256],
                           int method,
                           bool find_best_mh_params,
                           char mh_params_file_path[255],
                           bool _debug,
                           int size,
                           std::function<double(std::vector<double>)> decoder){
    _decoder = decoder;
    debug = _debug;
    strcpy(_mh_params_file_path, mh_params_file_path);

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

    // // define the total number of metaheuristics available
    // #define NUM_MH 9

    // // array of pointers to metaheuristic functions
    // void (*functions_MH[NUM_MH])(int, int) = {
    //     BRKGA,      // 0
    //     SA,         // 1
    //     GRASP,      // 2
    //     ILS,        // 3
    //     VNS,        // 4
    //     PSO,        // 5
    //     GA,         // 6
    //     LNS,        // 7
    //     BRKGA_CS    // 8
    // };

    // // declare and initialize the array of strings
    // const char *algorithms[] = {
    //     "BRKGA",     // 0
    //     "SA",        // 1
    //     "GRASP",     // 2
    //     "ILS",       // 3
    //     "VNS",       // 4
    //     "PSO",       // 5
    //     "GA",        // 6
    //     "LNS",       // 7
    //     "BRKGA_CS",  // 8
    //     // "RKO",       // 9
    //     // "MultiStart" // 10
    // };

    // if (method <= NUM_MH) strcpy(nameMH, algorithms[method]);

    // // get the associated function
    // void (*function_mh)(int, int) = functions_MH[method];

    // create initial solutions in the pool of solutions
    CreatePoolSolutions();

    // best solution found in this run
    bestSolution = pool[0];

    start_time = get_time_in_seconds();

    // calls the metaheuristic [method] function
    SA(method, find_best_mh_params);

    end_time = get_time_in_seconds();

    std::cout << std::endl;
    std::cout << "Terminou" << std::endl;

    return pool[0].rk;
}