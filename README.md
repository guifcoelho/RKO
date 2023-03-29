
# RKO - Random-key optimitization 

This is an implementation of the RKO to solve combinatorial optmization problems.

The C++ code of this algorithm has been designed to be easy of reuse. Users can only implement specific functions (read, decoder and local search). 

Here we have the RKO version 1.0 code.


## References

When using this algorithm in academic studies, please refer to the following works:

RKO version 1.0
[1] Chaves, A.A. (2023)
RKO
Available here in technical report form (as soon).

## Scope

This code has been designed to solve the Traveling Salesman Problem (TSP). Users need to configure only Problem.cpp file to solve other problem.


## Running the algorithm

* Enter the Program directory: `cd Program`
* Run the make command: `make`
* Run the RKO: `./runTest testScenario.csv K`, k is the number of metatheuristic used to solve the TSP

* Or compile via terminal: `g++ -std=c++11 -o runTest BRKGA_QL.cpp Problem.cpp -O3 -openmp` (or -fopenmp if Linux)


## Code structure

The code structure is documented in [1] and organized in the following manner:

* **SPECIFIC_CODE:**
    * **Problem.cpp**: Contains data structure of the problem, the input function, the decoders, and the local search heuristics.

* **GENERAL_CODE:**
    * **MH.cpp**: Contains all of the metaheuristic (MH) algorithm's mechanisms.
    * **Main.cpp: Contains the main function to start the algorithm and stores the global variables.
    * **Data.h**: Represents the data structures.
    * **Output.h**: Stores the outputs functions, including the best solution found and statistical analysis of the MH.

## File testScenario.csv is the input data problem and each line consists in:

- Instance Name
- Run mode (0 = debug, prints in the screen; 1 = run, prints in files)
- Number of implemented decoders
- Number of implemented local search heuristics (0 if local search is not available)
- Maximum rumtime (in seconds)
- Maximum number of runs
- Number of threads used in OpenMP
- Optimal or lower bound solution (if is known), 0 otherwise

Users need to create a folder named "Instances/ProblemName", where the instances must be; Users also need to create a folder named "Results" where the results files are writed.
