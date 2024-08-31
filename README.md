
# RKO - Random-key optimizer 

This is an implementation of the RKO to solve combinatorial optimization problems.

This algorithm's C++ code has been designed to be easy to reuse. Users can only implement specific functions (read and decoder). 

Here, we have the RKO version 1.0 code.


## References

When using this algorithm in academic studies, please refer to the following works:

RKO version 1.0
[1] Chaves, A.A. (2023)
RKO
Available here in technical report form (as soon as possible).

## Scope

This code has been designed to solve the Traveling Salesman Problem (TSP). To solve other problems, users need to configure only the Problem.cpp file.


## Running the algorithm

* Enter the Program directory: `cd Program`
* Run the make command: `make`
* Run the RKO: `./runTest testScenario.csv K`, k is the number of metaheuristics used to solve the TSP

* Or compile via terminal: `g++ -std=c++11 -o runTest BRKGA_QL.cpp Problem.cpp -O3 -fopenmp`


## Code structure

The code structure is documented in [1] and organized in the following manner:

* **SPECIFIC_CODE:**
    * **Problem.cpp**: Contains data structure of the problem, the input function, the decoders, and the local search heuristics.

* **GENERAL_CODE:**
    * **MH.cpp**: Contains all of the metaheuristic (MH) algorithm's mechanisms.
    * **Main.cpp**: Contains the main function to start the algorithm and stores the global variables.
    * **Data.h**: Represents the data structures.
    * **Output.h**: Stores the output functions, including the best solution found and statistical analysis of the MH.

## File testScenario.csv is the input data problem, and each line consists of:

- Instance Name
- Maximum runtime (in seconds)
- Maximum number of runs

Users need to create a folder named "Instances/ProblemName", where the instances must be; users also need to create a folder named "Results", where the results files are written.
