
# RKO - Random-key optimizer

This is an implementation of the RKO to solve combinatorial optimization problems. The code was prepared for Unix and Windows systems. However, the openMP paradigm needs to be enabled.

This algorithm's C++ code has been designed to be easy to reuse. Users can only implement specific functions (read and decoder).

Here, we have the RKO version 1.0 code.


## References

When using this algorithm in academic studies, please refer to the following works:

RKO version 1.0

[1] Chaves, A.A., Resende, M.G.C., Schuetz, M.J.A.,  Brubaker, J.K., Katzgraber, H.G., Arruda, E.F., Silva, R.M.A.
A Random-Key Optimizer for Combinatorial Optimization. This paper wiil be submitted to the Journal of Heuristics.

Available here in technical report form:
https://doi.org/10.48550/arXiv.2411.04293

## Scope

This code has been designed to solve the Traveling Salesman Problem (TSP). To solve other problems, users only need to configure the Problem.h file.


## Running the algorithm

* Enter the Program directory: `cd Program`
* Run the make command: `make rebuild`
* Run the RKO: `./runTest ../Instances/testScenario.csv K T D`, where k is the number of metaheuristics used to solve the problem, T is the tuning method (0 is offline and 1 is online), and D is the debug mode (1 to print information in the screen and 0 to save information in file)
* In Windows: runTest.exe ../Instances/testScenario.csv K T D
* Available metaheuristics:

"BRKGA"      = 0,

"SA"         = 1,

"GRASP"      = 2,

"ILS"        = 3,

"VNS"        = 4,

"PSO"        = 5,

"GA"         = 6,

"LNS"        = 7,

"BRKGA_CS"   = 8,

"RKO"        = 9,

"MultiStart" = 10.

* Or compile via terminal: `g++ -std=c++20 -o runTest main.cpp -O3 -fopenmp`


## Code structure

The code structure is documented in [1] and organized in the following manner:

* **SPECIFIC_CODE:**
    * **Problem.cpp**: Contains data structure of the problem, the input function, and the decoder.

* **GENERAL_CODE:**
    * **/MH**: Contains all of the metaheuristic (MH) algorithm's mechanisms.
    * **/Main**: Contains the main function to start the algorithm and stores the global variables.
    * **Data.h**: Represents the data structures.
    * **Output.h**: Stores the output functions, including the best solution found and statistical analysis of the MH.

## File testScenario.csv is the input data problem, and each line consists of:

- Instance Name

- Maximum runtime (in seconds)

- Maximum number of runs

Users need to create a folder named "Instances/ProblemName", where the instances must be; users also need to create a folder named "Results", where the results files are written.

## File Parameters:

Users can choose two parameter settings: parameter tuning (option 0) and parameter control (option 1).
 - For parameter tuning, users must inform the static configuration of each parameter in the file ParametersOffline.txt.
 - For parameter control, users must inform a set of possible values for each parameter in the file ParametersOnline.txt. We use the Q-Learning method to learn the best configuration for each metaheuristic during the search process.

## Python

To run from python, first install the RKOpy api via pip from the root directory:

```
$ pip install .
```

Then define the decoder function and call the RKOpy api.

```
import rkopy

def decoder(keys: list[float]):
    return sum(keys)

best_sol = rkopy.solve(
    "Example",
    1,                          # Will run with Simulated Annealing
    False,                      # Whether or not it should use the Q-Learning algorithm to find the best MH params
    "ParametersOffline.txt",    # Path to the MH params
    True,                       # Whether or not it should run o debug mode
    10,                         # Size of the random keys vector
    decoder                     # The decoder function
)
```
