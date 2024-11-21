#ifndef _ILS_H
#define _ILS_H

/************************************************************************************
 Method: ILS
 Description: search process of the Iterated Local Search
*************************************************************************************/
void ILS(int method, bool find_best_mh_params)
{
    int Iter = 0;                                   // count the number of iterations of the ILS
    int IterImprov = 0;                             // store the last iteration that improve the current solution

    TSol sBest,                                     // local optimal solution (current)
         sLine,                                     // neighborhood solution
         sBestLine;                                 // local optimal of the neighborhood solution

    double betaMin = 0.0;                           // minimum perturbation
    double betaMax = 0.0;                           // maximum perturbation

    float currentTime = 0;                          // computational time of the search process
    int improv = 0;                                 // improvement flag

    double start_timeMH = get_time_in_seconds();    // start computational time
    double end_timeMH = get_time_in_seconds();      // end computational time

    // Q-Learning parameters
    std::vector<TState> S;                      // finite state space
    int numPar = 0;                             // number of parameters
    int numStates = 0;                          // number of states
    int iCurr = 0;                              // current (initial) state
    double epsilon=0;                           // greed choice possibility
    double lf=0;                                // learning factor
    double df=0;                                // discount factor
    double R=0;                                 // reward
    std::vector <std::vector <TQ> > Q;          // Q-Table
    std::vector<int> ai;                        // actions
    float epsilon_max = 1.0;                    // maximum epsilon
    float epsilon_min = 0.1;                    // minimum epsilon
    int Ti = 1;                                 // number of epochs performed
    int restartEpsilon = 1;                     // number of restart epsilon
    int st = 0;                                 // current state
    int at = 0;                                 // current action

    // ** read file with parameter values
    numPar = 2;
    std::vector<std::vector<double>> parameters;
    parameters.resize(numPar);

    readParameters(method, parameters, numPar);

    // offline control
    if (!find_best_mh_params){
        betaMin = parameters[0][0];
        betaMax = parameters[1][0];
    }

    // online control
    else{
        // Q-Learning
        // create possible states of the Markov chain
        CreateStates(parameters, method, numStates, numPar, S);

        // number of restart epsilon
        restartEpsilon = 1;

        // maximum epsilon
        epsilon_max = 1.0;

        // current state
        iCurr = irandomico(0,numStates-1);

        // define parameters of ILS
        betaMin = S[iCurr].par[0];
        betaMax = S[iCurr].par[1];
    }

    // number of iterations
    Iter = 0;

    // create initial solution
    CreateInitialSolutions(sBest);
    sBest.ofv = decoder(sBest.rk);

    // apply local search
    RVND(sBest);

    // terminate the search process in MAXTIME
    end_timeMH = get_time_in_seconds();
    currentTime = end_timeMH - start_timeMH;

    // run the search process until stop criterion
    while (currentTime < MAXTIME)
    {
        if (stop_execution.load()) return;

        // increase the number of ILS iterations
        Iter++;

        // define the parameters considering the current state and evolve a new iteration of the ILS
        // Q-Learning
        if (find_best_mh_params){
            // set Q-Learning parameters
            SetQLParameter(currentTime, Ti, restartEpsilon, epsilon_max, epsilon_min, epsilon, lf, df);

            // choose a action at for current state st
            at = ChooseAction(S, st, epsilon);

            // execute action at of st
            iCurr = S[st].Ai[at];

            // define the parameters according of the current state
            betaMin = S[iCurr].par[0];
            betaMax = S[iCurr].par[1];
        }

        // new iteration of the ILS
        sLine = sBest;

        // Shake the current solution
        ShakeSolution(sLine, betaMin, betaMax);

        // calculate the OFV
        sLine.ofv = decoder(sLine.rk);

        //s*' <- local search (s')
        sBestLine = sLine;
        RVND(sBestLine);

        //s* <- acceptance criterion (s*, s*', historico)
        if (sBestLine.ofv < sBest.ofv)
        {
            sBest = sBestLine;
            IterImprov = Iter;
            improv = 1;
            // update the pool of solutions
            UpdatePoolSolutions(sBest, method);
        }
        // else
        // {
        //     // metropolis criterion
        //     if (randomico(0,1) < (exp(-(sBestLine.ofv - sBest.ofv)/(1000 - 1000*(currentTime / MAXTIME)))) )
        //     {
        //         sBest = sBestLine;
        //     }
        // }

        // history
        if (Iter - IterImprov > 1000)   //1000 iteracoes sem melhora
        {
            // restart the search process with a new random solution
            CreateInitialSolutions(sBest);
            sBest.ofv = decoder(sBest.rk);
            IterImprov = Iter;
        }

        // if (debug && method == 3) printf("\nIter: %d [%.3lf %.3lf] \t s'Best: %lf \t sBest: %lf", Iter, betaMin, betaMax, sBestLine.ofv, sBest.ofv);

        // Q-Learning 1
        if (find_best_mh_params){
            // The reward function is based on improvement of the current best fitness and binary reward
            if (improv){
                R = 1;
                improv = 0;
            }
            else{
                R = (sBest.ofv - sBestLine.ofv)/sBestLine.ofv;
            }

            // if (debug) printf("\t [%.4lf, %d, %d] ", R, st, at);

            // index of the next state
            int st_1 = S[st].Ai[at];

            // Update the Q-Table value
            S[st].Qa[at] = S[st].Qa[at] + lf*(R + df*S[st_1].maxQ - S[st].Qa[at]);

            if (S[st].Qa[at] > S[st].maxQ)
            {
                S[st].maxQ = S[st].Qa[at];
                S[st].maxA = at;
            }

            // Define the new current state st
            st = st_1;
        }

        // terminate the search process in MAXTIME
        end_timeMH = get_time_in_seconds();
        currentTime = end_timeMH - start_timeMH;
    }

    // print policy
    // if (debug and control == 1)
    //     PrintPolicy(S, st);
}

#endif