#ifndef _SA_H
#define _SA_H

/************************************************************************************
 Method: SA
 Description: search process of the Simulated Annealing
*************************************************************************************/
void SA(int method, int control)
{
    TSol s;                                  // current solution
    TSol sViz;                               // neighor solution
    TSol sBest;                              // best solution of SA

    double delta = 0;                        // difference between solutions
    double bestOFV;                          // value of the best solution in the current iteration

    double T0 = 0;                           // initial temperature
    double T = 0;                            // current temperature
    double alphaSA = 0;                      // cool rate
    int SAmax = 0;                           // number of iterations in a temperature T
    int IterT = 0;                           // iteracao corrente
    float betaMin = 0;                       // minimum perturbation
    float betaMax = 0;                       // maximum perturbation
    int reanneling = 0;                      // reanneling flag
    float currentTime = 0;                   // computational time of the search process

    struct timeval TstartMH, TendMH;         // computational time (unix systems)
    gettimeofday(&TstartMH, NULL);
    gettimeofday(&TendMH, NULL);

    // offline control
    if (control == 0){
        // define the parameters of the SA
        T0 = 1000000;
        SAmax = 100;
        alphaSA = 0.99;
        betaMin = 0.10;
        betaMax = 0.20;
    }

    // Create the initial solution with random keys
    CreateInitialSolutions(s); 
    s.ofv = Decoder(s);
    sBest = s;

    // run the search process until stop criterion
    while (currentTime < MAXTIME)
    {
        IterT = 0;
        if (!reanneling) T = T0;
        else T = T0*0.3;
        while (T > 0.0001 && currentTime < MAXTIME)
        {
            // offline
            if (control == 0)
            {
                if (T > 3000)
                {
                    SAmax = 300;
                    alphaSA = 0.99;
                    betaMin = 0.02;
                    betaMax = 0.05;
                }
                else if (T > 1000)
                {
                    SAmax = 200;
                    alphaSA = 0.98;
                    betaMin = 0.01;
                    betaMax = 0.03;
                }
                else 
                {
                    SAmax = 100;
                    alphaSA = 0.97;
                    betaMin = 0.01;
                    betaMax = 0.02;
                }
            }
            
            bestOFV = INFINITY;
            while (IterT < SAmax && currentTime < MAXTIME)
            {
                if (stop_execution.load()) return;      
                
                IterT++;

                // Shake the current solution
                sViz = s;
                ShakeSolution(sViz, betaMin, betaMax);

                // calculate the OFV
                sViz.ofv = Decoder(sViz);
                
                // value function is the best solution found in this iteration
                if (sViz.ofv < bestOFV)
                    bestOFV = sViz.ofv;
                
                // calculate the delta SA
                delta = sViz.ofv - s.ofv;

                // define from which solution to continue the search
                if (delta < 0)
                {
                    // update current solution
                    s = sViz;

                    // update the best solution found by SA
                    if (s.ofv < sBest.ofv)
                    {
                        sBest = s;

                        // update the pool of solutions
                        UpdatePoolSolutions(s, method);

                    }
                }
                else
                {
                    // metropolis criterion
                    double x = randomico(0,1);

                    if ( x < (exp(-delta/T)) )       
                        s = sViz;
                }
            } //End-SAmax

            // *************************************************************

            T = T * alphaSA;
            IterT = 0;

            // apply local search
            if (T < 1 && randomico(0,1) < 0.5)
            {
                RVND(s);
                // update the best solution found by SA
                if (s.ofv < sBest.ofv)
                {
                    sBest = s;

                    // update the pool of solutions
                    UpdatePoolSolutions(s, method);

                }
            }

            // if (debug && method == 1) printf("\nT: %lf \t current solution: %lf \t best solution: %lf", T, s.ofv, sBest.ofv);

            // terminate the search process in MAXTIME
            gettimeofday(&TendMH, NULL);
            currentTime = ((TendMH.tv_sec  - TstartMH.tv_sec) * 1000000u + TendMH.tv_usec - TstartMH.tv_usec) / 1.e6; 
        } //Fim-T

        // reanneling
        reanneling = 1;
    }
}

#endif