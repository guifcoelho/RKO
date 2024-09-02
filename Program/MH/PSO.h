#ifndef _PSO_H
#define _PSO_H

/************************************************************************************
 Method: PSO
 Description: search process of the Particle Swarm Optimization
*************************************************************************************/
void PSO(int method, int control)
{
    int tamPop = 0;                          // number of particles
    float c1 = 0.0;
    float c2 = 0.0;
    float w = 0.0;

    std::vector <TSol> X;                    // current solutions
    std::vector <TSol> Pbest;                // best solutions 
    std::vector <std::vector <float> > V;    // particle velocity

    TSol Gbest;                              // global best solution
    double bestOFcurrent = 0;                // best ofv found in the current generation

    // local variables
    int numGenerations = 0;                  // number of generations
    float currentTime = 0;                   // computational time of the search process
    int bestGeneration = 0;                  // number of generation that found the best solution

    struct timespec TstartMH, TendMH;               // computational time (unix systems)
    clock_gettime(CLOCK_MONOTONIC, &TstartMH);
    clock_gettime(CLOCK_MONOTONIC, &TendMH);

    // offline control
    if (control == 0){
        // define parameters of PSO (offline)
        tamPop = 50;
        c1 = 2.05;
        c2 = 2.05;
        w  = 0.73;
    }

    // initialize population
    X.clear(); 
    Pbest.clear();
    V.clear(); 

    X.resize(tamPop);
    Pbest.resize(tamPop);
    V.resize(tamPop, std::vector<float>(n));

    // Create the initial chromosomes with random keys 
    Gbest.ofv = INFINITY;
    for (int i=0; i<tamPop; i++)
    {
        // initialize X[i]
        CreateInitialSolutions(X[i]); 

        // fitness of X[i]
        X[i].ofv = Decoder(X[i]);

        // initialize V[i][j]
        for (int j=0; j<n; j++)
            V[i][j] = randomico(0,1);

        // initialize Gbest
        if (X[i].ofv < Gbest.ofv)
            Gbest = X[i];

        // initialize Pbest
        Pbest[i] = X[i];    
    }
    
    // run the evolutionary process until stop criterion
    while (currentTime < MAXTIME)
    {
    	// number of generations
        numGenerations++;

        bestOFcurrent = INFINITY;
        double media = 0;
        for (int i=0; i<tamPop; i++)
        {
            if (stop_execution.load()) return;      

            double probUpdate = 1.0; // randomico(0.1, 0.3);

            // update particles X[i]
            for (int j=0; j<n; j++)
            {
                float r1 = randomico(0,1);
                float r2 = randomico(0,1);

                // update v[i][j]
                // V[i][j] = (w * V[i][j]) + (c1 * r1 * (Pbest[i].vec[j].rk - X[i].vec[j].rk)) + 
                //                           (c2 * r2 * (Gbest.vec[j].rk - X[i].vec[j].rk));

                // update v[i][j] com fator de constricao
                V[i][j] = w * (V[i][j] + (c1 * r1 * (Pbest[i].rk[j] - X[i].rk[j])) + 
                                         (c2 * r2 * (Gbest.rk[j] - X[i].rk[j])));
                
                if (randomico(0,1) < probUpdate){
                    // update X[i][j]
                    double oldrk = X[i].rk[j];
                    X[i].rk[j] = X[i].rk[j] + V[i][j];  

                    if (X[i].rk[j] < 0.0 || X[i].rk[j] >= 1.0) {
                        X[i].rk[j] = oldrk; //randomico(0.0, 1.0); //oldrk;
                        V[i][j] = 0;
                    }
                }
            }

            // fitness
            X[i].ofv = Decoder(X[i]);

            // // local search
            // if (randomico(0,1)<0.98)
            //     NelderMeadSearch(X[i]);
            // else
            //     RVND(X[i]);

            // set the best ofv found in this generation
            if (X[i].ofv < bestOFcurrent){
                bestOFcurrent = X[i].ofv;
            }

            // set Pbest
            if (X[i].ofv < Pbest[i].ofv){
                Pbest[i] = X[i];
            }

            // set Gbest 
            if (X[i].ofv < Gbest.ofv){
                Gbest = X[i];
                bestGeneration = numGenerations;    
            }

            media += X[i].ofv;
        } 

        // local search
        double oldGbest = Gbest.ofv;
        int escolhido = irandomico(0,tamPop-1);
        RVND(X[escolhido]);
        if (X[escolhido].ofv < Gbest.ofv){
            Gbest = X[escolhido];
            bestGeneration = numGenerations;    
        }

        // if (debug && method == 5) printf("\nGen: %d \t sBest: %lf \t media: %lf", numGenerations, Gbest.ofv, media/tamPop);

        if (bestGeneration == numGenerations || Gbest.ofv < oldGbest){
            // update the pool of solutions
            UpdatePoolSolutions(Gbest, method);
        }

        // RESET
        if ((numGenerations - bestGeneration) > 1000) {
            X.clear(); 
            X.resize(tamPop);

            //Gbest.ofv = INFINITY;
            for (int i=0; i<tamPop; i++){
                if (stop_execution.load()) return;      
                
                CreateInitialSolutions(X[i]);
                X[i].ofv = Decoder(X[i]);

                if (X[i].ofv < Gbest.ofv)
                    Gbest = X[i];

                Pbest[i] = X[i];
            }

            bestGeneration = numGenerations;
        }
        
        // terminate the evolutionary process in MAXTIME
        clock_gettime(CLOCK_MONOTONIC, &TendMH);
        currentTime = (TendMH.tv_sec - TstartMH.tv_sec) + (TendMH.tv_nsec - TstartMH.tv_nsec) / 1e9;
    }

    // free memory of PSO components
    X.clear();
    V.clear();
    Pbest.clear();
}

#endif