#ifndef _ILS_H
#define _ILS_H

/************************************************************************************
 Method: ILS
 Description: search process of the Iterated Local Search
*************************************************************************************/
void ILS(int method, int control)
{
    int Iter = 0;                                   // count the number of iterations of the ILS
    int IterImprov = 0;                             // store the last iteration that improve the current solution

    TSol sBest,                                     // local optimal solution (current)
         sLine,                                     // neighborhood solution
         sBestLine;                                 // local optimal of the neighborhood solution

    double betaMin = 0.0;                           // minimum perturbation
    double betaMax = 0.0;                           // maximum perturbation

    float currentTime = 0;                          // computational time of the search process  

    struct timespec TstartMH, TendMH;               // computational time (unix systems)
    clock_gettime(CLOCK_MONOTONIC, &TstartMH);
    clock_gettime(CLOCK_MONOTONIC, &TendMH);

    // offline control
    if (control == 0){
        betaMin = 0.05;
        betaMax = 0.20;
    }

    // number of iterations
    Iter = 0;

    // create initial solution
    CreateInitialSolutions(sBest); 
    sBest.ofv = Decoder(sBest);

    // apply local search
    RVND(sBest);

    // run the search process until stop criterion
    while (currentTime < MAXTIME)
    {
        if (stop_execution.load()) return;      
        
        // increase the number of ILS iterations
        Iter++;
    
        // define the parameters considering the current state and evolve a new iteration of the ILS
        // offlline
        if (control == 0){
            betaMin = randomico(0.01, 0.05);
            betaMax = randomico(0.10, 0.20);
            // betaMin = 0.05;
            // betaMax = 0.20;
        }
    
        // new iteration of the ILS
        sLine = sBest;

        // Shake the current solution
        ShakeSolution(sLine, betaMin, betaMax);

        // calculate the OFV
        sLine.ofv = Decoder(sLine);

        //s*' <- local search (s')
        sBestLine = sLine;
        RVND(sBestLine);

        //s* <- acceptance criterion (s*, s*', historico)
        if (sBestLine.ofv < sBest.ofv)
        {
            sBest = sBestLine;
            IterImprov = Iter;
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

        //aplicar history
        if (Iter - IterImprov > 1000)   //1000 iteracoes sem melhora
        {
            //reiniciar a busca de uma solucao aleatoria
            CreateInitialSolutions(sBest); 
            sBest.ofv = Decoder(sBest);
            IterImprov = Iter;
        }

        // if (debug && method == 3) printf("\nIter: %d \t s'Best: %lf \t sBest: %lf", Iter, sBestLine.ofv, sBest.ofv);
        
        // terminate the evolutionary process in MAXTIME
        clock_gettime(CLOCK_MONOTONIC, &TendMH);
        currentTime = (TendMH.tv_sec - TstartMH.tv_sec) + (TendMH.tv_nsec - TstartMH.tv_nsec) / 1e9;
    }
}

#endif