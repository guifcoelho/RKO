#ifndef _VNS_H
#define _VNS_H

/************************************************************************************
 Method: VNS
 Description: search process of the Variable Neighborhood Search
*************************************************************************************/
void VNS(int method, int control)
{
    double beta = 0;                    // perturbation rate
    int Iter = 0;                       // current iteration
    int IterMelhora = 0;                // iteration that found the best solution
    int r = 0;                          // number of neighborhood strutures
    double betaMin = 0.0;               // intensity of perturbation

    TSol s,                             // current solution
         sBest,                         // best solution of VNS
         sLine,                         // neighborhood solution
         sBestLine;                     // best neighborhood solution

    
    float currentTime = 0;              // computational time of the search process

    struct timespec TstartMH, TendMH;   // computational time (unix systems)
    clock_gettime(CLOCK_MONOTONIC, &TstartMH);
    clock_gettime(CLOCK_MONOTONIC, &TendMH);

    // offline control
    if (control == 0){
        betaMin = 0.03;
        r = 6;
    }

    // Create the initial solution with random keys
    sBest.ofv = INFINITY;
    for (int i=0; i<1; i++)
    {
        CreateInitialSolutions(s); 
        s.ofv = Decoder(s);
        if (s.ofv < sBest.ofv)
            sBest = s;
    }

    // current solution
    s = sBest;

    // run the search process until stop criterion           
    while (currentTime < MAXTIME)
    {
        // current neighborhood
        int k = 1;
        while (k <= r && currentTime < MAXTIME)
        {
            if (stop_execution.load()) return;      
            
            Iter++;
            
            //s' <- perturb the best solution in the neighborhood k
            beta = randomico(k*betaMin,(k+1)*betaMin);

            // perturb the current solution (s)
            sLine = s;
            ShakeSolution(sLine, beta, beta);

            // calculate OFV
            sLine.ofv = Decoder(sLine);

            //s*' <- local search (s')
            sBestLine = sLine; 
            RVND(sBestLine);

            //s <- acceptance criterion (s,s*', historico)
            if (sBestLine.ofv < s.ofv)
            {
                s = sBestLine;
                
                // update the best solution found in VNS
                if (s.ofv < sBest.ofv){   
                    sBest = s;   

                    // return to the first neighborhood structure          
                    k = 1;
                    IterMelhora = Iter;

                    // update the pool of solutions
                    UpdatePoolSolutions(sBestLine, method);
                }
            }
            else
            {
                // next neighborhood structure
                k++; 
            }

            // if (debug && method == 4) printf("\nIter: %d \t s'Best: %lf \t sBest: %lf", Iter, sBestLine.ofv, sBest.ofv);

            //aplicar history
            if (Iter - IterMelhora > 1000)   
            {
                // restart the search process
                CreateInitialSolutions(s); 
                s.ofv = Decoder(s);
                IterMelhora = Iter;
            }

            // terminate the evolutionary process in MAXTIME
            clock_gettime(CLOCK_MONOTONIC, &TendMH);
            currentTime = (TendMH.tv_sec - TstartMH.tv_sec) + (TendMH.tv_nsec - TstartMH.tv_nsec) / 1e9;
        }
    }
}


#endif
