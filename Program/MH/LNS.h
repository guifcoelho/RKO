#ifndef _LNS_H
#define _LNS_H

/************************************************************************************
 Method: fareySequence
 Description: generate the farey sequence F of order num
*************************************************************************************/
void fareySequence(int num, std::vector<double> &F) {
    int a = 0, b = 1, c = 1, d = num;
    // printf("Farey Sequence of order %d:\n", num);
    // printf("%.4lf ", (double)a/b);
    F.push_back((double)a/b);

    while (c <= num) {
        int k = (num + b) / d;
        int temp_a = a, temp_b = b;
        a = c;
        b = d;
        c = k * c - temp_a;
        d = k * d - temp_b;
        if (a != b){
            // printf("%.4lf ", (double)a/b);
            F.push_back((double)a/b);
        }
        else{
            F.push_back(((double)a/b) - 0.00001);
        }   
    }
}

/************************************************************************************
 Method: LNS
 Description: search process of the Large Neighborhood Search
*************************************************************************************/
void LNS(int method, int control)
{
    double T0;                           // initial temperature
    double T;                            // current temperature
    double alphaLNS;                     // cool rate
    double betaMin = 0.0;                // minimum perturbation
    double betaMax = 0.0;                // maximum perturbation
    int reanneling = 0;                  // reanneling flag

    TSol sLine,                          // neighborhood solution
         sLineBest,                      // best neighborhood solution
         s,                              // current solution
         sBest;                          // best solution of LNS
        
    float currentTime = 0;               // computational time of the search process  

    struct timespec TstartMH, TendMH;    // computational time (unix systems)
    clock_gettime(CLOCK_MONOTONIC, &TstartMH);
    clock_gettime(CLOCK_MONOTONIC, &TendMH);

    // create a Farey sequence
    std::vector<double> F;
    fareySequence(7,F);

    // offline control
    if (control == 0){
        betaMin = 0.05;
        betaMax = 0.20;
        T0 = 10000;
        alphaLNS = 0.95;
    }

    // Create the initial solution with random keys
    CreateInitialSolutions(s); 
    s.ofv = Decoder(s);    
    sBest = s;

    // run the search process until stop criterion
    while (currentTime < MAXTIME)
    {
        if (!reanneling) T = T0;
        else T = T0*0.3;
        while (T > 0.0001 && currentTime < MAXTIME)
        {
            // Destroy the current solution s
            sLine = s;
            int intensity = irandomico(betaMin*n, betaMax*n);
            
            // define which rk will be deleted - Random Removal
            std::shuffle (RKorder.begin(), RKorder.end(),rng);

            // Repair the current solution
            sLineBest = sLine;
            sLineBest.ofv = INFINITY;

            // generate a new rk for each deleted position
            for (int k=0; k<intensity; k++)
            {
                int pos = RKorder[k];
                double OFVbest = INFINITY;
                double rkBest = 0;

                for (int j=0; j<(int)F.size()-1; j++)
                {
                    if (stop_execution.load()) return;      
                    
                    // generate a random value between two intervals of the Farey sequence
                    sLine.rk[pos] = randomico(F[j], F[j+1]);

                    sLine.ofv = Decoder(sLine);

                    if (sLine.ofv < OFVbest)
                    {
                        OFVbest = sLine.ofv;
                        rkBest = sLine.rk[pos];
                    }

                    // if (sLine.ofv < sLineBest.ofv)
                    //     sLineBest = sLine;
                }

                // continuar a busca a partir do melhor valor de random key encontrado para a rk pos
                sLine.ofv = OFVbest;
                sLine.rk[pos] = rkBest;
            }

            // local search
            sLineBest = sLine;
            if (randomico(0,1) < 0.2) 
                RVND(sLineBest);
            else
                NelderMeadSearch(sLineBest);
            
            // calculate delta
            double delta = sLineBest.ofv - s.ofv;

            // acceptance criterion
            if (delta <= 0)
            {
                s = sLineBest;

                if (s.ofv < sBest.ofv){
                    sBest = s;

                    // update the pool of solutions
                    UpdatePoolSolutions(s, method);
                }
            }
            else
            {
                double x = randomico(0,1);

                if ( x < (exp(-delta/T)) )       
                    s = sLineBest;
            }

            // if (debug && method == 7) printf("\nT: %lf \t s current: %lf \t sBest: %lf", T, s.ofv, sBest.ofv);
        
            // *************************************************************

            T = T * alphaLNS;

            // terminate the evolutionary process in MAXTIME
            clock_gettime(CLOCK_MONOTONIC, &TendMH);
            currentTime = (TendMH.tv_sec - TstartMH.tv_sec) + (TendMH.tv_nsec - TstartMH.tv_nsec) / 1e9;

        } //Fim-T

        // reanneling
        reanneling = 1;
    }
}

#endif