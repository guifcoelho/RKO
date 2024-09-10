#ifndef _GRASP_H
#define _GRASP_H

/************************************************************************************
 Method: LineSearch
 Description: Consider the line search in direction ei , where vector ei has zeros 
 in all components except the i-th, where it has value one. The objective function 
 is evaluated at points x + k·h·ei for k = 0,1,−1,2,−2,... such that li≤xi+k·h≤ui.
 Let k∗ the value of k that minimizes f(x+k·h·ei) subject to li ≤ xi +k·h ≤ ui
*************************************************************************************/
static void LineSearch(TSol s, float h, int i, double &bestZ, double &bestF)
{
    std::vector<double> rk;

    // find the best solution in line
    bestZ = 0;
    bestF = INFINITY;

    // generate k as possible random keys of position i, k = rk_i * tau | tau = {0, 1, -1, 2, -2, ...}
    double tau = 0;
    rk.push_back(s.rk[i] + tau * h);
    for (int j=0; j<(int)(1.0/h)+1; j+=2){
        tau++;

        if (((s.rk[i] + tau * h) >= 0) && ((s.rk[i] + tau * h) < 1))
            rk.push_back(s.rk[i] + tau * h);

        if (((s.rk[i] + (-1*tau) * h) >= 0) && ((s.rk[i] + (-1*tau) * h) < 1))
            rk.push_back(s.rk[i] + (-1*tau) * h);
    }

    // (sample greedy) This method is similar to the greedy algorithm, but instead of selecting the best among all possible options, 
    // it only considers q < m possible insertions (chosen uniformly at random) in each iteration. The most profitable among those is selected. 
    int q = ceil(log2((int)(1.0/h))) + 1; 

    if (q > (int)rk.size())
        q = rk.size();

    // choose the subset with q rks to calculate the decoder
    std::shuffle(rk.begin(), rk.end(), std::mt19937(std::random_device()()));

    // decoder solution s with rk j
    for (int j=0; j<q; j++)
    {  
        if (stop_execution.load()) return;      
        
        s.rk[i] = rk[j];   
        s.ofv = Decoder(s);

        if (s.ofv < bestF){
            bestZ = s.rk[i];
            bestF = s.ofv;
        }
    }

    rk.clear();
}

/************************************************************************************
 Method: ConstrutiveGreedyRandomized
 Description: It takes as input a solution vector s. Initially, the procedure allows 
 all coordinates of s to change (i.e. they are called unfixed). In turn, a discrete 
 line search is performed in each unfixed coordinate direction i of s with the other 
 n − 1 coordinates of s held at their current values.
*************************************************************************************/
static void ConstrutiveGreedyRandomized(TSol &s, float h, float alpha)
{
    // ** nao estou alterando o decoder
    // A fase construtiva inicia da solução corrente s e fará uma “perturbação” desta solução com uma intensidade entre Beta_min e Beta_max. 
    // Em cada iteração da fase construtiva escolhe-se aleatoriamente um conjunto de random keys que ainda não foram fixadas e para cada 
    // random key deste conjunto gera-se todos os possíveis valores de novas random keys baseado no valor de h. Escolhe-se aleatoriamente 
    // um destes valores e calcula o decoder (esta será a saída do line search). Com a lista de novas random keys geradas, construimos a RCL 
    // e selecionamos um valor aleatoriamente entre os melhores, fixando esta random key.

    std::vector<int> UnFixed(n);                // store random-keys not fixed
    std::vector<int> chosenRK;                  // store random-keys for the line search
    std::vector<int> RCL;                       // store the best candidate solutions
    std::vector<double> z(n);                   // store the best value of random-key i
    std::vector<double> g(n,INFINITY);          // store the value of the ofv for the solution s with random-key z_i

    double min, max;
    double betaMin = 0.5, 
           betaMax = 0.8;

    // initialize the points on the solution that can be changed
    for (int i = 0; i < n; i++) {
        UnFixed[i] = i;
    }

    // construct a solution by perturbing the current solution and choosing one of the best
    double intensity = randomico(betaMin, betaMax);
    
    // while (!UnFixed.empty())
    for (int j=0; j<n*intensity; j++)
    {
        // create a list of candidate solutions by perturbing a (not yet 'fixed') rk of the current solution
        min = INFINITY;
        max = -INFINITY;

        // choose the subset of random keys that will be searched
        chosenRK.clear();
        int kMax = UnFixed.size() * 0.1;
        if (kMax < 2) 
            kMax = UnFixed.size();

        chosenRK = UnFixed;
        std::shuffle(chosenRK.begin(), chosenRK.end(), rng);
        chosenRK.resize(kMax);

        // line search
        z.clear();
        z.resize(n);
        g.clear();
        g.resize(n,INFINITY);

        for (int k=0; k<kMax; k++) 
        {
            int i = chosenRK[k];

            TSol sAux = s;
            
            // linear search
            LineSearch(sAux, h, i, z[i], g[i]);
        }

        for(int i=0; i<n; i++)
        {
            if (min > g[i])
                min = g[i];

            if (max < g[i] && g[i] != INFINITY)
                max = g[i];
        }

        // create RCL
        RCL.clear();
        double threshold = min + alpha * (max - min);

        for (int i=0; i<n; i++)
        { 
            if (g[i] <= threshold)
            {
                RCL.push_back(i);
            }
        }

        // randomly select one of the best candidates to continue construction
        if (!RCL.empty()){
            int x = irandomico(0, (int)(RCL.size())-1);
            int kCurrent = RCL[x];  // index of the random key that will be fixed
            
            // update current solution
            s.rk[kCurrent] = z[kCurrent];
            s.ofv = g[kCurrent];

            // remove rk k from the set UnFixed
            for (int i=0; i<(int)UnFixed.size(); i++)
            {
                if (UnFixed[i] == kCurrent) 
                {
                    UnFixed.erase(UnFixed.begin()+i);
                    break;
                }
            }
        }
    }

    // update the solution found in the constructive phase
    s.ofv = Decoder(s);
}

/************************************************************************************
 Method: GRASP
 Description: Metaheurist Greedy Randomized Adpative Search Procedura.
*************************************************************************************/
void GRASP(int method, int control)
{
    // GRASP parameters
    TSol s;                              // current solution
    TSol sLine;                          // constructive solution
    TSol sLineBest;                      // local search solution
    TSol sBest;                          // best solution of GRASP

    float alphaGrasp = 0.1;              // greedy rate

    float currentTime = 0;               // computational time of the search process

    float h;                             // grid dense
    float hs = 0.12500;                  // start grid dense
    float he = 0.00012;                  // end grid dense

    struct timespec TstartMH, TendMH;               // computational time (unix systems)
    clock_gettime(CLOCK_MONOTONIC, &TstartMH);
    clock_gettime(CLOCK_MONOTONIC, &TendMH);

    // offline control
    if (control == 0){
        // define parameters of GRASP
        alphaGrasp = randomico(0.1, 0.9);
    }

    int iter = 0;

    // create an initial solution
    CreateInitialSolutions(s);
    s.ofv = Decoder(s);

    sBest = s;

    // run the search process until stop criterion (maxTime)
    while (currentTime < MAXTIME)
    {
        h = hs;
        int noImprov = 0;
        while (h >= he && currentTime < MAXTIME)
        {
            if (stop_execution.load()) return; 

            iter++;

            // offline control
            if (control == 0){
                alphaGrasp = randomico(0.1, 0.9);
            }

            // construct a greedy randomized solution
            sLine = s;
            ConstrutiveGreedyRandomized(sLine, h, alphaGrasp);
            
            // apply local search in current solution
            sLineBest = sLine;
            RVND(sLineBest);
            
            // update the best solution found by GRASP
            if (sLineBest.ofv < sBest.ofv){
                sBest = sLineBest;
                noImprov = 0;

                // update the pool of solutions
                UpdatePoolSolutions(sLineBest, method);
            }
            // make grid more dense
            else{
                noImprov++;

                if (noImprov > n/2)
                {
                    h = h/2;
                    noImprov = 0;
                }
            }

            // accept criterion
            if (sLineBest.ofv < s.ofv){
                s = sLineBest;
            }
            else{
                // Metropolis criterion
                if ( randomico(0.0,1.0) < (exp(-(sLineBest.ofv - s.ofv)/(1000 - 1000*(currentTime / MAXTIME)))) ){ 
                    s = sLineBest;
                }
            }

            // if (debug && method == 2) printf("\nIter: %d \t s': %lf \t s'best: %lf \t sBest: %lf", iter, sLine.ofv, sLineBest.ofv, sBest.ofv);

            // terminate the evolutionary process in MAXTIME
            clock_gettime(CLOCK_MONOTONIC, &TendMH);
            currentTime = (TendMH.tv_sec - TstartMH.tv_sec) + (TendMH.tv_nsec - TstartMH.tv_nsec) / 1e9;
        } 

        // restart the current solution 
        CreateInitialSolutions(s); 

        // decode the new solution
        s.ofv = Decoder(s);

        // terminate the search process in MAXTIME
        clock_gettime(CLOCK_MONOTONIC, &TendMH);
        currentTime = (TendMH.tv_sec - TstartMH.tv_sec) + (TendMH.tv_nsec - TstartMH.tv_nsec) / 1e9;
    }
}

#endif
