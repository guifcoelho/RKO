#ifndef _GA_H
#define _GA_H

/************************************************************************************
 Method: GA
 Description: search process of the Genetic Algorithm
*************************************************************************************/
void GA(int method, int control)
{
    int sizePop = 0;                    // population size
    float probCros = 0.0;               // crossover probability
    float probMut = 0.0;                // mutation probability

    std::vector <TSol> Pop;             // current population
    std::vector <TSol> PopInter;        // selected individuals to evolutionary process
    std::vector <TSol> PopNew;          // offsprings

    TSol bestInd;                       // best solution found in GA

    // useful local variables
    int numGenerations = 0;             // number of generations
    float currentTime = 0;              // computational time of the search process
    int bestGeneration = 0;             // generation in which found the best solution

    struct timeval TstartMH, TendMH;    // computational time (unix systems)
    gettimeofday(&TstartMH, NULL);
    gettimeofday(&TendMH, NULL);

     // offline control
    if (control == 0){
        // define parameters of GA
        sizePop = 500;
        probCros = 0.97;
        probMut = 0.005; 
    }

    // initialize population
    Pop.clear(); 
    PopNew.clear();
    PopInter.clear(); 

    Pop.resize(sizePop);
    PopNew.resize(sizePop);
    PopInter.resize(sizePop);

    // Create the initial chromosomes with random keys
    bestInd.ofv = 9999999999; 
    for (int i=0; i<sizePop; i++)
    {
        CreateInitialSolutions(Pop[i]); 
        Pop[i].ofv = Decoder(Pop[i]);  

        // set the best individual
        if (Pop[i].ofv < bestInd.ofv)
            bestInd = Pop[i];
    }
    
    // run the evolutionary process until stop criterion
    while (currentTime < MAXTIME)
    {
    	// number of generations
        numGenerations++;

        // Select parents using the TOURNAMENT selection method
        for (int j=0; j<sizePop; j++)
        {
            int p1, p2, p3, melhor;

            do{
                p1 = irandomico(0,sizePop-1);
                p2 = irandomico(0,sizePop-1);
                p3 = irandomico(0,sizePop-1);
            }while (p1 == p2 || p1 == p3 || p2 == p3);
            
            if ((Pop[p1].ofv < Pop[p2].ofv) && (Pop[p1].ofv < Pop[p3].ofv))
                melhor = p1;
            else
            if (Pop[p2].ofv < Pop[p3].ofv)
                melhor = p2;
            else
                melhor = p3;

            // insert the best individual into the intermediate population 
            PopInter[j] = Pop[melhor];
        }

        // perform uniform crossover between the selected pairs of parents
        double mediaFilhos = 0.0;
        for (int i=0; i<sizePop-1; i=i+2)
        {
            if (stop_execution.load()) return;      

            PopNew[i] = PopInter[i];
            PopNew[i+1] = PopInter[i+1];

            // check the probability of crossover
            if (randomico(0,1) < probCros)
            {
                for(int j = 0; j < n; j++)
                {
                    // swap the alleles of offspring i
                    if (randomico(0,1) < 0.5){
                        PopNew[i].rk[j] = PopInter[i+1].rk[j];
                    }

                    // swap the alleles of offspring i+1
                    if (randomico(0,1) < 0.5){
                        PopNew[i+1].rk[j] = PopInter[i].rk[j];
                    }

                    // check the probability of making a mutation in each gene of offspring i
                    if (randomico(0,1) <= probMut){
                        PopNew[i].rk[j] = randomico(0,1); 
                    }

                    // check the probability of making a mutation in each gene of offspring i+1
                    if (randomico(0,1) <= probMut){
                        PopNew[i+1].rk[j] = randomico(0,1); 
                    }
                }

                // calculate fitness
                PopNew[i].ofv   = Decoder(PopNew[i]);
                PopNew[i+1].ofv = Decoder(PopNew[i+1]);
            }

            mediaFilhos += PopNew[i].ofv + PopNew[i+1].ofv;

            // set the best individual
            if (PopNew[i].ofv < bestInd.ofv){
                bestInd = PopNew[i];
                bestGeneration = numGenerations;

                // update pool of solutions
                UpdatePoolSolutions(bestInd, method);
            }

            if (PopNew[i+1].ofv < bestInd.ofv){
                bestInd = PopNew[i+1];
                bestGeneration = numGenerations;

                // update pool of solutions
                UpdatePoolSolutions(bestInd, method);
            }
        }

        // Elitism - assigning the best individual in Pop
        int pos1 = irandomico(0,sizePop-1);
        if (randomico(0,1)<0.95)
            NelderMeadSearch(PopNew[pos1]);
        else 
            RVND(PopNew[pos1]);

        // set the best individual
        if (PopNew[pos1].ofv < bestInd.ofv){
            bestInd = PopNew[pos1];
            bestGeneration = numGenerations;

            // update pool of solutions
            UpdatePoolSolutions(bestInd, method);
        }
           
        // replace the population with offspring
        Pop = PopNew;

        // if (debug && method == 6) printf("\nGen: %d \t sBest: %lf \t average: %lf", numGenerations, bestInd.ofv, mediaFilhos/sizePop);

        // RESET
        if ((numGenerations - bestGeneration) > n*2) {
            Pop.clear(); 
            Pop.resize(sizePop);
            for (int i=0; i<sizePop; i++){
                if (stop_execution.load()) return;      
                
                CreateInitialSolutions(Pop[i]);
                Pop[i].ofv = Decoder(Pop[i]);
            }

            bestGeneration = numGenerations;
        }

        // terminate the evolutionary process in MAXTIME
        gettimeofday(&TendMH, NULL);
        currentTime = ((TendMH.tv_sec  - TstartMH.tv_sec) * 1000000u + TendMH.tv_usec - TstartMH.tv_usec) / 1.e6; 
    }

    // free memory of GA components
    Pop.clear();
    PopNew.clear();
    PopInter.clear();
}

#endif