#ifndef _BRKGA_H
#define _BRKGA_H

/************************************************************************************
 Method: PARAMETRICUNIFORMCROSSOVER
 Description: create a new offspring with parametric uniform crossover
*************************************************************************************/
static TSol ParametricUniformCrossover(int eliteSize, int popSize, double rhoe, std::vector <TSol> &Pop)
{	
	TSol s;

    int eliteParent = irandomico(0, eliteSize - 1);                 // one chromosome from elite set
    int nonEliteParent = irandomico(eliteSize, popSize-1);          // one chromosome from nonelite  population

    // best fit parent is parent elite
    if (Pop[eliteParent].ofv > Pop[nonEliteParent].ofv){
        int temp = eliteParent;
        eliteParent = nonEliteParent;
        nonEliteParent = temp;
    }

	// create a new offspring
	s.rk.resize(n);

    // Mate
    for(int j = 0; j < n; j++)
    {
        //copy alelos of top chromossom of the new generation
        if (randomico(0,1) < rhoe){
            s.rk[j] = Pop[eliteParent].rk[j];
        }
        else{
            s.rk[j] = Pop[nonEliteParent].rk[j];
        }
    }

    // set the flag of local search as zero
    s.flag = 0;

    return s;
}

/************************************************************************************
 Method: BRKGA
 Description: search process of the BRKGA
*************************************************************************************/
void BRKGA(int method, int control)
{
    int p = 0;          	                    // size of population
    double pe = 0.0;              	            // fraction of population to be the elite-set
    double pm = 0.0;          	                // fraction of population to be replaced by mutants
    double rhoe = 0.0;              	        // probability that offspring inherit an allele from elite parent
 
    std::vector <TSol> Pop;                     // current population
    std::vector <TSol> PopInter;                // intermediary population
    TSol bestInd;                               // best individual found in past generation

    int numGenerations = 0;                     // number of generations
    int bestGeneration = 0;                     // generation in which found the best solution
    float currentTime = 0;                      // computational time of the search process

    struct timeval TstartMH, TendMH;            // computational time (unix systems)
    gettimeofday(&TstartMH, NULL);
    gettimeofday(&TendMH, NULL);


    // offline control
    if (control == 0){
        p    = 610;
        pe   = 0.20;                                                         
        pm   = 0.15;                                                    
        rhoe = 0.70;
    }
    
    // initialize population
    Pop.clear(); 
    PopInter.clear(); 
    Pop.resize(p);
    PopInter.resize(p);

    // Create the initial chromosomes with random keys
    for (int i=0; i<p; i++)
    {
        CreateInitialSolutions(Pop[i]); 
        Pop[i].ofv = Decoder(Pop[i]);
        PopInter[i] = Pop[i];
    }
    
    // sort population in increase order of fitness
    sort(Pop.begin(), Pop.end(), sortByFitness);

    // save the best solution found
    bestInd = Pop[0];
    
    // run the evolutionary process until stop criterion
    while (currentTime < MAXTIME)
    {
    	// number of generations
        numGenerations++;

        // The 'Pe' best chromosomes are maintained, so we just copy these into PopInter:
        for (int i=0; i<(int)(p*pe); i++){
            // copy the chromosome for next generation
            PopInter[i] = Pop[i]; 
        }  

        // We'll mate 'p - pe - pm' pairs; initially, i = pe, so we need to iterate until i < p - pm:
        TSol bestOff;
        bestOff.ofv = INFINITY;
        for (int i = (int)(p*pe); i < p - (int)(p*pm); i++){
            if (stop_execution.load()) return;      

            // Parametric uniform crossover with mutation
            PopInter[i] = ParametricUniformCrossover((int)(p*pe), p, rhoe, Pop);
 
            // Calculate the fitness of new chromosomes
            PopInter[i].ofv = Decoder(PopInter[i]); 

            if (PopInter[i].ofv < bestOff.ofv){
                bestOff = PopInter[i];
            }
        }
        
        // We'll introduce 'pm' mutants:
        for (int i = p - (int)(p*pm) - (int)(p*pe); i < p; i++){
            if (stop_execution.load()) return;      
            
            CreateInitialSolutions(PopInter[i]);
            PopInter[i].ofv = Decoder(PopInter[i]);
            // PopInter[i].ofv = 999999999;
        }  
                
        // Update the current population
        Pop = PopInter;   

        // Sort population in increase order of fitness
        sort(Pop.begin(), Pop.end(), sortByFitness);

        // We improve the best fitness in the current population 
        if (Pop[0].ofv < bestInd.ofv){
                RVND(Pop[0]);

            bestGeneration = numGenerations;
            bestInd = Pop[0];
        
            // Update pool of solutions
            UpdatePoolSolutions(Pop[0], method);
        }

        // if (debug && method == 0) printf("\nGen: %d \t best ind: %lf \t off: %lf", numGenerations, bestInd.ofv, bestOff.ofv);

        // ******************************* SHAKING *****************************
        if ((numGenerations - bestGeneration) > 5*n) {
            if (stop_execution.load()) return;      

            // reset the number of generations without improvement
            bestGeneration = numGenerations;

            for (int i=0; i<p; i++){

                // Shake the elite set
                if (i < (int)(p*pe)){
                    ShakeSolution(Pop[i], rhoe, rhoe);
                }

                // reset the non-elite chromosomes
                else{
                    CreateInitialSolutions(Pop[i]);
                }
                
                // decode the new solution
                Pop[i].ofv = Decoder(Pop[i]);
            }

            sort(Pop.begin(), Pop.end(), sortByFitness);
        }
        // *********************************************************************

        // terminate the evolutionary process in MAXTIME
        gettimeofday(&TendMH, NULL);
        currentTime = ((TendMH.tv_sec  - TstartMH.tv_sec) * 1000000u + TendMH.tv_usec - TstartMH.tv_usec) / 1.e6; 
    }

    // free memory of BRKGA components
    Pop.clear();
    PopInter.clear();
}

#endif