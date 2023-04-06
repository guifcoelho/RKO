#include "BRKGA.h"

/************************************************************************************
			                  GENERAL FUNCTIONS
*************************************************************************************/

// Sort TSol by objective function
bool sortByFitnessBRKGA(const TSol &lhs, const TSol &rhs) { return lhs.ofv < rhs.ofv; }


void BRKGA()
{
    // free memory with problem data
    FreeMemoryProblem();

    // free memory of BRKGA components
    FreeMemory();

    //read data of the instance
    ReadData(nameTable);
    
    // initialize population
    Pop.clear(); 
    PopInter.clear(); 

    // define parameters for classic BRKGA                                    
    p       = 987;                                                          
    pe      = 0.20;                                                         
    pm      = 0.15;                                                         
    rhoe    = 0.70;                                                          

    Pop.resize(p);
    PopInter.resize(p);

    // Create the initial chromosomes with random keys
    for (int i=0; i<p; i++)
    {
        TSol ind = CreateInitialSolutions(); 
        Decoder(ind);
        Pop[i] = PopInter[i] = ind;
    }
    
    // sort population in increase order of fitness
    sort(Pop.begin(), Pop.end(), sortByFitnessBRKGA);

    // save the best solution found
    updateBestSolution(Pop[0]);
    
    // useful local variables
    int numGenerations = 0;             // number of generations
    int bestGeneration = 0;             // generation in which found the best solution
    double bestFitness = Pop[0].ofv;    // best fitness found in past generation
    float currentTime = 0;              // computational time of the search process
    int sumLS = 0;                      // number of local search applied in each generation
    int noImprov = 0;                   // number of generations without improvement in the best solution
    
    // run the evolutionary process until stop criterion
    while(1)
    {
    	// number of generations
        numGenerations++;

        // number of generations without improvement in the best solution
        noImprov++;

        // The 'Pe' best chromosomes are maintained, so we just copy these into PopInter:
        #pragma omp parallel for num_threads(MAX_THREADS)
        for (int i=0; i<(int)(p*pe); i++){
            // copy the chromosome for next generation
            PopInter[i] = Pop[i]; 
        }  

        // We'll mate 'p - pe - pm' pairs; initially, i = pe, so we need to iterate until i < p - pm:
        #pragma omp parallel for num_threads(MAX_THREADS)
        for (int i = (int)(p*pe); i < p - (int)(p*pm); i++){
            // Parametric uniform crossover with mutation
            PopInter[i] = ParametricUniformCrossover((int)(p*pe), p);
 
            // Calculate the fitness of new chromosomes
            Decoder(PopInter[i]); 
        }
        
        // We'll introduce 'pm' mutants:
        #pragma omp parallel for num_threads(MAX_THREADS)
        for (int i = p - (int)(p*pm) - (int)(p*pe); i < p; i++){
            PopInter[i] = CreateInitialSolutions();
            Decoder(PopInter[i]);
        }  
                
        // Update the current population
        Pop = PopInter;   

        // Sort population in increase order of fitness
        sort(Pop.begin(), Pop.end(), sortByFitnessBRKGA);
        updateBestSolution(Pop[0]);

        // We improve the best fitness in the current population 
        if (Pop[0].ofv < bestFitness){
            bestFitness = Pop[0].ofv;
            bestGeneration = numGenerations;
            noImprov = 0;
        }
        
        // ********************* LOCAL SEARCH IN COMMUNITIES *******************
        sumLS = 0;

        // Verify if there are local search heuristics 
        if (numLS > 0){   

            //apply local search when BRKGA found a new better solution or n*pe generations without improvements
            if (noImprov == 0 || noImprov > (int)n*pe){

                // restart the count of generations without improvements (local search)
                noImprov = 0;

                // Identify commuties in the Elite with Label Propagation method
	            IC();

	            std::vector <int> promisingSol; 
                promisingSol.clear();

	            for (int i=0; i < (int)(p*pe); i++) {
                    
                    // insert the individual index in the promising list
	                if (Pop[i].promising == 1){
	                	promisingSol.push_back(i);
	                }
                    
                    // generate caotic individual (crossover between one elite and one mutant)
                    else if (i > 0){
                        ChaoticInd(Pop[i]);
                        Decoder(Pop[i]);

                        // set flag as 0 to permit new local search
                        Pop[i].flag = 0;
                    }
	            }

	            #pragma omp parallel for num_threads(MAX_THREADS)
                for (unsigned int i=0; i < promisingSol.size(); i++){

                    // local search not influence the evolutionary process
                    TSol s = Pop[promisingSol[i]];
                    LocalSearch(s);
                    updateBestSolution(s);                    

                    // set flag as 1 to prevent new local search in the same solution
                    Pop[promisingSol[i]].flag = 1;
                }

                sumLS = promisingSol.size();
                promisingSol.clear();

                sort(Pop.begin(), Pop.end(), sortByFitnessBRKGA);
                updateBestSolution(Pop[0]);
	        }
	    }
        // *********************************************************************


        // ******************************* SHAKING *****************************
        if ((numGenerations - bestGeneration) > n) {
            
            if (debug) printf("\n\nShaking elite and reset non-elite...\n\n");

            // reset the number of generations without improvement
            bestGeneration = numGenerations;

            // Shake the elite set
            float shaking_type = 0.0;
            int intensity = n*rhoe;
            for(int e = 0; e < (int)(p*pe); e++) {
                for(int k = 0; k < intensity; k++) {
                    shaking_type = irandomico(1,4);
                    int i = irandomico(0, n - 1);
                    if(shaking_type == 1){
                        // Invert value
                        if (Pop[e].vec[i].rk > 0.0001)
                            Pop[e].vec[i].rk = 1.0 - Pop[e].vec[i].rk;
                        else
                            Pop[e].vec[i].rk = 0.9999;
                    }
                    else 
                    if (shaking_type == 2){
                        // Swap two random positions
                        int j = irandomico(0, n - 1);
                        double temp = Pop[e].vec[i].rk;
                        Pop[e].vec[i].rk = Pop[e].vec[j].rk;
                        Pop[e].vec[j].rk = temp;
                    }
                    else
                    if(shaking_type == 3){
                        // Change to random value
                        Pop[e].vec[i].rk = randomico(0,1);
                    }
                    i = irandomico(0, n - 2);
                    if(shaking_type == 4){
                        // Swap with neighbor
                        double temp = Pop[e].vec[i].rk;
                        Pop[e].vec[i].rk = Pop[e].vec[i+1].rk;
                        Pop[e].vec[i+1].rk = temp;
                    }
                }
                Decoder(Pop[e]);
            }

            // reset the non-elite chromosomes
            for (int i=(int)(p*pe); i<p; i++){
                Pop[i] = CreateInitialSolutions();
                Decoder(Pop[i]);
            }

            sort(Pop.begin(), Pop.end(), sortByFitnessBRKGA);
            updateBestSolution(Pop[0]);
            bestFitness = Pop[0].ofv;
        }
        // *********************************************************************

        // print screen 
        if (debug){
            printf("\nGeneration: %d [%d] \t %.2lf (%.2lf)  \t %.2lf",
                        numGenerations, sumLS, bestSolution.ofv, bestSolution.vec[n].rk, bestFitness);
        }

        // terminate the evolutionary process in MAXTIME
        gettimeofday(&Tend, NULL);
        currentTime = ((Tend.tv_sec  - Tstart.tv_sec) * 1000000u + Tend.tv_usec - Tstart.tv_usec) / 1.e6; 
        
        // stop criterium
        if (currentTime >= MAXTIME || bestSolution.ofv <= OPTIMAL){
            break;
        }
    }

    // free memory with problem data
    FreeMemoryProblem();

    // free memory of BRKGA-QL components
    FreeMemory();
}

void updateBestSolution(TSol s)
{
    // save the best solution found 
    if (s.ofv < bestSolution.ofv)
    {
        bestSolution = s;
        gettimeofday(&Tbest, NULL);
    }
}

TSol CreateInitialSolutions()
{
	TSol s;
	TVecSol aux;

    s.vec.clear();

	// create a random-key for each allelo (consider decoder type in the n-th random-key)
	for (int j = 0; j < n+1; j++)
	{
        aux.rk  = randomico(0,1);  // random value between [0,1[
        aux.sol = 0;
        s.vec.push_back(aux);
	}

    // flag to control the local search memory
    s.flag = 0;

	return s;
}

void ChaoticInd(TSol &s)
{
    // generate a caotic individual
    for (int k=0; k<n+1; k++)
    {      
        if (randomico(0,1) > rhoe)
           s.vec[k].rk = randomico(0,1);
    }

    // set the flag of local search as zero
    s.flag = 0;
}

TSol ParametricUniformCrossover(int eliteSize, int popSize)
{	
	TSol s;

    int eliteParent = irandomico(0, eliteSize - 1);                 // one chromosome from elite set
    int nonEliteParent = irandomico(0, popSize-1);                  // one chromosome from entire population

    // best fit parent is parent elite
    if (Pop[eliteParent].ofv > Pop[nonEliteParent].ofv){
        int temp = eliteParent;
        eliteParent = nonEliteParent;
        nonEliteParent = temp;
    }

	// create a new offspring
	s.vec.resize(n+1);

    // Mate: including decoder gene in the n-th rk 
    for(int j = 0; j < n+1; j++)
    {
        //copy alelos of top chromossom of the new generation
        if (randomico(0,1) < rhoe){
            s.vec[j].rk = Pop[eliteParent].vec[j].rk;
        }
        else{
            s.vec[j].rk = Pop[nonEliteParent].vec[j].rk;
        }
    }

    // set the flag of local search as zero
    s.flag = 0;

    return s;
}

double PearsonCorrelation(std::vector <TVecSol> X, std::vector <TVecSol> Y)
{
    double correlation = 0;
    double sumXY = 0;
    double sumX2 = 0;
    double sumY2 = 0;
    double sumX = 0;
    double sumY = 0;

    for(int j=0; j<n; j++)
    {
        sumX += X[j].rk;
        sumX2 += X[j].rk * X[j].rk;
        sumXY += X[j].rk * Y[j].rk;
        sumY += Y[j].rk;
        sumY2 += Y[j].rk * Y[j].rk;
    }

    //Pearson
    correlation= ((n*sumXY) - (sumX*sumY) ) / (sqrt( (n*sumX2 - sumX*sumX) * (n*sumY2 - sumY*sumY) ));
    return correlation;
}

void IC() 
{
    int Tpe = (int)p*pe;
    std::vector<std::vector<std::pair<int, double> > > listaArestas(Tpe, std::vector<std::pair<int, double> >());

	// create weighted (pearson correlation) graph
	int entrouAresta = 0;
	double pearson = 0.0;
	for (int i = 0; i < Tpe - 1; i++) {
		for (int j = i + 1; j < Tpe; j++)
		{
			pearson = PearsonCorrelation(Pop[i].vec, Pop[j].vec);
			if (pearson > 0.6) {
				entrouAresta++;
				listaArestas[i].push_back(std::make_pair(j, pearson));
				listaArestas[j].push_back(std::make_pair(i, pearson));
			}
		}
	}

	// apply clustering method
	LP(listaArestas);

	PromisingLP();
    listaArestas.clear();
}

void LP(std::vector<std::vector<std::pair<int, double> > > listaArestas)
{
    int nk = listaArestas.size();

	// Create vector with visit order
	std::vector<int> ordemVisita(nk);
	iota(ordemVisita.begin(), ordemVisita.end(), 0);

	// initialize each node with its own label
	for (int i = 0; i < nk; i++)
		Pop[i].label = i;

	int iteracao = 1;
	int labelVizinho, melhorLabel;
	double melhorPeso;
	std::map<int, double> totalLabels;
	std::map<int, double>::iterator it;

	int movimentos = 1;
	while (movimentos) {
		movimentos = 0;
		random_shuffle(ordemVisita.begin(), ordemVisita.end());
		for (auto idVertice : ordemVisita) {

			// Calculate the weigth of the labels
			totalLabels.clear();
			for (auto idVizinho : listaArestas[idVertice]) {
				labelVizinho = Pop[idVizinho.first].label;
				it = totalLabels.find(labelVizinho);
				if (it != totalLabels.end()) {
					it->second += idVizinho.second;
				}
				else {
					totalLabels[labelVizinho] = idVizinho.second;
				}
			}

			// Best label is itself initially
			melhorLabel = Pop[idVertice].label;
			melhorPeso = std::numeric_limits<double>::min();
			for (auto totais : totalLabels) {
				if (totais.second > melhorPeso) {
					melhorLabel = totais.first;
					melhorPeso = totais.second;
				}
			}

			if (melhorLabel != Pop[idVertice].label) {
				Pop[idVertice].label = melhorLabel;
				movimentos = 1;
			}
		}
		iteracao++;
	}

    ordemVisita.clear();
}

void PromisingLP()
{
    int Tpe = (int)p*pe;
    std::vector <int> grupos;
	int tamanhoGrupos = 0;

	// initialize promisings solutions
	for (int i = 0; i < Tpe; i++)
		Pop[i].promising = 0;

	// save labels defined by LP in groups
	int achei;

    for (int i = 0; i < Tpe; i++)
	{
		achei = 0;
		for (unsigned int j = 0; j < grupos.size(); j++)
		{
			if (Pop[i].label == grupos[j])
				achei = 1;
		}
		if (achei == 0)
		{
			tamanhoGrupos++;
			grupos.push_back(Pop[i].label);
		}
	}

	// find the best solution in the group (with flag = 0)
	for (unsigned int j = 0; j < grupos.size(); j++)
	{
		double menorFO = INFINITY;
		int localMenor = -1;
		int local = -1;
		for (int i = 0; i < Tpe; i++)
		{
			if (Pop[i].label == grupos[j])
			{
				// find the best solution of the group
				if (local == -1)
					local = i;

				// we not apply local search in this solution yet
                if (Pop[i].ofv < menorFO && Pop[i].flag == 0) 
				{
					menorFO = Pop[i].ofv;
					localMenor = i;
				}
			}
		}

		if (localMenor == -1)
			localMenor = local;

		if (Pop[localMenor].label != -1)
			Pop[localMenor].promising = 1;
	}
}

void FreeMemory()
{
    //methods
    Pop.clear();
    PopInter.clear();
}

double randomico(double min, double max)
{
	return ((double)(rand()%10000)/10000.0)*(max-min)+min;
}

int irandomico(int min, int max)
{
	return (int)randomico(0,max-min+1.0) + min;
}