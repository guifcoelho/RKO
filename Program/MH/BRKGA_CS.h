#ifndef _BRKGA_CS_H
#define _BRKGA_CS_H

/************************************************************************************
			                  DECLARATION
*************************************************************************************/

/************************************************************************************
 Method: BRKGA_QL
 Description: the evolutionary process of the BRKGA-QL
*************************************************************************************/
void BRKGA_CS(int method, int control);

/************************************************************************************
 Method: ChaoticInd
 Description: create a solution between a mutant and a elite individual
*************************************************************************************/
static void ChaoticInd(TSol &s, int rhoe);

/************************************************************************************
 Method: PARAMETRICUNIFORMCROSSOVER
 Description: create a new offspring with parametric uniform crossover
*************************************************************************************/
static TSol ParametricUniformCrossover(int elitesize, int popSize, double pm, double rhoe, std::vector <TSol> &Pop);

/************************************************************************************
 Method: PEARSON CORRELATION
 Description: calculate the Pearson correlation coefficient between two chromossoms
*************************************************************************************/
static double PearsonCorrelation(std::vector <double> s1, std::vector <double> s2);

/************************************************************************************
 Metodo: IC(TSol Pop)
 Description: apply clustering method to find promising solutions in the population
*************************************************************************************/
static void IC(int p, double pe, std::vector <TSol> &Pop);

/************************************************************************************
 Method: LP
 Description: Apply Label Propagation to find communities in the population
*************************************************************************************/
static void LP(std::vector<std::vector<std::pair<int, double> > > listaArestas, std::vector <TSol> &Pop);

/************************************************************************************
 Method: PROMISINGLP
 Description: Find the promising solutions to represent the communities
*************************************************************************************/
static void PromisingLP(int p, double pe, std::vector <TSol> &Pop);


/************************************************************************************
			                  IMPLEMENTATION
*************************************************************************************/

void BRKGA_CS(int method, int control)
{
    // BRKGA parameters
    int p = 1597;          	                    // size of population
    double pe = 0.20;              	            // fraction of population to be the elite-set
    double pm = 0.03;          	                    // fraction of population to be replaced by mutants
    double rhoe = 0.70;             	            // probability that offspring inherit an allele from elite parent

    // BRKGA variables
    std::vector <TSol> Pop;                         // current population
    std::vector <TSol> PopInter;                    // intermediary population
    TSol bestInd;                                   // best solution of the BRKGA-CS

    struct timespec TstartMH, TendMH;               // computational time (unix systems)
    clock_gettime(CLOCK_MONOTONIC, &TstartMH);
    clock_gettime(CLOCK_MONOTONIC, &TendMH);

    // Run
    int numGenerations = 0;                         // number of generations
    int bestGeneration = 0;                         // generation in which found the best solution
    float currentTime = 0;                          // computational time of the search process
    int sumLS = 0;                                  // number of local search applied in each generation
    int noImprovBRKGA = 0;                          // number of generations without improvement in the best solution

    // offline control
    if (control == 0)
    {
        p    = 1597;
        pe   = 0.30;                                                         
        pm   = 0.001;                                                    
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
    bestInd = Pop[0];

    // run the evolutionary process until stop criterion
    while (currentTime < MAXTIME)
    {
    	// number of generations
        numGenerations++;

        // number of generations without improvement in the best solution
        noImprovBRKGA++;

        // The 'Pe' best chromosomes are maintained, so we just copy these into PopInter:
        for (int i=0; i<(int)(p*pe); i++){
            // copy the chromosome for next generation
            PopInter[i] = Pop[i]; 
        }  

        // We'll mate 'P - Pe' pairs; initially, i = p*pe, so we need to iterate until i < p:
        double bestOFV = INFINITY;
        for (int i = (int)(p*pe); i < p; i++){            
            // Parametric uniform crossover with mutation
            PopInter[i] = ParametricUniformCrossover((int)(p*pe), p, pm, rhoe, Pop);
 
            // Calculate the fitness of new chromosomes
            PopInter[i].ofv = Decoder(PopInter[i]); 

            if (PopInter[i].ofv < bestOFV)
                bestOFV = PopInter[i].ofv;

            if (stop_execution.load()) return;  
        }
                
        // Update the current population
        Pop = PopInter;   

        // Sort population in increase order of fitness
        sort(Pop.begin(), Pop.end(), sortByFitness);

        // We improve the best fitness in the current population 
        if (Pop[0].ofv < bestInd.ofv){
            bestInd = Pop[0];
            bestGeneration = numGenerations;
            noImprovBRKGA = 0;

            // update the pool of solutions
            UpdatePoolSolutions(bestInd, method);
        }

        // if (debug) printf("\n%d: \tbestInd: %lf", numGenerations, bestInd.ofv);

        // ********************* LOCAL SEARCH IN COMMUNITIES *******************
        sumLS = 0;

        //apply local search when BRKGA found a new better solution or n*pe generations without improvements
        if (noImprovBRKGA == 0 || noImprovBRKGA > 10)
        {
            // restart the count of generations without improvements (local search)
            noImprovBRKGA = 0;

            // Identify commuties in the Elite with Label Propagation method
            IC(p, pe, Pop);

            std::vector <int> promisingSol; 
            promisingSol.clear();

            for (int i=0; i < (int)(p*pe); i++) {
                if (stop_execution.load()) return;  

                // insert the individual index in the promising list
                if (Pop[i].promising == 1){
                    promisingSol.push_back(i);
                }
                
                // generate caotic individual (crossover between one elite and one mutant)
                else {
                    ChaoticInd(Pop[i], rhoe);
                    Pop[i].ofv = Decoder(Pop[i]);

                    // set flag as 0 to permit new local search
                    Pop[i].flag = 0;
                }
            }

            for (unsigned int i=0; i < promisingSol.size(); i++){
                if (stop_execution.load()) return;      

                // local search
                if (i < 1)
                   RVND(Pop[promisingSol[i]]);
                else
                   NelderMeadSearch(Pop[promisingSol[i]]);

                // set flag as 1 to prevent new local search in the same solution
                Pop[promisingSol[i]].flag = 1;          

                // if (debug) printf("\n%d: \tBL: %lf", i, Pop[promisingSol[i]].ofv);
            }

            sumLS = promisingSol.size();
            promisingSol.clear();

            sort(Pop.begin(), Pop.end(), sortByFitness);
        }
        // *********************************************************************


        // ******************************* SHAKING *****************************
        if ((numGenerations - bestGeneration) > 5*n) {
            if (stop_execution.load()) return;  

            // reset the number of generations without improvement
            bestGeneration = numGenerations;

            // Shake the elite set
            for(int e = 0; e < (int)(p*pe); e++) {
                ShakeSolution(Pop[e],rhoe,rhoe);
                Pop[e].ofv = Decoder(Pop[e]);
            }

            // reset the non-elite chromosomes
            for (int i=(int)(p*pe); i<p; i++){
                CreateInitialSolutions(Pop[i]);
                Pop[i].ofv = Decoder(Pop[i]);
            }

            sort(Pop.begin(), Pop.end(), sortByFitness);
        }

        // ******************************* RESET *****************************

        // *********************************************************************

        // terminate the evolutionary process in MAXTIME
        clock_gettime(CLOCK_MONOTONIC, &TendMH);
        currentTime = (TendMH.tv_sec - TstartMH.tv_sec) + (TendMH.tv_nsec - TstartMH.tv_nsec) / 1e9;
    }

    // free memory of BRKGA-CS components
    Pop.clear();
    PopInter.clear();
}

static void ChaoticInd(TSol &s, int rhoe)
{
    // generate a caotic individual
    for (int k=0; k<n; k++)
    {      
        if (randomico(0,1) > rhoe)
           s.rk[k] = randomico(0,1);
    }

    // set the flag of local search as zero
    s.flag = 0;
}

static TSol ParametricUniformCrossover(int eliteSize, int popSize, double pm, double rhoe, std::vector <TSol> &Pop)
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
	s.rk.resize(n);

    // Mate
    for(int j = 0; j < n; j++)
    {
        // mutation
        if (randomico(0,1) < pm)
        {
            s.rk[j] = randomico(0,1);
        }
        else
        {
            //copy alelos of top chromossom of the new generation
            if (randomico(0,1) < rhoe){
                s.rk[j] = Pop[eliteParent].rk[j];
            }
            else{
                s.rk[j] = Pop[nonEliteParent].rk[j];
            }
        }
    }

    // set the flag of local search as zero
    s.flag = 0;

    return s;
}

static double PearsonCorrelation(std::vector <double> X, std::vector <double> Y)
{
    double correlation = 0;
    double sumXY = 0;
    double sumX2 = 0;
    double sumY2 = 0;
    double sumX = 0;
    double sumY = 0;

    for(int j=0; j<n; j++)
    {
        sumX += X[j];
        sumX2 += X[j] * X[j];
        sumXY += X[j] * Y[j];
        sumY += Y[j];
        sumY2 += Y[j] * Y[j];
    }

    //Pearson
    correlation= ((n*sumXY) - (sumX*sumY) ) / (sqrt( (n*sumX2 - sumX*sumX) * (n*sumY2 - sumY*sumY) ));
    return correlation;
}

static void IC(int p, double pe, std::vector <TSol> &Pop) 
{
    int Tpe = (int)p*pe;
    std::vector<std::vector<std::pair<int, double> > > listaArestas(Tpe, std::vector<std::pair<int, double> >());

    // pearson correlation factor
    double sigma = 0.6;                      

	// create weighted (pearson correlation) graph
	int entrouAresta = 0;
	double pearson = 0.0;
	for (int i = 0; i < Tpe - 1; i++) {
		for (int j = i + 1; j < Tpe; j++)
		{
			pearson = PearsonCorrelation(Pop[i].rk, Pop[j].rk);
			if (pearson > sigma) {
				entrouAresta++;
				listaArestas[i].push_back(std::make_pair(j, pearson));
				listaArestas[j].push_back(std::make_pair(i, pearson));
			}
            else{
                entrouAresta += 5;
            }
		}
	}

	// apply clustering method
	LP(listaArestas, Pop);

	PromisingLP(p, pe, Pop);
    listaArestas.clear();
}

static void LP(std::vector<std::vector<std::pair<int, double> > > listaArestas, std::vector <TSol> &Pop)
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
	while (movimentos) 
        {
		movimentos = 0;
		shuffle(ordemVisita.begin(), ordemVisita.end(),std::mt19937(std::random_device()()));
	        for (std::vector<int>::size_type idVertice=0; idVertice <ordemVisita.size(); idVertice++)
	        {
		    // Calculate the weigth of the labels
		    totalLabels.clear();
	            for (std::vector<std::pair<int, double> >::iterator itVizinho = listaArestas[idVertice].begin(); itVizinho != listaArestas[idVertice].end(); ++itVizinho) 
		    {
	                int idVizinho = itVizinho->first;
	                labelVizinho = Pop[idVizinho].label;
	                it = totalLabels.find(labelVizinho);
	                if (it != totalLabels.end()) {
	                    it->second += itVizinho->second;
	                } 
			else {
	                    totalLabels[labelVizinho] = itVizinho->second;
	                }
	            }

		    // Best label is itself initially
		    melhorLabel = Pop[idVertice].label;
		    melhorPeso = std::numeric_limits<double>::min();
	            for (std::map<int, double>::iterator itTotais = totalLabels.begin(); itTotais != totalLabels.end(); ++itTotais) {
	                if (itTotais->second > melhorPeso) {
	                    melhorLabel = itTotais->first;
	                    melhorPeso = itTotais->second;
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

static void PromisingLP(int p, double pe, std::vector <TSol> &Pop)
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


#endif
