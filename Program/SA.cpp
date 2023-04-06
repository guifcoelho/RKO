#include "SA.h"

/************************************************************************************
			                  GENERAL FUNCTIONS
*************************************************************************************/

void SA()
{
    // free memory with problem data
    FreeMemoryProblem();

    // free memory of BRKGA components
    FreeMemory();

    //read data of the instance
    ReadData(nameTable);

    float currentTime = 0;              // computational time of the search process

    // run the search process until stop criterion
    while(1)
    {
        // Create the initial solution with random keys
        sMelhor.ofv = 99999999999;
        for (int i=0; i<1000; i++)
        {
            TSol ind = CreateInitialSolutions(); 
            Decoder(ind);
            if (ind.ofv < sMelhor.ofv)
                sMelhor = ind;
        }
        
        // update the best solution found
        updateBestSolution(sMelhor);

    
        // define parameters
        T0 = 10000000;
        SAmax = 500;
        alfa = 0.98;

        // current solution
        s = sMelhor;
        IterT = 0;
        T = T0;

        while (T > 0.0001)
        {
            while (IterT < SAmax)
            {
                IterT++;

                // Shake the current solution
                sViz = s;
                int shaking_type = 0.0;
                int intensity = n*0.05;
                for(int k = 0; k < intensity; k++) 
                {
                    shaking_type = irandomico(1,4); //printf("\nType: %d\n",shaking_type);
                    int i = irandomico(0, n - 1);
                    if(shaking_type == 1){
                        // Invert value
                        if (sViz.vec[i].rk > 0.0001)
                            sViz.vec[i].rk = 1.0 - sViz.vec[i].rk;
                        else
                            sViz.vec[i].rk = 0.9999;
                    }
                    else 
                    if (shaking_type == 2){
                        // Swap two random positions
                        int j = irandomico(0, n - 1);
                        double temp = sViz.vec[i].rk;
                        sViz.vec[i].rk = sViz.vec[j].rk;
                        sViz.vec[j].rk = temp;
                    }
                    else
                    if(shaking_type == 3){
                        // Change to random value
                        sViz.vec[i].rk = randomico(0,1);
                    }
                    i = irandomico(0, n - 2);
                    if(shaking_type == 4){
                        // Swap with neighbor
                        double temp = sViz.vec[i].rk;
                        sViz.vec[i].rk = sViz.vec[i+1].rk;
                        sViz.vec[i+1].rk = temp;
                    }
                }

                // sortear um decoder
                // sViz.vec[n].rk = irandomico(0,1);

                // calcular a FO
                Decoder(sViz);
                
                // calcular o delta
                delta = sViz.ofv - s.ofv;

                //definir a partir de qual solucao continuar a busca
                if (delta < 0)
                {
                    s = sViz;
                    updateBestSolution(s);
                }
                else
                {
                    double x = randomico(0,1);

                    if ( x < (exp(-delta/T)) )       
                        s = sViz;
                }
            } //Fim-SAmax

            if (debug)
                printf("\nT = %.4f \t FO = %.2lf \t melhorFO = %.2lf", T, s.ofv, bestSolution.ofv);

            T = T * alfa;
            IterT = 0;

            // terminate the search process in MAXTIME
            gettimeofday(&Tend, NULL);
            currentTime = ((Tend.tv_sec  - Tstart.tv_sec) * 1000000u + Tend.tv_usec - Tstart.tv_usec) / 1.e6; 
            
            // stop criterium
            if (currentTime >= MAXTIME || bestSolution.ofv <= OPTIMAL){
                break;
            }

        } //Fim-T

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

// double PearsonCorrelation(std::vector <TVecSol> X, std::vector <TVecSol> Y)
// {
//     double correlation = 0;
//     double sumXY = 0;
//     double sumX2 = 0;
//     double sumY2 = 0;
//     double sumX = 0;
//     double sumY = 0;
//     for(int j=0; j<n; j++)
//     {
//         sumX += X[j].rk;
//         sumX2 += X[j].rk * X[j].rk;
//         sumXY += X[j].rk * Y[j].rk;
//         sumY += Y[j].rk;
//         sumY2 += Y[j].rk * Y[j].rk;
//     }
//     //Pearson
//     correlation= ((n*sumXY) - (sumX*sumY) ) / (sqrt( (n*sumX2 - sumX*sumX) * (n*sumY2 - sumY*sumY) ));
//     return correlation;
// }

// void IC() 
// {
//     int Tpe = (int)p*pe;
//     std::vector<std::vector<std::pair<int, double> > > listaArestas(Tpe, std::vector<std::pair<int, double> >());
// 	// create weighted (pearson correlation) graph
// 	int entrouAresta = 0;
// 	double pearson = 0.0;
// 	for (int i = 0; i < Tpe - 1; i++) {
// 		for (int j = i + 1; j < Tpe; j++)
// 		{
// 			pearson = PearsonCorrelation(Pop[i].vec, Pop[j].vec);
// 			if (pearson > 0.6) {
// 				entrouAresta++;
// 				listaArestas[i].push_back(std::make_pair(j, pearson));
// 				listaArestas[j].push_back(std::make_pair(i, pearson));
// 			}
// 		}
// 	}
// 	// apply clustering method
// 	LP(listaArestas);
// 	PromisingLP();
//     listaArestas.clear();
// }

// void LP(std::vector<std::vector<std::pair<int, double> > > listaArestas)
// {
    // int nk = listaArestas.size();
	// // Create vector with visit order
	// std::vector<int> ordemVisita(nk);
	// iota(ordemVisita.begin(), ordemVisita.end(), 0);
	// // initialize each node with its own label
	// for (int i = 0; i < nk; i++)
	// 	Pop[i].label = i;
	// int iteracao = 1;
	// int labelVizinho, melhorLabel;
	// double melhorPeso;
	// std::map<int, double> totalLabels;
	// std::map<int, double>::iterator it;
	// int movimentos = 1;
	// while (movimentos) {
	// 	movimentos = 0;
	// 	random_shuffle(ordemVisita.begin(), ordemVisita.end());
	// 	for (auto idVertice : ordemVisita) {
	// 		// Calculate the weigth of the labels
	// 		totalLabels.clear();
	// 		for (auto idVizinho : listaArestas[idVertice]) {
	// 			labelVizinho = Pop[idVizinho.first].label;
	// 			it = totalLabels.find(labelVizinho);
	// 			if (it != totalLabels.end()) {
	// 				it->second += idVizinho.second;
	// 			}
	// 			else {
	// 				totalLabels[labelVizinho] = idVizinho.second;
	// 			}
	// 		}
	// 		// Best label is itself initially
	// 		melhorLabel = Pop[idVertice].label;
	// 		melhorPeso = std::numeric_limits<double>::min();
	// 		for (auto totais : totalLabels) {
	// 			if (totais.second > melhorPeso) {
	// 				melhorLabel = totais.first;
	// 				melhorPeso = totais.second;
	// 			}
	// 		}
	// 		if (melhorLabel != Pop[idVertice].label) {
	// 			Pop[idVertice].label = melhorLabel;
	// 			movimentos = 1;
	// 		}
	// 	}
	// 	iteracao++;
	// }
    // ordemVisita.clear();
// }

// void PromisingLP()
// {
    // int Tpe = (int)p*pe;
    // std::vector <int> grupos;
	// int tamanhoGrupos = 0;
	// // initialize promisings solutions
	// for (int i = 0; i < Tpe; i++)
	// 	Pop[i].promising = 0;
	// // save labels defined by LP in groups
	// int achei;
    // for (int i = 0; i < Tpe; i++)
	// {
	// 	achei = 0;
	// 	for (unsigned int j = 0; j < grupos.size(); j++)
	// 	{
	// 		if (Pop[i].label == grupos[j])
	// 			achei = 1;
	// 	}
	// 	if (achei == 0)
	// 	{
	// 		tamanhoGrupos++;
	// 		grupos.push_back(Pop[i].label);
	// 	}
	// }
	// // find the best solution in the group (with flag = 0)
	// for (unsigned int j = 0; j < grupos.size(); j++)
	// {
	// 	double menorFO = INFINITY;
	// 	int localMenor = -1;
	// 	int local = -1;
	// 	for (int i = 0; i < Tpe; i++)
	// 	{
	// 		if (Pop[i].label == grupos[j])
	// 		{
	// 			// find the best solution of the group
	// 			if (local == -1)
	// 				local = i;
	// 			// we not apply local search in this solution yet
    //             if (Pop[i].ofv < menorFO && Pop[i].flag == 0) 
	// 			{
	// 				menorFO = Pop[i].ofv;
	// 				localMenor = i;
	// 			}
	// 		}
	// 	}
	// 	if (localMenor == -1)
	// 		localMenor = local;
	// 	if (Pop[localMenor].label != -1)
	// 		Pop[localMenor].promising = 1;
	// }
// }

void FreeMemory()
{
    //methods
    
}

double randomico(double min, double max)
{
	return ((double)(rand()%10000)/10000.0)*(max-min)+min;
}

int irandomico(int min, int max)
{
	return (int)randomico(0,max-min+1.0) + min;
}