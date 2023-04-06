#include "Problem.h"

// Sort TSol by random-keys
bool sortByRk(const TVecSol &lhs, const TVecSol &rhs) { return lhs.rk < rhs.rk; }

void ReadData(char nameTable[])
{ 
    char name[200] = "../Instances/";
    strcat(name,nameTable);

    FILE *arq;
    arq = fopen(name,"r");

    if (arq == NULL)
    {
        printf("\nERROR: File (%s) not found!\n",name);
        getchar();
        exit(1);
    }

    // => read data
    fscanf(arq, "%d", &nPontos);
    nHubs = 3;
    alfaH = 0.2;

    n = nPontos + (nPontos - nHubs) + ((nHubs * (nHubs-1))/2);

    // imprimir
    // printf("\nPontos: %d \nHubs: %d \nAlfa: %.1lf \nn: %d", nPontos, nHubs, alfaH, n);
    

    //  euclidean distance
    dist.clear();
    dist.resize(nPontos, std::vector<double>(nPontos));

    flow.clear();
    flow.resize(nPontos, std::vector<double>(nPontos));


    // read node informations
    while (!feof(arq))
    {
        int i, j;
        double d = 0.0, f = 0.0;
    	fscanf(arq, "%d %d %lf %lf", &i, &j, &f, &d);

        dist[i][j] = d/1000.0;
        flow[i][j] = f;
    }
    fclose(arq);

    // imprimir
    // for (int i=0; i<nPontos; i++)
    // {
    //     printf("\n");
    //     for (int j = 0; j < nPontos; j++)
    //     {
    //         printf("%.2lf\t",dist[i][j]);
    //     }
    // }

    // for (int i=0; i<nPontos; i++)
    // {
    //     printf("\n");
    //     for (int j = 0; j < nPontos; j++)
    //     {
    //         printf("%.2lf\t",flow[i][j]);
    //     }
    // }

    // getchar();
}

void Decoder(TSol &s)
{
    // copy the random-key sequence of current solution 
    TSol temp = s;

    // define decoder function based in the random-key of position n+1
    int dec = floor(s.vec[n].rk*numDecoders)+1;
    switch (dec)
    {
        case 1: 
            Dec1(s);
            break;

        default:
            break;
    }

    // return initial random-key sequence and maintain the solution sequence
    for (int i=0; i<n; i++){
        s.vec[i].rk = temp.vec[i].rk;
    }
}

void LocalSearch(TSol &s)
{
    // ***** we use a Random Variable Neighborhood Descent (RVND) as local search ****

    // current neighborhood
	int k = 1;

    // predefined number of neighborhood moves
    std::vector <int> NSL;
    std::vector <int> NSLAux;
    
    for (int i=1; i<=numLS; i++)
    {
        NSL.push_back(i);
        NSLAux.push_back(i);
    }

	while (!NSL.empty())
	{
        // current objective function
        double foCurrent = s.ofv;

        // randomly choose a neighborhood
        int pos = rand() % NSL.size();
        k = NSL[pos];

        switch (k)
        {
            case 1: 
                LS1(s); 
                break;

            default:
                break;
        }

        // we better the current solution
        if (s.ofv < foCurrent)
        {
            // refresh NSL
            NSL.clear();
            NSL = NSLAux;
        }
        else
        {
            // Remove N(n) from NSL
            NSL.erase(NSL.begin()+pos);
        }
	} //end while
}

double CalculateFitness(TSol s)
{
    // calculate objective function
    s.ofv = 0;
    for (int i=0; i<n; i++){
        s.ofv += dist[s.vec[i%n].sol][s.vec[(i+1)%n].sol];
    }

    return s.ofv;
}

void Dec1(TSol &s) // sort
{
    // create a initial solution of the problem
    s.ofv = 0;
    for (int j = 0; j < n+1; j++)
	{
        if (j < nPontos){
		    s.vec[j].sol = j;                    // parte 1
        }
        else{
            if (j < (nPontos + nPontos - nHubs)) // parte 2 
                s.vec[j].sol = 0;
            
            else 
                if (j < n)                       // parte 3
                    s.vec[j].sol = j - (nPontos + nPontos - nHubs);
                else
                    s.vec[j].sol = -1;
        }
	}

    // imprimir
    // printf("\ns: ");
    // for (int j = 0; j < n+1; j++)
	// {
    //     printf("%d (%.2lf) ",s.vec[j].sol, s.vec[j].rk);
    // }
    

    // inicializar as estruturas
    a.clear();
    a.resize(nPontos,0);

    h.clear();
    h.resize(nHubs,0);

    hP.clear();
    hP.resize(nPontos,-1);


    Tp.clear();
    Tp.resize(nHubs, std::vector<int>(nHubs));

    cT.clear();
    cT.resize(nHubs, std::vector<double>(nHubs));


    // sort random-key vector until nPontos
    sort(s.vec.begin(), s.vec.begin()+nPontos, sortByRk);   

    // ordenar as chaves dos arcos da arvore
    sort(s.vec.begin()+(nPontos + nPontos - nHubs), s.vec.end()-1, sortByRk);   

    // imprimir
    // printf("\ns: ");
    // for (int j = 0; j < n+1; j++)
	// {
    //     printf("%d (%.2lf) ",s.vec[j].sol, s.vec[j].rk);
    // }

    // definir os hubs
    for (int i = 0; i < nHubs; i++)
    {
        h[i] = s.vec[i].sol;
        hP[s.vec[i].sol] = i;
    }

    // printf("\nh[i]: ");
    // for (int i = 0; i < nHubs; i++)
    // {
    //     printf("%d ", h[i]);
    // }

    // printf("\nhP[i]: ");
    // for (int i = 0; i < nPontos; i++)
    // {
    //     printf("%d ", hP[i]);
    // }

    // alocar os pontos 
    int k = nPontos;
    for (int i = 0; i < nPontos; i++)
    {
        // o ponto eh um hub
        if (i < nHubs)
            a[s.vec[i].sol] = s.vec[i].sol;
        else
        {
            a[s.vec[i].sol] = h[floor(s.vec[k].rk * nHubs)];
            s.vec[k].sol = a[s.vec[i].sol];
            k++;
        }
    }

    // printf("\n\ns: ");
    // for (int j = 0; j < n+1; j++)
	// {
    //     printf("%d (%.2lf) ",s.vec[j].sol, s.vec[j].rk);
    // }

    // printf("\na[i]: ");
    // for (int i = 0; i < nPontos; i++)
    // {
    //     printf("[%d] %d, ", i, a[i]);
    // }
    
    
    // definir a arvore
    std::vector <std::vector <int> > tree;
    tree.clear();
    tree.resize((nHubs * (nHubs-1))/2, std::vector<int>(2));

    k = 0;
    for (int i=0; i<nHubs-1; i++)
    {
        for (int j=i+1; j<nHubs; j++)
        {
            tree[k][0] = h[i];
            tree[k][1] = h[j];
            k++;
        }
    }

    // adicionar arestas na arvore (substituir pela ideia do Kruskal)
    k = nPontos + nPontos - nHubs;
    for (int i=0; i<nHubs-1; i++)
    {
        int h1 = tree[s.vec[k].sol][0];
        int h2 = tree[s.vec[k].sol][1];

        Tp[hP[h1]][hP[h2]] = Tp[hP[h2]][hP[h1]] = 1;
        k++;
    }

    // printf("\n\ntree: ");
    // for (int i = 0; i < (nHubs * (nHubs-1))/2; i++)
    // {
    //     printf("\n%d \t%d", tree[i][0], tree[i][1]);
    // }

    // printf("\n\nT: ");
    // for (int i = 0; i < nHubs; i++)
    // {
    //     printf("\n");
    //     for (int j = 0; j < nHubs; j++)
    //     {
    //         printf("%d \t", Tp[i][j]);
    //     }
    // }
    

    // Calcular os custos entre hubs na árvore T
    for (int i=0; i<nHubs; i++)
    {
        for (int j=0; j<nHubs; j++)
        {
            cT[i][j] = 0;

            if (i != j)
            {
                int h1 = h[i];
                int h2 = h[j];

                // ha arco entre os hubs i e j
                if (Tp[i][j] == 1){
                    cT[i][j] = dist[h1][h2];
                }
                // nao ha arco entre os hubs i e j -> encontrar o caminho unico entre h1 e h2
                else{
                    int h3;
                    if ((i==0 && j==1) || (j==0 && i==1)) h3 = h[2];
                    if ((i==0 && j==2) || (j==0 && i==2)) h3 = h[1];
                    if ((i==1 && j==2) || (j==1 && i==2)) h3 = h[0];

                    if (i < j)
                        cT[i][j] = dist[h1][h3] + dist[h3][h2];
                    else
                        cT[i][j] = dist[h2][h3] + dist[h3][h1];
                }
            }
        }
    }

    // printf("\n\ncT: ");
    // for (int i = 0; i < nHubs; i++)
    // {
    //     printf("\n");
    //     for (int j = 0; j < nHubs; j++)
    //     {
    //         printf("%.2lf \t", cT[i][j]);
    //     }
    // }
    // getchar();


    // Calcular a função objetivo
    double custo = 0;
    for (int i=0; i<nPontos; i++){
        for (int j=0; j<nPontos; j++){
                //encontrar o custo do caminho entre i e j
                double firstMile = dist[i][a[i]];	                  // custo de i até seu hub
                double lastMile = dist[a[j]][j]; 	                  // custo do hub de j até j

                double middleMile = alfaH * cT[hP[a[i]]][hP[a[j]]];   // custo entre hubs de i e j

                custo += (firstMile + lastMile + middleMile) * flow[i][j];
        }
    }



    // s.ofv = CalculateFitness(s);
    s.ofv = custo;
}

void LS1(TSol &s) 
{
    
}

void FreeMemoryProblem()
{
    //specific problem
    dist.clear();
    flow.clear();
}
