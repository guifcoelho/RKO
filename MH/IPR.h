#ifndef _IPR_H
#define _IPR_H

// Variables declared in main.cpp
extern int n;                               // size of solution

/************************************************************************************
*** Metodo: PR()                                                                  ***
*** Funcao: Aplica a metodo Path Relinking                                        ***
*************************************************************************************/
/*static TSol IPR_All(TSol atual, TSol guia)
{
    TSol melhorCaminho = atual;     //melhor solucao obtida pelo PR
    TSol melhorIteracao = atual;    //melhor solucao em cada iteracao
    TSol sCorrente = atual;         //solucao corrente em cada iteracao
    TSol sViz = atual;              //solucao vizinha em cada iteracao

    //calcula a diferenca entre as solucoes
    int dist=0;
    std:: vector<int> diferenca(n,0);
    for (int i=0; i<n; i++)
    {
    	if (atual.vec[i].rk != guia.vec[i].rk){
            diferenca[i] = 1;
            dist++;
        }
    }

    //se houver diferenca entre as solucoes 'zerar' a melhor solucao do caminho
    if (dist > 0)
        melhorCaminho.ofv = INFINITO;

    //inicializa o processo de busca do PR
    int melhorindice=-1;
    int limiteCaminho = (int)(dist*0.5); //pesquisar somente em 50% do caminho
    while(dist > limiteCaminho)
    {
        melhorindice = -1;
        melhorIteracao.ofv = INFINITO;

        //examinar todas as trocas possiveis em uma iteracao do PR
        for(int i=0; i<n; i++)
        {
            if(diferenca[i] == 1)
            {
                //gerar o vizinho da solucao corrente
                sViz = sCorrente;
                sViz.vec[i].rk = guia.vec[i].rk;
                Decoder(sViz);

                //verficar se eh a melhor solucao da iteracao
                if(sViz.ofv < melhorIteracao.ofv){
                    melhorIteracao = sViz;
                    melhorindice = i;
                }
            }
        }

        //verificar se eh a melhor solucao do caminho
        if (melhorIteracao.ofv < melhorCaminho.ofv)
           melhorCaminho = melhorIteracao;

        //continuar a busca a partir da melhor solucao da iteracao
        sCorrente = melhorIteracao;

        //apagar o indice do cliente trocado nesta iteracao
        diferenca[melhorindice]=0;

        //diminuir o n�mero de movimentos poss�veis
        dist--;
    }

    return melhorCaminho;
}*/



/************************************************************************************
*** Metodo: IPR()                                                                 ***
*** Funcao: Aplica a metodo Implicit Path Relinking                               ***
*************************************************************************************/
static TSol IPR(TSol atual, TSol guia, double factor, int direction)
{
    TSol melhorCaminho = atual;     //melhor solucao obtida pelo PR
    TSol melhorIteracao = atual;    //melhor solucao em cada iteracao
    TSol sCorrente = atual;         //solucao corrente em cada iteracao
    TSol sViz = atual;              //solucao vizinha em cada iteracao

    //calcula a diferenca entre as solucoes
    // int dist=0;
    // std:: vector<int> diferenca(n,0);
    // for (int i=0; i<n; i++)
    // {
    // 	if (atual.vec[i].rk != guia.vec[i].rk){
    //         diferenca[i] = 1;
    //         dist++;
    //     }
    // }

    //  //se houver diferenca entre as solucoes 'zerar' a melhor solucao do caminho
    // if (dist > 0)
    //     melhorCaminho.ofv = INFINITO;

    double correlation = 0;
    // double sumXY = 0;
    // double sumX2 = 0;
    // double sumY2 = 0;
    // double sumX = 0;
    // double sumY = 0;

    // for(int j=0; j<n; j++)
    // {
    //     sumX += atual.vec[j].rk;
    //     sumX2 += atual.vec[j].rk * atual.vec[j].rk;
    //     sumXY += atual.vec[j].rk * guia.vec[j].rk;
    //     sumY += guia.vec[j].rk;
    //     sumY2 += guia.vec[j].rk * guia.vec[j].rk;
    // }

    // // Pearson
    // correlation= ((n*sumXY) - (sumX*sumY) ) / (sqrt( (n*sumX2 - sumX*sumX) * (n*sumY2 - sumY*sumY) ));

    // se houver diferenca entre as solucoes 'zerar' a melhor solucao do caminho e realizar a busca
    if (correlation < 0.8)
    {
        melhorCaminho.ofv = INFINITY;

        // inicializa o processo de busca do PR
        int melhorbloco=-1;

        int blockSize = 2;
        int numBlock = n/blockSize;
        std:: vector<int> fixedBlock(numBlock,0);

        int numIteracao = 0;
        // while(numIteracao < numBlock-1)
        while(numIteracao < numBlock * factor)
        {
            numIteracao++;
            melhorbloco=-1;
            melhorIteracao.ofv = INFINITY;

            // examinar todos os blocos possiveis em uma iteracao do PR
            for(int i=0; i<numBlock; i++)
            {
                if (fixedBlock[i] == 0){
                    int initialBlock = i*blockSize;
                    int finalBlock = initialBlock + blockSize;

                    // gerar o vizinho da solucao corrente
                    sViz = sCorrente;
                    for (int k=initialBlock; k<finalBlock && k<n; k++)
                    {
                        // internal PR
                        if (direction == 1)
                            sViz.rk[k] = guia.rk[k];
                        
                        // external PR
                        else
                        if (direction == -1)
                            sViz.rk[k] = 1.0 - guia.rk[k];

                    }
                    Decoder(sViz);

                    // verficar se eh a melhor solucao da iteracao
                    if(sViz.ofv < melhorIteracao.ofv){
                        melhorIteracao = sViz;
                        melhorbloco = i;
                    }
                }
            }

            // verificar se eh a melhor solucao do caminho
            if (melhorIteracao.ofv < melhorCaminho.ofv)
                melhorCaminho = melhorIteracao;

            // continuar a busca a partir da melhor solucao da iteracao
            sCorrente = melhorIteracao;

            // apagar o indice do bloco trocado nesta iteracao
            fixedBlock[melhorbloco] = 1;
        }
    }

    return melhorCaminho;
}

#endif