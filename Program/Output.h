/************************************************************************************
									IO Functions
*************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "Data.h"
 
void WriteSolutionScreen(TSol s, int n, float timeBest, float timeTotal, char instance[])
{
	printf("\n\n\n Instance: %s \nsol: ", instance);
	for (int i=0; i<n; i++)
		printf("%d ", s.vec[i].sol);

    printf("\nDecoder: %.2lf",s.vec[n].rk);
	printf("\nfo: %.5lf",s.ofv);
	printf("\nTotal time: %.3f",timeTotal);
	printf("\nBest time: %.3f\n\n",timeBest);
}

void WriteSolution(TSol s, int n, float timeBest, float timeTotal, char instance[])
{
	FILE *arq;
    arq = fopen("../Results/Solutions.txt","a");

	if (!arq)
	{
		printf("\n\nFile not found Solutions.txt!!!");
		getchar();
		exit(1);
	}

    fprintf(arq,"\n\nInstance: %s", instance);
	fprintf(arq,"\nSol: ");
	for (int i=0; i<n; i++)
		fprintf(arq,"%d ", s.vec[i].sol);
	fprintf(arq,"\nFO: %.5lf", s.ofv);
  	fprintf(arq,"\nBest time: %.3f",timeBest);
	fprintf(arq,"\nTotal time:%.3f \n",timeTotal);

	fclose(arq);
}

void WriteResults(double ofv, double ofvAverage, std::vector <double> ofvs, float timeBest, float timeTotal, char instance[])
{
	FILE *arq;
    arq = fopen("../Results/Results.csv","a");

	if (!arq)
	{
		printf("\n\nFile not found Results.xls!!!");
		getchar();
		exit(1);
	}

	fprintf(arq,"\n%s", instance);
    fprintf(arq,"\t%d", (int)ofvs.size());
    for (unsigned int i=0; i<ofvs.size(); i++){
        fprintf(arq,"\t%lf", ofvs[i]);
	}
	fprintf(arq,"\t%lf", ofv);
	fprintf(arq,"\t%lf", ofvAverage);
	fprintf(arq,"\t%.3f", timeBest);
	fprintf(arq,"\t%.3f", timeTotal);

	fclose(arq);
}