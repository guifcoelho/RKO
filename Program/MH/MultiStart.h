#ifndef _MultiStart_H
#define _MultiStart_H

void MultiStart(int method)
{
    //Multi Start
    static TSol s;                              // current solution
    static TSol sBest;                          // best solution

    static int IterT;                           // current iteration

    float currentTime = 0;                      // computational time of the search process
    IterT = 0;

    struct timespec TstartMH, TendMH;           // computational time (unix systems)
    clock_gettime(CLOCK_MONOTONIC, &TstartMH);
    clock_gettime(CLOCK_MONOTONIC, &TendMH);

    CreateInitialSolutions(sBest); 
    Decoder(sBest); 
    
    // run the search process until stop criterion
    while(currentTime < MAXTIME)
    {
        if (stop_execution.load()) return;      
        
        // Create the initial solution with random keys 
		CreateInitialSolutions(s); 
		Decoder(s); 

		if (s.ofv < sBest.ofv)
		{
			sBest = s;
		}
        
        // update the best solution found
        if (sBest.ofv < bestSolution.ofv) UpdatePoolSolutions(sBest, method);

        IterT++;

        // terminate the evolutionary process in MAXTIME
        clock_gettime(CLOCK_MONOTONIC, &TendMH);
        currentTime = (TendMH.tv_sec - TstartMH.tv_sec) + (TendMH.tv_nsec - TstartMH.tv_nsec) / 1e9;
    }
}

#endif