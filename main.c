#include <sys/time.h>
#include "tspevo.h"

#define MAXN 1001
#define VE
#define MULTI_THREAD 4

int main(int argc, char **argv){

    // Read configuration
    tspcfg config = argParse(argc, argv);

    //read input
    int i,j, N;
    double cost;
    scanf("%d", &N);

    int M = 50; // The number of vehicles in the fleet.

    double **costs = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        costs[i] = (double *)malloc(N * sizeof(double));
    }

    while (scanf("%d %d %lf", &i, &j, &cost ) != EOF) {
        costs[i][j] = cost;
    }

    // problem configuration
    config.indsize = N;
    config.fleetsize = M;
    config.costs = costs;

    // Measure time for statistics
    struct timeval tval_start, tval_finish, tval_whole_result;
    gettimeofday(&tval_start, NULL);
    
    individual best = tspevo(&config); // call tsp evo

    gettimeofday(&tval_finish, NULL);
    timersub(&tval_finish, &tval_start, &tval_whole_result);
    printf("Took: %ld.%06ld seconds\n", (long int)tval_whole_result.tv_sec, (long int)tval_whole_result.tv_usec);

    printf("Best: ");
    printIndividual(&best, N);
    free(best.perm);

    return 0;
}
