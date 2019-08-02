#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAXN 1001
#define V

// problem configuration
double costs[MAXN][MAXN];

// Data structures
typedef struct _individual{
    int *perm;
    double fitness;
} individual;

individual best;

// Function declarations
void magic(int npops, int popsize, int indsize);
individual *newIndividual(individual *new, int indsize);
void evalPopFitness(individual *population, int popsize, int indsize);
double evalIndFitness(individual *ind, int indsize);
void printIndividual(individual *ind, int indsize);

int main(){

    //read input
    int i,j, N;
    double cost;
    scanf("%d", &N);
    while( scanf("%d %d %lf", &i, &j, &cost ) != EOF ) {
        costs[i][j] = cost;
    }

    //parameters
    int popsize = 250;
    int npops = 50;

    // call tsp evo
    magic(npops, popsize, N);

    return 0;
}

/**
 * Evolutionary meta-heuristic algorithm for evolving solutions to the TSP problem.
 */
void magic(int npops, int popsize, int indsize) {

    srand(clock());

    clock_t start, finish;
    start = clock();

    individual *populations[npops];

    //generate initial population
    for (int p = 0; p < npops; p++) {
        individual *population = (individual *)malloc(popsize * sizeof(individual));

        for (int i = 0; i < popsize; i++) {
            newIndividual(population + i, indsize);
        }
        
        evalPopFitness(population, popsize, indsize);
        populations[p] = population;
    }

    best = populations[0][0];
    printIndividual(&best, indsize);

    for (int p = 0; p < npops; p++) {
        for (int i = 0; i < popsize; i++) {
            free(populations[p][i].perm);
        }
        free(populations[p]);
    }

    finish = clock();
    printf("Took %lf seconds\n", (double)((finish - start)/(double)CLOCKS_PER_SEC) );
}

/**
 * Generates a new random individual with the given size.
 * 
 * @param new The individual to populate with the generated configuration.
 * @param indsize The size of the individual.
 */
individual *newIndividual(individual *new, int indsize) {
    
    new->fitness = 0;
    new->perm = (int *)malloc(indsize * sizeof(int));

    for (int i = 0; i < indsize; i++) {
        new->perm[i] = i;
    }

    // Make 2*N swaps for randomization
    int tmp, i1, i2;
    for (int i = 0; i < 2 * indsize; i++) {
        i1 = rand() % indsize;
        i2 = rand() % indsize;
        tmp = new->perm[i1];
        new->perm[i1] = new->perm[i2];
        new->perm[i2] = tmp;
    }

    return new;
}

/**
 * Evaluates the fitness of all the individuals in a population.
 * 
 * @param population The population to evaluate
 * @param popsize The size of the population to evaluate
 * @param indsize The size of the individual.
 */
void evalPopFitness(individual *population, int popsize, int indsize) {
    int i;
    for (i = 0; i < popsize; i++) {
        population[i].fitness = evalIndFitness(&population[i], indsize);
    }
}

/**
 * Evaluates the fitness of an indiviual.
 * 
 * @param ind The individual to evaluate.
 * @param indsize The size of the individual.
 */
double evalIndFitness(individual *ind, int indsize) {
    int i=0;
    double fitness = 0;

    for (i = 0; i < indsize; i++) {
        if (i < indsize - 1)
            fitness += costs[ind->perm[i]][ind->perm[i+1]];
        else
            fitness += costs[ind->perm[i]][ind->perm[0]];
    }
    return fitness;
}

/**
 * Prints an individual's characteristics to stdout.
 * 
 * @param ind     The individual to display.
 * @param indsize The size of the individual.
 */
void printIndividual(individual *ind, int indsize) {

    printf("%.4lf - ", ind->fitness);

    int first = 0;
    for (int i = 0; i < indsize; i++) {
        if (ind->perm[i] < ind->perm[first])
            first = i;
    }

    for (int i = first; i < first + indsize; i++) {
        printf("%d ", ind->perm[i % indsize]);
    }
    printf("\n");
}
