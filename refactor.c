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

void magic(
        int npops, int popsize, int indsize, int ngens, 
        double swap_mutation_rate, double inversion_mutation_rate,
        int elitesize
);

individual *newIndividual(individual *new, int indsize);
void evalPopFitness(individual *population, int popsize, int indsize);
double evalIndFitness(individual *ind, int indsize);

void evolvePopulation(
    individual *population, int popsize, int indsize, 
    double swap_mutation_rate, double inversion_mutation_rate, 
    int elitesize
);

void swapMutation(individual *ind, int indsize);
void subsequenceInversionMutation(individual *ind, int indsize);
int fit_cmp(const void *ls, const void *rs);

void printIndividual(individual *ind, int indsize);


int main(){

    //read input
    int i,j, N;
    double cost;
    scanf("%d", &N);
    while( scanf("%d %d %lf", &i, &j, &cost ) != EOF ) {
        costs[i][j] = cost;
    }


    // ===== Parameterization =====

    // How much machinery for heavy lifting?
    int popsize = 250;
    int npops = 50;
    int ngens = 1000;

    // How much randomness?
    double swap_mutation_rate = 0.3;
    double inversion_mutation_rate = 0.3;

    // How much not randomness?
    int elitesize = 10;

    // ===== Parameterization =====

    // call tsp evo
    magic(npops, popsize, N, ngens, swap_mutation_rate, inversion_mutation_rate, elitesize);

    return 0;
}

/**
 * Evolutionary meta-heuristic algorithm for evolving solutions to the TSP problem.
 */
void magic(
        int npops, int popsize, int indsize, int ngens,
        double swap_mutation_rate, double inversion_mutation_rate,
        int elitesize
    ) {

    // Keep a reference of the best individual ever.
    individual best;

    // Get a good randomness for the algorithm
    srand(clock());

    // Measure time for statistics
    clock_t start, finish;
    start = clock();

    //generate initial populations
    individual **populations = (individual **)malloc(npops * sizeof(individual *));
    for (int p = 0; p < npops; p++) {
        individual *population = (individual *)malloc(popsize * sizeof(individual));
        for (int i = 0; i < popsize; i++) {
            newIndividual(population + i, indsize);
        }
        evalPopFitness(population, popsize, indsize);
        populations[p] = population;
    }
    best = populations[0][0];

    // generations loop
    for (int gen = 1; gen <= ngens; gen++) {

        // Find the best individual among all populations for reference.
        for(int p = 1; p < npops; p++){
            if(populations[p][0].fitness < best.fitness){
                best = populations[p][0];
            }
        }

        // Evolve the populations
        for (int p = 0; p < npops; p++) {
            evolvePopulation(populations[p], popsize, indsize, swap_mutation_rate, inversion_mutation_rate, elitesize);
        }

#ifdef V
        printf("Best: ");
        printIndividual(populations[0] + 5, indsize);

        printf("Best: ");
        printIndividual(&best, indsize);
#endif
    }

    // Free all the allocated memory for easier mem leak validation.
    for (int p = 0; p < npops; p++) {
        for (int i = 0; i < popsize; i++) {
            free(populations[p][i].perm);
        }
        free(populations[p]);
    }
    free(populations);

    // Time statistics
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
        population[i].fitness = evalIndFitness(&(population[i]), indsize);
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
 * Performs one generation of evolution on a given population.
 */
void evolvePopulation(
        individual *population, int popsize, int indsize, 
        double swap_mutation_rate, double inversion_mutation_rate,
        int elitesize
    ) {

    // Mutate the population
    for (int i = elitesize; i < popsize; i++) {
        
        if ((rand() % 10000)/10000.0 <= swap_mutation_rate) {
            swapMutation(&(population[i]), indsize);
        }

        if ((rand() % 10000)/10000.0 <= inversion_mutation_rate) {
            subsequenceInversionMutation(&(population[i]), indsize);
        }
    }

    evalPopFitness(population, popsize, indsize);

    // Sort individuals by fitness
    qsort(population, popsize, sizeof(individual), fit_cmp);
}

/**
 * Mutates an individual by swaping two of its positions.
 * 
 * @param ind     The individual to mutate.
 * @param indsize The size of the individual.
 */
void swapMutation(individual *ind, int indsize) {
    int i1 = rand() % indsize;
    
    int i2 = rand() % indsize;
    while (i2 == i1) {
        i2 = rand() % indsize;
    }

    int tmp = ind->perm[i1];
    ind->perm[i1] = ind->perm[i2];
    ind->perm[i2] = tmp;
}

/**
 * Mutates an individual by inverting a subsection of the permutation.
 */
void subsequenceInversionMutation(individual *ind, int indsize) {
    int range = rand() % indsize;
    int start = rand() % (indsize - range);

    int tmp, i, count = 0;
    for (i = start; i < start + ceil(range / 2.0); i++, count++) {
        int i1 = i;
        int i2 = start + range - count;
        
        tmp = ind->perm[i1];
        ind->perm[i1] = ind->perm[i2];
        ind->perm[i2] = tmp;
    }
}

/**
 * Compares the fitness of two individuals.
 */
int fit_cmp(const void *ls, const void *rs) {
    individual *a = (individual *)ls;
    individual *b = (individual *)rs;

    if (a->fitness < b->fitness)
        return -1;
    else if (a->fitness > b->fitness)
        return 1;
    else
        return 0;
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
