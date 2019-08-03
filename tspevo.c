/*
TODO:
 > Parameter tunning
 > Parallellization
 > Use edge flipping mutation operation
 > Check other crossover algorithms
*/

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
typedef struct _individual {
    int *perm;
    double fitness;
} individual;

typedef struct _couple {
    individual *a;
    individual *b;
} couple;

individual best;

// Function declarations

void magic(
        int npops, int popsize, int indsize, int ngens,
        int maxgrad0count,
        double swap_mutation_rate, double inversion_mutation_rate,
        int elitesize,
        int tournament_size, double crossover_rate,
        double population_migration_rate, double individual_migration_rate
);

individual **initializePopulations(int npops, int popsize, int indsize);
individual *newIndividual(individual *new, int indsize);
void evalPopFitness(individual *population, int popsize, int indsize);
double evalIndFitness(individual *ind, int indsize);

void evolvePopulation(
    individual *population, individual *nextpop, int popsize, int indsize, 
    double swap_mutation_rate, double inversion_mutation_rate, 
    int elitesize,
    int tournament_size, double crossover_rate
);

couple tournament(int tsize, individual *population, int popsize);
void crossOver(couple parents, individual *maria, individual *zezinho, int indsize);
void swapMutation(individual *ind, int indsize);
void subsequenceInversionMutation(individual *ind, int indsize);
int fit_cmp(const void *ls, const void *rs);
int fit_cmp_ptr(const void *ls, const void *rs);

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
    int ngens = 2000;

    // How much iterations without improvements before throwing the towel?
    int maxgrad0count = N*2;

    // How much randomness?
    double swap_mutation_rate = 0.15;
    double inversion_mutation_rate = 0.15;

    // How much not randomness?
    int elitesize = 3;
    int tournament_size = 5;

    // How much recombination of the individuals?
    double crossover_rate = 0.7;
    double population_migration_rate = 0.1;
    double individual_migration_rate = 0.01;

    // ===== Parameterization =====

    // call tsp evo
    magic(
        npops, popsize, N,
        ngens, maxgrad0count,
        swap_mutation_rate, inversion_mutation_rate,
        elitesize, tournament_size,
        crossover_rate, population_migration_rate, individual_migration_rate
    );

    return 0;
}

/**
 * Evolutionary meta-heuristic algorithm for evolving solutions to the TSP problem.
 */
void magic(
        int npops, int popsize, int indsize, int ngens,
        int maxgrad0count,
        double swap_mutation_rate, double inversion_mutation_rate,
        int elitesize,
        int tournament_size, double crossover_rate,
        double population_migration_rate, double individual_migration_rate
    ) {

    // Keep a reference of the best individual ever.
    individual best;

    // Get a good randomness for the algorithm
    srand(clock());

    // Measure time for statistics
    clock_t start, finish;
    start = clock();

    //generate initial populations
    individual **populations = initializePopulations(npops, popsize, indsize);
    individual **nextpops = initializePopulations(npops, popsize, indsize);

    best = populations[0][0];

    // generations loop
    int grad0count = 0;
    for (int gen = 1; gen <= ngens; gen++) {

#ifdef V
        clock_t gen_start, gen_finish;
        gen_start = clock();
#endif

        // Find the best individual among all populations for reference.
        double prev_best = best.fitness;
        for(int p = 1; p < npops; p++){
            if(populations[p][0].fitness < best.fitness){
                best = populations[p][0];
            }
        }

        if (best.fitness < prev_best) {
            grad0count = 0;
        } else {
            grad0count++;
            if (grad0count > maxgrad0count) {
                break;
            }
        }

        // Evolve the populations
        for (int p = 0; p < npops; p++) {
            evolvePopulation(
                populations[p], nextpops[p],
                popsize, indsize,
                swap_mutation_rate, inversion_mutation_rate,
                elitesize,
                tournament_size, crossover_rate
            );
        }

        //perform migrations between populations
        if ((rand() % 10000)/10000.0 <= population_migration_rate) {
            int p1 = rand() % npops;
            int p2 = rand() % npops;
            
            while(p2 == p1)
                p2 = rand() % npops;

            for (int k=0; k<popsize; k++) {
                if ((rand() % 10000)/10000.0 <= individual_migration_rate) {
                    int ind2 = rand() % popsize;
                    individual tmp = populations[p1][k];
                    populations[p1][k] = populations[p2][ind2];
                    populations[p2][ind2] = tmp;
                }
            }
        }

        // Swap nextpops with populations for double buffering (i.e. avoid memory allocation overhead).
        individual **tmp = populations;
        populations = nextpops;
        nextpops = tmp;

#ifdef V
        gen_finish = clock();
        printf("[%04d] Took %lf seconds\n", gen, (double)((gen_finish - gen_start)/(double)CLOCKS_PER_SEC) );
        printf("[%04d] Best: %lf\n", gen, best.fitness);
#endif
    }

#ifdef V
    printf("Best: ");
    printIndividual(&best, indsize);
#endif

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
 * Allocates memory for npops * popsize individuals initialized randomly.
 */
individual **initializePopulations(int npops, int popsize, int indsize) {
    individual **populations = (individual **)malloc(npops * sizeof(individual *));

    for (int p = 0; p < npops; p++) {
        individual *population = (individual *)malloc(popsize * sizeof(individual));

        for (int i = 0; i < popsize; i++) {
            newIndividual(&(population[i]), indsize);
        }

        evalPopFitness(population, popsize, indsize);
        populations[p] = population;
    }

    return populations;
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
        individual *population, individual *nextpop, int popsize, int indsize, 
        double swap_mutation_rate, double inversion_mutation_rate,
        int elitesize,
        int tournament_size, double crossover_rate
    ) {

    // Copy the elite from population[0 : elitesize] to nextpop[0 : elitesize]
    for (int i = 0; i < elitesize; i++) {
        nextpop[i].fitness = population[i].fitness;
        for (int j = 0; j < indsize; j++)
            nextpop[i].perm[j] = population[i].perm[j];
    }

    // Combine the population via sex (parent selection and crossover)
    for (int i = elitesize; i < popsize; i += 2) {
        couple parents = tournament(tournament_size, population, popsize);

        if ((rand() % 10000) / 10000.0 <= crossover_rate) {
            // Recombine the parents with a cross over strategy
            if (i + 1 < popsize) {
                crossOver(parents, &(nextpop[i]), &(nextpop[i + 1]), indsize);
            } else {
                crossOver(parents, &(nextpop[i]), NULL, indsize);
            }

        }else{
            // Copy the parents "as is"
            nextpop[i].fitness = parents.a->fitness;
            for (int j = 0; j < indsize; j++)
                nextpop[i].perm[j] = parents.a->perm[j];

            if (i + 1 < popsize) {
                nextpop[i+1].fitness = parents.a->fitness;
                for (int j = 0; j < indsize; j++)
                    nextpop[i+1].perm[j] = parents.a->perm[j];
            }
        }
    }

    // Mutate the next population
    for (int i = elitesize; i < popsize; i++) {
        
        if ((rand() % 10000)/10000.0 <= swap_mutation_rate) {
            swapMutation(&(nextpop[i]), indsize);
        }

        if ((rand() % 10000)/10000.0 <= inversion_mutation_rate) {
            subsequenceInversionMutation(&(nextpop[i]), indsize);
        }
    }

    // Re-evaluate the population's fitness and sort the individuals by fitness
    evalPopFitness(nextpop, popsize, indsize);
    qsort(nextpop, popsize, sizeof(individual), fit_cmp);
}

/**
 * Choose a couple of parents by k-way tournament selection.
 */
couple tournament(int tsize, individual *population, int popsize) {

    int used[tsize];
    memset(used, -1, sizeof(int) * tsize);

    individual *individuals[tsize];

    for (int i = 0; i < tsize; i++) {
        int idx = rand() % popsize;

        // Re-assign a random value untill its unique in this tournament
        for (int j = 0; j < tsize; j++) {
            if (used[j] == idx) {
                idx = rand() % popsize;
                j = 0;
            }
        }
        individuals[i] = &(population[idx]);
    }

    qsort(individuals, tsize, sizeof(individual *), fit_cmp_ptr);

    couple parents;
    parents.a = individuals[0];
    parents.b = individuals[1];

    return parents;
}

/**
 * Creates two new individuals based on the recombination of two parent individuals.
 * 
 * The crossover algorithm used is based on combining portions of each parent that have
 * references to the exact same positions. 
 * 
 * For example, for parents:
 * a = [5, 3, 2, 4, 1, 0]
 * b = [1, 2, 0, 5, 4, 3]
 * 
 * The following portions will be combined:
 * [5, -, -, 4, 1, -]
 * [-, 2, 0, -, -, 3]
 * and
 * [-, 3, 2, -, -, 0]
 * [1, -, -, 5, 4, -]
 */
void crossOver(couple parents, individual *maria, individual *zezinho, int indsize) {

    char used[indsize];
    memset(used, 0, sizeof(char) * indsize);

    int mothermap[indsize];

    for (int i = 0; i < indsize; i++) {
        mothermap[parents.b->perm[i]] = i;
    }

    individual *father = parents.a;
    individual *mother = parents.b;

    int count = 0;
    int idx = 0;
    int startIdx = 0;
    int cycle_count = 0;
    while (count < indsize ) {

        cycle_count++;

        for (int i = startIdx; i < indsize; i++) {
            if (!used[i]) {
                idx = i;
                startIdx = idx;
                break;
            }
        }

        // Decide whether each son will receive from one parent or the other
        // double r = (rand() % 10000) / 10000.0;
        
        while (!used[idx]) {
            if (cycle_count % 2) {
                maria->perm[idx] = mother->perm[idx];
                if (zezinho) zezinho->perm[idx] = father->perm[idx];
            } else {
                maria->perm[idx] = father->perm[idx];
                if (zezinho) zezinho->perm[idx] = mother->perm[idx];
            }

            count++;
            used[idx] = 1;
            idx = mothermap[father->perm[idx]];
        }
    }
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
 * Compares the fitness of two pointers to individuals.
 */
int fit_cmp_ptr(const void *ls, const void *rs) {
    individual *a = *(individual **)ls;
    individual *b = *(individual **)rs;

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
