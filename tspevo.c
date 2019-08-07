/*
TODO:
 > Parameter tunning -- Also make MIGRATIONS macro configurable
 > Use edge flipping mutation operation
 > Check other crossover algorithms
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "worker_pool.h"
#include <unistd.h>
#include "conf.h"

#define MAXN 1001
#define VE
#define MIGRATIONS 50

/**
 * The flag below makes the program run across multiple threads. However, the parallelization
 * makes it slower. I've added bulk jobs (i.e. 100 generations for a populations range per thread)
 * to avoid context switching and synchronization overheads, but to no avail. I've also ran tests
 * using small populations in order to rule out cache misses as the culprit, they don't seem to be
 * the culprit. This has been put on hold and we shall go on single threadedly.
 */
// #define MULTI_THREAD 4

volatile unsigned long popevo_count;

// problem configuration
double costs[MAXN][MAXN];

// Data structures
typedef struct _individual {
    unsigned char *perm;
    double fitness;
} individual;

typedef struct _couple {
    individual *a;
    individual *b;
} couple;

typedef struct _evo_job {
    individual **populations;
    individual **nextpops;
    int start;
    int end;
    tspcfg *cfg;
} evo_job;

individual best;

// Function declarations

individual magic(tspcfg *cfg);

individual **initializePopulations(int npops, int popsize, int indsize);
individual *newIndividual(individual *new, int indsize);
void shuffle(unsigned char *array, int size);
void evalPopFitness(individual *population, int popsize, int indsize);
double evalIndFitness(individual *ind, int indsize);

void popEvolutionThreadPoolWrapper(void *_job);
void evolvePopulation(individual *population, individual *nextpop, tspcfg *cfg);

couple tournament(int tsize, individual *population, int popsize);
void crossOver(couple parents, individual *maria, individual *zezinho, int indsize);
void swapMutation(individual *ind, int indsize);
void subsequenceInversionMutation(individual *ind, int indsize);
int fit_cmp(const void *ls, const void *rs);
int fit_cmp_ptr(const void *ls, const void *rs);

void free_the_world(individual **world, int npops, int popsize);
void printIndividual(individual *ind, int indsize);


int main(int argc, char **argv){

    // Read configuration
    tspcfg config = argParse(argc, argv);

    //read input
    int i,j, N;
    double cost;
    scanf("%d", &N);
    while( scanf("%d %d %lf", &i, &j, &cost ) != EOF ) {
        costs[i][j] = cost;
    }

    // Parameterization
    config.indsize = N;

    // Other helpder structures
#ifdef MULTI_THREAD
    config.thread_pool = pool_new(MULTI_THREAD, popEvolutionThreadPoolWrapper);
    pool_start(config.thread_pool);
#endif

    // call tsp evo
    individual best = magic(&config);

    // Free up the thread pool resources
#ifdef MULTI_THREAD
    pool_stop(config.thread_pool);
    pool_del(config.thread_pool);
#endif

    printf("Best: ");
    printIndividual(&best, N);
#ifdef V
    printf("Total evo calls: %lu\n", popevo_count);
#endif

    free(best.perm);

    return 0;
}

individual multi_thread_generations_loop(individual **populations, individual **nextpops, tspcfg *cfg);
individual single_thread_generations_loop(individual **populations, individual **nextpops, tspcfg *cfg);

/**
 * Evolutionary meta-heuristic algorithm for evolving solutions to the TSP problem.
 */
individual magic(tspcfg *cfg) {

    // Get a good randomness for the algorithm
    srand(clock());

    // Measure time for statistics
    struct timeval tval_start, tval_finish, tval_whole_result;
    gettimeofday(&tval_start, NULL);

    //generate initial populations
    individual **populations = initializePopulations(cfg->npops, cfg->popsize, cfg->indsize);
    individual **nextpops = initializePopulations(cfg->npops, cfg->popsize, cfg->indsize);

#ifdef MULTI_THREAD
    individual best = multi_thread_generations_loop(populations, nextpops, cfg);
#else
    individual best = single_thread_generations_loop(populations, nextpops, cfg);
#endif

    // Free all the allocated memory for easier mem leak validation.
    free_the_world(populations, cfg->npops, cfg->popsize);
    free_the_world(nextpops, cfg->npops, cfg->popsize);

    // Time statistics
    gettimeofday(&tval_finish, NULL);
    timersub(&tval_finish, &tval_start, &tval_whole_result);
    printf("Took: %ld.%06ld seconds\n", (long int)tval_whole_result.tv_sec, (long int)tval_whole_result.tv_usec);

    return best;
}

individual single_thread_generations_loop(individual **populations, individual **nextpops, tspcfg *cfg) {

    // Keep a reference of the best individual ever.
    individual best;
    best.perm = (unsigned char *)malloc(cfg->indsize * sizeof(unsigned char));
    best.fitness = populations[0][0].fitness;

    // generations loop
    int grad0count = 0;
    for (int gen = 1; gen <= cfg->ngens; gen++) {

#ifdef V
        struct timeval tval_gen_start, tval_gen_finish, tval_gen_result;
        gettimeofday(&tval_gen_start, NULL);
#endif

        // Sort the individuals in the populations by fitness
        for (int p = 0; p < cfg->npops; p++) {
            qsort(populations[p], cfg->popsize, sizeof(individual), fit_cmp);
        }

        // Find the best individual among all populations for reference.
        double prev_best = best.fitness;
        for (int p = 1; p < cfg->npops; p++) {
            if (populations[p][0].fitness < best.fitness) {
                best.fitness = populations[p][0].fitness;
                for (int j = 0; j < cfg->indsize; j++) {
                    best.perm[j] = populations[p][0].perm[j];
                }
            }
        }

        if (best.fitness < prev_best) {
            grad0count = 0;
        } else {
            grad0count++;
            if (grad0count > cfg->maxgrad0count) {
                break;
            }
        }

        // Evolve the populations
        for (int p = 0; p < cfg->npops; p++) {
            evolvePopulation(populations[p], nextpops[p], cfg);
        }

#ifdef MIGRATIONS
        // Have migrations every 100 generations
        // Migrations become more likely as the number of generations increases, proportional to the logarithm of the number of generations.
        if (gen % MIGRATIONS == 0) {
            for (int p1 = 0; p1 < cfg->npops; p1++) {

                if ((rand() % 10000)/10000.0 <= cfg->population_migration_rate * (log(gen) / log(2))) {
                    
                    int p2 = rand() % cfg->npops;
                    while(p2 == p1)
                        p2 = rand() % cfg->npops;

                    for (int k = 0; k < cfg->popsize; k++) {
                        if ((rand() % 10000)/10000.0 <= cfg->individual_migration_rate * (log(gen) / log(2))) {
                            int ind2 = rand() % cfg->popsize;
                            individual tmp = nextpops[p1][k];
                            nextpops[p1][k] = nextpops[p2][ind2];
                            nextpops[p2][ind2] = tmp;
                        }
                    }
                }
            }
        }
#endif

        // Swap nextpops with populations for double buffering (i.e. avoid memory allocation overhead).
        individual **tmp = populations;
        populations = nextpops;
        nextpops = tmp;

#ifdef V
        gettimeofday(&tval_gen_finish, NULL);
        timersub(&tval_gen_finish, &tval_gen_start, &tval_gen_result);
        printf("Took: %ld.%06ld seconds\n", (long int)tval_gen_result.tv_sec, (long int)tval_gen_result.tv_usec);
        printf("[%04d] Best: %lf\n", gen, best.fitness);
#endif
    }

    return best;
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
    new->perm = (unsigned char *)malloc(indsize * sizeof(unsigned char));

    for (int i = 0; i < indsize; i++) {
        new->perm[i] = i;
    }
    
    shuffle(new->perm, indsize);

    return new;
}

/**
 * Makes 2*N swaps on an array for randomization
 */
void shuffle(unsigned char *array, int size) {
    int tmp, i1, i2;
    for (int i = 0; i < 2 * size; i++) {
        i1 = rand() % size;
        i2 = rand() % size;
        tmp = array[i1];
        array[i1] = array[i2];
        array[i2] = tmp;
    }
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
    int i = 0;
    double fitness = 0;

    for (i = 0; i < indsize - 1; i++) {
        fitness += costs[ind->perm[i]][ind->perm[i+1]];
    }
    fitness += costs[ind->perm[indsize - 1]][ind->perm[0]];

    return fitness;
}

/**
 * Performs one generation of evolution on a given population.
 */
void evolvePopulation(individual *population, individual *nextpop, tspcfg *cfg) {

#ifdef V
    struct timeval tval_start, tval_finish, tval_result;
    gettimeofday(&tval_start, NULL);
#endif

    // Copy the elite from population[0 : elitesize] to nextpop[0 : elitesize]
    for (int i = 0; i < cfg->elitesize; i++) {
        nextpop[i].fitness = population[i].fitness;
        for (int j = 0; j < cfg->indsize; j++)
            nextpop[i].perm[j] = population[i].perm[j];
    }

    // Combine the population via sex (parent selection and crossover)
    for (int i = cfg->elitesize; i < cfg->popsize; i += 2) {
        couple parents = tournament(cfg->tournament_size, population, cfg->popsize);

        if ((rand() % 10000) / 10000.0 <= cfg->crossover_rate) {
            // Recombine the parents with a cross over strategy
            if (i + 1 < cfg->popsize) {
                crossOver(parents, &(nextpop[i]), &(nextpop[i + 1]), cfg->indsize);
            } else {
                crossOver(parents, &(nextpop[i]), NULL, cfg->indsize);
            }

        }else{
            // Copy the parents "as is"
            nextpop[i].fitness = parents.a->fitness;
            for (int j = 0; j < cfg->indsize; j++)
                nextpop[i].perm[j] = parents.a->perm[j];

            if (i + 1 < cfg->popsize) {
                nextpop[i+1].fitness = parents.a->fitness;
                for (int j = 0; j < cfg->indsize; j++)
                    nextpop[i+1].perm[j] = parents.a->perm[j];
            }
        }
    }

    // Mutate the next population
    for (int i = cfg->elitesize; i < cfg->popsize; i++) {
        
        if ((rand() % 10000)/10000.0 <= cfg->swap_mutation_rate) {
            swapMutation(&(nextpop[i]), cfg->indsize);
        }

        if ((rand() % 10000)/10000.0 <= cfg->inversion_mutation_rate) {
            subsequenceInversionMutation(&(nextpop[i]), cfg->indsize);
        }
    }

    // Replace some individuals randomly
    for (int i = cfg->elitesize; i < cfg->popsize; i++) {
        if ((rand() % 10000) / 10000.0 < cfg->individual_replacement_rate) {
            shuffle(nextpop[i].perm, cfg->indsize);
        }
    }

    // Eval the fitness of the population
    evalPopFitness(nextpop, cfg->popsize, cfg->indsize);

#ifdef V
    gettimeofday(&tval_finish, NULL);
    timersub(&tval_finish, &tval_start, &tval_result);
#endif

    popevo_count++;
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
        for (int j = 0; used[j] >= 0 && j < tsize; j++) {
            if (used[j] == idx) {
                idx = rand() % popsize;
                j = 0;
            }
        }
        individuals[i] = &(population[idx]);
    }

    double max_val = 0;
    for (int i = 1; i < tsize; i++) {
        max_val += individuals[i]->fitness;
    }

    int i1 = rand() % tsize;
    int i2 = rand() % tsize;
    while (i2 == i1) i2 = rand() % tsize;

    couple parents;
    parents.a = individuals[i1];
    parents.b = individuals[i2];

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

    unsigned char mothermap[indsize];

    for (int i = 0; i < indsize; i++) {
        mothermap[parents.b->perm[i]] = i;
    }

    individual *father = parents.a;
    individual *mother = parents.b;

    int count = 0;
    unsigned int idx = 0;
    unsigned int startIdx = 0;
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
    int range = rand() % (indsize / 2);
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
 * Liberates all the populations. Frees their allocated memory.
 */
void free_the_world(individual **world, int npops, int popsize) {
    for (int p = 0; p < npops; p++) {
        for (int i = 0; i < popsize; i++) {
            free(world[p][i].perm);
        }
        free(world[p]);
    }
    free(world);
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

#ifdef MULTI_THREAD
individual multi_thread_generations_loop(individual **populations, individual **nextpops, tspcfg *cfg) {
    
    int bulk_size = 5000;
    evo_job jobs[cfg->npops];

    int step = cfg->npops / cfg->thread_pool->num_threads / 2;
    step = step > 0 ? step : 1;

    // Keep a reference of the best individual ever.
    individual best;
    best.perm = (unsigned char *)malloc(cfg->indsize * sizeof(unsigned char));
    best.fitness = populations[0][0].fitness;

    // generations loop
    for (int gen = 1; gen <= cfg->ngens / bulk_size; gen++) {

#ifdef V
        struct timeval tval_bulk_start, tval_bulk_finish, tval_bulk_result;
        gettimeofday(&tval_bulk_start, NULL);
#endif

        // Find the best individual among all populations for reference.
        for (int p = 1; p < cfg->npops; p++) {
            if (populations[p][0].fitness < best.fitness) {
                best.fitness = populations[p][0].fitness;
                for (int j = 0; j < cfg->indsize; j++) {
                    best.perm[j] = populations[p][0].perm[j];
                }
            }
        }

        // Bulk insert the population evolution jobs to avoid thread context switching
        for (int p = 0; p < cfg->npops; p += step) {
            jobs[p].cfg = cfg;
            jobs[p].populations = populations;
            jobs[p].nextpops = nextpops;
            jobs[p].start = p;
            jobs[p].end = p + step;
            pool_add_work(cfg->thread_pool, &(jobs[p]));
        }

        pool_await_empty_queue(cfg->thread_pool);

#ifdef V
        gettimeofday(&tval_bulk_finish, NULL);
        timersub(&tval_bulk_finish, &tval_bulk_start, &tval_bulk_result);
        printf("Took: %ld.%06ld seconds\n", (long int)tval_bulk_result.tv_sec, (long int)tval_bulk_result.tv_usec);
        printf("[%04d] Best: %lf\n", gen, best.fitness);
#endif
    }

    return best;
}

void popEvolutionThreadPoolWrapper(void *_job) {

#ifdef V
    struct timeval tval_start, tval_finish, tval_result;
    gettimeofday(&tval_start, NULL);
#endif

    evo_job *job = (evo_job *)_job;
    int bulk_size = 5000;

    for (int bulk = 0; bulk < bulk_size; bulk++) {
        for (int i = job->start; i < job->end && i < job->cfg->npops; i++) {
            evolvePopulation(job->populations[i], job->nextpops[i], job->cfg);
            qsort(job->nextpops[i], job->cfg->popsize, sizeof(individual), fit_cmp);

            // Swap (double buffer)
            individual *tmp = job->populations[i];
            job->populations[i] = job->nextpops[i];
            job->nextpops[i] = tmp;
        }
    }

#ifdef V
    gettimeofday(&tval_finish, NULL);
    timersub(&tval_finish, &tval_start, &tval_result);
    printf("Bulk Job Took: %ld.%06ld seconds\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
#endif

}
#endif
