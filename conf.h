#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>

#include "worker_pool.h"

typedef struct _tspcfg {
    // What's the size of the problem?
    int indsize;
    int fleetsize;

    // How much machinery for heavy lifting?
    int nthreads;
    int popsize;
    int npops;
    int ngens;

    // How much iterations without improvements before throwing the towel?
    int maxgrad0count;

    // How much randomness?
    double swap_mutation_rate;
    double inversion_mutation_rate;
    double individual_replacement_rate;
    double optimal_mutation_rate;

    // How much not randomness?
    int elitesize;
    int tournament_size;

    // How much recombination of the individuals?
    double crossover_rate;
    int migrations;
    double population_migration_rate;
    double individual_migration_rate;

    // Other, helper structures
    worker_pool *thread_pool;
    unsigned int *rand_seeds;
    double **costs;
} tspcfg;


tspcfg argParse(int argc, char **argv);
