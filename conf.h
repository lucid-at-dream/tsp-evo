#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>

#ifdef MULTI_THREAD
#include "worker_pool.h"
#endif

typedef struct _tspcfg {
    // What's the size of the problem?
    int indsize;

    // How much machinery for heavy lifting?
    int popsize;
    int npops;
    int ngens;

    // How much iterations without improvements before throwing the towel?
    int maxgrad0count;

    // How much randomness?
    double swap_mutation_rate;
    double inversion_mutation_rate;
    double individual_replacement_rate;

    // How much not randomness?
    int elitesize;
    int tournament_size;

    // How much recombination of the individuals?
    double crossover_rate;
    double population_migration_rate;
    double individual_migration_rate;

    // Other, helper structures
#ifdef MULTI_THREAD
    worker_pool *thread_pool;
#endif
} tspcfg;


tspcfg argParse(int argc, char **argv);
