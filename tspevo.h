#pragma once

#include "conf.h"

typedef struct _individual {
    unsigned char *perm;
    double fitness;
} individual;

individual tspevo(tspcfg *cfg);
void printIndividual(individual *ind, int indsize);