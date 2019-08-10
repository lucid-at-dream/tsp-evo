#pragma once

#include "conf.h"

typedef struct _individual {
    unsigned short *perm;
    double fitness;
} individual;

individual tspevo(tspcfg *cfg);
void printIndividual(individual *ind, int indsize);