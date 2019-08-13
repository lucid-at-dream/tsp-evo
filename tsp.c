#include "tsp.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

typedef struct _solution {
    double cost;
    unsigned short *perm;
} solution;

typedef struct _problem {
    int N;
    unsigned short *nodes;
    double **costs;
    solution best;
} problem;

double *calc_lowerbounds(problem *cfg);
int cmp (const void *ls, const void* rs);

void solve(problem *cfg, double *lowerbounds);
void tsp(int idx, double currentcost, unsigned short *currentperm, unsigned short *visited, problem *cfg, double *lowerbounds);

void optimal_path(int N, unsigned short *nodes, double **costs) {

    problem cfg;
    cfg.N = N;
    cfg.nodes = nodes;
    cfg.costs = costs;
    cfg.best.cost = 1e20;
    cfg.best.perm = (unsigned short *)calloc(N, sizeof(unsigned short));

    double *lowerbounds = calc_lowerbounds(&cfg);

    // Solve the shortest path optimally
    solve(&cfg, lowerbounds);

    // Rearrange the original array according to the found solution
    unsigned short *tmp = (unsigned short *)malloc(N * sizeof(unsigned short));
    memcpy(tmp, nodes, N * sizeof(unsigned short));

    unsigned short *best = cfg.best.perm;

    for (int i = 0; i < N; i++) {
        nodes[i] = tmp[best[i]];
    }

    // Free allocated resources
    free(tmp);
    free(lowerbounds);
    free(cfg.best.perm);
}

double *calc_lowerbounds(problem *cfg) {
    double *mincost = malloc(cfg->N * sizeof(double));

    double mins[cfg->N];
    for (int i = 0; i < cfg->N; i++) {
        mins[i] = 1e20;
    }

    for (int i = 0; i < cfg->N; i++) {
        for (int j = 0; j < cfg->N; j++) {

            if (i == j) continue;

            int i1 = cfg->nodes[i], i2 = cfg->nodes[j];
            double cost = cfg->costs[i1][i2];

            if (cost < mins[i]) {
                mins[i] = cost;
            }
        }
    }
    qsort(mins, cfg->N, sizeof(double), cmp);

    mincost[0] = 0;
    for (int i = 1; i < cfg->N; i++) {
        mincost[i] = mincost[i-1] + mins[i-1];
    }

    return mincost;
}

int cmp (const void *ls, const void* rs){
    double a = *((double *) ls);
    double b = *((double *) rs);

    if( a < b )
        return -1;
    else if( a > b )
        return 1;
    else
        return 0;
}

void solve(problem *cfg, double *lowerbounds) {
    unsigned short *visited = (unsigned short *)calloc(cfg->N, sizeof(unsigned short));
    unsigned short *partialperm = (unsigned short *)calloc(cfg->N, sizeof(unsigned short));

    tsp(1, 0, partialperm, visited, cfg, lowerbounds);

    free(visited);
    free(partialperm);
}

void tsp(int idx, double currentcost, unsigned short *currentperm, unsigned short *visited, problem *cfg, double *lowerbounds) {

    int i;
    if (idx == cfg->N) {
        if (currentcost < cfg->best.cost) {
            cfg->best.cost = currentcost;
            memcpy(cfg->best.perm, currentperm, cfg->N * sizeof(unsigned short));
        }
        return;
    }

    if (currentcost + lowerbounds[cfg->N-idx-1] >= cfg->best.cost)
        return;

    for (i = 1; i < cfg->N; i++) {
        if (!visited[i]) {
            currentperm[idx] = i;

            int curnode = cfg->nodes[i];
            int nextnode = cfg->nodes[currentperm[idx-1]];

            visited[i] = 1;
            tsp(idx+1, currentcost + cfg->costs[curnode][nextnode], currentperm, visited, cfg, lowerbounds);
            visited[i] = 0;
        }
    }
}
