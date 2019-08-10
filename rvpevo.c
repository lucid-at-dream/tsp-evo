#include "rvpevo.h"

#include "conf.h"

void magic(tspcfg *cfg) {

    // Generate an initial population: A set of disjoint sets of nodes
    //  - Put a constraint on the maximum number of nodes that a vehicle can visit
    //  - What is the most adequate representation??

    // Fitness evaluation:
    //  - For each vehicle in a solution, compute the TSP using the set of nodes assigned to that vehicle.
    //  - The fitness of the solution is the sum of the solutions to each tsp problem.

    // Parent selection

    // Crossover algorithm:
    //  - ???

    // Mutation operators:
    //  - Reassign a node to another vehicle
    //  - More ?
}

