#!/usr/bin/python

from os import makedirs
from random import randint, random, choice, sample
from os import system
from numpy import mean
from time import time

prob_size = 100
inputs_folder = 'in'
binary = "./t"
ninputs = 4
run_count = 2

popsize = 100
npops = 10
ngens = 1000
tournament_size = 5
maxgrad0count = 100

def generate_input(N, file):
    with open(file, 'w') as f:
        f.write("%d\n" % N)
        for i in range(N):
            for j in range(N):
                if (i == j):
                    f.write("%d %d %d\n" % (i,j,0))
                else:
                    f.write("%d %d %.6lf\n" % (i,j,random()))

def generate_individual():
    ind = {
        "swap_mutation_rate": random(),
        "inversion_mutation_rate": random(),
        "individual_replacement_rate": random(),
        "elitesize": randint(1, min(20, popsize)),
        "crossover_rate": random(),
        "population_migration_rate": random(),
        "individual_migration_rate": random(),
        "#rundir": "runs/ind_%d" % randint(1, 100000000000000)
    }
    write_individual(ind)
    return ind

def crossover(ind1, ind2):
    newind = {}
    for k in ind1:
        newind[k] = ind1[k] if random() < 0.5 else ind2[k]
    newind["#rundir"] = "runs/ind_%d" % randint(1, 100000000000000)

    if "#fitness" in newind:
        newind.pop("#fitness")

    write_individual(newind)
    return newind

def mutate_individual(ind):
    newind = {}
    for k in ind:
        newind[k] = ind[k]

    param = choice([i for i in ind])
    while type(newind[param]) == str:
        param = choice([i for i in ind])

    if param in ['elitesize']:
        newind[param] += randint(-2, 2)
        newind[param] = min(newind[param], popsize)
        newind[param] = max(newind[param], 0)
    else:
        newind[param] += (random() -0.5) / 10
        newind[param] = min(newind[param], 1.0)
        newind[param] = max(newind[param], 0)
    newind["#rundir"] = "runs/ind_%d" % randint(1, 100000000000000)
    if "#fitness" in newind:
        newind.pop("#fitness")
    write_individual(newind)
    return newind

def write_individual(ind):
    makedirs(ind["#rundir"])
    with open("%s/config" % ind["#rundir"], 'w') as f:
        f.write("popsize %d\n" % popsize)
        f.write("npops %d\n" % npops)
        f.write("ngens %d\n" % ngens)
        f.write("tournament_size %d\n" % tournament_size)
        f.write("maxgrad0count %d\n" % maxgrad0count)
        f.write("\n".join(["{0} {1}".format(i, ind[i]) for i in ind]))

def evaluate_ind(ind):
    rundir = ind["#rundir"]
    scores = []
    for i in range(ninputs):
        for j in range(run_count):
            system("cat %s/%d_%d | %s -c %s/config | grep Best | awk '{print $2}' >> %s/out_%d" % (inputs_folder, prob_size, i, binary, rundir, rundir, i))

        with open("%s/out_%d" % (rundir, i), 'r') as f:
            score = [float(i) for i in f.read().split("\n") if i != ""]
            scores += [mean(score)]

    ind['#fitness'] = mean(scores)

def tournament(guys):
    return sorted(guys, key=lambda x: x['#fitness'])[:2]

# Generate input files to feed the individuals
makedirs(inputs_folder, exist_ok=True)
count=0
while count < ninputs:
    generate_input(prob_size, "%s/%d_%d" % (inputs_folder, prob_size, count))
    count += 1

# Parameterize this algorithm
_popsize = 50
_ngens = 1000
_elite = 2
_crossover_rate = 0.7
_tournament_size = 3
_mutation_rate = 0.3

# Generate initial population
start = time()

pop = [generate_individual() for i in range(_popsize)]
for ind in pop:
    evaluate_ind(ind)
pop = sorted(pop, key=lambda x: x['#fitness'])

stop = time()
print("Generation took %.2lf seconds" % (stop - start))

for gen in range(_ngens):

    start = time()

    print("[%.3d] Current best solution is %s with fitness %lf followed by %s" % (
        gen,
        pop[0]['#rundir'],
        pop[0]['#fitness'],
        str([i["#fitness"] for i in pop[1:6]])
    ))

    # Save the elite
    nextpop = pop[:_elite]

    # Generate offspring
    while len(nextpop) < _popsize:

        # Recombination
        if random() < _crossover_rate:
            parents = tournament(sample(pop, _tournament_size - 1) + [pop[len(nextpop) - 1]])
            nextpop += [crossover(*parents)]
        else:
            nextpop += [pop[len(nextpop) - 1]]

        # Mutation
        if random() < _mutation_rate:
            nextpop[-1] = mutate_individual(nextpop[-1])

    # Evalutation
    for ind in nextpop:
        if "#fitness" not in ind:
            evaluate_ind(ind)
    
    pop = nextpop
    nextpop = None

    pop = sorted(pop, key=lambda x: x['#fitness'])

    stop = time()
    print("Generation took %.2lf seconds" % (stop - start))




