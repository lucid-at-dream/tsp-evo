/*
TODO:

---Increase diversity of solutions---
Problem: big population size improves results, but many are equivalent.
Solution:
    > Imigrants
    > Multiple populations
    > Parameter tunning
---Increase diversity of solutions---


---Improve performance of the program---
    > re-use memory instead of freeing it
    > paralelization?
---Improve performance of the program---


---Use a edge flipping mutation---
Enhancement: There is a well known and good mutation that uses edge flipping.

Hints: It's currently one of the best solutions and you probably can find about it
       pretty easily.
---Use a edge flipping mutation---


---Check other crossover algorithms---
Enhancement: There must be other crossover algorithms that might be more interesting.

Hints: What about finding a path from one solution to the other and choosing something
       in between? How could you do that?
---Check other crossover algorithms---
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define merda printf("merda\n");
#define cona(x) printf("cona--%d\n", x);

#define MAXN 1001

//problem variables
int N;
double costs[MAXN][MAXN];
double costs_v[MAXN * MAXN];

typedef struct _individual{
    int *perm;
    double fitness;
}individual;

individual *best = NULL;

void printIndividual(individual *ind, int length);
int fit_cmp(const void *ls, const void *rs);
double fitness(individual *ind);
individual **survivors_elitism(individual **elite, individual **offspring, int eliteSize, int popsize);
individual **copyElite(individual **population, int length, int eliteSize, int popsize);
individual **tournamentSelection( int tsize, individual **population, int popsize );
individual *tournament(int tsize, individual **population, int popsize);
void evalPopFitness(individual **population, int popsize);
void mutate2(individual *ind, int length);
void mutate1(individual *ind, int length);
individual **crossOver(individual *father, individual *mother, int length);
int indexOf(int node, individual *ind, int length);
individual *newIndividual(int length);
individual *copyIndividual(individual *old, int length);
void deleteIndividual(individual *i);
individual **evolvePopulation(individual ** population, double mutation_rate, double crossover_rate, int popsize, int tournamentSize, int eliteSize);

int int_cmp(const void *a, const void *b) {
    int _a = *(int *)a;
    int _b = *(int *)b;
    return _a > _b ? 1 : _a == _b ? 0 : -1;
}

int main(){

    //read input
    int i,j;
    double cost;
    scanf("%d", &N);
    while( scanf("%d %d %lf", &i, &j, &cost ) != EOF ) {
        costs[i][j] = cost;
    }

    //parameters
    int ngens = 2000;
    double mutation_rate = 0.2;
    double crossover_rate = 0.7;
    double migration_rate = 0.01;
    double migrationProb = 0.2;
    double replacementProb = 0;
    double replacementInc = 0.1/(double)ngens;
    double replacementProbMax = 0.15;
    int popsize = 100;
    int tournamentSize = 3;
    int eliteSize = 1;
    int npops = 20;

    int maxgrad0count = 200;

    srand(clock());

    printf("starting the algorithm\n");
    clock_t start, finish;
    start = clock();

    individual **populations[npops];

    //generate initial population
    int p;
    for(p=0; p<npops; p++){
        individual **population = (individual **)malloc(sizeof(individual *) * popsize);
        for(i=0; i<popsize; i++){
            population[i] = newIndividual(N);
        }
        evalPopFitness(population, popsize);
        populations[p] = population;
    }

    best = populations[0][0];

    double prev = 0;
    int gen, k, grad0count=0;
    individual *aux;
    for( gen = 1; gen <= ngens; gen++ ){

        clock_t gen_start, gen_finish;
        gen_start = clock();

        if(replacementProb < replacementProbMax)
            replacementProb += replacementInc;
        //printf("%lf\n", replacementProb);

        //printf("gen - %d\n", gen);
        aux = populations[0][0];
        for(p=1; p<npops; p++){
            if(populations[p][0]->fitness < aux->fitness){
                aux = populations[p][0];
            }
        }

        if(aux->fitness < best->fitness) {
            double prev_best_fitness = best->fitness;
            best = aux;
            printf("[%d] best of best: %.4lf (-%.4lf)\n", best->fitness, prev_best_fitness);
        }

        if( fabs(prev - aux->fitness) <= 0.00001 )
            grad0count+=1;
        else
            grad0count=0;
        prev = aux->fitness;
        if(grad0count >= maxgrad0count)
            break;


        //perform the generation evolution
        for(p=0; p<npops; p++){
            populations[p] = evolvePopulation(populations[p], mutation_rate, crossover_rate, popsize, tournamentSize, eliteSize);
        }

        //perform replacements
        for(p=0; p<npops; p++){
            for(k=0; k<popsize; k++){
                if((rand() % 10000)/10000.0 <= replacementProb){
                    populations[p][k] = newIndividual(N);
                    evalPopFitness(&populations[p][k], 1);
                }
            }
        }

        //perform the migrations
        if( (rand() % 10000)/10000.0 <= mutation_rate ){
            int p1 = rand() % npops,
                p2 = rand() % npops;
            while(p1 == p2)
                p2 = rand() % npops;

            for(k=0; k<popsize; k++){
                if( (rand() % 10000)/10000.0 <= migrationProb ){
                    int ind2 = rand() % popsize;
                    individual *tmp = populations[p1][k];
                    populations[p1][k] = populations[p2][ind2];
                    populations[p2][ind2] = tmp;
                }
            }
        }

        gen_finish = clock();
        printf("[%d] %lf seconds\n", gen, (double)((gen_finish - gen_start)/(double)CLOCKS_PER_SEC) );
    }

    //fetch the best individuals
    for(p=0; p<npops; p++){
        qsort(populations[p], popsize, sizeof(individual *), fit_cmp);
    }

    finish = clock();

    printf("the algorithm ran in %lf seconds\n", (double)((finish - start)/(double)CLOCKS_PER_SEC) );

    int pos[npops];
    memset(pos, 0, sizeof(int)*npops);

    int minpop = 0;
    individual *min;

    for(i=0; i<5; i++){
        min = populations[0][pos[0]];
        for(p=0; p<npops; p++){
            if(populations[p][pos[p]]->fitness < min->fitness){
                min = populations[p][pos[p]];
                minpop = p;
            }
        }
        pos[minpop]++;
        printIndividual(min, N);
    }

    printf("best of best:\n");
    printIndividual(best,N);

    return 0;
}

individual **evolvePopulation(individual ** population,
                              double mutation_rate    ,
                              double crossover_rate   ,
                              int    popsize          ,
                              int    tournamentSize   ,
                              int    eliteSize        ){

    int i;

    //make a copy of the elite, so that one can mutate everything at will.
    individual **elite = copyElite(population, N, eliteSize, popsize);

    //select parents
    individual **parents = tournamentSelection(tournamentSize, population, popsize);

    //generate offspring
    individual **offspring = (individual **)malloc(sizeof(individual *)*popsize);
    for(i=0; i<popsize-1; i+=2){
        if( (rand() % 10000)/10000.0 > crossover_rate ){
            offspring[i] = copyIndividual( parents[i], N );
            offspring[i+1] = copyIndividual( parents[i+1], N );
        }else{
            individual **sons = crossOver(parents[i], parents[i+1], N);
            evalPopFitness(sons, 2);
            offspring[i] = sons[0];
            offspring[i+1] = sons[1];
        }
    }

    //mutate offspring
    for(i=0; i<popsize; i++){
        if( (rand() % 10000)/10000.0 > mutation_rate ){
            continue;
        }else{
            mutate1(offspring[i], N);
            offspring[i]->fitness = fitness(offspring[i]);
        }
    }

    //select survivors
    individual **nextPop = survivors_elitism(elite, offspring, eliteSize, popsize);

    /*clean unselected population*/
    /*clean crossovered parents*/
    /*clean unselected offspring*/

    return nextPop;
}

double fitness(individual *ind){
    int i=0;
    double cost = 0;
    for(i=0; i<N; i++){
        if(i<N-1)
            cost += costs[ind->perm[i]][ind->perm[i+1]];
        else
            cost += costs[ind->perm[i]][ind->perm[0]];
    }
    return cost;
}

individual **survivors_elitism(individual **elite, individual **offspring, int eliteSize, int popsize){
    individual **population = (individual **)malloc(sizeof(individual *) * popsize);

    qsort(offspring, popsize, sizeof(individual *), fit_cmp);

    int i;
    for(i=0; i<eliteSize; i++){
        population[i] = elite[i];
    }

    for(i=0; i<popsize-eliteSize; i++){
        population[eliteSize + i] = offspring[i];
    }

    return population;
}

int fit_cmp(const void *ls, const void *rs){
    individual *a = *((individual **)ls);
    individual *b = *((individual **)rs);

    if( a->fitness < b->fitness )
        return -1;
    else if( a->fitness > b->fitness )
        return 1;
    else
        return 0;
}

individual **copyElite(individual **population, int length, int eliteSize, int popsize){
    individual **elite = (individual **)malloc(sizeof(individual *) * eliteSize);

    qsort(population, popsize, sizeof(individual *), fit_cmp);

    int i;
    for(i=0; i<eliteSize; i++){
        elite[i] = copyIndividual( population[i], length );
    }

    return elite;
}

individual **tournamentSelection( int tsize, individual **population, int popsize ){
    individual **selected = (individual **)malloc(sizeof(individual *)*popsize);
    int i;
    for(i=0; i<popsize;i++){
        selected[i] = tournament(tsize, population, popsize);
    }
    return selected;
}

individual *tournament(int tsize, individual **population, int popsize){
    char used[popsize];
    memset(used, 0, sizeof(char)*popsize);

    individual *individuals[tsize];

    int i;
    for(i=0; i<tsize;i++){
        int idx = rand()%popsize;
        while(used[idx])
            idx = rand()%popsize;
        individuals[i] = population[idx];
    }

    individual *winner = individuals[0];
    for(i=1; i<tsize; i++){
        if(individuals[i]->fitness < winner->fitness){
            winner = individuals[i];
        }
    }

    return winner;
}

void evalPopFitness(individual **population, int popsize){
    int i;
    for(i=0; i<popsize; i++){
        population[i]->fitness = fitness(population[i]);
    }
}

void mutate2(individual *ind, int length){
    int range = rand()%5;
    int s = rand()%(length-range);

    int tmp, i;
    for(i=s; i<s+range-i; i++){
        tmp = ind->perm[i];
        ind->perm[i] = ind->perm[s+range-i];
        ind->perm[s+range-i] = tmp;
    }
}

void mutate1(individual *ind, int length){
    int i1 = rand()%length;
    int i2 = rand()%length;
    while(i2 == i1){
        i2 = rand()%length;
    }
    int tmp = ind->perm[i1];
    ind->perm[i1] = ind->perm[i2];
    ind->perm[i2] = tmp;
}

individual **crossOver(individual *father, individual *mother, int length){
    individual *son1 = (individual *)malloc(sizeof(individual));
    son1->perm = (int *)malloc(length * sizeof(int));

    individual *son2 = (individual *)malloc(sizeof(individual));
    son2->perm = (int *)malloc(length * sizeof(int));

    char used[length];
    memset(used, 0, sizeof(char)*length);

    int fathermap[length],
        mothermap[length];

    int i;
    for(i=0; i<length; i++){
        fathermap[ father->perm[i] ] = i;
        mothermap[ mother->perm[i] ] = i;
    }


    int count = 0;
    int idx, start_idx;
    while( count < length ){
        for(i=0; i<length; i++){
            if(!used[i]){
                idx = i;
                start_idx = idx;
                break;
            }
        }
        while(1){
            idx = mothermap[father->perm[idx]];
            son2->perm[idx] = mother->perm[idx];
            son1->perm[idx] = father->perm[idx];
            count++;
            used[idx] = 1;
            if( idx == start_idx )
                break;
        }
    }

    individual **sons = (individual **)malloc(sizeof(individual *)*2);
    sons[0] = son1;
    sons[1] = son2;
    return sons;
}

void printIndividual(individual *ind, int length){
    printf("%.4lf - ", ind->fitness);

    int i;
    int first = 0;

    for(i=0; i<length; i++){
        if (ind->perm[i] < ind->perm[first])
            first = i;
    }

    for(i=first; i<first + length; i++){
        printf("%d ", ind->perm[i % length]);
    }
    printf("\n");
}

int indexOf(int node, individual *ind, int length){
    int i;
    for(i=0; i<length; i++){
        if(ind->perm[i] == node)
            return i;
    }
}

individual *newIndividual(int length){
    individual *new = (individual *)malloc(sizeof(individual));
    new->perm = (int *)malloc(length * sizeof(int));
    new->fitness = 0;

    int i;
    for(i=0; i<length; i++){
        new->perm[i] = i;
    }

    int tmp, i1, i2;
    for(i=0; i<2 * length; i++){
        i1 = rand()%length;
        i2 = rand()%length;
        tmp = new->perm[i1];
        new->perm[i1] = new->perm[i2];
        new->perm[i2] = tmp;
    }

    return new;
}

individual *copyIndividual(individual *old, int length){
    individual *new = (individual *)malloc(sizeof(individual));
    new->perm = (int *)malloc(length * sizeof(int));
    new->fitness = old->fitness;
    memcpy(new->perm, old->perm, sizeof(int)*length);
    return new;
}

void deleteIndividual(individual *i){
    free(i->perm);
    free(i);
}
