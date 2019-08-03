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

#define merda printf("merda\n");
#define cona(x) printf("cona--%d\n", x);

#define MAXN 1001
#define VERBOSE

//problem variables
int N;
double costs[MAXN][MAXN];

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

void magic(
    int ngens,
    double mutation_rate, double crossover_rate, double migration_rate, double migrationProb,
    double replacementProb, double replacementInc, double replacementProbMax,
    int popsize, int tournamentSize, int eliteSize,
    int npops,
    int maxgrad0count
);

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
    double mutation_rate = 0.6;
    double crossover_rate = 0.3;
    double migration_rate = 0.2;
    double migrationProb = 0.1;
    double replacementProb = 0;
    double replacementInc = 0.0002;
    double replacementProbMax = 0.15;
    int popsize = 250;
    int tournamentSize = 5;
    int eliteSize = 2;
    int npops = 50;
    int maxgrad0count = 200;

    magic(
        ngens,
        mutation_rate, crossover_rate,
        migration_rate, migrationProb,
        replacementProb, replacementInc, replacementProbMax, 
        popsize, tournamentSize, eliteSize, 
        npops,
        maxgrad0count
    );

    return 0;
}

void magic(
    int ngens,
    double mutation_rate,
    double crossover_rate,
    double migration_rate,
    double migrationProb,
    double replacementProb,
    double replacementInc,
    double replacementProbMax,
    int popsize,
    int tournamentSize,
    int eliteSize,
    int npops,
    int maxgrad0count
) {

    int i,j;
    double cost;

    srand(clock());

#ifdef VERBOSE
    printf("starting the algorithm\n");
#endif

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
#ifdef VERBOSE
            printf("[%d] best of best: %.4lf (-%.4lf)\n", best->fitness, prev_best_fitness);
#endif
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
        if( (rand() % 10000)/10000.0 <= migration_rate ){
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

#ifdef VERBOSE
        printf("[%d] %lf seconds\n", gen, (double)((gen_finish - gen_start)/(double)CLOCKS_PER_SEC) );
#endif
    }

    //fetch the best individuals
    for(p=0; p<npops; p++){
        qsort(populations[p], popsize, sizeof(individual *), fit_cmp);
    }

    finish = clock();

#ifdef VERBOSE
    printf("the algorithm ran in %lf seconds\n", (double)((finish - start)/(double)CLOCKS_PER_SEC) );
#endif

    int pos[npops];
    memset(pos, 0, sizeof(int)*npops);

    int minpop = 0;
    individual *min;

#ifdef VERBOSE
    for(i=0; i<npops; i++){

        for(j=0; pos[j] < 0; j++);

        min = populations[0][pos[j]];
        for(p=0; p<npops; p++){
            if (pos[p] < 0)
                continue;
            if(populations[p][pos[p]]->fitness < min->fitness){
                min = populations[p][pos[p]];
                minpop = p;
            }
        }
        pos[minpop] = -1;
        printIndividual(min, N);
    }

    printf("best of best: ");
    printIndividual(best, N);
#endif
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

    qsort(population, popsize, sizeof(individual *), fit_cmp);

    return population;
}

