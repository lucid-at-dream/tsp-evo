#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#define MAXN 100

//problem variables
int N;
double costs[MAXN][MAXN];
double costs_v[MAXN * MAXN];

//best solution storage
double best = 1000;
int bestperm[MAXN];

//auxiliar variables
int perm[MAXN];
char visited[MAXN];

//bound related variables
double mincost[MAXN+1];

int cmp (const void *ls, const void* rs);
inline double lowerbound( double cost, int idx );
void solve();
void tsp(double cost, int idx);

int main(){
    int i,j;
    double cost;
    int count = 0;

    //read input
    scanf("%d", &N);

    double mins[MAXN];
    memset(mins, 100, N * sizeof(double));
    while( scanf("%d %d %lf", &i, &j, &cost ) != EOF ) {
        costs[i][j] = cost;
        if(cost < mins[i] && cost != 0){
            mins[i] = cost;
        }
    }

    //initialize lowerbounds
    qsort(mins, N, sizeof(double), cmp);

    mincost[0] = 0;
    for( i=1; i<N+1; i++ ){
        mincost[i] = mincost[i-1] + mins[i-1];
    }

    //initialize auxiliary vars
    memset(visited, 0, MAXN * sizeof(char));

    solve();

    //gen output
    char separator = ' ';
    for( i=0; i<N; i++ ){
        printf("%d%c", bestperm[i], i != N-1 ? separator : '\n');
    }
    printf("%lf\n", best);

    return 0;
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

void solve(){
    perm[0] = 0;
    visited[0] = 1;

    int i;
    for(i=1; i<N; i++) {
        perm[1] = i;
        visited[i] = 1;
        tsp(costs[0][i], 2);
        visited[i] = 0;
    }
    visited[0] = 0;
}

void tsp(double cost, int idx){

    if( cost + mincost[N-idx-1] >= best )
        return;

    int i;
    if( idx == N ){
        cost += costs[perm[idx-1]][0];
        if(cost < best){
            best = cost;
            memcpy(bestperm, perm, N*sizeof(int));
        }
        return;
    }

    for(i=1; i<N; i++){
        if( !visited[i] ){
            perm[idx] = i;
            visited[i] = 1;
            tsp( cost + costs[perm[idx-1]][i], idx+1 );
            visited[i] = 0;
        }
    }
}



