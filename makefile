
CFLAGS=-Wall -lm -lpthread
COPT=-O2 -mtune=native
#CDEB=-Og -g
CCFLAGS=${CFLAGS} ${COPT} ${CDEB}


all: conf.c conf.h queue.c queue.h worker_pool.c worker_pool.h tspevo.c
	gcc -c queue.c ${CCFLAGS}
	gcc -c worker_pool.c ${CCFLAGS}
	gcc -c conf.c ${CCFLAGS}
	gcc -c tsp.c ${CCFLAGS}
	gcc -c tspevo.c ${CCFLAGS}
	gcc *.o main.c -o t ${CCFLAGS}

clean:
	rm -rf *.o *.ghc cachegring.* t
