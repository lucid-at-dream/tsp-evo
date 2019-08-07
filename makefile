
CCFLAGS=-Wall -lm -lpthread -O2 -mtune=native

all: conf.c conf.h queue.c queue.h worker_pool.c worker_pool.h tspevo.c
	gcc -c queue.c ${CCFLAGS}
	gcc -c worker_pool.c ${CCFLAGS}
	gcc -c conf.c ${CCFLAGS}
	gcc *.o tspevo.c -o t ${CCFLAGS}

clean:
	rm -rf *.o *.ghc cachegring.* t
