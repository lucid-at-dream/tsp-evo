#pragma once

#include "queue.h"

#include <pthread.h>

typedef struct _worker_pool {
    int num_threads;
    void (*do_work)(void *);
    pthread_t *threads;
    queue *work_queue;
    volatile char stop;
    char idle;

    pthread_mutex_t queue_mutex;
    pthread_cond_t  await_idle_cond;
    pthread_cond_t  await_work_cond;
    pthread_cond_t  await_finish_cond;
} worker_pool;

worker_pool *pool_new(int num_threads, void (*)(void *));
void pool_del(worker_pool *pool);

void pool_start(worker_pool *pool);
void pool_stop(worker_pool *pool);

void pool_add_work(worker_pool *pool, void *job);
void pool_await_empty_queue(worker_pool *pool);