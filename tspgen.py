#!/bin/python

from random import random, randint

def generate(N):
    print(N)
    for i in range(N):
        for j in range(N):
            if (i == j):
                print(i,j,0)
            else:
                print(i,j,random())

if __name__ == '__main__':
    from sys import argv
    prob_size = int(argv[1])
    generate(prob_size)