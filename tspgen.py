#!/bin/python

from random import random, randint

from sys import argv

N = randint(5,20)
N = 1000
N = int(argv[1])

print(N)
for i in range(N):
	for j in range(N):
		if( i == j ):
			print(i,j,0)
		else:
			print(i,j,random())
