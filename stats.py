#!/usr/bin/python

from numpy import mean
from os import listdir

ninputs = 3

def get_ind(rundir):
    ind = {}

    with open("runs/%s/config" % rundir, "r") as f:
        for i in [j.strip() for j in f.read().split("\n") if j != "" and j[0] != "#"]:
            ind[i.split(" ")[0]] = float(i.split(" ")[-1])

    scores = []
    for i in range(ninputs):
        with open("runs/%s/out_%d" % (rundir, i), 'r') as f:
            score = [float(i) for i in f.read().split("\n") if i != ""]
            scores += [mean(score)]
    ind["#fitness"] = mean(scores)

    return ind

inds = []
for i in listdir('runs'):
    try:
        inds += [get_ind(i)]
    except:
        pass

from sys import argv
field = argv[1]

data = sorted([(x[field], x["#fitness"]) for x in inds])

x = [i[0] for i in data]
y = [i[1] for i in data]

from matplotlib import pyplot as p
p.title(field)
p.plot(x, y, '.')
p.show()

