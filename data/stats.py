from os import listdir
from numpy import std, mean
from sys import argv

d = argv[1]

best_files = [i for i in listdir(d) if i.endswith('.best')]

for fp in best_files:
  with open(d + "/" + fp, 'r') as f:
    data = [float(i) for i in f.read().split("\n") if i != ""]

  print("%s - %.4f (+%.4f)" % (fp.split(".")[0], mean(data), std(data)))

