#!/usr/bin/env python
from dp.gene import *
from dp.utils import parallel_map_dill
import os
import random
Gene.recalculateDists = True
f = GeneFactory(deserialize=True)
os.system("taskset -p 0xff %d" % os.getpid())

i = 1
def regen(genes):
    genes = list(genes)
    l = len(genes)
    random.shuffle(genes)
    def r(g):
        global i
        f.getGene(g)
        print("%d/%d" % (i,l))
        i += 1
    parallel_map_dill(2, r, (g.strip() for g in genes))

with open('data/dataset.txt','r') as genes:
    regen(genes)

regen((x.split(".")[0] for x in os.listdir("data/serialized_genes")))
