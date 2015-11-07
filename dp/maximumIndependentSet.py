#!/usr/bin/env python
from itertools import groupby, islice 
from copy import deepcopy
import random
import pickle
import time
import sys
import os

#ALL = 'data/blast/seqres.fasta'
#DISTS  = 'data/blast/blastdist-blastp.txt'
#THRESHOLD = 10e-20

class Grouper:
    def __init__(self, n=2):
        self.n = n
        self.c = 0
        self.g = True
    def __call__(self, a):
        self.c %= self.n
        if self.c == 0: self.g = not self.g
        self.c += 1
        return self.g

class Graph:
    def __init__(self, threshold, dfname, distsFile, all, graph = None, best = set()):
        self.graph = graph
        self.best = best
        if self.graph is None:
            self.graph = {}
            print("reading file "+all)
            with open(all) as seqres:
                for _,(desc, seq) in groupby(seqres, key=Grouper()):
                    name = desc.split(" ")[0][1:]
                    self.graph[name] = set()

            print("reading file "+distsFile)
            with open(distsFile) as dists:
                key = None
                for line in dists:
                    line = line.strip().split(' ')
                    if len(line) == 1:
                        key, = line
                    else:
                        dist = float(line[1])
                        if dist < threshold and key in self.graph:
                            self.graph[key].add(line[0])

            self.checkEdges()

        self.dfname = dfname
        print("graph initialised")

    def checkEdges(self):
        print("Checking edges")
        # Check if edges are both-way
        for v, ne in self.graph.items():
            for n in ne:
                if n in self.graph and v not in self.graph[n]:
                    #print("Adding edge %s -> %s!" % (n,v))
                    self.graph[n].add(v)

    @classmethod
    def load(cls, fname='graph.pickle'):
        print("loading graph from "+fname)
        with open(fname, 'rb') as graphF:
            data = pickle.load(graphF)
            return data
            #return cls(*data)

    def dump(self, fname='data/graph.pickle'):
        print("dumping graph to "+fname)
        with open(fname, 'wb') as graphFile:
            data = (self.graph, self.best)
            pickle.dump(graph, graphFile)
        with open(self.dfname, 'w') as dataset:
            dataset.write("\n".join(sorted(self.best)))

    def singleMaximalIndependentSet(self):
        #print("running singleMaximalIndependentSet()")
        graph = deepcopy(self.graph)
        independent = set()
        while graph:
            node = random.choice(tuple(graph.keys()))
            independent.add(node)
            for neighbor in graph.pop(node):
                try:
                    del graph[neighbor]
                except KeyError:
                    pass


        if len(independent) > len(self.best):
            self.best = independent
            print("Found better independent set of size %d!" % len(self.best))
            self.dump()

        return independent

    def randomMaximalIndependentSets(self, iters):
        bestL = len(self.best)
        for x in range(iters):
            try:
                self.singleMaximalIndependentSet()
            except KeyboardInterrupt:
                break

    def bestRatio(self):
        return len(self.best), len(self.graph)


if __name__ == "__main__":
    if len(sys.argv) not in (5,6):
        sys.stderr.write("Usage: %s FASTAFILE BLAST-DISTANCES THRESHOLD OUTPUT [ITERS]\n")
        sys.exit(1)
    if len(sys.argv) == 5: sys.argv.append(20)
    ALL, DISTS, THRESHOLD, DATASET, ITERS = sys.argv[1:] 

    graph = Graph(float(THRESHOLD), DATASET, DISTS, ALL)
    #graph = Graph.load()
    graph.randomMaximalIndependentSets(ITERS)

