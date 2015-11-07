#!/usr/bin/env python
from dp.gene import *
from dp.utils import parallel_map_dill
Gene.recalculateDists = True
f = GeneFactory(deserialize=True)
with open('data/dataset.txt','r') as genes:
	parallel_map_dill(20, f.getGene, (g.strip() for g in genes))
	#for gene in genes:
	#	f.getGene(gene.strip())
