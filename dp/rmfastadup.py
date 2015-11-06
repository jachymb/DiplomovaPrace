#!/usr/bin/python
"""Remove duplicates in a fasta file"""
import sys
from dp.associations import GeneAssociations
from dp.ontology import Ontology
from collections import Counter
seqs = set()
names = set()
fastafile = open(sys.argv[1])

MIN_SEQ_LEN = 32
MAX_SEQ_UNK = 0.1

asoc = GeneAssociations.fromFile(sys.argv[2])
ontology = Ontology(sys.argv[3])
ontology.setAssociations(asoc)
asoc.transitiveClosure()
associated = set()
for k,v in asoc.associations.items():
    associated.update({g.upper() for g in v})

#print(associated)

for l in fastafile:
    name, typ, *_ = l[1:].split(" ")
    name = name.upper()
    seq = next(fastafile)
    if typ != 'mol:protein' \
        or len(seq) < MIN_SEQ_LEN \
        or Counter(seq)['X']/len(seq) > MAX_SEQ_UNK \
        or name not in associated:
        continue
    if seq not in seqs and name not in names:
        sys.stdout.write(l)
        sys.stdout.write(seq)
        seqs.add(seq)
        names.add(name) # This is HACK because blastp is case insensitive

