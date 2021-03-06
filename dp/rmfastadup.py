#!/usr/bin/env python
"""Remove duplicates in a fasta file"""
import sys
from dp.associations import GeneAssociations
from dp.ontology import Ontology
from collections import Counter
from dp.utils import parseFasta
seqs = set()
names = set()
fastafile = open(sys.argv[1])

MIN_SEQ_LEN = 32
MAX_SEQ_UNK = 0.1

TAXONS_HOMMO_SAPIENS = {9606}
asoc = GeneAssociations.fromFile(sys.argv[2], taxons = TAXONS_HOMMO_SAPIENS)
ontology = Ontology(sys.argv[3])
ontology.setAssociations(asoc)
asoc.transitiveClosure()
associated = set()
for k,v in asoc.associations.items():
    associated.update({g.upper() for g in v})

ss = dict(parseFasta("data/ss.txt"))
#print(associated)

for l in fastafile:
    name, typ, *_ = l[1:].split(" ")
    name = name.upper()
    seq = next(fastafile)
    sskey = "%s:secstr" % name.replace("_",":")
    if typ != 'mol:protein' \
        or len(seq) < MIN_SEQ_LEN \
        or Counter(seq)['X']/len(seq) > MAX_SEQ_UNK \
        or name not in associated \
        or sskey not in ss \
        or seq in seqs \
        or name in names:
        continue

    sys.stdout.write(l)
    sys.stdout.write(seq)
    seqs.add(seq)
    names.add(name) # This is HACK because blastp is case insensitive

