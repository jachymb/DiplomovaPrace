#!/usr/bin/python
""" Remove proteins without secondary structure from dataset """

from dp.utils import parseFasta

dataset = "data/dataset.txt"
g = dict(parseFasta("data/ss.txt"))
with open(dataset, "r") as f:
    data = f.read().splitlines()


with open(dataset, "w") as f:
    r = []
    for g in data:
        name, strand = g.split("_")
        name = name.upper()
        if ("%s:%s:secstr" % (name, strand)) in g:
            r.append(g)
    f.write("\n".join(r))

