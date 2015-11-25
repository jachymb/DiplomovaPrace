ASSOCFILE_SER=data/associations/gene_association.goa_pdb.pickle
ASSOCFILE=data/associations/gene_association.goa_pdb
ONTOLOGY=data/obo/go-basic.obo
DATASET=data/dataset.txt
TREELIKER=TreeLiker.jar
BACKGROUND_KNOWLEDGE=data/backround_knowledge.txt
MAX=1298
MIN_ASOC=128
TEMPLATE="res(-a),res(+b),dist(+a,-b,#num),residue(+a, -r1),residue(+b, -r2),charge(+r1, #chrg),charge(+r2, #chrg),polarity(+r1, #pol),polarity(+r2, #pol),secstr(+a, #s),secstr(+b, #s)"
MEMORY=3G
PROCESSES=1

#Blastdbconfiguration
DBNAME=seqres
TASK=blastp #blastp-fast,blastp,blastp-short
VERBOSITY=3
BLAST_THRESHOLD="10e-20"
SAMPLING=64x2 # Sample size x num samples
