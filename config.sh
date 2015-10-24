ASSOCFILE_SER=data/associations/gene_association.goa_pdb.pickle
ASSOCFILE=data/associations/gene_association.goa_pdb
ONTOLOGY=data/obo/goslim_pir.obo
DATASET=data/dataset.txt
TREELIKER=TreeLiker.jar
BACKGROUND_KNOWLEDGE=data/backround_knowledge.txt
MAX_POSITIVE=50
MAX_NEGATIVE=100
MIN_ASOC=100
TEMPLATE="residue(+a, #resType), res(-a)"

# Blast db configuration
DBNAME=seqres
TASK=blastp # blastp-fast, blastp, blastp-short
