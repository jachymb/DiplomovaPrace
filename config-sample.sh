ASSOCFILE_SER=data/associations/gene_association.goa_pdb.pickle
ASSOCFILE=data/associations/gene_association.goa_pdb
#ONTOLOGY=data/obo/goslim_pir.obo
ONTOLOGY=data/obo/go-basic.obo
DATASET=data/dataset.txt
TREELIKER=TreeLiker.jar
BACKGROUND_KNOWLEDGE=data/backround_knowledge.txt
MAX=512
MIN_ASOC=128
#MAX=2048
#MIN_ASOC=1024
#TEMPLATE="residue(+a,#resType),res(-a)"
#TEMPLATE="residue(+a,-r1),residue(+a,#resType),res(-a),res(-b),residue(+b,#resType),residue(+b,-r2),dist(-a,!resType,-b,!resType,#num),next(-a,-b),positive(+r1),negative(+r1),neutral(+r1),nonpolar(+r1),polar(+r1),positive(+r2),negative(+r2),neutral(+r2),nonpolar(+r2),polar(+r2)"
TEMPLATE="res(-a),res(+b),dist(+a,-b,#num),residue(+a, -r1),residue(+b, -r2),charge(+r1, #chrg),charge(+r2, #chrg),polarity(+r1, #pol),polarity(+r2, #pol),secstr(+a, #s),secstr(+b, #s)"
#MEMORY=3584M
MEMORY=3G
PROCESSES=1
#Blastdbconfiguration
DBNAME=seqres
TASK=blastp #blastp-fast,blastp,blastp-short
VERBOSITY=3
BLAST_THRESHOLD="10e-20"
#RESERVED=0.1
