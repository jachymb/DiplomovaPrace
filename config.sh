ASSOCFILE_SER=data/associations/gene_association.goa_pdb.pickle
ASSOCFILE=data/associations/gene_association.goa_pdb
ONTOLOGY=data/obo/goslim_pir.obo
DATASET=data/dataset.txt
TREELIKER=TreeLiker.jar
BACKGROUND_KNOWLEDGE=data/backround_knowledge.txt
MAX_POSITIVE=500
MAX_NEGATIVE=1000
MIN_ASOC=100
#TEMPLATE="residue(+a,#resType),res(-a)"
TEMPLATE="residue(+a,-r1),residue(+a,#resType),res(-a),res(-b),residue(+b,#resType),residue(+b,-r2),dist(-a,!resType,-b,!resType,#num),next(-a,-b),positive(+r1),negative(+r1),neutral(+r1),nonpolar(+r1),polar(+r1),positive(+r2),negative(+r2),neutral(+r2),nonpolar(+r2),polar(+r2)"

#Blastdbconfiguration
DBNAME=seqres
TASK=blastp#blastp-fast,blastp,blastp-short
