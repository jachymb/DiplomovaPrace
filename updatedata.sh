#!/bin/sh
source ./config.sh

RSYNC=/usr/bin/rsync                             # location of local rsync
SERVER=rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated     # PDBe server name
PORT=873                                                      # port PDBe server is using

mkdir -p results
mkdir -p data/associations data/derived_data data/obo data/serialized_genes data/xml 
cd data

${RSYNC} -rlpt -v -z --delete --port=$PORT ${SERVER}/derived_data/ derived 2>/dev/null
${RSYNC} -rlpt -v -z --delete --port=$PORT ${SERVER}/data/structures/divided/XML/ xml 2>/dev/null

cd associations
wget -N http://geneontology.org/gene-associations/gene_association.goa_pdb.gz
yes | gunzip gene_association.goa_pdb.gz

cd ../obo
wget -N http://www.geneontology.org/ontology/subsets/goslim_{generic,plant,candida,pir,pombe,yeast,aspergillus,metagenomics,virus,chembl}.obo

cd ../.. 


./main.py \
  --dump-associations "$ASSOCFILE_SER" \
  --dataset "$DATASET" \
  "$ONTOLOGY" \
  "$ASSOCFILE"

