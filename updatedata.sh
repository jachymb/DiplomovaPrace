#!/bin/sh
if [ ! -f config.sh ]; then
  cp config-sample.sh config.sh
fi
source ./config.sh
                                                    # port PDBe server is using

mkdir -p results
mkdir -p data/associations data/derived_data data/obo data/serialized_genes data/xml 

cd data/associations
if [ ! -f gene_association.goa_pdb ]; then
  wget -N http://geneontology.org/gene-associations/gene_association.goa_pdb.gz
  yes | gunzip -f gene_association.goa_pdb.gz
fi
#cp gene_association.goa_pdb gene_association.goa_pdb_reserved

cd ../..

./main.py \
  --dump-associations "$ASSOCFILE_SER" \
  --dataset "$DATASET" \
  "$ONTOLOGY" \
  "$ASSOCFILE"

cd data

RSYNC=/usr/bin/rsync                             # location of local rsync
SERVER=rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated     # PDBe server name
PORT=873  
${RSYNC} -rlpt -v -z --delete --port=$PORT ${SERVER}/derived_data/ derived 2>/dev/null
${RSYNC} -rlpt -v -z --delete --port=$PORT ${SERVER}/data/structures/divided/XML/ xml 2>/dev/null

cd obo
wget -N http://www.geneontology.org/ontology/subsets/goslim_{generic,plant,candida,pir,pombe,yeast,aspergillus,metagenomics,virus,chembl}.obo

wget -N http://purl.obolibrary.org/obo/go/go-basic.obo
wget -N http://purl.obolibrary.org/obo/go.obo


