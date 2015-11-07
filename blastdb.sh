#!/bin/sh
source ./config.sh
export PYTHONPATH="$PWD"

# make a copy of original data and filter the data
mydir="data/blast"
mkdir -p $mydir
fasta=${DBNAME}.fasta
dp/rmfastadup.py data/derived/pdb_seqres.txt "$ASSOCFILE" "$ONTOLOGY"> "$mydir/$fasta"

# Remove old database
cd "$mydir"
rm ${DBNAME}{.phr,.pin,.psq,.psi,.psd,.pog}

# Re-create database
makeblastdb -in "$fasta" -dbtype prot -parse_seqids -out "$DBNAME"

# Compute blast distances and format output
dists=blastdist-${TASK}.txt
blastp -db "$DBNAME" -query "$fasta" -evalue 1.0 -task $TASK -outfmt '10 qseqid sseqid evalue' \
  | awk -F,  '{if(x!=$1) {delete list; x=$1; print x;} if($1 != $2 && !($2 in list)) {list[$2];print $2,$3;} }' \
  | tee "$dists" \
  | egrep '^[^ ]*$'

# Generate new independent datasets
cd -
dp/maximumIndependentSet.py "$mydir/$fasta" "$mydir/$dists" "$BLAST_THRESHOLD" "$DATASET"
