#!/bin/sh
# make a copy of original data
source ./config.sh
mydir="data/blast"
mkdir -p $mydir

./rmfastadup.py data/derived/pdb_seqres.txt "$ASSOCFILE" > $mydir/${DBNAME}.fasta

cd $mydir

# Remove old database
rm ${DBNAME}{.phr,.pin,.psq,.psi,.psd,.pog}

# Re-create database
makeblastdb -in ${DBNAME}.fasta -dbtype prot -parse_seqids -out $DBNAME

# Compute blast distances and format output
blastp -db $DBNAME -query ${DBNAME}.fasta -evalue 1.0 -task $TASK -outfmt '10 qseqid sseqid evalue' \
  | awk -F,  '{if(x!=$1) {delete list; x=$1; print x;} if($1 != $2 && !($2 in list)) {list[$2];print $2,$3;} }' \
  | tee blastdist-${TASK}.txt \
  | egrep '^[^ ]*$'
