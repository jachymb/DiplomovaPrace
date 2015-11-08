#!/bin/sh
#vmtouch -ld `awk '{print "GENES/"$1".pickle"}' < $DATASET`

export PYTHONWARNINGS=ignore

source ./config.sh

#taskset -c 1 ./main.py \
./main.py \
  --max "$MAX" \
  --min-associations "$MIN_ASOC" \
  --dataset "$DATASET" \
  --background-knowledge "$BACKGROUND_KNOWLEDGE" \
  --template "$TEMPLATE" \
  --treeliker "$TREELIKER" \
  --memory "$MEMORY" \
  --processes "$PROCESSES" \
  --verbosity "$VERBOSITY" \
  "$ONTOLOGY" \
  "$ASSOCFILE_SER"
