#!/bin/sh
#vmtouch -ld `awk '{print "GENES/"$1".pickle"}' < $DATASET`

export PYTHONWARNINGS=ignore

source ./config.sh

#taskset -c 1 ./main.py \
./main.py \
  --max-positive "$MAX_POSITIVE" \
  --max-negative "$MAX_NEGATIVE" \
  --min-associations "$MIN_ASOC" \
  --dataset "$DATASET" \
  --background-knowledge "$BACKGROUND_KNOWLEDGE" \
  --template "$TEMPLATE" \
  --treeliker "$TREELIKER" \
  --memory "$MEMORY" \
  --processes "$PROCESSES" \
  --verbosity "$VERBOSITY" \
  --recalculate-distances \
  --reserve "$RESERVED" \
  "$ONTOLOGY" \
  "$ASSOCFILE_SER"
