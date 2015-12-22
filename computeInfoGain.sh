#!/bin/sh

RESULTS=results-2k
for term in `ls $RESULTS|sort`; do
  gainall=$RESULTS/$term/gainall.txt
  if [ -f $gainall ]; then rm $gainall; fi
  for fold in `seq 0 7`; do
    outf=$RESULTS/$term/$fold/gainratio.txt
    if [ ! -f $outf ]; then
      java -Xmx2G -cp /usr/share/java/weka/weka.jar weka.attributeSelection.GainRatioAttributeEval \
            -M \
            -s "weka.attributeSelection.Ranker -N -1" \
            -i $RESULTS/$term/$fold/train.arff > $outf
    fi
    echo Done fold $fold on term $term
    head -22 $outf | tail -n+13 | sed 's/.* [0-9][0-9]* //' >> $gainall
  done
  sort $gainall | uniq -c | sort -rn > tex/selected-$term.txt
done

