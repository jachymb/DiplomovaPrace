#!/usr/bin/python -W ignore
import re
import csv
import numpy
import sys
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn import cross_validation
from sklearn.preprocessing import normalize, scale, StandardScaler
from sklearn.metrics import accuracy_score, recall_score, precision_score, confusion_matrix
from pathlib import Path
from collections import Counter

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

__all__ = ['readArff', 'learningTest']

POSTIVE_LABEL = 0
NEGATIVE_LABEL = 1
def readArff(filename):
    """ Parses the ARFF file """
    
    data = []
    labels = []

    def parseLine(line): # csv.reader could not do this.
        isopen = False
        current = ''
        for c in line:
            if c == "'":
                if isopen:
                    yield current
                    current = ''
                isopen = not isopen
            elif isopen:
                current += c

    with filename.open() as f:
        
        line = ''
        while line != '@data':
            line = f.readline().strip()
            if line.startswith("@attribute 'classification'"):
                line = line[line.find('{') + 1:line.find('}')]
                classes = {i:n for n,i in enumerate(parseLine(line))}

        for line in f.read().splitlines():
            record = list(parseLine(line))
            labels.append(classes[record[-1]])
            data.append([int(x) for x in record[:-1]])
    return numpy.array(data, dtype=float), numpy.array(labels), classes

def learningTest(cvdir):
    clasfifiers = (
            ("RBF SVM C=0.5", lambda : SVC(C=0.1)),
            ("RBF SVM C=1", SVC),
            ("RBF SVM C=2", lambda : SVC(C=10)),
            ("RBF SVM C=inf", lambda : SVC(C=numpy.inf)),
            ("Linear SVM C=1", lambda: SVC(kernel='linear')),
            ("Quadratic SVM C=1", lambda: SVC(kernel='poly', degree=2)),
            ("Cubic SVM C=1", lambda: SVC(kernel='poly', degree=3)),
            ("Random Forest", RandomForestClassifier),
            )

    scaler = StandardScaler()
    results = []
    with (cvdir / 'scores.txt').open('w') as output:
        print('Cross Validation results for file term %s:' % cvdir.name, file=output)
        print('Confusion Matrix: [[True positive, False postive], [False negative, True negative]]', file=output)
        counts = Counter(
                   (re.search(r'^"(.*?)"', l).groups()[0]
                    for l in (cvdir / 'dataset.txt').open()
                    if l.strip() != ''))
        print('Dataset entry counts: ' +
                ", ".join(
                    ("%s = %d" % item
                     for item
                     in sorted(counts.items()))),
            file = output)

        for name, Clf in clasfifiers:
            scores = []
            classifiers = []
            print('Testing the classifier %s:' % name, file=output)

            for train, test in zip(sorted(cvdir.glob('train_*.arff')), sorted(cvdir.glob('test_*.arff'))):
                clf = Clf()
                i = int(train.name[1+train.name.find('_'):train.name.find('.')])

                X_train, y_train, _ = readArff(train)
                X_test , y_test,  _ = readArff(test)

                X_train = scaler.fit_transform(X_train)
                X_test = scaler.transform(X_test, copy=True)
               
                pos = (y_train == POSTIVE_LABEL)
                neg = (y_train == NEGATIVE_LABEL)
                posWeight = numpy.sum(neg) / len(y_train)
                negWeight = numpy.sum(pos) / len(y_train)
                sample_weight = posWeight*pos + negWeight*neg

                clf.fit(X_train, y_train, sample_weight = sample_weight)
                y_pred = clf.predict(X_test)
                score = accuracy_score(y_test, y_pred)
                #prec = precision_score(y_test, y_pred, average='weighted')
                #reca = recall_score(y_test, y_pred, average='weighted')
                conf = confusion_matrix(y_test, y_pred)
                print("Fold %d score: %.2f, confusion matrix: %s" % (i, score*100.0, conf.tolist()), file=output)
                scores.append(score)
                clf.conf = conf
                clf.cvindex = i
                clf.name = name
                classifiers.append(clf)
                
            scores = numpy.array(scores)
            print("Average: %.2f%% (+/- %.2f%%)" % (scores.mean()*100.0, scores.std() * 2 * 100.0), file=output)
            print("Best: %.2f%%\n" % (scores.max()*100.0), file=output, flush=True)
            bestPerforming = numpy.argmax(scores)
            bestClassifier = classifiers[bestPerforming]
            bestClassifier.cvindex = bestPerforming

            with (cvdir / (name.replace(' ','_')+'.pickle')).open('wb') as bestSer:
                pickle.dump(bestClassifier, bestSer)

            #results.append(bestClassifier)
            yield bestClassifier

if __name__ == "__main__":
    cvdir = Path(sys.argv[1])
    for best in learningTest(cvdir):
        print(best.name, "done.")
