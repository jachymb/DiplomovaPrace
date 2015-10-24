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
from sklearn.metrics import accuracy_score, recall_score, precision_score
from pathlib import Path
from collections import Counter

FOLDS = 10
CVTESTS = 10
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
    clasfifiers = (("SVM", SVC), ("Random Forest", RandomForestClassifier))

    scaler = StandardScaler()
    results = []
    with (cvdir / 'scores.txt').open('w') as output:
        print('Cross Validation results for file term %s:' % cvdir.name, file=output)
        print('Dataset entry counts: ' +
                ", ".join(
                    ("%s - %d" % item
                        for item
                        in sorted(
                            Counter(
                                (re.search(r'^"(.*?)"', l).groups()[0]
                                 for l in (cvdir / 'dataset.txt').open()
                                 if l.strip() != '')
                            ).items()))),
            file = output)

        for name, Clf in clasfifiers:
            scores = []
            classifiers = []
            print('Testing the classifier %s:' % name, file=output)

            for train, test in zip(sorted(cvdir.glob('train_*.arff')), sorted(cvdir.glob('test_*.arff'))):
                clf = Clf()
                i = train.name[1+train.name.find('_'):train.name.find('.')]

                X_train, y_train, _ = readArff(train)
                X_test , y_test,  _ = readArff(test)

                X_train = scaler.fit_transform(X_train)
                X_test = scaler.transform(X_test, copy=True)

                clf.fit(X_train, y_train)
                y_pred = clf.predict(X_test)
                score = accuracy_score(y_test, y_pred)
                prec = precision_score(y_test, y_pred, average='weighted')
                reca = recall_score(y_test, y_pred, average='weighted')
                print("Fold %s score: %.2f, weighted precision: %.2f, weighted recall: %.2f" % (i, score*100.0, prec * 100.0, reca % 100.0), file=output)
                scores.append(score)
                classifiers.append(clf)
                
            scores = numpy.array(scores)
            print("Average: %.2f%% (+/- %.2f%%)" % (scores.mean()*100.0, scores.std() * 2 * 100.0), file=output)
            print("Best: %.2f%%\n" % (scores.max()*100.0), file=output)
            bestPerforming = numpy.argmax(scores)
            bestClassifier = classifiers[bestPerforming]
            bestClassifier.cvindex = bestPerforming

            with (cvdir / (name.replace(' ','_')+'.pickle')).open('wb') as bestSer:
                pickle.dump(bestClassifier, bestSer)

            #results.append(bestClassifier)
            yield bestClassifier
    #return results


if __name__ == "__main__":
    cvdir = Path(sys.argv[1])
    learningTest(cvdir)

#Counter(map(lambda l: re.search('^".*?"', l).group(), filter(lambda l: l.strip() != '', open('dataset.txt'))))
