#!/usr/bin/env python -W ignore
import re
import csv
import numpy
import sys
import pickle
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn import cross_validation
from sklearn.preprocessing import normalize, scale, StandardScaler
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score, recall_score, precision_score, confusion_matrix, precision_recall_curve
from pathlib import Path
from collections import Counter

from sklearn.decomposition import PCA
from sklearn.lda import LDA

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
            #TOOD, use third dict for params
            ("SGD", lambda: SGDClassifier(n_iter=100,alpha=0.01)),
            ("RBF SVM C=0.5", lambda : SVC(C=0.1)),
            ("RBF SVM C=1", lambda: SVC(shrinking=False, tol=1e-5)),
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
                pos = (y_test == POSTIVE_LABEL)
                neg = (y_test == NEGATIVE_LABEL)
                score = accuracy_score(y_test[pos], y_pred[pos])*posWeight + accuracy_score(y_test[neg], y_pred[neg])*negWeight
                #prec = precision_score(y_test, y_pred, average='weighted')
                #reca = recall_score(y_test, y_pred, average='weighted')
                conf = confusion_matrix(y_test, y_pred)
                print("Fold %d score: %.2f, confusion matrix: %s" % (i, score*100.0, conf.tolist()), file=output)
                scores.append(score)
                clf.conf = conf
                clf.cvindex = i
                clf.name = name
                classifiers.append(clf)

                #y_score = clf.decision_function(X_test)
                #precision, recall, _ = precision_recall_curve(y_test, y_score)
                #Plot Precision-Recall curve
                #plt.clf()
                #plt.plot(recall, precision, label='Precision-Recall curve')
                #plt.plot([1,0], [0,1], label='id')
                #plt.xlabel('Recall')
                #plt.ylabel('Precision')
                #plt.ylim([0.0, 1.05])
                #plt.xlim([0.0, 1.0])
                #plt.title('Precision-Recall '+name)
                #plt.legend(loc="lower left")
                #plt.show()
                #break

                #import matplotlib.pyplot as plt
                #target_names = ["pos","neg"]
                #pca = PCA(n_components=2)
                #X_r = pca.fit(X_train).transform(X_train)
                #print(X_r.shape)
                #plt.figure()
                #for c, i, target_name in zip("rg", [0, 1], target_names):
                #    plt.scatter(X_r[y_train == i, 0], X_r[y_train == i, 1], c=c, label=target_name)
                #plt.legend()
                #plt.title('PCA')

                #lda = LDA(n_components=10)
                #X_r2 = lda.fit(X_train, y_train).transform(X_test)
                #plt.figure()
                #print(X_r2.shape)
                #for c, i, target_name in zip("rg", [0, 1], target_names):
                #    plt.scatter(X_r2[y_test == i, 0], X_r2[y_test == i, 0], c=c, label=target_name)
                #plt.legend()
                #plt.title('LDA')

                #plt.show()
                
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
