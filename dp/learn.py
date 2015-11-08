#!/usr/bin/env python
import re
import csv
import numpy
import sys
import pickle
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn import cross_validation
from sklearn.preprocessing import normalize, scale, StandardScaler
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score, recall_score, precision_score, confusion_matrix, precision_recall_curve, roc_curve, auc, average_precision_score
from scipy import interp
from pathlib import Path
from collections import Counter
from dp.utils import NUM_FOLDS, in_directory

from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA

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
            ("AdaBoost-DecisionTree", AdaBoostClassifier),
            ("5-NN", lambda: KNeighborsClassifier(p=1, algorithm='kd_tree')),
            ("Random Forest", RandomForestClassifier),
            ("SGD", lambda: SGDClassifier(n_iter=100,alpha=0.01,loss="modified_huber")),
            ("RBF SVM C=1", lambda: SVC(shrinking=False, tol=1e-5,probability=True)),
            ("RBF SVM C=0.5", lambda : SVC(C=0.1,probability=True)),
            ("RBF SVM C=2", lambda : SVC(C=10,probability=True)),
            ("RBF SVM C=inf", lambda : SVC(C=numpy.inf,probability=True)),
            ("Linear SVM C=1", lambda: SVC(kernel='linear',probability=True)),
            ("Quadratic SVM C=1", lambda: SVC(kernel='poly', degree=2,probability=True)),
            ("Cubic SVM C=1", lambda: SVC(kernel='poly', degree=3,probability=True)),
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
            data = []
            print('Testing the classifier %s:' % name, file=output)
            clfdir = cvdir / name.replace(' ','_')
    
            for i in range(NUM_FOLDS): # NUM_FOLDS
                foldDir = cvdir / str(i)
                train = foldDir / 'train.arff'
                test = foldDir / 'test.arff'

                clf = Clf()

                X_train, y_train, _ = readArff(train)
                X_test , y_test,  _ = readArff(test)

                X_train = scaler.fit_transform(X_train)
                X_test = scaler.transform(X_test, copy=True)
               
                pos = (y_train == POSTIVE_LABEL)
                neg = (y_train == NEGATIVE_LABEL)
                posWeight = numpy.sum(neg) / len(y_train)
                negWeight = numpy.sum(pos) / len(y_train)
                sample_weight = posWeight*pos + negWeight*neg

                try:
                    clf.fit(X_train, y_train, sample_weight = sample_weight)
                except TypeError:
                    clf.fit(X_train, y_train)
                y_pred = clf.predict(X_test)
                pos = (y_test == POSTIVE_LABEL)
                neg = (y_test == NEGATIVE_LABEL)
                score = accuracy_score(y_test[pos], y_pred[pos])*posWeight + accuracy_score(y_test[neg], y_pred[neg])*negWeight
                #prec = precision_score(y_test, y_pred, average='weighted')
                #reca = recall_score(y_test, y_pred, average='weighted')
                conf = confusion_matrix(y_test, y_pred)
                print("Fold %d score: %.2f, confusion matrix: %s" % (i, score*100.0, conf.tolist()), file=output)
                clf.conf = conf
                clf.cvindex = i
                clf.name = name
                classifiers.append(clf)
                data.append((X_test, y_test))

            with in_directory(clfdir):
                legendprop = {'size': 10}
                #Plot Precision-Recall curve
                plt.clf()
                #plt.plot([0, 0], [1, 1], '--', color=(0.6, 0.6, 0.6))
                y_tests = []
                y_scores = []
                for i, (clf, (X_test, y_test)) in enumerate(zip(classifiers, data)):
                    try:
                        y_score = clf.decision_function(X_test)
                    except AttributeError:
                        y_score = clf.predict_proba(X_test)[:, 0]
                    precision, recall, _ = precision_recall_curve(y_test, y_score)
                    y_tests.extend(y_test)
                    y_scores.extend(y_score)

                    area = average_precision_score(y_test, y_score)
                    clf.pr_auc = area
                    plt.plot(recall, precision, label='Fold %d, AUC = %0.2f' % (i, area), lw=1)

                precision, recall, _ = precision_recall_curve(y_tests, y_scores)
                area = average_precision_score(y_tests, y_scores)

                plt.plot(recall, precision, 'k--', label='Mean, AUC = %0.2f' % (area), lw=2)
                plt.xlabel('Recall')
                plt.ylabel('Precision')
                plt.xlim([-0.05, 1.05])
                plt.ylim([-0.05, 1.05])
                plt.title('Precision-Recall: '+name)
                plt.legend(loc="lower center", prop=legendprop)
                plt.savefig('precision-recall.png')

                #Plot ROC curve
                plt.clf()
                plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6))
                mean_tpr = 0.0
                mean_fpr = numpy.linspace(0, 1, 100)
                all_tpr = []
                for i, (clf, (X_test, y_test)) in enumerate(zip(classifiers, data)):
                    probabs = clf.predict_proba(X_test)
                    fpr, tpr, _ = roc_curve(y_test, probabs[:, 1])
                    mean_tpr += interp(mean_fpr, fpr, tpr)
                    mean_tpr[0] = 0.0
                    roc_auc = auc(fpr, tpr)
                    clf.roc_auc = roc_auc
                    scores.append(roc_auc)
                    plt.plot(fpr, tpr, lw=1, label='Fold %d, AUC = %0.2f' % (i, roc_auc))


                mean_tpr /= NUM_FOLDS
                mean_tpr[-1] = 1.0
                mean_auc = auc(mean_fpr, mean_tpr)
                plt.plot(mean_fpr, mean_tpr, 'k--',
                        label='Mean, AUC = %0.2f' % mean_auc, lw=2)

                plt.xlim([-0.05, 1.05])
                plt.ylim([-0.05, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('Receiver operating characteristic: '+name)
                plt.legend(loc="lower right", prop=legendprop)
                plt.savefig('roc.png')


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
            #print("Average: %.2f%% (+/- %.2f%%)" % (scores.mean()*100.0, scores.std() * 2 * 100.0), file=output)
            print("Best: ROC AUC = %.2f%%\n" % (scores.max()), file=output, flush=True)
            bestPerforming = numpy.argmax(scores)
            bestClassifier = classifiers[bestPerforming]
            bestClassifier.cvindex = bestPerforming

            with (clfdir / 'best.pickle').open('wb') as bestSer:
                pickle.dump(bestClassifier, bestSer)



            #results.append(bestClassifier)
            yield bestClassifier

if __name__ == "__main__":
    cvdir = Path(sys.argv[1])
    for best in learningTest(cvdir):
        print(best.name, "done.")
