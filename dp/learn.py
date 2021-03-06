#!/usr/bin/env python

from collections import Counter, defaultdict
from dp.utils import NUM_FOLDS, debug, loadClf, rerun
from pathlib import Path
from scipy import interp
from itertools import product
from sklearn.covariance import EllipticEnvelope
from sklearn.dummy import DummyClassifier
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, BaggingClassifier
from sklearn import cross_validation
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score, recall_score, precision_score, confusion_matrix, precision_recall_curve, roc_curve, auc, average_precision_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler, Normalizer, MaxAbsScaler, RobustScaler
from sklearn.semi_supervised import LabelPropagation, LabelSpreading
from sklearn.svm import SVC
import bz2
import csv
import matplotlib.pyplot as plt
import numpy
import dill
import re
import sys

__all__ = ['readArff', 'learningTest']

VALIDATION_RATIO = 0.5
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

    #with filename.open() as f:
    with bz2.open(str(filename)+'.bz2', 'r') as f:
        
        line = ''
        while line != '@data':
            line = f.readline().decode().strip()
            if line.startswith("@attribute 'classification'"):
                line = line[line.find('{') + 1:line.find('}')]
                classes = {i:n for n,i in enumerate(parseLine(line))}

        for line in f.read().decode().splitlines():
            record = list(parseLine(line))
            labels.append(classes[record[-1]])
            data.append([int(x) for x in record[:-1]])
    return numpy.array(data, dtype=float), numpy.array(labels), classes

# sort of hack
def getGenes(cvdir):
    path = str(cvdir /'dataset.txt.bz2')
    with bz2.open(path, 'rb') as f:
        genes = [re.sub(r'.*?proteinName\((.*?)\).*',r'\1', x.decode().strip()) for x in f]
    fold = 0
    def getIndexes(s):
        return {int(x) for x in re.search("([0-9]+,)+[0-9]+",s).group().split(",")}
    with (cvdir / 'batch.treeliker').open() as f:
        for line in f:
            if line.startswith('set(output,'):
                trainIndices = getIndexes(next(f))
                testIndices = getIndexes(next(f))
                yield [x for i,x in enumerate(genes) if i in trainIndices], [x for i,x in enumerate(genes) if i in testIndices]
                fold += 1
                
legendprop = {'size': 10}
def plotRoc(term, clfName, title, clfs = None):

    mean_tpr = 0.0
    mean_fpr = numpy.linspace(0, 1, 100)
    all_tpr = []
    #plt.clf()
    if clfs is None:
        clfs = (loadClf(term, fold, clfName) for fold in range(NUM_FOLDS))
    for i,clf in enumerate(clfs):
        print("roc", clf)
    #for i, (clf, X_train, y_train, X_test, y_test, X_validation, y_validation,_,_,_) in enumerate(folds):

        probabs = clf.predict_proba(clf.X_test)
        try:
            fpr, tpr, _ = roc_curve(clf.y_test, probabs[:, 1])
        except IndexError:
            fpr, tpr, _ = roc_curve(clf.y_test, probabs[:, 0], pos_label=POSTIVE_LABEL)
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        fpr = numpy.nan_to_num(fpr)
        try:
            roc_auc = auc(fpr, tpr)
        except ValueError: # root node
            roc_auc = 1.

        clf.roc_auc = roc_auc
        plt.plot(fpr, tpr, lw=1, label='Fold %d, AUC = %0.2f' % (i, roc_auc))

    mean_tpr /= NUM_FOLDS
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    plt.plot(mean_fpr, mean_tpr, 'k--',
            label='Mean, AUC = %0.2f' % mean_auc, lw=2)

    plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6))
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic: '+title)
    plt.legend(loc="lower right", prop=legendprop)
    #plt.savefig(str(outdir/(clfName.replace(" ","_")+'_roc.png')))

def plotPrc(clfName, folds, outdir):
    y_tests = []
    y_scores = []
    plt.clf()
    for i, (clf, X_test, y_test, _, _, _, _,_,_,_) in enumerate(folds):
        try:
            y_score = clf.decision_function(X_test)
        except AttributeError:
            y_score = clf.predict_proba(X_test)[:, 0]
        precision, recall, _ = precision_recall_curve(y_test, y_score, pos_label=POSTIVE_LABEL)
        y_tests.extend(y_test)
        y_scores.extend(y_score)

        try:
            area = average_precision_score(y_test, y_score)
        except ValueError:
            area = 0.0
        clf.prc_auc = area
        plt.plot(recall, precision, label='Fold %d, AUC = %0.2f' % (i, area), lw=1)

    precision, recall, _ = precision_recall_curve(y_tests, y_scores, pos_label=POSTIVE_LABEL)
    try:
        area = average_precision_score(y_tests, y_scores)
    except ValueError:
        area = 0.0
    plt.plot(recall, precision, 'k--', label='Mean, AUC = %0.2f' % (area), lw=2)

    plt.title('Precision-Recall: %s\n%s'%(clfName,outdir.name.replace("_"," ")))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.legend(loc="lower center", prop=legendprop)

    plt.savefig(str(outdir/(clfName.replace(" ","_")+'_precision-recall.png')))

def plotPCA(X_train, y_train, X_test, y_test, outdir):
    #clf = loadClf(term, fold, clfName)
    #try:
    #    decision = clf.decision_function
    #    Vf = numpy.arange(-1.,1.1,0.1)
    #    V = (0.,)
    #except AttributeError:
    #    decision =  lambda x:clf.predict_proba(x)[:,0]
    #    Vf = numpy.arange(0.,1.05,0.05)
    #    V = (0.5,)
    scaler = MaxAbsScaler(copy=False)
    target_names = ("Positive","Negative")
    term = outdir.parent.name.replace("_", " ")
    pca = PCA(n_components=2)
    pca.fit(X_train)
    scaler.fit(pca.transform(X_train))
    #delta = 0.025
    #a=numpy.arange(-1., 1., delta)
    #b=numpy.arange(-1., 1., delta)
    #A,B = numpy.meshgrid(a,b)
    #C=numpy.empty(A.shape)
    for X, y, n in ((X_train, y_train, 'training'), (X_test, y_test, 'testing')):
        X_r = scaler.transform(pca.transform(X))
        inlier = (numpy.abs(X_r[:,0]) <= 1) & (numpy.abs(X_r[:,1]) <= 1)
        #print(X_r)
        plt.clf()

        #for k,l in product(range(len(a)),range(len(b))):
        #    C[k][l] = decision(pca.inverse_transform(scaler.inverse_transform(((A[k][l],B[k][l]),))))
        #print(C)
        #cfp = plt.contourf(A,B,C,Vf,cmap=plt.cm.bone)
        #cfp.cmap.set_under('black')
        #cfp.cmap.set_over('white')
        #plt.contour(A,B,C,V,colors=("b",))
        #y=clf.predict(X)
        for c, i, target_name in zip("rg", (0, 1), target_names):
            plt.scatter(X_r[(y == i) & inlier, 0], X_r[(y == i) & inlier, 1],
                    c = c,
                    label = target_name,
                    marker = ",",
                    s = 1,#0.8,#1/numpy.sqrt(2),
                    #edgecolors='none',
                    linewidth = 0,
                    alpha = 0.7)
        plt.legend()
        plt.title('PCA for %s on %s data' % (term, n))
        plt.savefig(str(outdir/('pca-%s.png' % (n,))))
        plt.savefig(str(outdir/('pca-%s.ps' % (n,))))

#def plotLDA(X_train, X_test, y_train, y_test, outdir):
    #target_names = ["pos","neg"]
    #lda = LDA(n_components=2)
    #X_r2 = lda.fit(X_train, y_train).transform(X_test)
    #plt.clf()
    #for c, i, target_name in zip("rg", [0, 1], target_names):
    #    plt.scatter(X_r2[y_test == i, 0], X_r2[y_test == i, 0], c=c, label=target_name)
    #plt.legend()
    #plt.title('LDA')
    #plt.savefig(str(outdir/'lda.png'))

def learningTest(cvdir):
    debug("Starting learning in node %s." % cvdir.name) 
    clasfifiers = (
            #TOOD, use third dict for params
            #("Bagged SVM", lambda: BaggingClassifier(SVC(C=0.1,probability=True))),
            #("LabelPropagation RBF", LabelPropagation),
            #("LabelSpreading RBF", LabelSpreading),
            #("LabelSpreading-7nn", lambda: LabelSpreading(kernel='knn')),
            #("LabelPropagation-7nn", lambda: LabelPropagation(kernel='knn')),
            #!("AdaBoost-DecisionTree", AdaBoostClassifier),
            #("5-NN", lambda: KNeighborsClassifier(p=1, algorithm='kd_tree')),
            #!("Random Forest", RandomForestClassifier),
            #("SGD", lambda: SGDClassifier(n_iter=100,alpha=0.01,loss="modified_huber")),
            ("RBF SVM C=1", lambda: SVC(shrinking=False, tol=1e-5,probability=True)),
            #("RBF SVM C=0.5", lambda : SVC(C=0.1,probability=True)),
            #("RBF SVM C=2", lambda : SVC(C=10,probability=True)),
            #("RBF SVM C=inf", lambda : SVC(C=numpy.inf,probability=True)),
            #("Linear SVM C=1", lambda: SVC(kernel='linear',probability=True)),
            #("Quadratic SVM C=1", lambda: SVC(kernel='poly', degree=2,probability=True)),
            #("Cubic SVM C=1", lambda: SVC(kernel='poly', degree=3,probability=True)),
            )

    scaler = StandardScaler(copy=False)
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

        #alldata = defaultdict(list)
        alldata = []
        for i, (g_train, g_test) in zip(range(NUM_FOLDS), getGenes(cvdir)) :
            foldDir = cvdir / str(i)
            train = foldDir / 'train.arff'
            test = foldDir / 'test.arff'


            X_train, y_train, _ = readArff(train)
            X_test , y_test,  _ = readArff(test)
            assert len(g_train) == len(y_train) and len(g_test) == len(y_test)

            # Preprocess
            #envelope = EllipticEnvelope(contamination=0.05)
            #envelope.fit(X_train)
            #inliers = envelope.predict(X_train) == 1

            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)
            norms = [numpy.linalg.norm(x) for x in X_train]
            #print(numpy.mean(norms))
            #print(numpy.median(norms))
            #print(numpy.max(norms))
            #print(numpy.min(norms))
            #print(numpy.sqrt(numpy.cov(norms)))

            plotPCA(X_train, y_train, X_test, y_test, foldDir)
            splitIndex = round(len(y_test)*VALIDATION_RATIO)
            X_validation, y_validation, g_validation = X_test[:splitIndex], y_test[:splitIndex], g_test[:splitIndex]
            X_test, y_test, g_test = X_test[splitIndex:], y_test[splitIndex:], g_test[splitIndex:]
            assert len(g_train) == len(y_train) and len(g_test) == len(y_test)

            #plotLDA(X_train, X_test, y_train, y_test, foldDir)

            for name, Clf in clasfifiers:
                alldata.append(( name, i))
                serFile = foldDir / (name + ".pickle.bz2")
                if serFile.is_file() and not rerun:
                    continue

                debug("Fitting clasifier %s for fold %d of %d in node %s." % (name, i+1, NUM_FOLDS, cvdir.name))
                if cvdir.name != 'molecular_function':
                    clf = Clf()
                else:
                    clf = DummyClassifier(strategy='constant', constant=POSTIVE_LABEL)
                    clf.decision_function = lambda a: [1.0]*len(a)

                print('Testing the classifier %s:' % name, file=output)
       
                if clf.__module__.startswith('sklearn.semi_supervised'):
                    y_train = - y_train
                    y_test = - y_test

                pos = (y_train == POSTIVE_LABEL)
                neg = (y_train == NEGATIVE_LABEL)
                posWeight = numpy.sum(neg) / len(y_train)
                negWeight = numpy.sum(pos) / len(y_train)
                sample_weight = posWeight*pos + negWeight*neg

                try:
                    #clf.fit(X_train, y_train)
                    clf.fit(X_train, y_train, sample_weight = sample_weight)
                except TypeError:
                    clf.fit(X_train, y_train)

                y_pred = clf.predict(X_test)
                #compare = numpy.empty((len(y_pred),2))
                #compare[:,0] = y_test
                #compare[:,1] = y_pred
                #print(compare)
                pos = (y_test == POSTIVE_LABEL)
                neg = (y_test == NEGATIVE_LABEL)
                #prec = precision_score(y_test, y_pred, average='weighted')
                #reca = recall_score(y_test, y_pred, average='weighted')
                score = accuracy_score(y_test[pos], y_pred[pos])*posWeight + accuracy_score(y_test[neg], y_pred[neg])*negWeight
                conf = confusion_matrix(y_test, y_pred)
                if len(conf) == 1:
                    conf = numpy.array([[conf[0][0], 0],[0,0]])
                print("Fold %d score: %.2f, confusion matrix: %s" % (i, score*100.0, conf.tolist()), file=output)
                clf.conf = conf
                clf.fold = i
                clf.name = name
                clf.cvdir = cvdir
                clf.X_train = X_train
                clf.X_test = X_test
                clf.X_validation = X_validation
                clf.y_train = y_train
                clf.y_test = y_test
                clf.y_validation = y_validation
                clf.g_train = g_train
                clf.g_test = g_test
                clf.g_validation = g_validation

                #alldata[name].append((clf, X_train, y_train, X_test, y_test, X_validation, y_validation, g_train, g_test, g_validation))

                with bz2.open(str(serFile), 'wb') as ser:
                    ser.write(dill.dumps(clf))
                
                debug("Finished fitting clasifier %s for fold %d of %d in node %s." % (name, i+1, NUM_FOLDS, cvdir.name))

       
        #for clfName, folds in alldata.items():
        #    plotRoc(clfName, folds, cvdir)
        #    plotPrc(clfName, folds, cvdir)

            # Musí to být tady, protože funkce výše mohou objekt ještě upravovat
        #    for i, (clf,_,_,_,_,_,_,_,_,_) in enumerate(folds):
        #        foldDir = cvdir / str(i)
        #        with (foldDir / (name.replace(' ','')+'.pickle')).open('wb') as ser:
        #            pickle.dump(clf, ser)
    debug("Finished learning in node %s." % cvdir.name) 
    return alldata
            
if __name__ == "__main__":
    cvdir = Path(sys.argv[1])
    learningTest(cvdir)
