from bayespy.nodes import Categorical, Mixture, Gaussian
from bayespy.inference.vmp.nodes.stochastic import Stochastic
from bayespy.inference import VB
from dp.learn import POSTIVE_LABEL, NEGATIVE_LABEL
from dp.utils import loadClf, getTermPath
from itertools import product
import numpy
import copy
import sys
import random

PRIOR = 1

def _or(p_false, p_true):
    """
    Build probability table for OR-operation of two parents

    p_false: Probability table to use if both are FALSE

    p_true: Probability table to use if one or both is TRUE
    """
    return np.take([p_false, p_true], [[FALSE, TRUE], [TRUE, TRUE]], axis=0)

class BayesNet:
    def __init__(self, fold, clfName, ontology):
        self.ontology = ontology
        self.fold = fold
        self.clfName = clfName
        #self.topsort = ontology.termsByDepth()
        self.lenValidation = None
        self.allobserved = {}
        self.allhidden = {}
        self.extranodes = set()

    def generateCPD(self, term): #, X_train, y_train, X_test, y_test, X_validation, y_validation, g_train, g_test, g_validation):
        
        clf = loadClf(self.ontology[term]['name'], self.fold, self.clfName)
        posTrain = sum(clf.y_train == POSTIVE_LABEL)
        negTrain = sum(clf.y_train == NEGATIVE_LABEL)
        totalTrain = posTrain + negTrain

        children = sorted(self.ontology[term]['children'])
        parents = sorted(self.ontology[term]['parents'])


        labels = {l : PRIOR for l in product(*((POSTIVE_LABEL,NEGATIVE_LABEL),)*(len(children)+1))}
        if children:
            childNodes = [self.ontology[child]['node'][self.fold][self.clfName] for child in children]
            for gene,y in zip(clf.g_train, clf.y_train):
                event = []
                for child in children:
                    event.append(POSTIVE_LABEL if gene in self.ontology.associations[child] else NEGATIVE_LABEL)
                event.append(POSTIVE_LABEL if gene in self.ontology.associations[term] else NEGATIVE_LABEL)
                assert (gene in self.ontology.associations[term]) == (y == POSTIVE_LABEL)
                event = tuple(event)

                labels[event] += 1
            def countBoth(event):
                return labels[event[:-1]+(POSTIVE_LABEL,)] + labels[event[:-1]+(NEGATIVE_LABEL,)]
            cprior = PRIOR * (2 ** len(children))
            
            types = [Mixture]*(len(children)-1) + [Categorical]
            mixparams = [i for s in zip(childNodes, types) for i in s]
            cpd = numpy.empty((2,)*(len(children)+1))
            for event, counted in labels.items():
                v=cpd
                for b in event[:-1]:
                    v = v[b]
                    
                hid = event[-1]
                print("Event: ", event)
                if POSTIVE_LABEL not in event[:-1]: # Všichni potomci označeni "ne"
                    v[hid] = counted/countBoth(event)
                    print("Stored %d / %d" % (counted,countBoth(event)))
                else:
                    v[hid] = {POSTIVE_LABEL: 0.99, NEGATIVE_LABEL:0.01}[hid]
                    print("Stored %d : %d" % (hid, v[hid]))


            #print(term,"<-",",".join(children))
            print(cpd)
            #print(labels)

            hidden = Mixture(*mixparams, cpd)
            hidden.params = cpd
            

        else: #No children
            #hidden = DiscreteDistribution({'0': posTrain / totalTrain, '1': negTrain / totalTrain})
            params = (posTrain / totalTrain, negTrain / totalTrain)
            hidden = Categorical(params)
            hidden.params = params

        #print("Hidden node %s:" % term)
        #print(repr(hidden))
        #print([p for p in hidden.parents if isinstance(p, Stochastic)])
        #print(hidden.get_moments())

        conf = clf.conf + PRIOR
        #posTest, negTest = numpy.sum(conf, 1) 
        posTest, negTest = numpy.sum(conf, 0) 
        #print("Confusion matrix:")
        #print(conf)
       
        try:
            assert term != self.ontology.root
            pos_decisions = clf.decision_function(clf.X_test[clf.y_test==POSTIVE_LABEL])
            neg_decisions = clf.decision_function(clf.X_test[clf.y_test==NEGATIVE_LABEL])
            means = [numpy.mean(pos_decisions)], [numpy.mean(neg_decisions)]
            maxprec = 100.0
            precs = [[numpy.min((1/numpy.var(pos_decisions), maxprec))]], [[numpy.min((1/numpy.var(neg_decisions), maxprec))]]
        #else:
        except (ValueError, AssertionError):
            means = [-1.], [1.]
            precs = [[1.]], [[1.]]
        print("Gaussian params:", term, self.ontology[term]['name'], means, precs)
        observed = Mixture(hidden, Gaussian, means, precs)
        #observed = ConditionalProbabilityTable([
        #        ['0', '0', conf[0][0] / posTest], # if term != root else 1.],
        #        ['0', '1', conf[0][1] / posTest], # if term != root else 0.],
        #        ['1', '0', conf[1][0] / negTest], # if term != root else 0.],
        #        ['1', '1', conf[1][1] / negTest]], #if term != root else 1.]],
        #    [hidden.distribution])
        #print("Observed node %s - %s:" % (term, self.ontology[term]['name']))
        #print(repr(observed))
        #print([p for p in observed.parents if isinstance(p, Stochastic)])

        self.ontology[term]['node'][self.fold][self.clfName] = hidden
        #self.ontology[term]['clf'][self.fold][self.clfName] = clf, X_validation, y_validation, g_validation
        assert self.lenValidation is None or self.lenValidation == len(clf.y_validation)
        self.lenValidation = len(clf.y_validation)
        self.allobserved[term] = observed
        self.allhidden[term] = hidden
        self.extranodes.update((p for p in hidden.parents if isinstance(p, Stochastic)))
    
    def bake(self):
        pass

    def getCopy(self):
        #return copy.deepcopy(self.allhidden), copy.deepcopy(self.allobserved), copy.deepcopy(self.extranodes)
        return copy.deepcopy((self.allhidden, self.allobserved, self.extranodes))

    def predict(self):
        #print(self.network.graph)
        self.predictions = {term : numpy.empty((self.lenValidation, 2), dtype=float) for term in self.ontology.ontology}

        classifiers = {
                term : lambda: loadClf(self.ontology[term]['name'], self.fold, self.clfName) #self.ontology[term]['clf'][self.fold][self.clfName]
                for term
                in self.ontology.ontology
                }
        #for term, (clf,X,y,g) in classifiers.items():
        #    print(term, ":", repr(self.clfName), repr(clf.name), self.fold, clf.fold)

        observations = {
                term : clf.decision_function(clf.X_validation) if term != self.ontology.root else numpy.array([-1.]*len(clf.X_validation))
                for term, clff #(clf, X, y, g)
                in classifiers.items()
                for clf in (clff(),)
                }
        #print("observations:")
        #print(observations)
        gt = {
                term : clf().y_validation
                for term, clf
                in classifiers.items()}
        #print("gt:")
        #print(gt)

        for i in range(self.lenValidation):
            observation = {term : pred[i] for term, pred in observations.items()}
            #print("Observation for gene %d" % i)
            #print(observation)
            #print(self.network.forward_backward(observation))
            hidden, observed, extra = self.getCopy()
            #print(hidden)
            #print(observed)
            #print(observation)
            #for term, node in hidden.items():
            #    print(i, term, node.get_moments()[0])
            for k,v in observation.items():
                observed[k].observe((v,))
            #    print("%s observes %s" % (k, v))
            allv = (*hidden.values(), *observed.values(), *extra)
            #print([(x, [p for p in x.parents if isinstance(p, Stochastic)]) for x in [*hidden.values(), *observed.values()]])
            Q = VB(*allv)
            Q.update(*allv, tol=1e-7, repeat=1000, verbose=True)
            #print("---")

            for term, node in hidden.items():
                #print(i, term, node.get_moments()[0])
                self.predictions[term][i,:] = node.get_moments()[0]

        #print("predictions:")
        #print(self.predictions)
        for term in observations:
                compare = numpy.empty((len(gt[term]),4), dtype=float)
                compare[:,0] = gt[term]
                compare[:,1] = observations[term]
                compare[:,2] = self.predictions[term][:,1]
                compare[:,3] = numpy.round(self.predictions[term][:,1])
                print(term, self.ontology[term]['name'])
                print(compare)

    def nodeAsClf(self, term):
        clf = loadClf(self.ontology[term]['name'], self.fold, self.clfName)
        #clf.predict_proba = lambda *a: self.predictions[term] #HACK
        #clf.__dict__['predict_proba'] = lambda *a: self.predictions[term] #HACK
        class FakeClassifier:
            def predict_proba(*a): 
                return self.predictions[term]
        for a in ('X_test','y_test',):
        #    if a != 'predict_proba':
            setattr(FakeClassifier, a, getattr(clf, a))
        return FakeClassifier()
                
