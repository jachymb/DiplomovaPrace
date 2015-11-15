from pomegranate import DiscreteDistribution, ConditionalProbabilityTable, JointProbabilityTable, State, BayesianNetwork
from dp.learn import POSTIVE_LABEL, NEGATIVE_LABEL
import numpy
from itertools import product
import copy

PRIOR = 1

class BayesNet:
    def __init__(self, fold, clfName, ontology):
        self.ontology = ontology
        self.fold = fold
        self.clfName = clfName
        #self.topsort = ontology.termsByDepth()
        self.network = BayesianNetwork(clfName.replace(" ","_"))
        self.lenValidation = None

    @staticmethod
    def chname(name):
        return name
        #return name.replace(" ","_").replace(":","")

    def generateCPD(self, term, clf, X_train, y_train, X_test, y_test, X_validation, y_validation, g_train, g_test, g_validation):
        root = self.ontology.root
        
        posTrain = sum(y_train == POSTIVE_LABEL)
        negTrain = sum(y_train == NEGATIVE_LABEL)
        totalTrain = posTrain + negTrain

        children = sorted(self.ontology[term]['children'])
        parents = sorted(self.ontology[term]['parents'])

        childNodes = [self.ontology[child]['node'][self.fold][self.clfName] for child in children]

        labels = {l : PRIOR for l in product(*(('0','1'),)*(len(children)+1))}
        if children:
            for gene,y in zip(g_train, y_train):
                event = []
                for child in children:
                    event.append('0' if gene in self.ontology.associations[child] else '1')
                event.append('0' if gene in self.ontology.associations[term] else '1')
                assert (gene in self.ontology.associations[term]) == (y == POSTIVE_LABEL)
                event = tuple(event)

                labels[event] += 1
            cprior = PRIOR * (2 ** len(children))
            hidden = ConditionalProbabilityTable(
                        [
                            list(event) + [counted/(cprior+(posTrain if event[-1] == '0' else negTrain))] # if term != root else {'0':1.,'1':0.}[event[0]]]
                            #list(event) + [counted/total] # if term != root else {'0':1.,'1':0.}[event[0]]]
                            for event, counted in labels.items()],
                        [node.distribution for node in childNodes])

        else: #No children
            hidden = DiscreteDistribution({'0': posTrain / totalTrain, '1': negTrain / totalTrain})

        #print("Hidden node %s:" % term)
        #print(hidden)
        #hidden.freeze()
        hidden = State(hidden, name=self.chname("H"+term))

        clf.conf += PRIOR
        posTest, negTest = numpy.sum(clf.conf, 1) 
        
        observed = ConditionalProbabilityTable([
                ['0', '0', clf.conf[0][0] / posTest], # if term != root else 1.],
                ['1', '0', clf.conf[0][1] / posTest], # if term != root else 0.],
                ['0', '1', clf.conf[1][0] / negTest], # if term != root else 0.],
                ['1', '1', clf.conf[1][1] / negTest]], #if term != root else 1.]],
            [hidden.distribution])
        #observed.freeze()
        #print("Observed node %s:" % term)
        #print(observed)
        observed = State(observed, name=self.chname(term))

        self.network.add_states([hidden, observed])
        self.network.add_transition(hidden, observed)
        for child in children:
            childNode = self.ontology[child]['node'][self.fold][self.clfName]
            self.network.add_transition(childNode, hidden)


        self.ontology[term]['node'][self.fold][self.clfName] = hidden
        self.ontology[term]['clf'][self.fold][self.clfName] = clf, X_validation, y_validation, g_validation
        assert self.lenValidation is None or self.lenValidation == len(y_validation)
        self.lenValidation = len(y_validation)
    
    def bake(self):
        self.network.bake()

    def predict(self):
        #print(self.network.graph)
        self.predictions = {term : numpy.empty((self.lenValidation, 2), dtype=float) for term in self.ontology.ontology}

        classifiers = {term : self.ontology[term]['clf'][self.fold][self.clfName] for term in self.ontology.ontology}
        #for term, (clf,X,y,g) in classifiers.items():
        #    print(term, ":", repr(self.clfName), repr(clf.name), self.fold, clf.fold)

        observations = {
                self.chname(term) : clf.predict(X)
                for term, (clf, X, y, g)
                in classifiers.items()}
        print("observations:")
        print(observations)
        gt = {
                self.chname(term) : y
                for term, (clf, X, y, g)
                in classifiers.items()}
        print("gt:")
        print(gt)

        for i in range(self.lenValidation):
            observation = {term : str(pred[i]) for term,pred in observations.items()}
            #print("Observation for gene %d" % i)
            #print(observation)
            #print(self.network.forward_backward(observation))
            beliefs = self.network.predict_proba(observation) 
            for state, belief in zip(self.network.states, beliefs):
                if not state.name.startswith('H'): continue

                term = state.name[1:]
                distribution = belief.parameters
                #print(distribution)
                self.predictions[term][i][0] = distribution[0]['0']
                self.predictions[term][i][1] = distribution[0]['1']
    
        print("predictions:")
        print(self.predictions)

    def nodeAsClf(self, term):
        class BayesCLF:
            def predict_proba(X_test):
                classifier = self.ontology[term]['clf'][self.fold][self.clfName]
                res = numpy.empty((self.lenValidation, 2))
                for i in range(self.lenValidation):
                    pass

        return BayesCLF()
                
