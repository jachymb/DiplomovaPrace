from pomegranate import DiscreteDistribution, ConditionalProbabilityTable, State, BayesianNetwork
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
        return name.replace(" ","_").replace(":","")

    def generateCPD(self, term, clf, X_train, y_train, X_test, y_test, X_validation, y_validation, g_train, g_test, g_validation):
        root = self.ontology.root
        
        posTrain = sum(y_train == POSTIVE_LABEL)
        negTrain = sum(y_train == NEGATIVE_LABEL)
        totalTrain = posTrain + negTrain

        children = sorted(self.ontology[term]['children'])

        childNodes = [self.ontology[child]['node'][self.clfName] for child in children]

        labels = {l : PRIOR for l in product(*(('0','1'),)*(len(children)+1))}
        for gene in g_train:
            event = []
            for child in children:
                event.append('0' if child in self.ontology.associations[term] else '1')
            event.append('0' if term in self.ontology.associations[term] else '1')
            event = tuple(event)

            labels[event] += 1

        total = len(g_train) + len(labels)*PRIOR

        hidden = ConditionalProbabilityTable(
                    [
                        list(event) + [counted/total if term != root else {'0':1.,'1':0.}[event[0]]]
                        for event, counted in labels.items()],
                    [node.distribution for node in childNodes])


        #else: #No children
        #    hidden = DiscreteDistribution({'0': posTrain / totalTrain, '1': negTrain / totalTrain})
        hidden.freeze()
        hidden = State(hidden, name=self.chname(term+" hidden"))

        for child in children:
            childNode = self.ontology[child]['node'][self.clfName]
            self.network.add_transition(hidden, childNode)
            #self.network.add_transition(childNode, hidden)

        posTest, negTest = numpy.sum(clf.conf + PRIOR, 1) 
        
        observed = ConditionalProbabilityTable([
                ['0', '0', clf.conf[0][0] / posTest if term != root else 1.],
                ['0', '1', clf.conf[0][1] / posTest if term != root else 0.],
                ['1', '0', clf.conf[1][0] / negTest if term != root else 0.],
                ['1', '1', clf.conf[1][1] / negTest if term != root else 1.]],
            [hidden.distribution])
        observed.freeze()
        observed = State(observed, name=self.chname(term))

        self.network.add_states([hidden, observed])
        self.network.add_transition(hidden, observed)
        #self.network.add_transition(observed, hidden)

        self.ontology[term]['node'][self.clfName] = hidden
        self.ontology[term]['clf'][self.clfName] = clf, X_validation, y_validation, g_validation
        assert self.lenValidation is None or self.lenValidation == len(y_validation)
        self.lenValidation = len(y_validation)
    
    def bake(self):
        self.network.bake()

    def predict(self):
        print(self.network.graph)

        classifiers = {term : self.ontology[term]['clf'][self.clfName] for term in self.ontology.ontology}
        for term, (clf,X,y,g) in classifiers.items():
            print(term, ":", repr(self.clfName), repr(clf.name), self.fold, clf.fold)

        observations = {
                self.chname(term) : clf.predict(X)
                for term, (clf, X, y, g)
                in classifiers.items()}

        for i in range(self.lenValidation):
            observation = {term : str(pred[i]) for term,pred in observations.items()}
            print("Observation for gene %d" % i)
            print(observation)
            print(self.network.forward_backward(observation))
            #print(self.network.predict_proba(observation))
            #print(self.network.predict(observation))
            print(y)
