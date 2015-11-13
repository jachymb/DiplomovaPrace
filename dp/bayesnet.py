from pomegranate import DiscreteDistribution, ConditionalProbabilityTable, State, BayesianNetwork
from dp.learn import POSTIVE_LABEL, NEGATIVE_LABEL
import numpy
from itertools import product
import copy

PRIOR = 1

class BayesNet:
    def __init__(self, clfName, ontology):
        self.ontology = ontology
        self.clfName = clfName
        #self.topsort = ontology.termsByDepth()
        self.network = BayesianNetwork(clfName.replace(" ","_"))

    def generateCPD(self, term, clf, X_train, y_train, X_test, y_test):
        root = self.ontology.root
        
        posTrain = sum(y_train == POSTIVE_LABEL)
        negTrain = sum(y_train == NEGATIVE_LABEL)
        totalTrain = posTrain + negTrain

        children = sorted(self.ontology[term]['children'])
        if children is not None:
            childNodes = [self.ontology[child]['node'][self.clfName] for child in children]

            labels = {l : PRIOR for l in product(*(('0','1'),)*(len(children)+1))}
            for gene in self.ontology.genes:
                event = []
                for child in children:
                    event.append('0' if child in self.ontology.associations[term] else '1')
                event.append('0' if term in self.ontology.associations[term] else '1')
                event = tuple(event)

                labels[event] += 1

            total = len(self.ontology.genes) + len(labels)*PRIOR

            hidden = ConditionalProbabilityTable(
                        [
                            list(event) + [counted/total if term != root else {'0':1.,'1':0.}[event[0]]]
                            for event, counted in labels.items()],
                        [node.distribution for node in childNodes])


        else: #No children
            hidden = DiscreteDistribution({'0': posTrain / totalTrain, '1': negTrain / totalTrain})
        hidden.freeze()
        hidden = State(hidden, name=term.replace(" ","_").replace(":",""))
        print("Hidden node:", hidden)

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
        observed = State(observed, name=(term+"_prediction").replace(" ","_").replace(":",""))
        print("Observed node:", observed)

        self.network.add_states([hidden, observed])
        self.network.add_transition(hidden, observed)
        #self.network.add_transition(observed, hidden)

        self.ontology[term]['node'][self.clfName] = hidden
    
    def bake(self):
        self.network.bake()

    def setClassifiers(self,clfs):
        # assign clasifiers to nodes
        pass

    def predict(self):
        pass
