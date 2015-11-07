import sys
#import progressbar
import random
import traceback
import json
from dp.gene import GeneFactory
from dp.utils import debug, parallel_map_dill
from dp.treeliker import TreeLikerWrapper
from collections import defaultdict
from itertools import groupby, product
from pathlib import Path
from pomegranate import DiscreteDistribution, ConditionalProbabilityTable, State, BayesianNetwork
import dp.utils

__all__ = ["Onotology"]

class Ontology:
    """ Class representing the Gene Onotology graph."""

    def __init__(self, inputFileName, namespace = 'molecular_function'):
        """Constructor, reads and parses the ontology OBO file."""
        debug("Reading ontology file %s... " % inputFileName)
        self.root = None
        self.namespace = namespace
        ontology = defaultdict(lambda: defaultdict(list))
        with open(inputFileName) as go:
            terms = groupby(go.read().splitlines(), lambda x: x != '')

            for b, term in terms:
                term = list(term)
                if not b or term[0] != '[Term]': continue
                nonlists = ('id', 'def', 'name', 'namespace', 'is_obsolete')
                # Do some deciphering here...
                term = defaultdict(list, [
                        (a, [y[1] for y in b][0 if a in nonlists else slice(None)])
                        for a,b in groupby(
                            [x.split(': ', 1) for x in term[1:]],
                            lambda x: x[0])])

                # Filter terms by namespace, discard obsolete terms
                if term['namespace'] != namespace or term['is_obsolete'] == 'true':
                    continue

                # Decide root node
                if term['name'] == namespace:
                    assert self.root is None
                    self.root = term['id']

                # Save the term to ontology
                ontology[term['id']]['name'] = term['name']
                for ref in term['is_a']:
                    refid, refname = ref.split(' ! ')
                    ontology[refid]['children'].append(term['id'])
                    ontology[term['id']]['parents'].append(refid)

        self.ontology = ontology

        self.associations = None
        self.geneFactory = GeneFactory()
        debug("Initialized ontology for file %s... " % inputFileName)

    def genTranslations(self):
        return {self.ontology[term]['name'] : term for term in self.ontology}

    def setAssociations(self, associations, attrname = 'associations'):
        associations.translations = self.genTranslations()
        associations.ontology = self
        self.__setattr__(attrname, associations)

    def _termDepth(self, term):
        #if term == self.root: return 0
        parents = self.ontology[term]['parents']
        if len(parents) == 0: return 0
        return 1 + min(map(self._termDepth, parents))

    def deleteSmallTerms(self, lb):
        """ Delete terms which don't have enough associations.
        Do not call this prior to calling transitiveClosure()"""
        notEnough = [term for term in self.ontology if len(self.associations[term]) < lb]
        notEnough.sort(key = lambda x: -self._termDepth(x)) # Yes, it's necessary.
        for term in notEnough:
            for parent in self.ontology[term]['parents']:
                if term in self.ontology[parent]['children']:
                    self.ontology[parent]['children'].remove(term)
            del self.ontology[term]
            del self.associations.associations[term]
            if hasattr(self, 'reserved'):
                del self.reserved.associations[term]
        debug("Deleted terms with not enough associations. Left %d terms." % len(self.ontology))

    def jsonExport(self, output = sys.stdout):
        """Export generated onotology in JSON"""
        json.dump(self.ontology, output)

    def dotExport(self, direction = "parents", output = sys.stdout):
        """ Export as a graph to the DOT format for graphical presentation."""
        print("digraph ontology {", file=output)
        def nodename(t): return t.replace(":", "")
        for term, props in self.ontology.items():
            print('%s [label="%s"]' % (nodename(term), props["name"]), file=output)

        for term, props in self.ontology.items():
            for related in props[direction]:
                print('%s -> %s' % (nodename(term), nodename(related)), file=output)
        print("}", file=output)

    def __iter__(self):
        """Iterate over the underlying dict."""
        return self.ontology.__iter__()

    def __getitem__(self, item):
        """An instance of this class can be indexed by GO terms or their names."""
        if item in self.ontology:
            return self.ontology[item]

        term = self.getTermByName(item)
        if term in self.ontology:
            return self.ontology[term]

        raise KeyError(item)

    def display(self, lb = 0):
        for term in self.ontology:
            numTerms = len(self.associations[term])
            children = ", ".join(self.ontology[term]['children'])
            if numTerms >= lb:
                #print(term, self.ontology[term]['name'], numTerms, "->", children)
                print(term, self.ontology[term]['name'], numTerms)

    def generateExamples(self, term, output, maxAssociations = None):
        """Generates examples from genes associated with the term in logical reprentation in the pseudo-prolog syntax. """
        def getRecord(geneName):
            try:
                gene = self.geneFactory.getGene(geneName)
                e = '"%s" %s' % (term, ", ".join(gene.logicalRepresentation()))
                assert e != ''
                print(e, file=output, flush=True)
                #if pbar is not None:
                #    pbar.update(pbar.currval + 1)
            except Exception as exc:
                traceback.print_exc()

        genes = sorted(self.associations[term])
        genes = list(genes)
        random.shuffle(genes)
        genes = genes[:maxAssociations]
        for gene in genes:
            getRecord(gene)
        #with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            #executor.map(getRecord, genes)
            #for future in concurrent.futures.as_completed(s):
            #    debug("Future completed. " + str(future))
            #    if future.exception() is not None:
            #        debug("Exception: " + future.exception())

    def generateDataset(self, term, output, maxPositive = None, maxNegative = None, associations = None):
        """Generate whole dataset directly usable for learning. The terms are used as learning classes."""
        # FIXME: Remove maxNegative parameter
        print(len(self.associations.associations))
        print(term)
        print(len(self.associations[term]))
        if associations is None:
            associations = self.associations

        totalPos = len(associations[term]) if maxPositive is None else min(maxPositive, len(associations[term]))
        maxNegative = round(totalPos / associations.getRatio(term))
        totalNeg = len(associations['~'+term]) if maxPositive is None else min(maxNegative, len(associations['~'+term]))
        total = totalPos + totalNeg
        debug("Generating dataset for term: %s. Using %d postive and %d negative examples." % (self[term]['name'], totalPos, totalNeg))
        #if dp.utils.verbosity >= 2:
            #pbar = progressbar.ProgressBar(maxval=total, widgets = (
            #    progressbar.Bar(), ' ', progressbar.Counter(), '/'+str(total), ' =', progressbar.Percentage()))
            #pbar.start()
        #else:
            #pbar = None
        self.generateExamples(    term, output, maxPositive)
        self.generateExamples('~'+term, output, maxNegative)
        #if dp.utils.verbosity >= 2:
        #    pbar.finish()
        debug("Finished generating dataset for term: %s" % (self[term]['name']))

    def getTermByName(self, name):
        """Returns human-readable name of the given GO term."""
        for k,v in self.ontology.items():
            if v['name'] == name:
                return k
        raise KeyError('Name %s is not in the ontology!' % name)

    def completeTest(self, maxPositive, maxNegative, treelikerPath, template, processes = 1):
        bestClassifiers = []
        terms = sorted(self.ontology.keys(), key = lambda x: (-self._termDepth(x), x))[1:3] # This sorting is needed later in bnet learning
        treeliker = TreeLikerWrapper(self, treelikerPath, template)
        treeliker.runValidation(self.reserved)
        sys.exit()
        def processTerm(term):
            return treeliker.runTermTest(term, maxPositive, maxNegative)
            
        parallel_map_dill(processes, processTerm, terms)

        #for term in terms:
        #    if term == self.root:
        #        continue

        #    best = treeliker.runTermTest(term, maxPositive, maxNegative)
        #    bestClassifiers.append(best)
        #    print("Best:",best)
            
        #bnet = self._toBayessNet(bestClassifiers, terms)

    def _toBayessNet(self, classifiers, terms):
        assert hasattr(self, 'reserved')

        bnet = BayesianNetwork(self.ontology[self.root]['name'])
        hidden_nodes = {}
        observed_nodes = {}
        classes = [('1','0')]
        for term in terms:
            children = sorted(self.ontology[term]['children'])
            cpd_hidden = [[*x, 1.0] for x in product(*classes*(len(children)+2))]
            node_hidden = ConditionalProbabilityTable(cpd_hidden, [hidden_nodes[child] for child in children])

            cpd_observed = [[*x, 1.0] for x in product(*classes*2)]
            node_observed = ConditionalProbabilityTable(cpd_observed, [node_hidden])
            hidden_nodes[term] = node_hidden
            observed_nodes[term] = node_observed
            print(term, "observed:")
            print(node_observed)
            print(term, "hidden:")
            print(node_hidden)
                
            #state = State(node, name = self.ontology[term]['name'])
            #bnet.add_state(state)



        return
        genes = [*self.reserved.dataset]
        #terms = [*self.ontology.keys()]
        labels = [
                [(t in assocTerms) for t in terms]
                for gene, assocTerms
                in ((gene, tuple(self.reserved.inverseAssoc(gene))) for gene in  genes)]
        labels = np.array(labels, dtype=bool)
        for i in range(8):
            X_train, X_test, y_train, y_test = cross_validation.train_test_split(genes, labels)
            clfs = zip(*classifiers)
            for clftype in clfs:
                pass
                # Potřebujeme další validační množinu označenou treelikerem
