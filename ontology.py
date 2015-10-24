__all__ = ["Onotology"]

import sys
import progressbar
import random
import traceback
import json
from gene import GeneFactory
from utils import debug
from collections import defaultdict
from itertools import groupby
from pathlib import Path
from treeliker import TreeLikerWrapper


class Ontology:
    """ Class representing the Gene Onotology graph."""

    def __init__(self, inputFileName, namespace):
        """Constructor, reads and parses the ontology OBO file."""
        debug("Reading ontology file %s... " % inputFileName, False)
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
        debug("Done.")

        self.associations = None
        self.geneFactory = GeneFactory()

    def genTranslations(self):
        return {self.ontology[term]['name'] : term for term in self.ontology}

    def setAssociations(self, associations, attrname = 'associations'):
        associations.translations = self.genTranslations()
        associations.ontology = self
        self.__setattr__(attrname, associations)

    def _termDepth(self, term):
        if term == self.root: return 0
        return 1 + min(map(self._termDepth, self.ontology[term]['parents']))

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

    def generateExamples(self, term, pbar, output, maxAssociations = None):
        """Generates examples from genes associated with the term in logical reprentation in the pseudo-prolog syntax. """
        def getRecord(geneName):
            try:
                gene = self.geneFactory.getGene(geneName)
                e = '"%s" %s' % (term, ", ".join(gene.logicalRepresentation()))
                assert e != ''
                print(e, file=output, flush=True)
                pbar.update(pbar.currval + 1)
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

    def generateDataset(self, term, output, maxPositive = None, maxNegative = None, testRatio = 0.1):
        """Generate whole dataset directly usable for learning. The terms are used as learning classes."""
        debug("Generating dataset for term: %s" % (term))
        totalPos = len(self.associations[term]) if maxPositive is None else min(maxPositive, len(self.associations[term]))
        totalNeg = len(self.associations[term]) if maxPositive is None else min(maxNegative, len(self.associations['~'+term]))
        total = totalPos + totalNeg
        debug("We use %d postive and %d negative examples." % (totalPos, totalNeg))
        pbar = progressbar.ProgressBar(maxval=total, widgets = (
            progressbar.Bar(), ' ', progressbar.Counter(), '/'+str(total), ' =', progressbar.Percentage()))

        sampleCounts = []
        pbar.start()
        self.generateExamples(    term, pbar, output, maxPositive)
        self.generateExamples('~'+term, pbar, output, maxNegative)
        pbar.finish()

        # Calculate train set and test set indices
        #posTest = round(totalPos * testRatio)
        #negTest = round(totalNeg * testRatio)
        #testSet = [*range(posTest)] + [*range(totalPos, totalPos+negTest)]
        #trainSet = [*range(posTest, totalPos)] + [*range(totalPos+negTest, totalPos+totalNeg)]

        debug("Finished generating dataset.")

    def getTermByName(self, name):
        """Returns human-readable name of the given GO term."""
        for k,v in self.ontology.items():
            if v['name'] == name:
                return k
        raise KeyError('Name %s is not in the ontology!' % name)

    def completeTest(self, maxPositive, maxNegative, treelikerPath, template):
        bestClassifiers = []
        terms = sorted(self.ontology.keys())[:3]
        treeliker = TreeLikerWrapper(self, treelikerPath, template)
        for term in terms:
            if term == self.root:
                continue

            best = treeliker.runTermTest(term, maxPositive, maxNegative)
            bestClassifiers.append(best)
            
        #bnet = self._toBayessNet(bestClassifiers, terms)

    def _toBayessNet(self, classifiers, terms):
        assert hasattr(self, 'reserved')

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
