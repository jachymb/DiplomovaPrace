from collections import defaultdict
from contextlib import ExitStack
from dp.bayesnet import BayesNet
from dp.gene import GeneFactory
from dp.treeliker import TreeLikerWrapper
from dp.utils import debug, parallel_map_dill, getTermPath
from itertools import groupby, product
from pathlib import Path
import dp.utils
import json
#import progressbar
import random
import subprocess
import sys
import traceback

__all__ = ["Onotology"]

class Ontology:
    """ Class representing the Gene Onotology graph."""

    def __init__(self, inputFileName, namespace = 'molecular_function'):
        """Constructor, reads and parses the ontology OBO file."""
        debug("Reading ontology file %s... " % inputFileName)
        self.root = None
        self.namespace = namespace
        ontology = defaultdict(lambda: defaultdict(list))
        self.inputFileName = Path(inputFileName)
        with self.inputFileName.open() as go:
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
                ontology[term['id']]['name'] = term['name'].replace('_', ' ') # FIXME KDYBY BLBLO, ODEBRAT replace
                for ref in term['is_a']:
                    refid, refname = ref.split(' ! ')
                    ontology[refid]['children'].append(term['id'])
                    ontology[term['id']]['parents'].append(refid)
                # This is used by Bayes nets
                ontology[term['id']]['node'] = {} # clfName : node
                ontology[term['id']]['clf'] = {} # clfName : Classifier

        self.ontology = {**ontology}

        self.associations = None
        self.geneFactory = GeneFactory()
        debug("Initialized ontology for file %s... " % inputFileName)

    def genTranslations(self):
        return {self.ontology[term]['name'] : term for term in self.ontology}

    def setAssociations(self, associations, attrname = 'associations'):
        associations.translations = self.genTranslations()
        associations.ontology = self
        self.__setattr__(attrname, associations)
        associations.transitiveClosure()

    def _termDepth(self, term):
        #if term == self.root: return 0
        parents = self.ontology[term]['parents']
        if len(parents) == 0: return 0
        return 1 + min(map(self._termDepth, parents))

    def deleteSmallTerms(self, lb):
        """ Delete terms which don't have enough associations.
        Do not call this prior to calling transitiveClosure()"""
        notEnough = self.termsByDepth(terms = (term for term in self.ontology if len(self.associations[term]) < lb))
        for term in notEnough:
            if term in self.associations.associations:
                del self.associations.associations[term]
            if hasattr(self, 'reserved') and term in self.reserved.associations:
                del self.reserved.associations[term]
                
            for d1, d2 in (('parents', 'children'),('children','parents')):
                for parent in self.ontology[term][d1]:
                    if term in self.ontology[parent][d2]:
                        self.ontology[parent][d2].remove(term)
            del self.ontology[term]
        #debug("Deleted %d terms with not enough associations. Left %d terms." % (len(notEnough), len(self.ontology)))

    def jsonExport(self):
        """Export generated onotology in JSON"""
        with (self.inputFileName.parent / (self.inputFileName.stem + '.json')).open('w') as output:
            json.dump(self.ontology, output)

    def dotExport(self):
        """ Export as a graph to the DOT format for graphical presentation."""
        debug("Exporting ontology to dot.")
        def nodename(t): return t.replace(":", "")
        for direction in ("parents", "children"):
            dotFile = self.inputFileName.parent / (self.inputFileName.stem +"_"+direction+ '.dot')
            with dotFile.open('w') as output:
                print("digraph ontology {", file=output)
                print('overlap="false";', file=output)
                print('root="%s";' % nodename(self.root), file=output)
                for term, props in self.ontology.items():
                    name = props['name']
                    if not name:
                        continue
                    if ',' in name:
                        name = name.replace(', ', r',\n')
                    else:
                        name = name.replace('binding', r'\nbinding').replace('activity',r'\nactivity').replace(r' \n',r'\n')
                    print('%s [fontsize=8,label="%s"]' % (nodename(term), name), file=output)

                for term, props in self.ontology.items():
                    for related in props[direction]:
                        print('%s -> %s' % (nodename(term), nodename(related)), file=output)
                print("}", file=output)
            for fmt in ('ps', 'png'):
                outFile = dotFile.parent / (dotFile.stem + '.' + fmt)
                #print(" ".join(['dot', '-T'+fmt, str(dotFile), '-o', str(outFile)]))
                try:
                    subprocess.Popen(['dot', '-T'+fmt, str(dotFile), '-o', str(outFile).replace(".","_dot.")]) 
                    subprocess.Popen(['twopi', '-T'+fmt, str(dotFile), '-o', str(outFile).replace(".","_twopi.")]) 
                except IOError:
                    pass
        debug("Finished dot export.")

    def overView(self):
        terms = self.termsByDepth(False)
        total = len(self.associations[self.root])
        for term in terms:
            name = self.ontology[term]['name']
            asoc = len(self.associations[term])
            ratio = asoc / total
            children = ", ".join((repr(self[ch]['name']) for ch in self.ontology[term]['children']))
            depth = self._termDepth(term)

            print("%s\t%d\t%.2f\t%d\t%s\t%s"% (term, asoc, ratio, depth, name, children), file=sys.stderr)

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

    def termsByDepth(self, leavesFirst = True, terms = None):
        if terms is None:
            terms = self.ontology.keys()
        m = -1 if leavesFirst else 1
        return sorted(terms, key = lambda x: (m*self._termDepth(x), x))

    def generateExamplesUnified(self):
        debug("Generating unified datasets.")
        terms = self.termsByDepth(False)
        #rootname = self.ontology[self.root]['name']
        with ExitStack() as stack: # Closes all files when exited
            files = [(term, stack.enter_context((getTermPath(term) / 'dataset.txt').open('w')))
                    for term
                    in (self[t]['name'] for t in self.ontology.keys())
                    ]#if term != rootname]
            #for i, geneName in enumerate(self.genes):
            for geneName in self.genes:
                #debug("%d. Writing gene %s." % (i, geneName))
                gene = self.geneFactory.getGene(geneName)
                repg = ", ".join(gene.logicalRepresentation())
                for term, output in files:
                    if geneName not in self.associations[term]:
                        term = '~'+term
                    e = '"%s" %s' % (term, repg)
                    print(e, file=output)

    #def generateExamples(self, term, output, associations, maxAssociations = None):
     #    """Generates examples from genes associated with the term in logical reprentation in the pseudo-prolog syntax. """
     #   def getRecord(geneName):
     #       try:
     #           gene = self.geneFactory.getGene(geneName)
     #           e = '"%s" %s' % (term, ", ".join(gene.logicalRepresentation()))
     #           print(e, file=output, flush=True)
     #           #if pbar is not None:
     #           #    pbar.update(pbar.currval + 1)
     #       except Exception as exc:
     #           traceback.print_exc()

     #   genes = sorted(associations[term])
     #   #genes = list(genes)
     #   random.shuffle(genes)
     #   genes = genes[:maxAssociations]
     #   for i,gene in enumerate(genes):
     #       getRecord(gene)
     #       #debug("%s %d/%d" %(self[term]['name'] if not term.startswith('~') else term,i,maxAssociations))
     #   #with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
     #       #executor.map(getRecord, genes)
            #for future in concurrent.futures.as_completed(s):
            #    debug("Future completed. " + str(future))
            #    if future.exception() is not None:
            #        debug("Exception: " + future.exception())

    #def generateDataset(self, term, output, maxPositive = None, maxNegative = None, associations = None):
    #    """Generate whole dataset directly usable for learning. The terms are used as learning classes."""
    #    # FIXME: Remove maxNegative parameter
    #    if associations is None:
    #        associations = self.associations

    #    totalPos = len(associations[term]) if maxPositive is None else min(maxPositive, len(associations[term]))
    #    maxNegative = round(totalPos / associations.getRatio(term))
    #    totalNeg = len(associations['~'+term]) if maxPositive is None else min(maxNegative, len(associations['~'+term]))
    #    total = totalPos + totalNeg
    #    termname = self[term]['name'] if not term.startswith('~') else termname
    #    debug("Generating dataset for term: %s. Using %d postive and %d negative examples." % (termname, totalPos, totalNeg))
        #if dp.utils.verbosity >= 2:
            #pbar = progressbar.ProgressBar(maxval=total, widgets = (
            #    progressbar.Bar(), ' ', progressbar.Counter(), '/'+str(total), ' =', progressbar.Percentage()))
            #pbar.start()
        #else:
            #pbar = None
    #    self.generateExamples(    term, output, associations, maxPositive)
    #    self.generateExamples('~'+term, output, associations, maxNegative)
        #if dp.utils.verbosity >= 2:
        #    pbar.finish()
    #    debug("Finished generating dataset for term: %s" % (termname))

    def getTermByName(self, name):
        """Returns human-readable name of the given GO term."""
        for k,v in self.ontology.items():
            if v['name'] == name:
                return k
        raise KeyError('Name %s is not in the ontology!' % name)

    def completeTest(self, treelikerPath, template, processes = 1):
        self.generateExamplesUnified()
        bestClassifiers = []
        terms = self.termsByDepth() # This sorting is needed later in bnet learning
        treeliker = TreeLikerWrapper(self, treelikerPath, template)
        #treeliker.runValidation(self.reserved)
        def processTerm(term):
            return term, treeliker.runTermTest(term)
        
        nets = defaultdict(dict)
        allresults = parallel_map_dill(processes, processTerm, terms)
        for term, learned in allresults:
            for clfName, folds in learned.items():

                for clf, X_train, y_train, X_test, y_test, X_validation, y_validation, g_test, g_train, g_validation in folds:
                    i = clf.fold
                    if clfName in nets[i]:
                        net = nets[i][clfName]
                    else:
                        net = BayesNet(i, clfName, self)
                        nets[i][clfName] = net
                    print(i,net.fold)

                    net.generateCPD(term, clf, X_train, y_train, X_test, y_test, X_validation, y_validation, g_train, g_test, g_validation) 

        for i, byClf in sorted(nets.items()):
            for clfName, net in byClf.items():
                net.bake()
                net.predict()
        debug("Finished complete test.")

