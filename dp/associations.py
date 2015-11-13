__all__ = ["GeneAssociations"]

import pickle
import random
from dp.gene import Gene
from dp.utils import debug
from collections import defaultdict

class GeneAssociations:
    """This class associates GO terms with genes."""

    def __init__(self, associations, alltaxons, dataset = None):
        self.associations = associations
        self.alltaxons = alltaxons
        self.translations = {}
        self.ontology = None
        self.dataset = dataset

    @classmethod
    def fromFile(cls, inputFileName, taxons = None, dataset = None):
        """Decides file type and reads relevant data."""
        debug("Reading gene associations file %s...%s" % (inputFileName, ("" if dataset is None else " Dataset size is %d." % len(dataset))))
        #open = gzip.open if inputFileName.endswith(".gz") else __builtins__.open

        if inputFileName.endswith('.pickle') or inputFileName.endswith('.pickle_reserved'):
            # Serialized data = much faster
            with open(inputFileName, 'rb') as f:
                associations, alltaxons = pickle.load(f)
        else:
            associations = defaultdict(set)
            alltaxons = set()

            with open(inputFileName, 'rb') as associationFile:
                for line in associationFile.read().decode('utf8').splitlines():
                    if line.startswith('!'): continue
                    line = line.split('\t')
                    taxon = {int(x.split(':')[1]) for x in line[12].split('|')}
                    alltaxons.update(taxon)
                    gene = Gene.canonicalName(line[2])
                    term = line[4]
                    if (taxons is None or taxons.intersection(taxon)) and \
                       (dataset is None or gene in dataset):
                        associations[term].add(gene)
        debug("Finished reading gene associations file %s... " % inputFileName)
        #if dataset is not None:
        #    d = dataset.difference(allgenes)
        #    if d:
        #        debug("Missing genes: %s!!!" % ", ".join(d))
        return cls(associations, alltaxons, dataset)

    def transitiveClosure(self):
        """Transitive closure of associations makes genes to be associated to all parents of nodes they are currently associated to."""
        #debug("Calculating transitive closure... ", False)
        def getChildGenes(term):
            for child in self.ontology[term]['children']:
                self.associations[term].update(getChildGenes(child))
            return self.associations[term]
        getChildGenes(self.ontology.root)

        # Remove from associations terms not in ontology
        termsToDelete = [term for term in self.associations if term not in self.ontology.ontology]
        for term in termsToDelete:
            del self.associations[term]

    def __getitem__(self, item):
        """An instance of this class can be indexed by GO terms."""
        if item in self.associations:
            return self.associations[item]
        elif item in self.translations:
            term = self.translations[item]
            if term in self.associations:
                return self.associations[term]
        elif item.startswith("GO:"):
            return set()
        elif item.startswith("~"):
            return self.complement(item[1:])

        raise KeyError(item)

    def complement(self, term):
        """ Returns genes that are NOT associated with the term."""
        if term not in self.associations and term not in self.translations:
            raise KeyError(term)
        return self[self.ontology.root].difference(self[term])

    def serialize(self, fName):
        """Serializes data to a file = faster future use."""
        debug("Serializing gene associations to file %s..." % fName)
        data = (self.associations, self.alltaxons)
        with open(fName, 'wb') as f:
            pickle.dump(data, f)
        debug("Finished serializing gene associations to file %s..." % fName)

    def inverseAssoc(self, gene):
        """Returns terms associated with gene"""

        gene = Gene.canonicalName(gene)
        for term in self.associations:
            if gene in self[term]:
                yield term

    def getRatio(self, term):
        """Returns relative number of genes associated with terms compared to all genes"""
        return len(self[term]) / len(self.associations[self.ontology.root])

    def delgene(self, gene):
        terms = self.inverseAssoc(gene)
        for t in terms:
            self.associations[t].remove(gene)

    def shrink(self, toSize, minTermAssociations):
        random.seed(0)
        debug("Shrinking associations")
        allgenes = sorted(self.associations[self.ontology.root])
        size = len(allgenes)
        while size > toSize:
            todel = random.choice(allgenes)
            allgenes.remove(todel)
            self.delgene(todel)
            self.ontology.deleteSmallTerms(minTermAssociations)

            allgenes = sorted(self.associations[self.ontology.root])
            size = len(allgenes)

        self.ontology.genes = allgenes

        debug("Finished shrinking associations. Left with %d genes." % (size))
