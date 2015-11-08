#!/usr/bin/env python
import sys
import os

import random
random.seed(0)

from optparse import OptionParser
from dp.ontology import Ontology
from dp.associations import GeneAssociations
from dp.gene import Gene
from dp.treeliker import TreeLikerWrapper
import dp.utils
#TAXONS_SOLANCEAE = {13442,243963,24663,4071,4074,4078,4081,4085,4101,4111,4113,448155}
#TAXONS_SACCHAROMYCES_CEREVISIAE = {559292}
#TAXONS_SACCHAROMYCES_CEREVISIAE = {4932,11008,1182966,1182967,1182968,1247190,1293430,1294303,1294306,1294310,1294311,1294312,1294313,1294314,1294315,1294316,1294317,1294318,1294319,1294320,1294321,1294322,1294323,1294325,1294326,1294327,1294330,1294332,1294333,1294334,1294335,1294336,1294337,1294338,1294339,1294342,1294343,1294344,1294345,1294346,1294348,1294349,1294350,1294351,1294354,1294355,1294356,1294358,1294360,1294361,1294363,1294364,1294369,1294371,1294373,1294374,1294377,1294378,1294379,1294380,1294381,1294382,1294383,1294384,1294387,1294388,1296266,468558,502869,574961,580240,643680,721032,764097,764101,889517,947042,947044}
TAXONS_HOMMO_SAPIENS = {9606}
#TAXONOMY_NAMES = 'names.dmp'
TAXONS = None
TAXONS = TAXONS_HOMMO_SAPIENS

os.system("taskset -p 0xff %d" % os.getpid())
if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.usage = "%prog [options] OBO-FILE ASSOCIATIONS-FILE"
    parser.add_option("-a", "--dump-associations", dest="associationsDump", help="Dump serialized associations to file specified by this parameter and exit")
    parser.add_option("-c", "--processes", dest="processes", type="int", help="Maximum spawned subprocesses.", default=1)
    parser.add_option("-b", "--background-knowledge", dest="backgroundKnowledge", help="File with background knowledge appended to each gene.")
    parser.add_option("-d", "--dataset", dest="dataset", help="Set of genes to use.")
    parser.add_option("-e", "--memory", dest="memory", help="Maximum TreeLiker memory.")
    parser.add_option("-i", "--recalculate-distances", dest="recalcDists", help="Reacalculate distances (if max distance or discreatisation were changed. Coordinates are used the same.)", action="store_true")
    parser.add_option("-l", "--min-associations", dest="lb", type="int", help="Minimum number of posive examples. (If there are not enough, term is ignored)", default=1)
    #parser.add_option("-m", "--max-positive", dest="max_positive", type="int", help="Maximum positive samples.")
    #parser.add_option("-n", "--max-negative", dest="max_negative", type="int", help="Maximum negative samples.")
    parser.add_option("-m", "--max", dest="max", type="int", default=2048, help="Maximum dataset size.")
    parser.add_option("-p", "--template", dest="template", help="The template for TreeLiker. Use with -f.")
    parser.add_option("-r", "--reserve", dest="reserve", type="float", help="Ratio of genes (or absolute number if >= 1) to reserve for Bayessian learning.", default = 1024)
    parser.add_option("-s", "--no-deserialize", dest="deserialize", help="Don't deserialize stored gene data.", action="store_false")
    parser.add_option("-v", "--verbosity", dest="verbosity", type="int", default=2, help="0 = Silent, 1 = Hide dynamic elements, 2  = Show everything")

    parser.add_option("-x", "--treeliker", dest="treeliker", help="TreeLiker jar binary.")

    options, args = parser.parse_args()

    if len(args) != 2:
        parser.error("Incorect number of arguments!")

    oboFileName, associationsFileName = args

    dp.utils.verbosity = options.verbosity

    if options.backgroundKnowledge:
        with open(options.backgroundKnowledge) as bk:
            Gene.backgroundKnowledge = bk.read().splitlines()

    ontology = Ontology(oboFileName)
    ontology.geneFactory.deserialize = True if options.deserialize is None else False
    #associations = GeneAssociations(associationsFileName, TAXONS_SACCHAROMYCES_CEREVISIAE)

    dataset = None
    if options.dataset:
        # FIXME: When dataset is changed, serialized associations need to be regenerated. This is serious bug if we don't seed random
        dataset = [*open(options.dataset).read().splitlines()]
        random.shuffle(dataset)
        assert options.reserve > 0.0
        if options.reserve < 1.0: # Use ratio
            splitIndex = int(options.reserve * len(dataset))
        else:
            splitIndex = int(options.reserve)
        reserved = set(dataset[:splitIndex])
        #dataset = set(dataset[splitIndex:])
        dataset = set(dataset)

    associations = GeneAssociations.fromFile(associationsFileName, taxons = TAXONS, dataset = dataset)
    #reservedAssociations = GeneAssociations.fromFile(associationsFileName+"_reserved", dataset = reserved)
    ontology.setAssociations(associations)
    #ontology.setAssociations(reservedAssociations, 'reserved')
 
    if options.associationsDump:
        associations.serialize(options.associationsDump)
        #reservedAssociations.serialize(options.associationsDump+"_reserved")
        sys.exit()

    ontology.deleteSmallTerms(options.lb)
    associations.shrink(options.max, options.lb)
    
    ontology.overView()
    ontology.dotExport()


    #ontology._toBayessNet(None,None)
    TreeLikerWrapper.maxMemory = options.memory
    Gene.recalculateDists = bool(options.recalcDists)

    ontology.completeTest(options.treeliker, options.template, options.processes)
