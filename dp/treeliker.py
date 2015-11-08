import subprocess
import os
from pathlib import Path
from sklearn import cross_validation
from dp.learn import learningTest
from dp.utils import in_directory, debug, RESULTS, getTermPath, NUM_FOLDS, TEST_SIZE
import dp
import sys


__all__ = ["TreeLikerWrapper"]

class TreeLikerWrapper:
    maxMemory = None
    def __init__(self, ontology, treeliker, template):
        self.ontology = ontology
        self.treeliker = str(Path(treeliker).resolve())
        self.template = template

    def runValidation(self, dataset):
        validationPath = RESULTS / 'validation.txt'
        with validationPath.open('w') as output:
            self.ontology.generateDataset(self.ontology.root, output, associations=dataset, maxPositive = len(dataset[self.ontology.root]))

        arff = 'validation.arff'

        batchPath = RESULTS / 'validation.treeliker'
        batchFile = "set(algorithm, relf_grounding_counting)\n" \
                    "set(output_type, single)\n" \
                    "set(output, '%s')\n" \
                    "set(examples, '%s')\n" \
                    "set(template, [%s])\n" \
                    "work(yes)\n"
                    # covered_class shouldn't be used here!
 
        batchFile %= (
                arff,
                validationPath.name,
                self.template)

        with batchPath.open('w') as bf:
            bf.write(batchFile)

        self._runTreeLiker(RESULTS, batchPath)

        return RESULTS / arff

    def _runTreeLiker(self, resultPath, batchPath):
        cmd = ["java", "-cp", self.treeliker, "ida.ilp.treeLiker.TreeLikerMain", "-batch", batchPath.name]
        if self.maxMemory is not None:
            cmd.insert(1, '-Xmx'+self.maxMemory)

        debug("Starting treeliker for "+resultPath.name)
        with in_directory(resultPath), subprocess.Popen(cmd, stdout = subprocess.PIPE, bufsize = 1, universal_newlines=True) as treelikerProc:
            prev = 0
            i = 1
            for _line in treelikerProc.stdout:
                line = '\r%d : %s' % (i, _line.rstrip())
                if _line.startswith('Fold') and dp.utils.verbosity == 1:
                    debug("%s: %s" % (batchPath.name, _line))
                elif dp.utils.verbosity >= 2:
                    debug(line.ljust(prev), end=_line.startswith('Fold'))
                prev = len(line)
                if _line.startswith('Processing'):
                    i+=1
        if dp.utils.verbosity >= 2:
            sys.stderr.write("\n")

        debug("Finished treeliker for "+resultPath.name)

    def runTermTest(self, term):
        term = self.ontology[term]['name']
        debug("Preparing for TreeLiker on term %s." % term)

        resultPath = getTermPath(term)
        batchPath = resultPath / 'batch.treeliker'

        datasetPath = resultPath / 'dataset.txt'

        batchFile = "set(algorithm, relf_grounding_counting)\n" \
                    "set(output_type, train_test)\n" \
                    "set(examples, '%s')\n" \
                    "set(template, [%s])\n" \
                    "set(covered_class, '%s')\n\n" % (
                        datasetPath.name,
                        self.template,
                        term)

        with datasetPath.open() as ds:
            dataSetLen = len([*ds]) # Counts lines

        for i, (train, test) in enumerate(cross_validation.KFold(dataSetLen, NUM_FOLDS)):
            path = resultPath / str(i)
            if not path.is_dir():
                path.mkdir()
                
            batchFile += "set(output, '%s')\n" \
                         "set(train_set, [%s])\n" \
                         "set(test_set, [%s])\n" \
                         "work(yes)\n" % (
                             path.name,
                             ",".join(map(str,train)),
                             ",".join(map(str,test)))

        with batchPath.open('w') as bf:
            bf.write(batchFile)

        self._runTreeLiker(resultPath, batchPath)

        return tuple(learningTest(resultPath))
