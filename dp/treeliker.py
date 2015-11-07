import subprocess
import os
from pathlib import Path
from dp.learn import learningTest
from dp.utils import in_directory, debug
import dp
import sys

RESULTS = Path('results')
NUM_FOLDS = 8

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

        debug("Starting treeliker for "+batchPath.name)
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

        debug("Finished treeliker for "+batchPath.name)

    def runTermTest(self, term, maxPositive, maxNegative):
        term = self.ontology[term]['name']

        resultPath = RESULTS / term.replace(' ','_')
        if not resultPath.is_dir():
            resultPath.mkdir()
        batchPath = resultPath / 'batch.treeliker'

        datasetPath = resultPath / 'dataset.txt'
        with datasetPath.open('w') as output:
            self.ontology.generateDataset(term, output, maxPositive, maxNegative)

        batchFile = "set(algorithm, relf_grounding_counting)\n" \
                    "set(output_type, cv)\n" \
                    "set(output, '%s')\n" \
                    "set(examples, '%s')\n" \
                    "set(template, [%s])\n" \
                    "set(num_folds, %d)\n" \
                    "set(covered_class, '%s')\n" \
                    "work(yes)\n"

        batchFile %= (
                os.curdir,
                datasetPath.name,
                self.template,
                NUM_FOLDS,
                term)

        with batchPath.open('w') as bf:
            bf.write(batchFile)

        self._runTreeLiker(resultPath, batchPath)

        return tuple(learningTest(resultPath))
