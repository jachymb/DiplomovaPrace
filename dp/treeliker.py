import subprocess
import os
from pathlib import Path
from dp.learn import learningTest

RESULTS = Path('results')
NUM_FOLDS = 8

class TreeLikerWrapper:
    def __init__(self, ontology, treeliker, template):
        self.ontology = ontology
        self.treeliker = str(Path(treeliker).resolve())
        self.template = template

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

        pwd = os.path.abspath(os.curdir)
        os.chdir(str(resultPath))

        with subprocess.Popen(
                ["java", "-Xmx1G", "-cp", self.treeliker, "ida.ilp.treeLiker.TreeLikerMain", "-batch", batchPath.name],
                stdout = subprocess.PIPE, bufsize = 1, universal_newlines=True) as treelikerProc:
            prev = 0
            i = 1
            for _line in treelikerProc.stdout:
                line = '\r%d : %s' % (i, _line.rstrip())
                print(line.ljust(prev), end="\n" if _line.startswith('Fold') else "", flush=True)
                prev = len(line)
                if _line.startswith('Processing'):
                    i+=1
        print()
        os.chdir(pwd)

        return learningTest(resultPath)
