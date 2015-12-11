from collections import defaultdict
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
import bz2
import dill
import multiprocessing
import os
import sys

NUM_FOLDS = 8
TEST_SIZE = 0.4
RESULTS = Path('results') # I don't like this to in configuration so it's here
rerun = False

def getTermPath(term):
    # Expects term name, not GO id
    path = RESULTS / term.replace(' ','_').replace('/', '-')
    if not path.is_dir():
        path.mkdir()
    return path

def loadClf(term, fold, clfName):
    path = str(getTermPath(term) / str(fold) / (clfName + '.pickle.bz2'))
    with bz2.open(path, 'rb') as ser:
        return dill.load(ser)

verbosity = 2
def debug(s, end=True):
    s = datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' ' + s
    if verbosity > 0:
        print(s, file=sys.stderr, end = "\n" if end else "", flush=True)

# Inspired by http://stackoverflow.com/questions/8804830/python-multiprocessing-pickling-error
def apply_packed_function_for_map(args):
    dumped_function, item, args, kwargs = args
    target_function = dill.loads(dumped_function)
    res = target_function(item, *args, **kwargs)
    return res

def pack_function_for_map(target_function, items, *args, **kwargs):
    dumped_function = dill.dumps(target_function)
    dumped_items = [(dumped_function, item, args, kwargs) for item in items]
    return apply_packed_function_for_map, dumped_items

def parallel_map_dill(workers, function, iterable):
    pool = multiprocessing.Pool(processes=workers)
    return pool.imap(*pack_function_for_map(function, iterable))

def parseFasta(path):
    buf=None
    if isinstance(path, str):
        path = Path(path)
    with path.open() as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if buf is not None:
                    yield name, buf
                name=line[1:]
                buf=""
            else:
                buf += line


