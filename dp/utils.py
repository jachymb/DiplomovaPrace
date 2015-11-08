import sys
from collections import defaultdict
import os
from contextlib import contextmanager
import dill
import multiprocessing
from datetime import datetime
from pathlib import Path

NUM_FOLDS = 4
TEST_SIZE = 0.4
RESULTS = Path('results') # I don't like this to in configuration so it's here

def getTermPath(term):
    # Expects term name, not GO id
    path = RESULTS / term.replace(' ','_')
    if not path.is_dir():
        path.mkdir()
    return path

verbosity = 2
def debug(s, end=True):
    s = datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' ' + s
    if verbosity > 0:
        print(s, file=sys.stderr, end = "\n" if end else "", flush=True)

@contextmanager
def in_directory(path):
    pwd = os.path.abspath(os.curdir)
    path = Path(path)
    if not path.is_dir():
        path.mkdir()
    os.chdir(str(path))
    yield
    os.chdir(pwd)

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
    return pool.map(*pack_function_for_map(function, iterable))

def flipdict(d):
    new = defaultdict(dict)
    for k, v in d.items():
        for s, r in v.items():
            new[s][k] = r
    return new

