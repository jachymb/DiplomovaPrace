import sys
import os
from contextlib import contextmanager
import dill
import multiprocessing

verbosity = 2
def debug(s, end=True):
    if verbosity > 0:
        print(s, file=sys.stderr, end = "\n" if end else "", flush=True)


@contextmanager
def in_directory(path):
    pwd = os.path.abspath(os.curdir)
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

