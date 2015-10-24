import sys

def debug(s, end=True):
    print(s, file=sys.stderr, end = "\n" if end else "", flush=True)

