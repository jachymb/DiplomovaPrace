# OBOSLETE STUFF

class Taxonomy:
    """This class can be used to get taxon IDs of organisms. """

    def __init__(self, inputFileName):
        self.names = {int(e[0]) : e[1] for e in ([x.strip() for x in line.split("|")] for line in open(inputFileName).read().splitlines()) if e[3] == 'scientific name'}

    def __getitem__(self, item):
        try:
            return self.names[item]
        except KeyError:
            return None


